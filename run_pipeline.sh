#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

CORES=20
PROJECT_NR=""
SAMPLES=""
SAMPLES_GIVEN=0
PROJECT_NR_GIVEN=0
DIAMOND_DB="/path/to/DiamondDB.dmnd"
VERBOSE=0
DRYRUN=0
BINNING_BACKEND="mycobinr"
GNN_INPUT_ROOT="PacBio_output/gnn_inputs"
GNN_PYTHON="python3"
GENE_PREDICTION="busco"
TIARA_MODE="off"
TIARA_MODE_SET=0
TIARA_MODE_SET=0
TIARA_TSV=""
METAEUK_REF=""
METAEUK_SIF="bin/singularity/metaeuk.sif"
PRODIGAL_SIF="bin/singularity/prodigal.sif"
TIARA_SIF="bin/singularity/tiara.sif"
CONTAINER_RUNTIME="auto"
CONTAINER_IMAGES="bin/containers/images.json"
SETUP_CONTAINERS=1
declare -a READS=() SAMPLE_IDS=() NAMED_SAMPLES=()

usage() {
  cat <<'EOF'
Usage: run_pipeline.sh [options]

Input choices (do not mix):
  --reads FILE [FILE ...]     FASTQ paths; sample IDs derived from filenames
  --sample-id ID             Explicit ID for each --reads file, same order
  --sample ID=FILE           Explicit ID and FASTQ (repeat for many samples)
  --samples LIST             Legacy IDs, e.g. 001,002 (looks in PacBio_data/)
  --project_nr PREFIX        Output prefix (default: empty with direct reads;
                             nr_ with legacy --samples)

Other options:
  --cores N                  Number of CPU cores (default: 20)
  --db PATH                  DIAMOND .dmnd file (required)
  --gene-prediction MODE     busco (default) or standalone
  --metaeuk-ref PATH         MetaEuk reference FASTA (standalone)
  --metaeuk-sif PATH         Optional override for MetaEuk SIF
  --prodigal-sif PATH        Optional override for Prodigal SIF
  --tiara-mode MODE          off (default), run, or file
  --tiara-tsv PATH           Tiara TSV if mode=file (one sample only)
  --tiara-sif PATH           Optional override for Tiara SIF if mode=run
  --container-runtime MODE  auto (default), docker, apptainer, singularity, none
  --binning-backend MODE   mycobinr (default), gnn, both
  --gnn-input-root DIR    Existing ML-compatible input root (per sample/input)
  --gnn-python PATH       Python executable with ML dependencies
                            (GNN enables Tiara by default; use --tiara-mode off to override)
                            (GNN enables Tiara by default; use --tiara-mode off to override)
  --container-images PATH   JSON manifest for images (default bin/containers/images.json)
  --no-container-setup      Do not provision images before actual execution
  -n, --dry-run              Snakemake dry-run
  -v, --verbose              Verbose shell output
  -h, --help                 Show this help

Examples:
  bash run_pipeline.sh --reads /data/pr_197_001.hifireads.fastq.gz --db /db/nr.dmnd -n
  bash run_pipeline.sh --sample lichenA=/data/x.fastq.gz --sample lichenB=/data/y.fq.gz --db /db/nr.dmnd -n
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cores)        CORES="${2:?}"; shift 2;;
    --project_nr)   PROJECT_NR="${2-}"; PROJECT_NR_GIVEN=1; shift 2;;
    --samples)      SAMPLES="${2-}"; SAMPLES_GIVEN=1; shift 2;;
    --reads)
      shift
      if [[ $# -eq 0 || "$1" == -* ]]; then
        echo 'Error: --reads requires at least one FASTQ path' >&2; exit 2
      fi
      while [[ $# -gt 0 && "$1" != -* ]]; do
        READS+=("$1"); shift
      done;;
    --sample-id)    SAMPLE_IDS+=("${2:?}"); shift 2;;
    --sample)       NAMED_SAMPLES+=("${2:?}"); shift 2;;
    --db|--diamond_db) DIAMOND_DB="${2:?}"; shift 2;;
    --gene-prediction) GENE_PREDICTION="${2:?}"; shift 2;;
    --metaeuk-ref) METAEUK_REF="${2:?}"; shift 2;;
    --metaeuk-sif) METAEUK_SIF="${2:?}"; shift 2;;
    --prodigal-sif) PRODIGAL_SIF="${2:?}"; shift 2;;
    --tiara-mode) TIARA_MODE="${2:?}"; TIARA_MODE_SET=1; shift 2;;
    --tiara-tsv) TIARA_TSV="${2:?}"; shift 2;;
    --tiara-sif) TIARA_SIF="${2:?}"; shift 2;;
    --container-runtime|--runtime) CONTAINER_RUNTIME="${2:?}"; shift 2;;
    --binning-backend) BINNING_BACKEND="${2:?}"; shift 2;;
    --gnn-input-root) GNN_INPUT_ROOT="${2:?}"; shift 2;;
    --gnn-python) GNN_PYTHON="${2:?}"; shift 2;;
    --container-images) CONTAINER_IMAGES="${2:?}"; shift 2;;
    --no-container-setup) SETUP_CONTAINERS=0; shift;;
    -n|--dry-run)   DRYRUN=1; shift;;
    -v|--verbose)   VERBOSE=1; shift;;
    -h|--help)      usage; exit 0;;
    --) shift; break;;
    *) echo "Unknown parameter: $1" >&2; usage; exit 2;;
  esac
done
[[ "$VERBOSE" -eq 1 ]] && set -x

[[ "$CORES" =~ ^[0-9]+$ && "$CORES" -ge 1 ]] || { echo 'Error: --cores must be a positive integer' >&2; exit 2; }
[[ -f "$DIAMOND_DB" ]] || { echo "Error: DIAMOND DB not found: $DIAMOND_DB" >&2; exit 2; }
DIAMOND_DB="$(realpath "$DIAMOND_DB")"
[[ "${DIAMOND_DB##*.}" == dmnd ]] || echo "Warning: not a .dmnd file: $DIAMOND_DB" >&2
[[ "$GENE_PREDICTION" == busco || "$GENE_PREDICTION" == standalone ]] || { echo 'Error: --gene-prediction must be busco or standalone' >&2; exit 2; }
[[ "$TIARA_MODE" == off || "$TIARA_MODE" == run || "$TIARA_MODE" == file ]] || { echo 'Error: --tiara-mode must be off, run or file' >&2; exit 2; }

if [[ "$TIARA_MODE" == file ]]; then
  [[ -f "$TIARA_TSV" ]] || { echo "Error: Tiara TSV not found: $TIARA_TSV" >&2; exit 2; }
  TIARA_TSV="$(realpath "$TIARA_TSV")"
fi
if [[ "$TIARA_MODE" == file ]]; then
  [[ -f "$TIARA_TSV" ]] || { echo "Error: Tiara TSV not found: $TIARA_TSV" >&2; exit 2; }
  TIARA_TSV="$(realpath "$TIARA_TSV")"
fi
if [[ "$GENE_PREDICTION" == standalone ]]; then
  [[ -f "$METAEUK_REF" ]] || { echo 'Error: --metaeuk-ref must name an existing reference FASTA' >&2; exit 2; }
  METAEUK_REF="$(realpath "$METAEUK_REF")"
  if [[ "$TIARA_MODE" == file ]]; then
    [[ -f "$TIARA_TSV" ]] || { echo 'Error: --tiara-tsv must exist' >&2; exit 2; }
    TIARA_TSV="$(realpath "$TIARA_TSV")"
  fi
fi

# Resolve the requested container engine (Snakemake itself has no Docker
# --use-docker flag; Docker rule commands use docker_exec.py).
[[ "$BINNING_BACKEND" == mycobinr || "$BINNING_BACKEND" == gnn || "$BINNING_BACKEND" == both ]] || { echo "Invalid --binning-backend" >&2; exit 2; }
# Original GNN uses Tiara as a feature: provide it unless explicitly disabled.
if [[ "$TIARA_MODE_SET" -eq 0 && "$BINNING_BACKEND" != mycobinr ]]; then
  TIARA_MODE=run
fi
# Original GNN uses Tiara as a feature: provide it unless explicitly disabled.
if [[ "$TIARA_MODE_SET" -eq 0 && "$BINNING_BACKEND" != mycobinr ]]; then
  TIARA_MODE=run
fi
case "$CONTAINER_RUNTIME" in
  auto)
    if command -v apptainer >/dev/null 2>&1; then
      CONTAINER_RUNTIME=apptainer
    elif command -v singularity >/dev/null 2>&1; then
      CONTAINER_RUNTIME=singularity
    elif command -v docker >/dev/null 2>&1; then
      CONTAINER_RUNTIME=docker
    else
      echo 'Error: no Docker, Apptainer or Singularity executable in PATH' >&2; exit 2
    fi;;
  docker|apptainer|singularity|none) ;;
  *) echo 'Error: --container-runtime must be auto, docker, apptainer, singularity or none' >&2; exit 2;;
esac
[[ -f "$CONTAINER_IMAGES" ]] || { echo "Error: image manifest not found: $CONTAINER_IMAGES" >&2; exit 2; }

# Preserve the historical prefix only for legacy PacBio_data lookup.
if [[ "$PROJECT_NR_GIVEN" -eq 0 ]]; then
  if [[ "${#READS[@]}" -eq 0 && "${#NAMED_SAMPLES[@]}" -eq 0 ]]; then
    PROJECT_NR='nr_'
  else
    PROJECT_NR=''
  fi
fi

# A real JSON file avoids Snakemake 9 YAML coercion of sample 001 to integer 1.
INPUT_CFG="$(mktemp "${TMPDIR:-/tmp}/schwammerl-reads.XXXXXXXX.json")"
trap 'rm -f "$INPUT_CFG"' EXIT
sample_args=(--output "$INPUT_CFG" --project-nr "$PROJECT_NR")
if [[ "$SAMPLES_GIVEN" -eq 1 ]]; then
  sample_args+=(--samples "$SAMPLES" --legacy-samples-given)
fi
for f in "${READS[@]}"; do sample_args+=(--reads "$f"); done
for id in "${SAMPLE_IDS[@]}"; do sample_args+=(--sample-id "$id"); done
for pair in "${NAMED_SAMPLES[@]}"; do sample_args+=(--sample "$pair"); done

echo 'Input samples:'
python3 bin/scripts/prepare_samples.py "${sample_args[@]}"

if [[ "$TIARA_MODE" == file ]]; then
  sample_count="$(python3 -c 'import json,sys; print(len(json.load(open(sys.argv[1]))["samples"]))' "$INPUT_CFG")"
  [[ "$sample_count" -eq 1 ]] || { echo 'Error: --tiara-mode=file currently supports only one sample' >&2; exit 2; }
fi


# Bind the original external FASTQs and databases, without copying or linking.
# realpath() in prepare_samples.py also handles files reached through symlinks.
declare -a BIND_DIRS=("$(dirname "$DIAMOND_DB")")
while IFS= read -r d; do
  [[ -n "$d" ]] && BIND_DIRS+=("$d")
done < <(python3 - "$INPUT_CFG" <<'PY'
import json, os, sys
cfg = json.load(open(sys.argv[1]))
for path in cfg.get('reads', {}).values():
    print(os.path.dirname(path))
PY
)
if [[ "$GENE_PREDICTION" == standalone ]]; then
  BIND_DIRS+=("$(dirname "$METAEUK_REF")")
fi
[[ "$TIARA_MODE" == file ]] && BIND_DIRS+=("$(dirname "$TIARA_TSV")")
# Deduplicate bind roots without changing their order.
declare -a UNIQUE_DIRS=()
declare -A SEEN_BIND=()
for d in "${BIND_DIRS[@]}"; do
  if [[ ! -v SEEN_BIND["$d"] ]]; then
    UNIQUE_DIRS+=("$d")
    SEEN_BIND["$d"]=1
  fi
done
BINDS="$(IFS=,; echo "${UNIQUE_DIRS[*]}")"

declare -a SNK_ARGS=(--cores "$CORES" --snakefile ./bin/Snakefile
                     --printshellcmds --rerun-incomplete --latency-wait 60
                     --configfile "$INPUT_CFG")
if [[ "$CONTAINER_RUNTIME" == apptainer || "$CONTAINER_RUNTIME" == singularity ]]; then
  SNK_ARGS+=(--use-apptainer --apptainer-args "--bind $BINDS")
fi
declare -a CFG=("binning_backend=$BINNING_BACKEND" "gnn_input_root=$GNN_INPUT_ROOT" "gnn_python=$GNN_PYTHON" "diamond_db=$DIAMOND_DB" "gene_prediction=$GENE_PREDICTION" "tiara_mode=$TIARA_MODE"
                "container_runtime=$CONTAINER_RUNTIME" "container_images=$CONTAINER_IMAGES")
if [[ "$GENE_PREDICTION" == standalone ]]; then
  CFG+=("metaeuk_ref=$METAEUK_REF" "metaeuk_sif=$METAEUK_SIF" "prodigal_sif=$PRODIGAL_SIF" "tiara_sif=$TIARA_SIF")
fi
[[ "$TIARA_MODE" == file ]] && CFG+=("tiara_tsv=$TIARA_TSV")

printf 'Running Snakemake: cores=%s, prediction=%s, Tiara=%s, runtime=%s, output prefix=%q\n' \
       "$CORES" "$GENE_PREDICTION" "$TIARA_MODE" "$CONTAINER_RUNTIME" "$PROJECT_NR"
if [[ "$DRYRUN" -eq 0 && "$SETUP_CONTAINERS" -eq 1 && "$CONTAINER_RUNTIME" != none ]]; then
  # Build/pull only images used in the selected workflow. Existing ones are cached.
  required_tools=(longqc flye coverm busco diamond telfinder report)
  if [[ "$GENE_PREDICTION" == standalone ]]; then
    required_tools+=(metaeuk prodigal)
  fi
  [[ "$TIARA_MODE" == run ]] && required_tools+=(tiara)
  bash bin/setup.sh --runtime "$CONTAINER_RUNTIME" --images "$CONTAINER_IMAGES" \
    --tools "${required_tools[@]}"
fi
[[ "$DRYRUN" -eq 1 ]] && SNK_ARGS+=(-n)
snakemake "${SNK_ARGS[@]}" --config "${CFG[@]}"

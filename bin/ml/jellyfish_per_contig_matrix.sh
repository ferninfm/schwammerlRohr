#!/usr/bin/env bash
set -euo pipefail

# -------------------------
# SETTINGS (anpassen)
# -------------------------
SIF="/media/HDD4/fernando/01_maximilian/05_apptainer/container/jellyfish.sif"
IN_FASTA="/media/HDD1/fernando/01_Paul/Tests/Test_final/schwammerlRohr/PacBio_output/flye_output/pr_197_002_flye/assembly.fasta"
OUTDIR="/media/HDD4/fernando/01_maximilian/09_jellyfish/output/pr_197_002"
K=4
THREADS=4
HASH_SIZE="1M"       # für k=4 reicht 1M locker
CANONICAL="no"       # "yes" würde -C setzen (dann NICHT mehr 256 Spalten)

# Apptainer cache/tmp (optional)
export APPTAINER_CACHEDIR="/media/HDD4/fernando/01_maximilian/05_apptainer/apptainer_cache"
export APPTAINER_TMPDIR="/media/HDD4/fernando/01_maximilian/05_apptainer/apptainer_tmp"

# -------------------------
# CHECKS
# -------------------------
[[ -f "$SIF" ]] || { echo "ERROR: SIF not found: $SIF" >&2; exit 1; }
[[ -f "$IN_FASTA" ]] || { echo "ERROR: FASTA not found: $IN_FASTA" >&2; exit 1; }

mkdir -p "$OUTDIR"
WORKDIR="$OUTDIR/work"
CONTIGDIR="$WORKDIR/contigs"
mkdir -p "$CONTIGDIR"

MATRIX_CSV="$OUTDIR/assembly.k${K}.per_contig_matrix.csv"
KMERS_TXT="$WORKDIR/kmers_k${K}.txt"

# -------------------------
# 1) Generate all kmers (4^K) in fixed order
# -------------------------
gen_kmers() {
  local k="$1"
  if [[ "$k" -eq 0 ]]; then
    echo ""
    return
  fi
  local smaller
  while IFS= read -r smaller; do
    for b in A C G T; do
      echo "${b}${smaller}"
    done
  done < <(gen_kmers $((k-1)))
}

echo "Generating k-mer list (k=$K) ..."
gen_kmers "$K" > "$KMERS_TXT"

# -------------------------
# 2) Split FASTA into single-contig FASTA files
# -------------------------
echo "Splitting FASTA into per-contig files ..."
awk -v out="$CONTIGDIR" '
  BEGIN { file=""; }
  /^>/ {
    name=$0; sub(/^>/,"",name);
    gsub(/[ \t\r\n]/,"_",name);
    gsub(/[^A-Za-z0-9_.-]/,"_",name);
    file=out "/" name ".fa";
    print $0 > file;
    next;
  }
  {
    if (file != "") print $0 > file;
  }
' "$IN_FASTA"

# -------------------------
# 3) Write CSV header
# -------------------------
{
  printf "contig"
  while IFS= read -r kmer; do
    printf ",%s" "$kmer"
  done < "$KMERS_TXT"
  printf "\n"
} > "$MATRIX_CSV"

# -------------------------
# 4) For each contig: jellyfish count -> dump -> row in matrix
# -------------------------
echo "Counting per-contig k-mers with Jellyfish ..."
COUNT_ARGS=(jellyfish count -m "$K" -s "$HASH_SIZE" -t "$THREADS" -o)
if [[ "$CANONICAL" == "yes" ]]; then
  COUNT_ARGS+=( -C )
fi

shopt -s nullglob
n=0

for fa in "$CONTIGDIR"/*.fa; do
  n=$((n+1))
  contig="$(basename "$fa" .fa)"
  jf="$WORKDIR/${contig}.jf"
  dump="$WORKDIR/${contig}.dump.tsv"

  # count
  apptainer exec --bind /media:/media "$SIF" \
    "${COUNT_ARGS[@]}" "$jf" "$fa" >/dev/null

  # dump counts (kmer \t count)
  apptainer exec --bind /media:/media "$SIF" \
    jellyfish dump -c "$jf" > "$dump"

  # build one CSV row in fixed kmer order
  awk -v contig="$contig" -v kmers="$KMERS_TXT" '
    BEGIN{
      # read kmers order
      while((getline k < kmers) > 0){
        order[++m]=k
      }
      close(kmers)
    }
    {
      cnt[$1]=$2
    }
    END{
      printf "%s", contig
      for(i=1;i<=m;i++){
        k=order[i]
        printf ",%s", (k in cnt ? cnt[k] : 0)
      }
      printf "\n"
    }
  ' "$dump" >> "$MATRIX_CSV"

  # optional cleanup per contig (spart Platz)
  rm -f "$jf" "$dump"

  if (( n % 50 == 0 )); then
    echo "  processed $n contigs ..."
  fi
done

echo "Done. Matrix written to:"
echo "  $MATRIX_CSV"

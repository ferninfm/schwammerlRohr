#!/bin/bash
set -euo pipefail

# ---- Settings ----
SAMPLE="pr_197_002"
BIN="complete"

# Input assembly (Canu contigs)
IN_FASTA="/media/HDD1/fernando/01_Paul/Tests/Test_final/schwammerlRohr/PacBio_output/flye_output/pr_197_002_flye/assembly.fasta"

# Output folder for Tiara + split
OUTDIR="/media/HDD4/fernando/01_maximilian/06_tiara/output/${SAMPLE}/${BIN}"

# Apptainer + Tiara container
TIARA_SIF="/media/HDD4/fernando/01_maximilian/05_apptainer/container/tiara.sif"

# Apptainer cache/tmp (gegen Home-Quota)
export APPTAINER_CACHEDIR="/media/HDD4/fernando/01_maximilian/05_apptainer/apptainer_cache"
export APPTAINER_TMPDIR="/media/HDD4/fernando/01_maximilian/05_apptainer/apptainer_tmp"


if [ -d "$OUTDIR" ]; then
  echo "Removing old Tiara output: $OUTDIR"
  rm -rf "$OUTDIR"
fi


mkdir -p "$OUTDIR"

echo "SAMPLE:   $SAMPLE"
echo "BIN:      $BIN"
echo "IN_FASTA: $IN_FASTA"
echo "OUTDIR:   $OUTDIR"
echo "TIARA:    $TIARA_SIF"
echo


# ---- Run Tiara ----

OUTFILE="$OUTDIR/${BIN}.classification.tsv"

echo "[1/2] Running Tiara..."
apptainer exec \
  --bind /media:/media \
  "$TIARA_SIF" \
  tiara \
    -i "$IN_FASTA" \
    -o "$OUTFILE"

echo "[1/2] Tiara finished."
echo

# ---- Find classification file ----
CLASS_FILE="$OUTFILE"
echo "Classification file: $CLASS_FILE"
echo

# ---- Split IDs ----
echo "[2/2] Splitting contigs into eukaryota / prokaryota..."
awk 'BEGIN{FS="\t"} $2=="eukarya"{split($1,a," "); print a[1]}' "$CLASS_FILE" > "$OUTDIR/eukarya.ids"
awk 'BEGIN{FS="\t"} ($2=="bacteria" || $2=="archaea" || $2=="unknown"){split($1,a," "); print a[1]}' "$CLASS_FILE" > "$OUTDIR/prokaryota.ids"

# ---- Split FASTA (samtools) ----
samtools faidx "$IN_FASTA"

samtools faidx "$IN_FASTA" $(cat "$OUTDIR/eukarya.ids") \
  > "$OUTDIR/${BIN}_euk.fasta" || true

samtools faidx "$IN_FASTA" $(cat "$OUTDIR/prokaryota.ids") \
  > "$OUTDIR/${BIN}_prok.fasta" || true

echo "Done."
echo "Euk FASTA:  $OUTDIR/${BIN}_euk.fasta"
echo "Prok FASTA: $OUTDIR/${BIN}_prok.fasta"
echo

echo "Counts:"
echo -n "  euk contigs:  "; grep -c "^>" "$OUTDIR/${BIN}_euk.fasta" 2>/dev/null || echo 0
echo -n "  prok contigs: "; grep -c "^>" "$OUTDIR/${BIN}_prok.fasta" 2>/dev/null || echo 0

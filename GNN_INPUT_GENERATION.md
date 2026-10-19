# Generating missing GNN inputs (preview update)

This updates the **unified preview** rather than the historical repository checkout. It does not alter the GNN's Python implementation.

## Additions
- Tiara is automatically selected when `--binning-backend gnn|both`, unless `--tiara-mode off` is given explicitly. Its **existing** Snakemake rule generates the classification TSV; also available in BUSCO prediction mode.
- `gnn_kmer_counts` generates the same full 256 noncanonical 4-mer count columns expected by the original ML loader, without depending on the former hardcoded Jellyfish shell script.
- `gnn_bam_links` uses the existing `bin/ml/bam_links.py` with the CoverM BAM cache. It **fails clearly** if a single BAM cannot be found; this step has not yet been validated with a real CoverM run.
- `gnn_prepare_inputs` stages assembly, coverage, kmers, links, Tiara, and original BUSCO-DIAMOND taxonomic evidence using symlinks. `original_gnn` then invokes unchanged GNN source.

## Remaining limitations
- `bam_links.py` requires pysam in the Python environment running Snakemake rules. Its existing link-detection algorithm requires suitable supplementary alignments.
- CoverM may not retain the BAM required by that algorithm. If absent, implement an explicit, persistent all-read alignment with minimap2 + samtools; do **not** silently substitute an unsuitable BAM.
- BUSCO-DIAMOND is still the taxonomy source for the first GNN integration. Standalone MetaEuk/Prodigal classifications are not yet merged into it.
- The GNN Python environment, containerization, `rule all` selection, and end-to-end execution have not been fully tested.
- The current `rule all` still includes R reporting even in GNN-only mode.

## Command example (dry-run)
```bash
bash run_pipeline.sh --cores 16 --reads /path/to/pr_197_001.hifireads.fastq.gz \
  --db /path/to/nr.dmnd --container-runtime docker \
  --binning-backend gnn --gnn-python /path/to/gnn-env/bin/python --dry-run
```

Do not merge into `main` until the cache-BAM/link rule and full DAG pass on your system.

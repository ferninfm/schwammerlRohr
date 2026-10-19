# Unified schwammerlRohr preview (from uploaded source snapshots)

This directory includes the **unchanged original** `Cluster_ML/*.py` model files within `bin/ml/Cluster_ML/`, plus original auxiliary ML scripts. It uses a small wrapper and an optional `gnn.smk` rule. No separate Metagenome-Binning-ML checkout is required.

**Important limit**: The current GNN rule consumes a previously prepared `--gnn-input-root/<sample>/input/` containing `assembly.fasta`, `mean_coverage.tsv`, `assembly.k4.per_contig_matrix.csv`, and `links.tsv` (and, when needed, `complete.classification.tsv` and `<sample>_taxprot`). It does **not yet generate** all these from the raw read stage. This preview is a single repository and Snakemake workflow, but it is **not yet a validated raw-FASTQ-to-GNN pipeline**. No GNN model parameters or taxonomy features were deliberately switched off.

Command (once input folder has been prepared):

```bash
bash run_pipeline.sh --reads /path/to/sample.hifireads.fastq.gz \
 --db /path/to/nr.dmnd --container-runtime docker \
 --binning-backend both --gnn-input-root /path/to/cluster_data \
 --gnn-python /path/to/gnn-env/bin/python --dry-run
```

`--binning-backend gnn` currently **adds** the GNN output; the existing report is not suppressed. For strict backend selection, the original report targets in `rule all` need conditional handling. The scientific model remains exactly as uploaded.

### Integration work remaining
1. Connect k-mer matrix and links as explicitly declared Snakemake outputs or use existing upstream files.
2. Ensure retained BAMs if links are built from alignments.
3. Ensure original classification TSV formats/paths are supplied; do not silently run without taxonomy.
4. Test the real Snakemake DAG, GNN Python environment, and a full sample before merge into main.

This bundle is a *source preview*, not a Git push or an end-to-end validated release. Copy/review changes against your live `feat/standalone-gene-prediction` branch, which may differ from the uploaded source snapshot.

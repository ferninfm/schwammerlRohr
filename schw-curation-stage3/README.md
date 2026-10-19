# schwammerlRohr optional post-GNN curation and reassembly — prototype

This package does **not** modify your ML GNN, existing Snakemake rules, or containers. It builds on the earlier offline curation v2 viewer. It is deliberately a separable first implementation, not a validated end-to-end release.

## Stage 1: Prepare a complete contig evidence table

Use your original GNN `*__binning_results.tsv`. Pass additional TSVs keyed by `contig_id`. For numerical flagging, select the **actual columns** present in your files (prefixed as `filename.column` by the joiner):

```bash
python3 prepare_evidence.py \
  --assignments sample__binning_results.tsv \
  --evidence coverages.tsv --evidence composition.tsv \
  --gc-column composition.gc_content \
  --coverage-columns coverages.mean_coverage \
  --kmer-columns composition.AAAA,composition.AAAC,composition.AAAG \
  --output contig_evidence.tsv
```

`--kmer-columns` should cover the **full numerical k-mer vector** from your actual feature file, not merely the three illustrative features shown above. Automatic detection recognizes fields containing `tetra`, `kmer`, or `4mer`. Robust modified z-scores use the median and MAD of the **original** bin. Bins with too few contigs or zero MAD get no flag for that metric; this is *not* a biological confidence estimate. K-mer distance is relative to the bin's robust per-feature centroid, scaled dimension-wise; it is a heuristic that can be replaced later.

NCBI lineages can optionally be resolved against a **local** taxonomy dump:

```bash
python3 taxonomy_paths.py --input contig_evidence.tsv --output contig_evidence_taxonomy.tsv \
  --taxdump /path/to/ncbi_taxdump --taxid-column taxonomy.taxid
```

The taxonomy table's taxid column must exist and be identified explicitly. The lineage function does not call NCBI or require an API key. Download `nodes.dmp` and `names.dmp` from the official NCBI taxonomy dump separately. No key is hardcoded.

## Stage 2: Interactive review outside cluster

```bash
python3 build_curation.py --assignments sample__binning_results.tsv \
  --qc-table contig_evidence_taxonomy.tsv --output sample_curation.html
```

If taxonomy was not prepared, supply `contig_evidence.tsv`. Open the self-contained HTML locally. Filter flagged contigs; merge/split bins; tag possible haplotigs, alternative alleles, repeat/paralog copies, organellar variants, assembly artifacts and a related contig ID; annotate bin taxonomy. **Export both TSV and JSON**. The JSON is the complete snapshot of manual edits; the TSV is the input to the next stage. Taxid edits are not automatically resolved in the offline browser: resolve/export a consistent lineage later from the taxonomy snapshot or add an optional online service in a later release.

## Stage 3: Explicit second Snakemake run

The files `curation_export.py`, `recruit_reads.py` and `reassembly.smk` create curated FASTAs, recruited FASTQs, and optionally hifiasm GFAs, without rerunning the GNN. Keep this standalone Snakefile separate from the primary workflow for now:

```bash
snakemake -s reassembly.smk --cores 16 --config \
  curation_assembly=/shared/assembly.fasta \
  curation_assignments=/shared/curated_assignments.tsv \
  curation_reads=/shared/reads.fastq.gz \
  curation_bam=/shared/all_reads.sorted.bam \
  curation_script_dir=/shared/curation_scripts \
  curation_out=/shared/curation_run1 \
  curated_reassembly_inputs
```

The BAM must be the **full read-to-original-assembly alignment**, using exactly the same contig IDs; it is not automatically created or verified by this prototype. `samtools` must be on PATH. Recruitment currently filters unmapped, secondary, and supplementary alignments, uses alignment span and MAPQ, and flags competitive mappings between bins; it does **not** solve repeats, chimeras, subread/HiFi identifier differences, or all alignment ambiguity. Test against real data before trusting the reads as bin-specific. In the reassembly invocation, specify `--use-apptainer` and `--config hifiasm_image=docker://...` for an image you have verified, and request `.../hifiasm/bin_0001.bp.p_ctg.gfa` as target to run hifiasm. The runtime and image are **not validated here**.

### Important workflow limitations

- The reassembly module is a separate *second Snakemake invocation* for now; this is intentional so that GUI editing is never run on the HPC nodes.
- A changed curated TSV is a normal Snakemake input and triggers downstream reruns.
- The first exporter checks that every curated contig matches the FASTA. By default it also requires all FASTA contigs to be represented in the curation table.
- Per-bin FASTA file stems (`bin_0001`, etc.) are deterministic for each curation snapshot, with a manifest retaining curated bin IDs.
- The reassembler currently emits hifiasm GFA; graph-to-FASTA and post-reassembly plotting/alternative-contig interpretation are **not implemented yet**.
- The binning TSV and evidence joiner's actual field names need verification on a real sample; this package has synthetic-only tests.

# The original ML source lives at bin/ml/Cluster_ML and is never modified.
# GNN uses the existing Flye assembly and CoverM results. New feature rules
# only produce inputs not yet available from the main workflow.
GNN_INPUT_ROOT = config.get("gnn_input_root", "PacBio_output/gnn_inputs")
GNN_ML_DIR = "bin/ml/Cluster_ML"
GNN_PYTHON = config.get("gnn_python", "python3")
GNN_TIARA = config.get("gnn_use_tiara", True)

rule gnn_kmer_counts:
    input:
        fasta=f"{OUTPUT_DIR}/flye_output/{PROJECT_NR}" + "{sample}_flye/assembly.fasta"
    output:
        matrix=f"{OUTPUT_DIR}/gnn_features/" + "{sample}/assembly.k4.per_contig_matrix.csv"
    shell:
        "python3 bin/scripts/gnn_kmers.py --fasta {input.fasta:q} --output {output.matrix:q}"

rule gnn_bam_links:
    input:
        fasta=f"{OUTPUT_DIR}/flye_output/{PROJECT_NR}" + "{sample}_flye/assembly.fasta",
        coverage_dir=f"{OUTPUT_DIR}/coverm_output/{PROJECT_NR}" + "{sample}_mean_cov"
    output:
        links=f"{OUTPUT_DIR}/gnn_features/" + "{sample}/links.tsv"
    threads: min(THREADS, 8)
    shell:
        """
        python3 bin/scripts/gnn_links_from_cache.py \
          --coverage-dir {input.coverage_dir:q} --assembly {input.fasta:q} \
          --output {output.links:q} --threads {threads}
        """

rule gnn_prepare_inputs:
    input:
        fasta=f"{OUTPUT_DIR}/flye_output/{PROJECT_NR}" + "{sample}_flye/assembly.fasta",
        coverage_dir=f"{OUTPUT_DIR}/coverm_output/{PROJECT_NR}" + "{sample}_mean_cov",
        kmer=f"{OUTPUT_DIR}/gnn_features/" + "{sample}/assembly.k4.per_contig_matrix.csv",
        links=f"{OUTPUT_DIR}/gnn_features/" + "{sample}/links.tsv",
        tiara=lambda wc: [tiara_path(wc)] if GNN_TIARA and TIARA_MODE != "off" else [],
        tax=lambda wc: f"{OUTPUT_DIR}/diamond_output/{PROJECT_NR}{wc.sample}_taxonomy_prots/{PROJECT_NR}{wc.sample}_taxprot"
    output:
        prepared=directory(GNN_INPUT_ROOT + "/{sample}/input")
    params:
        tiara_arg=lambda wc: "--tiara " + tiara_path(wc) if GNN_TIARA and TIARA_MODE != "off" else "",
        tax_arg=lambda wc: f"--tax {OUTPUT_DIR}/diamond_output/{PROJECT_NR}{wc.sample}_taxonomy_prots/{PROJECT_NR}{wc.sample}_taxprot"
    shell:
        """
        python3 bin/scripts/gnn_stage_inputs.py \
          --assembly {input.fasta:q} \
          --coverage {input.coverage_dir:q}/mean_coverage.tsv \
          --kmers {input.kmer:q} --links {input.links:q} \
          {params.tiara_arg} {params.tax_arg} \
          --sample {wildcards.sample:q} --outdir {output.prepared:q}
        """

rule original_gnn:
    input:
        prepared=GNN_INPUT_ROOT + "/{sample}/input"
    output:
        done="PacBio_output/gnn/{sample}/.complete"
    params:
        ml=GNN_ML_DIR,
        python=GNN_PYTHON,
        input_dir=lambda wc: f"{GNN_INPUT_ROOT}/{wc.sample}/input",
        output_dir=lambda wc: f"PacBio_output/gnn/{wc.sample}"
    shell:
        """
        python3 bin/scripts/run_original_gnn.py \
          --ml-dir {params.ml:q} --sample {wildcards.sample:q} \
          --input-dir {params.input_dir:q} --output-dir {params.output_dir:q} \
          --python {params.python:q}
        """

"""Optional second Snakemake invocation after offline curation.

Invoke with --config curation_assembly=... curation_assignments=...
          curation_reads=... curation_bam=... curation_out=...
No GNN or upstream pipeline rules are modified.
"""
import csv
from pathlib import Path

CURATION_OUT=config.get('curation_out','curated_reassembly')
CURATION_ASSEMBLY=config['curation_assembly']
CURATION_ASSIGNMENTS=config['curation_assignments']
CURATION_READS=config['curation_reads']
CURATION_BAM=config['curation_bam']
CURATION_SCRIPT_DIR=config.get('curation_script_dir','bin/scripts')
HIFIASM_IMAGE=config.get('hifiasm_image')

# Read the curated snapshot to determine dynamic bin IDs at DAG construction.
with open(CURATION_ASSIGNMENTS) as f:
    _curated=list(csv.DictReader(f,delimiter='\t'))
BIN_COUNT=len({r['curated_bin'] for r in _curated})
STEMS=[f'bin_{i:04d}' for i in range(1,BIN_COUNT+1)]

rule curated_bin_export:
    input: assembly=CURATION_ASSEMBLY,assignments=CURATION_ASSIGNMENTS
    output: manifest=f'{CURATION_OUT}/bins/bins_manifest.tsv',mapping=f'{CURATION_OUT}/bins/contig_to_bin.tsv',fastas=expand(f'{CURATION_OUT}/bins/{{stem}}.fasta',stem=STEMS)
    shell:
        'python {CURATION_SCRIPT_DIR}/curation_export.py --assembly {input.assembly:q} --assignments {input.assignments:q} --output {CURATION_OUT:q}' + '/bins'

rule curated_read_recruitment:
    input: bam=CURATION_BAM,reads=CURATION_READS,mapping=rules.curated_bin_export.output.mapping
    output: manifest=f'{CURATION_OUT}/reads/recruitment_manifest.tsv',summary=f'{CURATION_OUT}/reads/read_recruitment.tsv',fastqs=expand(f'{CURATION_OUT}/reads/{{stem}}.fastq.gz',stem=STEMS),ambiguous=f'{CURATION_OUT}/reads/ambiguous.fastq.gz'
    shell:
        'python {CURATION_SCRIPT_DIR}/recruit_reads.py --bam {input.bam:q} --fastq {input.reads:q} --mapping {input.mapping:q} --output {CURATION_OUT:q}/reads'

# Hifiasm is optional: call explicit GFA targets; requires --config hifiasm_image=...
# Container paths are configured by the user; image provisioning belongs to the chosen runtime.
rule curated_hifiasm:
    input: reads=f'{CURATION_OUT}/reads/{{bin}}.fastq.gz',manifest=rules.curated_read_recruitment.output.manifest
    output: gfa=f'{CURATION_OUT}/hifiasm/{{bin}}.bp.p_ctg.gfa'
    threads: 16
    container: HIFIASM_IMAGE
    shell:
        'mkdir -p {CURATION_OUT:q}/hifiasm && hifiasm -t {threads} -l0 -o {CURATION_OUT:q}/hifiasm/{wildcards.bin} {input.reads:q}'

rule curated_reassembly_inputs:
    input:
        rules.curated_bin_export.output.manifest,
        rules.curated_read_recruitment.output.manifest
    output: touch(f'{CURATION_OUT}/inputs.ready')

# Run explicit hifiasm GFA targets per bin after checking recruitment manifest.

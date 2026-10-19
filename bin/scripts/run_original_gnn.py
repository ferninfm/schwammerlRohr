#!/usr/bin/env python3
"""Execute original Cluster_ML/main.py for a supplied sample without editing upstream."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

REQUIRED = ('assembly.fasta','mean_coverage.tsv','assembly.k4.per_contig_matrix.csv','links.tsv')
OPTIONAL = ('complete.classification.tsv',)

def run(args):
    ml_dir = Path(args.ml_dir).resolve()
    if not (ml_dir / 'main.py').is_file():
        raise SystemExit(f'Missing ML main.py: {ml_dir / "main.py"}')
    original = (ml_dir / 'main.py').read_text(encoding='utf-8')
    # Only substitute path configuration in the ephemeral entry point.
    pattern = r'(?m)^(\s*SAMPLE\s*=\s*DATA_DIR\s*/\s*)["\']pr_197_002["\'](\s*)$'
    matches = list(re.finditer(pattern, original))
    if len(matches) != 1:
        raise SystemExit('Expected one original SAMPLE = DATA_DIR / "pr_197_002" in main.py; upstream changed. No model execution.')
    changed = original[:matches[0].start()] + matches[0].group(1) + repr(args.sample) + matches[0].group(2) + original[matches[0].end():]
    src=Path(args.input_dir).resolve()
    missing=[name for name in REQUIRED if not (src/name).is_file()]
    if missing:
        raise SystemExit('Missing preprocessed ML inputs: ' + ', '.join(missing) + f' in {src}')
    dest=Path(args.output_dir).resolve()
    if dest.exists() and any(dest.iterdir()):
        raise SystemExit(f'Output directory must be empty or absent to prevent overwriting: {dest}')
    with tempfile.TemporaryDirectory(prefix='schwammerl-gnn-') as tmp:
        tmp=Path(tmp)
        input_dir=tmp/'cluster_data'/args.sample/'input'; input_dir.mkdir(parents=True)
        for name in REQUIRED+OPTIONAL:
            candidate=src/name
            if candidate.is_file(): (input_dir/name).symlink_to(candidate)
        # Optional preexisting classification inputs, using original filename conventions.
        for pat in (f'{args.sample}_taxprot',f'{args.sample}_diamond_blastx.tsv',
                    f'{args.sample}_classification_centrifuger.tsv',
                    f'{args.sample}_pf.kraken2.tsv',f'{args.sample}.kraken2.pluspfp.tsv'):
            candidate=src/pat
            if candidate.is_file(): (input_dir/pat).symlink_to(candidate)
        generated=tmp/'entrypoint.py'; generated.write_text(changed,encoding='utf-8')
        env=os.environ.copy()
        env['PYTHONPATH']=str(ml_dir)+os.pathsep+env.get('PYTHONPATH','')
        subprocess.run([args.python, str(generated)],cwd=tmp,env=env,check=True)
        result=tmp/'cluster_data'/args.sample/'output'
        if not result.is_dir() or not list(result.rglob('*binning_results.tsv')):
            raise SystemExit('GNN finished without its expected binning_results.tsv')
        dest.mkdir(parents=True,exist_ok=True)
        import shutil
        shutil.copytree(result,dest,dirs_exist_ok=True)
        (dest/'.complete').write_text('ok\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ml-dir',required=True,help='Directory with original main.py')
    p.add_argument('--sample',required=True)
    p.add_argument('--input-dir',required=True,help='Directory containing original ML inputs')
    p.add_argument('--output-dir',required=True)
    p.add_argument('--python',default=sys.executable)
    args=p.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',args.sample) or args.sample in ('.','..'):
        p.error('invalid sample ID')
    run(args)
if __name__=='__main__': main()

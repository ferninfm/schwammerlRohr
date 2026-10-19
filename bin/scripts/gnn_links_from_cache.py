#!/usr/bin/env python3
"""Use the pre-existing CoverM BAM cache for the original BAM links method."""
import argparse, subprocess, sys
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--coverage-dir',required=True);p.add_argument('--assembly',required=True);p.add_argument('--output',required=True);p.add_argument('--threads',default='4')
    a=p.parse_args();cache=Path(a.coverage_dir)/'bam_cache'
    bams=list(cache.rglob('*.bam')) if cache.exists() else []
    if len(bams)!=1:
        p.error(f'Expected exactly one retained BAM in {cache}, found {len(bams)}. Do not silently regenerate links or use a different BAM; check CoverM cache retention.')
    script=Path(__file__).resolve().parents[1]/'ml'/'bam_links.py'
    subprocess.run([sys.executable,str(script),'--bam',str(bams[0]),'--fasta',a.assembly,'--output',a.output,'--threads',a.threads],check=True)
if __name__=='__main__':main()

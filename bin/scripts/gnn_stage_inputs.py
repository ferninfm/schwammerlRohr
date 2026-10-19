#!/usr/bin/env python3
"""Stage ML-compatible source files without changing the upstream GNN code."""
import argparse, shutil
from pathlib import Path

def run(a):
    out=Path(a.outdir)
    if out.exists(): shutil.rmtree(out)
    out.mkdir(parents=True)
    items={'assembly.fasta':a.assembly,'mean_coverage.tsv':a.coverage,'assembly.k4.per_contig_matrix.csv':a.kmers,'links.tsv':a.links}
    if a.tiara: items['complete.classification.tsv']=a.tiara
    if a.tax: items[a.sample+'_taxprot']=a.tax
    for filename,src in items.items():
        src=Path(src).resolve()
        if not src.is_file(): raise FileNotFoundError(src)
        (out/filename).symlink_to(src)

def main():
    p=argparse.ArgumentParser()
    for x in ('assembly','coverage','kmers','links','sample','outdir'):p.add_argument('--'+x,required=True)
    p.add_argument('--tiara');p.add_argument('--tax');a=p.parse_args();run(a)
if __name__=='__main__':main()

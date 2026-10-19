#!/usr/bin/env python3
"""Generate the original 256-column Jellyfish-compatible noncanonical 4-mer counts."""
import argparse, csv
from collections import Counter
from itertools import product
from pathlib import Path

KMERS = [''.join(k) for k in product('ACGT', repeat=4)]

def sequences(path):
    name = None; chunks = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(chunks)
                name = line[1:].split()[0]; chunks=[]
            else: chunks.append(line.strip().upper())
    if name is not None: yield name, ''.join(chunks)

def run(src, dest):
    dest=Path(dest); dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['contig',*KMERS])
        for name,seq in sequences(src):
            counts=Counter(seq[i:i+4] for i in range(max(0,len(seq)-3)))
            w.writerow([name,*[counts[k] for k in KMERS]])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();run(a.fasta,a.output)

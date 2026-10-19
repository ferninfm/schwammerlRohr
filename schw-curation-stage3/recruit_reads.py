#!/usr/bin/env python3
"""Recruit original FASTQ reads using primary BAM alignments, competitively by aligned bases."""
import argparse,csv,gzip,subprocess
from collections import defaultdict
from pathlib import Path

def open_text(p,mode='rt'):
    return gzip.open(p,mode) if str(p).endswith('.gz') else open(p,mode)

def recruit(bam,fastq,mapping,output,min_mapq=10,min_aligned=500,ambiguity_ratio=.9):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    with open(mapping) as f:
        records=list(csv.DictReader(f,delimiter='\t'))
    contigs={r['contig_id']:r['bin_id'] for r in records}
    if len(contigs)!=len(records):raise ValueError('Duplicate contig ID in mapping')
    cmd=['samtools','view','-F','0x904',str(bam)] # unmapped, secondary, supplementary excluded
    strengths=defaultdict(lambda:defaultdict(int))
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True)
    assert proc.stdout
    try:
        for line in proc.stdout:
            f=line.rstrip('\n').split('\t')
            if len(f)<11:continue
            read,ref=f[0],f[2]
            if ref not in contigs or int(f[4])<min_mapq:continue
            # CIGAR query alignment measure: M/= /X, not clipping or insertions
            import re
            aligned=sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]) if op in 'M=X')
            if aligned>=min_aligned:strengths[read][contigs[ref]]+=aligned
    finally:
        proc.stdout.close()
    if proc.wait()!=0:raise RuntimeError('samtools view failed')
    selection={};ambiguous=set()
    for read,scores in strengths.items():
        top=sorted(scores.items(),key=lambda a:(-a[1],a[0]));best=top[0]
        if len(top)>1 and top[1][1]>=best[1]*ambiguity_ratio:ambiguous.add(read)
        else:selection[read]=best[0]
    names=sorted(set(contigs.values()));handles={b:gzip.open(out/(f'bin_{i:04d}.fastq.gz'),'wt') for i,b in enumerate(names,1)}
    handles['__ambiguous__']=gzip.open(out/'ambiguous.fastq.gz','wt')
    seen=set()
    try:
        with open_text(fastq) as inp:
            while True:
                h=inp.readline()
                if not h:break
                seq=inp.readline();plus=inp.readline();qual=inp.readline()
                if not (h.startswith('@') and plus.startswith('+') and seq and qual):raise ValueError('Invalid 4-line FASTQ')
                read=h[1:].split()[0]
                if read in seen:raise ValueError('Duplicate FASTQ read identifier '+read)
                seen.add(read)
                target='__ambiguous__' if read in ambiguous else selection.get(read)
                if target:handles[target].write(h+seq+plus+qual)
    finally:
        for h in handles.values():h.close()
    with open(out/'recruitment_manifest.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['bin_id','fastq'])
        for i,b in enumerate(names,1):w.writerow([b,str((out/f'bin_{i:04d}.fastq.gz').resolve())])
    with open(out/'read_recruitment.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['read_id','assigned_bin','ambiguous','best_aligned_bases','competing_bins'])
        for read,scores in strengths.items():
            top=sorted(scores.items(),key=lambda a:-a[1]);w.writerow([read,'' if read in ambiguous else selection[read],read in ambiguous,top[0][1],','.join(b for b,_ in top[1:])])
    print(f'Recruited {len(selection)} unambiguous and {len(ambiguous)} ambiguous reads; {len(seen)} original FASTQ records examined')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['bam','fastq','mapping','output']:p.add_argument('--'+k,required=True)
    p.add_argument('--min-mapq',type=int,default=10);p.add_argument('--min-aligned',type=int,default=500);p.add_argument('--ambiguity-ratio',type=float,default=.9)
    a=p.parse_args();recruit(a.bam,a.fastq,a.mapping,a.output,a.min_mapq,a.min_aligned,a.ambiguity_ratio)
if __name__=='__main__':main()

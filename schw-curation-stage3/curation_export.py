#!/usr/bin/env python3
"""Validate curated contig assignments and export per-bin FASTA without altering original assembly."""
import argparse,csv,gzip,json,re
from collections import defaultdict
from pathlib import Path

def read_fasta(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as src:
        name=None; parts=[]
        for line in src:
            if line.startswith('>'):
                if name is not None:yield name,''.join(parts)
                name=line[1:].split()[0];parts=[]
            elif name is not None:parts.append(line.strip())
        if name is not None:yield name,''.join(parts)

def export(assembly, assignments, output, allow_missing=False):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    with open(assignments,newline='') as f:
        rd=csv.DictReader(f,delimiter='\t')
        required={'contig_id','curated_bin'}
        if not required.issubset(rd.fieldnames or []):raise ValueError('Expected contig_id and curated_bin')
        rows=list(rd)
    ids=[r['contig_id'] for r in rows]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate contig in curated assignments')
    bin_for={r['contig_id']:r['curated_bin'] for r in rows}
    if any(not v for v in bin_for.values()):raise ValueError('Blank curated bin')
    # Filesystem safe, stable file names separate from original bin labels.
    names={b:f'bin_{i:04d}' for i,b in enumerate(sorted(set(bin_for.values())),1)}
    seen=set();counts=defaultdict(int);handles={}
    try:
        for cid,seq in read_fasta(assembly):
            if cid in seen:raise ValueError(f'Duplicate FASTA identifier {cid}')
            seen.add(cid)
            if cid not in bin_for:continue
            b=bin_for[cid]; counts[b]+=1
            if b not in handles:handles[b]=open(out/(names[b]+'.fasta'),'w')
            h=handles[b];h.write(f'>{cid}\n')
            for pos in range(0,len(seq),80):h.write(seq[pos:pos+80]+'\n')
    finally:
        for h in handles.values():h.close()
    missing=set(bin_for)-seen
    extra=seen-set(bin_for)
    if missing:raise ValueError(f'{len(missing)} curated contig(s) missing from FASTA; examples {sorted(missing)[:5]}')
    if extra and not allow_missing:raise ValueError(f'{len(extra)} FASTA contig(s) absent from curation; examples {sorted(extra)[:5]}')
    with open(out/'bins_manifest.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['bin_id','file_stem','fasta','contig_count'])
        for b in names:w.writerow([b,names[b],str((out/(names[b]+'.fasta')).resolve()),counts[b]])
    with open(out/'contig_to_bin.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['contig_id','bin_id'])
        for cid,b in bin_for.items():w.writerow([cid,b])
    print(f'Exported {len(names)} bins, {sum(counts.values())} contigs to {out}')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly',required=True);p.add_argument('--assignments',required=True);p.add_argument('--output',required=True)
    p.add_argument('--allow-unassigned-assembly-contigs',action='store_true')
    a=p.parse_args();export(a.assembly,a.assignments,a.output,a.allow_unassigned_assembly_contigs)
if __name__=='__main__':main()

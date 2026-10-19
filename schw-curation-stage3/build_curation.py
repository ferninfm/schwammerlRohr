#!/usr/bin/env python3
"""Build a self-contained, offline HTML curation viewer from the original GNN result TSV.
No model code changes and no third-party Python dependencies.
"""
import argparse
import csv
import gzip
import html
import json
from pathlib import Path
import re

KEY='contig_id'

def open_text(p):
    return gzip.open(p, 'rt', encoding='utf-8-sig', newline='') if str(p).endswith('.gz') else open(p, encoding='utf-8-sig', newline='')

def read_tsv(path, sectioned=False):
    with open_text(path) as f:
        header=None; records=[]
        for line in f:
            line=line.rstrip('\n\r')
            if not line.strip():continue
            if sectioned and (line.startswith('===') or line.startswith('---') or (not '\t' in line and not line.startswith(KEY+'\t'))):continue
            vals=next(csv.reader([line],delimiter='\t'))
            if header is None:
                if KEY not in vals:raise ValueError(f'{path}: no contig_id header')
                header=vals;continue
            if len(vals)!=len(header):
                if sectioned:continue
                raise ValueError(f'{path}: row with {len(vals)} columns, expected {len(header)}')
            records.append(dict(zip(header,vals)))
    if header is None:raise ValueError(f'{path}: no header found')
    return header,records

def make_data(assignment, evidence):
    headers,rows=read_tsv(assignment,True)
    if not {'contig_id','bin_id'}.issubset(headers):raise ValueError('Expected contig_id and bin_id in ML output')
    seen=set()
    for r in rows:
        cid=r[KEY]
        if cid in seen:raise ValueError(f'Duplicate contig ID in assignments: {cid}')
        seen.add(cid)
        r['original_bin']=r['bin_label'] if r.get('bin_label') else ('Noise' if r['bin_id']=='-1' else 'Organelle' if r['bin_id']=='-2' else 'Bin '+r['bin_id'])
        r['curated_bin']=r['original_bin']
    columns=list(headers)
    for p in evidence:
        h,rr=read_tsv(p)
        if len({r[KEY] for r in rr})!=len(rr):raise ValueError(f'{p}: duplicate contig_id (aggregate gene hits first)')
        idx={r[KEY]:r for r in rr}
        prefix=re.sub('[^A-Za-z0-9_]+','_',Path(p).name.replace('.tsv','').replace('.gz',''))
        for c in h:
            if c==KEY:continue
            name=f'{prefix}.{c}'
            columns.append(name)
            for r in rows:r[name]=idx.get(r[KEY],{}).get(c,'')
    return rows,columns

def build(assignment, output, evidence):
    rows,cols=make_data(assignment,evidence)
    template=Path(__file__).with_name('viewer.html').read_text(encoding='utf-8')
    safe_json=json.dumps({'rows':rows,'columns':cols,'source':Path(assignment).name},ensure_ascii=False,separators=(',',':')).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    result=template.replace('/*__DATA__*/null',safe_json)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(result,encoding='utf-8')
    print(f'Wrote {output} ({len(rows)} contigs, {len(cols)} source fields, {output.stat().st_size/1048576:.1f} MB)')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assignments',required=True,type=Path,help='Original *__binning_results.tsv from ML')
    p.add_argument('--evidence',action='append',default=[],type=Path,help='Optional contig_id-keyed TSV; repeat')
    p.add_argument('--output',required=True,type=Path,help='Standalone HTML')
    p.add_argument('--qc-table',type=Path,help='Prepared full contig table, already containing qc_* columns')
    a=p.parse_args();build(a.assignments,a.output,([a.qc_table] if a.qc_table else []) + a.evidence)
if __name__=='__main__':main()

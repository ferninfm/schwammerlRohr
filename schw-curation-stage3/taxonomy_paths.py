#!/usr/bin/env python3
"""Resolve taxids to local NCBI taxdump lineages; never needs or embeds an API key.
Requires nodes.dmp and names.dmp downloaded separately from NCBI.
"""
import argparse,csv
from pathlib import Path

def parse_dump(path):
    with open(path) as f:
        for line in f:yield [field.strip() for field in line.split('|')]

def load_taxdump(folder):
    parents={};ranks={};names={}
    for fields in parse_dump(Path(folder)/'nodes.dmp'):
        parents[fields[0]]=fields[1];ranks[fields[0]]=fields[2]
    for fields in parse_dump(Path(folder)/'names.dmp'):
        if fields[3]=='scientific name':names[fields[0]]=fields[1]
    return parents,ranks,names

def resolve(taxid,parents,ranks,names):
    seen=set();path=[];taxid=str(taxid).strip();cur=taxid
    while cur in parents and cur not in seen:
        seen.add(cur);path.append((cur,names.get(cur,''),ranks.get(cur,'')))
        nxt=parents[cur]
        if nxt==cur:break
        cur=nxt
    if not path:return '', '', 'unresolved'
    if cur not in parents:return '', '', 'broken_lineage'
    return names.get(taxid,''),'; '.join(f'{n} ({i})' for i,n,_ in reversed(path)),'resolved'

def annotate(inp,out,taxdump,col):
    parents,ranks,names=load_taxdump(taxdump)
    with open(inp,newline='') as f:reader=csv.DictReader(f,delimiter='\t');rows=list(reader);fields=reader.fieldnames or []
    if col not in fields:raise ValueError('Missing taxid column: '+col)
    for r in rows:
        name,path,status=resolve(r.get(col,''),parents,ranks,names)
        r['ncbi_taxon_name']=name;r['ncbi_taxon_path']=path;r['ncbi_resolution_status']=status
    with open(out,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields+['ncbi_taxon_name','ncbi_taxon_path','ncbi_resolution_status'],delimiter='\t');w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--taxdump',required=True);p.add_argument('--taxid-column',default='taxid')
    a=p.parse_args();annotate(a.input,a.output,a.taxdump,a.taxid_column)
if __name__=='__main__':main()

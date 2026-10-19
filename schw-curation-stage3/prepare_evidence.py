#!/usr/bin/env python3
"""Build contig-level metrics and robust within-bin outlier flags (stdlib only)."""
import argparse,csv,math,statistics
from pathlib import Path
from build_curation import make_data

def number(value):
    try:
        x=float(value)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None

def median_mad(v):
    med=statistics.median(v)
    return med,statistics.median(abs(x-med) for x in v)

def metric_columns(rows,choice,patterns):
    if choice:return [x.strip() for x in choice.split(',') if x.strip()]
    return [k for k in rows[0] if any(p in k.lower() for p in patterns) and sum(number(r.get(k)) is not None for r in rows)>=max(2,len(rows)//3)] if rows else []

def calculate(rows, gc=None, coverage=None, kmers=None, threshold=3.5, min_group=5):
    if not rows:return rows
    gc_cols=metric_columns(rows,gc,['gc_content','gc_percent','gc_pct'])
    cov_cols=metric_columns(rows,coverage,['mean_coverage','coverage','depth'])
    k_cols=metric_columns(rows,kmers,['tetra','kmer','4mer'])
    if not (gc_cols or cov_cols or k_cols):raise ValueError('No metrics found: specify --gc-column, --coverage-columns or --kmer-columns, with --evidence tables')
    bybin={}
    for r in rows:bybin.setdefault(r['original_bin'],[]).append(r)
    for group in bybin.values():
        for r in group:
            flags=[]; detail=[]
            for kind,cols in [('gc',gc_cols),('coverage',cov_cols)]:
                for col in cols:
                    v=[number(x.get(col)) for x in group];v=[x for x in v if x is not None]
                    x=number(r.get(col))
                    if x is None or len(v)<min_group:continue
                    center,mad=median_mad(v)
                    # zero MAD means no reliable spread estimate: do not overflag
                    if mad<=0:continue
                    score=.67448975*(x-center)/mad
                    r[f'qc_{kind}_z_{col}']=f'{score:.3f}'
                    if abs(score)>threshold: flags.append(kind);detail.append(f'{col}:{score:.2f}')
            if k_cols:
                # Euclidean distance in robust standardized composition space; no inferred k-mer normalization
                norms={}
                for col in k_cols:
                    vals=[number(x.get(col)) for x in group]; vals=[v for v in vals if v is not None]
                    if len(vals)>=min_group:
                        med,mad=median_mad(vals)
                        if mad>0:norms[col]=(med,mad)
                if len(norms)>=2 and len(group)>=min_group:
                    distances=[]
                    for candidate in group:
                        components=[((number(candidate.get(c))-med)/(1.4826*mad))**2 for c,(med,mad) in norms.items() if number(candidate.get(c)) is not None]
                        distances.append(math.sqrt(sum(components)/len(components)) if len(components)>=2 else None)
                    finite=[v for v in distances if v is not None]
                    if len(finite)>=min_group:
                        center,mad=median_mad(finite)
                        idx=group.index(r); dist=distances[idx]
                        if dist is not None and mad>0:
                            z=.67448975*(dist-center)/mad
                            r['qc_kmer_distance']=f'{dist:.4f}';r['qc_kmer_z']=f'{z:.3f}'
                            if z>threshold:flags.append('kmer');detail.append(f'kmer:{z:.2f}')
            r['qc_outlier_flags']=','.join(sorted(set(flags)))
            r['qc_outlier_detail']=';'.join(detail)
            r['qc_bin_size']=str(len(group))
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assignments',required=True);p.add_argument('--evidence',action='append',default=[])
    p.add_argument('--gc-column');p.add_argument('--coverage-columns');p.add_argument('--kmer-columns')
    p.add_argument('--threshold',type=float,default=3.5);p.add_argument('--min-bin-size',type=int,default=5)
    p.add_argument('--output',required=True)
    a=p.parse_args();rows,_=make_data(a.assignments,a.evidence)
    calculate(rows,a.gc_column,a.coverage_columns,a.kmer_columns,a.threshold,a.min_bin_size)
    columns=list(dict.fromkeys(k for r in rows for k in r))
    Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    with open(a.output,'w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
    print(f'{len(rows)} contigs; {sum(bool(r["qc_outlier_flags"]) for r in rows)} flagged')
if __name__=='__main__':main()

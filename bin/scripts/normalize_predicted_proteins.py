#!/usr/bin/env python3
"""Normalize MetaEuk/Prodigal protein IDs; retain lossless original-ID map.

Uses MetaEuk C_acc (pipe header field 2) and Prodigal GFF ID-to-seqid
mapping. Strictly fail on unrecognized headers instead of misassigning genes.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from route_tiara import iter_fasta

PRODIGAL_ID = re.compile(r'(?:^|[;\s])ID=([^;\s]+)')


def fasta_ids(path):
    ids = set()
    for header, _ in iter_fasta(path):
        seq_id = header.split()[0]
        if seq_id in ids:
            raise ValueError(f'Duplicate contig ID in assembly: {seq_id}')
        ids.add(seq_id)
    return ids


def prodigal_gff_map(gff):
    mapping = {}
    with open(gff, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) < 9 or row[2] != 'CDS':
                continue
            match = PRODIGAL_ID.search(row[8])
            if match:
                gene_id = match.group(1)
                if gene_id in mapping:
                    raise ValueError(f'Duplicate Prodigal GFF ID {gene_id}')
                mapping[gene_id] = row[0]
    return mapping


def normalize(assembly, metaeuk, prodigal, gff, proteins_out, mapping_out):
    contigs = fasta_ids(assembly)
    prodigal_map = prodigal_gff_map(gff)
    Path(proteins_out).parent.mkdir(parents=True, exist_ok=True)
    Path(mapping_out).parent.mkdir(parents=True, exist_ok=True)
    records = 0
    with open(proteins_out, 'w') as proteins, open(mapping_out, 'w', newline='') as mapping:
        writer = csv.writer(mapping, delimiter='\t', lineterminator='\n')
        writer.writerow(['protein_id', 'contig_id', 'predictor', 'original_protein_header'])
        for source, path in [('metaeuk', metaeuk), ('prodigal', prodigal)]:
            for original_header, sequence in iter_fasta(path):
                if source == 'metaeuk':
                    header_fields = original_header.split()[0].split('|')
                    if len(header_fields) < 8:
                        raise ValueError(f'MetaEuk protein header missing expected fields: {original_header}')
                    contig = header_fields[1]
                else:
                    id_match = PRODIGAL_ID.search(original_header)
                    if not id_match:
                        raise ValueError(f'Prodigal protein has no ID= tag: {original_header}')
                    raw_gene_id = id_match.group(1)
                    if raw_gene_id not in prodigal_map:
                        raise ValueError(f'Prodigal ID {raw_gene_id!r} missing from GFF')
                    contig = prodigal_map[raw_gene_id]
                if contig not in contigs:
                    raise ValueError(f'Predicted protein mapped to unknown contig: {contig!r}')
                if not sequence:
                    raise ValueError(f'Empty predicted protein: {original_header}')
                records += 1
                identifier = f"{'ME' if source == 'metaeuk' else 'PG'}_{records:08d}"
                proteins.write('>' + identifier + '\n')
                for i in range(0, len(sequence), 80):
                    proteins.write(sequence[i:i + 80] + '\n')
                writer.writerow([identifier, contig, source, original_header])
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly', required=True)
    p.add_argument('--metaeuk', required=True)
    p.add_argument('--prodigal', required=True)
    p.add_argument('--prodigal-gff', required=True)
    p.add_argument('--proteins', required=True)
    p.add_argument('--mapping', required=True)
    a = p.parse_args()
    print(f"Normalized {normalize(a.assembly, a.metaeuk, a.prodigal, a.prodigal_gff, a.proteins, a.mapping)} proteins")


if __name__ == '__main__':
    main()

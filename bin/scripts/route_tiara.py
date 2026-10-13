#!/usr/bin/env python3
"""Route assembly contigs to gene predictors using Tiara's two-stage TSV.

Conservative routing: unknown and organellar contigs enter BOTH branches.
The TSV is complete even for short contigs missing from Tiara classification.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

EUK = {'eukarya', 'unknown', 'organelle', 'mitochondria', 'plastid'}
PROK = {'bacteria', 'archaea', 'prokarya', 'unknown', 'organelle', 'mitochondria', 'plastid'}
LABELS = EUK | PROK


def iter_fasta(path):
    header = None
    chunks = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header = line[1:].rstrip('\n\r')
                chunks = []
            elif line.strip():
                if header is None:
                    raise ValueError(f'FASTA sequence before header in {path}')
                chunks.append(line.strip())
    if header is not None:
        yield header, ''.join(chunks)


def read_tiara(path):
    labels = {}
    with open(path, encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for line_num, row in enumerate(reader, 1):
            if not row or not row[0].strip() or row[0].startswith('#'):
                continue
            raw_id = row[0].strip()
            if raw_id.lower() in {'sequence id', 'sequence_id', 'seqid', 'contig_id', 'sequence_id'}:
                continue
            contig_id = raw_id.split()[0]
            if len(row) < 2:
                raise ValueError(f'Tiara {path}:{line_num}: expected at least two columns')
            label = row[1].strip().lower()
            if label not in LABELS:
                raise ValueError(f'Tiara {path}:{line_num}: unexpected classification {label!r}')
            if contig_id in labels:
                raise ValueError(f'Duplicate Tiara contig ID {contig_id!r}')
            stage2 = row[2].strip() if len(row) > 2 else ''
            labels[contig_id] = (label, stage2)
    return labels


def write_fasta_record(f, header, sequence):
    f.write('>' + header + '\n')
    for i in range(0, len(sequence), 80):
        f.write(sequence[i:i + 80] + '\n')


def route(assembly, classification, euk_out, prok_out, route_out):
    labels = read_tiara(classification) if classification else {}
    seen = set()
    counts = {'euk': 0, 'prok': 0, 'both': 0, 'missing_tiara': 0}
    for path in (euk_out, prok_out, route_out):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(euk_out, 'w') as euk, open(prok_out, 'w') as prok, open(route_out, 'w', newline='') as report:
        writer = csv.writer(report, delimiter='\t', lineterminator='\n')
        writer.writerow(['contig_id', 'tiara_first_stage', 'tiara_second_stage', 'metaeuk', 'prodigal', 'length_bp'])
        for header, seq in iter_fasta(assembly):
            cid = header.split()[0]
            if cid in seen:
                raise ValueError(f'Duplicate assembly contig identifier {cid!r}')
            seen.add(cid)
            label, stage2 = labels.get(cid, ('unknown', ''))
            if cid not in labels:
                counts['missing_tiara'] += 1
            is_euk, is_prok = label in EUK, label in PROK
            if is_euk:
                write_fasta_record(euk, header, seq)
                counts['euk'] += 1
            if is_prok:
                write_fasta_record(prok, header, seq)
                counts['prok'] += 1
            if is_prok and is_euk:
                counts['both'] += 1
            if not (is_euk or is_prok):
                raise AssertionError(f'Lost contig {cid!r} ({label})')
            writer.writerow([cid, label, stage2, int(is_euk), int(is_prok), len(seq)])
    if not seen:
        raise ValueError('Assembly FASTA contains no contigs')
    extras = labels.keys() - seen
    if extras:
        raise ValueError(f'Tiara contains {len(extras)} IDs absent from assembly; e.g. {sorted(extras)[0]}')
    return counts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly', required=True)
    p.add_argument('--tiara', default=None,
                   help='Tiara TSV (omit to send every contig to both predictors)')
    p.add_argument('--euk', required=True)
    p.add_argument('--prok', required=True)
    p.add_argument('--report', required=True)
    args = p.parse_args()
    print(json.dumps(route(args.assembly, args.tiara, args.euk, args.prok, args.report)))


if __name__ == '__main__':
    main()

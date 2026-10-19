#!/usr/bin/env python3
"""
bam_links.py

Derives contig-contig links from supplementary BAM alignments,
excluding telomeric ends.

Equivalent to incorporate_bam_bubbles() from report_kiruna_v2.Rmd.
The taxonomy-reassignment part of that function is intentionally omitted —
it will be handled inside the GNN model.

Contig lengths and telomere detection are derived directly from the FASTA,
so no separate lengths file or telomere TSV is needed.

Telomere detection: scans the first/last --scan-window bp of each contig for
TTAGGG or CCCTAA repeats — the same two motifs hardcoded in incorporate_bam_bubbles().

Usage:
    python bam_links.py \\
        --bam    <file.bam> \\
        --fasta  <assembly.fasta> \\
        --output <links.tsv> \\
        [--flank 10000] \\
        [--threshold 0.7] \\
        [--scan-window 500] \\
        [--threads 4]
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pysam

# ── Constants ────────────────────────────────────────────────────────────────

# Motifs checked in incorporate_bam_bubbles() — hardcoded in the R code
TELOMERE_MOTIFS = {"TTAGGG", "CCCTAA"}

FLANK_DEFAULT       = 10_000
THRESHOLD_DEFAULT   = 0.5
SCAN_WINDOW_DEFAULT = 500


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bam",         required=True,
                   help="BAM file (sorted, with supplementary alignments)")
    p.add_argument("--fasta",       required=True,
                   help="Assembly FASTA (assembly.fasta from Flye)")
    p.add_argument("--output",      required=True,
                   help="Output TSV: contig_a  contig_b  jaccard")
    p.add_argument("--flank",       type=int,   default=FLANK_DEFAULT,
                   help=f"Flank size in bp (default: {FLANK_DEFAULT})")
    p.add_argument("--threshold",   type=float, default=THRESHOLD_DEFAULT,
                   help=f"Min Jaccard to keep edge (default: {THRESHOLD_DEFAULT})")
    p.add_argument("--scan-window", type=int,   default=SCAN_WINDOW_DEFAULT,
                   help=f"Bases scanned at each contig end for telomere motifs "
                        f"(default: {SCAN_WINDOW_DEFAULT})")
    p.add_argument("--threads",     type=int,   default=4,
                   help="Threads for BAM reading (default: 4)")
    return p.parse_args()


# ── FASTA: lengths + telomere detection ──────────────────────────────────────

def read_fasta(
    fasta_path:  str,
    scan_window: int,
) -> tuple[dict[str, int], set, set]:
    """
    Parse the assembly FASTA in one pass.

    Returns:
        lengths   — {contig_id: length}
        left_tel  — contigs with TTAGGG/CCCTAA in the first scan_window bp
        right_tel — contigs with TTAGGG/CCCTAA in the last  scan_window bp

    Telomere detection mirrors the two motifs hardcoded in
    incorporate_bam_bubbles():
        telomeres$leftmotif  %in% c("TTAGGG", "CCCTAA")
        telomeres$rightmotif %in% c("TTAGGG", "CCCTAA")
    """
    lengths:   dict[str, int] = {}
    left_tel:  set            = set()
    right_tel: set            = set()

    current_name: str       = ""
    current_seq:  list[str] = []

    def _process(name: str, seq_parts: list[str]) -> None:
        seq = "".join(seq_parts).upper()
        lengths[name] = len(seq)

        left_end  = seq[:scan_window]
        right_end = seq[-scan_window:]

        if any(m in left_end  for m in TELOMERE_MOTIFS):
            left_tel.add(name)
        if any(m in right_end for m in TELOMERE_MOTIFS):
            right_tel.add(name)

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    _process(current_name, current_seq)
                # contig name is everything up to the first space
                current_name = line[1:].split()[0]
                current_seq  = []
            else:
                current_seq.append(line)

    if current_name:
        _process(current_name, current_seq)

    return lengths, left_tel, right_tel


# ── Include-range logic (mirrors R exactly) ───────────────────────────────────

def build_include_ranges(
    lengths:   dict[str, int],
    left_tel:  set,
    right_tel: set,
    flank:     int,
) -> dict[str, tuple[int, int, int, int]]:
    """
    For each contig build (sl, el, sr, er) — the two flank intervals to use.
    Coordinates are 0-based, half-open: [sl, el) and [sr, er).
    A deactivated flank is stored as equal start/end (sl==el or sr==er).

    Mirrors the R logic in incorporate_bam_bubbles():
      sl=1, el=10000, sr=len-10000, er=len   (R 1-based)
    → sl=0, el=10000, sr=len-10000, er=len   (Python 0-based)
    """
    ranges: dict[str, tuple[int, int, int, int]] = {}

    for contig, length in lengths.items():

        # Initial flanks
        sl = 0
        el = min(flank, length)
        sr = max(0, length - flank)
        er = length

        # R: include_ranges$sr[include_ranges$sr < 1] <- 1
        if sr < 1:
            sr = 1

        # Telomere exclusion
        if contig in left_tel:
            sl, el = 0, 0
        if contig in right_tel:
            sr, er = 0, 0

        # R: Take out duplicated (sr==1 → right flank overlaps left start)
        if sr == 1:
            sr, er = 0, 0

        # Collapse overlapping flanks into single range
        # R: ir <- el > sr & sr != 0 → el=er, sr=0, er=0
        if sr != 0 and el > sr:
            el = er
            sr, er = 0, 0

        ranges[contig] = (sl, el, sr, er)

    return ranges


def _pos_in_range(pos: int, sl: int, el: int, sr: int, er: int) -> bool:
    """True if 0-based start position falls in either active flank interval."""
    in_left  = el > sl and sl <= pos < el
    in_right = er > sr and sr <= pos < er
    return in_left or in_right


def _assign_flank(pos: int, length: int, contig: str) -> str:
    """Mirrors R: pos >= tlen/2 → _r, else → _l."""
    if pos >= length / 2:
        return f"{contig}_r"
    return f"{contig}_l"


# ── BAM reading ───────────────────────────────────────────────────────────────

def collect_flank_reads(
    bam_path: str,
    ranges:   dict[str, tuple[int, int, int, int]],
    lengths:  dict[str, int],
    threads:  int,
) -> dict[str, set]:
    """
    Return {flank_id: set_of_read_qnames}.

    Mirrors R logic:
      1. Find qnames of reads with flag 2048 or 2064 (supplementary).
      2. Keep ALL alignments (primary + supplementary) of those reads.
      3. Filter to positions within include_ranges.
      4. Label each read as contig_l or contig_r.
    """
    # Pass 1: supplementary qnames — use pysam native flag instead of exact values
    print("[bam_links] Pass 1 — collecting supplementary read names ...", flush=True)
    supp_qnames: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_supplementary:
                supp_qnames.add(read.query_name)

    print(f"[bam_links]   {len(supp_qnames):,} reads with supplementary alignments",
          flush=True)

    # Pass 2: assign to flanks
    print("[bam_links] Pass 2 — assigning reads to flanks ...", flush=True)
    flank_reads: dict[str, set] = defaultdict(set)
    skipped_range     = 0
    skipped_no_contig = 0

    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
        for read in bam.fetch(until_eof=True):
            if read.query_name not in supp_qnames:
                continue

            contig = read.reference_name
            if contig is None or contig not in ranges:
                skipped_no_contig += 1
                continue

            pos = read.reference_start  # 0-based
            sl, el, sr, er = ranges[contig]

            if not _pos_in_range(pos, sl, el, sr, er):
                skipped_range += 1
                continue

            flank_id = _assign_flank(pos, lengths[contig], contig)
            flank_reads[flank_id].add(read.query_name)

    total = sum(len(v) for v in flank_reads.values())
    print(
        f"[bam_links]   {total:,} read-flank assignments across {len(flank_reads)} flanks",
        flush=True,
    )
    if skipped_range or skipped_no_contig:
        print(
            f"[bam_links]   skipped: {skipped_range:,} out-of-range, "
            f"{skipped_no_contig:,} contig not in FASTA",
            flush=True,
        )

    return dict(flank_reads)


# ── Jaccard computation ───────────────────────────────────────────────────────

def compute_jaccard_edges(
    flank_reads: dict[str, set],
    threshold:   float,
) -> list[tuple[str, str, float]]:
    """
    Compute Jaccard for every flank pair.

    Mirrors R:
        sum(foo_net[,fa]!=0 & foo_net[,fb]!=0) /
        sum(foo_net[,fa]!=0 | foo_net[,fb]!=0)
    which is the Jaccard index on the read sets.
    """
    flanks  = list(flank_reads.keys())
    n       = len(flanks)
    n_pairs = n * (n - 1) // 2

    print(
        f"[bam_links] Computing Jaccard for {n} flanks ({n_pairs:,} pairs) ...",
        flush=True,
    )

    edges: list[tuple[str, str, float]] = []

    for i, fa in enumerate(flanks):
        sa = flank_reads[fa]
        for fb in flanks[i + 1:]:
            sb = flank_reads[fb]

            intersection = len(sa & sb)
            if intersection == 0:
                continue

            jaccard = intersection / len(sa | sb)
            if jaccard >= threshold:
                edges.append((fa, fb, jaccard))

    print(
        f"[bam_links]   {len(edges):,} flank-flank edges with Jaccard >= {threshold}",
        flush=True,
    )
    return edges


# ── Flank → Contig edges ──────────────────────────────────────────────────────

def flanks_to_contig_edges(
    flank_edges: list[tuple[str, str, float]],
) -> list[tuple[str, str, float]]:
    """
    Convert flank-level edges to contig-level edges.
    Same-contig pairs (l↔r) are dropped.
    If multiple flank pairs link the same contig pair, keep max Jaccard.
    """
    best: dict[tuple[str, str], float] = {}

    for fa, fb, j in flank_edges:
        ca = fa.rsplit("_", 1)[0]
        cb = fb.rsplit("_", 1)[0]

        if ca == cb:
            continue

        key = (min(ca, cb), max(ca, cb))
        if j > best.get(key, 0.0):
            best[key] = j

    result = [(a, b, j) for (a, b), j in sorted(best.items())]
    print(f"[bam_links]   {len(result):,} contig-contig edges", flush=True)
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    # --- FASTA: lengths + telomere detection ---
    print("[bam_links] Reading FASTA ...", flush=True)
    lengths, left_tel, right_tel = read_fasta(args.fasta, args.scan_window)
    print(
        f"[bam_links]   {len(lengths):,} contigs | "
        f"{len(left_tel)} left-telomere ends | "
        f"{len(right_tel)} right-telomere ends",
        flush=True,
    )

    # --- Include ranges ---
    ranges = build_include_ranges(lengths, left_tel, right_tel, args.flank)

    fully_excluded = sum(1 for sl, el, sr, er in ranges.values() if sl == el == sr == er == 0)
    if fully_excluded:
        print(
            f"[bam_links]   {fully_excluded} contigs fully excluded "
            f"(both ends telomeric or too short)",
            flush=True,
        )

    # --- BAM ---
    flank_reads = collect_flank_reads(args.bam, ranges, lengths, args.threads)

    if not flank_reads:
        print(
            "[bam_links] ERROR: no flank reads collected. "
            "Check that the BAM contains supplementary alignments and that "
            "contig names match between BAM and FASTA.",
            file=sys.stderr,
        )
        sys.exit(1)

    # --- Jaccard → contig edges ---
    flank_edges  = compute_jaccard_edges(flank_reads, args.threshold)
    contig_edges = flanks_to_contig_edges(flank_edges)

    # --- Write output ---
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w") as fh:
        fh.write("contig_a\tcontig_b\tjaccard\n")
        for a, b, j in contig_edges:
            fh.write(f"{a}\t{b}\t{j:.6f}\n")

    print(f"[bam_links] Done → {out}", flush=True)


if __name__ == "__main__":
    main()

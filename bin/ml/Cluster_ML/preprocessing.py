from pathlib import Path
from typing import Dict, Union

import numpy as np
import pandas as pd


# ── GC-Content + Länge ────────────────────────────────────────────────────────

def gc_content_length_table(fasta_path) -> pd.DataFrame:
    """Liest FASTA → DataFrame(contig, gc_content, length)."""
    results = []
    contig  = None
    a = c = g = t = length = 0

    def flush():
        nonlocal contig, a, c, g, t, length
        if contig is None:
            return
        n = a + c + g + t
        results.append((contig, (g + c) / n if n else 0.0, length))
        a = c = g = t = length = 0

    with Path(fasta_path).open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                contig = line[1:].split()[0]
            else:
                seq = line.upper()
                length += len(seq)
                a += seq.count("A")
                c += seq.count("C")
                g += seq.count("G")
                t += seq.count("T")
    flush()
    return pd.DataFrame(results, columns=["contig", "gc_content", "length"])


# ── Coverage ──────────────────────────────────────────────────────────────────

def coverage_table(path) -> pd.DataFrame:
    """Liest CoverM TSV → DataFrame(contig_id, coverage)."""
    df = pd.read_csv(path, sep="\t").iloc[:, :2]
    df.columns = ["contig_id", "coverage"]
    return df


# ── k-mer Verarbeitung ────────────────────────────────────────────────────────

def _reverse_complement(kmer: str) -> str:
    return kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def reduce_kmer_to_canonical(df: pd.DataFrame) -> pd.DataFrame:
    """Reduziert k-mer-Spalten auf kanonische Form (min(kmer, revcomp))."""
    canon = {k: min(k, _reverse_complement(k)) for k in df.columns}
    return df.rename(columns=canon).T.groupby(level=0).sum().T


def normalize_kmer_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Normiert absolute k-mer-Counts auf relative Häufigkeiten (Zeilensumme = 1)."""
    return df.div(df.sum(axis=1), axis=0)


def process_kmer_matrix(path: Union[str, Path], canonical: bool = False) -> pd.DataFrame:
    """Liest k-mer-Matrix CSV, optional kanonisch reduziert, gibt normierte Frequenzen zurück."""
    df = pd.read_csv(Path(path), index_col="contig")
    df.index.name = "contig_id"
    if canonical:
        df = reduce_kmer_to_canonical(df)
    return normalize_kmer_counts(df)


# ── Tiara ─────────────────────────────────────────────────────────────────────

def _tiara_domain_group(cls: str) -> str:
    """Mappt Tiara-Klasse → Domain-Label."""
    cls = str(cls).lower().strip()
    if cls == "eukarya":
        return "eukarya"
    if cls in ("bacteria", "archaea", "prokarya"):
        return "prokaryote"
    if cls == "organelle_mito":
        return "organelle_mito"
    if cls == "organelle_plastid":
        return "organelle_plastid"
    return "unknown"


def load_tiara(tiara_path: Union[str, Path]) -> pd.DataFrame:
    """Liest Tiara TSV → DataFrame(contig_id, tiara_class).

    Löst Organellen-Subtypen (mitochondrion/plastid) über zweite Tiara-Stufe auf.
    """
    df = pd.read_csv(Path(tiara_path), sep="\t")
    id_col  = df.columns[0]
    fst_col = next(
        (c for c in ("class_fst_tier", "class_fst_stage", "class") if c in df.columns),
        df.columns[1],
    )
    snd_col = next(
        (c for c in ("class_snd_tier", "class_snd_stage") if c in df.columns),
        None,
    )

    def _classify_row(row) -> str:
        fst_raw = str(row[fst_col]).lower().strip()
        domain  = _tiara_domain_group(fst_raw)
        if domain == "unknown" and fst_raw == "organelle" and snd_col is not None:
            snd = str(row[snd_col]).lower().strip()
            if snd in ("mitochondrion", "mitochondria"):
                return "organelle_mito"
            if snd == "plastid":
                return "organelle_plastid"
            return "organelle_mito"
        return domain

    result = pd.DataFrame({
        "contig_id":   df[id_col].astype(str),
        "tiara_class": df.apply(_classify_row, axis=1),
    }).set_index("contig_id")
    print(f"[tiara] {len(result)} Contigs: {result['tiara_class'].value_counts().to_dict()}")
    return result


def tiara_node_features(contig_ids: pd.Index, tiara_df: pd.DataFrame) -> np.ndarray:
    """Gibt (n, 4) binäres Array zurück: [is_eukarya, is_prokaryote, is_organelle_mito, is_organelle_plastid]."""
    group_to_col: Dict[str, int] = {
        "eukarya": 0, "prokaryote": 1, "organelle_mito": 2, "organelle_plastid": 3,
    }
    feat = np.zeros((len(contig_ids), 4), dtype=np.float32)
    for i, cid in enumerate(contig_ids):
        if cid in tiara_df.index:
            col = group_to_col.get(tiara_df.loc[cid, "tiara_class"])
            if col is not None:
                feat[i, col] = 1.0
    print(f"[tiara] {int((feat.sum(1) > 0).sum())}/{len(contig_ids)} Contigs klassifiziert")
    return feat


def get_domain_groups(tiara_df: pd.DataFrame, contig_ids: list) -> dict:
    """Mappt Contigs auf Tiara-Domain → {domain: int_indices_array}.

    Domains: "eukarya" | "prokaryote" | "organelle" | "unknown"
    """
    groups: dict[str, list[int]] = {}
    for i, cid in enumerate(contig_ids):
        if cid in tiara_df.index:
            domain = str(tiara_df.loc[cid, "tiara_class"])
            if domain in ("organelle_mito", "organelle_plastid"):
                domain = "organelle"
        else:
            domain = "unknown"
        groups.setdefault(domain, []).append(i)
    return {k: np.array(v, dtype=np.int64) for k, v in groups.items()}

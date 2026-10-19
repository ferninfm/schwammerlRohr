from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from preprocessing import gc_content_length_table, coverage_table, process_kmer_matrix, tiara_node_features


def build_feature_table(
    fasta_path: Union[str, Path],
    coverage_tsv: Union[str, Path],
    kmer_path: Union[str, Path],
    use_canonical: bool = False,
    kmer_pca_components: int = 0,
    tiara_df: pd.DataFrame | None = None,
    use_length: bool = False,
) -> tuple[pd.DataFrame, np.ndarray, int, int]:
    """Builds X = [coverage(1), length(0|1), kmer(n), tiara(0|4)].
    Returns (feature_df, X, n_tiara, n_length).
    GC-content is kept in feature_df for plotting only — not used as model input.
    k-mer downweighting happens in gnn.py after standardization.
    Length is log1p-transformed and standardized normally (continuous feature).
    """
    # ── Coverage + GC-Content + Länge ───────────────────────────────────────
    gc_df = gc_content_length_table(fasta_path).rename(columns={"contig": "contig_id"})
    cov_df = coverage_table(coverage_tsv)
    feat = gc_df.merge(cov_df, on="contig_id", how="left")
    feat["coverage"] = np.log1p(feat["coverage"].fillna(0.0).astype(float))
    feat["length_log"] = np.log1p(feat["length"].fillna(0.0).astype(float))
    feat = feat[["contig_id", "gc_content", "coverage", "length_log"]].set_index("contig_id")

    # ── k-mer ────────────────────────────────────────────────────────────────
    kmer_df = process_kmer_matrix(kmer_path, canonical=use_canonical)
    n_raw = kmer_df.shape[1]
    mode = "kanonisch (136)" if use_canonical else "alle (256)"

    if kmer_pca_components > 0:
        pca = PCA(n_components=kmer_pca_components, random_state=42)
        arr = pca.fit_transform(kmer_df.values)
        print(f"[feature] k-mer PCA: {n_raw} ({mode}) → {kmer_pca_components} Komp. "
              f"({pca.explained_variance_ratio_.sum():.1%} Varianz)")
        kmer_df = pd.DataFrame(arr, index=kmer_df.index,
                               columns=[f"kmer_pca_{i}" for i in range(kmer_pca_components)])
    else:
        print(f"[feature] k-mer: {n_raw} Spalten ({mode})")

    feat = feat.join(kmer_df, how="left")
    feat[kmer_df.columns] = feat[kmer_df.columns].fillna(0.0)

    # ── Basis-X: coverage [+ length] + kmer ─────────────────────────────────
    n_length = 0
    base_cols = ["coverage"]
    if use_length:
        base_cols.append("length_log")
        n_length = 1
        print("[feature] Länge als Feature aktiv (log1p-transformiert)")
    X = feat[base_cols + list(kmer_df.columns)].to_numpy(dtype=np.float32)

    # ── Tiara (optional) ─────────────────────────────────────────────────────
    n_tiara = 0
    if tiara_df is not None:
        tiara_feat = tiara_node_features(feat.index, tiara_df)
        X = np.hstack([X, tiara_feat])
        n_tiara = tiara_feat.shape[1]

    return feat, X, n_tiara, n_length

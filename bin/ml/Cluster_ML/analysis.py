import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import hdbscan
import umap


# ── UMAP ──────────────────────────────────────────────────────────────────────

def run_umap(
    X: np.ndarray,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    metric: str = "euclidean",
) -> np.ndarray:
    """Standardisiert X und reduziert mit UMAP. min_dist=0.0 vor HDBSCAN empfohlen."""
    X_scaled = StandardScaler().fit_transform(X)
    reducer = umap.UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=42,
    )
    X_umap = reducer.fit_transform(X_scaled)
    print(f"[UMAP] {X.shape} → {X_umap.shape}")
    return X_umap


# ── HDBSCAN ───────────────────────────────────────────────────────────────────

def run_hdbscan(
    X_any_dim: np.ndarray,
    min_cluster_size: int = 10,
    min_samples: int = None,
    cluster_selection_method: str = "eom",
    cluster_selection_epsilon: float = 0.0,
    prescaled: bool = False,
    plot: bool = True,
    title: str = "",
    coords_df=None,
    save_path=None,
    label_offset: int = 0,
) -> np.ndarray:
    """HDBSCAN-Clustering. prescaled=True überspringt Skalierung (nach UMAP).

    cluster_selection_epsilon: Cluster näher als epsilon werden gemergt (0.0 = aus).
    coords_df: feature_df mit gc_content/coverage für zweiten Plot.
    label_offset: verschiebt Bin-IDs im Plot für Domain-Split.
    """
    X_scaled = X_any_dim if prescaled else StandardScaler().fit_transform(X_any_dim)
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method=cluster_selection_method,
        cluster_selection_epsilon=cluster_selection_epsilon,
    )
    labels    = clusterer.fit_predict(X_scaled)
    n_cluster = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise   = int(np.sum(labels == -1))
    print(f"HDBSCAN: {n_cluster} Cluster | {n_noise} Noise-Punkte")

    if not plot:
        return labels

    # Globale Bin-IDs für Plot (lokale labels + offset)
    display_labels = np.where(labels == -1, -1, labels + label_offset)
    unique_labels  = sorted(set(display_labels))
    cluster_ids    = [lbl for lbl in unique_labels if lbl != -1]

    cmap      = plt.cm.tab20
    color_map = {lbl: cmap(i / max(len(cluster_ids), 1)) for i, lbl in enumerate(cluster_ids)}
    color_map[-1] = (0.4, 0.4, 0.4, 1.0)

    legend_kwargs = dict(title="Cluster-ID", markerscale=1.5, fontsize=8,
                         title_fontsize=9, framealpha=0.8)

    # Punktgröße aus Contig-Länge (log1p → [8, 220] pt²)
    pt_sizes = None
    if coords_df is not None and "length_log" in coords_df.columns:
        len_vals = coords_df["length_log"].values.astype(float)
        lo, hi   = len_vals.min(), len_vals.max()
        pt_sizes = np.interp(len_vals, [lo, hi], [8, 220]) if hi > lo else np.full(len(len_vals), 30.0)

    def _scatter(ax, x_vals, y_vals, xlabel, ylabel, ax_title):
        for lbl in unique_labels:
            mask = display_labels == lbl
            s    = pt_sizes[mask] if pt_sizes is not None else 12
            ax.scatter(x_vals[mask], y_vals[mask], c=[color_map[lbl]], s=s,
                       alpha=0.55, linewidths=0,
                       label="Noise" if lbl == -1 else f"Bin {lbl}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(ax_title)

    # Embedding-Achsen: 2D → UMAP, sonst PCA zur Visualisierung
    if X_scaled.shape[1] == 2:
        emb_x, emb_y = X_scaled[:, 0], X_scaled[:, 1]
        emb_labels   = ("UMAP 1", "UMAP 2")
    else:
        X_2d         = PCA(n_components=2).fit_transform(X_scaled)
        emb_x, emb_y = X_2d[:, 0], X_2d[:, 1]
        emb_labels   = ("PC1 (viz)", "PC2 (viz)")

    if coords_df is not None:
        fig, (ax1, ax_leg, ax2) = plt.subplots(
            1, 3, figsize=(24, 8),
            gridspec_kw={"width_ratios": [5, 1.8, 5]},
        )
        fig.suptitle(f"HDBSCAN {title}".strip(), fontsize=13)
        _scatter(ax1, coords_df["gc_content"].values, coords_df["coverage"].values,
                 "GC Content", "log(Coverage)", "GC Content vs Coverage")
        _scatter(ax2, emb_x, emb_y, emb_labels[0], emb_labels[1],
                 f"{emb_labels[0]} vs {emb_labels[1]}")
        ax_leg.axis("off")
        handles, leg_labels = ax1.get_legend_handles_labels()
        ax_leg.legend(handles, leg_labels, loc="center", **legend_kwargs)
    else:
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.suptitle(f"HDBSCAN {title}".strip(), fontsize=13)
        _scatter(ax, emb_x, emb_y, emb_labels[0], emb_labels[1],
                 f"{emb_labels[0]} vs {emb_labels[1]}")
        ax.legend(loc="best", **legend_kwargs)

    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
        print(f"[plot] gespeichert: {save_path}")
    plt.show()
    plt.close(fig)
    return labels

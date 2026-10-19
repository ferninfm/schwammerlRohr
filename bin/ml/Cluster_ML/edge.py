from pathlib import Path
from typing import Dict, List, Optional, Set

import numpy as np
import pandas as pd
from taxoniq import Taxon


# ── Broad-Group via Taxoniq (für Cannot-Link) ─────────────────────────────────

_BROAD_GROUP_TAXIDS: Dict[int, str] = {
    4751:  "fungi",         # Fungi
    33090: "viridiplantae", # Viridiplantae
    2:     "bacteria",      # Bacteria
    2157:  "archaea",       # Archaea
}


def _broad_group_taxon(taxon: Taxon) -> str:
    """Bestimmt Breitgruppe via ranked_lineage (6 Knoten statt ~25 bei full lineage)."""
    try:
        for t in taxon.ranked_lineage:
            group = _BROAD_GROUP_TAXIDS.get(t.tax_id)
            if group:
                return group
    except Exception:
        pass
    return "unknown"


# ── Graph-Subset (für Domain-Split) ──────────────────────────────────────────

def subset_graph(
    X: np.ndarray,
    edge_index: np.ndarray,
    edge_weight: np.ndarray,
    global_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Schneidet X, edge_index und edge_weight auf eine Teilmenge von Knoten zu.

    global_indices: ganzzahlige Positionen in X (Form [n_sub]).
    Gibt (X_sub, edge_index_sub, edge_weight_sub) zurück, wobei edge_index_sub
    lokale Indizes (0 … n_sub-1) verwendet. Nur Kanten deren BEIDE Endpunkte
    in der Teilmenge liegen werden behalten.
    """
    n_total = X.shape[0]
    n_sub   = len(global_indices)

    g2l = np.full(n_total, -1, dtype=np.int64)
    g2l[global_indices] = np.arange(n_sub, dtype=np.int64)

    X_sub = X[global_indices]

    if edge_index.shape[1] == 0:
        return X_sub, np.zeros((2, 0), dtype=np.int64), np.zeros(0, dtype=np.float32)

    src_local = g2l[edge_index[0]]
    dst_local = g2l[edge_index[1]]
    keep = (src_local >= 0) & (dst_local >= 0)

    ei_sub = np.vstack([src_local[keep], dst_local[keep]])
    ew_sub = edge_weight[keep].astype(np.float32)
    return X_sub, ei_sub, ew_sub


# ── BAM-Links → Kanten ────────────────────────────────────────────────────────

def build_edge_index_from_links(
    links_tsv: str | Path,
    contig2idx: Dict[str, int],
    undirected: bool = True,
    contig_taxonomy: Optional[Dict[str, Taxon]] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Liest BAM-Links TSV (contig_a, contig_b, jaccard) → (edge_index, edge_weight).

    Wenn contig_taxonomy gegeben: Kanten zwischen verschiedenen Breitgruppen
    werden entfernt (Cannot-Link via ranked_lineage).
    """
    df = pd.read_csv(Path(links_tsv), sep="\t")
    src_list: List[int]   = []
    dst_list: List[int]   = []
    wgt_list: List[float] = []
    skipped    = 0
    cannot_lnk = 0

    # Breitgruppen-Cache: Taxon.tax_id → group-string (vermeidet Mehrfachlookups)
    group_cache: Dict[int, str] = {}

    def _cached_group(taxon: Taxon) -> str:
        tid = taxon.tax_id
        if tid not in group_cache:
            group_cache[tid] = _broad_group_taxon(taxon)
        return group_cache[tid]

    for _, row in df.iterrows():
        ca = str(row["contig_a"])
        cb = str(row["contig_b"])
        j  = float(row["jaccard"])

        if ca not in contig2idx or cb not in contig2idx:
            skipped += 1
            continue

        if contig_taxonomy is not None:
            tax_a = contig_taxonomy.get(ca)
            tax_b = contig_taxonomy.get(cb)
            if tax_a is not None and tax_b is not None:
                ga = _cached_group(tax_a)
                gb = _cached_group(tax_b)
                if ga != "unknown" and gb != "unknown" and ga != gb:
                    cannot_lnk += 1
                    continue

        src_list.append(contig2idx[ca])
        dst_list.append(contig2idx[cb])
        wgt_list.append(j)

    if skipped:
        print(f"[edge] {skipped} Links übersprungen (Contig nicht in Feature-Matrix)")
    if cannot_lnk:
        print(f"[edge] {cannot_lnk} Cannot-Link Kanten entfernt")

    if not src_list:
        print("[edge] WARNUNG: Keine BAM-Links geladen")
        return np.zeros((2, 0), dtype=np.int64), np.zeros(0, dtype=np.float32)

    src = np.array(src_list, dtype=np.int64)
    dst = np.array(dst_list, dtype=np.int64)
    wgt = np.array(wgt_list, dtype=np.float32)
    edge_index = np.vstack([src, dst])
    edge_weight = wgt

    if undirected:
        edge_index  = np.hstack([edge_index, edge_index[::-1, :]])
        edge_weight = np.concatenate([edge_weight, edge_weight])

    print(f"[edge] BAM-Links: {edge_index.shape[1]} Kanten ({len(src_list)} Links)")
    return edge_index, edge_weight


# ── Taxonomie Must-Link Kanten ────────────────────────────────────────────────

def build_must_link_edges(
    contig_taxonomy: Dict[str, Taxon],
    contig2idx: Dict[str, int],
    must_link_min_depth: int = 8,
    undirected: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Must-Link Kanten für Contig-Paare mit ausreichend tiefem gemeinsamem Vorfahren.

    Effizienz-Strategie:
      1. Contigs nach Konsens-Taxon (tax_id) gruppieren.
      2. Lineage jedes eindeutigen Taxons EINMAL als Set vorberechnen.
      3. LCA-Tiefe = len(b_lineage) - k  (Position im lineage-Array, kein Extra-Lookup).
      4. Kanten nur für Contig-Paare deren Taxon-Paar die Tiefenschwelle erfüllt.
      → O(n_unique_taxa²) LCA-Berechnungen statt O(n_contigs²).

    Tiefe = len(lca.lineage): größer = spezifischer (Taxoniq: spezifisch→root).
      depth ≥ 20 → Gewicht 1.00  (Gattungs-/Artebene)
      depth ≥ 15 → Gewicht 0.85  (Familien-/Ordnungsebene)
      depth ≥ min_depth → Gewicht 0.60
    """
    # ── Schritt 1: Contigs nach Konsens-Taxon gruppieren ─────────────────────
    taxon_groups: Dict[int, List[int]] = {}  # tax_id → [node_idx, ...]
    taxon_cache:  Dict[int, Taxon]     = {}  # tax_id → Taxon

    for cid, taxon in contig_taxonomy.items():
        if cid not in contig2idx:
            continue
        tid = taxon.tax_id
        taxon_groups.setdefault(tid, []).append(contig2idx[cid])
        taxon_cache[tid] = taxon

    unique_tids = list(taxon_groups.keys())
    n_unique    = len(unique_tids)

    # ── Schritt 2: Lineage EINMAL pro Taxon vorberechnen ─────────────────────
    # lineage_sets: schnelles Membership-Lookup für LCA
    # lineage_lists: für positionsbasierte Tiefenberechnung ohne Extra-Lookup
    # taxon_depth: Eigentiefe für Same-Group-Fall
    lineage_sets:  Dict[int, Set[int]]    = {}
    lineage_lists: Dict[int, List[Taxon]] = {}
    taxon_depth:   Dict[int, int]         = {}

    for tid, taxon in taxon_cache.items():
        lin = list(taxon.lineage)
        lineage_lists[tid] = lin
        lineage_sets[tid]  = {t.tax_id for t in lin}
        taxon_depth[tid]   = len(lin)

    # ── Schritt 3: LCA-Tiefe pro Taxon-Paar + Kanten generieren ──────────────
    src_list: List[int]   = []
    dst_list: List[int]   = []
    wgt_list: List[float] = []

    for i in range(n_unique):
        tid_a   = unique_tids[i]
        ids_a   = lineage_sets[tid_a]
        idx_a   = taxon_groups[tid_a]

        for j in range(i, n_unique):
            tid_b = unique_tids[j]
            idx_b = taxon_groups[tid_b]

            if i == j:
                # LCA = Taxon selbst → Tiefe direkt aus vorberechneter Tiefe
                lca_depth = taxon_depth[tid_a]
            else:
                # LCA-Tiefe: Position k in b_lineage wo b's Vorfahre in a_ids auftaucht
                # → len(b_lineage) - k  (identisch mit len(lca_node.lineage), kein Extra-Lookup)
                b_lin     = lineage_lists[tid_b]
                lca_depth = 1  # root-Fallback
                for k, node in enumerate(b_lin):
                    if node.tax_id in ids_a:
                        lca_depth = len(b_lin) - k
                        break

            if lca_depth < must_link_min_depth:
                continue

            weight = 1.0 if lca_depth >= 20 else (0.85 if lca_depth >= 15 else 0.60)

            if i == j:
                # Innerhalb einer Gruppe: alle Paare (p < q)
                for p in range(len(idx_a)):
                    for q in range(p + 1, len(idx_a)):
                        src_list.append(idx_a[p])
                        dst_list.append(idx_a[q])
                        wgt_list.append(weight)
            else:
                # Zwischen zwei Gruppen: vollständiges Kreuzprodukt
                for na in idx_a:
                    for nb in idx_b:
                        src_list.append(na)
                        dst_list.append(nb)
                        wgt_list.append(weight)

    if not src_list:
        print(f"[edge] Taxonomy Must-Link: keine Kanten (min_depth={must_link_min_depth})")
        return np.zeros((2, 0), dtype=np.int64), np.zeros(0, dtype=np.float32)

    src = np.array(src_list, dtype=np.int64)
    dst = np.array(dst_list, dtype=np.int64)
    wgt = np.array(wgt_list, dtype=np.float32)
    edge_index = np.vstack([src, dst])
    edge_weight = wgt

    if undirected:
        edge_index  = np.hstack([edge_index, edge_index[::-1, :]])
        edge_weight = np.concatenate([edge_weight, edge_weight])

    print(f"[edge] Taxonomy Must-Link: {edge_index.shape[1]} Kanten "
          f"({len(src_list)} Paare | {n_unique} unique Taxa | "
          f"min_depth={must_link_min_depth})")
    return edge_index, edge_weight

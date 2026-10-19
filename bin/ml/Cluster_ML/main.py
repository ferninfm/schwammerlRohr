import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from feature import build_feature_table
from analysis import run_umap, run_hdbscan
from gnn import train_gae_embeddings
from postprocessing_taxonomy import (load_taxonomy_hits, get_contig_consensus_taxonomy,
                                     assign_bin_taxonomy)
from preprocessing import load_tiara, get_domain_groups
from edge import build_edge_index_from_links, build_must_link_edges, subset_graph


def main():

    # ── Pfade ─────────────────────────────────────────────────────────────────
    DATA_DIR = Path("cluster_data")
    SAMPLE   = DATA_DIR / "pr_197_002"
    INPUT    = SAMPLE / "input"
    FASTA    = INPUT / "assembly.fasta"
    COVERAGE = INPUT / "mean_coverage.tsv"
    KMER     = INPUT / "assembly.k4.per_contig_matrix.csv"
    LINKS    = INPUT / "links.tsv"
    TIARA    = INPUT / "complete.classification.tsv"

    # Taxonomie-Eingabepfade pro Tool (automatisch aus Sample-Name aufgebaut)
    _S = SAMPLE.name
    TAX_PATHS = {
        "diamond_blastp":  INPUT / f"{_S}_taxprot",
        "diamond_blastx":  INPUT / f"{_S}_diamond_blastx.tsv",
        "centrifuger":     INPUT / f"{_S}_classification_centrifuger.tsv",
        "kraken2_pf":      INPUT / f"{_S}_pf.kraken2.tsv",
        "kraken2_pluspfp": INPUT / f"{_S}.kraken2.pluspfp.tsv",
    }

    # ── Konfiguration ─────────────────────────────────────────────────────────
    USE_CANONICAL_KMERS  = True
    KMER_PCA_COMPONENTS  = 25
    KMER_WEIGHT          = 0.4
    MUST_LINK_MIN_DEPTH  = 8
    CONFIDENCE_MIN_DEPTH = 8   # Min. gemeinsame Vorfahren-Tiefe für Confidence

    USE_TIARA            = True
    USE_LENGTH           = None
    USE_TAX_CANNOT_LINK  = True
    USE_TAX_MUST_LINK    = True
    USE_TAX_BIN_TAXONOMY = True
    USE_DOMAIN_SPLIT     = None
    MIN_DOMAIN_SIZE      = 30

    # "diamond_blastp" | "diamond_blastx" | "centrifuger" | "kraken2_pf" | "kraken2_pluspfp"
    TAX_TOOL = "diamond_blastp"

    # HDBSCAN: epsilon > 0 mergt nahe Cluster; min_samples > 1 macht Clustering konservativer
    HDBSCAN_EPSILON           = 0.7
    HDBSCAN_MIN_SAMPLES       = 1
    HDBSCAN_MIN_CLUSTER_SIZE  = 24 

    # ── Output-Ordner (pro Tool + Domain-Split-Modus ein Unterordner) ───────────
    _MODE   = "domain_split" if USE_DOMAIN_SPLIT else "no_domain_split"
    _TAG    = f"{_S}__{TAX_TOOL}__{_MODE}"
    OUTPUT  = SAMPLE / "output" / TAX_TOOL / _MODE
    PLOTS   = OUTPUT / "plots"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)

    # ── Console-Output in Datei spiegeln (mit Zeitstempel pro Zeile) ─────────
    class _Tee:
        def __init__(self, fh):
            self._fh        = fh
            self._out       = sys.__stdout__
            self._bol       = True   # beginning-of-line flag
        def write(self, msg):
            self._out.write(msg)
            if not msg:
                return
            stamped = []
            for ch in msg:
                if self._bol and ch != "\n":
                    stamped.append(f"[{datetime.now().strftime('%H:%M:%S')}] ")
                    self._bol = False
                stamped.append(ch)
                if ch == "\n":
                    self._bol = True
            self._fh.write("".join(stamped))
        def flush(self):
            self._out.flush()
            self._fh.flush()

    _log_path = OUTPUT / f"{_TAG}__console_output.txt"
    _log_fh   = open(_log_path, "w", encoding="utf-8")
    sys.stdout = _Tee(_log_fh)

    print(f"[output] {OUTPUT}")

    # ── Tiara laden ───────────────────────────────────────────────────────────
    tiara_df = None
    if USE_TIARA:
        if TIARA.exists():
            print("\n── Tiara ────────────────────────────────────────────────")
            tiara_df = load_tiara(TIARA)
        else:
            print(f"[tiara] WARNUNG: {TIARA} nicht gefunden")

    # ── Taxonomie laden ───────────────────────────────────────────────────────
    tax_hits        = None
    contig_taxonomy = None
    need_taxonomy   = USE_TAX_CANNOT_LINK or USE_TAX_MUST_LINK or USE_TAX_BIN_TAXONOMY

    if need_taxonomy:
        TAX_INPUT = TAX_PATHS.get(TAX_TOOL.lower())
        if TAX_INPUT is None:
            print(f"[taxonomy] WARNUNG: Unbekanntes TAX_TOOL='{TAX_TOOL}'")
        elif not TAX_INPUT.exists():
            print(f"[taxonomy] WARNUNG: {TAX_INPUT} nicht gefunden")
        else:
            print(f"\n── Taxonomie ({TAX_TOOL.upper()}) ───────────────────────────────")
            tax_hits = load_taxonomy_hits(TAX_INPUT, tool=TAX_TOOL)
            if USE_TAX_CANNOT_LINK or USE_TAX_MUST_LINK:
                contig_taxonomy = get_contig_consensus_taxonomy(tax_hits)

    # ── Feature-Matrix ────────────────────────────────────────────────────────
    print("\n── Feature-Matrix ───────────────────────────────────────")
    feature_df, X, n_tiara, n_length = build_feature_table(
        fasta_path=FASTA,
        coverage_tsv=COVERAGE,
        kmer_path=KMER,
        use_canonical=USE_CANONICAL_KMERS,
        kmer_pca_components=KMER_PCA_COMPONENTS,
        tiara_df=tiara_df,
        use_length=USE_LENGTH,
    )
    print(f"X: {X.shape[0]} Contigs × {X.shape[1]} Features "
          f"({n_length} Länge, {n_tiara} Tiara)")

    # ── Organellen herausfiltern ──────────────────────────────────────────────
    organelle_feature_df = pd.DataFrame()
    organelle_ids: set   = set()

    if tiara_df is not None:
        for cid in feature_df.index:
            if cid in tiara_df.index:
                if str(tiara_df.loc[cid, "tiara_class"]) in ("organelle_mito", "organelle_plastid"):
                    organelle_ids.add(cid)

    if organelle_ids:
        organelle_feature_df = feature_df.loc[sorted(organelle_ids)].copy()
        keep_mask  = np.array([cid not in organelle_ids for cid in feature_df.index])
        feature_df = feature_df[keep_mask]
        X          = X[keep_mask]
        print(f"\n[organelle] {len(organelle_ids)} entfernt → {len(feature_df)} für GNN")

    # ── Graph ─────────────────────────────────────────────────────────────────
    print("\n── Graph ────────────────────────────────────────────────")
    contig2idx = {cid: i for i, cid in enumerate(feature_df.index)}

    edge_index, edge_weight = build_edge_index_from_links(
        links_tsv=LINKS,
        contig2idx=contig2idx,
        contig_taxonomy=contig_taxonomy if USE_TAX_CANNOT_LINK else None,
    )

    if USE_TAX_MUST_LINK and contig_taxonomy is not None:
        di_index, di_weight = build_must_link_edges(
            contig_taxonomy=contig_taxonomy,
            contig2idx=contig2idx,
            must_link_min_depth=MUST_LINK_MIN_DEPTH,
        )
        edge_index  = np.hstack([edge_index, di_index])
        edge_weight = np.concatenate([edge_weight, di_weight])
        print(f"[edge] Gesamt: {edge_index.shape[1]} Kanten")

    contig_ids = list(feature_df.index)

    # ── GNN + UMAP + HDBSCAN ─────────────────────────────────────────────────
    def _run_gnn_umap_hdbscan(X_in, ei_in, ew_in, n_nbrs, min_cs,
                               coords, title_tag, save_path=None, label_offset=0):
        """GNN → UMAP → HDBSCAN. Gibt (labels, Z_umap) zurück."""
        print(f"\n── GNN {title_tag} ──────────────────────────────────────────")
        Z_out = train_gae_embeddings(
            X=X_in, edge_index=ei_in, edge_weight=ew_in,
            out_dim=32, hidden=64, epochs=1000, lr=5e-4,
            alpha=0.6, dropout=0.1,
            kmer_weight=KMER_WEIGHT,
            n_tiara_features=n_tiara,
            n_length_features=n_length,
        )
        Z_u = run_umap(Z_out, n_components=2, n_neighbors=n_nbrs, min_dist=0.0)
        print(f"\n── HDBSCAN {title_tag} ─────────────────────────────────────")
        lbl = run_hdbscan(
            Z_u, min_cluster_size=min_cs, min_samples=HDBSCAN_MIN_SAMPLES,
            cluster_selection_epsilon=HDBSCAN_EPSILON,
            cluster_selection_method="eom", title=title_tag,
            coords_df=coords, save_path=save_path, label_offset=label_offset,
            prescaled=True,
        )
        return lbl, Z_u

    # ── Binning ───────────────────────────────────────────────────────────────
    domain_groups: dict | None = None

    if USE_DOMAIN_SPLIT and tiara_df is not None:
        print("\n── Domain-Split Binning ─────────────────────────────────")
        domain_groups = get_domain_groups(tiara_df, contig_ids)
        domain_groups = {k: v for k, v in domain_groups.items()
                         if k not in ("organelle", "organelle_mito", "organelle_plastid")}

        # Taxonomie-Reklassifikation: Unknown-Contigs mit Treffer → richtiger Domain
        if contig_taxonomy is not None and "unknown" in domain_groups:
            _BACTERIA_ID  = 2
            _ARCHAEA_ID   = 2157
            _EUKARYOTA_ID = 2759
            _reclassified = {"prokaryote": [], "eukarya": []}
            _still_unknown = []

            for _idx in domain_groups["unknown"]:
                _cid = contig_ids[_idx]
                if _cid in contig_taxonomy:
                    _lin_ids = {t.tax_id for t in contig_taxonomy[_cid].lineage}
                    if _BACTERIA_ID in _lin_ids or _ARCHAEA_ID in _lin_ids:
                        _reclassified["prokaryote"].append(_idx)
                    elif _EUKARYOTA_ID in _lin_ids:
                        _reclassified["eukarya"].append(_idx)
                    else:
                        _still_unknown.append(_idx)
                else:
                    _still_unknown.append(_idx)

            for _dom, _idx_list in _reclassified.items():
                if _idx_list:
                    domain_groups[_dom] = list(domain_groups.get(_dom, [])) + _idx_list
                    print(f"[domain] {len(_idx_list)} Unknown → {_dom} (Taxonomie)")

            if _still_unknown:
                domain_groups["unknown"] = _still_unknown
                print(f"[domain] {len(_still_unknown)} Contigs verbleiben als Unknown")
            else:
                del domain_groups["unknown"]
                print("[domain] Alle Unknown reklassifiziert")

        labels     = np.full(len(contig_ids), -1, dtype=np.int64)
        next_label = 0

        for domain_name, dom_idx in domain_groups.items():
            n_dom = len(dom_idx)
            if n_dom < MIN_DOMAIN_SIZE:
                print(f"[domain] {domain_name}: {n_dom} Contigs → übersprungen")
                continue

            X_sub, ei_sub, ew_sub = subset_graph(X, edge_index, edge_weight, dom_idx)
            n_nbrs     = min(30, n_dom - 1)
            min_cs     = max(5, min(24, n_dom // 10))
            safe_tag   = domain_name.replace(" ", "_")
            plot_file  = PLOTS / f"{_TAG}__hdbscan_{safe_tag}.png"
            plot_title = f"{_S}  |  {TAX_TOOL}  |  {domain_name}"

            sub_lbl, _ = _run_gnn_umap_hdbscan(
                X_sub, ei_sub, ew_sub, n_nbrs, min_cs,
                feature_df.iloc[dom_idx], plot_title,
                save_path=plot_file, label_offset=next_label,
            )
            max_sub = int(np.max(sub_lbl)) if np.any(sub_lbl >= 0) else -1
            for local_i, global_i in enumerate(dom_idx):
                if sub_lbl[local_i] != -1:
                    labels[global_i] = sub_lbl[local_i] + next_label
            next_label += max_sub + 1 if max_sub >= 0 else 0

        print(f"\n── Domain-Split: {next_label} Bins gesamt ───────────────────")
        print(pd.Series(labels).value_counts())

    else:
        plot_file  = PLOTS / f"{_TAG}__hdbscan.png"
        plot_title = f"{_S}  |  {TAX_TOOL}  |  Gesamt"
        labels, _  = _run_gnn_umap_hdbscan(
            X, edge_index, edge_weight, n_nbrs=30, min_cs=HDBSCAN_MIN_CLUSTER_SIZE,
            coords=feature_df, title_tag=plot_title, save_path=plot_file,
        )
        print(pd.Series(labels).value_counts())

    # ── Contig-Domain-Mapping ─────────────────────────────────────────────────
    _contig_domain = np.full(len(contig_ids), "unknown", dtype=object)
    if USE_DOMAIN_SPLIT and domain_groups is not None:
        for _dom, _idx in domain_groups.items():
            _contig_domain[_idx] = _dom

    # ── Bin-Taxonomie ─────────────────────────────────────────────────────────
    tax_df = pd.DataFrame()
    if USE_TAX_BIN_TAXONOMY and tax_hits is not None:
        print(f"\n── Bin-Taxonomie ({TAX_TOOL.upper()}) ───────────────────────────────")
        tax_df = assign_bin_taxonomy(
            labels=labels, contig_ids=contig_ids, tax_hits=tax_hits,
            n_total_contigs=len(contig_ids),
            confidence_min_depth=CONFIDENCE_MIN_DEPTH,
        )

        _lbl_arr = np.array(labels, dtype=np.int64)

        def _bin_majority_domain(bin_label_str: str) -> str:
            if bin_label_str == "Noise":
                mask = (_lbl_arr == -1)
            else:
                try:
                    mask = (_lbl_arr == int(bin_label_str.split()[-1]))
                except ValueError:
                    return "unknown"
            if not mask.any():
                return "unknown"
            vals, counts = np.unique(_contig_domain[mask], return_counts=True)
            return str(vals[counts.argmax()])

        tax_df["domain"] = tax_df["bin_label"].apply(_bin_majority_domain)

        _lines = []
        if USE_DOMAIN_SPLIT:
            _cols = ["bin_label", "domain", "n_contigs", "n_contigs_with_hits",
                     "taxonomy_short", "taxonomy_rank", "taxonomy_depth", "confidence"]
            for _dom_key, _header in [("prokaryote", "── Prokaryoten"),
                                       ("eukarya",    "── Eukaryoten"),
                                       ("unknown",    "── Unknown / Andere")]:
                _sub = tax_df[
                    (tax_df["domain"] == _dom_key) &
                    (tax_df["bin_label"] != "Noise")
                ][_cols].sort_values("confidence", ascending=False)
                if not _sub.empty:
                    _hdr = f"\n{_header} ──────────────────────────────────────────"
                    _tbl = _sub.to_string(index=False)
                    print(_hdr)
                    print(_tbl)
                    _lines += [_hdr, _tbl]
            _noise = tax_df[tax_df["bin_label"] == "Noise"][_cols]
            if not _noise.empty:
                _hdr = "\n── Noise ──────────────────────────────────────────────"
                _tbl = _noise.to_string(index=False)
                print(_hdr)
                print(_tbl)
                _lines += [_hdr, _tbl]
        else:
            _cols  = ["bin_label", "n_contigs", "n_contigs_with_hits",
                      "taxonomy_short", "taxonomy_rank", "taxonomy_depth", "confidence"]
            tax_df = tax_df.drop(columns=["domain"])
            _hdr = "\n── Gesamt ──────────────────────────────────────────────"
            _tbl = tax_df[_cols].sort_values("confidence", ascending=False).to_string(index=False)
            print(_hdr)
            print(_tbl)
            _lines += [_hdr, _tbl]



    # ── Bin-Taxonomie-Lookup ──────────────────────────────────────────────────
    bin_tax_short_map: dict = {}
    if not tax_df.empty:
        for _, row in tax_df.iterrows():
            bin_tax_short_map[int(row["bin_id"])] = row["taxonomy_short"]

    # ── Taxonomie-Helfer pro Contig ───────────────────────────────────────────
    def _tax_path(cid: str) -> str:
        if contig_taxonomy and cid in contig_taxonomy:
            return "; ".join(t.scientific_name for t in reversed(contig_taxonomy[cid].lineage))
        return ""

    def _tax_short(cid: str) -> str:
        if contig_taxonomy and cid in contig_taxonomy:
            return contig_taxonomy[cid].scientific_name
        return ""

    # ── Binning Results aufbauen ──────────────────────────────────────────────
    rows = []
    for i, cid in enumerate(contig_ids):
        lbl = int(labels[i])
        rows.append({
            "contig_id":     cid,
            "bin_id":        lbl,
            "bin_label":     "Noise" if lbl == -1 else f"Bin {lbl}",
            "domain":        _contig_domain[i] if USE_DOMAIN_SPLIT else "",
            "bin_tax_short": bin_tax_short_map.get(lbl, ""),
            "tax_path":      _tax_path(cid),
            "tax_short":     _tax_short(cid),
        })

    for cid in organelle_feature_df.index:
        org_cls = str(tiara_df.loc[cid, "tiara_class"]) if (
            tiara_df is not None and cid in tiara_df.index) else "organelle"
        rows.append({
            "contig_id":     cid,
            "bin_id":        -2,
            "bin_label":     "Organelle",
            "domain":        org_cls,
            "bin_tax_short": "",
            "tax_path":      _tax_path(cid),
            "tax_short":     _tax_short(cid),
        })

    results_df = pd.DataFrame(rows)

    # ── Sortierung ────────────────────────────────────────────────────────────
    if USE_DOMAIN_SPLIT:
        _DOMAIN_ORDER = {"prokaryote": 0, "eukarya": 1, "unknown": 2,
                         "organelle_mito": 3, "organelle_plastid": 3}
        results_df["_sort"] = results_df["domain"].map(lambda d: _DOMAIN_ORDER.get(d, 2))
        results_df = (results_df
                      .sort_values(["_sort", "bin_id", "contig_id"])
                      .drop(columns=["_sort"])
                      .reset_index(drop=True))
    else:
        results_df = results_df.sort_values(["bin_id", "contig_id"]).reset_index(drop=True)

    # ── Binning Results schreiben ─────────────────────────────────────────────
    cols = ["contig_id", "bin_id", "bin_label", "domain", "bin_tax_short", "tax_path", "tax_short"]
    if not USE_DOMAIN_SPLIT:
        cols = [c for c in cols if c != "domain"]

    def _sec(title: str) -> str:
        return f"\n{'=' * 60}\n  {title}\n{'=' * 60}\n"

    def _bin(title: str) -> str:
        return f"\n{'-' * 50} {title} {'-' * 50}\n"

    out_path = OUTPUT / f"{_TAG}__binning_results.tsv"
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")

        if USE_DOMAIN_SPLIT:
            written_org = False
            for dom_key, dom_title in [("prokaryote", "PROKARYOTEN"),
                                        ("eukarya",    "EUKARYOTEN"),
                                        ("unknown",    "UNKNOWN"),
                                        ("organelle_mito",   "ORGANELLEN"),
                                        ("organelle_plastid", "ORGANELLEN")]:
                dom_df = results_df[results_df["domain"] == dom_key]
                if dom_df.empty:
                    continue
                if dom_key in ("organelle_mito", "organelle_plastid"):
                    if not written_org:
                        fh.write(_sec("ORGANELLEN"))
                        written_org = True
                else:
                    fh.write(_sec(dom_title))
                for bin_lbl, bin_grp in dom_df.groupby("bin_label", sort=False):
                    fh.write(_bin(bin_lbl))
                    for _, r in bin_grp.iterrows():
                        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        else:
            noise_df = results_df[results_df["bin_label"] == "Noise"]
            bins_df  = results_df[results_df["bin_label"].str.startswith("Bin")]
            org_df   = results_df[results_df["bin_label"] == "Organelle"]
            for bin_lbl, bin_grp in bins_df.groupby("bin_id"):
                fh.write(_bin(f"Bin {bin_lbl}"))
                for _, r in bin_grp.iterrows():
                    fh.write("\t".join(str(r[c]) for c in cols) + "\n")
            if not noise_df.empty:
                fh.write(_bin("Noise"))
                for _, r in noise_df.iterrows():
                    fh.write("\t".join(str(r[c]) for c in cols) + "\n")
            if not org_df.empty:
                fh.write(_sec("ORGANELLEN"))
                fh.write(_bin("Organellen"))
                for _, r in org_df.iterrows():
                    fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"→ {out_path}")
    print(f"→ Plots: {PLOTS}/")

    # ── Console-Log schließen ─────────────────────────────────────────────────
    sys.stdout = sys.__stdout__
    _log_fh.close()
    print(f"→ Console-Log: {_log_path}")


if __name__ == "__main__":
    main()

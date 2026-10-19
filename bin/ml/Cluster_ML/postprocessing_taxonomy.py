from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from taxoniq import Taxon


# ── Taxonomie-Konsens ─────────────────────────────────────────────────────────

def _majority_consensus(
    tax_ids: List[int],
    min_fraction: float = 0.5,
    taxon_cache: Optional[Dict[int, Taxon]] = None,
) -> Optional[Taxon]:
    """Tiefster NCBI-Knoten den >= min_fraction der Hits teilen (Majority-LCA).

    taxon_cache: gemeinsamer Cache über mehrere Aufrufe → ein DB-Lookup pro tax_id.
    Tiefe wird inline aus Lineage-Position berechnet.
    """
    if taxon_cache is None:
        taxon_cache = {}

    valid: List[Taxon] = []
    for tid in tax_ids:
        if not tid:
            continue
        if tid not in taxon_cache:
            try:
                taxon_cache[tid] = Taxon(int(tid))
            except Exception:
                continue
        valid.append(taxon_cache[tid])

    if not valid:
        return None

    threshold    = len(valid) * min_fraction
    votes:       Dict[int, int]   = {}
    anc_obj:     Dict[int, Taxon] = {}
    depth_cache: Dict[int, int]   = {}

    for taxon in valid:
        lin = taxon.lineage
        n   = len(lin)
        for k, ancestor in enumerate(lin):
            tid = ancestor.tax_id
            votes[tid] = votes.get(tid, 0) + 1
            if tid not in depth_cache:
                depth_cache[tid] = n - k
                anc_obj[tid]     = ancestor

    best: Optional[Taxon] = None
    best_depth = -1
    for tid, count in votes.items():
        if count < threshold:
            continue
        depth = depth_cache[tid]
        if depth > best_depth:
            best_depth = depth
            best = anc_obj[tid]

    return best


# ── Tool-spezifische Loader ───────────────────────────────────────────────────

def _parse_contig_id(query: str) -> Optional[str]:
    """Extrahiert contig_id aus Diamond Query-Format: ...|contig_XXX:start-end|..."""
    for part in query.split("|"):
        if part.startswith("contig_"):
            return part.split(":")[0]
    return None


def load_diamond(
    path: str | Path,
    evalue_threshold: float = 1e-10,
) -> Dict[str, List[int]]:
    """Liest Diamond BLASTP Output → {contig_id: [tax_id, ...]}.

    Format (Tab-separiert): query | tax_id | evalue | taxonomy_string
    query enthält die contig_id im Format: ...|contig_XXX:start-end|...
    """
    df = pd.read_csv(Path(path), sep="\t", header=None,
                     names=["query", "tax_id", "evalue", "taxonomy"])
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce").fillna(1.0)
    df["tax_id"] = pd.to_numeric(df["tax_id"], errors="coerce").fillna(0).astype(int)

    n_before = len(df)
    df = df[(df["evalue"] <= evalue_threshold) & (df["tax_id"] > 0)]
    print(f"[taxonomy] {len(df)}/{n_before} Hits nach E-value + tax_id Filter")

    df["contig_id"] = df["query"].apply(_parse_contig_id)
    df = df.dropna(subset=["contig_id"])

    result: Dict[str, List[int]] = {}
    for contig_id, group in df.groupby("contig_id"):
        ids = group["tax_id"].tolist()
        if ids:
            result[str(contig_id)] = ids

    print(f"[taxonomy] {len(result)} Contigs mit Hits geladen")
    return result


def _fetch_accession_taxids(
    accessions: List[str],
    batch_size: int = 200,
    sleep_sec: float = 0.35,
) -> Dict[str, int]:
    """Fragt NCBI eutils nach tax_ids für Protein-Accession-Nummern.

    Gibt {accession: tax_id} zurück. Accessions ohne Treffer werden
    weggelassen.
    """
    import time
    try:
        import requests
    except ImportError:
        print("[taxonomy] WARNUNG: 'requests' nicht installiert — "
              "pip install requests. Kein tax_id-Lookup möglich.")
        return {}

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    acc2taxid: Dict[str, int] = {}

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        try:
            resp = requests.get(
                f"{base}/esummary.fcgi",
                params={"db": "protein", "id": ",".join(batch), "retmode": "json"},
                timeout=30,
            )
            resp.raise_for_status()
            data = resp.json().get("result", {})
            for uid, doc in data.items():
                if uid == "uids" or not isinstance(doc, dict):
                    continue
                taxid = int(doc.get("taxid", 0))
                acc_ver = doc.get("accessionversion", "")
                if taxid > 0 and acc_ver:
                    acc2taxid[acc_ver] = taxid
                    acc2taxid[acc_ver.split(".")[0]] = taxid
        except Exception as e:
            print(f"[taxonomy] NCBI-Anfrage fehlgeschlagen (Batch {i}): {e}")
        time.sleep(sleep_sec)

    return acc2taxid


def load_diamond_blastx(
    path: str | Path,
    evalue_threshold: float = 1e-10,
    cache: bool = True,
) -> Dict[str, List[int]]:
    """Liest Diamond BLASTX Output (outfmt 6) → {contig_id: [tax_id, ...]}.

    Format (12 Spalten): qseqid sseqid pident length mismatch gapopen
    qstart qend sstart send evalue bitscore

    Da outfmt 6 keine tax_id enthält, werden die Protein-Accessions über
    die NCBI eutils API zu tax_ids aufgelöst. Das Ergebnis wird als
    JSON-Cache neben der Eingabedatei gespeichert, damit spätere Runs
    die API nicht erneut aufrufen.
    """
    import json

    path = Path(path)
    cache_path = path.with_suffix(".taxid_cache.json")

    df = pd.read_csv(path, sep="\t", header=None,
                     names=["contig_id", "sseqid", "pident", "length",
                             "mismatch", "gapopen", "qstart", "qend",
                             "sstart", "send", "evalue", "bitscore"])
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce").fillna(1.0)

    n_before = len(df)
    df = df[df["evalue"] <= evalue_threshold]
    print(f"[taxonomy] {len(df)}/{n_before} Hits nach E-value Filter")

    unique_accs = df["sseqid"].dropna().unique().tolist()

    # Cache laden oder NCBI abfragen
    if cache and cache_path.exists():
        with open(cache_path) as f:
            acc2taxid = json.load(f)
        print(f"[taxonomy] tax_id Cache geladen: {cache_path.name} "
              f"({len(acc2taxid)} Einträge)")
    else:
        print(f"[taxonomy] Frage NCBI nach tax_ids für "
              f"{len(unique_accs)} Accessions …")
        acc2taxid = _fetch_accession_taxids(unique_accs)
        print(f"[taxonomy] {len(acc2taxid)} tax_ids erhalten")
        if cache:
            with open(cache_path, "w") as f:
                json.dump(acc2taxid, f)
            print(f"[taxonomy] Cache gespeichert: {cache_path.name}")

    df["tax_id"] = df["sseqid"].map(acc2taxid).fillna(0).astype(int)
    df = df[df["tax_id"] > 0]

    result: Dict[str, List[int]] = {}
    for contig_id, group in df.groupby("contig_id"):
        ids = group["tax_id"].tolist()
        if ids:
            result[str(contig_id)] = ids

    print(f"[taxonomy] {len(result)} Contigs mit tax_id-Hits geladen")
    return result


def load_centrifuger(
    path: str | Path,
    min_hit_length: int = 100,
) -> Dict[str, List[int]]:
    """Liest Centrifuger TSV Output → {contig_id: [tax_id]}.

    Format: readID seqID taxID score 2ndBestScore hitLength queryLength numMatches
    Filter: taxID > 1, taxID != 32630 (synthetic), hitLength >= min_hit_length
    """
    df = pd.read_csv(Path(path), sep="\t")
    df.columns = [c.strip() for c in df.columns]

    n_before = len(df)
    df = df[
        (df["taxID"] > 1) &
        (df["taxID"] != 32630) &
        (df["hitLength"] >= min_hit_length)
    ]
    print(f"[taxonomy] {len(df)}/{n_before} Hits nach taxID + hitLength Filter")

    result: Dict[str, List[int]] = {}
    for _, row in df.iterrows():
        result[str(row["readID"])] = [int(row["taxID"])]

    print(f"[taxonomy] {len(result)} Contigs mit Hits geladen")
    return result


def load_kraken2(
    path: str | Path,
) -> Dict[str, List[int]]:
    """Liest Kraken2 Output → {contig_id: [tax_id]}.

    Format: C/U  sequence_id  tax_id  length  kmer_hits
    Unklassifizierte (U / taxID 0) werden ausgeschlossen — gleich wie Diamond/Centrifuger.
    Kein weiteres Filtern: rohe Klassifikation bleibt erhalten (inkl. Homo sapiens).
    """
    result: Dict[str, List[int]] = {}
    n_total      = 0
    n_classified = 0

    with open(Path(path)) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            n_total += 1
            status = parts[0]
            cid    = str(parts[1])
            try:
                tid = int(parts[2])
            except ValueError:
                continue
            if status == "U" or tid == 0:
                continue
            n_classified += 1
            result[cid] = [tid]

    print(f"[taxonomy] Kraken2: {n_classified} klassifiziert | "
          f"{n_total - n_classified} unklassifiziert von {n_total} Contigs")
    return result


def load_taxonomy_hits(
    path: str | Path,
    tool: str = "diamond",
    **kwargs,
) -> Dict[str, List[int]]:
    """Dispatcher → wählt Loader anhand von tool.

    tool: "diamond" | "diamond_blastp" | "diamond_blastx"
          "centrifuger"
          "kraken2" | "kraken2_pf" | "kraken2_pluspfp"
    """
    tool = tool.lower().strip()
    if tool in ("diamond", "diamond_blastp"):
        return load_diamond(path, **kwargs)
    if tool == "diamond_blastx":
        return load_diamond_blastx(path, **kwargs)
    if tool == "centrifuger":
        return load_centrifuger(path, **kwargs)
    if tool in ("kraken2", "kraken2_pf", "kraken2_pluspfp"):
        return load_kraken2(path, **kwargs)
    raise ValueError(f"Unbekanntes Tool: '{tool}'. "
                     f"Erlaubt: diamond, centrifuger, kraken2_pf, kraken2_pluspfp")


# ── Konsens-Taxonomie pro Contig ──────────────────────────────────────────────

def get_contig_consensus_taxonomy(
    tax_hits: Dict[str, List[int]],
    min_fraction: float = 0.5,
) -> Dict[str, Taxon]:
    """Majority-LCA pro Contig → {contig_id: Taxon}.

    Gemeinsamer taxon_cache: jede tax_id nur einmal in Taxoniq aufgelöst.
    """
    shared_cache: Dict[int, Taxon] = {}
    result: Dict[str, Taxon] = {}
    for cid, tax_ids in tax_hits.items():
        taxon = _majority_consensus(tax_ids, min_fraction, taxon_cache=shared_cache)
        if taxon is not None:
            result[cid] = taxon
    print(f"[taxonomy] Konsens: {len(result)}/{len(tax_hits)} Contigs aufgelöst")
    return result


# ── Bin-Level Taxonomie ───────────────────────────────────────────────────────

def assign_bin_taxonomy(
    labels: np.ndarray,
    contig_ids: List[str],
    tax_hits: Dict[str, List[int]],
    min_fraction: float = 0.5,
    n_total_contigs: Optional[int] = None,
    confidence_min_depth: int = 8,
) -> pd.DataFrame:
    """Majority-LCA pro Bin → DataFrame mit Taxonomie-Zuweisung.

    Spalten: bin_id, bin_label, n_contigs, n_contigs_with_hits, n_hits,
             taxonomy, taxonomy_short, taxonomy_rank, taxonomy_depth, confidence

    confidence = Anteil der Contigs im Bin deren Hit mit der Bin-Taxonomie
                 einen gemeinsamen Vorfahren bei Tiefe >= confidence_min_depth
                 teilt, geteilt durch alle Contigs im Bin.
    """

    bin_map: Dict[int, List[str]] = {}
    for cid, lbl in zip(contig_ids, labels):
        bin_map.setdefault(int(lbl), []).append(cid)

    shared_cache: Dict[int, Taxon] = {}
    rows = []
    for bin_id in sorted(bin_map):
        contigs = bin_map[bin_id]
        all_tax_ids: List[int] = []
        n_with_hits = 0
        for cid in contigs:
            hits = tax_hits.get(cid, [])
            if hits:
                n_with_hits += 1
                all_tax_ids.extend(hits)

        if all_tax_ids:
            lca_taxon = _majority_consensus(all_tax_ids, min_fraction,
                                            taxon_cache=shared_cache)
            if lca_taxon is not None:
                lineage   = lca_taxon.lineage
                tax_full  = "; ".join(t.scientific_name for t in reversed(lineage))
                tax_short = lca_taxon.scientific_name
                depth     = len(lineage)
                rank      = lca_taxon.rank.name
                # lineage[0] = Taxon selbst (tiefste Ebene), lineage[-1] = root
                # Tiefe eines Knotens = len(lineage) - position_in_lineage
                lca_depth_map = {
                    t.tax_id: len(lineage) - i
                    for i, t in enumerate(lineage)
                }

                # Confidence: Anteil Contigs im Bin deren Hit mit der
                # Bin-Taxonomie einen gemeinsamen Vorfahren bei Tiefe
                # >= confidence_min_depth teilt.
                n_matching = 0
                for cid in contigs:
                    for tid in tax_hits.get(cid, []):
                        try:
                            if tid not in shared_cache:
                                shared_cache[tid] = Taxon(int(tid))
                            hit_lin_ids = {t.tax_id for t in shared_cache[tid].lineage}
                            common = lca_depth_map.keys() & hit_lin_ids
                            if common:
                                deepest = max(lca_depth_map[a] for a in common)
                                if deepest >= confidence_min_depth:
                                    n_matching += 1
                                    break
                        except Exception:
                            continue
                confidence = round(n_matching / len(contigs), 3)
            else:
                tax_full = tax_short = "unknown"
                depth = 0
                rank  = "unknown"
                confidence = 0.0
        else:
            tax_full = tax_short = "unknown"
            depth = 0
            rank  = "unknown"
            confidence = 0.0

        rows.append({
            "bin_id":              bin_id,
            "bin_label":           "Noise" if bin_id == -1 else f"Bin {bin_id}",
            "n_contigs":           len(contigs),
            "n_contigs_with_hits": n_with_hits,
            "n_hits":              len(all_tax_ids),
            "taxonomy":            tax_full,
            "taxonomy_short":      tax_short,
            "taxonomy_rank":       rank,
            "taxonomy_depth":      depth,
            "confidence":          confidence,
        })

    return pd.DataFrame(rows).sort_values("bin_id").reset_index(drop=True)

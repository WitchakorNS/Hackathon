"""
Layer 6 — Pathway Annotation (DETERMINISTIC).

Hypergeometric (over-representation) enrichment of the significant metabolites
against the local KEGG-like pathway table — chosen over GSEA-style permutation
because the metabolite counts are tiny (Devflow 4.1). Also emits the hub x pathway
cross-link (Devflow 4.2): a metabolite that is BOTH a stable network hub AND sits
in an enriched pathway.
"""
from __future__ import annotations

from typing import Dict, List

from scipy import stats

ENRICH_ALPHA = 0.10


def enrich(
    significant_chebis: List[str],
    measured_chebis: List[str],
    pathways: List[dict],
    hub_chebis: List[str],
    name_of: Dict[str, str],
) -> tuple[List[dict], Dict[str, bool]]:
    measured = set(measured_chebis)
    hits = set(significant_chebis)
    M = len(measured)          # background size
    n = len(hits)              # drawn (significant)

    results: List[dict] = []
    in_enriched: Dict[str, bool] = {c: False for c in measured_chebis}

    for pw in pathways:
        members = measured & set(pw["chebi_members"])
        K = len(members)       # successes in background
        if K == 0:
            continue
        overlap = members & hits
        k = len(overlap)       # observed successes
        # P(X >= k) hypergeometric survival (upper tail, inclusive)
        pval = float(stats.hypergeom.sf(k - 1, M, K, n)) if k > 0 else 1.0
        enriched = pval < ENRICH_ALPHA and k > 0
        if enriched:
            for c in overlap:
                in_enriched[c] = True
        results.append({
            "pathway_id": pw["pathway_id"],
            "name": pw["name"],
            "k_hits": k,
            "k_members": K,
            "pval": round(pval, 4),
            "enriched": enriched,
            "matched_metabolites": sorted(name_of.get(c, c) for c in overlap),
        })

    results.sort(key=lambda r: r["pval"])

    # Hub x pathway cross-link (explicit inner join on chebi_id).
    hub_set = set(hub_chebis)
    for r in results:
        r["hub_crosslinked"] = any(
            in_enriched.get(c, False) for c in measured
            if name_of.get(c, c) in r["matched_metabolites"] and c in hub_set
        )
    return results, in_enriched

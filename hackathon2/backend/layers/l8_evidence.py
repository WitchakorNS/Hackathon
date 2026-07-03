"""
Layer 8 — Evidence Integration (DETERMINISTIC).

Combines the five evidence layers into one convergence score in [0,1]:
  * differential  : -log10(adj p)
  * effect size   : |Cohen's d|
  * network       : degree centrality
  * pathway       : membership in an enriched pathway (0/1)
  * ml            : SHAP-proxy importance

Each is min-max scaled to [0,1] BEFORE combining (Devflow 6.1 — otherwise an
unbounded layer would dominate), then equal-weighted (0.2 each, Devflow 6.2).

Tier is derived from the score, but capped at "Moderate" when n < LOW_POWER_N
regardless of how high the score is (Devflow Risk Register, Stage 6).
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np

from schema import LOW_POWER_N, MetaboliteResult
from layers.l1_resolution import Reference
from layers.l3_differential import DiffRow
from layers.l4_network import NetworkResult


def _minmax(values: Dict[str, float]) -> Dict[str, float]:
    v = np.array(list(values.values()), dtype=float)
    lo, hi = np.nanmin(v), np.nanmax(v)
    rng = hi - lo
    if rng == 0:
        return {k: 0.0 for k in values}
    return {k: float((x - lo) / rng) for k, x in values.items()}


def _tier(score: float, n_samples: int) -> str:
    if score >= 0.70:
        tier = "High (Tier 1)"
    elif score >= 0.40:
        tier = "Moderate (Tier 2)"
    else:
        tier = "Low (Tier 3)"
    # Low-power cap: cannot claim High on thin data.
    if n_samples < LOW_POWER_N and tier == "High (Tier 1)":
        tier = "Moderate (Tier 2)"
    return tier


def integrate(
    diff: Dict[str, DiffRow],
    network: NetworkResult,
    modules: Dict[str, str],
    shap: Dict[str, float],
    in_enriched: Dict[str, bool],
    n_disease: int,
    reference: Reference,
) -> Dict[str, MetaboliteResult]:
    chebis = list(diff.keys())

    # Raw per-layer evidence signals.
    ev_diff = {c: -np.log10(max(diff[c].pval_adj, 1e-12)) for c in chebis}
    ev_effect = {c: abs(diff[c].effect) for c in chebis}
    ev_net = {c: float(network.degree.get(c, 0)) for c in chebis}
    ev_path = {c: (1.0 if in_enriched.get(c, False) else 0.0) for c in chebis}
    ev_ml = {c: float(shap.get(c, 0.0)) for c in chebis}

    s_diff, s_effect, s_net, s_ml = _minmax(ev_diff), _minmax(ev_effect), _minmax(ev_net), _minmax(ev_ml)

    results: Dict[str, MetaboliteResult] = {}
    for c in chebis:
        score = 0.2 * (s_diff[c] + s_effect[c] + s_net[c] + ev_path[c] + s_ml[c])
        ref = reference.by_chebi(c) or {}
        row = diff[c]
        results[ref.get("hmdb_id", c)] = MetaboliteResult(
            id=ref.get("hmdb_id", c),
            chebi_id=c,
            kegg=ref.get("kegg_id", ""),
            name=ref.get("name", c),
            dir=row.direction,
            fc=round(row.log2fc, 2),
            pval=round(row.pval_adj, 3),
            effect=round(row.effect, 2),
            shap=round(ev_ml[c], 2),
            tier=_tier(score, n_disease),
            score=round(score, 2),
            role=network.role.get(c, "Peripheral"),
            degree=int(network.degree.get(c, 0)),
            module=modules.get(c, "M1"),
            low_power_warning=row.low_power_warning,
        )
    return results

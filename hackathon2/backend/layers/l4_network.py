"""
Layer 4 — Correlation Network (DETERMINISTIC).

Builds a metabolite-metabolite graph from the disease cohort. Pearson AND Spearman
are both computed (Devflow 3.1 — some metabolic relationships aren't linear); the
edge weight uses the stronger-of-the-two by magnitude.

Node role is assigned from degree + edge stability, where "stability" is measured
deterministically by splitting the samples into two fixed halves and checking how
many edges persist (a cheap stand-in for the bootstrap in Devflow 3.2).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

EDGE_THRESHOLD = 0.35   # |r| below this is not an edge


@dataclass
class NetworkResult:
    degree: Dict[str, int]
    role: Dict[str, str]
    stable: Dict[str, bool]
    edges: List[Tuple[str, str, float]]


def _corr_edges(z: pd.DataFrame, chebis: List[str]) -> Dict[Tuple[str, str], float]:
    pear = z.corr(method="pearson")
    spear = z.corr(method="spearman")
    edges: Dict[Tuple[str, str], float] = {}
    for i, a in enumerate(chebis):
        for b in chebis[i + 1:]:
            rp = pear.loc[a, b]
            rs = spear.loc[a, b]
            r = rp if abs(rp) >= abs(rs) else rs
            if np.isfinite(r) and abs(r) >= EDGE_THRESHOLD:
                edges[(a, b)] = float(r)
    return edges


def build(disease_z: pd.DataFrame) -> NetworkResult:
    chebis = list(disease_z.columns)
    edges_full = _corr_edges(disease_z, chebis)

    # Deterministic stability: split rows in half, keep edges present in both.
    n = disease_z.shape[0]
    half = n // 2
    e1 = _corr_edges(disease_z.iloc[:half], chebis)
    e2 = _corr_edges(disease_z.iloc[half:], chebis)
    stable_edges = set(e1) & set(e2)

    degree: Dict[str, int] = {c: 0 for c in chebis}
    stable_degree: Dict[str, int] = {c: 0 for c in chebis}
    for (a, b) in edges_full:
        degree[a] += 1
        degree[b] += 1
        if (a, b) in stable_edges:
            stable_degree[a] += 1
            stable_degree[b] += 1

    role: Dict[str, str] = {}
    stable: Dict[str, bool] = {}
    for c in chebis:
        deg = degree[c]
        is_stable = stable_degree[c] >= max(2, deg // 2)
        stable[c] = is_stable
        if deg >= 5 and is_stable:
            role[c] = "Stable Hub"
        elif deg >= 3:
            role[c] = "Hub"
        else:
            role[c] = "Peripheral"

    edge_list = [(a, b, w) for (a, b), w in edges_full.items()]
    return NetworkResult(degree, role, stable, edge_list)

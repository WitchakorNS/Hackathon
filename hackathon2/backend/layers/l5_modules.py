"""
Layer 5 — Module Discovery (DETERMINISTIC).

The Devflow design calls for Louvain community detection; here we use a
deterministic spectral 2-way partition of the (absolute) correlation graph, which
recovers the same energy-vs-nitrogen split without a random Louvain seed. When the
full pipeline is enabled, swap this for python-louvain with a resolution sweep.

Outputs per-metabolite module labels plus module-level metadata:
  * stability_ari : agreement of the partition across two fixed sample halves.
  * eigenmetabolite: mean disease z-score across the module's members (signed).
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd

from schema import ModuleResult

# Anchor: the cluster containing L-Lactic acid is labelled M1 (energy).
_M1_ANCHOR = "CHEBI:422"       # L-Lactic acid
_MODULE_LABELS = {
    "M1": ("Module M1 (Energy)", "orange"),
    "M2": ("Module M2 (Nitrogen Shunt)", "gray"),
}


def _spectral_bisect(z: pd.DataFrame, chebis: List[str]) -> np.ndarray:
    """Return a {0,1} label per chebi via the Fiedler vector of the affinity graph."""
    corr = z.corr(method="pearson").reindex(index=chebis, columns=chebis).to_numpy()
    corr = np.nan_to_num(corr, nan=0.0)
    W = np.abs(corr)
    np.fill_diagonal(W, 0.0)
    d = W.sum(axis=1)
    d_safe = np.where(d > 0, d, 1.0)
    # Normalized Laplacian L = I - D^-1/2 W D^-1/2
    Dinv = np.diag(1.0 / np.sqrt(d_safe))
    L = np.eye(len(chebis)) - Dinv @ W @ Dinv
    vals, vecs = np.linalg.eigh(L)
    fiedler = vecs[:, 1]           # second-smallest eigenvector
    return (fiedler >= 0).astype(int)


def _ari(a: np.ndarray, b: np.ndarray) -> float:
    """Adjusted Rand Index between two label vectors (small-n, pure numpy)."""
    from itertools import combinations
    n = len(a)
    if n < 2:
        return 1.0
    same_a = np.array([a[i] == a[j] for i, j in combinations(range(n), 2)])
    same_b = np.array([b[i] == b[j] for i, j in combinations(range(n), 2)])
    tp = np.sum(same_a & same_b)
    tn = np.sum(~same_a & ~same_b)
    total = len(same_a)
    agree = (tp + tn) / total
    # rescale [0.5,1] agreement -> [0,1]-ish stability; clamp
    return float(max(0.0, min(1.0, 2 * agree - 1)))


def detect(disease_z: pd.DataFrame, name_of: Dict[str, str]) -> tuple[Dict[str, str], Dict[str, ModuleResult]]:
    chebis = list(disease_z.columns)
    labels = _spectral_bisect(disease_z, chebis)

    # Map {0,1} -> {M1,M2} so the cluster with the anchor becomes M1.
    anchor_idx = chebis.index(_M1_ANCHOR) if _M1_ANCHOR in chebis else 0
    m1_label = labels[anchor_idx]
    assign = {c: ("M1" if labels[i] == m1_label else "M2") for i, c in enumerate(chebis)}

    # Stability via two fixed halves.
    n = disease_z.shape[0]
    half = n // 2
    lab1 = _spectral_bisect(disease_z.iloc[:half], chebis)
    lab2 = _spectral_bisect(disease_z.iloc[half:], chebis)
    ari = _ari(lab1, lab2)

    modules: Dict[str, ModuleResult] = {}
    for key in ("M1", "M2"):
        members = [c for c in chebis if assign[c] == key]
        if not members:
            continue
        eigen = float(disease_z[members].mean(axis=1).mean())  # mean member z in disease
        label, color = _MODULE_LABELS[key]
        modules[key] = ModuleResult(
            key=key,
            label=label,
            members=[name_of.get(c, c) for c in members],
            stability_ari=round(ari, 2),
            eigenmetabolite=round(eigen, 2),
            color=color,
        )
    return assign, modules

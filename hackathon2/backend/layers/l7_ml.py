"""
Layer 7 — Feature Importance (DETERMINISTIC proxy for TreeSHAP).

The full pipeline trains XGBoost + nested CV + Optuna and explains it with TreeSHAP.
Until those deps are enabled we use a deterministic, leakage-free importance proxy:
each metabolite's univariate separability between healthy and disease, measured as
|AUC - 0.5| * 2 (a rank statistic, so it's stable and needs no random model fit).

This keeps the same *shape* of output (a global importance per metabolite in [0,1])
that index.html's SHAP panel consumes, so nothing on the frontend changes when the
real TreeSHAP is swapped in later.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd
from scipy import stats


def _auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """AUC via the Mann-Whitney U statistic (rank-based, deterministic)."""
    n1, n2 = len(pos), len(neg)
    if n1 == 0 or n2 == 0:
        return 0.5
    ranks = stats.rankdata(np.concatenate([pos, neg]))
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2.0
    return float(u1 / (n1 * n2))


def importance(healthy_z: pd.DataFrame, disease_z: pd.DataFrame) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for c in disease_z.columns:
        d = disease_z[c].dropna().to_numpy()
        h = healthy_z[c].dropna().to_numpy()
        auc = _auc(d, h)
        out[c] = round(abs(auc - 0.5) * 2.0, 4)   # 0 = no signal, 1 = perfect
    return out

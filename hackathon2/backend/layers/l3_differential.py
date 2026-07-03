"""
Layer 3 — Differential Analysis (REAL).

Per metabolite:
  1. Shapiro-Wilk normality test decides t-test vs Mann-Whitney (don't hardcode).
  2. log2 fold change  = mean(disease log2) - mean(healthy log2).
  3. Cohen's d effect size (pooled SD).
  4. Benjamini-Hochberg FDR across the metabolites of THIS pipeline only.
  5. low_power_warning flag when n < LOW_POWER_N (Devflow: fatty_liver n=18).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats

from schema import LOW_POWER_N


@dataclass
class DiffRow:
    chebi_id: str
    log2fc: float
    direction: str
    effect: float
    pval_raw: float
    pval_adj: float
    test: str
    low_power_warning: bool


def _cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = len(a), len(b)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    pooled = np.sqrt(((na - 1) * va + (nb - 1) * vb) / max(na + nb - 2, 1))
    if pooled == 0:
        return 0.0
    return float((a.mean() - b.mean()) / pooled)


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    return adj


def differential(healthy_log2: pd.DataFrame, disease_log2: pd.DataFrame) -> Dict[str, DiffRow]:
    chebis = list(healthy_log2.columns)
    n_disease = disease_log2.shape[0]
    low_power = n_disease < LOW_POWER_N

    rows: List[DiffRow] = []
    raw_p: List[float] = []
    for c in chebis:
        h = healthy_log2[c].dropna().to_numpy()
        d = disease_log2[c].dropna().to_numpy()
        if len(h) < 3 or len(d) < 3:
            continue

        # Normality on the disease arm (small-sample safe) -> pick test.
        try:
            normal = stats.shapiro(d).pvalue > 0.05 and stats.shapiro(h).pvalue > 0.05
        except Exception:
            normal = False
        if normal:
            stat, p = stats.ttest_ind(d, h, equal_var=False)
            test = "welch_t"
        else:
            stat, p = stats.mannwhitneyu(d, h, alternative="two-sided")
            test = "mann_whitney"

        log2fc = float(d.mean() - h.mean())
        rows.append(DiffRow(
            chebi_id=c,
            log2fc=log2fc,
            direction="Up" if log2fc >= 0 else "Down",
            effect=_cohens_d(d, h),
            pval_raw=float(p),
            pval_adj=float("nan"),  # filled after BH
            test=test,
            low_power_warning=low_power,
        ))
        raw_p.append(float(p))

    adj = _bh_fdr(np.array(raw_p)) if raw_p else np.array([])
    for row, a in zip(rows, adj):
        row.pval_adj = float(a)

    return {r.chebi_id: r for r in rows}

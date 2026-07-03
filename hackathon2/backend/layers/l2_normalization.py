"""
Layer 2 — Normalization (REAL).  The silent-failure danger zone (Devflow Stage 1).

Pipeline: PQN  ->  log2(x + 1)  ->  z-score vs Healthy baseline.
Guarded by a QC checkpoint on L-Alanine's healthy z-score median (must be ~0).

Why these choices (from the doc):
  * PQN over total-sum: robust to a single spiking metabolite.
  * log2(x+1): +1 pseudocount keeps the NMR floor (0.02) non-negative, matches lit.
  * z-score vs Healthy only: preserves the between-group signal (independent
    per-group z-scoring would force both groups to mean 0 and erase the contrast).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import pandas as pd

from schema import GATE_CHEBI


@dataclass
class NormResult:
    chebis: List[str]
    healthy_log2: pd.DataFrame   # samples x chebi
    disease_log2: pd.DataFrame
    healthy_z: pd.DataFrame
    disease_z: pd.DataFrame
    qc: dict


def _to_wide(long_df: pd.DataFrame, chebis: List[str]) -> pd.DataFrame:
    wide = long_df.pivot_table(
        index="sample_id", columns="chebi_id", values="value", aggfunc="mean"
    )
    return wide.reindex(columns=chebis)


def _pqn(wide: pd.DataFrame, reference_spectrum: pd.Series) -> pd.DataFrame:
    """Probabilistic Quotient Normalization against a fixed reference spectrum."""
    quotients = wide.divide(reference_spectrum, axis=1)
    dilution = quotients.median(axis=1)             # per-sample dilution factor
    dilution = dilution.replace(0, np.nan).fillna(1.0)
    return wide.divide(dilution, axis=0)


def normalize(healthy_long: pd.DataFrame, disease_long: pd.DataFrame) -> NormResult:
    # Align on metabolites measured in BOTH groups (the resolvable overlap).
    common = sorted(set(healthy_long["chebi_id"]) & set(disease_long["chebi_id"]))
    h_wide = _to_wide(healthy_long, common)
    d_wide = _to_wide(disease_long, common)

    # PQN reference spectrum = median healthy sample (computed from Healthy only,
    # then applied to Disease too — same reference for a fair comparison).
    ref_spectrum = h_wide.median(axis=0)
    h_pqn = _pqn(h_wide, ref_spectrum)
    d_pqn = _pqn(d_wide, ref_spectrum)

    # log2(x + 1)
    h_log2 = np.log2(h_pqn + 1.0)
    d_log2 = np.log2(d_pqn + 1.0)

    # z-score using Healthy mean/SD (apply identical transform to Disease)
    mu = h_log2.mean(axis=0)
    sd = h_log2.std(axis=0, ddof=1).replace(0, np.nan)
    h_z = (h_log2 - mu) / sd
    d_z = (d_log2 - mu) / sd

    # --- QC checkpoint: Alanine healthy z-score median must sit near 0 ---
    ala_median = float(h_z[GATE_CHEBI].median()) if GATE_CHEBI in h_z else float("nan")
    # Tolerance is the *normalization-bug* threshold, not a distributional test:
    # a real bug (wrong reference axis, un-centred z, wrong group) sends this to
    # 0.5+ or NaN. A skewed-but-correct metabolite median sits well inside 0.20.
    TOL = 0.20
    qc = {
        "checkpoint_metabolite": GATE_CHEBI,
        "alanine_healthy_zscore_median": round(ala_median, 4),
        "tolerance": TOL,
        "passed": bool(abs(ala_median) < TOL),
    }

    return NormResult(common, h_log2, d_log2, h_z, d_z, qc)

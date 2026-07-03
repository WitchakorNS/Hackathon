"""
Synthetic NMR data source — used *now*, before real files exist.

It fabricates biologically-shaped raw data so the Layer 1-3 statistics run for
real (not hardcoded). Two latent module factors (energy vs nitrogen) are baked
into the covariance so that Layer 4-5 genuinely *recovers* the M1/M2 split
rather than being told about it.

Design points that mirror the Devflow document:
  * L-Alanine is emitted in ALL five sources (gate check).
  * breast_cancer is the largest cohort (n≈699) — the "best case" primary.
  * fatty_liver is a tiny edge case (n=18) with few overlapping metabolites.
  * a couple of junk feature labels are injected to exercise orphan detection.
Everything is seeded => fully reproducible.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
from typing import Dict, List

from .base import DataSource, RawCohort

LN2 = np.log(2.0)

# Per-metabolite generation spec for the PRIMARY contrast (healthy vs breast_cancer).
# label            : raw feature label as it would appear in the source file
# baseline         : healthy geometric-mean concentration (linear NMR units)
# log2fc           : intended disease shift in log2 units (sign = direction)
# module           : latent block driving correlation structure
# loading          : how strongly this metabolite loads on its module factor
_PRIMARY_SPEC: List[dict] = [
    {"label": "lactate",       "baseline": 1.80, "log2fc":  0.95, "module": "M1", "loading": 0.72},
    {"label": "citrate",       "baseline": 1.20, "log2fc":  0.62, "module": "M1", "loading": 0.68},
    {"label": "Malic acid",    "baseline": 0.90, "log2fc":  0.50, "module": "M1", "loading": 0.66},
    {"label": "alanine",       "baseline": 1.05, "log2fc":  0.30, "module": "M1", "loading": 0.62},
    {"label": "pyruvate",      "baseline": 0.70, "log2fc":  0.55, "module": "M1", "loading": 0.66},
    {"label": "D-Glucose",     "baseline": 5.20, "log2fc": -0.30, "module": "M2", "loading": 0.58},
    {"label": "succinate",     "baseline": 0.60, "log2fc": -0.55, "module": "M2", "loading": 0.66},
    {"label": "glutamine",     "baseline": 2.40, "log2fc": -0.80, "module": "M2", "loading": 0.72},
    {"label": "glutamate",     "baseline": 1.10, "log2fc": -0.68, "module": "M2", "loading": 0.68},
]

# Slight anti-correlation between the two module factors sharpens the split so
# the deterministic spectral partition is stable across sample halves.
_FACTOR_ANTICORR = 0.30
_NOISE_SIGMA = 0.30

# Which metabolites each *other* source happens to measure (for overlap realism).
# Alanine appears everywhere (gate check).
_SOURCE_COVERAGE: Dict[str, List[str]] = {
    "healthy_m":     [s["label"] for s in _PRIMARY_SPEC],
    "healthy_f":     [s["label"] for s in _PRIMARY_SPEC],
    "obesity":       ["alanine", "lactate", "glutamine", "citrate", "D-Glucose", "pyruvate"],
    "fatty_liver":   ["alanine", "lactate", "glutamine", "citrate", "succinate", "Malic acid"],
    "breast_cancer": [s["label"] for s in _PRIMARY_SPEC],
}

_SOURCE_N = {
    "healthy_m": 110, "healthy_f": 108, "obesity": 60,
    "fatty_liver": 18, "breast_cancer": 699,
}

# Junk labels to prove orphan/anti-join detection works (won't resolve to CHEBI).
_ORPHANS = ["unknown_peak_7.32ppm", "contaminant_edta"]


class SyntheticDataSource(DataSource):
    def __init__(self, seed: int = 20260703):
        self.seed = seed

    def describe(self) -> str:
        return f"SyntheticDataSource(seed={self.seed}) - fabricated NMR, reproducible"

    def _seed_for(self, src: str) -> int:
        # Stable across processes — Python's built-in hash() is randomized by
        # PYTHONHASHSEED, which would make "reproducible" a lie between restarts.
        digest = hashlib.sha256(f"{self.seed}:{src}".encode()).hexdigest()
        return int(digest[:8], 16)

    def _gen_matrix(self, rng, labels, n, group, disease_frac):
        """Generate a long-form frame for one source with module-correlated noise."""
        spec = {s["label"]: s for s in _PRIMARY_SPEC}
        # Shared latent factors per sample => within-module correlation.
        # fM2 is mildly anti-correlated with fM1 to separate the two modules.
        fM1 = rng.normal(0, 1, n)
        indep = rng.normal(0, 1, n)
        fM2 = -_FACTOR_ANTICORR * fM1 + np.sqrt(1 - _FACTOR_ANTICORR**2) * indep
        is_disease = (rng.random(n) < disease_frac).astype(float) if group == "disease" else np.zeros(n)
        rows = []
        sample_ids = [f"{group[:3]}_{i:04d}" for i in range(n)]
        for lab in labels:
            if lab not in spec:
                continue
            s = spec[lab]
            latent = fM1 if s["module"] == "M1" else fM2
            mu = np.log(s["baseline"]) + s["log2fc"] * LN2 * is_disease
            noise = rng.normal(0, _NOISE_SIGMA, n)
            raw = np.exp(mu + noise + s["loading"] * latent)
            raw = np.clip(raw, 0.02, None)  # NMR floor seen in obesity (Devflow 1.2)
            for sid, v in zip(sample_ids, raw):
                rows.append((lab, sid, float(v), group))
        return rows

    def load(self) -> Dict[str, RawCohort]:
        cohorts: Dict[str, RawCohort] = {}
        for src, labels in _SOURCE_COVERAGE.items():
            rng = np.random.default_rng(self._seed_for(src))
            n = _SOURCE_N[src]
            group = "healthy" if src.startswith("healthy") else "disease"
            # For "healthy" sources everything is baseline; disease sources shift.
            disease_frac = 1.0 if group == "disease" else 0.0
            rows = self._gen_matrix(rng, labels, n, group, disease_frac)
            # Inject an orphan feature into a couple of sources.
            if src in ("obesity", "breast_cancer"):
                orphan = _ORPHANS[0] if src == "obesity" else _ORPHANS[1]
                for sid in sorted({r[1] for r in rows}):
                    rows.append((orphan, sid, float(rng.random()), group))
            frame = pd.DataFrame(rows, columns=["feature_label", "sample_id", "value", "group"])
            cohorts[src] = RawCohort(pipeline_id=src, frame=frame)
        return cohorts

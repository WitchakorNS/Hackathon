"""
Layer 0 — Ingestion interface.

``DataSource`` is the seam that lets us "build now, map data later". Everything
downstream consumes ``RawCohort`` objects; it does not care whether they were
synthesised or read from real NMR exports.

To plug in real data later: implement a new ``DataSource`` (see nmr_files.py)
and point the pipeline config at it. No other file changes.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List

import pandas as pd


@dataclass
class RawCohort:
    """
    One source's raw measurements in tidy long form.

    ``frame`` columns:
        feature_label : str   raw identifier as it appears in the source
                              (name / synonym / kegg id — NOT yet resolved)
        sample_id     : str
        value         : float raw concentration (linear scale)
        group         : str   "healthy" | "disease"
    """
    pipeline_id: str
    frame: pd.DataFrame

    @property
    def n_samples(self) -> int:
        return self.frame["sample_id"].nunique()

    @property
    def feature_labels(self) -> List[str]:
        return sorted(self.frame["feature_label"].unique().tolist())


class DataSource(ABC):
    """A provider of raw cohorts. Implementations: synthetic now, NMR files later."""

    @abstractmethod
    def load(self) -> Dict[str, RawCohort]:
        """Return {pipeline_id: RawCohort} for every source this provider knows."""
        raise NotImplementedError

    @abstractmethod
    def describe(self) -> str:
        """Human-readable one-liner about where the data came from (for logs/UI)."""
        raise NotImplementedError

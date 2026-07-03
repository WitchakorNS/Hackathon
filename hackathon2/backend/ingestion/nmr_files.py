"""
Real NMR file source — the "map data later" hook. STUB.

When you have the actual exports, implement ``load()`` here to read them into the
same ``RawCohort`` shape the synthetic source produces. Nothing downstream changes.

Expected on-disk layout (adjust to your real format):
    data/
      healthy_m.csv      wide: rows = features, cols = samples (+ 'feature' col)
      healthy_f.csv
      obesity.csv
      fatty_liver.csv
      breast_cancer.csv

Typical steps (Devflow Stage 0.2 — pandas wide-to-long melt + anti-join):
    df = pd.read_csv(path)
    long = df.melt(id_vars="feature", var_name="sample_id", value_name="value")
    long["feature_label"] = long["feature"]
    long["group"] = "healthy" if name.startswith("healthy") else "disease"
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict

from .base import DataSource, RawCohort


class NmrFileDataSource(DataSource):
    def __init__(self, data_dir: str | Path):
        self.data_dir = Path(data_dir)

    def describe(self) -> str:
        return f"NmrFileDataSource(data_dir={self.data_dir}) — NOT YET IMPLEMENTED"

    def load(self) -> Dict[str, RawCohort]:
        raise NotImplementedError(
            "Real NMR ingestion is not wired yet. Provide files under "
            f"'{self.data_dir}' and implement NmrFileDataSource.load(). "
            "The rest of the pipeline is ready to consume them unchanged."
        )

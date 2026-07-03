"""Layer 0 — Ingestion. Swappable data sources behind one interface."""
from .base import DataSource, RawCohort
from .synthetic import SyntheticDataSource

__all__ = ["DataSource", "RawCohort", "SyntheticDataSource"]

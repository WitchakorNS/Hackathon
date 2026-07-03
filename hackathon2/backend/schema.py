"""
schema.py — Central data contract (source of truth for column naming).

This is the *first* thing built (Devflow Stage 0.1). Every layer downstream keys
off ``chebi_id`` / ``pipeline_id`` / ``sample_id``, so raw data can come from
synthetic generators today and real NMR files later WITHOUT touching Layer 1-8.

Nothing here depends on where the data came from — that is the whole point.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional

# --- Canonical key columns (the "convention กลาง" from Devflow Stage 0.1) ---
COL_CHEBI = "chebi_id"        # canonical molecule key — every layer joins on this
COL_PIPELINE = "pipeline_id"  # which disease cohort / source pipeline
COL_SAMPLE = "sample_id"      # one biological sample (row) within a source
COL_VALUE = "value"           # measured concentration (raw NMR, linear scale)
COL_GROUP = "group"           # "healthy" | "disease" (for differential contrast)

# --- The five data sources named in the Devflow document ---
SOURCES = ["healthy_m", "healthy_f", "obesity", "fatty_liver", "breast_cancer"]

# Alanine is the "gatekeeper" metabolite — present in ALL five sources.
# Used as the Layer 1 gate check and the Layer 2 QC checkpoint.
GATE_CHEBI = "CHEBI:16977"  # L-Alanine

# n below this => confidence_tier is capped at "Moderate" (Devflow Risk Register).
LOW_POWER_N = 30


@dataclass
class MetaboliteResult:
    """
    Per-metabolite output row. Field names mirror exactly what index.html's
    ``metData`` object consumes, so the frontend can ``fetch()`` this as-is.
    """
    id: str            # HMDB accession (frontend key)
    chebi_id: str      # canonical key used internally
    kegg: str
    name: str
    dir: str           # "Up" | "Down"
    fc: float          # log2 fold change (disease / healthy)
    pval: float        # adjusted p-value (BH-FDR)
    effect: float      # signed effect size (Cohen's d)
    shap: float        # Layer 7 global importance (deterministic proxy)
    tier: str          # "High (Tier 1)" | "Moderate (Tier 2)" | "Low (Tier 3)"
    score: float       # Layer 8 convergence score [0,1]
    role: str          # "Stable Hub" | "Hub" | "Peripheral"
    degree: int        # Layer 4 network degree
    module: str        # "M1" | "M2" (Layer 5)
    low_power_warning: bool = False

    def to_frontend(self) -> dict:
        """Shape used by index.html metData (drops internal-only fields)."""
        d = asdict(self)
        d.pop("chebi_id", None)
        d.pop("low_power_warning", None)
        return d


@dataclass
class ModuleResult:
    key: str                       # "M1" | "M2"
    label: str                     # "Module M1 (Energy)"
    members: List[str]             # metabolite names
    stability_ari: float           # cluster stability [0,1]
    eigenmetabolite: float         # signed module-level shift
    color: str                     # "orange" | "gray" (frontend styling hint)


@dataclass
class PipelineOutput:
    generated_at: str
    primary_cohort: str
    n_samples: Dict[str, int]
    ingestion: dict
    resolution: dict
    qc: dict
    metabolites: Dict[str, MetaboliteResult]
    modules: Dict[str, ModuleResult]
    pathways: List[dict]
    log: List[str] = field(default_factory=list)

    def to_json(self) -> dict:
        return {
            "generated_at": self.generated_at,
            "primary_cohort": self.primary_cohort,
            "n_samples": self.n_samples,
            "ingestion": self.ingestion,
            "resolution": self.resolution,
            "qc": self.qc,
            "metabolites": {k: v.to_frontend() for k, v in self.metabolites.items()},
            "modules": {k: asdict(v) for k, v in self.modules.items()},
            "pathways": self.pathways,
            "log": self.log,
        }

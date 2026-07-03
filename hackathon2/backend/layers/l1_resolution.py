"""
Layer 1 — Resolution (REAL).

Resolves raw ``feature_label`` -> canonical ``chebi_id`` against the local
reference cache. Everything downstream keys off chebi_id, so this is the layer
that MUST be correct — hence the Alanine gate check (Devflow Stage 0.3 / Gate 0):
if L-Alanine does not resolve in all five sources, the whole run is halted.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import pandas as pd

from schema import GATE_CHEBI, SOURCES
from ingestion.base import RawCohort

_REF_PATH = Path(__file__).resolve().parent.parent / "reference" / "metabolite_reference.json"


class Reference:
    """Local canonical lookup (CHEBI cache). Loaded once from JSON."""

    def __init__(self, path: Path = _REF_PATH):
        data = json.loads(path.read_text(encoding="utf-8"))
        self.metabolites = data["metabolites"]
        self.pathways = data["kegg_pathways"]
        self._index: Dict[str, dict] = {}
        for m in self.metabolites:
            keys = {m["name"].lower(), m["kegg_id"].lower(), m["chebi_id"].lower()}
            keys |= {s.lower() for s in m.get("synonyms", [])}
            for k in keys:
                self._index[k] = m

    def resolve_label(self, label: str) -> dict | None:
        return self._index.get(str(label).strip().lower())

    def by_chebi(self, chebi_id: str) -> dict | None:
        for m in self.metabolites:
            if m["chebi_id"] == chebi_id:
                return m
        return None


@dataclass
class ResolutionResult:
    resolved: Dict[str, pd.DataFrame]          # pipeline_id -> long frame w/ chebi_id
    present_in_sources: Dict[str, List[str]]   # chebi_id -> [source, ...]
    orphans: Dict[str, List[str]]              # pipeline_id -> [unresolved labels]
    gate_passed: bool
    report: dict = field(default_factory=dict)


def resolve(cohorts: Dict[str, RawCohort], reference: Reference) -> ResolutionResult:
    resolved: Dict[str, pd.DataFrame] = {}
    present: Dict[str, set] = {}
    orphans: Dict[str, List[str]] = {}

    for src, cohort in cohorts.items():
        df = cohort.frame.copy()
        mapped = df["feature_label"].map(
            lambda lab: (reference.resolve_label(lab) or {}).get("chebi_id")
        )
        df["chebi_id"] = mapped
        # anti-join: labels that did not resolve are orphans
        orphan_labels = sorted(df.loc[df["chebi_id"].isna(), "feature_label"].unique().tolist())
        orphans[src] = orphan_labels
        clean = df.dropna(subset=["chebi_id"]).reset_index(drop=True)
        resolved[src] = clean
        for c in clean["chebi_id"].unique():
            present.setdefault(c, set()).add(src)

    present_lists = {c: sorted(list(s)) for c, s in present.items()}

    # --- GATE CHECK: Alanine must be present in ALL five sources ---
    alanine_sources = set(present_lists.get(GATE_CHEBI, []))
    gate_passed = alanine_sources == set(SOURCES)

    report = {
        "sources": list(cohorts.keys()),
        "n_resolved_features": {s: int(df["chebi_id"].nunique()) for s, df in resolved.items()},
        "orphans": orphans,
        "gate": {
            "metabolite": GATE_CHEBI,
            "required_sources": SOURCES,
            "found_in": sorted(alanine_sources),
            "passed": gate_passed,
        },
    }
    return ResolutionResult(resolved, present_lists, orphans, gate_passed, report)

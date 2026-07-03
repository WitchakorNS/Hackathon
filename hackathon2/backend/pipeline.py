"""
Pipeline orchestrator — wires Layer 0 (ingestion) through Layer 8 (evidence).

Run order here is the RUNTIME order (not the Devflow *development* order). The
primary differential contrast is Healthy (healthy_m + healthy_f) vs breast_cancer,
which is exactly the single metabolite set index.html renders.

    from pipeline import run_pipeline
    out = run_pipeline()          # -> PipelineOutput
    out.to_json()                 # -> dict the API/frontend consumes
"""
from __future__ import annotations

import datetime as _dt
from typing import Dict

import pandas as pd

from schema import PipelineOutput, SOURCES
from ingestion.base import DataSource
from ingestion.synthetic import SyntheticDataSource
from layers.l1_resolution import Reference, resolve
from layers import l2_normalization as l2
from layers import l3_differential as l3
from layers import l4_network as l4
from layers import l5_modules as l5
from layers import l6_pathway as l6
from layers import l7_ml as l7
from layers import l8_evidence as l8

PRIMARY_COHORT = "breast_cancer"
HEALTHY_SOURCES = ["healthy_m", "healthy_f"]


class GateError(RuntimeError):
    """Raised when the Layer 1 Alanine gate check fails — halt, don't propagate."""


def run_pipeline(source: DataSource | None = None) -> PipelineOutput:
    log: list[str] = []
    source = source or SyntheticDataSource()
    reference = Reference()

    # --- Layer 0: Ingestion ---
    cohorts = source.load()
    log.append(f"[L0] Ingestion via {source.describe()}")
    log.append(f"[L0] Sources loaded: {', '.join(cohorts)}  "
               f"(n = {{{', '.join(f'{k}:{v.n_samples}' for k, v in cohorts.items())}}})")

    # --- Layer 1: Resolution + gate check ---
    res = resolve(cohorts, reference)
    log.append(f"[L1] Orphans dropped: "
               f"{sum(len(v) for v in res.orphans.values())} labels "
               f"({res.orphans})")
    if not res.gate_passed:
        raise GateError(
            f"Alanine gate FAILED — resolved only in {res.report['gate']['found_in']}, "
            f"expected all of {SOURCES}"
        )
    log.append(f"[L1] GATE PASSED — Alanine resolves in all {len(SOURCES)} sources")

    # Build the primary contrast: combined healthy vs breast_cancer.
    healthy_long = pd.concat([res.resolved[s] for s in HEALTHY_SOURCES], ignore_index=True)
    disease_long = res.resolved[PRIMARY_COHORT]

    # --- Layer 2: Normalization + QC ---
    norm = l2.normalize(healthy_long, disease_long)
    log.append(f"[L2] PQN + log2(x+1) + z-score on {len(norm.chebis)} metabolites")
    log.append(f"[L2] QC Alanine healthy z-median = "
               f"{norm.qc['alanine_healthy_zscore_median']} "
               f"({'PASS' if norm.qc['passed'] else 'FAIL'})")

    # --- Layer 3: Differential ---
    diff = l3.differential(norm.healthy_log2, norm.disease_log2)
    n_disease = norm.disease_z.shape[0]
    n_sig = sum(1 for r in diff.values() if r.pval_adj < 0.05)
    log.append(f"[L3] Differential ({PRIMARY_COHORT}, n={n_disease}): "
               f"{n_sig}/{len(diff)} significant at FDR<0.05")

    # --- Layer 4: Network ---
    net = l4.build(norm.disease_z)
    log.append(f"[L4] Network: {len(net.edges)} edges, "
               f"{sum(1 for r in net.role.values() if 'Hub' in r)} hubs")

    # --- Layer 5: Modules ---
    name_of = {c: (reference.by_chebi(c) or {}).get("name", c) for c in norm.chebis}
    module_of, modules = l5.detect(norm.disease_z, name_of)
    log.append(f"[L5] Modules: " + ", ".join(
        f"{k}={len(m.members)} (ARI {m.stability_ari})" for k, m in modules.items()))

    # --- Layer 6: Pathway ---
    # "Signature hits" = FDR-significant AND strong effect. Using only FDR<0.05 is
    # degenerate here (every metabolite in this small panel is significant), so the
    # over-representation test needs the stronger contrast to be meaningful.
    sig_chebis = [c for c, r in diff.items() if r.pval_adj < 0.05 and abs(r.effect) >= 0.8]
    hub_chebis = [c for c, r in net.role.items() if "Hub" in r]
    pathways, in_enriched = l6.enrich(sig_chebis, norm.chebis, reference.pathways, hub_chebis, name_of)
    log.append(f"[L6] Pathway: {sum(1 for p in pathways if p['enriched'])} enriched "
               f"of {len(pathways)} tested")

    # --- Layer 7: ML importance ---
    shap = l7.importance(norm.healthy_z, norm.disease_z)
    log.append(f"[L7] SHAP-proxy importance computed for {len(shap)} metabolites")

    # --- Layer 8: Evidence integration ---
    metabolites = l8.integrate(diff, net, module_of, shap, in_enriched, n_disease, reference)
    n_high = sum(1 for m in metabolites.values() if "High" in m.tier)
    log.append(f"[L8] Convergence scored: {n_high} High-tier signatures")

    return PipelineOutput(
        generated_at=_dt.datetime.now().isoformat(timespec="seconds"),
        primary_cohort=PRIMARY_COHORT,
        n_samples={k: v.n_samples for k, v in cohorts.items()},
        ingestion={
            "source": source.describe(),
            "orphans": res.orphans,
            "n_resolved_features": res.report["n_resolved_features"],
        },
        resolution={
            "gate": res.report["gate"],
            "present_in_sources": res.present_in_sources,
        },
        qc=norm.qc,
        metabolites=metabolites,
        modules=modules,
        pathways=pathways,
        log=log,
    )


if __name__ == "__main__":
    import json
    out = run_pipeline()
    print(json.dumps(out.to_json(), indent=2, ensure_ascii=False))

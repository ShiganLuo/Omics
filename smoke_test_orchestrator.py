"""Smoke test for LLM-orchestrator architecture in mode_auto.

Loads an existing advanced h5ad (already has leiden, cell_type, rank_genes_groups,
X_umap from a previous run), then exercises:
  1. _collect_cluster_state — full perception snapshot
  2. _parse_orchestrator_decision — JSON schema validation
  3. _apply_orchestrator_decision — simulated LLM decisions applied to adata

No real LLM calls. We feed pre-canned decision dicts to verify execution path.
"""

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np

sys.path.insert(0, "/home/luosg/Data/genomeStability/workflow/Omics/modules/scanpy/bin")
import scRNAseq as scr


def main():
    adata_path = Path("/home/luosg/Data/genomeStability/output/luancao/scRNAseq/common/5_combine_h5ad/ovaries/ovaries_cellranger_advanced.h5ad")
    print(f"[1/4] Loading {adata_path.name}...")
    adata = ad.read_h5ad(adata_path)
    print(f"      n_obs={adata.n_obs} n_vars={adata.n_vars}")
    print(f"      clusters={sorted(adata.obs['leiden'].unique(), key=lambda x: int(x))}")
    print(f"      cell_types={sorted(adata.obs['cell_type'].unique())}")

    # Synthesize a minimal annotations dict from the existing cell_type column.
    # In real mode_auto, _run_one_annotation_pass produces this; here we fake it
    # by reading existing cell_type + adding shape-matching defaults.
    annotations = {}
    for cl in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
        mask = adata.obs["leiden"] == cl
        ct = adata.obs.loc[mask, "cell_type"].iloc[0]
        annotations[cl] = {
            "cell_type": ct,
            "cell_type_raw": ct,
            "key_markers": [],
            "canonical_markers": [],
            "reasoning": f"smoke-test synthetic reasoning for cluster {cl}",
            "confidence": "high",
            "quality_flag": None,
            "is_subcluster": False,
            "parent_cluster": None,
            "should_merge": False,
            "references": {},
        }

    # Synthesize quality_reports matching the existing cell_type column.
    quality_reports = []
    for cl in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
        mask = adata.obs["leiden"] == cl
        n_cells = int(mask.sum())
        # Make cluster 11 (smallest) flagged as low_quality for the test.
        should_filter = (cl == "11")
        quality_reports.append({
            "cluster": cl,
            "n_cells": n_cells,
            "mean_genes": float(adata.obs.loc[mask, "n_genes_by_counts"].mean()),
            "mean_counts": float(adata.obs.loc[mask, "total_counts"].mean()),
            "pct_mt": float(adata.obs.loc[mask, "pct_counts_mt"].mean()),
            "ai_cell_type": annotations[cl]["cell_type"],
            "ai_confidence": "high",
            "ai_quality_flag": "low_quality" if should_filter else None,
            "flags": ["low_quality"] if should_filter else [],
            "should_filter": should_filter,
            "filter_mode": "cell",
        })
    flagged = [r for r in quality_reports if r["should_filter"]]
    print(f"[2/4] Synthetic flagged clusters: {[r['cluster'] for r in flagged]}")

    # ── Test 1: _collect_cluster_state builds a full perception snapshot ──
    print("\n[3/4] Running _collect_cluster_state...")
    separation_diag = {"separated_types": [], "misannotated": [], "cluster_distances": {}}
    continuity_diag = {"discontinuous_clusters": [], "details": {}, "n_discontinuous": 0}
    state = scr._collect_cluster_state(
        adata=adata,
        annotations=annotations,
        quality_reports=quality_reports,
        separation_diag=separation_diag,
        continuity_diag=continuity_diag,
        tissue_cell_types={},
    )
    print(f"      state['n_cells']={state['n_cells']}")
    print(f"      state['n_clusters']={state['n_clusters']}")
    print(f"      state['clusters'] has {len(state['clusters'])} cluster records")
    sample_cl = sorted(state['clusters'].keys(), key=lambda x: int(x))[0]
    sample = state['clusters'][sample_cl]
    print(f"      sample cluster '{sample_cl}': keys={sorted(sample.keys())}")
    print(f"        cell_type={sample['cell_type']}")
    print(f"        umap_centroid={sample['umap_centroid']}")
    print(f"        umap_distances count={len(sample['umap_distances'])}")
    print(f"        should_filter={sample['should_filter']}")

    # Verify prompt builds
    prompt = scr._build_orchestrator_prompt(state, iteration=1, max_iterations=3)
    print(f"      prompt length: {len(prompt)} chars (head: {prompt[:120]!r})")

    # ── Test 2: _parse_orchestrator_decision validates simulated LLM JSON ──
    print("\n[4/4] Testing _parse_orchestrator_decision + _apply_orchestrator_decision...")

    # Decision A: filter whole cluster 11 + correct cluster 5's cell_type
    simulated_raw = {
        "action": "filter_and_recluster",
        "reasoning": "Cluster 11 is small + low quality; cluster 5 looks mislabeled.",
        "filter_plan": {
            "whole_cluster_removals": ["11"],
            "cell_level_clusters": [],
        },
        "new_resolution": 0.6,
        "annotation_corrections": {
            "5": {"new_cell_type": "Stromal cells", "merge_to_cluster": None,
                  "reasoning": "looks like stromal based on synthetic test"},
        },
    }
    parsed = scr._parse_orchestrator_decision(simulated_raw)
    print(f"      Decision A parsed: action={parsed['action']}")
    print(f"        whole_cluster_removals={parsed['filter_plan']['whole_cluster_removals']}")
    print(f"        new_resolution={parsed['new_resolution']}")
    print(f"        corrections for {list(parsed['annotation_corrections'].keys())}")

    # Apply it
    adata_test = adata.copy()
    adata_test.obs["cell_type"] = adata_test.obs["leiden"].map(
        {cl: ann["cell_type"] for cl, ann in annotations.items()}
    ).astype("category")
    n_cells_before = adata_test.n_obs
    adata_after, summary = scr._apply_orchestrator_decision(
        adata=adata_test,
        decision=parsed,
        annotations=annotations,
        quality_reports=quality_reports,
    )
    print(f"      Applied: action={summary['action']}")
    print(f"      n_cells: {n_cells_before} -> {adata_after.n_obs}")
    print(f"      applied_changes:")
    for c in summary["applied_changes"]:
        print(f"        - {c}")

    # Decision B: invalid action -> fallback to accept
    parsed_bad = scr._parse_orchestrator_decision({"action": "delete_everything"})
    assert parsed_bad["action"] == "accept", "invalid action must fallback to accept"
    print(f"      Decision B (invalid action) correctly fell back to: {parsed_bad['action']}")

    # Decision C: correct_annotations with merge_to_cluster
    simulated_merge = {
        "action": "correct_annotations",
        "reasoning": "merge cluster 6 into 0 (same cell_type, low marker divergence)",
        "annotation_corrections": {
            "6": {"new_cell_type": "X", "merge_to_cluster": "0",
                  "reasoning": "merge test"},
        },
    }
    parsed_c = scr._parse_orchestrator_decision(simulated_merge)
    adata_test2 = adata.copy()
    adata_test2.obs["cell_type"] = adata_test2.obs["leiden"].map(
        {cl: ann["cell_type"] for cl, ann in annotations.items()}
    ).astype("category")
    _, summary_c = scr._apply_orchestrator_decision(
        adata=adata_test2,
        decision=parsed_c,
        annotations=annotations,
        quality_reports=quality_reports,
    )
    print(f"      Decision C merge applied: {summary_c['applied_changes']}")
    print(f"        cluster 6 cell_type is now: {annotations['6']['cell_type']}")
    print(f"        annotations['6'].get('merged_to'): {annotations['6'].get('merged_to')}")

    print("\n✅ ALL SMOKE TESTS PASSED")


if __name__ == "__main__":
    main()

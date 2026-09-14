#!/usr/bin/env python3
"""Fix orchestrator prompt: use spatial metrics instead of cluster counting."""
import sys

path = sys.argv[1]
with open(path, "r", encoding="utf-8") as f:
    content = f.read()

changes = 0

# ── Fix 1: accept criteria — use separation_diag metrics, not cluster count ──
old1 = '''1. "accept" — translations AND structure are good enough. Pipeline saves and exits.
   ONLY when ALL of these are true:
   - Most clusters have trust=high or trust=medium
   - NO same cell_type appears in >2 clusters (no over-clustering)
   - NO synonym pairs (e.g. Fibroblasts + Stromal, Endothelial + Lymphatic_endothelial)
   - Any trust=low clusters are genuinely uncertain (Unknown, Unverified_TE) with no alternative'''

new1 = '''1. "accept" — translations AND structure are good enough. Pipeline saves and exits.
   ONLY when ALL of these are true:
   - Most clusters have trust=high or trust=medium
   - separation_diag shows NO over_clustering entries (same-type clusters with
     high Jaccard >0.3 are over-fragmented and need merging)
   - Same-type clusters (if any) are spatially close in UMAP (small max_distance)
     with distinct marker profiles (low Jaccard = legitimate subtypes)
   - Any trust=low clusters are genuinely uncertain (Unknown, Unverified_TE) with no alternative'''

assert old1 in content, "Fix 1: old string not found"
content = content.replace(old1, new1, 1)
changes += 1
print("Fix 1: accept criteria now uses separation_diag metrics")

# ── Fix 2: adjust_resolution criteria — use Jaccard, not cluster count ──
old2 = '''4. "adjust_resolution_and_recluster" — keep all cells, change Leiden
   resolution, re-cluster. PREFERRED when:
   - Same cell_type appears in 3+ clusters (over-fragmentation)
   - separation_diag shows high Jaccard (>0.3) between same-type clusters
   - Many trust=low clusters without clear alternatives
   Lower resolution by ~0.15-0.2 each round. Current resolution is in the state.
   REQUIRES: new_resolution."""

rules_block = """## Rules

- TWO kinds of problems need fixing: annotation problems (wrong label) and
  structural problems (over-fragmentation, synonyms). The LLM must check BOTH.
- For each trust=low cluster: if alternative is set, consider
  annotation_corrections; if not, the cluster may need reclustering.
- STRUCTURAL CHECK (mandatory before accepting):
  1. Count how many distinct cell_types appear. For each, count clusters.
     If any cell_type has 3+ clusters -> prefer adjust_resolution_and_recluster.
  2. Look for synonym pairs (Fibroblasts/Stromal, Endothelial/Lymphatic_endothelial,
     Smooth_muscle/Myocyte, etc.). If found -> correct_annotations to unify.
  3. Check separation_diag.over_clustering and separation_diag.misannotated.
     If any exist -> prefer adjust_resolution_and_recluster or correct_annotations.
- should_drop=true means the cluster is a clear artifact
  (low_quality/te_dominated/unknown + trust=low). To drop it, add
  to whole_cluster_removals. To keep it but clean bad cells, add to
  cell_level_clusters.
- "annotation_corrections" only rewrites the cell_type STRING column.
  It does NOT change leiden cluster IDs. To physically merge clusters,
  use adjust_resolution_and_recluster.
- If the same problem persists across rounds, choose a more aggressive
  action next time (lower resolution, filter instead of accept, etc.).
- NEVER accept on the first round if there are structural issues. The
  first round is for discovery; accept only when the structure is clean."""'''

new2 = '''4. "adjust_resolution_and_recluster" — keep all cells, change Leiden
   resolution, re-cluster. PREFERRED when:
   - separation_diag.over_clustering shows same-type clusters with Jaccard >0.3
     (same markers split across clusters = genuine over-fragmentation)
   - Many trust=low clusters without clear alternatives
   Lower resolution by ~0.15-0.2 each round. Current resolution is in the state.
   REQUIRES: new_resolution."""

rules_block = """## Rules

- TWO kinds of problems need fixing: annotation problems (wrong label) and
  structural problems (over-fragmentation, synonyms). The LLM must check BOTH.
- For each trust=low cluster: if alternative is set, consider
  annotation_corrections; if not, the cluster may need reclustering.
- STRUCTURAL CHECK (mandatory before accepting):
  1. Check separation_diag.over_clustering: these are cell_types whose clusters
     have HIGH Jaccard (>0.3) = same markers split across clusters = must merge
     via adjust_resolution_and_recluster.
  2. Check separation_diag.misannotated: these are same-name clusters with LOW
     Jaccard (<=0.3) = different markers under same label = likely misannotation.
     Use correct_annotations to rename the flagged cluster.
  3. Same-type clusters with LOW Jaccard and SMALL max_distance are legitimate
     biological subtypes (e.g. proliferating vs quiescent) — this is NORMAL,
     do NOT treat as over-fragmentation.
  4. Look for synonym pairs (Fibroblasts/Stromal, Endothelial/Lymphatic_endothelial,
     Smooth_muscle/Myocyte, etc.). If found -> correct_annotations to unify.
- should_drop=true means the cluster is a clear artifact
  (low_quality/te_dominated/unknown + trust=low). To drop it, add
  to whole_cluster_removals. To keep it but clean bad cells, add to
  cell_level_clusters.
- "annotation_corrections" only rewrites the cell_type STRING column.
  It does NOT change leiden cluster IDs. To physically merge clusters,
  use adjust_resolution_and_recluster.
- If the same problem persists across rounds, choose a more aggressive
  action next time (lower resolution, filter instead of accept, etc.).
- NEVER accept on the first round if there are structural issues. The
  first round is for discovery; accept only when the structure is clean."""'''

assert old2 in content, "Fix 2: old string not found"
content = content.replace(old2, new2, 1)
changes += 1
print("Fix 2: adjust_resolution + rules_block now use Jaccard-based logic")

# Write back
with open(path, "w", encoding="utf-8") as f:
    f.write(content)

print(f"Total changes: {changes}")
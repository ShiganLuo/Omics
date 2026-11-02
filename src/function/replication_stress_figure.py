#!/usr/bin/env python
"""Replication-load / fork-protection / damage-response overview figure.

Builds a two-part figure from raw TEcount gene+TE count matrices:

- Panel A: module-score heatmap (per-gene z of log2-CPM, module mean).
  Modules are authoritative gene sets only (Reactome + MSigDB Hallmark,
  human/mouse matched, see replication_stress_{human,mouse}.gmt):
  replication-stress response arm ATR activation -> Fanconi/fork
  protection -> translesion synthesis -> HRR/fork restart; damage arm
  DSB response + p53 pathway; proliferation load (E2F targets) as
  confound control. Z-scores are computed within species.
- Panels B/C: TE reads fraction (%) vs ATR-activation z (B) and vs
  p53-pathway z (C).

Inputs are raw artifacts only (TEcount_name matrices + the raw-file row mask
that separates gene rows from TE rows); all derived values are computed here.
Set RS_QUANT=tpm to score modules on the precomputed gene TPM tables
(function/{species}/gsva/expression_tpm.tsv) instead of CPM; output then
lands in replication_stress/tpm_version/.

Notes
-----
TLSC2 is excluded (RNA-degradation library artifact, long-gene bias), see
thesis figure S1.
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.text
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path("/data/pub/zhousha/Totipotent20251031/output/pluripotent2totipotent/RNAseq")
COUNTS_DIR = ROOT / "counts" / "TEtranscripts"
TPM_DIR = ROOT / "function"
OUT_DIR = ROOT / "function" / "replication_stress"
# Gene-level quantification for module scores: 'cpm' (from raw counts) or
# 'tpm' (precomputed function/{species}/gsva/expression_tpm.tsv). TE reads
# fraction always comes from the counts matrix (TPM tables have no TE rows).
QUANT = os.environ.get("RS_QUANT", "cpm")
if QUANT == "tpm":
    OUT_DIR = OUT_DIR / "tpm_version"
EXCLUDE = {
    "TLSC2",  # RNA-degradation library artifact, see thesis figure S1
    # GSE175615 fibroblasts are total-RNA (Ribo-Zero) libraries, not poly(A)-
    # comparable; thesis fibroblast controls are HDF n=4 (GSE307377) and
    # WT MEF n=3 (GSE266445), both poly(A).
    "hFibroblasts1", "hFibroblasts2", "hFibroblasts3",
    "mFibroblasts1", "mFibroblasts2",
}

# Gene sets live in the shared geneset assets; replication_stress_{human,mouse}.gmt
# holds exactly the seven authoritative modules used here (Reactome + MSigDB
# Hallmark, see file header rows for source ids). Figure name == GMT row name.
GENESET_DIR = Path("/data/pub/zhousha/Totipotent20251031/workflow/Omics/assests/geneset")
GENESET_FILES = {"GRCm39": "replication_stress_mouse.gmt", "GRCh38": "replication_stress_human.gmt"}
# Causal chain: TE pressure -> fork stalling -> (sensing/protection/
# restart = stress response) -> DSB -> p53 output. Proliferation = confound control.
MODULE_ORDER = [
    # replication-stress response arm (sensing -> protection -> tolerance -> restart)
    "ATR activation (Reactome)",
    "Fork protection / Fanconi anemia (Reactome)",
    "Translesion synthesis (Reactome)",
    "Fork restart / HRR (Reactome)",
    # damage arm (lesion readouts)
    "DNA double strand break response (Reactome)",
    "p53 pathway (Hallmark)",
    # confound control
    "Proliferation load (E2F/G2M, Hallmark)",
]
B_MODULE = "ATR activation (Reactome)"
C_MODULE = "p53 pathway (Hallmark)"
PANEL_MODULES = [B_MODULE, C_MODULE]
# Groups expected to increase in this order (fibroblast -> pluripotent ->
# totipotent-like); used only for the response-module ordering check.
ORDERINGS = {"GRCm39": ["mFibroblasts", "mESC", "ciTotiSC"],
             "GRCh38": ["hFibroblasts", "hESC", "hTBLC"]}

# Group resolution: raw.json first, fibroblasts fall back to name prefix.
GROUP_COLORS = {
    "mFibroblasts": "#7f7f7f", "hFibroblasts": "#bdbdbd",
    "mESC": "#1f77b4", "hESC": "#1f77b4", "prEpiSC": "#9ecae1",
    "TLSC": "#2ca02c", "ciTotiSC": "#d62728",
    "ci8CLC": "#ff7f0e", "hTBLC": "#d62728",
}
# Column order: fibroblast -> pluripotent -> intermediate/8C-like -> totipotent-like.
SPECIES_LAYOUT = {
    "GRCm39": {
        "order": ["mFibroblasts", "mESC", "TLSC", "ciTotiSC"],
        "blocks": ["mFibroblasts1", "mFibroblasts2", "mFibroblasts3", "mFibroblasts4",
                   "mFibroblasts5", "mESC", "mESC1", "mESC2", "TLSC1"]
                 + [f"ciTotiSC{i}-{j}" for i in range(1, 5) for j in (1, 2)],
    },
    "GRCh38": {
        "order": ["hFibroblasts", "hESC", "prEpiSC", "ci8CLC", "hTBLC"],
        "blocks": [f"hFibroblasts{i}" for i in range(1, 8)]
                  + ["WIBR3_E8_hESC1", "WIBR3_E8_hESC2", "WIBR3_E8_hESC3",
                     "prEpiSC_r1", "prEpiSC_r2"]
                  + [f"ci8CLC_r{i}" for i in range(1, 5)]
                  + [f"WIBR3_E8_hTBLC{i}" for i in range(1, 6)],
    },
}


def load_counts(species: str) -> pd.DataFrame:
    """Load the gene-symbol TEcount matrix with TE rows kept as-is.

    Parameters
    ----------
    species : str
        ``GRCm39`` or ``GRCh38``.

    Returns
    -------
    pandas.DataFrame
        Raw counts, index = gene symbol or TE id (``name:family:class``),
        columns = samples.
    """
    name_path = COUNTS_DIR / species / f"{species}_TEcount_name.tsv"
    raw_path = COUNTS_DIR / species / f"{species}_TEcount.tsv"
    df = pd.read_csv(name_path, sep="\t", index_col=0)
    mask = pd.read_csv(raw_path, sep="\t", usecols=[0]).iloc[:, 0].astype(str)
    if len(mask) != len(df):
        raise ValueError(f"row-count mismatch: {name_path} vs {raw_path}")
    # The name file collapses TE ids to their first segment; restore the raw id
    # where the raw file marks a TE row (contains ':'), keeping the row order.
    df.index = [
        raw_id if ":" in raw_id else sym
        for raw_id, sym in zip(mask, df.index.astype(str))
    ]
    log.info(f"{species}: {df.shape[0]} rows x {df.shape[1]} samples")
    return df


def sample_groups() -> dict[str, str]:
    """Resolve sample -> group from raw.json, fibroblasts via name prefix."""
    raw = json.loads((ROOT / "raw.json").read_text())
    groups: dict[str, str] = {}
    for grp, samples in raw["Params"]["StringTie"]["sample_groups"].items():
        for s in samples:
            groups[s] = grp
    return groups


def te_reads_fraction(counts: pd.DataFrame) -> pd.Series:
    """Per-sample fraction of counted reads overlapping TEs (percent)."""
    is_te = counts.index.astype(str).str.contains(":")
    return counts.loc[is_te].sum() / counts.sum() * 100.0


def load_modules(species: str) -> dict[str, list[str]]:
    """Read module gene lists from the geneset GMT file.

    Parameters
    ----------
    species : str
        ``GRCm39`` or ``GRCh38``.

    Returns
    -------
    dict of str to list of str
        Module name -> member gene symbols, in file order.
    """
    path = GENESET_DIR / GENESET_FILES[species]
    modules: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        name, _desc, *genes = line.split("\t")
        modules[name] = genes
    if not modules:
        raise ValueError(f"no modules parsed from {path}")
    missing = [m for m in MODULE_ORDER if m not in modules]
    if missing:
        raise ValueError(f"{species}: modules missing from GMT files: {missing}")
    log.info(f"{species}: {len(MODULE_ORDER)} modules from {path.name}")
    return {m: modules[m] for m in MODULE_ORDER}


def gene_z(counts: pd.DataFrame, expr: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per-gene z-scores of log2 expression, within one species.

    Parameters
    ----------
    counts : pandas.DataFrame
        Raw counts, genes x samples. Used directly when ``expr`` is None
        (log2 CPM+1 computed here).
    expr : pandas.DataFrame, optional
        Precomputed gene-level expression table (genes x samples, e.g. TPM);
        log2(x+1) is applied here. Rows with ':' (TE ids) are dropped.

    Returns
    -------
    pandas.DataFrame
        Samples x genes z-scores.
    """
    if expr is None:
        cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
        logexpr = np.log2(cpm + 1)
    else:
        logexpr = np.log2(expr.astype(float) + 1)
    return logexpr.T.apply(lambda s: (s - s.mean()) / s.std(ddof=0), axis=0)


def load_tpm(species: str, samples: list[str]) -> pd.DataFrame:
    """Load the precomputed gene-level TPM table, restricted to ``samples``.

    Parameters
    ----------
    species : str
        ``GRCm39`` or ``GRCh38``.
    samples : list[str]
        Sample names to keep (columns).

    Returns
    -------
    pandas.DataFrame
        Genes x samples TPM matrix.
    """
    path = TPM_DIR / species / "gsva" / "expression_tpm.tsv"
    tpm = pd.read_csv(path, sep="\t", index_col=0)
    missing = [s for s in samples if s not in tpm.columns]
    if missing:
        raise ValueError(f"{species}: samples missing from {path.name}: {missing}")
    return tpm[samples]


def resolve_genes(columns: pd.Index, genes: list[str]) -> list[str]:
    """Map requested symbols to matrix rows (case-insensitive, genes only)."""
    lookup: dict[str, str] = {}
    for sym in columns:
        sym = str(sym)
        if ":" not in sym:
            lookup.setdefault(sym.upper(), sym)
    return [lookup[g.upper()] for g in genes if g.upper() in lookup]


def module_scores(counts: pd.DataFrame, species: str, modules: dict[str, list[str]],
                  expr: pd.DataFrame | None = None) -> pd.DataFrame:
    """Mean per-gene z-score (log2 expression, within-species) per module.

    Parameters
    ----------
    counts : pandas.DataFrame
        Raw counts for one species.
    species : str
        Species tag for logging.
    modules : dict of str to list of str
        Module name -> gene symbols (from :func:`load_modules`).

    Returns
    -------
    pandas.DataFrame
        Modules x samples mean z-scores.
    """
    z = gene_z(counts, expr)
    rows = {}
    for module, genes in modules.items():
        present = resolve_genes(z.columns, genes)
        if len(present) < max(3, 0.6 * len(genes)):
            raise ValueError(f"{species} {module}: only {len(present)}/{len(genes)} genes found")
        rows[module] = z[present].mean(axis=1)
        log.info(f"{species} {module}: {len(present)}/{len(genes)} genes")
    return pd.DataFrame(rows).T


def robustness(counts: pd.DataFrame, te_frac: pd.Series, groups: pd.Series,
               species: str, module: str, genes: list[str],
               ordering: list[str] | None = None, n_half: int = 50,
               expr: pd.DataFrame | None = None) -> dict:
    """Leave-one-out and random half-split stability of the TE%-module link.

    Parameters
    ----------
    counts : pandas.DataFrame
        Raw counts for one species.
    te_frac : pandas.Series
        TE reads fraction per sample.
    groups : pandas.Series
        Sample -> group.
    species, module : str
        Reporting tags.
    genes : list of str
        Module member symbols.
    ordering : list of str or None
        Groups expected to increase in this order; ``None`` skips the check.
    n_half : int
        Number of random half-split resamples (two halves each).

    Returns
    -------
    dict
        Spearman r range over leave-one-out and half-splits, sign flips, and
        group-ordering violation counts.
    """
    z = gene_z(counts, expr)
    present = resolve_genes(z.columns, genes)
    if len(present) < 3:
        raise ValueError(f"{species} {module}: only {len(present)} genes resolved")

    def spear(gs: list[str]) -> float:
        return float(te_frac.corr(z[gs].mean(axis=1), method="spearman"))

    r_base = spear(present)
    loo = [spear([g for g in present if g != drop]) for drop in present]
    rng = np.random.default_rng(0)
    halves: list[float] = []
    resamples = [[g for g in present if g != drop] for drop in present]
    for _ in range(n_half):
        perm = list(rng.permutation(present))
        half_a, half_b = perm[: len(perm) // 2], perm[len(perm) // 2:]
        halves += [spear(half_a), spear(half_b)]
        resamples += [half_a, half_b]
    halves_arr = np.asarray(halves)

    def ordering_ok(gs: list[str]) -> bool:
        means = z[gs].mean(axis=1).groupby(groups).mean()
        vals = np.array([means.get(g, np.nan) for g in ordering], dtype=float)
        return bool(np.all(np.diff(vals) > 0))

    viol = "NA"
    if ordering is not None:
        n_bad = sum(0 if ordering_ok(gs) else 1 for gs in resamples)
        viol = f"{n_bad}/{len(resamples)}"
    row = {
        "species": species, "module": module, "n_genes": len(present),
        "r_base": round(r_base, 3),
        "r_loo_min": round(min(loo), 3), "r_loo_max": round(max(loo), 3),
        "r_half_p2.5": round(float(np.percentile(halves_arr, 2.5)), 3),
        "r_half_p97.5": round(float(np.percentile(halves_arr, 97.5)), 3),
        "loo_sign_flips": sum(1 for r in loo if np.sign(r) != np.sign(r_base)),
        "ordering_violations": viol,
    }
    log.info(f"{species} {module}: r={row['r_base']} "
             f"loo[{row['r_loo_min']},{row['r_loo_max']}] "
             f"half[{row['r_half_p2.5']},{row['r_half_p97.5']}] "
             f"sign_flips={row['loo_sign_flips']} order_viol={viol}")
    return row


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    groups = sample_groups()
    data = {}
    rob_rows = []
    for species, layout in SPECIES_LAYOUT.items():
        counts = load_counts(species)
        counts = counts.loc[:, [c for c in counts.columns if c not in EXCLUDE]]
        mods = load_modules(species)
        expr = load_tpm(species, list(counts.columns)) if QUANT == "tpm" else None
        scores = module_scores(counts, species, mods, expr)
        te_frac = te_reads_fraction(counts)
        grp = pd.Series({s: groups.get(s, s.rstrip("0123456789")) for s in counts.columns})
        cols = [c for c in layout["blocks"] if c in counts.columns]
        missing = set(counts.columns) - set(cols)
        if missing:
            raise ValueError(f"{species}: samples not placed in layout: {sorted(missing)}")
        data[species] = {"scores": scores[cols], "te_frac": te_frac[cols], "group": grp[cols]}
        for g in grp.unique():
            vals = te_frac[grp[grp == g].index]
            log.info(f"{species} {g}: TE reads% = {vals.min():.2f}-{vals.max():.2f}")
        for module in PANEL_MODULES:
            ordering = ORDERINGS[species] if module == B_MODULE else None
            rob_rows.append(robustness(counts, te_frac, grp, species, module,
                                       mods[module], ordering, expr=expr))
    pd.DataFrame(rob_rows).to_csv(OUT_DIR / "robustness.tsv", sep="\t", index=False)

    # Export derived values for transparency.
    out_tbl = []
    for species, d in data.items():
        t = d["scores"].T.join(d["te_frac"].rename("te_reads_pct")).join(d["group"].rename("group"))
        t.insert(0, "species", species)
        out_tbl.append(t.reset_index(names="sample"))
    pd.concat(out_tbl, ignore_index=True).to_csv(OUT_DIR / "module_scores.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.5, 8.0))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.6, 0.5, 1.0], wspace=0.05, hspace=0.32)
    gs_left = gs[:, 0].subgridspec(2, 1, height_ratios=[14, 0.7], hspace=0.0)
    ax_hm = fig.add_subplot(gs_left[0])
    ax_strip = fig.add_subplot(gs_left[1])  # no sharex: shared axis re-draws ticks
    # Column 1 is a spacer; the colorbar hugs the heatmap's right edge.
    pos = ax_hm.get_position()
    ax_cb = fig.add_axes([pos.x1 + 0.006, pos.y0, 0.012, pos.height])
    ax_b = fig.add_subplot(gs[0, 2])
    ax_c = fig.add_subplot(gs[1, 2])

    # ---- Panel A: module heatmap (z-score, one scale for the whole figure).
    hm = pd.concat([data[s]["scores"] for s in SPECIES_LAYOUT], axis=1)
    vmax = 2.0
    im = ax_hm.imshow(np.clip(hm.values, -vmax, vmax), aspect="auto",
                      cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax_hm.set_yticks(range(hm.shape[0]))
    ax_hm.set_yticklabels(hm.index, fontsize=9)
    ax_hm.set_xticks([])  # no sample ticks/labels; the group strip labels replace them
    ax_strip.set_xticks([])
    ax_hm.set_ylabel("module score (z)", fontsize=10)
    for spine in ax_hm.spines.values():
        spine.set_visible(False)
    # Species divider. Labels in axes fraction (imshow inverts the y axis).
    n_mouse = data["GRCm39"]["scores"].shape[1]
    ax_hm.axvline(n_mouse - 0.5, color="black", lw=1.2)
    frac_mouse = n_mouse / hm.shape[1]
    ax_hm.text(frac_mouse / 2, 1.03, "Mus musculus", transform=ax_hm.transAxes,
               ha="center", va="bottom", fontsize=10, fontweight="bold", clip_on=False)
    ax_hm.text((frac_mouse + 1) / 2, 1.03, "Homo sapiens", transform=ax_hm.transAxes,
               ha="center", va="bottom", fontsize=10, fontweight="bold", clip_on=False)
    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("module score (z, per-species)", fontsize=9)

    # Group colour strip under the heatmap.
    ax_strip.set_ylim(0, 1)
    ax_strip.set_yticks([])
    ax_strip.set_xticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)
    groups_all = pd.concat([data[s]["group"] for s in SPECIES_LAYOUT])
    for x, g in enumerate(groups_all):
        ax_strip.add_patch(Rectangle((x - 0.5, 0.15), 1.0, 0.7,
                                     facecolor=GROUP_COLORS.get(g, "#cccccc"), edgecolor="none"))
    # One horizontal label centred under each group block.
    start = 0
    for x in range(1, len(groups_all) + 1):
        if x == len(groups_all) or groups_all.iloc[x] != groups_all.iloc[start]:
            ax_strip.text((start + x - 1) / 2, -0.12, groups_all.iloc[start],
                          ha="center", va="top", fontsize=6.5, clip_on=False)
            start = x
    ax_strip.set_xlim(-0.5, hm.shape[1] - 0.5)

    # ---- Panels B/C: TE reads fraction vs module z.
    def scatter(ax, module: str, ylabel: str) -> None:
        r_text = []
        for species, d in data.items():
            marker = "o" if species == "GRCm39" else "^"
            for g, idx in d["group"].groupby(d["group"]).groups.items():
                ax.scatter(d["te_frac"].loc[idx], d["scores"].loc[module, idx],
                           c=GROUP_COLORS.get(g, "#cccccc"), marker=marker,
                           s=34, edgecolors="white", linewidths=0.5, label=g, zorder=3)
            r = d["te_frac"].corr(d["scores"].loc[module], method="spearman")
            tag = "mouse" if species == "GRCm39" else "human"
            r_text.append(f"{tag} r={r:.2f}")
            log.info(f"{species} TE% vs {module}: Spearman r = {r:.2f}")
        ax.text(0.97, 0.05, "\n".join(r_text), transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="#333333")
        ax.axhline(0, color="#999999", lw=0.8, ls="--", zorder=1)
        ax.set_xlabel("TE reads fraction (%)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)

    scatter(ax_b, B_MODULE, "ATR / replication-stress response (z)")
    scatter(ax_c, C_MODULE, "p53 pathway (Hallmark, z)")
    ax_b.set_title("Replication-stress response tracks TE load", fontsize=10)
    ax_c.set_title("DNA damage response vs TE load", fontsize=10)

    # One legend for both scatter panels; marker colours equal point colours.
    handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor="white", label=g)
               for g in ["mFibroblasts", "hFibroblasts", "mESC", "hESC", "prEpiSC",
                         "TLSC", "ci8CLC", "ciTotiSC", "hTBLC"]]
    handles += [plt.Line2D([], [], marker="o", color="#333333", ls="", label="mouse"),
                plt.Line2D([], [], marker="^", color="#333333", ls="", label="human")]
    ax_b.legend(handles=handles, fontsize=7, ncol=1, frameon=False,
                loc="upper left", bbox_to_anchor=(1.22, 1.02))

    if os.environ.get("RS_DEBUG"):
        fig.canvas.draw()
        bb_hm = ax_hm.get_window_extent()
        print(f"AX_HM x{bb_hm.x0:.0f}-{bb_hm.x1:.0f} y{bb_hm.y0:.0f}-{bb_hm.y1:.0f}")
        for a in fig.findobj():
            try:
                bb = a.get_window_extent()
            except Exception:
                continue
            if bb.y0 < bb_hm.y0 + 20 and bb.y1 > bb_hm.y0 - 20 and bb.width > 0:
                print(f"NEAR {type(a).__name__} x{bb.x0:.0f}-{bb.x1:.0f} y{bb.y0:.0f}-{bb.y1:.0f} "
                      f"vis={a.get_visible()} parent={getattr(a, 'axes', None)}")
    fig.savefig(OUT_DIR / "replication_stress_overview.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "replication_stress_overview.pdf", bbox_inches="tight")
    log.info(f"written: {OUT_DIR / 'replication_stress_overview.png'}")


if __name__ == "__main__":
    main()

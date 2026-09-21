#!/usr/bin/env python3
"""Plot pseudo-bulk DEG + GO results.

Inputs:
  - {analysis_dir}/deg_pseudobulk/deg_ovaries/*.csv
  - {analysis_dir}/deg_pseudobulk/deg_uterus/*.csv
  - {analysis_dir}/go_macaque/*.csv (clusterProfiler results)

Outputs (saved to {analysis_dir}/figures_pseudobulk/):
  - ovaries_DEG_count_heatmap.png
  - uterus_DEG_count_heatmap.png
  - volcano_<tissue>_<cell_type>.png
  - heatmap_top_<tissue>_<cell_type>.png
  - go_dotplot_<cell_type>_<comparison>.png
  - uterus_DEG_overlap.png

Usage:
  python plot_pseudobulk.py \
    --analysis-dir /home/luosg/Data/genomeStability/output/luancao/scRNAseq/analysis/results
"""
import argparse, os, sys, warnings
warnings.filterwarnings("ignore")
os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

# Uterus comparison label map (clusterProfiler output uses 'a_vs_b' format)
UTERUS_COMP_MAP = {
    "zigong-21310-10XSC3_vs_ZIGONG-21224-10XSC3": "Sham_vs_OV-POI",
    "zigong-11238-10XSC3_vs_zigong-1411200-10XSC3": "Intact_vs_OV-only",
    "zigong-1411200-10XSC3_vs_zigong-2200671-10XSC3": "OV-only_vs_TRA",
}
OVARIES_COMP_MAP = {
    "luanchao-11238-10XSC3_vs_luanchao-21310-10XSC3": "Aging_vs_Youth",
}


def load_all_deg(deg_dir, label_prefix):
    rows = []
    if not os.path.isdir(deg_dir):
        return pd.DataFrame()
    for f in os.listdir(deg_dir):
        if f.startswith(label_prefix) and f.endswith("_DEG.csv"):
            df = pd.read_csv(f"{deg_dir}/{f}")
            rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def deg_count_heatmap(deg_df, comp_map, title, output_path):
    if len(deg_df) == 0:
        return
    # Normalize comparison labels - use short description
    deg_df = deg_df.copy()
    deg_df['comp_short'] = deg_df['comparison'].map(
        lambda x: comp_map.get(x, x.replace('_vs_', ' vs '))
    )
    pv = deg_df.pivot_table(
        index='cell_type', columns='comp_short', values='significant',
        aggfunc='sum', fill_value=0
    )
    fig, ax = plt.subplots(figsize=(max(8, len(pv.columns) * 3), max(6, len(pv) * 0.4)))
    sns.heatmap(pv, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'Significant DEG count'})
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def volcano_plot(deg_df, cell_type, comparison, output_path, title=None,
                log2fc_clip=5, nlogp_clip=18):
    """Volcano plot with biologically meaningful log2FC range.

    log2FC is clipped to [-log2fc_clip, +log2fc_clip] for display only.
    """
    sub = deg_df[(deg_df['cell_type'] == cell_type) &
                 (deg_df['comparison'] == comparison)].copy()
    if len(sub) < 3:
        return False

    # Clip log2FC for display (original values preserved in CSV)
    sub['log2FC_plot'] = sub['log2FC'].clip(-log2fc_clip, log2fc_clip)
    sub['nlogp'] = -np.log10(sub['pval'].clip(1e-300))
    sub['nlogp'] = sub['nlogp'].clip(upper=nlogp_clip)
    sub['sig'] = sub['significant']

    fig, ax = plt.subplots(figsize=(8, 6))
    not_sig = sub[~sub['sig']]
    sig_df = sub[sub['sig']]
    up = sig_df[sig_df['log2FC'] > 0]
    dn = sig_df[sig_df['log2FC'] < 0]

    # Only plot NS points that are NOT in the up/down significant region
    # (filter out the bulk of NS in the central area to reduce clutter)
    not_sig_outside = not_sig[(not_sig['log2FC_plot'].abs() > 1) |
                              (not_sig['nlogp'] > -np.log10(0.05))]

    if len(not_sig_outside) > 0:
        ax.scatter(not_sig_outside['log2FC_plot'], not_sig_outside['nlogp'],
                   s=2, c='#cccccc', alpha=0.3, edgecolor='none',
                   label=f'NS ({len(not_sig)})', zorder=1)
    # Up (red) - larger and on top
    if len(up) > 0:
        ax.scatter(up['log2FC_plot'], up['nlogp'],
                   s=25, c='#d62728', alpha=0.9, edgecolor='darkred',
                   linewidth=0.3, label=f'Up ({len(up)})', zorder=4)
    # Down (blue)
    if len(dn) > 0:
        ax.scatter(dn['log2FC_plot'], dn['nlogp'],
                   s=25, c='#1f77b4', alpha=0.9, edgecolor='darkblue',
                   linewidth=0.3, label=f'Down ({len(dn)})', zorder=4)

    # Label top genes: equal number per direction (5 up + 5 down)
    n_labels = 5
    top_up = up.nlargest(n_labels, 'nlogp')
    top_dn = dn.nlargest(n_labels, 'nlogp')
    for _, r in pd.concat([top_up, top_dn]).iterrows():
        ax.annotate(r['gene'], (r['log2FC_plot'], r['nlogp']),
                    fontsize=8, ha='left', va='bottom', color='black',
                    zorder=4)

    # Threshold lines
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.7, alpha=0.6)
    ax.axvline(1, color='black', linestyle='--', linewidth=0.7, alpha=0.6)
    ax.axvline(-1, color='black', linestyle='--', linewidth=0.7, alpha=0.6)

    # Label clipping
    ax.text(-log2fc_clip * 1.05, nlogp_clip * 0.98,
            f'log2FC clipped at ±{log2fc_clip}', ha='left', va='top',
            fontsize=8, color='gray', style='italic')

    ax.set_xlim(-log2fc_clip * 1.1, log2fc_clip * 1.1)
    ax.set_ylim(0, nlogp_clip)
    ax.set_xlabel('log2FC (clipped at ±' + str(log2fc_clip) + ')', fontsize=11)
    ax.set_ylabel('-log10(pvalue)', fontsize=11)
    if title is None:
        title = f'{cell_type}: {comparison.replace("_vs_", " vs ")}'
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right', framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return True


def top_deg_heatmap(adata, deg_df, cell_type, comparison, output_path, top_n=20):
    sub = deg_df[(deg_df['cell_type'] == cell_type) &
                 (deg_df['comparison'] == comparison) &
                 (deg_df['significant'])]
    if len(sub) < 3:
        return False
    sub = sub.sort_values('pval').head(top_n)
    genes = sub['gene'].tolist()

    sub_adata = adata[adata.obs['cell_type'] == cell_type]
    raw = sub_adata.raw.to_adata()
    gene_idx = [i for i, g in enumerate(raw.var_names) if g in genes]
    found_genes = [raw.var_names[i] for i in gene_idx]
    if not gene_idx:
        return False

    rows = []
    for sample in sorted(sub_adata.obs['sample_id'].unique()):
        mask = sub_adata.obs['sample_id'] == sample
        sample_cells = raw[mask]
        if sample_cells.X.shape[0] == 0:
            continue
        means = np.array(sample_cells.X[:, gene_idx].mean(axis=0)).flatten()
        rows.append({'sample': sample, **dict(zip(found_genes, means))})
    expr = pd.DataFrame(rows).set_index('sample')
    if expr.shape[0] < 2 or expr.shape[1] < 2:
        return False

    expr_z = (expr - expr.mean()) / (expr.std() + 1e-10)

    fig, ax = plt.subplots(figsize=(10, max(4, len(expr_z) * 0.5)))
    sns.heatmap(expr_z, cmap='RdBu_r', center=0, annot=False,
                linewidths=0.3, ax=ax, cbar_kws={'label': 'Z-score'})
    ax.set_title(
        f'{cell_type} top {len(found_genes)} DEGs ({comparison.replace("_vs_", " vs ")})',
        fontsize=11, fontweight='bold'
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return True


def go_dotplot(go_file, output_path, title, top_n=12):
    """GO dotplot with size legend and compact labels.

    size  = -log10(pvalue)  (bigger = more significant)
    color = pvalue         (more red = more significant)
    x     = gene ratio (DE genes in term / all DE genes)
    """
    if not os.path.exists(go_file):
        return False
    df = pd.read_csv(go_file)
    if len(df) == 0:
        return False
    df = df.sort_values('pvalue').head(top_n).copy()
    if 'GeneRatio' in df.columns:
        parts = df['GeneRatio'].str.split('/', expand=True)
        df['gene_ratio_num'] = parts[0].astype(float) / parts[1].astype(float)
    else:
        df['gene_ratio_num'] = df['Count'] / 100
    df['nlogp'] = -np.log10(df['pvalue'].clip(1e-300))
    df['Description'] = df['Description'].astype(str).str[:60]

    # Color by -log10(p), bigger = more red
    df['_color'] = df['nlogp']
    vmin, vmax = df['_color'].min(), df['_color'].max()
    if vmin == vmax:
        vmin, vmax = 0, max(1, vmax)

    fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.4)))
    sc_plot = ax.scatter(
        df['gene_ratio_num'], range(len(df)),
        s=df['nlogp'] * 25,
        c=df['_color'],
        cmap='YlOrRd',
        vmin=vmin, vmax=vmax,
        alpha=0.85,
        edgecolor='black', linewidth=0.5,
    )
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['Description'], fontsize=9)
    ax.set_xlabel('Gene Ratio (DE in term / total DE)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.invert_yaxis()
    cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('-log10(pvalue)', fontsize=10)

    # Size legend
    nlogp_legend_vals = [2, 5, 10]
    if df['nlogp'].max() < 10:
        nlogp_legend_vals = [v for v in nlogp_legend_vals if v <= df['nlogp'].max()]
    if nlogp_legend_vals:
        # Place size legend on the right
        x_legend = df['gene_ratio_num'].max() * 1.15
        for i, nlp in enumerate(nlogp_legend_vals):
            ax.scatter([x_legend], [-(i + 0.5)],
                       s=nlp * 25, c='gray', alpha=0.5,
                       edgecolor='black', linewidth=0.5)
            ax.text(x_legend + df['gene_ratio_num'].max() * 0.05,
                    -(i + 0.5), f'-log10(p)={nlp:.0f}',
                    va='center', fontsize=8)
        ax.text(x_legend - df['gene_ratio_num'].max() * 0.02,
                0.5, 'Size:',
                ha='right', va='center', fontsize=9, fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return True


def deg_overlap_bar(deg_df, sets_to_compare, output_path, title):
    """Bar plot showing DEG overlap between comparisons."""
    sig_per_set = {}
    for label, comp in sets_to_compare:
        sub = deg_df[(deg_df['comparison'] == comp) & (deg_df['significant'])]
        sig_per_set[label] = set(sub['gene'])

    if not all(sig_per_set.values()):
        return

    keys = list(sig_per_set.keys())
    n = len(keys)

    overlap_data = {}
    for i, k1 in enumerate(keys):
        s1 = sig_per_set[k1]
        for j, k2 in enumerate(keys):
            if j < i:
                continue
            s2 = sig_per_set[k2]
            if i == j:
                overlap_data[f'Only {k1}'] = len(s1 - sig_per_set[keys[0]] - sig_per_set[keys[-1]])
            else:
                common = s1 & s2
                if n == 3 and (i == 0 and j == 2):
                    overlap_data[f'All three'] = len(common & sig_per_set[keys[1]])
                else:
                    other = [sig_per_set[keys[k]] for k in range(n)
                             if k != i and k != j]
                    only_this_pair = common.copy()
                    for o in other:
                        only_this_pair -= o
                    overlap_data[f'{k1} & {k2}'] = len(only_this_pair)

    # Sort by count
    overlap_data = dict(sorted(overlap_data.items(), key=lambda x: -x[1]))

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = sns.color_palette('tab10', len(overlap_data))
    bars = ax.barh(list(overlap_data.keys()), list(overlap_data.values()),
                   color=colors)
    ax.set_xlabel('Number of significant DEGs', fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold')
    for bar, val in zip(bars, overlap_data.values()):
        ax.text(val + 20, bar.get_y() + bar.get_height() / 2,
                str(val), va='center', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-dir', required=True,
                        help='Analysis results directory')
    parser.add_argument('--ovaries-h5ad', required=True)
    parser.add_argument('--uterus-h5ad', required=True)
    parser.add_argument('--top-cell-types', nargs='+', default=None,
                        help='Cell types to plot volcano/heatmap for')
    args = parser.parse_args()

    deg_dir = f"{args.analysis_dir}/deg_pseudobulk"
    go_dir = f"{args.analysis_dir}/go_macaque"
    fig_dir = f"{args.analysis_dir}/figures_pseudobulk"
    os.makedirs(fig_dir, exist_ok=True)

    # Load DEG results
    ov_deg = load_all_deg(f"{deg_dir}/deg_ovaries", "ovaries_")
    ut_deg = load_all_deg(f"{deg_dir}/deg_uterus", "uterus_")
    print(f"Ovaries DEG rows: {len(ov_deg)}")
    print(f"Uterus DEG rows: {len(ut_deg)}")

    # Save combined
    all_deg = pd.concat([ov_deg, ut_deg], ignore_index=True)
    all_deg.to_csv(f"{deg_dir}/all_pseudobulk_DEG.csv", index=False)
    print(f"Combined DEG: {len(all_deg)} rows -> {deg_dir}/all_pseudobulk_DEG.csv")

    # Load h5ad for expression heatmap
    import scanpy as sc
    print("Loading h5ad...")
    adata_ov = sc.read_h5ad(args.ovaries_h5ad)
    adata_ut = sc.read_h5ad(args.uterus_h5ad)

    # 2. Ovaries
    short_label = OVARIES_COMP_MAP.get('luanchao-11238-10XSC3_vs_luanchao-21310-10XSC3', 'Aging_vs_Youth')
    fig, ax = plt.subplots(figsize=(8, 6))
    deg_count_heatmap(ov_deg, OVARIES_COMP_MAP,
                      f'Ovaries: Significant DEG count per cell type ({short_label.replace("_", " ")})',
                      f"{fig_dir}/ovaries_DEG_count_heatmap.png")
    plt.close()  # prevent second plot from overlaying
    deg_count_heatmap(ut_deg, UTERUS_COMP_MAP,
                      'Uterus: Significant DEG count per cell type (pseudo-bulk edgeR)',
                      f"{fig_dir}/uterus_DEG_count_heatmap.png")
    plt.close()

    # 2. Volcano + heatmap for top cell types
    print("\n2. Volcano + heatmap for top cell types...")

    # Ovaries top cell types
    if args.top_cell_types:
        ov_top = args.top_cell_types
    else:
        ov_sig_count = ov_deg[ov_deg['significant']].groupby('cell_type').size()
        ov_top = ov_sig_count.sort_values(ascending=False).head(3).index.tolist()
    print(f"  Ovaries top: {ov_top}")
    for ct in ov_top:
        comp = 'luanchao-11238-10XSC3_vs_luanchao-21310-10XSC3'
        ct_safe = ct.replace(' ', '_')
        # Short title using comp map
        comp_short = OVARIES_COMP_MAP.get(comp, comp).replace('_', ' ')
        title = f'{ct} ({comp_short})'
        ok = volcano_plot(ov_deg, ct, comp,
                          f"{fig_dir}/volcano_ovaries_{ct_safe}.png",
                          title=title)
        if ok:
            print(f"    volcano: {ct}")
        ok = top_deg_heatmap(adata_ov, ov_deg, ct, comp,
                             f"{fig_dir}/heatmap_top_ovaries_{ct_safe}.png")
        if ok:
            print(f"    heatmap: {ct}")

    # Uterus top cell types (per comparison)
    ut_sig_count = ut_deg[ut_deg['significant']].groupby(['cell_type', 'comparison']).size()
    ut_top = ut_sig_count.reset_index(name='n').sort_values('n', ascending=False).head(6)
    print(f"  Uterus top:")
    for _, row in ut_top.iterrows():
        ct = row['cell_type']
        comp = row['comparison']
        ct_safe = ct.replace(' ', '_')
        comp_short = UTERUS_COMP_MAP.get(comp, comp).replace('_', ' ')
        title = f'{ct} ({comp_short})'
        ok = volcano_plot(ut_deg, ct, comp,
                          f"{fig_dir}/volcano_uterus_{ct_safe}_{comp_short.replace(' vs ', '-')}.png",
                          title=title)
        if ok:
            print(f"    volcano: {ct} | {comp_short}")
        ok = top_deg_heatmap(adata_ut, ut_deg, ct, comp,
                             f"{fig_dir}/heatmap_top_uterus_{ct_safe}_{comp_short.replace(' vs ', '-')}.png")
        if ok:
            print(f"    heatmap: {ct} | {comp_short}")

    # 3. GO dotplots (top results)
    print("\n3. GO dotplots...")
    if os.path.isdir(go_dir):
        for go_file in sorted(os.listdir(go_dir)):
            if not go_file.endswith('.csv'):
                continue
            if '_go_up.csv' not in go_file and '_go_down.csv' not in go_file:
                continue
            # Skip KEGG (mostly empty for macaque)
            if '_kegg_' in go_file:
                continue
            # Build short title
            base = go_file.replace('.csv', '')
            if 'ovaries_' in base:
                tissue = 'Ovaries'
                rest = base.replace('ovaries_', '').replace('luanchao-21310-10XSC3_vs_luanchao-11238-10XSC3_', '')
            elif 'uterus_' in base:
                tissue = 'Uterus'
                # Extract comparison and cell type
                for comp_key, comp_short in UTERUS_COMP_MAP.items():
                    if comp_key in base:
                        rest = base.replace('uterus_', '').replace(comp_key + '_', '')
                        break
                else:
                    rest = base.replace('uterus_', '')
            else:
                tissue = ''
                rest = base
            direction = 'UP' if base.endswith('_go_up') else 'DOWN'
            ct_clean = rest.replace('_go_up', '').replace('_go_down', '').replace('_', ' ')
            title = f'{tissue} {ct_clean} GO {direction}'
            out = f"{fig_dir}/go_dotplot_{base}.png"
            ok = go_dotplot(f"{go_dir}/{go_file}", out, title, top_n=12)
            if ok:
                print(f"    GO: {base}")

    # 4. Uterus DEG overlap
    print("\n4. Uterus DEG overlap...")
    deg_overlap_bar(
        ut_deg,
        [
            ('Sham vs OV-POI', 'zigong-21310-10XSC3_vs_ZIGONG-21224-10XSC3'),
            ('Intact vs OV-only', 'zigong-11238-10XSC3_vs_zigong-1411200-10XSC3'),
            ('OV-only vs TRA', 'zigong-1411200-10XSC3_vs_zigong-2200671-10XSC3'),
        ],
        f"{fig_dir}/uterus_DEG_overlap.png",
        'Uterus: DEG overlap across 3 comparisons'
    )

    # Summary
    print(f"\n=== Generated figures in {fig_dir}/ ===")
    for f in sorted(os.listdir(fig_dir)):
        if f.endswith('.png'):
            size = os.path.getsize(f"{fig_dir}/{f}")
            print(f"  {f}: {size/1024:.0f}KB")

    print("\nDone.")


if __name__ == "__main__":
    main()
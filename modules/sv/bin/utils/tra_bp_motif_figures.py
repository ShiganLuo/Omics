#!/usr/bin/env python3
"""Figures for tra_bp_motif outputs (stage 4-7, standalone from the TSVs).

tra_bp_motif_figures.py --indir OUT --outdir OUT/figs \
    [--telocen-bed Mmul_T2T.telocen_regions.bed --fai genome.fa.fai]
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'figure.dpi': 150, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
})
C_BP, C_MATCH, C_UNIF = '#4C72B0', '#C44E52', '#999999'
C_YOUTH, C_AGED, C_SHARED = '#4C72B0', '#C44E52', '#55A868'
C_ENR, C_DEP, C_NS = '#C44E52', '#4C72B0', '#BBBBBB'


def load(indir):
    d = {}
    d['comp'] = pd.read_csv(os.path.join(indir, 'stage4_composition.tsv'), sep='\t')
    d['motif'] = pd.read_csv(os.path.join(indir, 'stage5_motif_enrichment.tsv'), sep='\t')
    d['kmer'] = pd.read_csv(os.path.join(indir, 'stage6_kmer.tsv'), sep='\t')
    d['mh'] = pd.read_csv(os.path.join(indir, 'stage7_junction_mh.tsv'), sep='\t')
    d['bps'] = pd.read_csv(os.path.join(indir, 'windows_features.tsv'), sep='\t')
    bp_path = os.path.join(indir, 'breakpoints.tsv')
    d['breakpoints'] = pd.read_csv(bp_path, sep='\t') if os.path.exists(bp_path) else None
    return d


def fig5_motif_enrichment(motif, outdir):
    """log2 OR (matched background) for structural motifs and TF PWMs."""
    struct = motif[motif.kind == 'structural'].copy()
    tf = motif[motif.kind == 'TF_PWM'].copy()

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    fig.subplots_adjust(wspace=0.9)
    for ax, df, title in ((axes[0], struct, 'structural motifs'),
                          (axes[1], tf, 'TF motifs (HOMER PWM)')):
        df = df.sort_values('or_bg_matched')
        y = np.arange(len(df))
        orr = df.or_bg_matched.values
        lo = df.ci_low_matched.values
        hi = df.ci_high_matched.values
        fdr = df.fdr_bg_matched.values
        colors = [C_ENR if (o > 1 and f < 0.05) else
                  C_DEP if (o < 1 and f < 0.05) else C_NS
                  for o, f in zip(orr, fdr)]
        ax.axvline(0, color='#666666', lw=0.8, ls='--')
        ax.scatter(np.log2(np.clip(orr, 1e-3, None)), y, c=colors, s=22, zorder=3)
        for yi, o, l, h in zip(y, orr, lo, hi):
            ax.plot([np.log2(max(l, 1e-3)), np.log2(max(h, 1e-3))],
                    [yi, yi], color=colors[yi], lw=1.2, alpha=0.7)
        ax.set_yticks(y)
        # keep the motif identity only: "Name(Family)/source/Homer" -> "Name (Family)"
        labels = [m.split('/')[0].replace('(', ' (') for m in df.motif]
        ax.set_yticklabels(labels, fontsize=7.5)
        ax.set_xlabel('log2 odds ratio vs matched background')
        ax.set_title(title)
    fig.suptitle('fig5 breakpoint motif enrichment (501bp windows; FDR<0.05 colored)',
                 fontsize=10, y=1.02)
    fig.savefig(os.path.join(outdir, 'fig5_motif_enrichment.png'))
    plt.close(fig)


def fig6_composition_kmer(comp, kmer, outdir):
    """Left: composition features vs backgrounds. Right: top enriched k-mers."""
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4),
                             gridspec_kw={'width_ratios': [1, 1.3]})

    ax = axes[0]
    feats = comp.feature.tolist()
    x = np.arange(len(feats))
    w = 0.27
    ax.bar(x - w, comp.mean_bp, w, label='breakpoint', color=C_BP)
    ax.bar(x, comp.mean_bg_matched, w, label='bg matched', color=C_MATCH)
    ax.bar(x + w, comp.mean_bg_uniform, w, label='bg uniform', color=C_UNIF)
    ax.set_xticks(x)
    ax.set_xticklabels(feats, rotation=35, ha='right', fontsize=7)
    ax.set_title('window composition')
    ax.legend(fontsize=7, frameon=False)

    ax = axes[1]
    top = kmer.sort_values(['fdr_bg_matched', 'or_bg_matched'],
                           ascending=[True, False]).head(15)
    top = top.iloc[::-1]
    y = np.arange(len(top))
    colors = [C_ENR if f < 0.05 else C_NS for f in top.fdr_bg_matched]
    ax.barh(y, np.log2(top.or_bg_matched.values), color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(top.kmer.values, fontsize=7, family='monospace',
                       rotation=0)
    ax.axvline(0, color='#666666', lw=0.8, ls='--')
    ax.set_xlabel('log2 OR vs matched background')
    ax.set_title('top enriched canonical 6-mers (ranked by FDR)')
    fig.suptitle('fig6 composition and k-mer background', fontsize=10, y=1.02)
    fig.savefig(os.path.join(outdir, 'fig6_composition_kmer.png'))
    plt.close(fig)


def fig7_mechanism(mh, outdir):
    """Left: mechanism x group stacked bars. Right: junction MH vs null."""
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2))

    ax = axes[0]
    order = ['TE-NAHR', 'microhomology-altEJ', 'short-MH', 'blunt-NHEJ-like']
    ct = pd.crosstab(mh.mechanism, mh.grp).reindex(order).fillna(0)
    bottom = np.zeros(len(ct))
    for grp, color in (('youth', C_YOUTH), ('aged', C_AGED), ('shared', C_SHARED)):
        if grp in ct.columns:
            ax.bar(ct.index, ct[grp].values, bottom=bottom, label=grp, color=color)
            bottom += ct[grp].values
    ax.set_xticks(range(len(ct.index)))
    ax.set_xticklabels(ct.index, rotation=20, ha='right', fontsize=7)
    ax.set_ylabel('events')
    ax.set_title('mechanism call x group')
    ax.legend(fontsize=7, frameon=False)

    ax = axes[1]
    bins = np.arange(-0.5, 11.5, 1)
    ax.hist(mh.mh_called.dropna(), bins=bins, alpha=0.75, color=C_BP,
            label=f'called join (n={mh.mh_called.notna().sum()})')
    ax.hist(mh.null_mh_max.dropna(), bins=bins, alpha=0.6, color=C_UNIF,
            label=f'null pairs (n={mh.null_mh_max.notna().sum()})')
    ax.set_xlabel('junction microhomology (bp)')
    ax.set_ylabel('events')
    ax.set_title('reference-flank microhomology')
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle('fig7 mechanism classification', fontsize=10, y=1.02)
    fig.savefig(os.path.join(outdir, 'fig7_mechanism.png'))
    plt.close(fig)


def fig8_context(bps_feats, breakpoints, telocen_bed, fai, outdir):
    """Left: structure motif presence vs backgrounds. Right: telomere/cen context."""
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2))

    ax = axes[0]
    keys = [('g4', 'G4'), ('zdna', 'Z-DNA'), ('ir_hits', 'inverted rep'),
            ('dr_hits', 'direct rep'), ('ag_mirror', 'AG mirror'),
            ('telo_any', 'telomere rep'), ('polyat', 'poly-A/T'),
            ('l1en_relax', 'L1 EN')]
    labels = [k[1] for k in keys]
    rates = []
    for label in ('breakpoint', 'bg_matched', 'bg_uniform'):
        sub = bps_feats[bps_feats.label == label]
        rates.append([(sub[k[0]] > 0).mean() for k in keys])
    x = np.arange(len(keys))
    w = 0.27
    ax.bar(x - w, rates[0], w, label='breakpoint', color=C_BP)
    ax.bar(x, rates[1], w, label='bg matched', color=C_MATCH)
    ax.bar(x + w, rates[2], w, label='bg uniform', color=C_UNIF)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha='right', fontsize=7)
    ax.set_ylabel('fraction of windows')
    ax.set_title('structural motif presence')
    ax.legend(fontsize=7, frameon=False)

    ax = axes[1]
    if breakpoints is not None:
        tags = ['telomere_core', 'subtelomere', 'cen_core', 'cen_flank', 'none']
        obs = breakpoints.telocen_tag.value_counts().reindex(tags).fillna(0)
        # genome expectation from region lengths
        genome = 0
        if fai and telocen_bed and os.path.exists(fai) and os.path.exists(telocen_bed):
            with open(fai) as fh:
                genome = sum(int(line.split('\t')[1]) for line in fh)
            region = {t: 0 for t in tags}
            with open(telocen_bed) as fh:
                for line in fh:
                    f = line.rstrip('\n').split('\t')
                    if len(f) >= 4 and f[3] in region:
                        region[f[3]] += int(f[2]) - int(f[1])
            exp = pd.Series({t: region[t] / genome * len(breakpoints)
                             for t in tags})
        else:
            exp = pd.Series({t: np.nan for t in tags})
        x = np.arange(len(tags))
        ax.bar(x, obs.values, 0.42, label='observed', color=C_BP)
        ax.bar(x + 0.42, exp.values, 0.42, label='genome expectation', color=C_UNIF)
        ax.set_xticks(x + 0.21)
        ax.set_xticklabels(tags, rotation=25, ha='right', fontsize=7)
        ax.set_ylabel('breakpoints')
        ax.set_title('telomere/centromere context')
        ax.legend(fontsize=7, frameon=False)
    else:
        ax.text(0.5, 0.5, 'breakpoints.tsv not found', ha='center')
    fig.suptitle('fig8 breakpoint context', fontsize=10, y=1.02)
    fig.savefig(os.path.join(outdir, 'fig8_context.png'))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--indir', required=True, help='tra_bp_motif output dir')
    ap.add_argument('-o', '--outdir', default='', help='figure dir (default <indir>/figs)')
    ap.add_argument('-c', '--telocen-bed', default='')
    ap.add_argument('-x', '--fai', default='')
    args = ap.parse_args()
    outdir = args.outdir or os.path.join(args.indir, 'figs')
    os.makedirs(outdir, exist_ok=True)
    d = load(args.indir)
    fig5_motif_enrichment(d['motif'], outdir)
    fig6_composition_kmer(d['comp'], d['kmer'], outdir)
    fig7_mechanism(d['mh'], outdir)
    fig8_context(d['bps'], d['breakpoints'], args.telocen_bed, args.fai, outdir)
    print('[done] figures in', outdir)


if __name__ == '__main__':
    main()

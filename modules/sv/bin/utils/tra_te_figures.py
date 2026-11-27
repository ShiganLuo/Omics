#!/usr/bin/env python3
"""Figures for tra_te_analysis outputs (stage4, standalone from the 5 TSVs).

tra_te_figures.py --indir OUT --outdir OUT/figs
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

plt.rcParams.update({
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'figure.dpi': 150, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
})
C_YOUTH, C_AGED, C_SHARED = '#4C72B0', '#C44E52', '#55A868'
GRAY = '#999999'


def load(indir):
    d = {}
    d['ev'] = pd.read_csv(os.path.join(indir, 'tra_events_annotated.tsv'), sep='\t')
    d['s1'] = pd.read_csv(os.path.join(indir, 'stage1_enrichment.tsv'),
                          sep='\t', index_col=0, names=['k', 'v'], skiprows=1)
    d['s2'] = pd.read_csv(os.path.join(indir, 'stage2_family_pairs.tsv'), sep='\t')
    d['s3'] = pd.read_csv(os.path.join(indir, 'stage3_homology.tsv'), sep='\t')
    return d


def parse_key(k):
    m = re.findall(r"\('([^']+)', (\d+)\)", str(k))
    return (m[0][0], int(m[0][1])), (m[1][0], int(m[1][1]))


def fig1_enrichment(s1, outdir, n_ev=None, n_rep=2000, seed=42):
    """Observed TE-hit rates vs background / permutation null.

    The right panel's null histogram is a binomial reconstruction
    parameterised by the script's permutation outputs (n events,
    p = perm_single^2) — the pipeline stores only summary values.
    """
    bg = float(s1.loc['bg', 'v'])
    perm_single = float(s1.loc['perm_single', 'v'])
    perm_both = float(s1.loc['perm_both', 'v'])
    obs_bp1 = float(s1.loc['obs_bp1', 'v'])
    obs_bp2 = float(s1.loc['obs_bp2', 'v'])
    obs_both = float(s1.loc['obs_both', 'v'])
    n_ev = int(n_ev or s1.loc['n', 'v'])

    rng = np.random.default_rng(seed)
    p_single = perm_single
    null_both = rng.binomial(n=n_ev, p=p_single ** 2, size=n_rep) / n_ev

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    ax = axes[0]
    bars = [('obs bp1', obs_bp1), ('obs bp2', obs_bp2),
            ('random single', perm_single), ('genome bg', bg)]
    x = np.arange(len(bars))
    ax.bar(x, [b[1] for b in bars],
           color=[C_YOUTH, C_AGED, GRAY, 'lightgray'], width=.62)
    for i, (_, v) in enumerate(bars):
        ax.text(i, v + .01, f'{v:.2f}', ha='center', fontsize=8)
    ax.set_xticks(x, [b[0] for b in bars], rotation=20, ha='right')
    ax.set_ylabel('fraction of breakpoints in TE')
    ax.set_title('Single breakpoint in TE')
    ax.set_ylim(0, max(obs_bp1, obs_bp2) * 1.25)

    ax = axes[1]
    ax.hist(null_both, bins=30, color='lightgray', edgecolor='white')
    ax.axvline(perm_both, color=GRAY, ls='--', lw=1,
               label=f'null mean {perm_both:.2f}')
    ax.axvline(obs_both, color=C_AGED, lw=1.6,
               label=f'observed {obs_both:.2f}')
    p = float(s1.loc['p_binom_both', 'v'])
    ax.set_title(f'Both breakpoints in TE (binom p={p:.4f})')
    ax.set_xlabel('fraction of events with both breakpoints in TE')
    ax.set_ylabel('permutation count')
    ax.legend(frameon=False, fontsize=8)

    fig.savefig(os.path.join(outdir, 'fig1_stage1_enrichment.png'))
    plt.close(fig)


def fig2_pairs(s2, outdir):
    """Family pair obs vs shuffle expectation (only pairs with exp>=1,
    per method doc 6.4: exp<1 ratios are unstable, reference only)."""
    d = s2[s2['exp_under_shuffle'] >= 1].sort_values(
        'obs', ascending=False).head(10).copy()
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.barh(y + .18, d['obs'], height=.34, color=C_YOUTH, label='observed')
    ax.barh(y - .18, d['exp_under_shuffle'], height=.34, color=GRAY,
            label='shuffle expectation')
    ax.set_yticks(y, d['pair'])
    ax.axvline(1, color='lightgray', lw=.8)
    for yi, (o, e) in zip(y, zip(d['obs'], d['exp_under_shuffle'])):
        if e > 0:
            r = o / e
            tag, tc = '', 'k'
            if e < 1:
                tag, tc = f'{r:.2f} (exp<1, ref only)', GRAY
            elif r > 1.5:
                tag, tc = f'{r:.2f} ENRICHED', C_AGED
            elif r < 0.7:
                tag, tc = f'{r:.2f} DEPLETED', C_YOUTH
            else:
                tag = f'{r:.2f}'
            ax.text(max(o, e) + .3, yi, tag,
                    va='center', fontsize=7.5, color=tc)
    ax.set_xlabel('both-TE events')
    ax.set_title('TE family pairing (obs vs expectation, exp>=1 only)')
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    fig.savefig(os.path.join(outdir, 'fig2_stage2_family_pairs.png'))
    plt.close(fig)


def fig3_homology(s3, outdir):
    """Tiers by group + ident vs hsp_len scatter."""
    h = s3.copy()
    h['nahr'] = (h.ident >= 90) & (h.hsp_len >= 100)
    order = ['youth', 'aged', 'shared']
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.2),
                             gridspec_kw={'width_ratios': [1, 1.4]})

    ax = axes[0]
    tiers = ['NAHR-grade', '80-90%', 'no HSP']
    cmap = {'youth': C_YOUTH, 'aged': C_AGED, 'shared': C_SHARED}
    bottom = np.zeros(3)
    for g in order:
        vals = np.array([
            ((h.grp == g) & h.nahr).sum(),
            ((h.grp == g) & ~h.nahr & (h.ident >= 80)).sum(),
            ((h.grp == g) & (h.ident == 0)).sum()])
        ax.bar(np.arange(3), vals, bottom=bottom, color=cmap[g], label=g)
        bottom += vals
    for i, t in enumerate(bottom):
        if t:
            ax.text(i, t + .4, int(t), ha='center', fontsize=8)
    ax.set_xticks(np.arange(3), tiers, fontsize=8)
    ax.set_ylim(0, bottom.max() * 1.3)
    ax.set_ylabel('both-TE events')
    ax.set_title('Homology tiers (n=46)\nNAHR-grade: >=90% & >=100bp')
    ax.legend(frameon=False, fontsize=8, loc='upper left',
              bbox_to_anchor=(0.0, 1.0))

    ax = axes[1]
    for g in order:
        sub = h[(h.grp == g) & (h.ident > 0)]
        ax.scatter(sub.ident, sub.hsp_len, s=26, alpha=.85,
                   color=cmap.get(g, GRAY), label=g, edgecolor='white', lw=.4)
    ax.axvline(90, color='k', ls=':', lw=.9)
    ax.axhline(100, color='k', ls=':', lw=.9)
    ax.text(90.3, ax.get_ylim()[1] * .5, 'NAHR\nthreshold', fontsize=7,
            color='k', ha='left')
    ax.set_xlabel('best HSP identity (%)')
    ax.set_ylabel('best HSP length (bp)')
    ax.set_yscale('log')
    ax.set_title('TE-copy homology per event')
    ax.legend(frameon=False, fontsize=8)

    fig.savefig(os.path.join(outdir, 'fig3_stage3_homology.png'))
    plt.close(fig)


def fig4_panorama(ev, s3, outdir, fai=None):
    """Circos-style panorama with cubic Bezier links (circlize-like).

    Chromosomes sized by real length (from .fai when given), arc length
    proportional to length; links are cubic Bezier curves whose control
    points pull toward the center -- near arcs stay shallow, far arcs
    bend deep, exactly like circos.genomicLink in circlize.
    """
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches

    def lab(c):
        try:
            n = int(c.split('_')[1].split('.')[0][3:])
        except (ValueError, IndexError):
            return c
        if 406 <= n <= 425:
            return f'chr{n - 405}'
        return 'chrX' if n == 426 else 'chrY'

    length = {}
    if fai and os.path.exists(fai):
        for line in open(fai):
            f_ = line.split('\t')
            length[f_[0]] = int(f_[1])
    chrs = sorted(set(ev.chrom1) | set(ev.chrom2),
                  key=lambda c: -sum(ev.chrom1 == c) - sum(ev.chrom2 == c))
    if not length:
        for c in chrs:
            length[c] = (ev[ev.chrom1 == c].pos1.max()
                         if (ev.chrom1 == c).any() else 1) * 1.1
    total = sum(length[c] for c in chrs)
    n_chr = len(chrs)
    cmap = plt.cm.tab20(np.linspace(0, 1, n_chr))

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    ax.set_aspect('equal')
    ax.axis('off')

    gap = 0.012
    scale = (2 * np.pi - gap * n_chr) / total
    theta_start, edges = 0, {}
    for c in chrs:
        w = length[c] * scale
        edges[c] = (theta_start, theta_start + w)
        theta_start += w + gap

    def to_xy(c, pos, r):
        t = edges[c][0] + (edges[c][1] - edges[c][0]) * pos / length[c]
        return r * np.cos(t), r * np.sin(t)

    R = 1.0
    s3n = s3.copy()
    s3n['key'] = s3n['key'].astype(str)
    ev_key = ev['key'].astype(str)
    for _, r in ev.iterrows():
        t1 = (edges[r.chrom1][0] + (edges[r.chrom1][1] - edges[r.chrom1][0])
              * r.pos1 / length[r.chrom1])
        t2 = (edges[r.chrom2][0] + (edges[r.chrom2][1] - edges[r.chrom2][0])
              * r.pos2 / length[r.chrom2])
        dt = abs((t1 - t2 + np.pi) % (2 * np.pi) - np.pi)
        depth = R * (1 - 0.55 * min(dt / np.pi, 1.0))
        x1, y1 = to_xy(r.chrom1, r.pos1, R)
        x2, y2 = to_xy(r.chrom2, r.pos2, R)
        cm_ = (t1 + t2) / 2 if dt < np.pi else (t1 + t2) / 2 + np.pi
        cx, cy = depth * np.cos(cm_), depth * np.sin(cm_)
        is_nahr = False
        hit = s3n[s3n['key'] == ev_key.loc[r.name] if hasattr(r, 'name')
                  else s3n['key'] == str(r.key)]
        if len(hit):
            row = hit.iloc[0]
            is_nahr = row.ident >= 90 and row.hsp_len >= 100
        if is_nahr:
            col, lw, al, z = '#D62728', 1.4, .95, 3
        else:
            col, lw, al, z = '#7F7F7F', .6, .35, 2
        c1 = (x1 * .55 + cx * .45, y1 * .55 + cy * .45)
        c2 = (x2 * .55 + cx * .45, y2 * .55 + cy * .45)
        path = mpath.Path([(x1, y1), c1, c2, (x2, y2)],
                          [mpath.Path.MOVETO, mpath.Path.CURVE4,
                           mpath.Path.CURVE4, mpath.Path.CURVE4])
        ax.add_patch(mpatches.PathPatch(
            path, fill=False, edgecolor=col, lw=lw, alpha=al, zorder=z))

    lab_angles = np.linspace(0, 2 * np.pi, n_chr, endpoint=False) + np.pi / n_chr
    arc_mid = {c: sum(edges[c]) / 2 for c in chrs}
    order_by_angle = sorted(chrs, key=lambda c: arc_mid[c])
    slot_of = {c: lab_angles[i] for i, c in enumerate(order_by_angle)}

    for i, c in enumerate(chrs):
        t0, t1_ = edges[c]
        tt = np.linspace(t0, t1_, 64)
        ax.fill(R * np.cos(tt), R * np.sin(tt), color=cmap[i], alpha=.4,
                lw=0, zorder=4)
        ax.plot(R * np.cos(tt), R * np.sin(tt), color=cmap[i], lw=2.2,
                alpha=.9, zorder=5, solid_capstyle='butt')
        n_bp = int(((ev.chrom1 == c).sum() + (ev.chrom2 == c).sum()))
        tm = arc_mid[c]
        la = slot_of[c]
        ax.plot([1.02 * np.cos(tm), 1.2 * np.cos(la)],
                [1.02 * np.sin(tm), 1.2 * np.sin(la)],
                color=cmap[i], lw=.7, alpha=.8, zorder=1)
        ang = np.degrees(la)
        rot = ang - 90 if ang <= 180 else ang + 90
        ax.text(1.23 * np.cos(la), 1.23 * np.sin(la),
                f'{lab(c)} ({n_bp})', rotation=rot, rotation_mode='anchor',
                ha='center', va='center', fontsize=8, color=cmap[i])

    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_title('TRA panorama: inter-chromosomal events\n'
                 f'(n={len(ev)} links; red = NAHR-grade homology)',
                 fontsize=11, pad=18)
    fig.savefig(os.path.join(outdir, 'fig4_tra_panorama.png'))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', required=True)
    ap.add_argument('--outdir', default=None)
    ap.add_argument('--fai', default=None,
                    help='reference .fai for real chromosome lengths in fig4')
    args = ap.parse_args()
    outdir = args.outdir or os.path.join(args.indir, 'figs')
    os.makedirs(outdir, exist_ok=True)

    d = load(args.indir)
    fig1_enrichment(d['s1'], outdir)
    fig2_pairs(d['s2'], outdir)
    fig3_homology(d['s3'], outdir)
    fig4_panorama(d['ev'], d['s3'], outdir, fai=args.fai)
    print('[done] figures in', outdir)


if __name__ == '__main__':
    main()

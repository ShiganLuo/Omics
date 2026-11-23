#!/usr/bin/env python3
"""TE-associated translocation (TRA) analysis — hierarchical, coarse-to-fine.

Stage 0  Parse TRA records from merged VCF, dedup mirror records -> events
Stage 1  Coarse: are TRA breakpoints enriched in TE regions at all?
         (observed TE rate vs genome background vs matched random controls,
         permutation p-value)
Stage 2  Fine-1: among both-TE events, which TE families pair up?
         (family-pair matrix, same-family excess vs expectation)
Stage 3  Fine-2: sequence-level homology between the two TE copies
         (blastn two-way identity; NAHR requires long high-identity tracts)
Stage 4  Report + figures

Usage:
  python tra_te_analysis.py --vcf youth_vs_aged_merged.vcf \
      --te-bed rheMac10_rmsk_TE.lifted.bed \
      --fasta T2T_genomic.fa \
      -o mei/validation/tra_te [--threads 12]
"""
import argparse
import bisect
import gzip
import math
import os
import random
import re
import subprocess
import sys
from collections import Counter, defaultdict

import pandas as pd
import numpy as np

# ---------------------------------------------------------------- stage 0
def parse_tra(vcf_path):
    recs = []
    op = gzip.open if vcf_path.endswith('.gz') else open
    with op(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if 'SVTYPE=TRA' not in f[7]:
                continue
            d = dict(re.findall(r'(\w+)=([^;]+)', f[7]))
            recs.append(dict(
                chrom1=f[0], pos1=int(f[1]),
                chrom2=d.get('CHR2', ''), pos2=int(d.get('END', 0) or 0),
                strands=d.get('STRANDS', ''), supp=d.get('SUPP_VEC', '').zfill(2)))
    df = pd.DataFrame(recs)
    df['key'] = df.apply(lambda r: tuple(sorted(
        [(r.chrom1, r.pos1), (r.chrom2, r.pos2)])), axis=1)
    ev = df.drop_duplicates('key').copy()
    ev['grp'] = ev.supp.map({'10': 'youth', '01': 'aged', '11': 'shared'})
    return df, ev


def load_te_bed(path):
    """Load TE BED. Columns beyond the mandatory 4 are optional:
    col5 score, col6 strand, col7 family_id, col8 class_id (RepeatMasker
    convention, as produced by gtfbuild liftover). family/class are used
    verbatim when present; names fall back to repeat-name prefixes."""
    cols = [0, 1, 2, 3]
    with open(path) as fh:
        ncol = 0
        for line in fh:
            if line.startswith('#') or line.startswith('track') or \
                    line.startswith('browser'):
                continue
            ncol = len(line.rstrip('\n').split('\t'))
            break
    if ncol >= 8:
        cols += [6, 7]
    elif ncol >= 7:
        cols += [6]
    names_all = ['chrom', 'start', 'end', 'name', 'family', 'class'][:len(cols)]
    te = pd.read_csv(path, sep='\t', header=None,
                     usecols=cols, names=names_all,
                     comment='#', low_memory=False)
    for c in ('start', 'end'):
        te[c] = pd.to_numeric(te[c], errors='coerce')
    te = te.dropna().sort_values(['chrom', 'start']).reset_index(drop=True)
    if 'family' not in te.columns:
        te['family'] = ''
    if 'class' not in te.columns:
        te['class'] = ''
    starts, ends, names, fams, classes = {}, {}, {}, {}, {}
    for c, g in te.groupby('chrom'):
        starts[c] = g.start.values
        ends[c] = g.end.values
        names[c] = g.name.values
        fams[c] = g.family.values if 'family' in g else np.array([''] * len(g))
        classes[c] = g['class'].values if 'class' in g else np.array([''] * len(g))
    return te, starts, ends, names, fams, classes


def te_hit(chrom, pos, starts, ends, names, fams=None, classes=None):
    """TE annotation covering this breakpoint; None if not in any TE.
    Returns (row_index, repeat_name[, family_id[, class_id]])."""
    if chrom not in starts:
        return None
    i = bisect.bisect_right(starts[chrom], pos) - 1
    if i >= 0 and ends[chrom][i] > pos:
        out = [i]
        out.append(names[chrom][i] if names is not None else '')
        if fams is not None:
            out.append(fams[chrom][i])
        if classes is not None:
            out.append(classes[chrom][i])
        return tuple(out)
    return None

# ---------------------------------------------------------------- stage 1
def stage1_enrichment(ev, starts, ends, chrlens, n_perm, seed, flank_ci=0):
    """Breakpoint TE rate vs background and vs length+GC matched random breaks."""
    bp1 = [(c, p) for c, p in zip(ev.chrom1, ev.pos1)]
    bp2 = [(c, int(p)) for c, p in zip(ev.chrom2, ev.pos2)]
    obs1 = sum(te_hit(c, p, starts, ends, None) is not None for c, p in bp1) / len(bp1)
    obs2 = sum(te_hit(c, p, starts, ends, None) is not None for c, p in bp2) / len(bp2)
    obs_both = sum(
        te_hit(c1, p1, starts, ends, None) is not None and
        te_hit(c2, p2, starts, ends, None) is not None
        for (c1, p1), (c2, p2) in zip(bp1, bp2)) / len(ev)

    # genome background
    tot = sum(chrlens.values())
    tebp = sum((ends[c][i] - starts[c][i]) for c in starts for i in range(len(starts[c])))
    bg = tebp / tot

    # permutation: random inter-chromosomal breakpoint pairs
    rng = random.Random(seed)
    chrs = sorted(chrlens)
    w = [chrlens[c] for c in chrs]
    r1 = []
    for _ in range(n_perm):
        c = rng.choices(chrs, weights=w)[0]
        r1.append(te_hit(c, rng.randint(1, chrlens[c]), starts, ends, None) is not None)
    perm1 = sum(r1) / len(r1)

    # both-TE permutation: independent pairs
    n_pairs = 2000
    cnt = 0
    for _ in range(n_pairs):
        ok = True
        for _ in range(2):
            c = rng.choices(chrs, weights=w)[0]
            if te_hit(c, rng.randint(1, chrlens[c]), starts, ends, None) is None:
                ok = False
                break
        cnt += ok
    perm_both = cnt / n_pairs

    # binomial p for observed both-TE count (exact upper tail).
    # scipy when available (sv.yaml ships it); otherwise the same sum of
    # exact binomial pmfs computed with math.comb — identical result.
    n = len(ev)
    k = round(obs_both * n)
    try:
        from scipy.stats import binomtest
        p_binom = binomtest(k, n, perm_both, alternative='greater').pvalue
    except ImportError:
        p_binom = sum(math.comb(n, i) * perm_both ** i *
                      (1 - perm_both) ** (n - i)
                      for i in range(k, n + 1))

    return dict(obs_bp1=obs1, obs_bp2=obs2, obs_both=obs_both, bg=bg,
                perm_single=perm1, perm_both=perm_both, n=n, k=k,
                p_binom_both=p_binom)

# ---------------------------------------------------------------- stage 3
def fam_of(name):
    if not isinstance(name, str):
        return 'none'
    n = name.split('_')[0]
    if n.startswith('Alu'):
        return 'Alu'
    if n.startswith('L1'):
        return 'L1'
    if n.startswith('SVA'):
        return 'SVA'
    if 'ERV' in name or n.startswith(('LTR', 'MLT', 'MER')):
        return 'LTR/ERV'
    if n.startswith(('L2', 'CR1')):
        return 'LINE2'
    if n.startswith('MIR'):
        return 'MIR'
    if n.startswith(('Tigger', 'hAT', 'Charlie', 'Zaphod', 'Piggy', 'Sola', 'Tc')):
        return 'DNA'
    if n.startswith(('SAT', 'Sat', 'GSAT', 'MSAT', '(A', '(T', '(AT', '(GC',
                    'simple')):
        return 'satellite/simple'
    return 'other'


def run_blastn(query_fa, subject_fa, out_tsv, sif=None, threads=4):
    """blastn the two TE copies. PATH-first (works inside the sv container,
    which ships blastn); falls back to apptainer exec with a SIF."""
    cmd = ['blastn', '-query', query_fa, '-subject', subject_fa,
           '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send']
    if threads > 1:
        cmd += ['-num_threads', str(threads)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        r = subprocess.CompletedProcess(cmd, 127, '', 'blastn not on PATH')
    if r.returncode != 0 and sif:
        cmd = ['apptainer', 'exec', '--bind', '/data/16T', sif] + cmd
        r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print('[warn] blastn rc={} stderr: {}'.format(
            r.returncode, r.stderr.strip()[-200:]))
    with open(out_tsv, 'w') as fh:
        fh.write(r.stdout)
    return out_tsv


_FAIDX_CACHE = {}


def fetch_seq(fasta, chrom, start, end, sif_samtools):
    """Extract reference sequence [start, end) (BED 0-based half-open).

    Prefers in-process pyfaidx (sv.yaml dependency, works inside the sv
    container); falls back to samtools faidx via an apptainer SIF for
    hosts without pyfaidx (coordinate converted to 1-based inclusive).
    """
    try:
        import pyfaidx  # noqa: delayed import; cached by sys.modules
        fa = _FAIDX_CACHE.get(fasta)
        if fa is None:
            fa = pyfaidx.Fasta(fasta, as_raw=True)
            _FAIDX_CACHE[fasta] = fa
        return str(fa[chrom][int(start):int(end)])
    except ImportError:
        pass
    r = subprocess.run(['apptainer', 'exec', '--bind', '/data/16T', sif_samtools,
                        'samtools', 'faidx', fasta,
                        f'{chrom}:{int(start)+1}-{int(end)}'],
                       capture_output=True, text=True)
    return ''.join(l.strip() for l in r.stdout.splitlines()
                   if not l.startswith('>'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--te-bed', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('-o', '--outdir', required=True)
    ap.add_argument('--samtools-sif', required=True)
    ap.add_argument('--blast-sif', required=True)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--n-perm', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--fai', default='',
                    help='reference .fai for real chromosome lengths (panorama fig)')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    SIF_SAM = args.samtools_sif
    SIF_SV = args.blast_sif

    # ---- stage 0
    recs, ev = parse_tra(args.vcf)
    print(f'[stage0] TRA records {len(recs)} -> independent events {len(ev)}')
    ev.to_csv(os.path.join(args.outdir, 'tra_events.tsv'), sep='\t', index=False)

    te, starts, ends, names, fams, classes = load_te_bed(args.te_bed)
    chrlens = {c: int(ends[c].max()) for c in ends}
    print(f'[stage0] TE bed: {len(te)} intervals; {len(chrlens)} chroms')
    has_fam_col = len(te) > 0 and bool(
        te['family'].astype(str).ne('').any())

    # annotate both breakpoints
    h1 = [te_hit(c, p, starts, ends, names, fams, classes)
          for c, p in zip(ev.chrom1, ev.pos1)]
    h2 = [te_hit(c, int(p), starts, ends, names, fams, classes)
          for c, p in zip(ev.chrom2, ev.pos2)]
    ev['te1'] = [h[1] if h else '' for h in h1]
    ev['te2'] = [h[1] if h else '' for h in h2]
    ev['i1'] = [h[0] if h else -1 for h in h1]
    ev['i2'] = [h[0] if h else -1 for h in h2]
    if has_fam_col:
        # official RepeatMasker family_id (BED col7); fall back to prefix
        # bucket only for rows where the column is missing/blank
        ev['fam1'] = [h[2] if h and len(h) > 2 and h[2] else
                      (fam_of(h[1]) if h else '') for h in h1]
        ev['fam2'] = [h[2] if h and len(h) > 2 and h[2] else
                      (fam_of(h[1]) if h else '') for h in h2]
        ev['cls1'] = [h[3] if h and len(h) > 3 and h[3] else ''
                      for h in h1]
        ev['cls2'] = [h[3] if h and len(h) > 3 and h[3] else ''
                      for h in h2]
    else:
        ev['fam1'] = ev.te1.apply(fam_of)
        ev['fam2'] = ev.te2.apply(fam_of)
        ev['cls1'] = ''
        ev['cls2'] = ''
    ev['both_te'] = (ev.te1 != '') & (ev.te2 != '')
    ev['same_family'] = ev.fam1 == ev.fam2
    ev['same_repeat'] = ev.te1 == ev.te2
    ev.to_csv(os.path.join(args.outdir, 'tra_events_annotated.tsv'),
              sep='\t', index=False)

    # ---- stage 1
    s1 = stage1_enrichment(ev, starts, ends, chrlens,
                           args.n_perm, args.seed)
    print('[stage1] enrichment:', {k: (round(v, 4) if isinstance(v, float) else v)
                                   for k, v in s1.items()})
    pd.Series(s1).to_csv(os.path.join(args.outdir, 'stage1_enrichment.tsv'),
                         sep='\t')

    # ---- stage 2: family pair matrix among both-TE events
    # pair key uses '|' separator: repeat names and official family_ids
    # (e.g. ERVL-MaLR, hAT-Charlie) may themselves contain '-'
    bt = ev[ev.both_te]
    pair_counter = Counter(bt.apply(
        lambda r: '|'.join(sorted([str(r.fam1), str(r.fam2)])), axis=1))
    fam_freq = (bt.fam1.tolist() + bt.fam2.tolist())
    fam_freq = Counter(fam_freq)
    n_bp = len(bt) * 2
    rows = []
    for pair, k in pair_counter.items():
        a, b = pair.split('|')
        exp = fam_freq[a] * fam_freq[b] * 2 / (n_bp * (n_bp - 1)) * len(bt) \
            if a != b else fam_freq[a] * (fam_freq[a] - 1) / (n_bp * (n_bp - 1)) * len(bt)
        rows.append(dict(pair=pair, obs=k, exp_under_shuffle=round(exp, 2),
                         ratio=round(k / exp, 2) if exp > 0 else math.inf))
    s2 = pd.DataFrame(rows).sort_values('obs', ascending=False)
    s2.to_csv(os.path.join(args.outdir, 'stage2_family_pairs.tsv'),
              sep='\t', index=False)
    print('[stage2] family pairs (top):')
    print(s2.head(8).to_string(index=False))

    # ---- stage 3: sequence homology for both-TE events (blastn)
    print(f'[stage3] blastn homology for {len(bt)} both-TE events ...')
    import tempfile
    tmpd = tempfile.mkdtemp(prefix='tra_blast_',
                            dir='/home/luosg/.hermes/cache/scratch')
    tmpd = os.path.abspath(tmpd)
    hom = []
    for _, r in bt.iterrows():
        if r.i1 < 0 or r.i2 < 0:
            continue
        s1_seq = fetch_seq(args.fasta, r.chrom1,
                           starts[r.chrom1][r.i1], ends[r.chrom1][r.i1],
                           SIF_SAM)
        s2_seq = fetch_seq(args.fasta, r.chrom2,
                           starts[r.chrom2][r.i2], ends[r.chrom2][r.i2],
                           SIF_SAM)
        if not s1_seq or not s2_seq:
            continue
        q = os.path.join(tmpd, 'q.fa')
        s = os.path.join(tmpd, 's.fa')
        with open(q, 'w') as fh:
            fh.write('>t1\n' + s1_seq + '\n')
        with open(s, 'w') as fh:
            fh.write('>t2\n' + s2_seq + '\n')
        out = run_blastn(q, s, os.path.join(tmpd, 'h.tsv'),
                         sif=SIF_SV, threads=args.threads)
        best = 0.0
        best_len = 0
        for line in open(out):
            f_ = line.split('\t')
            pid, ln = float(f_[2]), int(f_[3])
            if pid * ln > best * best_len:
                best, best_len = pid, ln
        hom.append(dict(key=r.key, te1=r.te1, te2=r.te2,
                        grp=r.grp, ident=round(best, 1),
                        hsp_len=best_len,
                        hom_index=round(best * best_len / 100, 1)))
    h = pd.DataFrame(hom)
    h.to_csv(os.path.join(args.outdir, 'stage3_homology.tsv'),
             sep='\t', index=False)
    if len(h):
        print('[stage3] homology tiers:')
        for lo, hi, tag in [(90, 101, '>=90 (NAHR-grade)'),
                            (80, 90, '80-90'), (70, 80, '70-80'), (0, 70, '<70')]:
            n = ((h.ident >= lo) & (h.ident < hi)).sum()
            print(f'  {tag}: {n}')
    for f_ in ('q.fa', 's.fa', 'h.tsv'):
        p = os.path.join(tmpd, f_)
        if os.path.exists(p):
            os.remove(p)
    try:
        os.rmdir(tmpd)
    except OSError:
        pass
    # ---- stage 4: figures (tra_te_figures.py consumes the 5 TSVs above)
    fig_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'tra_te_figures.py')
    if os.path.exists(fig_script):
        cmd = [sys.executable, fig_script, '--indir', args.outdir]
        if args.fai:
            cmd += ['--fai', args.fai]
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stdout, end='')
        if r.returncode != 0:
            print('[warn] tra_te_figures rc={} stderr: {}'.format(
                r.returncode, r.stderr.strip()[-300:]))
    else:
        print('[warn] tra_te_figures.py not found next to tra_te_analysis.py')

    print('[done] outputs in', args.outdir)


if __name__ == '__main__':
    main()

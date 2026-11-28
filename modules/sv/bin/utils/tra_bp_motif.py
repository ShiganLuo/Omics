#!/usr/bin/env python3
"""Breakpoint-sequence motif mining for TRA events (stage 4-5, continues tra_te_analysis).

For every TRA breakpoint (tra_events_annotated.tsv from tra_te_analysis.py) this
script extracts reference windows, builds TWO background sets (uniform random
and GC+TE-context matched), and asks three questions:

  1. Background composition  — GC, CpG O/E, homopolymer, tandem-repeat content
     of breakpoint windows vs both backgrounds (Mann-Whitney).
  2. Motif enrichment        — curated structural motifs (G-quadruplex,
     L1 endonuclease consensus, Z-DNA, inverted/direct repeats, polypurine /
     H-DNA tracts, AT-rich MAR-like tracts, telomeric repeats, poly-A/T),
     canonical k-mer enrichment, and TF PWMs (HOMER .motif format) scanned on
     both strands; Fisher exact + BH-FDR against BOTH backgrounds
     (unmatched OR = total enrichment incl. TE-mediated; matched OR = residual
     enrichment beyond TE/GC composition).
  3. Junction microhomology  — reference-flank microhomology of the called join
     (form recovered from per-sample BND ALT when available) and the maximum
     over all four BND join geometries; null = matched background position
     pairs. Combined with the stage-3 TE homology tier this yields a per-event
     mechanism call (TE-NAHR / microhomology-alt-EJ / blunt-NHEJ-like).

Outputs (TSVs + window FASTAs) are consumed by tra_bp_motif_figures.py and by
the HOMER de novo step in sv_tra_bp_motif.sh.

Usage:
  python tra_bp_motif.py \
      --events tra_te/tra_events_annotated.tsv \
      --homology tra_te/stage3_homology.tsv \
      --fasta GCF_049350105.2_T2T-MMU8v2.0_genomic.fa \
      --te-bed rheMac10_rmsk_TE.lifted.bed \
      --telocen-bed Mmul_T2T.telocen_regions.bed \
      --motif-dir tra_te/bp_motif/motif_lib \
      --sample-vcf luanchao-11236.sv.vcf.gz --sample-vcf luanchao-21310.sv.vcf.gz \
      -o tra_te/bp_motif [--seed 42] [--max-events 0]
"""
import argparse
import bisect
import gzip
import math
import os
import random
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

RC_TABLE = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(seq):
    """Reverse complement."""
    return seq.translate(RC_TABLE)[::-1]


# ------------------------------------------------------------ FASTA access
class FaiFasta:
    """Random-access FASTA reader driven by the samtools .fai index."""

    def __init__(self, fasta, fai=None):
        self.path = fasta
        self.fai = fai or (fasta + '.fai')
        self.info = {}
        self.order = []
        with open(self.fai) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                self.info[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
                self.order.append(f[0])
        self._fh = open(fasta, 'rb')

    def chrom_len(self, chrom):
        return self.info[chrom][0]

    def fetch(self, chrom, start, end):
        """Fetch [start, end) 0-based half-open; clamps to chrom bounds."""
        length, offset, linebases, linewidth = self.info[chrom]
        start = max(0, start)
        end = min(length, end)
        if end <= start:
            return ''
        s_off = offset + (start // linebases) * linewidth + (start % linebases)
        e_off = offset + ((end - 1) // linebases) * linewidth + ((end - 1) % linebases)
        self._fh.seek(s_off)
        raw = self._fh.read(e_off - s_off + 1)
        return raw.replace(b'\n', b'').replace(b'\r', b'').decode().upper()


# ------------------------------------------------------- TE / telocen lookup
def load_te_bed(path):
    """Load TE BED (chrom, start, end, name[, ... family, class])."""
    cols = [0, 1, 2, 3]
    with open(path) as fh:
        ncol = 0
        for line in fh:
            if line.startswith(('#', 'track', 'browser')):
                continue
            ncol = len(line.rstrip('\n').split('\t'))
            break
    if ncol >= 8:
        cols += [6, 7]
    elif ncol >= 7:
        cols += [6]
    names = ['chrom', 'start', 'end', 'name', 'family', 'class'][:len(cols)]
    te = pd.read_csv(path, sep='\t', header=None, usecols=cols, names=names,
                     comment='#', low_memory=False)
    for c in ('start', 'end'):
        te[c] = pd.to_numeric(te[c], errors='coerce')
    te = te.dropna().sort_values(['chrom', 'start']).reset_index(drop=True)
    for c in ('family', 'class'):
        if c not in te.columns:
            te[c] = ''
    idx = {}
    for c, g in te.groupby('chrom'):
        idx[c] = (g.start.values, g.end.values, g.name.values,
                  g.family.values, g['class'].values)
    return idx


def te_hit(te_idx, chrom, pos):
    """(name, family, class) of the TE covering pos, or None."""
    if chrom not in te_idx:
        return None
    starts, ends, names, fams, classes = te_idx[chrom]
    i = bisect.bisect_right(starts, pos) - 1
    if i >= 0 and ends[i] > pos:
        return (names[i], fams[i], classes[i])
    return None


def load_telocen(path):
    """chrom -> list of (start, end, tag) for telomere/centromere regions."""
    out = defaultdict(list)
    if not path:
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            out[f[0]].append((int(f[1]), int(f[2]), f[3]))
    return out


def telocen_context(tc, chrom, pos):
    """(region_tag, dist_to_region) for a 1-based position."""
    best_tag, best_dist = 'none', None
    for start, end, tag in tc.get(chrom, []):
        if start <= pos <= end:
            return tag, 0
        d = min(abs(pos - start), abs(pos - end))
        if best_dist is None or d < best_dist:
            best_tag, best_dist = tag, d
    return 'none', (best_dist if best_dist is not None else -1)


# ------------------------------------------------------ structural motifs
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'W': '[AT]', 'S': '[CG]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}


def iupac_regex(pattern):
    return ''.join(IUPAC[c] for c in pattern.upper())


def max_run_simple(seq):
    """(max homopolymer run, max tandem-repeat run) for mono/di/tri/tetra units.

    Tandem scan walks one phase at a time (O(n) per phase) instead of
    extending from every position (O(n^2) on homopolymers).
    """
    max_homo = 0
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1] and seq[i] in 'ACGT':
            run += 1
            max_homo = max(max_homo, run)
        else:
            run = 1
    max_tandem = 0
    n = len(seq)
    for unit_len in (2, 3, 4):
        for phase in range(unit_len):
            j = phase
            while j + unit_len <= n:
                unit = seq[j:j + unit_len]
                if 'N' in unit:
                    j += unit_len
                    continue
                run = 1
                k = j + unit_len
                while k + unit_len <= n and seq[k:k + unit_len] == unit:
                    run += 1
                    k += unit_len
                if run > 1:
                    max_tandem = max(max_tandem, run * unit_len)
                j = k if run > 1 else j + unit_len
    return max_homo, max_tandem


def at_rich_frac(seq, win=30):
    """Max AT fraction in any `win`-bp sliding window (window clipped at edges)."""
    if len(seq) < win:
        win = len(seq)
    if win == 0:
        return 0.0
    best = 0.0
    for i in range(0, len(seq) - win + 1):
        sub = seq[i:i + win]
        frac = sum(sub.count(b) for b in 'AT') / win
        best = max(best, frac)
    return best


def purine_tract(seq, min_len=15):
    """Longest purine (A/G) run and whether it qualifies as an H-DNA substrate."""
    best = 0
    run = 0
    for ch in seq:
        if ch in 'AG':
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


COMP_BASE = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def scan_inverted_repeats(seq, arm=8, max_spacer=12):
    """Count inverted repeats (palindromic stems) and best stem length.

    Seed: exact `arm`-mer whose downstream partner equals its reverse
    complement with 0..max_spacer bp between the arms; each seed is extended
    outwards one base at a time while the stem stays perfectly complementary.
    """
    n = len(seq)
    pos_of = defaultdict(list)
    for i in range(n - arm + 1):
        pos_of[seq[i:i + arm]].append(i)
    best_stem = 0
    hits = 0
    seen = set()
    for i in range(n - arm + 1):
        target = rc(seq[i:i + arm])
        for j in pos_of.get(target, ()):
            spacer = j - (i + arm)
            if spacer < 0 or spacer > max_spacer:
                continue
            x, y = i, j + arm          # left arm start / right arm end (exclusive)
            while x - 1 >= 0 and y < n and seq[x - 1] == COMP_BASE.get(seq[y], ''):
                x -= 1
                y += 1
            if (x, y) in seen:
                continue
            seen.add((x, y))
            hits += 1
            best_stem = max(best_stem, y - j + i + arm - x)
    return hits, best_stem


def scan_direct_repeats(seq, arm=8, max_spacer=20):
    """Count direct repeats: exact `arm`-mer duplicated with 1..max_spacer gap."""
    n = len(seq)
    pos_of = defaultdict(list)
    for i in range(n - arm + 1):
        pos_of[seq[i:i + arm]].append(i)
    best = 0
    hits = 0
    seen = set()
    for i in range(n - arm + 1):
        for j in pos_of.get(seq[i:i + arm], ()):
            spacer = j - (i + arm)
            if spacer < 1 or spacer > max_spacer:
                continue
            a, b = i, j               # start of unit1 / unit2
            while a - 1 >= 0 and b - 1 >= i + arm and seq[a - 1] == seq[b - 1]:
                a -= 1
                b -= 1
            c, d = i + arm, j + arm   # end of unit1 / unit2 (exclusive)
            while c < b and d < n and seq[c] == seq[d]:
                c += 1
                d += 1
            if (a, b) in seen:
                continue
            seen.add((a, b))
            hits += 1
            best = max(best, c - a)
    return hits, best


RE_G4 = re.compile(r'(?:G{3,}[ACGT]{1,7}){3}G{3,}')
RE_G4C = re.compile(r'(?:C{3,}[ACGT]{1,7}){3}C{3,}')
RE_L1EN_STRICT = re.compile(r'TTTTAA|TTAAAA')
RE_L1EN_RELAX = re.compile(r'TTTTA|TAAAA')
RE_ZDNA = re.compile(r'(?:CG){4,}|(?:CA){5,}|(?:TG){5,}')
RE_POLYAT = re.compile(r'A{10,}|T{10,}')
RE_AG_MIRROR = re.compile(r'(?:AG){6,}|(?:GA){6,}')
RE_TELO = re.compile(r'(?:TTAGGG){2,}|(?:CCCTAA){2,}')
RE_TELO_ANY = re.compile(r'TTAGGG|CCCTAA')
RE_ATRICH = None  # computed via at_rich_frac


def window_features(seq):
    """Structural motif feature dict for one window."""
    gc = (seq.count('G') + seq.count('C')) / max(1, len(seq))
    cpg = seq.count('CG')
    cpg_oe = (cpg * len(seq)) / max(1, (seq.count('C') * seq.count('G')))
    max_homo, max_tandem = max_run_simple(seq)
    pur = purine_tract(seq)
    ir_hits, ir_stem = scan_inverted_repeats(seq)
    dr_hits, dr_span = scan_direct_repeats(seq)
    at_rich = at_rich_frac(seq)
    return {
        'gc': gc,
        'cpg_oe': cpg_oe if (seq.count('C') * seq.count('G')) > 0 else float('nan'),
        'max_homopolymer': max_homo,
        'max_tandem': max_tandem,
        'purine_tract': pur,
        'ir_hits': ir_hits,
        'ir_stem': ir_stem,
        'dr_hits': dr_hits,
        'dr_span': dr_span,
        'at_rich_max': at_rich,
        'g4': len(RE_G4.findall(seq)),
        'g4c': len(RE_G4C.findall(seq)),
        'l1en_strict': len(RE_L1EN_STRICT.findall(seq)),
        'l1en_relax': len(RE_L1EN_RELAX.findall(seq)),
        'zdna': len(RE_ZDNA.findall(seq)),
        'polyat': len(RE_POLYAT.findall(seq)),
        'ag_mirror': len(RE_AG_MIRROR.findall(seq)),
        'telo': len(RE_TELO.findall(seq)),
        'telo_any': len(RE_TELO_ANY.findall(seq)),
    }


# ------------------------------------------------------------- PWM scanning
def parse_homer_motifs(motif_dir):
    """Parse HOMER .motif files -> {name: (L, 4) probability matrix (A,C,G,T)}."""
    out = {}
    files = sorted(f for f in os.listdir(motif_dir) if f.endswith('.motif'))
    for fn in files:
        rows = []
        name = fn[:-6]
        with open(os.path.join(motif_dir, fn)) as fh:
            for line in fh:
                if line.startswith('>'):
                    parts = line[1:].strip().split('\t')
                    name = parts[1] if len(parts) > 1 and parts[1] else parts[0]
                    continue
                vals = line.split()
                if len(vals) == 4:
                    try:
                        rows.append([float(v) for v in vals])
                    except ValueError:
                        continue
        if rows:
            out[name] = np.array(rows, dtype=float)
    return out


def pwm_logodds(mat, pseudocount=1e-3):
    """log2 odds matrix vs uniform background, plus max attainable score."""
    p = np.clip(mat, pseudocount, 1.0)
    lo = np.log2(p / 0.25)
    return lo, float(lo.max(axis=1).sum())


def scan_pwm(seq, lo, max_score, threshold=0.8):
    """Best log-odds hit (both strands) and its distance to the window center.

    Sliding-window scores are accumulated per motif position with numpy slice
    additions (vectorised over the window axis).
    """
    cut = threshold * max_score
    L = lo.shape[0]
    n = len(seq)
    if n < L:
        return -math.inf, None
    idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    codes = np.array([idx.get(ch, -1) for ch in seq], dtype=int)
    comp = np.array([idx.get(ch, -1) for ch in rc(seq)], dtype=int)

    def best_one(cod):
        nw = n - L + 1
        scores = np.zeros(nw)
        for j in range(L):
            window = cod[j:j + nw]
            ok = window >= 0
            contrib = np.where(ok, lo[j, np.where(ok, window, 0)], -10.0)
            scores += contrib
        i = int(np.argmax(scores))
        return float(scores[i]), i

    s1, i1 = best_one(codes)
    s2, i2 = best_one(comp)
    best = max(s1, s2)
    if best < cut:
        return best, None
    i = i1 if s1 >= s2 else i2
    center = n // 2
    dist = abs(i + L // 2 - center)
    return best, dist


# ------------------------------------------------------------- k-mer stats
def kmer_counts(seqs, k):
    """Canonical (strand-merged) k-mer occurrence counts."""
    cnt = Counter()
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if 'N' in kmer:
                continue
            cnt[min(kmer, rc(kmer))] += 1
    return cnt


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR over a list of p-values."""
    n = len(pvals)
    order = np.argsort(pvals)
    fdr = np.empty(n)
    prev = 1.0
    for rank_from_end, i in enumerate(order[::-1]):
        rank = n - rank_from_end
        val = min(prev, pvals[i] * n / rank)
        fdr[i] = val
        prev = val
    return fdr.tolist()


def fisher_or(a_hit, a_tot, b_hit, b_tot):
    """Odds ratio + Fisher exact p (two-sided) for hit rates a vs b."""
    from scipy.stats import fisher_exact
    t = [[a_hit, a_tot - a_hit], [b_hit, b_tot - b_hit]]
    orr, p = fisher_exact(t, alternative='two-sided')
    # Woolf CI on log OR
    cells = [a_hit + 0.5, a_tot - a_hit + 0.5, b_hit + 0.5, b_tot - b_hit + 0.5]
    logor = math.log((cells[0] * cells[3]) / (cells[1] * cells[2]))
    se = math.sqrt(sum(1.0 / c for c in cells))
    ci_lo, ci_hi = math.exp(logor - 1.96 * se), math.exp(logor + 1.96 * se)
    return orr, p, ci_lo, ci_hi


# -------------------------------------------------- junction microhomology
def mh_len(left, right, max_k=30):
    """Longest suffix(left)==prefix(right) overlap up to max_k."""
    best = 0
    for k in range(1, min(max_k, len(left), len(right)) + 1):
        if left[-k:] == right[:k]:
            best = k
        else:
            break
    return best


def join_microhomology(fa, chrom1, pos1, chrom2, pos2, flank=40, max_k=30):
    """Microhomology of the four BND join geometries.

    x = left flank of bp1 (ends at pos1), y = right flank of bp1 (starts at
    pos1+1), u = left flank of bp2, v = right flank of bp2 (all 1-based VCF
    POS semantics). Geometry A = t[p[, B = t]p], C = ]p]t, D = [p[t.
    """
    x = fa.fetch(chrom1, pos1 - flank, pos1)
    y = fa.fetch(chrom1, pos1, pos1 + flank)
    u = fa.fetch(chrom2, pos2 - flank, pos2)
    v = fa.fetch(chrom2, pos2, pos2 + flank)
    rcu, rcv = rc(u), rc(v)
    return {
        'A': mh_len(x, v, max_k),        # ...x | v...
        'B': mh_len(x, rcu, max_k),      # ...x | rc(u)...
        'C': mh_len(u, y, max_k),        # ...u | y...
        'D': mh_len(rcv, y, max_k),      # ...rc(v) | y...
    }


FORM_TO_GEOM = {'t[p[': 'A', 't]p]': 'B', ']p]t': 'C', '[p[t': 'D'}


def parse_sample_bnd(vcf_path):
    """Mate coordinates + ALT join form for every BND record in a sample VCF."""
    recs = []
    op = gzip.open if vcf_path.endswith('.gz') else open
    with op(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if 'SVTYPE=BND' not in f[7] and 'SVTYPE=TRA' not in f[7]:
                continue
            m = re.search(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]', f[4])
            if not m:
                continue
            # normalise join form from bracket placement in the ALT string
            alt = f[4]
            if re.match(r'^[ACGTNacgtn]*\[\S+\[\S*$', alt):
                form = 't[p['
            elif re.match(r'^[ACGTNacgtn]*\]\S+\]\S*$', alt):
                form = 't]p]'
            elif re.match(r'^\[\S+\[[ACGTNacgtn]*$', alt):
                form = '[p[t'
            elif re.match(r'^\]\S+\][ACGTNacgtn]*$', alt):
                form = ']p]t'
            else:
                form = 't[p['
            recs.append((f[0], int(f[1]), m.group(1), int(m.group(2)), form))
    return recs


def recover_form(events, sample_recs, tol=200):
    """Per-event join form from per-sample BND records (coordinate-matched)."""
    by_chrom = defaultdict(list)
    for chrom1, pos1, chrom2, pos2, form in sample_recs:
        by_chrom[chrom1].append((pos1, chrom2, pos2, form))
        by_chrom[chrom2].append((pos2, chrom1, pos1, form))
    forms, sources = [], []
    for _, r in events.iterrows():
        got = None
        for pos, mchrom, mpos, form in by_chrom.get(r.chrom1, ()):
            if mchrom == r.chrom2 and abs(pos - r.pos1) <= tol \
                    and abs(mpos - r.pos2) <= tol:
                got = form
                break
        if got is None:
            for pos, mchrom, mpos, form in by_chrom.get(r.chrom2, ()):
                if mchrom == r.chrom1 and abs(pos - r.pos2) <= tol \
                        and abs(mpos - r.pos1) <= tol:
                    got = form
                    break
        forms.append(got or 't[p[')
        sources.append('sample_vcf' if got else 'merged_default')
    return forms, sources


# ----------------------------------------------------------- background set
def sample_background(fa, te_idx, breakpoints, k, rng, matched, exclude,
                      gc_tol=0.05, max_tries=100):
    """k random windows per breakpoint; matched = same TE-context + GC bin.

    `exclude` is a list of (chrom, pos) breakpoint coordinates; sampled
    positions must stay >=10kb away from all of them.
    """
    out = []  # (chrom, pos, source_bp_index)
    chrom_lens = {c: fa.chrom_len(c) for c in fa.info}
    chroms = sorted(chrom_lens)
    excl = defaultdict(list)
    for ec, ep in exclude:
        excl[ec].append(ep)
    for c in excl:
        excl[c].sort()

    def near_breakpoint(c, p, radius=10000):
        arr = excl.get(c)
        if not arr:
            return False
        i = bisect.bisect_left(arr, p)
        for j in (i - 1, i):
            if 0 <= j < len(arr) and abs(arr[j] - p) < radius:
                return True
        return False

    for bi, (chrom, pos, _te) in enumerate(breakpoints):
        te_here = te_hit(te_idx, chrom, pos) is not None
        gc_here = None
        if matched:
            w = fa.fetch(chrom, pos - 250, pos + 251)
            if w:
                gc_here = (w.count('G') + w.count('C')) / len(w)
        got = 0
        best = None
        for _ in range(max_tries):
            c = chrom if matched else rng.choice(chroms)
            p = rng.randint(1000, chrom_lens[c] - 1000)
            if near_breakpoint(c, p):
                continue
            if matched:
                if (te_hit(te_idx, c, p) is not None) != te_here:
                    continue
                w = fa.fetch(c, p - 250, p + 251)
                if not w:
                    continue
                gc_cand = (w.count('G') + w.count('C')) / len(w)
                if gc_here is not None:
                    dev = abs(gc_cand - gc_here)
                    if dev > gc_tol:
                        if best is None or dev < best[0]:
                            best = (dev, c, p)
                        continue
            out.append((c, p, bi))
            got += 1
            if got >= k:
                break
        if got < k and matched and best is not None:
            out.append((best[1], best[2], bi))
    return out


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-e', '--events', required=True,
                    help='tra_events_annotated.tsv from tra_te_analysis.py')
    ap.add_argument('-H', '--homology', default='',
                    help='stage3_homology.tsv (NAHR tier per event)')
    ap.add_argument('-f', '--fasta', required=True, help='reference genome FASTA')
    ap.add_argument('-x', '--fai', default='', help='FASTA .fai (default <fasta>.fai)')
    ap.add_argument('-t', '--te-bed', required=True, help='TE annotation BED')
    ap.add_argument('-c', '--telocen-bed', default='',
                    help='telomere/centromere region BED (optional)')
    ap.add_argument('-m', '--motif-dir', default='',
                    help='directory of HOMER .motif files for TF PWM scanning')
    ap.add_argument('-s', '--sample-vcf', action='append', default=[],
                    help='per-sample SV VCF (repeatable) to recover join form')
    ap.add_argument('-o', '--outdir', required=True)
    ap.add_argument('-w', '--window', type=int, default=250,
                    help='half-window around the breakpoint (default 250 -> 501bp)')
    ap.add_argument('-k', '--bg-per-bp', type=int, default=10,
                    help='background windows per breakpoint per background set')
    ap.add_argument('-K', '--kmer', type=int, default=6,
                    help='k for k-mer enrichment (default 6)')
    ap.add_argument('-r', '--seed', type=int, default=42)
    ap.add_argument('-n', '--max-events', type=int, default=0,
                    help='smoke test: analyse only the first N events')
    ap.add_argument('--pwm-threshold', type=float, default=0.8,
                    help='PWM hit threshold x max_score (default 0.8)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)
    fa = FaiFasta(args.fasta, args.fai or None)
    te_idx = load_te_bed(args.te_bed)
    tc = load_telocen(args.telocen_bed)

    ev = pd.read_csv(args.events, sep='\t')
    if args.max_events:
        ev = ev.head(args.max_events).copy()
    ev['event_id'] = np.arange(len(ev))
    print(f'[stage0] events {len(ev)} breakpoints {2 * len(ev)}')

    # ---- breakpoint table -------------------------------------------------
    rows = []
    for _, r in ev.iterrows():
        for side, chrom, pos in ((1, r.chrom1, int(r.pos1)), (2, r.chrom2, int(r.pos2))):
            hit = te_hit(te_idx, chrom, pos)
            tag, dist = telocen_context(tc, chrom, pos)
            rows.append(dict(event_id=r.event_id, side=side, chrom=chrom, pos=pos,
                             grp=r.grp, te_name=hit[0] if hit else '',
                             te_fam=hit[1] if hit else '',
                             te_cls=hit[2] if hit else '',
                             te_hit=hit is not None,
                             telocen_tag=tag, telocen_dist=dist))
    bps = pd.DataFrame(rows)

    # ---- cluster near-duplicate breakpoints (500bp) for independent stats ---
    bps['cluster'] = -1
    cid = 0
    for chrom, g in bps.groupby('chrom'):
        last = None
        for idx_, pos in sorted(zip(g.index, g.pos), key=lambda t: t[1]):
            if last is None or pos - last > 500:
                cid += 1
                last = pos
            bps.loc[idx_, 'cluster'] = cid
    reps = bps.drop_duplicates('cluster')[['chrom', 'pos', 'te_hit', 'cluster']]
    print(f'[stage0] breakpoint clusters (<=500bp merged): {len(reps)}')

    # ---- windows ----------------------------------------------------------
    half = args.window
    rep_list = [(r.chrom, r.pos, r.te_hit) for r in reps.itertuples()]
    exclude = [(r.chrom, r.pos) for r in reps.itertuples()]
    cluster_to_bi = {int(c): bi for bi, c in enumerate(reps.cluster)}

    def make_windows(positions, label):
        seqs, beds, kept = [], [], []
        for pi, (chrom, pos, _te) in enumerate(positions):
            s = fa.fetch(chrom, pos - half, pos + half + 1)
            if len(s) < 2 * half * 0.6:
                continue
            seqs.append(s)
            beds.append((chrom, pos - half, pos + half + 1, label, pos))
            kept.append(pi)
        return seqs, beds, kept

    bp_seqs, bp_bed, bp_kept = make_windows(rep_list, 'breakpoint')
    print(f'[stage4] breakpoint windows: {len(bp_seqs)} x {2 * half + 1}bp')

    bg_m = sample_background(fa, te_idx, rep_list, args.bg_per_bp, rng,
                             matched=True, exclude=exclude)
    bg_u = sample_background(fa, te_idx, rep_list, args.bg_per_bp, rng,
                             matched=False, exclude=exclude)
    bgm_seqs, bgm_bed, bgm_kept = make_windows(bg_m, 'bg_matched')
    bgu_seqs, bgu_bed, bgu_kept = make_windows(bg_u, 'bg_uniform')
    print(f'[stage4] background windows: matched {len(bgm_seqs)} uniform {len(bgu_seqs)}')

    with open(os.path.join(args.outdir, 'windows.bed'), 'w') as fh:
        for bed in (bp_bed, bgm_bed, bgu_bed):
            for chrom, s, e, label, pos in bed:
                fh.write(f'{chrom}\t{s}\t{e}\t{label}|{chrom}:{pos}\n')
    for fname, seqs in (('windows_bp.fa', bp_seqs),
                        ('windows_bg_matched.fa', bgm_seqs),
                        ('windows_bg_uniform.fa', bgu_seqs)):
        with open(os.path.join(args.outdir, fname), 'w') as fh:
            for i, s in enumerate(seqs):
                fh.write(f'>{fname[:-3]}_{i}\n{s}\n')

    # ---- stage 4: composition --------------------------------------------
    comp_keys = ['gc', 'cpg_oe', 'max_homopolymer', 'max_tandem', 'purine_tract',
                 'at_rich_max']
    feats_bp = pd.DataFrame([window_features(s) for s in bp_seqs])
    feats_m = pd.DataFrame([window_features(s) for s in bgm_seqs])
    feats_u = pd.DataFrame([window_features(s) for s in bgu_seqs])
    feats_bp.insert(0, 'label', 'breakpoint')
    feats_m.insert(0, 'label', 'bg_matched')
    feats_u.insert(0, 'label', 'bg_uniform')
    pd.concat([feats_bp, feats_m, feats_u]).to_csv(
        os.path.join(args.outdir, 'windows_features.tsv'), sep='\t', index=False)

    from scipy.stats import mannwhitneyu
    comp_rows = []
    for key in comp_keys:
        a = feats_bp[key].dropna().values
        b = feats_m[key].dropna().values
        c = feats_u[key].dropna().values
        p_m = mannwhitneyu(a, b).pvalue if len(a) > 1 and len(b) > 1 else float('nan')
        p_u = mannwhitneyu(a, c).pvalue if len(a) > 1 and len(c) > 1 else float('nan')
        comp_rows.append(dict(feature=key, mean_bp=a.mean(), mean_bg_matched=b.mean(),
                              mean_bg_uniform=c.mean(), p_matched=p_m, p_uniform=p_u))
    pd.DataFrame(comp_rows).to_csv(
        os.path.join(args.outdir, 'stage4_composition.tsv'), sep='\t', index=False)

    # ---- stage 5a: structural motif enrichment ---------------------------
    struct_keys = ['g4', 'g4c', 'l1en_strict', 'l1en_relax', 'zdna', 'polyat',
                   'ag_mirror', 'telo', 'telo_any', 'ir_hits', 'dr_hits']
    kind_map = {
        'g4': 'G4 (G-rich)', 'g4c': 'G4 (C-rich)', 'l1en_strict': 'L1 EN core',
        'l1en_relax': 'L1 EN relaxed', 'zdna': 'Z-DNA', 'polyat': 'poly-A/T',
        'ag_mirror': 'AG mirror (H-DNA)', 'telo': 'telomere repeat x2',
        'telo_any': 'telomere repeat', 'ir_hits': 'inverted repeat',
        'dr_hits': 'direct repeat',
    }
    n_bp, n_m, n_u = len(feats_bp), len(feats_m), len(feats_u)
    rows5 = []
    for key in struct_keys:
        a_hit = int((feats_bp[key] > 0).sum())
        m_hit = int((feats_m[key] > 0).sum())
        u_hit = int((feats_u[key] > 0).sum())
        or_m, p_m, lo_m, hi_m = fisher_or(a_hit, n_bp, m_hit, n_m)
        or_u, p_u, lo_u, hi_u = fisher_or(a_hit, n_bp, u_hit, n_u)
        rows5.append(dict(motif=kind_map[key], kind='structural', feature=key,
                          hits_bp=a_hit, rate_bp=a_hit / n_bp,
                          hits_bg_matched=m_hit, rate_bg_matched=m_hit / n_m,
                          or_bg_matched=or_m, ci_low_matched=lo_m,
                          ci_high_matched=hi_m, p_bg_matched=p_m,
                          hits_bg_uniform=u_hit, rate_bg_uniform=u_hit / n_u,
                          or_bg_uniform=or_u, p_bg_uniform=p_u,
                          median_dist_bp=''))

    # ---- stage 5b: TF PWM enrichment ------------------------------------
    tf_rows = []
    tf_hitvec = {'bp': {}, 'm': {}, 'u': {}}
    if args.motif_dir:
        motifs = parse_homer_motifs(args.motif_dir)
        print(f'[stage5] TF PWMs: {len(motifs)} motifs from {args.motif_dir}')
        for name, mat in sorted(motifs.items()):
            lo, max_score = pwm_logodds(mat)
            out = {}
            for tag, seqs in (('bp', bp_seqs), ('m', bgm_seqs), ('u', bgu_seqs)):
                hits, dists, vec = 0, [], []
                for s in seqs:
                    best, dist = scan_pwm(s, lo, max_score, args.pwm_threshold)
                    hit = dist is not None
                    vec.append(hit)
                    if hit:
                        hits += 1
                        dists.append(dist)
                out[tag] = (hits, dists)
                tf_hitvec[tag][name] = np.array(vec, dtype=bool)
            a_hit, dists = out['bp']
            m_hit = out['m'][0]
            u_hit = out['u'][0]
            or_m, p_m, lo_m, hi_m = fisher_or(a_hit, n_bp, m_hit, n_m)
            or_u, p_u, lo_u, hi_u = fisher_or(a_hit, n_bp, u_hit, n_u)
            med = float(np.median(dists)) if dists else float('nan')
            tf_rows.append(dict(motif=name, kind='TF_PWM', feature=name,
                                hits_bp=a_hit, rate_bp=a_hit / n_bp,
                                hits_bg_matched=m_hit, rate_bg_matched=m_hit / n_m,
                                or_bg_matched=or_m, ci_low_matched=lo_m,
                                ci_high_matched=hi_m, p_bg_matched=p_m,
                                hits_bg_uniform=u_hit, rate_bg_uniform=u_hit / n_u,
                                or_bg_uniform=or_u, p_bg_uniform=p_u,
                                median_dist_bp=med))
    t5 = pd.DataFrame(rows5 + tf_rows)
    t5['fdr_bg_matched'] = bh_fdr(t5['p_bg_matched'].tolist())
    t5['fdr_bg_uniform'] = bh_fdr(t5['p_bg_uniform'].tolist())
    t5.sort_values('p_bg_matched').to_csv(
        os.path.join(args.outdir, 'stage5_motif_enrichment.tsv'), sep='\t', index=False)

    # ---- augmented per-window table: features + TE class + TF hits -------
    pos_lists = ([(rep_list[pi][0], rep_list[pi][1]) for pi in bp_kept],
                 [(bg_m[pi][0], bg_m[pi][1]) for pi in bgm_kept],
                 [(bg_u[pi][0], bg_u[pi][1]) for pi in bgu_kept])
    all_pos = pos_lists[0] + pos_lists[1] + pos_lists[2]
    all_feats = pd.concat([feats_bp, feats_m, feats_u], ignore_index=True)
    all_feats.insert(1, 'chrom', [p[0] for p in all_pos])
    all_feats.insert(2, 'pos', [p[1] for p in all_pos])

    def _te_class(cp):
        h = te_hit(te_idx, cp[0], cp[1])
        return (h[2] or 'nonTE') if h else 'nonTE'

    all_feats['te_cls'] = [_te_class(p) for p in all_pos]
    for name in sorted(tf_hitvec['bp']):
        col = 'TF_' + re.sub(r'[^A-Za-z0-9]+', '_', name)[:24]
        all_feats[col] = (list(tf_hitvec['bp'][name]) +
                          list(tf_hitvec['m'].get(name, [None] * len(bgm_seqs))) +
                          list(tf_hitvec['u'].get(name, [None] * len(bgu_seqs))))
    all_feats.to_csv(os.path.join(args.outdir, 'windows_all.tsv'), sep='\t',
                     index=False)

    # ---- stage 5c: k-mer enrichment -------------------------------------
    k = args.kmer
    cnt_bp = kmer_counts(bp_seqs, k)
    cnt_m = kmer_counts(bgm_seqs, k)
    cnt_u = kmer_counts(bgu_seqs, k)
    tot_bp = sum(cnt_bp.values())
    tot_m = sum(cnt_m.values())
    tot_u = sum(cnt_u.values())
    krows = []
    for kmer in sorted(set(cnt_bp) | set(cnt_m) | set(cnt_u)):
        a, b, c = cnt_bp.get(kmer, 0), cnt_m.get(kmer, 0), cnt_u.get(kmer, 0)
        if a + b + c < 5:
            continue
        or_m, p_m, lo_m, hi_m = fisher_or(a, tot_bp, b, tot_m)
        or_u, p_u, lo_u, hi_u = fisher_or(a, tot_bp, c, tot_u)
        krows.append(dict(kmer=kmer, count_bp=a, rate_bp=a / tot_bp,
                          count_bg_matched=b, rate_bg_matched=b / tot_m,
                          or_bg_matched=or_m, p_bg_matched=p_m,
                          count_bg_uniform=c, rate_bg_uniform=c / tot_u,
                          or_bg_uniform=or_u, p_bg_uniform=p_u))
    t6 = pd.DataFrame(krows)
    if len(t6):
        t6['fdr_bg_matched'] = bh_fdr(t6['p_bg_matched'].tolist())
        t6 = t6.sort_values(['fdr_bg_matched', 'or_bg_matched'],
                            ascending=[True, False])
    t6.to_csv(os.path.join(args.outdir, 'stage6_kmer.tsv'), sep='\t', index=False)
    print(f'[stage5] k-mer table: {len(t6)} canonical {k}-mers')

    # ---- stage 7: junction microhomology + mechanism ---------------------
    hom = {}
    if args.homology and os.path.exists(args.homology):
        h = pd.read_csv(args.homology, sep='\t')
        for _, r in h.iterrows():
            hom[str(r['key'])] = (float(r.ident), float(r.hsp_len),
                                  float(r.hom_index))
    sample_recs = []
    for vcf in args.sample_vcf:
        sample_recs.extend(parse_sample_bnd(vcf))
    forms, sources = recover_form(ev, sample_recs)

    # matched-background null pairs (one per event, from that bp's matched set)
    mgmt_by_src = defaultdict(list)
    for c, p, bi in bg_m:
        mgmt_by_src[bi].append((c, p))
    clu1 = bps[bps.side == 1].set_index('event_id')['cluster'].to_dict()
    clu2 = bps[bps.side == 2].set_index('event_id')['cluster'].to_dict()

    mh_rows = []
    for j, (_, r) in enumerate(ev.iterrows()):
        geom = join_microhomology(fa, r.chrom1, int(r.pos1), r.chrom2, int(r.pos2))
        form = forms[j]
        called = geom[FORM_TO_GEOM.get(form, 'A')]
        max_geom = max(geom, key=geom.get)
        # null pair drawn from the two breakpoints' matched backgrounds
        null = None
        cand1 = mgmt_by_src.get(cluster_to_bi.get(int(clu1.get(r.event_id, -1)), -1), [])
        cand2 = mgmt_by_src.get(cluster_to_bi.get(int(clu2.get(r.event_id, -1)), -1), [])
        if cand1 and cand2:
            c1, p1 = rng.choice(cand1)
            c2, p2 = rng.choice(cand2)
            ng = join_microhomology(fa, c1, p1, c2, p2)
            null = max(ng.values())
        hrow = hom.get(str(r.get('key', '')), (float('nan'), float('nan'),
                                               float('nan')))
        mh_rows.append(dict(
            event_id=r.event_id, key=r.get('key', ''), grp=r.grp,
            chrom1=r.chrom1, pos1=r.pos1, chrom2=r.chrom2, pos2=r.pos2,
            join_form=form, form_source=sources[j],
            mh_A=geom['A'], mh_B=geom['B'], mh_C=geom['C'], mh_D=geom['D'],
            mh_called=called, mh_max=geom[max_geom], mh_max_geom=max_geom,
            te_ident=hrow[0], te_hsp_len=hrow[1], te_hom_index=hrow[2],
            null_mh_max=null))
    mh = pd.DataFrame(mh_rows)

    # mechanism call
    def mech(r):
        nahr = (not math.isnan(r.te_ident)) and r.te_ident >= 90 and r.te_hsp_len >= 100
        if nahr:
            return 'TE-NAHR'
        if r.mh_called >= 4:
            return 'microhomology-altEJ'
        if r.mh_called <= 1:
            return 'blunt-NHEJ-like'
        return 'short-MH'

    mh['nahr'] = (~mh.te_ident.isna()) & (mh.te_ident >= 90) & (mh.te_hsp_len >= 100)
    mh['mechanism'] = mh.apply(mech, axis=1)

    # structural tags per breakpoint side
    tag_map = {}
    for _, b in bps.iterrows():
        w = fa.fetch(b.chrom, b.pos - half, b.pos + half + 1)
        feats = window_features(w)
        tag_map[(b.event_id, b.side)] = feats
    for key in ('g4', 'zdna', 'ag_mirror', 'telo_any'):
        mh[f'{key}_bp1'] = [tag_map[(r.event_id, 1)][key] > 0 for r in mh.itertuples()]
        mh[f'{key}_bp2'] = [tag_map[(r.event_id, 2)][key] > 0 for r in mh.itertuples()]
    mh['structure_tag'] = (mh.g4_bp1 | mh.g4_bp2 | mh.zdna_bp1 | mh.zdna_bp2 |
                           mh.ag_mirror_bp1 | mh.ag_mirror_bp2)
    mh['telo_tag'] = mh.telo_any_bp1 | mh.telo_any_bp2
    mh['l1en_tag'] = [
        tag_map[(r.event_id, 1)]['l1en_relax'] > 0 or
        tag_map[(r.event_id, 2)]['l1en_relax'] > 0 for r in mh.itertuples()]

    # ---- stage 5d: mechanism-stratified motif presence --------------------
    # bp window i (bp_kept order) <-> reps row bp_kept[i] -> cluster -> mechanism
    prio = {'TE-NAHR': 3, 'microhomology-altEJ': 2, 'short-MH': 1,
            'blunt-NHEJ-like': 0}
    mech_by_cluster = {}
    for row in mh.itertuples():
        for clu in (clu1.get(row.event_id), clu2.get(row.event_id)):
            if clu is None:
                continue
            clu = int(clu)
            cur = mech_by_cluster.get(clu)
            if cur is None or prio[row.mechanism] > prio[cur]:
                mech_by_cluster[clu] = row.mechanism
    clu_of_window = [int(reps.iloc[pi].cluster) for pi in bp_kept]
    mech_of_window = np.array([mech_by_cluster.get(c, 'unknown')
                               for c in clu_of_window])
    nahr_mask = mech_of_window == 'TE-NAHR'
    other_mask = np.isin(mech_of_window, list(prio)) & ~nahr_mask
    strat_sets = {key: (feats_bp[key].values > 0) for key in struct_keys}
    strat_sets.update({name: vec for name, vec in tf_hitvec['bp'].items()})
    strat_rows = []
    for key, vec in strat_sets.items():
        a, na = int(vec[nahr_mask].sum()), int(nahr_mask.sum())
        b, nb = int(vec[other_mask].sum()), int(other_mask.sum())
        if na and nb:
            orr, p, lo, hi = fisher_or(a, na, b, nb)
        else:
            orr, p, lo, hi = (float('nan'),) * 4
        strat_rows.append(dict(motif=key, hits_nahr=a, n_nahr=na,
                               rate_nahr=a / max(1, na),
                               hits_non_nahr=b, n_non_nahr=nb,
                               rate_non_nahr=b / max(1, nb),
                               or_nahr_vs_non=orr, ci_low=lo, ci_high=hi, p=p))
    t5s = pd.DataFrame(strat_rows)
    if len(t5s):
        t5s['fdr'] = bh_fdr(t5s['p'].fillna(1.0).tolist())
        t5s = t5s.sort_values('p')
    t5s.to_csv(os.path.join(args.outdir, 'stage5_stratified.tsv'), sep='\t',
               index=False)

    # ---- breakpoints.tsv (per-breakpoint features + context) -------------
    feat_df = pd.DataFrame([dict(event_id=eid, side=side, **feats)
                            for (eid, side), feats in sorted(tag_map.items())])
    bpdf = bps.merge(feat_df, on=['event_id', 'side'])
    bpdf['cluster_mechanism'] = bpdf.cluster.map(mech_by_cluster)
    bpdf.to_csv(os.path.join(args.outdir, 'breakpoints.tsv'), sep='\t', index=False)

    mh.to_csv(os.path.join(args.outdir, 'stage7_junction_mh.tsv'), sep='\t',
              index=False)
    ev_out = ev.merge(
        mh.drop(columns=['grp', 'chrom1', 'pos1', 'chrom2', 'pos2', 'key']),
        on='event_id', how='left')
    ev_out.to_csv(os.path.join(args.outdir, 'events_mechanism.tsv'), sep='\t',
                  index=False)

    # ---- summary ---------------------------------------------------------
    n_null = mh.null_mh_max.notna().sum()
    null_ge4 = int((mh.null_mh_max.dropna() >= 4).sum())
    obs_ge4 = int((mh.mh_called >= 4).sum())
    if n_null:
        _, p_mh, _, _ = fisher_or(obs_ge4, len(mh), null_ge4, int(n_null))
    else:
        p_mh = float('nan')
    with open(os.path.join(args.outdir, 'summary.txt'), 'w') as fh:
        fh.write(f'events {len(ev)}\nbreakpoints {len(bps)}\n'
                 f'clusters {len(reps)}\n'
                 f'bg windows matched {len(bgm_seqs)} uniform {len(bgu_seqs)}\n')
        fh.write('\nmechanism counts:\n')
        fh.write(mh.mechanism.value_counts().to_string() + '\n')
        fh.write('\nmechanism x group:\n')
        fh.write(pd.crosstab(mh.mechanism, mh.grp).to_string() + '\n')
        fh.write(f'\njunction MH >=4: obs {obs_ge4}/{len(mh)} '
                 f'null {null_ge4}/{int(n_null)} fisher p={p_mh:.3g}\n')
        fh.write(f'join forms: {mh.join_form.value_counts().to_dict()}\n')
        fh.write(f'structure-tagged events: {int(mh.structure_tag.sum())}\n')
        fh.write(f'telomere-repeat-tagged events: {int(mh.telo_tag.sum())}\n')
        fh.write(f'L1-EN-tagged events: {int(mh.l1en_tag.sum())}\n')
    print('[done] outputs in', args.outdir)


if __name__ == '__main__':
    main()

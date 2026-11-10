#!/usr/bin/env python3
"""Microhomology at eccDNA junction breakpoints.

Sequence-level repeat check at breakpoints: for every decoil junction the two
joined breakpoint ends each contribute the --flank bases immediately inside
their fragment (in molecule orientation, upstream fragment read 3'-ward,
downstream fragment 5'-ward).  The microhomology length is the largest k such
that the last k bases of the upstream flank equal the first k bases of the
downstream flank — i.e. a short direct repeat straddling the junction whose
placement is ambiguous (the MMEJ/alt-EJ footprint).

Statistics use junctions collapsed across samples (both breakpoint sides
within --tolerance bp); the null model redraws both breakpoints at random
positions on the same chromosomes with the same flank orientation and
recomputes the same statistic (--n-draws).

Outputs: per-junction TSV, threshold summary TSV, stdout summary.

Parameters
----------
--ecdna-bed : ecDNA fragment BED as ``SAMPLE=PATH`` (repeatable; the file may
    live anywhere; typically ``reconstruct.ecDNA.filtered.bed``)
--ecdna-links : junction links TSV as ``SAMPLE=PATH`` (repeatable; any file
    location; typically ``reconstruct.links.ecDNA.filtered.txt``)
--genome-fasta : reference FASTA with .fai index
--flank : flank length in bp (default 50)
--tolerance : junction collapse distance in bp (default 200)
--n-draws : null draws (default 1000)
--seed : RNG seed (default 0)
--out-dir : output directory
"""
from __future__ import annotations

import argparse
import logging
import os
import random
import statistics
import sys
from dataclasses import dataclass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repeat_utils import (FastaReader, Junction, empirical_p,  # noqa: E402
                          load_fragments, load_junctions, revcomp,
                          sample_map, sample_path_spec)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("ecdna_junction_microhomology")

THRESHOLDS = (1, 2, 3, 5, 8, 10)

# both junction breakpoints bucketed to --tolerance resolution
JunctionKey = tuple[tuple[str, int], tuple[str, int]]


@dataclass
class JuncRecord:
    """One junction's microhomology call."""
    sample: str
    circ: str
    up_chrom: str
    up_pos: int
    down_chrom: str
    down_pos: int
    mh_len: int
    flank_up: str
    flank_down: str
    key: JunctionKey


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--ecdna-bed", action="append", required=True,
                        type=sample_path_spec, metavar="SAMPLE=PATH",
                        help="ecDNA fragment BED (repeatable; any file location)")
    parser.add_argument("--ecdna-links", action="append", required=True,
                        type=sample_path_spec, metavar="SAMPLE=PATH",
                        help="Junction links TSV (repeatable; any file location)")
    parser.add_argument("--genome-fasta", required=True,
                        help="Reference FASTA with .fai index")
    parser.add_argument("--flank", type=int, default=50,
                        help="Flank length in bp (default: 50)")
    parser.add_argument("--tolerance", type=int, default=200,
                        help="Junction collapse distance (default: 200)")
    parser.add_argument("--n-draws", type=int, default=1000,
                        help="Null draws (default: 1000)")
    parser.add_argument("--seed", type=int, default=0,
                        help="RNG seed (default: 0)")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    args = parser.parse_args(argv)
    try:
        args.ecdna_bed = sample_map(args.ecdna_bed)
        args.ecdna_links = sample_map(args.ecdna_links)
    except ValueError as exc:
        parser.error(str(exc))
    return args


def junction_geometry(junc: Junction) -> tuple[str, str]:
    """Role of each joined end: ``('3p', '5p')`` up/downstream at junction.

    side 1 = end-coordinate side, side 0 = start-coordinate side; for a ``+``
    fragment the molecule runs start->end so side 1 is its 3' terminus and
    side 0 its 5' terminus (flipped for ``-`` strand).
    """
    def terminus(frag, side: int) -> str:
        plus = frag.strand == "+"
        return "3p" if (side == 1) == plus else "5p"
    return terminus(junc.frag1, junc.side1), terminus(junc.frag2, junc.side2)


def flank_3p(frag, fa: FastaReader, w: int) -> str:
    """Last ``w`` bases of the fragment in molecule orientation."""
    if frag.strand == "+":
        return fa.fetch(frag.chrom, frag.end - w, frag.end)
    return revcomp(fa.fetch(frag.chrom, frag.start, frag.start + w))


def flank_5p(frag, fa: FastaReader, w: int) -> str:
    """First ``w`` bases of the fragment in molecule orientation."""
    if frag.strand == "+":
        return fa.fetch(frag.chrom, frag.start, frag.start + w)
    return revcomp(fa.fetch(frag.chrom, frag.end - w, frag.end))


def microhomology_len(flank_a: str, flank_b: str) -> int:
    """Largest k with ``flank_a[-k:] == flank_b[:k]`` (0 when none)."""
    max_k = min(len(flank_a), len(flank_b))
    for k in range(max_k, 0, -1):
        if flank_a[-k:] == flank_b[:k]:
            return k
    return 0


def junction_key(junc: Junction, tolerance: int) -> JunctionKey:
    """Collapse key: both breakpoint sides rounded to ``tolerance`` buckets."""
    g1, g2 = junction_geometry(junc)
    b1 = (junc.frag1.chrom,
          (junc.frag1.oriented_end_coord if g1 == "3p"
           else junc.frag1.oriented_start_coord))
    b2 = (junc.frag2.chrom,
          (junc.frag2.oriented_start_coord if g2 == "5p"
           else junc.frag2.oriented_end_coord))
    pair = tuple(sorted([(c, p // tolerance) for c, p in (b1, b2)]))
    return pair  # type: ignore[return-value]


def random_flank(rng: random.Random, chrom: str, sizes: dict[str, int],
                 fa: FastaReader, w: int, role: str) -> str:
    """Random flank with the same orientation convention as ``role``."""
    size = sizes[chrom]
    if role == "3p":  # sequence ending at the drawn breakpoint
        p = rng.randrange(w, size)
        return fa.fetch(chrom, p - w, p)
    p = rng.randrange(0, size - w)
    return fa.fetch(chrom, p, p + w)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)
    fa = FastaReader(args.genome_fasta)
    sizes = fa.chrom_sizes()

    fragments = load_fragments(args.ecdna_bed)
    only_bed = set(fragments) - set(args.ecdna_links)
    if only_bed:
        logger.warning(f"no links for samples {sorted(only_bed)}; "
                       f"their junctions are skipped")
    junctions = load_junctions(args.ecdna_links, fragments)
    logger.info(f"{len(junctions)} junctions from "
                f"{len(args.ecdna_links)} samples")

    # ---------------- observed microhomology ---------------- #
    records: list[JuncRecord] = []
    seen: dict[JunctionKey, JuncRecord] = {}
    for junc in junctions:
        t1, t2 = junction_geometry(junc)
        if t1 == "3p" and t2 == "5p":
            up, down = junc.frag1, junc.frag2
            up_t, down_t = t1, t2
        elif t1 == "5p" and t2 == "3p":
            up, down = junc.frag2, junc.frag1
            up_t, down_t = t2, t1
        else:
            logger.warning(f"{junc.sample} circ{junc.circ_id}: ambiguous "
                           f"junction geometry ({t1}, {t2}); treating n1 as "
                           f"upstream")
            up, down = junc.frag1, junc.frag2
            up_t, down_t = t1, t2
        flank_up = flank_3p(up, fa, args.flank) if up_t == "3p" \
            else flank_5p(up, fa, args.flank)
        flank_down = flank_5p(down, fa, args.flank) if down_t == "5p" \
            else flank_3p(down, fa, args.flank)
        mh = microhomology_len(flank_up, flank_down)
        rec = JuncRecord(sample=junc.sample, circ=junc.circ_id,
                         up_chrom=up.chrom, up_pos=up.oriented_end_coord,
                         down_chrom=down.chrom,
                         down_pos=down.oriented_start_coord,
                         mh_len=mh, flank_up=flank_up, flank_down=flank_down,
                         key=junction_key(junc, args.tolerance))
        records.append(rec)
        if rec.key not in seen or mh > seen[rec.key].mh_len:
            seen[rec.key] = rec
    uniq = list(seen.values())
    logger.info(f"{len(uniq)} unique junctions after "
                f"{args.tolerance} bp collapse")

    junc_path = os.path.join(args.out_dir, "junction_microhomology.tsv")
    with open(junc_path, "w") as fh:
        fh.write("sample\tcirc\tup_chrom\tup_pos\tdown_chrom\tdown_pos\t"
                 "mh_len\tup_flank_tail\tdown_flank_head\tunique_representative\n")
        for rec in records:
            rep = "yes" if seen[rec.key] is rec else "no"
            fh.write(f"{rec.sample}\t{rec.circ}\t{rec.up_chrom}\t"
                     f"{rec.up_pos}\t{rec.down_chrom}\t{rec.down_pos}\t"
                     f"{rec.mh_len}\t{rec.flank_up[-12:]}\t"
                     f"{rec.flank_down[:12]}\t{rep}\n")
    logger.info(f"Wrote {junc_path}")

    # ---------------- null draws ---------------- #
    obs_mh = [r.mh_len for r in uniq]
    null_counts = {t: [] for t in THRESHOLDS}
    null_mean = []
    for draw in range(args.n_draws):
        mh_null = []
        for rec in uniq:
            fa_up = random_flank(rng, rec.up_chrom, sizes, fa,
                                 args.flank, "3p")
            fa_down = random_flank(rng, rec.down_chrom, sizes, fa,
                                   args.flank, "5p")
            mh_null.append(microhomology_len(fa_up, fa_down))
        null_mean.append(statistics.mean(mh_null))
        for t in THRESHOLDS:
            null_counts[t].append(sum(1 for m in mh_null if m >= t))
        if (draw + 1) % 200 == 0:
            logger.info(f"null draws: {draw + 1}/{args.n_draws}")

    sum_path = os.path.join(args.out_dir, "microhomology_summary.tsv")
    with open(sum_path, "w") as fh:
        fh.write("threshold\tobs_n\tnull_mean_n\tenrichment\tp_emp\n")
        for t in THRESHOLDS:
            obs_n = sum(1 for m in obs_mh if m >= t)
            mean_n = statistics.mean(null_counts[t])
            enr = obs_n / mean_n if mean_n > 0 else float("inf")
            p = empirical_p(float(obs_n),
                            [float(c) for c in null_counts[t]])
            fh.write(f"mh>={t}\t{obs_n}\t{mean_n:.2f}\t{enr:.2f}\t{p:.4f}\n")
    logger.info(f"Wrote {sum_path}")

    # ---------------- stdout summary ---------------- #
    n = len(uniq)
    print(f"\n=== eccDNA junction microhomology "
          f"({n} unique junctions, flank {args.flank} bp) ===")
    for t in THRESHOLDS:
        obs_n = sum(1 for m in obs_mh if m >= t)
        mean_n = statistics.mean(null_counts[t])
        p = empirical_p(float(obs_n), [float(c) for c in null_counts[t]])
        print(f"  mh>={t:<2d} obs {obs_n}/{n}  null mean {mean_n:.1f}  "
              f"p={p:.4f}")
    print(f"mean microhomology: obs {statistics.mean(obs_mh):.2f} bp, "
          f"null {statistics.mean(null_mean):.2f} bp")
    return 0


if __name__ == "__main__":
    sys.exit(main())

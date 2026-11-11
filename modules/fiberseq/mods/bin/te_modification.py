#!/usr/bin/env python3
"""Per-locus TE modification quantification with family-level collapse.

Two-stage engine for TE (transposable element) methylation analysis:

Stage 1 (locus level): for every TE locus in a RepeatMasker-style BED and
every sample, computes modification calls per covering fiber (m6A or CpG,
same call-block logic as extract_modification_longtable). A locus row is
kept per sample only when >= min_fibers fibers cover it.

Stage 2 (collapse): groups loci by subfamily (BED col4) or family (col7,
configurable) and computes the coverage-weighted mean per group x sample —
loci weighted by their covering-fiber count, uncovered loci excluded from
that sample's mean. Output is the sigtest long table (feature x sample x
group x value) at the FAMILY/SUBFAMILY level, plus a per-locus value
table (all samples, annotation columns kept) consumed downstream to rank
top loci within significant families.

BED format (rheMac10_rmsk_TE.lifted.bed convention):
    chrom start end subfamily score strand family class
col4 = subfamily (AluY / L1PA2 / MIRb ...), col7 = family (Alu / L1 ...),
col8 = class (SINE / LINE / LTR / DNA). Input MUST be sorted
`sort -k1,1 -k2,2n` (ChromSweep contract).

Usage:
    python te_modification.py \
        --mods m6a:s1.m6a.bed.gz,m6a:s2.m6a.bed.gz \
        --msp s1.msp.bed.gz,s2.msp.bed.gz \
        --te te.sorted.bed --design design.tsv \
        --te-level subfamily --min-fibers 3 \
        --out-long te_m6a_subfamily_long.tsv \
        --out-loci te_m6a_loci.tsv
"""
from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO, Tuple

# Self-contained module bin scripts import the module-private helper
# (verbatim src/common copies) — modules.md 脚本撰写规范 #0.
from _fiberlib import ChromSweep, setup_logger

logger = setup_logger("TeModification", stream=sys.stdout)


def parse_group_samples(spec: str) -> Dict[str, str]:
    """Parse 'gA:s1+s2,gB:s3' (node.py group_samples) into sample -> group.

    Raises ValueError naming the malformed token so a config typo surfaces
    at DAG/run time with the offending element.
    """
    sample_group: Dict[str, str] = {}
    for tok in spec.split(","):
        if ":" not in tok:
            raise ValueError(f"malformed group spec token {tok!r} "
                             "(expected group:sample[+sample])")
        group, s_str = tok.split(":", 1)
        for s in s_str.split("+"):
            if s.strip():
                sample_group[s.strip()] = group.strip()
    return sample_group


class TeLocus:
    """One TE locus row with its subfamily/family/class annotation."""

    __slots__ = ("chrom", "start", "end", "subfamily", "family", "cls",
                 "calls", "fibers")

    def __init__(self, chrom: str, start: int, end: int,
                 subfamily: str, family: str, cls: str) -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.subfamily = subfamily
        self.family = family
        self.cls = cls
        self.calls: int = 0
        self.fibers: int = 0

    @property
    def key(self) -> str:
        """Stable locus id: coordinates + subfamily (dedup key)."""
        return f"{self.chrom}:{self.start}-{self.end}|{self.subfamily}"


def read_te_bed(path: str) -> List[TeLocus]:
    """Parse the TE BED (col4 subfamily, col7 family, col8 class).

    Comment (#) / track / browser lines are skipped. Missing annotation
    columns degrade to "." so malformed rows never abort a 4M-row file.
    """
    loci: List[TeLocus] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            sub = f[3].strip() if len(f) > 3 and f[3].strip() else "."
            fam = f[6].strip() if len(f) > 6 and f[6].strip() else "."
            cls = f[7].strip() if len(f) > 7 and f[7].strip() else "."
            loci.append(TeLocus(f[0], int(f[1]), int(f[2]), sub, fam, cls))
    return loci


def quantify_sample(loci: List[TeLocus], mod_bed: str, msp_bed: str) -> None:
    """Fill calls/fibers on every locus for one sample (two sweeps).

    Sweep 1 (modification BED12): counts call blocks overlapping each
    locus interval. Sweep 2 (MSP): distinct fiber names covering the
    locus — the quant_matrix n_fibers convention. Loci absent from a
    chromosome keep zeros (reported downstream as uncovered).
    """
    anchors = [(t.chrom, t.start, t.end) for t in loci]
    sweep = ChromSweep([mod_bed], anchors)
    for k, (_c, _s, e), ov in sweep.advance():
        n = 0
        t = loci[k]
        for rec in ov[0]:
            for b_start, b_end in rec.iter_blocks():
                if b_start < t.end and b_end > t.start:
                    n += 1
        t.calls = n
    sweep = ChromSweep([msp_bed], anchors)
    for k, (_c, _s, _e), ov in sweep.advance():
        loci[k].fibers = len({r.name for r in ov[0]})


def collapse_to_level(
    loci: List[TeLocus],
    per_sample: Dict[str, List[TeLocus]],
    sample_ids: Sequence[str],
    min_fibers: int,
    level: str,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Tuple[str, str]], List[str]]:
    """Coverage-weighted mean per te-level group x sample.

    Loci with fibers < min_fibers in a sample are excluded from THAT
    sample's mean (they carry no reliable signal), not zero-filled —
    a group's n naturally reflects usable loci. Returns (values, level ->
    (family, class) annotation, level order).
    """
    attr = "subfamily" if level == "subfamily" else "family"
    sums: Dict[str, Dict[str, float]] = {}
    weights: Dict[str, Dict[str, float]] = {}
    annot: Dict[str, Tuple[str, str]] = {}
    order: List[str] = []
    for sid in sample_ids:
        for t in per_sample[sid]:
            if t.fibers < min_fibers:
                continue
            name = getattr(t, attr)
            if name == ".":
                continue
            if name not in sums:
                sums[name] = {}
                weights[name] = {}
                annot[name] = (t.family, t.cls)
                order.append(name)
            val = t.calls / t.fibers
            sums[name][sid] = sums[name].get(sid, 0.0) + val * t.fibers
            weights[name][sid] = weights[name].get(sid, 0.0) + t.fibers
    values: Dict[str, Dict[str, float]] = {}
    for name in order:
        values[name] = {}
        for sid in sample_ids:
            w = weights[name].get(sid, 0.0)
            if w > 0:
                values[name][sid] = sums[name][sid] / w
    return values, annot, order


def sort_te_bed(te_path: str, tmp_dir: str) -> str:
    """Sort a TE BED `k1,1 -k2,2n` into tmp scratch (ChromSweep contract).

    Returns the original path when it is already sorted (cheap adjacent
    check); otherwise an external `LC_ALL=C sort` writes the sorted copy
    under <tmp_dir>/te_modification/ so multi-GB inputs never load into
    memory. The sorted copy is reused across runs when it exists and is
    newer than the source.
    """
    last_chrom, last_start, unsorted = None, -1, 0
    with open(te_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t", 3)
            if len(f) < 3:
                continue
            start = int(f[1])
            if f[0] == last_chrom and start < last_start:
                unsorted += 1
                if unsorted > 100:
                    break
            last_chrom, last_start = f[0], start
    if unsorted == 0:
        return te_path
    out_dir = Path(tmp_dir) / "te_modification"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / (Path(te_path).name + ".sorted.bed")
    src_mtime = Path(te_path).stat().st_mtime
    if not (out_path.is_file() and out_path.stat().st_mtime >= src_mtime):
        logger.info(f"TE bed not coordinate-sorted; sorting -> {out_path}")
        subprocess.run(
            ["bash", "-c",
             f"LC_ALL=C sort -k1,1 -k2,2n {shlex.quote(te_path)} > {shlex.quote(str(out_path))}"],
            check=True,
        )
    return str(out_path)


def run(args: argparse.Namespace) -> int:
    """Stage 0 (sort) + stage 1 (per-locus quantify) + stage 2 (collapse)."""
    entries = [e for e in args.mods.split(",") if e.strip()]
    sample_ids = [e.split(":")[0] for e in entries]
    mod_beds = [e.split(":", 1)[1] for e in entries]
    msp_beds = [m for m in args.msp.split(",") if m.strip()]
    if len(msp_beds) != len(sample_ids):
        logger.error("--msp must list one bed per sample, same order as --mods")
        return 1
    sample_group = parse_group_samples(args.group_samples)
    missing = [s for s in sample_ids if s not in sample_group]
    if missing:
        logger.error(f"samples absent from design TSV: {missing}")
        return 1

    te_path = sort_te_bed(args.te, args.tmp_dir)
    loci = read_te_bed(te_path)
    logger.info(
        f"TE bed: {len(loci)} loci; level={args.te_level} "
        f"min_fibers={args.min_fibers}; samples={sample_ids}"
    )
    # per-locus per-sample tables (fresh copies: values are sample-specific)
    per_sample: Dict[str, List[TeLocus]] = {}
    for sid, mod_bed, msp_bed in zip(sample_ids, mod_beds, msp_beds):
        sample_loci = read_te_bed(args.te)
        quantify_sample(sample_loci, mod_bed, msp_bed)
        per_sample[sid] = sample_loci
        n_cov = sum(1 for t in sample_loci if t.fibers >= args.min_fibers)
        logger.info(f"sample {sid}: {n_cov}/{len(sample_loci)} loci covered "
                    f"(>= {args.min_fibers} fibers)")

    values, annot, order = collapse_to_level(
        loci, per_sample, sample_ids, args.min_fibers, args.te_level
    )
    # family-level long table for sigtest
    out = Path(args.out_long)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("feature\tsample\tgroup\tvalue\n")
        for name in order:
            for sid in sample_ids:
                if sid in values[name]:
                    fh.write(
                        f"{name}\t{sid}\t{sample_group[sid]}\t"
                        f"{values[name][sid]:.6f}\n"
                    )
    logger.info(f"Wrote family-level long table: {out} ({len(order)} {args.te_level}s)")

    # per-locus table (all samples wide + annotation) for top-loci ranking
    loci_out = Path(args.out_loci)
    with open(loci_out, "w") as fh:
        fh.write("chrom\tstart\tend\tsubfamily\tfamily\tclass\t"
                 + "\t".join(f"{s}.m6a_per_fiber" for s in sample_ids)
                 + "\t" + "\t".join(f"{s}.fibers" for s in sample_ids) + "\n")
        for i, t in enumerate(loci):
            vals = []
            for sid in sample_ids:
                st = per_sample[sid][i]
                vals.append(f"{st.calls / st.fibers:.6f}" if st.fibers else "NA")
            fibs = [str(per_sample[sid][i].fibers) for sid in sample_ids]
            fh.write("\t".join([t.chrom, str(t.start), str(t.end),
                                t.subfamily, t.family, t.cls]
                               + vals + fibs) + "\n")
    logger.info(f"Wrote per-locus table: {loci_out}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="te_modification.py",
        description="Per-locus TE modification quantification + family-level collapse.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--mods", required=True,
                   help="comma list of sample:mod_bed.bed.gz (m6a or cpg)")
    p.add_argument("--msp", required=True,
                   help="comma list of per-sample msp.bed.gz (fiber denominator)")
    p.add_argument("--te", required=True,
                   help="TE BED sorted k1,1 -k2,2n (col4 subfamily, col7 family, col8 class)")
    p.add_argument("--group-samples", required=True,
                   help="group spec 'gA:s1+s2,gB:s3' (sample -> group mapping)")
    p.add_argument("--te-level", choices=["subfamily", "family"],
                   default="subfamily",
                   help="collapse/test level (subfamily = col4, family = col7)")
    p.add_argument("--min-fibers", type=int, default=3,
                   help="per-sample locus coverage gate for collapsing")
    p.add_argument("--out-long", required=True,
                   help="family-level sigtest long TSV")
    p.add_argument("--out-loci", required=True,
                   help="per-locus value TSV (top-loci ranking input)")
    p.add_argument("--tmp-dir", default="/tmp",
                   help="scratch dir for the sorted TE copy when input is unsorted")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    if args.log:
        logger = setup_logger("TeModification", log_file=str(args.log),
                              stream=sys.stdout)
    logger.info(f"Command: {' '.join(sys.argv)}")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())

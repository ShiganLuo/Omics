#!/usr/bin/env python3
"""Rank top loci within significant TE subfamilies/families.

Consumes the family-level sigtest result (significance_test.tsv whose
features are subfamily or family names) plus the per-locus value table
from te_modification.py, and emits per significant group the loci with
the largest between-group mean difference — a BROWSING shortlist for
regionplot, NOT locus-level statistics (family-level p-values are the
only statistical assertion in this design).

Ranking rule per locus: |mean(groupA) - mean(groupB)| over the two
design groups, computed from per-locus per-sample values; loci below
min_fibers in ANY sample of both groups are excluded (reliability gate
mirrors the collapse step).

Usage:
    python te_top_loci.py \
        --sig significance_test.tsv --loci te_m6a_loci.tsv \
        --design design.tsv --groups aged,youth \
        --top-n 5 --min-fibers 3 \
        --out-tsv top_loci.tsv --out-bed top_loci.bed
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO, Tuple

# Copied verbatim from src/common/util/LogUtil.py (module bin scripts are
# self-contained — see modules.md 脚本撰写规范 #0).


def setup_logger(
    logger_name: str,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    stream: Optional[TextIO] = None,
) -> logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


logger = setup_logger("TeTopLoci", stream=sys.stdout)


def read_significant(sig_tsv: str) -> List[str]:
    """Feature names (subfamily/family) flagged significant=True."""
    names: List[str] = []
    with open(sig_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        low = [h.strip().lower() for h in header]
        f_idx, s_idx = low.index("feature"), low.index("significant")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > s_idx and f[s_idx].strip().lower() in ("true", "1", "yes"):
                names.append(f[f_idx])
    return names


def read_loci(loci_tsv: str, groups: Dict[str, List[str]],
              value_suffix: str) -> Tuple[List[dict], List[str]]:
    """Parse the per-locus table; returns loci rows and sample column order.

    Each row carries per-sample value/fiber fields keyed by sample id
    (columns <sample>.m6a_per_fiber / <sample>.fibers from te_modification).
    """
    with open(loci_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        val_col = {h[: -len("." + value_suffix)]: i for i, h in enumerate(header)
                   if h.endswith("." + value_suffix)}
        fib_col = {h[: -len(".fibers")]: i for i, h in enumerate(header)
                   if h.endswith(".fibers")}
        rows: List[dict] = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            rows.append({
                "chrom": f[0], "start": f[1], "end": f[2],
                "subfamily": f[3], "family": f[4], "class": f[5],
                "vals": {s: f[i] for s, i in val_col.items()},
                "fibers": {s: f[i] for s, i in fib_col.items()},
            })
    return rows, list(val_col)


def locus_diff(row: dict, group_samples: Dict[str, List[str]],
               min_fibers: int) -> Optional[float]:
    """|mean(groupA) - mean(groupB)| for one locus; None when excluded.

    Excluded when either group has no sample whose fiber count clears
    min_fibers with a numeric value (reliability gate mirroring the
    collapse step in te_modification.py).
    """
    means: List[float] = []
    for samples in group_samples.values():
        vals = []
        for s in samples:
            try:
                if int(row["fibers"].get(s, "0")) >= min_fibers:
                    v = float(row["vals"][s])
                    vals.append(v)
            except (ValueError, KeyError):
                continue
        if not vals:
            return None
        means.append(sum(vals) / len(vals))
    return abs(means[0] - means[1])


def run(args: argparse.Namespace) -> int:
    """Write per-significant-group top loci as TSV + regionplot BED."""
    groups = {}
    for pair in args.groups.split(","):
        g, s_str = pair.split(":", 1)
        groups[g] = [s.strip() for s in s_str.split("+") if s.strip()]
    if len(groups) != 2:
        logger.error(f"--groups needs exactly two group:sample[+sample] entries, got {groups}")
        return 1

    sig_names = set(read_significant(args.sig))
    logger.info(f"significant groups from {args.sig}: {len(sig_names)}")

    rows, _samples = read_loci(args.loci, groups, args.value_suffix)
    # bucket loci by BOTH subfamily and family so either tested level matches
    by_level: Dict[str, List[dict]] = {}
    for row in rows:
        for key in (row["subfamily"], row["family"]):
            if key in sig_names:
                by_level.setdefault(key, []).append(row)

    out_tsv, out_bed = Path(args.out_tsv), Path(args.out_bed)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    n_total = 0
    with open(out_tsv, "w") as tfh, open(out_bed, "w") as bfh:
        tfh.write("group\tchrom\tstart\tend\tsubfamily\tfamily\tclass\t"
                  "abs_group_mean_diff\tn_samples_used\n")
        bfh.write("#chrom\tstart\tend\tname\tscore\tstrand\tgroup\tabs_diff\n")
        for name in sorted(sig_names):
            cands = by_level.get(name, [])
            scored = []
            for row in cands:
                d = locus_diff(row, groups, args.min_fibers)
                if d is not None:
                    scored.append((d, row))
            scored.sort(key=lambda x: (-x[0], x[1]["chrom"], int(x[1]["start"])))
            for d, row in scored[: args.top_n]:
                locus_name = f"{row['subfamily']}@{row['chrom']}_{row['start']}_{row['end']}"
                tfh.write("\t".join([name, row["chrom"], row["start"], row["end"],
                                     row["subfamily"], row["family"], row["class"],
                                     f"{d:.6f}", str(len(groups))]) + "\n")
                bfh.write("\t".join([row["chrom"], row["start"], row["end"],
                                     locus_name, f"{d:.6f}", ".", name,
                                     f"{d:.6f}"]) + "\n")
                n_total += 1
    logger.info(f"Wrote {n_total} top loci across {len(sig_names)} significant "
                f"groups (top {args.top_n} each) -> {out_tsv}, {out_bed}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="te_top_loci.py",
        description="Top loci within significant TE subfamilies (regionplot shortlist).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sig", required=True, help="family-level significance_test.tsv")
    p.add_argument("--loci", required=True, help="per-locus table from te_modification.py")
    p.add_argument("--groups", required=True,
                   help="two group specs 'gA:s1+s2,gB:s3+s4' (design groups -> samples)")
    p.add_argument("--value-suffix", default="m6a_per_fiber",
                   help="value column suffix in the loci table")
    p.add_argument("--top-n", type=int, default=5,
                   help="loci per significant group (largest |group mean diff|)")
    p.add_argument("--min-fibers", type=int, default=3,
                   help="per-sample locus coverage gate for ranking")
    p.add_argument("--out-tsv", required=True, help="ranked loci TSV")
    p.add_argument("--out-bed", required=True, help="regionplot-ready BED")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    if args.log:
        logger = setup_logger("TeTopLoci", log_file=str(args.log), stream=sys.stdout)
    logger.info(f"Command: {' '.join(sys.argv)}")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())

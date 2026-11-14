#!/usr/bin/env python3
"""Split sigtest results into significant vs non-significant object files.

Reads one ``significance_test.tsv`` (sigtest.py output), writes the
significant subset as BOTH a TSV (with stats) and a BED (chrom/start/end/
name/score=p), ready to feed regionplot's Params.regionplot.regions or
any interval tool. Non-significant and skipped features go to a separate
TSV so the split is auditable (row counts add up).

Feature ids follow the mods long-table convention
``<chrom>:<start>-<end>|<name>`` (quant_matrix region_key); anything that
does not parse is reported and dropped from the BED (kept in the TSV).

Usage:
    python split_significant.py -i significance_test.tsv \
        --out-tsv significant.tsv --out-bed significant.bed \
        [--min-abs-effect 0] [--alpha 0.05]
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import List, Optional, Sequence, TextIO

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


logger = setup_logger("SplitSignificant", stream=sys.stdout)

# feature id: <chrom>:<start>-<end>|<name>
_REGION_RE = re.compile(r"^(\S+):(\d+)-(\d+)\|(.+)$")


def split_results(args: argparse.Namespace) -> int:
    """Write significant TSV+BED (BED optionally top-N by p) and the rest."""
    in_path = Path(args.input)
    if not in_path.is_file():
        logger.error(f"input not found: {in_path}")
        return 1

    header: List[str] = []
    sig_rows: List[List[str]] = []
    rest_rows: List[List[str]] = []
    n_bad = 0
    with open(in_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        low = [h.strip().lower() for h in header]
        f_idx = low.index("feature")
        p_idx = low.index("p_adjusted") if "p_adjusted" in low else low.index("p_value")
        sig_idx = low.index("significant")
        eff_idx = low.index("effect_size") if "effect_size" in low else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            is_sig = f[sig_idx].strip().lower() in ("true", "1", "yes")
            if is_sig and args.min_abs_effect > 0 and eff_idx is not None:
                try:
                    if abs(float(f[eff_idx])) < args.min_abs_effect:
                        is_sig = False
                except ValueError:
                    pass
            (sig_rows if is_sig else rest_rows).append(f)

    # BED = top-N by p (ascending); TSV keeps ALL significant rows so the
    # full set stays auditable. top_n<=0 (or > available) keeps everything.
    # rank only rows that yield a BED interval, so top-N means N BED rows
    def _parseable(row: List[str]) -> bool:
        return _REGION_RE.match(row[f_idx]) is not None

    bed_pool = [r for r in sig_rows if _parseable(r)]
    bed_rows = bed_pool
    n_trimmed = 0
    if args.top_n and args.top_n > 0 and len(bed_pool) > args.top_n:
        def _p(row: List[str]) -> float:
            try:
                return float(row[p_idx])
            except (ValueError, IndexError):
                return float("inf")

        ranked = sorted(bed_pool, key=_p)
        bed_rows = ranked[: args.top_n]
        n_trimmed = len(bed_pool) - args.top_n

    for path, rows in ((args.out_tsv, sig_rows), (args.out_rest_tsv, rest_rows)):
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")

    bed = Path(args.out_bed) if args.out_bed else None
    if bed is None:
        # no BED requested (TE axis: features are family names, no coords)
        logger.info(" BED skipped (--out-bed absent)")
    else:
        bed.parent.mkdir(parents=True, exist_ok=True)
        with open(bed, "w") as fh:
            fh.write("#chrom\tstart\tend\tname\tscore\tstrand\tp_adjusted\tstars\n")
            for r in bed_rows:
                m = _REGION_RE.match(r[f_idx])
                if not m:
                    n_bad += 1
                    continue
                chrom, start, end, name = m.groups()
                p = r[p_idx]
                stars = r[low.index("stars")] if "stars" in low else ""
                fh.write(f"{chrom}\t{start}\t{end}\t{name}\t{p}\t.\t{p}\t{stars}\n")

    bed_desc = str(bed) if bed is not None else "no BED"
    logger.info(
        f"significant={len(sig_rows)} non-significant/skipped={len(rest_rows)} "
        f"-> {args.out_tsv}, {bed_desc}"
        + (f"; BED limited to top {len(bed_rows)} by p ({n_trimmed} trimmed)"
           if n_trimmed else "")
        + (f" ({n_bad} unparseable feature ids stayed TSV-only)" if n_bad else "")
    )
    if n_bad:
        logger.warning(f"{n_bad} feature ids did not match chrom:start-end|name")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="split_significant.py",
        description="Split sigtest results into significant (TSV+BED) and rest.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, help="sigtest significance_test.tsv")
    p.add_argument("--out-tsv", required=True, help="significant subset TSV (all columns)")
    p.add_argument("--out-bed", default=None,
                   help="significant subset BED (regionplot-ready)")
    p.add_argument("--out-rest-tsv", required=True, help="non-significant remainder TSV")
    p.add_argument("--min-abs-effect", type=float, default=0.0,
                   help="also require |effect_size| >= this for the BED")
    p.add_argument("--top-n", type=int, default=0,
                   help="BED keeps only the N smallest-p significant regions "
                        "(0 = keep all; the TSV always keeps all)")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    if args.log:
        logger = setup_logger("SplitSignificant", log_file=str(args.log), stream=sys.stdout)
    logger.info(f"Command: {' '.join(sys.argv)}")
    return split_results(args)


if __name__ == "__main__":
    sys.exit(main())

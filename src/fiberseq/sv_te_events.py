#!/usr/bin/env python3
"""TE-associated SV / translocation events from an SV VCF + TE annotation.

Classifies each SV (pbsv VCF: DEL/INS/DUP/INV/BND) by its overlap with TE
annotation and reports:

  - per-event table: SV type, size, breakpoints, overlapping TE subfamily/
    family/class (either breakpoint or the whole span)
  - event counts by SV type x TE class/family (the "number and kind" table)
  - inter-chromosomal (translocation-like) events: BND across chromosomes
    plus same-type pairs of BNDs with reciprocal partners
  - TE-family enrichment: observed TE-overlapping events vs expectation from
    shuffled breakpoints (genome bp matched, --n-perm permutations)

Usage:
  python sv_te_events.py --sv-vcf sample.sv.vcf.gz \
      --te-bed GRCm39_GENCODE_rmsk_TE.lifted.bed --fai genome.fa.fai \
      -o sv_te_out/ --sample DMSO
"""
from __future__ import annotations

import argparse
import gzip
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402
from src.common.util.FiberlibUtil import write_tsv  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("SvTeEvents", log_file=str(log_file))
    else:
        lg = setup_logger("SvTeEvents")
    for h in lg.handlers:
        if type(h) is logging.StreamHandler:
            h.setStream(sys.stdout)
    return lg


logger = _get_logger()


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_svs(path: str) -> List[dict]:
    """Parse pbsv-style VCF into event dicts with both endpoints."""
    events = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            svtype = info.get("SVTYPE", "BND")
            chrom, pos = f[0], int(f[1])
            chrom2, pos2 = chrom, int(info.get("END", pos))
            if svtype == "BND":
                alt = f[4] if len(f) > 4 else ""
                for tok in alt.replace("]", "[").split("["):
                    if ":" in tok:
                        try:
                            c2, p2 = tok.rsplit(":", 1)
                            chrom2, pos2 = c2.strip(), int(p2)
                            break
                        except ValueError:
                            continue
            events.append({
                "chrom": chrom, "pos": pos, "chrom2": chrom2, "pos2": pos2,
                "svtype": svtype,
                "svlen": abs(int(info.get("SVLEN", pos2 - pos)) or pos2 - pos),
                "filter": f[6] if len(f) > 6 else ".",
            })
    return events


def load_te_bed(path: str):
    """TE BED sorted per chrom with parallel arrays for bisect queries."""
    per_chrom: Dict[str, List[Tuple[int, int, str, str, str]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            per_chrom[f[0]].append((s, e,
                                    f[3] if len(f) > 3 else ".",
                                    f[6] if len(f) > 6 else ".",
                                    f[7] if len(f) > 7 else "."))
    for c in per_chrom:
        per_chrom[c].sort()
    return per_chrom


def overlap_tes(per_chrom, chrom: str, start: int, end: int, slop: int):
    """TEs overlapping [start-slop, end+slop) on chrom."""
    hits = []
    tes = per_chrom.get(chrom)
    if not tes:
        return hits
    import bisect
    lo = bisect.bisect_left([t[1] for t in tes], start - slop)
    for i in range(lo, len(tes)):
        s, e, sub, fam, cls = tes[i]
        if s > end + slop:
            break
        if s < end + slop and e > start - slop:
            hits.append((sub, fam, cls))
    return hits


def annotate_events(events, per_chrom, slop: int) -> List[dict]:
    """Attach TE hits at both breakpoints and over the span."""
    rows = []
    for i, ev in enumerate(events):
        bp1 = overlap_tes(per_chrom, ev["chrom"], ev["pos"] - slop, ev["pos"] + slop, 0)
        bp2 = overlap_tes(per_chrom, ev["chrom2"], ev["pos2"] - slop, ev["pos2"] + slop, 0)
        span = (overlap_tes(per_chrom, ev["chrom"], ev["pos"], ev["pos2"], 0)
                if ev["chrom"] == ev["chrom2"] else [])
        subs = {h[0] for h in bp1 + bp2 + span}
        fams = {h[1] for h in bp1 + bp2 + span}
        clss = {h[2] for h in bp1 + bp2 + span}
        rows.append({
            **ev,
            "interchrom": ev["chrom"] != ev["chrom2"],
            "n_bp1_te": len(bp1), "n_bp2_te": len(bp2),
            "bp1_te_subfamilies": ";".join(sorted(subs)) if subs else "",
            "bp1_te_families": ";".join(sorted(fams)) if fams else "",
            "bp1_te_classes": ";".join(sorted(clss)) if clss else "",
            "te_associated": bool(subs),
        })
    return rows


def count_matrix(rows: List[dict], key: str) -> List[dict]:
    """SV type x TE key count matrix."""
    field = {"class": "bp1_te_classes", "family": "bp1_te_families",
             "subfamily": "bp1_te_subfamilies"}[key]
    mat: Dict[Tuple[str, str], int] = defaultdict(int)
    types = set()
    for r in rows:
        t = r["svtype"]
        types.add(t)
        vals = [v for v in r.get(field, "").split(";") if v]
        if vals:
            for v in set(vals):
                mat[(t, v)] += 1
        else:
            mat[(t, "(none)")] += 1
    keys = sorted({k[1] for k in mat})
    out = []
    for t in sorted(types):
        base = {"svtype": t, "total": sum(v for (tt, _), v in mat.items() if tt == t)}
        for k in keys:
            base[k] = mat.get((t, k), 0)
        out.append(base)
    return out


def perm_enrichment(rows: List[dict], fai_path: str, n_perm: int, seed: int,
                    te_bed: str = "") -> List[dict]:
    """Per-TE-class enrichment of breakpoints vs shuffled positions.

    Uses prebuilt per-chrom sorted interval arrays + np.searchsorted for the
    null queries (the observed events were annotated the same way above).
    """
    lengths = []
    with open(fai_path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2 and int(f[1]) > 5_000_000:
                lengths.append((f[0], int(f[1])))
    total_bp = sum(ln for _, ln in lengths)
    len_map = dict(lengths)
    rng = np.random.default_rng(seed)

    # prebuild per-chrom TE arrays once (chrom, starts[], ends[], classes[])
    per_chrom_starts: Dict[str, np.ndarray] = {}
    per_chrom_ends: Dict[str, np.ndarray] = {}
    per_chrom_cls: Dict[str, List[str]] = {}
    tmp = defaultdict(list)
    with open(te_bed) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                tmp[f[0]].append((int(f[1]), int(f[2]), f[7] if len(f) > 7 else "."))
            except ValueError:
                continue
    for c, items in tmp.items():
        items.sort()
        per_chrom_starts[c] = np.array([x[0] for x in items], dtype=np.int64)
        per_chrom_ends[c] = np.array([x[1] for x in items], dtype=np.int64)
        per_chrom_cls[c] = [x[2] for x in items]
    del tmp

    classes = sorted({c for r in rows for c in r.get("bp1_te_classes", "").split(";") if c})
    if not classes:
        return []
    obs = {c: sum(1 for r in rows if c in r.get("bp1_te_classes", "").split(";"))
           for c in classes}

    chrom_names = [c for c, _ in lengths]
    chrom_probs = np.array([ln / total_bp for _, ln in lengths])

    perm_hits = {c: np.zeros(n_perm, dtype=np.int64) for c in classes}
    n_events = len(rows)
    for p_i in range(n_perm):
        chrom_pick = rng.choice(chrom_names, size=n_events, p=chrom_probs)
        pos_pick = rng.integers(0, max(2, int(max(len_map.values()))), size=n_events)
        for c, pos in zip(chrom_pick, pos_pick):
            starts = per_chrom_starts.get(c)
            if starts is None or pos >= len_map.get(c, 1):
                continue
            lo = np.searchsorted(starts, pos + 100, side="right")
            hi = np.searchsorted(starts, pos - 100, side="left")
            for j in range(hi, lo):
                if per_chrom_ends[c][j] > pos - 100 and per_chrom_cls[c][j] in perm_hits:
                    perm_hits[per_chrom_cls[c][j]][p_i] += 1
                    break
    out = []
    for c in classes:
        null = perm_hits[c].astype(float)
        p = float((1 + (null >= obs[c]).sum()) / (1 + n_perm))
        out.append({"te_class": c, "observed_events": obs[c],
                    "perm_mean": float(null.mean()), "perm_p": p,
                    "fold": obs[c] / max(float(null.mean()), 1e-9)})
    return out


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="sv_te_events.py",
        description="TE-associated SV/translocation events from SV VCF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sv-vcf", required=True)
    p.add_argument("--te-bed", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("--slop", type=int, default=100,
                   help="bp slack around a breakpoint for TE overlap")
    p.add_argument("--n-perm", type=int, default=1000,
                   help="permutations for TE-class enrichment")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--min-svlen", type=int, default=50)
    p.add_argument("-l", "--log", type=Path, default=None)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    events = [e for e in load_svs(args.sv_vcf) if e["svlen"] >= args.min_svlen]
    logger.info(f"SVs loaded: {len(events)}")
    per_chrom = load_te_bed(args.te_bed)
    logger.info(f"TE contigs: {len(per_chrom)}")

    rows = annotate_events(events, per_chrom, args.slop)
    header = ["chrom", "pos", "chrom2", "pos2", "svtype", "svlen", "filter",
              "interchrom", "n_bp1_te", "n_bp2_te", "bp1_te_subfamilies",
              "bp1_te_families", "bp1_te_classes", "te_associated"]
    write_tsv(str(args.outdir / f"{args.sample}.sv_te_events.tsv"), header,
              [[r.get(c, "") for c in header] for r in rows])

    for key, fname in (("class", "svtype_x_class"), ("family", "svtype_x_family")):
        mat = count_matrix(rows, key)
        if mat:
            cols = ["svtype", "total"] + sorted({c for r in mat for c in r if c not in ("svtype", "total")})
            write_tsv(str(args.outdir / f"{args.sample}.{fname}.tsv"), cols,
                      [[r.get(c, 0) for c in cols] for r in mat])

    inter = [r for r in rows if r["interchrom"]]
    write_tsv(str(args.outdir / f"{args.sample}.interchrom_events.tsv"),
              header, [[r.get(c, "") for c in header] for r in inter])
    logger.info(f"interchromosomal events: {len(inter)}; "
                f"TE-associated: {sum(1 for r in rows if r['te_associated'])}/{len(rows)}")

    if args.n_perm > 0:
        enr = perm_enrichment(rows, args.fai, args.n_perm, args.seed, args.te_bed)
        if enr:
            write_tsv(str(args.outdir / f"{args.sample}.te_class_enrichment.tsv"),
                      ["te_class", "observed_events", "perm_mean", "perm_p", "fold"],
                      [[e[c] for c in ("te_class", "observed_events",
                                       "perm_mean", "perm_p", "fold")] for e in enr])
    logger.info(f"done -> {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

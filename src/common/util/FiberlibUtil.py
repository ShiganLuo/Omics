#!/usr/bin/env python3
"""Shared helpers for Fiber-seq downstream analysis scripts.

Provides BED/BED12 parsing (gzip-aware), per-chrom interval indexing,
windowed per-base signal aggregation, and small stats/plot utilities used
by the analysis scripts in this directory.

Element semantics follow `ft extract` BED12 output: one row per fiber,
name column is the fiber (read) name, and each BED12 block is one element
(m6A call, nucleosome, or MSP). Plain BED3-BED6 files are also accepted and
treated as a single block spanning the full interval.
"""
from __future__ import annotations

import gzip
import os
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]


@dataclass
class BedRecord:
    """One BED/BED12 row with its element blocks as absolute intervals."""

    chrom: str
    start: int
    end: int
    name: str = "."
    score: float = 0.0
    strand: str = "."
    blocks: List[Tuple[int, int]] = field(default_factory=list)

    def iter_blocks(self) -> Iterator[Tuple[int, int]]:
        for b in self.blocks:
            yield b


def open_maybe_gzip(path: str):
    """Open plain or gzip-compressed text file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_bed(path: str) -> List[BedRecord]:
    """Parse a BED3-BED12 file (optionally gzipped) into BedRecord list.

    Parameters
    ----------
    path : str
        Input BED path (plain or .gz). Comment/track lines are skipped.

    Returns
    -------
    list of BedRecord
        Each record carries absolute block intervals; BED3-BED6 rows get a
        single block spanning the full interval.
    """
    records: List[BedRecord] = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue  # skip header lines (e.g. TSV headers in .bed files)
            chrom = f[0]
            name = f[3] if len(f) > 3 else "."
            score = float(f[4]) if len(f) > 4 and f[4] not in (".", "") else 0.0
            strand = f[5] if len(f) > 5 and f[5] else "."
            blocks: List[Tuple[int, int]] = []
            if len(f) >= 12:
                n = int(f[9])
                sizes = [int(x) for x in f[10].rstrip(",").split(",") if x != ""]
                starts = [int(x) for x in f[11].rstrip(",").split(",") if x != ""]
                for i in range(min(n, len(sizes), len(starts))):
                    bs = start + starts[i]
                    blocks.append((bs, bs + sizes[i]))
            if not blocks:
                blocks = [(start, end)]
            records.append(BedRecord(chrom, start, end, name, score, strand, blocks))
    return records


class RegionIndex:
    """Per-chromosome interval index for overlap queries.

    Intervals are indexed by their first block; `query` returns indices of
    records whose indexed interval overlaps [start, end).
    """

    def __init__(self, records: Sequence[BedRecord]) -> None:
        self.by_chrom: Dict[str, Tuple[List[int], List[int], List[int]]] = {}
        tmp: Dict[str, List[Tuple[int, int, int]]] = {}
        for i, r in enumerate(records):
            tmp.setdefault(r.chrom, []).append((r.start, r.end, i))
        for chrom, items in tmp.items():
            items.sort()
            self.by_chrom[chrom] = (
                [x[0] for x in items],
                [x[1] for x in items],
                [x[2] for x in items],
            )

    def query(self, chrom: str, start: int, end: int) -> List[int]:
        """Return record indices on *chrom* with interval overlap [start, end)."""
        if chrom not in self.by_chrom:
            return []
        starts, ends, idxs = self.by_chrom[chrom]
        hi = bisect_left(starts, end)
        out = []
        for i in range(hi):
            if ends[i] > start:
                out.append(idxs[i])
        return out


def _merged_intervals(
    records: Sequence[BedRecord], chrom: str, win_start: int, win_end: int
) -> List[Tuple[int, int]]:
    """Merge all blocks per fiber clipped to the window; return union intervals."""
    per_fiber: Dict[str, List[Tuple[int, int]]] = {}
    for r in records:
        if r.chrom != chrom:
            continue
        for bs, be in r.blocks:
            s = max(bs, win_start)
            e = min(be, win_end)
            if s < e:
                per_fiber.setdefault(r.name, []).append((s, e))
    merged: List[Tuple[int, int]] = []
    for spans in per_fiber.values():
        spans.sort()
        cs, ce = spans[0]
        for s, e in spans[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                merged.append((cs, ce))
                cs, ce = s, e
        merged.append((cs, ce))
    return merged


def block_signal(
    records: Sequence[BedRecord],
    chrom: str,
    win_start: int,
    win_end: int,
    uniq_fibers: bool = True,
) -> np.ndarray:
    """Per-base element signal inside [win_start, win_end).

    With ``uniq_fibers=True`` each fiber contributes at most 1 per base
    (fiber coverage); with False every block adds one (mark density).
    """
    length = win_end - win_start
    diff = np.zeros(length + 1, dtype=np.int64)
    if uniq_fibers:
        spans = _merged_intervals(records, chrom, win_start, win_end)
    else:
        spans = []
        for r in records:
            if r.chrom != chrom:
                continue
            for bs, be in r.blocks:
                s = max(bs, win_start)
                e = min(be, win_end)
                if s < e:
                    spans.append((s, e))
    for s, e in spans:
        diff[s - win_start] += 1
        diff[e - win_start] -= 1
    return np.cumsum(diff[:length])


def coverage_signal(
    records_a: Sequence[BedRecord],
    records_b: Sequence[BedRecord],
    chrom: str,
    win_start: int,
    win_end: int,
) -> np.ndarray:
    """Per-base union fiber coverage from two tiling element sets (e.g. nuc+msp)."""
    return block_signal(list(records_a) + list(records_b), chrom, win_start, win_end)


def bh_fdr(pvals: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg FDR correction; returns q-values in input order."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank_from_end, i in enumerate(reversed(order)):
        rank = n - rank_from_end
        val = min(prev, pvals[i] * n / rank)
        q[i] = val
        prev = val
    return q


def atomic_write(path: str, write_fn) -> None:
    """Write via *write_fn(file_handle)* to a temp file, then os.replace."""
    tmp = f"{path}.tmp"
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(tmp, "w", encoding="utf-8") as fh:
        write_fn(fh)
    os.replace(tmp, path)


def write_tsv(path: str, header: Sequence[str], rows: Sequence[Sequence]) -> None:
    """Atomically write a TSV table."""

    def _w(fh):
        fh.write("\t".join(str(x) for x in header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    atomic_write(path, _w)


def safe_div(a: float, b: float) -> float:
    """Division that returns 0.0 when the denominator is zero."""
    return a / b if b else 0.0

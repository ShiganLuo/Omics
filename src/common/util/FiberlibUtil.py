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


def iter_bed_chrom(path: str) -> Iterator[Tuple[str, List[BedRecord]]]:
    """Stream a BED file grouped per chromosome, in file order.

    Each yielded chunk holds only the records of ONE chromosome, so peak
    memory is bounded by the largest chromosome instead of the whole file
    (the whole-file `read_bed` list OOMs on multi-GB Fiber-seq BEDs).
    Requires chromosome-grouped input (sort -k1,1); a chromosome that
    reappears later yields a second chunk for it.
    """
    chunk_chrom: Optional[str] = None
    chunk: List[BedRecord] = []
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
                continue
            chrom = f[0]
            if chrom != chunk_chrom:
                if chunk:
                    assert chunk_chrom is not None
                    yield chunk_chrom, chunk
                chunk_chrom = chrom
                chunk = []
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
            chunk.append(BedRecord(chrom, start, end, name, score, strand, blocks))
    if chunk:
        assert chunk_chrom is not None
        yield chunk_chrom, chunk


def stream_bed_records(path: str) -> Iterator[BedRecord]:
    """Stream a BED file record by record (O(1) memory, file order).

    Same parsing rules as `read_bed`; callers that only need records near
    a few anchors should retain selectively instead of materializing lists.
    """
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
                continue
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
            yield BedRecord(chrom, start, end, name, score, strand, blocks)


class ChromSweep:
    """Position-sweep join between big BEDs and a small set of anchors.

    Chromosome processing order follows the files themselves (all paths
    must be chromosome-grouped in the same order); within a chromosome,
    anchors are evaluated in start order. Only records overlapping the
    current anchor window are retained, so memory stays proportional to
    the anchor neighborhood, not the file. Anchors on chromosomes absent
    from the files yield empty per-path lists (same as zero coverage).
    """

    def __init__(self, paths: Sequence[str],
                 anchors: Sequence[Tuple[str, int, int]]) -> None:
        self._paths = list(paths)
        self._streams = [stream_bed_records(p) for p in self._paths]
        self._pending: List[Optional[BedRecord]] = [
            next(s, None) for s in self._streams
        ]
        self._active: List[List[BedRecord]] = [[] for _ in self._paths]
        self._last: List[Optional[str]] = [None] * len(self._paths)
        self._left: List[set] = [set() for _ in self._paths]
        self._anchors = list(anchors)

    def _advance_path(self, i: int, chrom: str) -> None:
        """Position stream *i* on chromosome *chrom* (or exhaust it)."""
        p = self._pending[i]
        prev = self._last[i]
        exited: set = set()
        while p is not None and p.chrom != chrom:
            x = p.chrom
            if x != prev:
                if prev is not None:
                    exited.add(prev)
                prev = x
            p = next(self._streams[i], None)
            self._pending[i] = p
        if p is not None:
            if prev is not None and prev != chrom:
                exited.add(prev)
                if chrom in exited or chrom in self._left[i]:
                    raise ValueError(
                        f"{self._paths[i]}: chromosome {chrom} re-appears after "
                        f"later chromosomes; BED not chromosome-grouped "
                        f"(run `sort -k1,1` on the input)")
            self._last[i] = chrom
        else:
            if prev is not None:
                exited.add(prev)
                if chrom in exited:
                    raise ValueError(
                        f"{self._paths[i]}: records for chromosome {chrom} were "
                        f"already passed; BED chromosome order diverges between "
                        f"files (run `sort -k1,1` on every input)")
            self._last[i] = None
        self._left[i] |= exited

    def _exit_block(self, i: int, chrom: str) -> None:
        """Move stream *i* past the remainder of its current *chrom* block."""
        p = self._pending[i]
        while p is not None and p.chrom == chrom:
            p = next(self._streams[i], None)
            self._pending[i] = p
        self._left[i].add(chrom)
        if p is not None:
            nxt = p.chrom
            if nxt in self._left[i]:
                raise ValueError(
                    f"{self._paths[i]}: chromosome {nxt} re-appears after "
                    f"later chromosomes; BED not chromosome-grouped "
                    f"(run `sort -k1,1` on the input)")
            self._last[i] = nxt
        else:
            self._last[i] = None

    def _anchors_by_chrom(self) -> Dict[str, List[Tuple[int, int, int]]]:
        buckets: Dict[str, List[Tuple[int, int, int]]] = {}
        for k, (c, s, e) in enumerate(self._anchors):
            buckets.setdefault(c, []).append((s, e, k))
        for c in buckets:
            buckets[c].sort()
        return buckets

    def advance(self) -> Iterator[Tuple[int, Tuple[str, int, int], List[List[BedRecord]]]]:
        """Yield (original_anchor_index, anchor, per-path overlapping records)."""
        buckets = self._anchors_by_chrom()
        seen: set = set()
        n = len(self._paths)
        while True:
            head = self._pending[0]
            if head is None:
                break
            chrom = head.chrom
            for i in range(n):
                self._advance_path(i, chrom)
                self._active[i] = []
            for (s, e, k) in buckets.get(chrom, []):
                seen.add(k)
                per_path: List[List[BedRecord]] = []
                for i in range(n):
                    p = self._pending[i]
                    while p is not None and p.chrom == chrom and p.start < e:
                        if self._active[i] and p.start < self._active[i][-1].start:
                            raise ValueError(
                                f"{self._paths[i]}: BED not sorted by start within "
                                f"chromosome {chrom} at position {p.start}; "
                                f"run `sort -k1,1 -k2,2n` on the input")
                        self._active[i].append(p)
                        p = next(self._streams[i], None)
                        self._pending[i] = p
                    self._active[i] = [r for r in self._active[i] if r.end > s]
                    per_path.append(
                        [r for r in self._active[i] if r.start < e and r.end > s])
                yield k, (chrom, s, e), per_path
            self._exit_block(0, chrom)
        for k, (c, s, e) in enumerate(self._anchors):
            if k not in seen:
                yield k, (c, s, e), [[] for _ in range(n)]


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
    Accumulates the difference array directly from records -- no
    intermediate span list, so memory is O(window) regardless of input.
    """
    length = win_end - win_start
    if length <= 0:
        return np.zeros(0, dtype=np.int64)
    diff = np.zeros(length + 1, dtype=np.int64)
    if uniq_fibers:
        # Union per fiber NAME across ALL passed records (e.g. the msp+nuc
        # concatenation of track: one fiber's MSP and nucleosome blocks must
        # merge into a single coverage contribution). Merging happens per
        # fiber and emits straight into the difference array -- the only
        # state is the window-local span list of the current fiber.
        per_fiber: Dict[str, List[Tuple[int, int]]] = {}
        for r in records:
            if r.chrom != chrom:
                continue
            for bs, be in r.blocks:
                s = bs if bs > win_start else win_start
                e = be if be < win_end else win_end
                if s < e:
                    spans = per_fiber.get(r.name)
                    if spans is None:
                        per_fiber[r.name] = [(s, e)]
                    else:
                        spans.append((s, e))
        for spans in per_fiber.values():
            if len(spans) == 1:
                s, e = spans[0]
                diff[s - win_start] += 1
                diff[e - win_start] -= 1
                continue
            spans.sort()
            cs, ce = spans[0]
            for s, e in spans[1:]:
                if s <= ce:
                    if e > ce:
                        ce = e
                else:
                    diff[cs - win_start] += 1
                    diff[ce - win_start] -= 1
                    cs, ce = s, e
            diff[cs - win_start] += 1
            diff[ce - win_start] -= 1
    else:
        for r in records:
            if r.chrom != chrom:
                continue
            for bs, be in r.blocks:
                s = bs if bs > win_start else win_start
                e = be if be < win_end else win_end
                if s < e:
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

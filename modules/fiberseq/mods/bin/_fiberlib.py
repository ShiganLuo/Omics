"""Module-private helpers for the fiberseq mods bin scripts.

Verbatim copies from src/common/util/ (modules.md 脚本撰写规范 #0:
module bin scripts are self-contained — no imports from outside the
module directory, helpers are copied in with provenance). Contains only
the pieces the mods scripts use: BED/BED12 parsing, the ChromSweep
position-sweep join, and the shared logger factory.

- setup_logger            <- src/common/util/LogUtil.py
- BedRecord, open_maybe_gzip, read_bed, stream_bed_records, ChromSweep
                          <- src/common/util/FiberlibUtil.py
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Sequence, Tuple, TextIO
import logging


# ---------------------------------------------------------------------------
# Copied verbatim from src/common/util/LogUtil.py.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Copied verbatim from src/common/util/FiberlibUtil.py.
# ---------------------------------------------------------------------------
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

#!/usr/bin/env python3
"""Shared multiprocess prefetching for the Fiber-seq profiling scripts.

Parsing the multi-GB element BEDs (gzip) is single-thread bound. This
helper runs ONE background parser process per BED that walks the file once
and ships per-chromosome data to the parent as compact numpy arrays
(fiber id + start + end per block). Fiber identity across FILES is
preserved by an IdMapper in the parent: the parser sends the names of
newly-seen fibers once per chromosome, and the parent maps every file's
parser-local ids into one shared namespace (a read seen in msp, nuc and
m6a collapses to one id), which the union-coverage / classification code
requires. The parent only does array math while the next chromosome is
already being parsed -- no giant BedRecord pickles.
"""
from __future__ import annotations

import multiprocessing as mp
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np

from src.common.util.FiberlibUtil import stream_bed_records  # noqa: E402

_PARSER_ERROR = "__parser_error__"


def _parser_worker(path: str, conn) -> None:
    """Parse BED into per-chrom compact block arrays and ship them over a pipe."""
    try:
        names: Dict[str, int] = {}
        cur: Optional[str] = None
        fids: List[int] = []
        starts: List[int] = []
        ends: List[int] = []
        new_pairs: List[Tuple[int, str]] = []

        def flush():
            if cur is None:
                return
            conn.send((
                cur,
                np.array(fids, dtype=np.int64),
                np.array(starts, dtype=np.int64),
                np.array(ends, dtype=np.int64),
                new_pairs,
            ))

        for r in stream_bed_records(path):
            if r.chrom != cur:
                flush()
                cur = r.chrom
                fids, starts, ends, new_pairs = [], [], [], []
            fid = names.get(r.name, -1)
            if fid < 0:
                fid = len(names)
                names[r.name] = fid
                new_pairs.append((fid, r.name))
            for bs, be in r.blocks:
                if bs < be:
                    fids.append(fid)
                    starts.append(bs)
                    ends.append(be)
        flush()
        conn.send(None)  # EOF
    except Exception as e:  # propagate to the parent
        try:
            conn.send((_PARSER_ERROR, e))
        except Exception:
            pass
    finally:
        conn.close()


class IdMapper:
    """Maps per-file parser fiber ids into one shared unified namespace.

    One IdMapper instance may serve multiple parser streams (e.g. the msp
    and the nuc BED): each stream keeps its own dense parser-id -> unified
    id lookup (``extend(stream, pairs)`` / ``arr(stream)``), while identical
    fiber names collapse to the same unified id across streams.
    """

    def __init__(self) -> None:
        self._names: Dict[str, int] = {}
        self._next = 0
        self._maps: Dict[str, List[int]] = {}

    def extend(self, stream: str, new_pairs: List[Tuple[int, str]]) -> None:
        m = self._maps.setdefault(stream, [])
        for pfid, name in new_pairs:
            if pfid != len(m):
                raise ValueError(
                    f"non-contiguous parser fiber id {pfid} for stream "
                    f"'{stream}' (expected {len(m)}); stream corrupted")
            u = self._names.get(name, -1)
            if u < 0:
                u = self._next
                self._next += 1
                self._names[name] = u
            m.append(u)

    def arr(self, stream: str) -> np.ndarray:
        return np.asarray(self._maps[stream], dtype=np.int64)


class ChromBlockStream:
    """Iterator of (chrom, fid, start, end, new_pairs) block arrays for a BED.

    Yields chromosomes in file order (chromosome-grouped input), each once.
    ``new_pairs`` lists (parser_fid, name) of fibers first seen in that
    chromosome; feed them to the shared IdMapper before using ``fid``.
    """

    def __init__(self, path: str):
        ctx = mp.get_context("fork")
        parent_conn, child_conn = ctx.Pipe(duplex=False)
        self._proc = ctx.Process(target=_parser_worker, args=(path, child_conn),
                                 daemon=True)
        self._proc.start()
        child_conn.close()
        self._conn = parent_conn
        self._path = path
        self._done = False

    def __iter__(self) -> Iterator[Tuple[str, np.ndarray, np.ndarray, np.ndarray, list]]:
        while True:
            item = self._conn.recv()
            if item is None:
                self._done = True
                self._proc.join()
                break
            if item[0] == _PARSER_ERROR:
                raise RuntimeError(f"{self._path}: parser failed: {item[1]}")
            yield item

    def close(self) -> None:
        if not self._done and self._proc.is_alive():
            self._proc.terminate()
        self._conn.close()

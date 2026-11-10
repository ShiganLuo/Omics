#!/usr/bin/env python3
"""Shared utilities for eccDNA repeat-content analyses.

Provides three data access layers used by the sibling scripts
(`breakpoint_repeat.py`, `ecdna_repeat_enrichment.py`,
`junction_microhomology.py`):

- :class:`RmskIndex` — UCSC rmsk GTF (GENCODE TE track) parsed into per-chrom
  element lists with merged interval views for fast overlap-bp queries.
- :class:`FastaReader` — .fai-indexed random access to the reference genome
  (stdlib only; no pysam dependency).
- circle loaders — decoil `reconstruct.ecDNA.filtered.bed` fragments and
  `reconstruct.links.ecDNA.filtered.txt` junctions.

Notes
-----
Node ids in the links file are 1-based data-row indices of the matching
`reconstruct.ecDNA.filtered.bed`.  A link ``n1.side1 <-> n2.side0`` joins the
right end of fragment ``n1`` to the left end of fragment ``n2`` (side 1 =
end-coordinate side, side 0 = start-coordinate side).  All current fragments
are ``+`` strand; ``-`` strand is handled defensively via molecule
orientation (3'-terminus of the upstream fragment meets the 5'-terminus of
the downstream fragment).
"""
from __future__ import annotations

import argparse
import gzip
import logging
import os
import pickle
import re
from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import dataclass, field

logger = logging.getLogger("ecdna_repeat_utils")

# key for merged-interval queries: 'any', ('fam', idx) or ('cls', idx)
RmskKey = str | tuple[str, int]

_ATTR_RE = {
    "family_id": re.compile(r'family_id "([^"]*)"'),
    "class_id": re.compile(r'class_id "([^"]*)"'),
}


# --------------------------------------------------------------------------- #
# RepeatMasker annotation
# --------------------------------------------------------------------------- #
class RmskIndex:
    """RepeatMasker element index over a UCSC rmsk GTF.

    Parameters
    ----------
    gtf_path : str
        Plain or gzipped GTF with ``family_id`` and ``class_id`` attributes
        (e.g. GRCm39_GENCODE_rmsk_TE.gtf[.gz]).
    cache_dir : str, optional
        Directory for a pickle cache of the parsed index (keyed by source
        path + size + mtime).  ``None`` disables caching.

    Attributes
    ----------
    families : list of str
        Family names indexed by ``fam_idx``.
    classes : list of str
        Class names indexed by ``cls_idx``.
    family_class : list of int
        ``cls_idx`` of each family (``family_class[fam_idx]``).
    elements : dict of str to list of tuple
        ``chrom -> [(start, end, fam_idx, cls_idx), ...]`` sorted by start
        (0-based half-open).
    """

    def __init__(self, gtf_path: str, cache_dir: str | None = None):
        self.gtf_path = gtf_path
        cached = self._load_cache(cache_dir) if cache_dir else None
        if cached is not None:
            (self.elements, self.families,
             self.classes, self.family_class) = cached
        else:
            (self.elements, self.families,
             self.classes, self.family_class) = self._parse_gtf()
            if cache_dir:
                self._save_cache(cache_dir)
        # merged-interval caches: key -> chrom -> (intervals, starts, ends)
        self._merged: dict[RmskKey, dict[str, tuple]] = {}
        self._built: set[str] = set()
        self._elem_starts: dict[str, list[int]] = {}
        self.max_elem_len = max(
            (e - s for els in self.elements.values() for s, e, _, _ in els),
            default=0,
        )

    # -- parsing ----------------------------------------------------------- #
    def _parse_gtf(self) -> tuple[dict[str, list[tuple[int, int, int, int]]],
                                  list[str], list[str], list[int]]:
        elements: dict[str, list[tuple[int, int, int, int]]] = {}
        family_names: list[str] = []
        class_names: list[str] = []
        family_class: list[int] = []
        fam_idx: dict[str, int] = {}
        cls_idx: dict[str, int] = {}
        opener = gzip.open if self.gtf_path.endswith(".gz") else open
        n = 0
        with opener(self.gtf_path, "rt") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                attrs = f[8]
                m = _ATTR_RE["family_id"].search(attrs)
                fam = m.group(1) if m else "Unknown"
                m = _ATTR_RE["class_id"].search(attrs)
                cls = m.group(1) if m else "Unknown"
                if cls not in cls_idx:
                    cls_idx[cls] = len(class_names)
                    class_names.append(cls)
                if fam not in fam_idx:
                    fam_idx[fam] = len(family_names)
                    family_names.append(fam)
                    family_class.append(cls_idx[cls])
                elements.setdefault(f[0], []).append(
                    (int(f[3]) - 1, int(f[4]), fam_idx[fam], cls_idx[cls]))
                n += 1
        for els in elements.values():
            els.sort()
        logger.info(f"Loaded {n} rmsk elements "
                    f"({len(family_names)} families, {len(class_names)} classes) "
                    f"from {self.gtf_path}")
        return elements, family_names, class_names, family_class

    # -- pickle cache ------------------------------------------------------ #
    def _cache_path(self, cache_dir: str) -> str:
        st = os.stat(self.gtf_path)
        key = (f"v2.{os.path.basename(self.gtf_path)}."
               f"{st.st_size}.{st.st_mtime_ns}")
        return os.path.join(cache_dir, f"rmsk_index.{key}.pkl")

    def _load_cache(self, cache_dir: str):
        path = self._cache_path(cache_dir)
        if not os.path.exists(path):
            return None
        try:
            with open(path, "rb") as fh:
                data = pickle.load(fh)
            logger.info(f"Loaded rmsk index cache {path}")
            return data
        except Exception as exc:  # noqa: BLE001 — cache miss is not fatal
            logger.warning(f"Cache {path} unreadable ({exc}); re-parsing GTF")
            return None

    def _save_cache(self, cache_dir: str) -> None:
        os.makedirs(cache_dir, exist_ok=True)
        path = self._cache_path(cache_dir)
        with open(path, "wb") as fh:
            pickle.dump((self.elements, self.families,
                         self.classes, self.family_class), fh, protocol=4)
        logger.info(f"Saved rmsk index cache {path}")

    # -- merged interval views --------------------------------------------- #
    @staticmethod
    def _merge(intervals: list[tuple[int, int]]) -> tuple[
            list[tuple[int, int]], list[int], list[int]]:
        out: list[tuple[int, int]] = []
        for s, e in sorted(intervals):
            if out and s <= out[-1][1]:
                if e > out[-1][1]:
                    out[-1] = (out[-1][0], e)
            else:
                out.append((s, e))
        starts = [s for s, _ in out]
        ends = [e for _, e in out]
        return out, starts, ends

    def _view(self, chrom: str, key: RmskKey) -> tuple[
            list[tuple[int, int]], list[int], list[int]]:
        """Return ``(intervals, starts, ends)`` merged for ``key`` on chrom.

        All keys of a chromosome are built in one pass the first time the
        chromosome is queried (``_build_chrom``); missing keys are empty.
        """
        if chrom not in self._built:
            self._build_chrom(chrom)
        empty: tuple[list[tuple[int, int]], list[int], list[int]] = (
            [], [], [])
        return self._merged.get(key, {}).get(chrom, empty)

    def _build_chrom(self, chrom: str) -> None:
        """Merge every key of ``chrom`` in a single scan over its elements."""
        buckets: dict[RmskKey, list[tuple[int, int]]] = defaultdict(list)
        for s, e, fi, ci in self.elements.get(chrom, []):
            buckets["any"].append((s, e))
            buckets[("fam", fi)].append((s, e))
            buckets[("cls", ci)].append((s, e))
        for key, iv in buckets.items():
            self._merged.setdefault(key, {})[chrom] = self._merge(iv)
        self._built.add(chrom)

    def merged(self, chrom: str, key: RmskKey = "any") -> list[tuple[int, int]]:
        """Merged disjoint intervals on ``chrom`` for key ``'any'``,
        ``('fam', idx)`` or ``('cls', idx)``."""
        return self._view(chrom, key)[0]

    # -- queries ------------------------------------------------------------ #
    def overlap_bp(self, chrom: str, start: int, end: int,
                   key: RmskKey = "any") -> int:
        """Bases of ``[start, end)`` covered by merged intervals of ``key``."""
        _, starts, ends = self._view(chrom, key)
        i = bisect_right(ends, start)
        total = 0
        while i < len(starts) and starts[i] < end:
            total += min(end, ends[i]) - max(start, starts[i])
            i += 1
        return total

    def contains(self, chrom: str, pos: int, key: RmskKey = "any") -> bool:
        """True when point ``pos`` lies inside a merged interval of ``key``."""
        _, starts, ends = self._view(chrom, key)
        i = bisect_left(ends, pos + 1)
        return i < len(starts) and starts[i] <= pos

    def nearest_distance(self, chrom: str, pos: int,
                         key: RmskKey = "any") -> int | None:
        """Distance from ``pos`` to the nearest interval of ``key`` (0 when
        inside); ``None`` when the chromosome has no such interval."""
        _, starts, ends = self._view(chrom, key)
        if not starts:
            return None
        i = bisect_left(ends, pos + 1)
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(starts):
                s, e = starts[j], ends[j]
                d = 0 if s <= pos < e else min(abs(pos - s), abs(pos - (e - 1)))
                best = d if best is None else min(best, d)
        return best

    def elements_overlapping(self, chrom: str, start: int, end: int,
                             ) -> list[tuple[int, int, int, int]]:
        """Raw elements overlapping ``[start, end)`` (may nest/overlap)."""
        els = self.elements.get(chrom, [])
        if chrom not in self._elem_starts:
            self._elem_starts[chrom] = [s for s, _, _, _ in els]
        starts = self._elem_starts[chrom]
        i = bisect_right(starts, end - 1)
        out = []
        j = i - 1
        while j >= 0:
            s, e, fi, ci = els[j]
            if e > start and s < end:
                out.append(els[j])
            if s + self.max_elem_len < start:
                break
            j -= 1
        return out

    def family_bp(self, key: RmskKey) -> int:
        """Total genome bases covered by merged intervals of ``key``."""
        total = 0
        for chrom, els in self.elements.items():
            if key == "any" or any(
                    (("fam", fi) == key or ("cls", ci) == key)
                    for _, _, fi, ci in els):
                total += sum(e - s for s, e in self.merged(chrom, key))
        return total


# --------------------------------------------------------------------------- #
# Reference genome
# --------------------------------------------------------------------------- #
class FastaReader:
    """Minimal .fai-indexed FASTA random access (stdlib only)."""

    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self._fh = open(fasta_path, "rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        fai = fasta_path + ".fai"
        with open(fai) as fh:
            for line in fh:
                f = line.split("\t")
                self.index[f[0]] = (int(f[1]), int(f[2]),
                                    int(f[3]), int(f[4]))
        logger.info(f"FASTA index: {len(self.index)} sequences from {fai}")

    def chrom_sizes(self) -> dict[str, int]:
        return {c: v[0] for c, v in self.index.items()}

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Uppercase sequence for ``[start, end)`` (0-based half-open)."""
        if chrom not in self.index:
            raise KeyError(f"{chrom} not in {self.fasta_path}.fai")
        length, offset, line_bases, line_width = self.index[chrom]
        start = max(0, start)
        end = min(end, length)
        if end <= start:
            return ""
        byte_start = offset + (start // line_bases) * line_width + start % line_bases
        last = end - 1
        byte_end = offset + (last // line_bases) * line_width + last % line_bases
        self._fh.seek(byte_start)
        raw = self._fh.read(byte_end - byte_start + 1).decode("ascii")
        return raw.replace("\n", "").replace("\r", "").upper()

    def close(self) -> None:
        self._fh.close()


_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    """Reverse complement of an uppercase nucleotide string."""
    return seq.translate(_COMP)[::-1]


# --------------------------------------------------------------------------- #
# decoil circle structures
# --------------------------------------------------------------------------- #
@dataclass
class Fragment:
    """One merged fragment row of ``reconstruct.ecDNA.filtered.bed``."""
    row: int          # 1-based data-row index (= links node id)
    chrom: str
    start: int        # 0-based half-open
    end: int
    circ_id: str
    fragment_id: str
    strand: str
    sample: str

    @property
    def oriented_start_coord(self) -> int:
        """Coordinate of the fragment's 5' terminus in molecule orientation."""
        return self.start if self.strand == "+" else self.end

    @property
    def oriented_end_coord(self) -> int:
        """Coordinate of the fragment's 3' terminus in molecule orientation."""
        return self.end if self.strand == "+" else self.start


@dataclass
class Junction:
    """One row of ``reconstruct.links.ecDNA.filtered.txt``."""
    sample: str
    circ_id: str
    n1: int
    n2: int
    side1: int
    side2: int
    frag1: Fragment = field(repr=False)
    frag2: Fragment = field(repr=False)


def sample_path_spec(spec: str) -> tuple[str, str]:
    """argparse type for ``SAMPLE=PATH`` CLI values.

    Parameters
    ----------
    spec : str
        ``SAMPLE=PATH`` where PATH is an existing file (any location).

    Returns
    -------
    (sample, path) : tuple of str

    Raises
    ------
    argparse.ArgumentTypeError
        If the spec is not ``SAMPLE=PATH`` or PATH does not exist.
    """
    if "=" not in spec:
        raise argparse.ArgumentTypeError(
            f"expected SAMPLE=PATH, got {spec!r}")
    sample, path = spec.split("=", 1)
    if not sample or not path:
        raise argparse.ArgumentTypeError(
            f"expected SAMPLE=PATH, got {spec!r}")
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"input file not found: {path}")
    return sample, path


def sample_map(pairs: list[tuple[str, str]]) -> dict[str, str]:
    """Build ``{sample: path}`` from parsed specs, rejecting duplicates.

    Parameters
    ----------
    pairs : list of (sample, path)

    Returns
    -------
    dict of str to str

    Raises
    ------
    ValueError
        If a sample name appears more than once.
    """
    out: dict[str, str] = {}
    for sample, path in pairs:
        if sample in out:
            raise ValueError(f"duplicate sample {sample!r}: "
                             f"{out[sample]} and {path}")
        out[sample] = path
    return out


def load_fragments(beds: dict[str, str]) -> dict[str, list[Fragment]]:
    """Load per-sample fragment rows of ``reconstruct.ecDNA.filtered.bed``.

    Parameters
    ----------
    beds : dict of str to str
        ``{sample: bed_path}`` — paths may be anywhere on disk.

    Returns
    -------
    fragments : dict of str to list of Fragment
        Per sample, rows in file order (``Fragment.row`` = links node id).
    """
    fragments: dict[str, list[Fragment]] = {}
    for sample, bed in beds.items():
        rows: list[Fragment] = []
        with open(bed) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 6:
                    continue
                rows.append(Fragment(row=len(rows) + 1, chrom=f[0],
                                     start=int(f[1]), end=int(f[2]),
                                     circ_id=f[3], fragment_id=f[4],
                                     strand=f[5], sample=sample))
        fragments[sample] = rows
        logger.info(f"{sample}: {len(rows)} fragments from {bed}")
    return fragments


def load_junctions(links: dict[str, str],
                   fragments: dict[str, list[Fragment]]) -> list[Junction]:
    """Load ``reconstruct.links.ecDNA.filtered.txt`` rows resolved to fragments.

    Parameters
    ----------
    links : dict of str to str
        ``{sample: links_path}``; every sample must also be a key of
        ``fragments`` (node ids index that sample's fragment rows).
    fragments : dict of str to list of Fragment
        Output of :func:`load_fragments`.

    Returns
    -------
    junctions : list of Junction
    """
    missing = set(links) - set(fragments)
    if missing:
        raise KeyError(f"links given without fragments for: {sorted(missing)}")
    junctions: list[Junction] = []
    for sample, path in links.items():
        rows = fragments[sample]
        n_junc = 0
        with open(path) as fh:
            fh.readline()  # header
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                n1, n2, s1, s2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
                if not (1 <= n1 <= len(rows) and 1 <= n2 <= len(rows)):
                    logger.warning(f"{sample}: link node out of range, "
                                   f"skipped ({line.strip()})")
                    continue
                junctions.append(Junction(sample=sample, circ_id=f[4],
                                          n1=n1, n2=n2, side1=s1, side2=s2,
                                          frag1=rows[n1 - 1],
                                          frag2=rows[n2 - 1]))
                n_junc += 1
        logger.info(f"{sample}: {n_junc} junctions from {path}")
    return junctions


# --------------------------------------------------------------------------- #
# statistics helpers
# --------------------------------------------------------------------------- #
def empirical_p(obs: float, draws: list[float], greater: bool = True) -> float:
    """One-sided empirical p-value ``(1 + #{draws as/extreme}) / (1 + n)``."""
    if not draws:
        return float("nan")
    n_ext = sum(1 for d in draws if (d >= obs if greater else d <= obs))
    return (1 + n_ext) / (1 + len(draws))


def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR adjusted p-values (input order preserved)."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    prev = 1.0
    for rank_from_end, i in enumerate(reversed(order)):
        rank = n - rank_from_end
        val = min(prev, pvals[i] * n / rank)
        adj[i] = val
        prev = val
    return adj

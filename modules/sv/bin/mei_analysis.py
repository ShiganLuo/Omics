#!/usr/bin/env python
"""MEI (mobile element insertion) analysis for SURVIVOR-merged SV VCFs.

Consumes the merged VCF produced by the Omics sv module (SURVIVOR merge of
pbsv calls). Insertion ALT alleles carry the inserted sequence (pbsv AAL
kept by SURVIVOR), which is the primary evidence for TE classification.

Pipeline (per INS event):
  1. Extract inserted sequence from ALT (fallback: none -> dropped with count)
  2. Classify TE family by RepeatMasker (subprocess, -no_is, species lib);
     fall back to k-mer matching against an rmsk-derived TE consensus FASTA
     when RepeatMasker is unavailable
  3. MEI fingerprints: poly(A) tail in inserted seq, TSD (target site
     duplication) by comparing REF flank in sample RAL vs ALT junctions,
     5'-truncation proxy (coverage of the TE consensus by the hit)
  4. Context annotation: gene overlap (GTF), insertion into pre-existing TE
     regions (rmsk TE bed) -> exon / intron / intergenic / TE-into-TE
  5. Full-length L1 detection: >=6000 bp + L1 family + full-length-ish hit

Stages (the entry runs analysis + figures; --skip-plot to skip figures):
  analysis (this file) -> 6-figure suite (utils/mei_plot.py).
TE insertion dating (utils/te_age_analysis.py) is a SEPARATE analysis,
run standalone with --mei-dir afterwards.

Outputs (all under --outdir):
  mei_events.tsv        per-event table (all INS, classified or not)
  mei_summary.json      counts per (group, family), full-length L1 list
  mei_family_by_group.tsv   family x SUPP_VEC matrix
  mei_context.tsv       context (exon/intron/intergenic/TE) x group x family
  mei_full_length_L1.tsv    candidate active L1 insertions
  mei_log.txt           run log
"""

import argparse
import gzip
import hashlib
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from typing import (TYPE_CHECKING, Any, Callable, Dict, List, Optional, Tuple,
                    Union)

import pandas as pd

from utils.common import setup_logger

if TYPE_CHECKING:
    import pysam  # pyright: ignore[reportMissingImports] — run env has pysam

# record dict produced by load_merged_vcf and enriched downstream
Rec = Dict[str, Any]
# TE bucket spec entry (see the bucket schema below)
Bucket = Dict[str, Any]

logger = setup_logger("mei_analysis")


# ---------------------------------------------------------------- VCF parsing
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_alt_insertion(ref: Optional[str], alt: Optional[str]) -> Optional[str]:
    """Return the pure inserted sequence of a simple insertion ALT.

    Handles symbolic ``<INS>`` alleles, plain-sequence ALTs and
    left-aligned padded ALTs sharing a prefix (and, rarely, a suffix)
    with REF — the common pbsv/SURVIVOR layout ``REF=A,
    ALT=A<inserted>``.

    Parameters
    ----------
    ref : str or None
        REF allele from the VCF record.
    alt : str or None
        ALT allele from the VCF record.

    Returns
    -------
    str or None
        The inserted sequence with shared REF prefix/suffix trimmed, or
        None for symbolic alleles (``<INS>``, ``.``), ALT == REF, or
        when nothing remains after trimming.

    Examples
    --------
    >>> parse_alt_insertion("A", "AACGT")
    'CGT'
    """
    if alt is None or alt.startswith("<") or alt == ".":
        return None
    ref = ref or ""  # None REF: trim against an empty anchor
    if alt == ref:
        return None
    # shared-prefix trim (REF=A ALT=AACGT... -> ACGT...)
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    ins = alt[i:]
    # also trim a shared suffix (rare, keeps middle)
    if len(ins) > 0 and len(ref) - i > 0:
        j = 0
        while (j < len(ins) and j < len(ref) - i
               and ins[len(ins) - 1 - j] == ref[len(ref) - 1 - j]):
            j += 1
        if j:
            ins = ins[: len(ins) - j]
    return ins if len(ins) > 0 else None


def parse_info(info: str) -> Dict[str, Union[str, bool]]:
    """Parse a VCF INFO column into a key/value dict.

    Parameters
    ----------
    info : str
        The INFO column of a VCF record (``;``-separated ``key=value``
        or flag entries).

    Returns
    -------
    dict
        Mapping of INFO keys to their string values; flag entries map
        to True.
    """
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def load_merged_vcf(path: str) -> Tuple[List[Rec], List[str]]:
    """Load INS records with inserted sequence from a merged SV VCF.

    Only ``SVTYPE=INS`` records whose ALT carries an actual sequence are
    kept (symbolic/no-sequence ALTs are counted and dropped). The
    inserted sequence is the primary TE-classification evidence.

    Parameters
    ----------
    path : str
        SURVIVOR-merged SV VCF path (``.gz`` or plain text).

    Returns
    -------
    recs : list of dict
        One dict per INS record with keys ``chrom``, ``pos``, ``id``,
        ``svlen`` (INFO SVLEN, falling back to sequence length),
        ``ins_len``, ``seq`` (inserted sequence), ``supp_vec``,
        ``supp``, ``samples``.
    samples : list of str
        Sample column names from the ``#CHROM`` header line.
    """
    opener = gzip.open if path.endswith(".gz") else open
    recs = []
    dropped_symbolic = 0
    with opener(path, "rt") as fh:
        samples: List[str] = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = parse_info(f[7])
            if info.get("SVTYPE") != "INS":
                continue
            ins = parse_alt_insertion(f[3], f[4])
            if ins is None:
                dropped_symbolic += 1
                continue
            svlen = 0
            try:
                svlen = int(info.get("SVLEN", len(ins)))
            except ValueError:
                svlen = len(ins)
            dr = []
            for s in f[9:]:
                fmt = s.split(":")
                # DR is Number=2 "ref,variant support reads"; pbsv/SURVIVOR
                # writes a single combined value pre-merge -> take max field
                m = re.search(r"DR=([^;:\t]+)", s)  # not standard; fallback
                dr.append(fmt)
            recs.append({
                "chrom": f[0], "pos": int(f[1]), "id": f[2],
                "svlen": svlen, "ins_len": len(ins), "seq": ins,
                "supp_vec": info.get("SUPP_VEC", "."),
                "supp": info.get("SUPP", "."),
                "samples": samples,
            })
    logger.info(f"parsed INS records with sequence: {len(recs)} "
        f"(dropped symbolic/no-seq: {dropped_symbolic})")
    return recs, samples


# ------------------------------------------------------- RepeatMasker classify
def find_rm_home() -> Optional[str]:
    """RepeatMasker share dir (holds famdb.py + Libraries/famdb)."""
    rm = shutil.which("RepeatMasker")
    if not rm:
        return None
    home = os.path.normpath(
        os.path.join(os.path.dirname(rm), "..", "share", "RepeatMasker"))
    return home if os.path.isdir(home) else None


def build_famdb_dir(h5_paths: List[str], workdir: str) -> str:
    """Assemble a FamDB directory RepeatMasker can consume.

    A FamDB 'database' is a directory of partition h5 files plus a
    ``rmlib.config`` manifest, and it needs exactly one ROOT file (the
    one carrying Taxonomy/NamesCache — shipped with RepeatMasker as
    e.g. dfam39_full.0.h5). We combine RepeatMasker's own root/partition
    files with the user-provided partitions (symlinks) and rewrite the
    manifest. Returns the ``--libdir`` parent (its ``famdb/`` subdir is
    the database).
    """
    parent = os.path.join(workdir, "famdb_lib")
    famdb = os.path.join(parent, "famdb")
    os.makedirs(famdb, exist_ok=True)

    files: Dict[str, str] = {}  # basename -> source path
    home = find_rm_home()
    if home:
        bundled = os.path.join(home, "Libraries", "famdb")
        if os.path.isdir(bundled):
            for fn in sorted(os.listdir(bundled)):
                if fn.endswith(".h5"):
                    files[fn] = os.path.join(bundled, fn)
    for p in h5_paths:
        files[os.path.basename(p)] = p  # user partitions override same-name
    if not files:
        raise RuntimeError("no FamDB h5 files found (RepeatMasker bundled "
                           "root missing and no --dfam-h5 given)")

    for fn, src in files.items():
        dst = os.path.join(famdb, fn)
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(os.path.abspath(src), dst)

    entries = ",\n".join(
        f"    '{fn}' => {{\n      'size' => {os.path.getsize(os.path.join(famdb, fn))}\n    }}"
        for fn in sorted(files))
    listing = ", ".join(f"'{fn}'" for fn in sorted(files))
    with open(os.path.join(famdb, "rmlib.config"), "w") as fh:
        fh.write("{\n  'famdb_files' => {\n%s\n  },\n  'files' => [\n    %s\n  ]\n}\n"
                 % (entries, listing))
    logger.info(f"famdb dir ready: {famdb} ({len(files)} h5 file(s))")
    return parent


def export_species_lib(famdb_parent: str, species: str, out_fa: str) -> str:
    """Export the curated species consensus library via RepeatMasker's own
    famdb.py (stage 80 = all search stages, curated, ancestors+descendants).

    This is the official Dfam curation path: species-specific thresholds
    are baked into the ``[S:...]`` headers so RepeatMasker enforces them.
    """
    home = find_rm_home()
    famdb_py = os.path.join(home, "famdb.py") if home else None
    if not famdb_py or not os.path.exists(famdb_py):
        raise RuntimeError("RepeatMasker's famdb.py not found — cannot export "
                           "the species library from the Dfam h5")
    cmd = [sys.executable, famdb_py, "-i", os.path.join(famdb_parent, "famdb"),
           "families", species, "--curated", "--ancestors", "--descendants",
           "--include-class-in-name", "--stage", "80", "-f", "fasta_name"]
    with open(out_fa, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    n = sum(1 for line in open(out_fa) if line.startswith(">"))
    if r.returncode != 0 or n == 0:
        raise RuntimeError(
            f"famdb.py export failed ({r.returncode}, {n} seqs): "
            f"{r.stderr[:500]}")
    logger.info(f"exported {n} curated '{species}' consensuses -> {out_fa}")
    return out_fa


def run_repeatmasker(seqs: List[str], ids: List[str], workdir: str,
                     species_lib: Optional[str], threads: int = 4,
                     libdir: Optional[str] = None,
                     rm_species: Optional[str] = None,
                     config_stamp: Optional[Dict[str, Any]] = None,
                     force: bool = False,
                     ) -> Dict[str, Dict[str, Any]]:
    """Classify inserted sequences with RepeatMasker against a TE library.

    Writes ``ins.fa``, runs RepeatMasker on it and parses the ``.out``
    into per-record best hits. An existing ``ins.fa.out`` is reused
    (cache hit on reruns; the FASTA is only written when the ``.out``
    is absent). The ``.out`` layout is variable-width, so columns are
    anchored from the right edge (ID -> rep end/begin -> class/family ->
    repeat name -> strand); repeat names may contain ``/`` (e.g.
    ``ALR/Alpha``), which rules out left-anchored parsing.

    Parameters
    ----------
    seqs : list of str
        Inserted sequences, aligned with ``ids``.
    ids : list of str
        Record identifiers used as FASTA headers.
    workdir : str
        Working directory holding ``ins.fa`` / ``ins.fa.out``.
    species_lib : str
        Species TE consensus FASTA for RepeatMasker ``-lib``.
    threads : int, optional
        RepeatMasker ``-pa`` parallel threads (default 4).

    Returns
    -------
    dict
        ``{id: hit}`` with the longest-covered-span hit per record. Each
        hit has ``begin``, ``end``, ``repeat`` (RM repeat NAME),
        ``class`` (RM class/family, e.g. ``SINE/Alu``), ``family``
        (suffix after ``/`` or empty), ``sw``, ``strand``, ``n_hits``.
        Empty dict when RepeatMasker fails or produces no ``.out``.
    """
    fa = os.path.join(workdir, "ins.fa")
    out = os.path.join(workdir, "ins.fa.out")
    stamp_path = os.path.join(workdir, "rm_config.stamp.json")
    # cache key: classification config + the exact query set (content hash)
    ins_md5 = hashlib.md5(
        "\n".join(f"{i}\n{s}" for i, s in zip(ids, seqs)).encode()).hexdigest()
    stamp = json.dumps({**(config_stamp or {}), "ins_md5": ins_md5},
                       sort_keys=True)
    if force and os.path.exists(out):
        logger.info("--force-rerun: discarding cached ins.fa.out")
        os.remove(out)
    prev = None
    if os.path.exists(stamp_path):
        try:
            with open(stamp_path) as fh:
                prev = fh.read().strip()
        except OSError:
            prev = None
    if os.path.exists(out) and prev != stamp:
        logger.info("classification inputs changed since cached ins.fa.out "
                    "(library/mode/species/query set) — re-running RepeatMasker")
        os.remove(out)
    if not os.path.exists(out):
        with open(fa, "w") as fh:
            for i, (sid, s) in enumerate(zip(ids, seqs)):
                fh.write(f">{sid}\n")
                for j in range(0, len(s), 60):
                    fh.write(s[j : j + 60] + "\n")
        cmd = ["RepeatMasker", "-no_is", "-s", "-pa", str(threads)]
        if libdir and rm_species:
            cmd += ["-libdir", libdir, "-species", rm_species]
        else:
            cmd += ["-lib", species_lib]
        cmd += ["-dir", workdir, fa]
        logger.info(f"running RepeatMasker on {len(seqs)} sequences (this can take a while)")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            logger.info(f"RepeatMasker failed: {r.stderr[:500]}")
            return {}
        if not os.path.exists(out):
            logger.info("RepeatMasker .out missing")
            return {}
    else:
        logger.info(f"reusing cached {out}")
    with open(stamp_path, "w") as fh:
        fh.write(stamp)
    hits = defaultdict(list)
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split()
            if len(f) < 15 or not f[0].isdigit():
                continue
            # RepeatMasker .out layout is variable-width; use stable anchors:
            #   f[0]=SW score, f[4]=query id, f[5]/f[6]=query begin/end,
            #   then [strand], repeat name, class/family, rep-begin/end...
            # Parse from the RIGHT: last col = ID(*), before it rep end,
            # before rep begin, before class/family, before repeat name.
            try:
                sw = int(f[0])
                sid = f[4]
                begin, end = int(f[5]), int(f[6])
            except (ValueError, IndexError):
                continue
            # find the class/family token by scanning from the RIGHT: the last
            # three tokens are ID, and (in the 15-col layout) rep-end/rep-begin
            # sit left of the ID; class/family and repeat name are fixed
            # anchors relative to them. Repeat names may themselves contain
            # '/' (e.g. ALR/Alpha in the macaque lib), so anchoring on the
            # FIRST '/'-token from the left mis-parses satellite rows into
            # strand '+' / 'C' as the "repeat name".
            rep = cls = fam = ""
            strand = "+"
            # Right-anchored walk: skip the trailing '*' column if present,
            # land on the ID, then walk left across the whole repeat-position
            # block (numbers and (numbers), in either order), and the first
            # non-numeric token left of the block is the class/family; one
            # further left the repeat name (may itself contain '/', e.g.
            # ALR/Alpha), one further the strand.
            def _num(t):
                return t.isdigit()

            def _pnum(t):
                # (123) or (-1): the (left) column can be negative
                body = t.strip("()")
                return t.startswith("(") and body.lstrip("-").isdigit()

            k = len(f) - 1
            if f[k] == "*":
                k -= 1
            if k < 8 or not f[k].rstrip("*").isdigit():
                continue
            j = k - 1
            while j >= 8 and (_num(f[j]) or _pnum(f[j])):
                j -= 1
            if j < 8:
                continue
            cls = f[j]
            rep = f[j - 1].rstrip("?")
            strand = f[j - 2] if j - 2 >= 8 and f[j - 2] in ("+", "-", "C") else "+"
            fam = ""
            if "/" in cls:
                fam = cls.split("/", 1)[1]
            hits[sid].append({
                "begin": begin, "end": end, "repeat": rep,
                "class": cls, "family": fam, "sw": sw, "strand": strand,
            })
    # best hit per record: longest covered span
    best = {}
    for sid, hs in hits.items():
        hs.sort(key=lambda h: (h["end"] - h["begin"]), reverse=True)
        h = hs[0]
        h["n_hits"] = len(hs)
        best[sid] = h
    return best


# ---------------------------------------------------------------- fingerprints
def poly_a_tail(seq: str, min_len: int = 12, min_frac: float = 0.9) -> bool:
    """Detect an A-rich 3' tail (L1/Alu/SVA retrotransposition hallmark).

    Scans the last 60 bases for the longest uninterrupted A-run and the
    overall A fraction; a tail is called when either criterion passes,
    so interrupted tails still count.

    Parameters
    ----------
    seq : str
        Inserted sequence (3' end is the string end).
    min_len : int, optional
        Minimum uninterrupted A-run length to call a tail (default 12).
    min_frac : float, optional
        Minimum A fraction over the terminal window (default 0.9).

    Returns
    -------
    bool
        True when an A-rich tail is detected.
    """
    tail_win = seq[-60:] if len(seq) >= 60 else seq
    best_run = 0
    cur = 0
    for ch in tail_win:
        cur = cur + 1 if ch in "Aa" else 0
        best_run = max(best_run, cur)
    a_frac = (tail_win.count("A") + tail_win.count("a")) / max(1, len(tail_win))
    return best_run >= min_len or a_frac >= min_frac


def microhomology_tsd(ref_fa_idx: "pysam.FastaFile", chrom: str, pos: int,
                      ins_seq: str, flank: int = 30,
                      ) -> Optional[Dict[str, Any]]:
    """Detect the target site duplication (TSD) at an insertion breakpoint.

    Looks for the longest k-mer (7-20 bp) occurring immediately 5' AND
    3' of the breakpoint — the duplicated sequence hallmark of
    retrotransposition. VCF POS is the 1-based anchor base; the
    insertion sits between POS and POS+1 (pysam fetch is 0-based
    half-open).

    Parameters
    ----------
    ref_fa_idx : pysam.FastaFile
        Indexed reference FASTA.
    chrom : str
        Chromosome of the insertion.
    pos : int
        1-based VCF POS of the anchor base.
    ins_seq : str
        Inserted sequence (unused; kept for call symmetry).
    flank : int, optional
        Reference bases fetched on each side (default 30).

    Returns
    -------
    dict or None
        ``{"tsd": str, "len": int}`` for the longest duplicated k-mer,
        or None when no 7-20 bp microhomology is found or the fetch
        fails.
    """
    try:
        up = ref_fa_idx.fetch(chrom, max(0, pos - 1 - flank), pos)  # ends at POS
        down = ref_fa_idx.fetch(chrom, pos, pos + 1 + flank)  # starts at POS+1
    except Exception:
        return None
    up, down = up.upper(), down.upper()
    for k in range(20, 6, -1):
        if len(up) >= k and len(down) >= k:
            a = up[-k:]
            b = down[:k]
            if a == b and set(a) != {"N"}:
                return {"tsd": a, "len": k}
    return None


# ------------------------------------------------- family bucket profiles
# A profile is an ORDERED list of buckets; the first bucket whose match
# level hits claims the event (so put specific buckets before generic
# ones, e.g. IAP/ETn before a bare LTR bucket).
#
# Bucket schema (all matching is EXACT, case-insensitive):
#   name       display/bucket name (never used for matching)
#   subfamily  match against the RM repeat NAME (AluY / B1_Mm / IAP1 ...)
#   family     match against the RM class/family suffix ("SINE/Alu" -> Alu;
#              entries without a slash match the whole string)
#   class      match against the RM class/family prefix ("SINE", "LTR" ...)
# Level probing order is subfamily -> family -> class; ONLY the first
# non-empty level list of a bucket participates in matching.
#   catch_all  "classified"   RM hit but no bucket matched
#              "unclassified" no RM hit at all
#   retro_fingerprint  true -> include this bucket in the fingerprint figure
#   zero_slot   false -> omit the bucket from figures when it has zero
#               events (default true: keep the slot, labelled "0")

def parse_te_buckets(spec: Union[str, List[Bucket]]) -> List[Bucket]:
    """Validate the --te-buckets value into an ordered bucket list.

    ONE input form: the bucket list itself (schema documented above),
    transported as a JSON string on the CLI. The flow config carries
    the list natively (Params.sv.mei.te_buckets is the per-species list;
    sv.smk json.dumps it here).

    Parameters
    ----------
    spec : str or list
        JSON string encoding the bucket list, or the already-parsed
        list.

    Returns
    -------
    list of dict
        The validated ordered bucket list.

    Raises
    ------
    ValueError
        If the JSON does not parse, the value is not a non-empty list,
        a bucket lacks ``name``, a bucket has no match level and no
        ``catch_all``, ``catch_all`` is not ``"classified"`` /
        ``"unclassified"``, or a match-level value is not a list of
        strings.
    """
    if isinstance(spec, str):
        try:
            spec = json.loads(spec)
        except json.JSONDecodeError as e:
            raise ValueError(
                f"--te-buckets must be a JSON list of buckets: {e}") from e
    if not isinstance(spec, list) or not spec:
        raise ValueError(
            "TE buckets must be a non-empty list; Params.sv.mei.te_buckets "
            "is empty or missing — regenerate the run config with run.py "
            "(it fills the species list from genome.default)")
    for i, b in enumerate(spec):
        if not isinstance(b, dict) or "name" not in b:
            raise ValueError(f"bucket #{i} must be a dict with a 'name'")
        levels = [k for k in ("subfamily", "family", "class")
                  if b.get(k) not in (None, [])]
        ca = b.get("catch_all")
        if not levels and not ca:
            raise ValueError(
                f"bucket '{b['name']}' has none of subfamily/family/class/"
                "catch_all")
        if ca and ca not in ("classified", "unclassified"):
            raise ValueError(
                f"bucket '{b['name']}': catch_all must be 'classified' or "
                "'unclassified'")
        for k in levels:
            if not (isinstance(b[k], list)
                    and all(isinstance(x, str) for x in b[k])):
                raise ValueError(
                    f"bucket '{b['name']}': '{k}' must be a list of strings")
    return spec


def make_bucket_fn(buckets: List[Bucket]) -> Callable[[str, str], str]:
    """Compile a bucket list into a classifier closure.

    Each bucket contributes at most one match level: the FIRST non-empty
    list among ``subfamily`` / ``family`` / ``class`` (probe order).
    Matching is exact and case-insensitive; events with a RepeatMasker
    hit that no bucket claims go to the ``"classified"`` catch-all
    (default name ``"other"``), events without a hit go to the
    ``"unclassified"`` catch-all.

    Parameters
    ----------
    buckets : list of dict
        Ordered bucket list (see the schema above); first hit wins.

    Returns
    -------
    callable
        ``classify(te_repeat, te_class) -> str`` mapping one event's
        RepeatMasker repeat NAME and class/family field to a bucket
        name.
    """
    compiled = []  # (name, level, {lowercased values})
    catch = {}
    for b in buckets:
        for level in ("subfamily", "family", "class"):
            vals = b.get(level)
            if vals:
                compiled.append((b["name"], level, {v.lower() for v in vals}))
                break  # only the FIRST non-empty level participates
        if b.get("catch_all"):
            catch[b["catch_all"]] = b["name"]
    def classify(rep: str, cls: str) -> str:
        rep = (rep or "").strip()
        cls = (cls or "").strip()
        if not rep or rep.lower() == "none":
            return catch.get("unclassified", "unclassified")
        rep_l, cls_l = rep.lower(), cls.lower()
        family = cls_l.split("/")[-1] if "/" in cls_l else cls_l
        klass = cls_l.split("/")[0] if "/" in cls_l else cls_l
        field = {"subfamily": rep_l, "family": family, "class": klass}
        for name, level, vals in compiled:
            if field[level] in vals:
                return name
        return catch.get("classified", "other")
    return classify


# ------------------------------------------------------------------ annotation
def load_gtf_genes(gtf_path: str) -> List[Tuple[str, int, int, str]]:
    """Load gene spans from a GTF annotation.

    Parameters
    ----------
    gtf_path : str
        GTF path (plain or gzipped).

    Returns
    -------
    list of tuple
        ``(chrom, start, end, gene_name)`` per ``gene`` feature, 1-based
        inclusive coordinates; ``gene_name`` falls back to ``gene_id``.
    """
    genes = []
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            genes.append((f[0], int(f[3]), int(f[4]),
                          attrs.get("gene_name") or attrs.get("gene_id", "")))
    return genes


def load_exons(gtf_path: str) -> List[Tuple[str, int, int]]:
    """Load exon spans from a GTF annotation.

    Parameters
    ----------
    gtf_path : str
        GTF path (plain or gzipped).

    Returns
    -------
    list of tuple
        ``(chrom, start, end)`` per ``exon`` feature, 1-based inclusive
        coordinates.
    """
    exons = []
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            exons.append((f[0], int(f[3]), int(f[4])))
    return exons


def load_te_bed(bed_path: str) -> List[Tuple[str, int, int, str]]:
    """Load reference TE regions from a BED file.

    Parameters
    ----------
    bed_path : str
        BED path (plain or gzipped), e.g. a RepeatMasker TE BED.

    Returns
    -------
    list of tuple
        ``(chrom, start, end, name)`` per record; ``name`` defaults to
        ``"TE"`` when the BED has no name column. Rows with
        non-integer coordinates are skipped.
    """
    tes = []
    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                tes.append((f[0], int(f[1]), int(f[2]),
                            f[3] if len(f) > 3 else "TE"))
            except ValueError:
                continue
    return tes


def annotate_context(pos_rec: Rec, genes: List[Tuple[str, int, int, str]],
                      exons: List[Tuple[str, int, int]],
                      tes: List[Tuple[str, int, int, str]],
                      binsearch: Tuple[List[int], List[int], List[int]],
                      ) -> Tuple[str, str, str]:
    """Annotate the genomic context of an insertion breakpoint.

    Walks back from the bisect insertion point to find the containing
    feature: gene span first (intron/exon), then TE region.

    Parameters
    ----------
    pos_rec : dict
        Record with ``chrom`` and ``pos`` keys (1-based).
    genes : list of tuple
        ``load_gtf_genes`` output, sorted by start.
    exons : list of tuple
        ``load_exons`` output, sorted by start.
    tes : list of tuple
        ``load_te_bed`` output, sorted by start.
    binsearch : tuple of list
        ``(gene_starts, exon_starts, te_starts)`` — the start-column
        indexes used for bisect lookups.

    Returns
    -------
    context : str
        ``"exon"``, ``"intron"`` or ``"intergenic"``.
    gene : str
        Name of the containing gene, or empty.
    te_hit : str
        Name of the overlapping reference TE region, or empty.
    """
    chrom, pos = pos_rec["chrom"], pos_rec["pos"]
    gidx, eidx, tidx = binsearch
    ctx = "intergenic"
    gene = ""
    import bisect
    # genes: any gene span containing pos
    starts = gidx
    i = bisect.bisect_right(starts, pos) - 1
    while i >= 0 and genes_sorted_check(genes, i, chrom):
        c, s, e, gn = genes[i]
        if c == chrom and s <= pos <= e:
            gene = gn
            ctx = "intron"
            break
        i -= 1
    if ctx == "intron":
        estarts = eidx
        j = bisect.bisect_right(estarts, pos) - 1
        while j >= 0:
            c, s, e = exons[j]
            if c == chrom and s <= pos <= e:
                ctx = "exon"
                break
            if c != chrom:
                break
            j -= 1
    te_hit = ""
    tstarts = tidx
    k = bisect.bisect_right(tstarts, pos) - 1
    while k >= 0:
        c, s, e, tn = tes[k]
        if c == chrom and s <= pos <= e:
            te_hit = tn
            break
        if c != chrom:
            break
        k -= 1
    return ctx, gene, te_hit


def genes_sorted_check(genes: List[Tuple[str, int, int, str]], i: int,
                       chrom: str) -> bool:
    """Guard for the backward gene-span walk in `annotate_context`.

    Parameters
    ----------
    genes : list of tuple
        Gene spans (unused; kept for call symmetry).
    i : int
        Current index into the gene list.
    chrom : str
        Chromosome (unused; kept for call symmetry).

    Returns
    -------
    bool
        True while the index is in range.
    """
    return i >= 0


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True,
                    help="SURVIVOR-merged SV VCF (pbsv source, INS with ALT seq)")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--te-lib", default=None,
                    help="TE consensus FASTA for RepeatMasker -lib (legacy path; "
                         "optional when --dfam-h5 is given)")
    ap.add_argument("--dfam-h5", action="append", default=None,
                    help="Dfam FamDB h5 partition file(s); combined with the "
                         "RepeatMasker-bundled root file into a famdb database")
    ap.add_argument("--rm-species", default=None,
                    help="species for the Dfam famdb export / -species search "
                         "(required with --dfam-h5; the flow fills it from the "
                         "genome via node.py GENOME_RM_SPECIES)")
    ap.add_argument("--dfam-mode", choices=["export", "species"], default="export",
                    help="'export': curated consensus export via famdb.py then RM "
                         "-lib (recommended: global best hit + species thresholds); "
                         "'species': RM -species staged search on the famdb dir")
    ap.add_argument("--gtf", default=None, help="gene annotation GTF(.gz)")
    ap.add_argument("--exon-gtf", default=None,
                    help="optional separate exon GTF; default: same as --gtf")
    ap.add_argument("--te-bed", default=None,
                    help="reference rmsk TE regions BED (TE-into-TE annotation)")
    ap.add_argument("--fasta", default=None,
                    help="reference FASTA (indexed) for TSD detection")
    ap.add_argument("--min-len", type=int, default=100,
                    help="minimum inserted-sequence length to classify (default 100)")
    ap.add_argument("--full-length-l1", type=int, default=6000,
                    help="min SVLEN for full-length L1 candidacy (default 6000)")
    ap.add_argument("--full-length-family", default="L1",
                    help="family considered for full-length candidacy "
                         "(default L1; set to e.g. ERVK for other clades)")
    ap.add_argument("--te-buckets", required=True,
                    help="ordered TE bucket list as JSON (the flow fills "
                         "Params.sv.mei.te_buckets from genome.default). "
                         "Buckets are ordered; matching is exact and "
                         "case-insensitive at the first non-empty level of "
                         "subfamily/family/class")
    ap.add_argument("--ctrl-name", default=None,
                    help="label for the control branch (default: ctrl)")
    ap.add_argument("--exp-name", default=None,
                    help="label for the experimental branch (default: exp)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--max-records", type=int, default=0,
                    help="debug: cap number of INS records (0 = all)")
    ap.add_argument("--title", default=None,
                    help="figure/report title (default: parent dir of --outdir)")
    ap.add_argument("--skip-plot", action="store_true",
                    help="skip the figure stage (utils/mei_plot.py)")
    ap.add_argument("--force-rerun", action="store_true",
                    help="recompute even when the input fingerprint matches "
                         "the cached run stamp (default: reuse cache)")

    return ap


def _run_stamp(args: argparse.Namespace, buckets: List[Bucket]) -> Dict[str, Any]:
    """Input fingerprint for whole-run caching (see main).

    Covers everything that can change an output: the VCF (path/size/
    mtime), every analysis parameter, the bucket profile and the
    reference files.
    """
    st = os.stat(args.input)
    return {
        "input": os.path.abspath(args.input),
        "input_size": st.st_size,
        "input_mtime": int(st.st_mtime),
        "min_len": args.min_len,
        "max_records": args.max_records,
        "full_length_l1": args.full_length_l1,
        "full_length_family": args.full_length_family,
        "ctrl_name": args.ctrl_name,
        "exp_name": args.exp_name,
        "title": args.title,
        "te_buckets": buckets,
        "te_lib": args.te_lib,
        "dfam_h5": args.dfam_h5,
        "dfam_mode": args.dfam_mode,
        "rm_species": args.rm_species,
        "gtf": args.gtf,
        "exon_gtf": args.exon_gtf,
        "te_bed": args.te_bed,
        "fasta": args.fasta,
    }


def load_and_label_records(args: argparse.Namespace,
                           ) -> Tuple[List[Rec], str, str]:
    """Load the merged VCF and label each record by SUPP_VEC group.

    Assumes control-first merge order. Custom branch labels keep
    figures/tables readable for any design (e.g. --ctrl-name serum).
    """
    recs, samples = load_merged_vcf(args.input)
    if args.max_records:
        recs = recs[: args.max_records]
    logger.info(f"sample columns: {samples}")
    cb = args.ctrl_name or "ctrl"
    eb = args.exp_name or "exp"
    vec = {"10": f"{cb}_specific", "01": f"{eb}_specific", "11": "shared"}
    for r in recs:
        r["group"] = vec.get(r["supp_vec"], f"other({r['supp_vec']})")
    return recs, cb, eb


def classify_insertions(recs: List[Rec], args: argparse.Namespace,
                        buckets: List[Bucket]) -> None:
    """RepeatMasker classification + family bucket assignment."""
    classify = make_bucket_fn(buckets)
    # ---------------- classification (RepeatMasker) ----------------

    to_classify = [r for r in recs if r["ins_len"] >= args.min_len]
    logger.info(f"classifying {len(to_classify)} INS >= {args.min_len} bp "
        f"(shorter kept as 'short_unclassified')")

    rm_dir = os.path.join(args.outdir, "repeatmasker")
    os.makedirs(rm_dir, exist_ok=True)
    ids = [f"ins{i}" for i in range(len(to_classify))]
    species_lib = args.te_lib
    libdir = None
    rm_species = None
    if args.dfam_h5:
        if not args.rm_species:
            raise SystemExit("--rm-species is required with --dfam-h5 "
                             "(e.g. --rm-species 'Mus musculus')")
        famdb_parent = build_famdb_dir(args.dfam_h5, rm_dir)
        if args.dfam_mode == "export":
            species_lib = export_species_lib(
                famdb_parent, args.rm_species,
                os.path.join(rm_dir, "famdb_export.fa"))
        else:
            species_lib = None
            libdir, rm_species = famdb_parent, args.rm_species
    elif not species_lib:
        raise SystemExit("one of --te-lib / --dfam-h5 is required")
    best = run_repeatmasker([r["seq"] for r in to_classify], ids, rm_dir,
                            species_lib, threads=args.threads,
                            libdir=libdir, rm_species=rm_species,
                            config_stamp={"te_lib": args.te_lib,
                                          "dfam_h5": args.dfam_h5,
                                          "dfam_mode": args.dfam_mode,
                                          "rm_species": args.rm_species},
                            force=args.force_rerun)

    for r, sid in zip(to_classify, ids):
        h = best.get(sid)
        r["rm_id"] = sid  # ins.fa query id: joins mei_events.tsv to ins.fa.out rows
        if h:
            cov = (h["end"] - h["begin"] + 1) / r["ins_len"]
            r["te_repeat"] = h["repeat"]
            r["te_class"] = h["class"]
            r["te_family"] = h["family"]
            r["te_cov"] = round(cov, 3)
            r["te_strand"] = h["strand"]
            r["te_sw"] = h["sw"]
        else:
            r["te_repeat"] = "none"
            r["te_class"] = ""
            r["te_family"] = ""
            r["te_cov"] = 0.0
            r["te_strand"] = ""
            r["te_sw"] = 0

    # family bucket: profile-driven (see the bucket schema / --te-buckets)
    for r in recs:
        if r.get("te_repeat"):
            r["family"] = classify(r["te_repeat"], r["te_class"])
        else:
            r["family"] = ("short_unclassified"
                           if r["ins_len"] < args.min_len else "unclassified")

    return


def annotate_fingerprints(recs: List[Rec], args: argparse.Namespace,
                          buckets: List[Bucket]) -> int:
    """poly(A) tail / TSD fingerprints; returns the poly(A) hit count."""
    fingerprint_fams = {b["name"] for b in buckets
                        if b.get("retro_fingerprint")}
    # ---------------- fingerprints ----------------
    fa_idx = None
    if args.fasta:
        try:
            import pysam
            fa_idx = pysam.FastaFile(args.fasta)
        except ImportError:
            logger.info("pysam not available; TSD detection skipped")
    n_polya = 0
    for r in recs:
        r["polyA"] = bool(poly_a_tail(r["seq"])) if r["ins_len"] >= 200 else None
        if r["polyA"]:
            n_polya += 1
        if fa_idx is not None and r["family"] in fingerprint_fams:
            r["tsd"] = microhomology_tsd(fa_idx, r["chrom"], r["pos"], r["seq"])
        else:
            r["tsd"] = None

    return n_polya


def annotate_genomic_context(recs: List[Rec], args: argparse.Namespace) -> None:
    """exon/intron/intergenic and TE-into-TE context fields."""
    # ---------------- context annotation ----------------
    genes = exons = tes = []
    gidx = eidx = tidx = []
    if args.gtf:
        genes = sorted(load_gtf_genes(args.gtf), key=lambda g: (g[0], g[1]))
        exons = sorted(load_exons(args.exon_gtf or args.gtf),
                       key=lambda e: (e[0], e[1]))
        gidx = [g[1] for g in genes]
        eidx = [e[1] for e in exons]
    if args.te_bed:
        tes = sorted(load_te_bed(args.te_bed), key=lambda t: (t[0], t[1]))
        tidx = [t[1] for t in tes]
    for r in recs:
        if genes:
            ctx, gene, te_hit = annotate_context(r, genes, exons, tes,
                                                 (gidx, eidx, tidx))
            r["context"] = ctx
            r["gene"] = gene
            r["in_te_region"] = te_hit
        else:
            r["context"] = r["gene"] = r["in_te_region"] = ""

def mark_full_length(recs: List[Rec], args: argparse.Namespace) -> List[Rec]:
    """Flag full-length candidates and return them."""
    # full-length L1 candidates
    fl_fam = args.full_length_family
    fl = [r for r in recs
          if r["family"] == fl_fam and r["svlen"] >= args.full_length_l1]
    for r in fl:
        r["full_length_L1"] = True
    for r in recs:
        r.setdefault("full_length_L1", False)
    return fl


def _bh_fdr(pvals: List[float]) -> List[float]:
    """Benjamini-Hochberg FDR with monotonicity enforcement."""
    n = len(pvals)
    order = sorted(range(n), key=pvals.__getitem__)
    adj = [1.0] * n
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        prev = min(prev, pvals[idx] * n / (i + 1))
        adj[idx] = prev
    return adj


def _side_metrics(sub: "pd.DataFrame", branch: str) -> Dict[str, float]:
    """Per-branch metrics for one family: its specific calls + shared ones.

    Shared calls count toward both branches, so each branch sees every
    insertion that is present in its biological state.
    """
    side = sub[sub["group"].isin([f"{branch}_specific", "shared"])]
    polya = side["polyA"].dropna().astype(bool)
    in_te = side["in_te_region"].fillna("").astype(str) != ""
    return {
        "polyA_pct": round(100.0 * float(polya.mean()), 1) if len(polya) else float("nan"),
        "in_te_pct": round(100.0 * float(in_te.mean()), 1) if len(side) else float("nan"),
        "median_ins_len": float(side["ins_len"].median()) if len(side) else float("nan"),
    }


def diff_by_group(recs: List[Rec], cb: str, eb: str, key: str = "family",
                  ) -> "pd.DataFrame":
    """Rank families (or individual TE repeats) by between-branch difference.

    ``key`` selects the grouping column: ``"family"`` for the TE buckets,
    ``"te_repeat"`` for individual RepeatMasker repeat names (records
    without a repeat hit are dropped at that level). Per group: branch-
    specific call counts, Haldane-corrected log2 fold change of the exp vs
    ctrl call share, two-sided Fisher exact p on the call-count table and
    BH-FDR, plus per-branch poly(A) rate, in-TE rate and median insertion
    length. Rows are sorted by FDR then p.

    The ``focus`` column marks triage candidates (NOT significance):
    |log2FC| >= 1 (2x share), p < 0.05 and at least 10 specific calls in
    total.

    Notes
    -----
    No biological replicates: the p values describe CALL-COUNT imbalance
    (coverage and calling noise included), not population statistics —
    read this as an exploratory ranking of "most different", not as
    significance.
    """
    from scipy.stats import fisher_exact

    df = pd.DataFrame([{k: r.get(k) for k in
                        ("family", "te_repeat", "group", "ins_len", "polyA",
                         "in_te_region")} for r in recs])
    if key == "te_repeat":
        df = df[~df["te_repeat"].fillna("").astype(str).str.lower()
                  .isin(["", "none"])]
    ctrl, exp = f"{cb}_specific", f"{eb}_specific"
    n_ctrl = int((df["group"] == ctrl).sum())
    n_exp = int((df["group"] == exp).sum())

    rows: List[Dict[str, Any]] = []
    for name, sub in df.groupby(key, sort=False):
        a = int((sub["group"] == ctrl).sum())
        b = int((sub["group"] == exp).sum())
        if a + b == 0:
            continue
        _, p = fisher_exact([[a, max(n_ctrl - a, 0)],
                             [b, max(n_exp - b, 0)]])
        cm = _side_metrics(sub, cb)
        em = _side_metrics(sub, eb)
        rows.append({
            key: name,
            "n_ctrl_specific": a,
            "n_exp_specific": b,
            "n_shared": int((sub["group"] == "shared").sum()),
            "log2fc_exp_vs_ctrl": math.log2(((b + 0.5) / (n_exp + 1))
                                            / ((a + 0.5) / (n_ctrl + 1))),
            "p_value": float(p),
            "ctrl_polyA_pct": cm["polyA_pct"], "exp_polyA_pct": em["polyA_pct"],
            "ctrl_in_te_pct": cm["in_te_pct"], "exp_in_te_pct": em["in_te_pct"],
            "ctrl_median_ins_len": cm["median_ins_len"],
            "exp_median_ins_len": em["median_ins_len"],
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr"] = _bh_fdr(out["p_value"].tolist())
    # focus rule (triage, not significance): 2x share, p < 0.05, and
    # enough specific calls to trust the ratio
    out["focus"] = ((out["log2fc_exp_vs_ctrl"].abs() >= 1.0)
                    & (out["p_value"] < 0.05)
                    & (out["n_ctrl_specific"] + out["n_exp_specific"] >= 10))
    return out.sort_values(["fdr", "p_value"], kind="mergesort").reset_index(drop=True)


def write_outputs(recs: List[Rec], args: argparse.Namespace,
                  buckets: List[Bucket], cb: str, eb: str, n_polya: int,
                  ) -> Dict[str, Any]:
    """Write the tsv tables and mei_summary.json; return the summary."""
    fl_fam = args.full_length_family
    # ---------------- outputs ----------------
    cols = ["chrom", "pos", "id", "svlen", "ins_len", "supp_vec", "group",
            "family", "te_repeat", "te_class", "te_family", "te_cov",
            "te_strand", "polyA", "tsd", "context", "gene", "in_te_region",
            "full_length_L1", "rm_id"]
    df = pd.DataFrame([{c: r.get(c) for c in cols} for r in recs])
    df["tsd"] = df["tsd"].apply(lambda t: t["tsd"] if t else "")
    df.to_csv(os.path.join(args.outdir, "mei_events.tsv"), sep="\t", index=False)

    fam_group = df.groupby(["family", "group"]).size().unstack(fill_value=0)
    fam_group.to_csv(os.path.join(args.outdir, "mei_family_by_group.tsv"),
                     sep="\t")

    # context table: TE buckets only (exclude repeat-noise classes —
    # simple/low-complexity/satellite — so the table stays interpretable)
    noise_classes = {"Simple_repeat", "Satellite", "Low_complexity"}
    te_bucket_names = [b["name"] for b in buckets
                       if not b.get("catch_all")
                       and not (set(b.get("class") or []) & noise_classes)]
    ctx_tbl = (df[df["family"].isin(te_bucket_names)]
               .groupby(["family", "group", "context"]).size()
               .unstack(fill_value=0))
    ctx_tbl.to_csv(os.path.join(args.outdir, "mei_context.tsv"), sep="\t")

    fl_df = df[df["full_length_L1"]]
    fl_df.to_csv(os.path.join(args.outdir, "mei_full_length_L1.tsv"),
                 sep="\t", index=False)

    diff_df = diff_by_group(recs, cb, eb, key="family")
    diff_df.to_csv(os.path.join(args.outdir, "mei_diff.tsv"), sep="\t",
                   index=False)
    diff_rep_df = diff_by_group(recs, cb, eb, key="te_repeat")
    diff_rep_df.to_csv(os.path.join(args.outdir, "mei_diff_by_repeat.tsv"),
                       sep="\t", index=False)

    summary = {
        "input": args.input,
        "branches": {"control": cb, "experiment": eb},
        "te_buckets": buckets,
        "full_length_family": fl_fam,
        "n_ins_total": len(recs),
        "n_classified": int((df["family"] != "unclassified").sum()),
        "n_with_polyA": int(n_polya),
        "full_length_L1_candidates": len(fl_df),
        "family_counts_overall": df["family"].value_counts().to_dict(),
        "family_by_group": {fam: row.to_dict()
                            for fam, row in fam_group.iterrows()},
        "context_counts": df["context"].value_counts().to_dict(),
        # exploratory ranking only — no biological replicates
        "diff_top": json.loads(diff_df.head(10).to_json(orient="records")),
        "diff_top_repeat": json.loads(
            diff_rep_df.head(10).to_json(orient="records")),
    }
    with open(os.path.join(args.outdir, "mei_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def run_figures(args: argparse.Namespace) -> None:
    """Stage 2: figures via utils/mei_plot.py."""
    if not args.skip_plot:
        utils_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils")
        title = args.title or os.path.basename(
            os.path.dirname(os.path.abspath(args.outdir.rstrip("/"))))
        plot_cmd = [sys.executable, os.path.join(utils_dir, "mei_plot.py"),
                    "--mei-dir", args.outdir, "--title", title]
        if args.ctrl_name:
            plot_cmd += ["--ctrl-name", args.ctrl_name]
        if args.exp_name:
            plot_cmd += ["--exp-name", args.exp_name]
        logger.info("figures -> " + " ".join(plot_cmd))
        subprocess.run(plot_cmd, check=True)


def main() -> None:
    args = build_arg_parser().parse_args()
    buckets = parse_te_buckets(args.te_buckets)
    logger.info(f"TE buckets ({len(buckets)}): "
        + ", ".join(b["name"] for b in buckets))
    os.makedirs(args.outdir, exist_ok=True)

    # whole-run cache: identical inputs + existing outputs -> skip
    run_stamp_path = os.path.join(args.outdir, "mei_run.stamp.json")
    run_stamp = json.dumps(_run_stamp(args, buckets), sort_keys=True)
    if (not args.force_rerun
            and os.path.exists(run_stamp_path)
            and os.path.exists(os.path.join(args.outdir, "mei_events.tsv"))
            and os.path.exists(os.path.join(args.outdir, "mei_summary.json"))):
        with open(run_stamp_path) as fh:
            if fh.read().strip() == run_stamp:
                logger.info("inputs unchanged — reusing cached outputs "
                            "(--force-rerun to recompute)")
                return

    recs, cb, eb = load_and_label_records(args)
    classify_insertions(recs, args, buckets)
    n_polya = annotate_fingerprints(recs, args, buckets)
    annotate_genomic_context(recs, args)
    mark_full_length(recs, args)
    summary = write_outputs(recs, args, buckets, cb, eb, n_polya)
    run_figures(args)
    with open(run_stamp_path, "w") as fh:
        fh.write(run_stamp)
    logger.info(f"done. outputs in {args.outdir}")
    logger.info(f"families: {summary['family_counts_overall']}")


if __name__ == "__main__":
    main()

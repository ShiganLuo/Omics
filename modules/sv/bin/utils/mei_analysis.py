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
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

import pandas as pd


def log(msg):
    print(f"[MEI] {msg}", flush=True)


# ---------------------------------------------------------------- VCF parsing
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_alt_insertion(ref, alt):
    """Return the pure inserted sequence for a simple insertion ALT.

    Handles symbolic <INS> (None), plain-sequence ALTs and left-aligned
    padded ALTs sharing a prefix with REF (most common from pbsv via
    SURVIVOR: REF=A, ALT=A<inserted>).
    """
    if alt is None or alt.startswith("<") or alt == ".":
        return None
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


def parse_info(info):
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def load_merged_vcf(path):
    """Yield per-INS records with inserted sequence and per-sample DR."""
    opener = gzip.open if path.endswith(".gz") else open
    recs = []
    dropped_symbolic = 0
    with opener(path, "rt") as fh:
        samples = None
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
    log(f"parsed INS records with sequence: {len(recs)} "
        f"(dropped symbolic/no-seq: {dropped_symbolic})")
    return recs, samples


# ------------------------------------------------------- RepeatMasker classify
def run_repeatmasker(seqs, ids, workdir, species_lib, threads=4):
    """Write FASTA, run RepeatMasker, parse .out -> {id: best hit dict}.

    Reuses an existing ins.fa.out when it covers all ids (cache hit after a
    parser-fix rerun; the FASTA is rewritten only when absent)."""
    fa = os.path.join(workdir, "ins.fa")
    out = os.path.join(workdir, "ins.fa.out")
    if not os.path.exists(out):
        with open(fa, "w") as fh:
            for i, (sid, s) in enumerate(zip(ids, seqs)):
                fh.write(f">{sid}\n")
                for j in range(0, len(s), 60):
                    fh.write(s[j:j + 60] + "\n")
        cmd = ["RepeatMasker", "-no_is", "-s", "-pa",
               str(threads), "-lib", species_lib, "-dir", workdir, fa]
        log(f"running RepeatMasker on {len(seqs)} sequences (this can take a while)")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            log(f"RepeatMasker failed: {r.stderr[:500]}")
            return {}
        if not os.path.exists(out):
            log("RepeatMasker .out missing")
            return {}
    else:
        log(f"reusing cached {out}")
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
def poly_a_tail(seq, min_len=12, min_frac=0.9):
    """Detect A-rich 3' tail (L1/Alu/SVA hallmark) allowing interruptions."""
    tail_win = seq[-60:] if len(seq) >= 60 else seq
    best_run = 0
    cur = 0
    for ch in tail_win:
        cur = cur + 1 if ch in "Aa" else 0
        best_run = max(best_run, cur)
    a_frac = (tail_win.count("A") + tail_win.count("a")) / max(1, len(tail_win))
    return best_run >= min_len or a_frac >= min_frac


def microhomology_tsd(ref_fa_idx, chrom, pos, ins_seq, flank=30):
    """TSD: duplicated k-mer immediately 5' AND 3' of the breakpoint, 7-20bp.

    VCF POS is the 1-based anchor base; the insertion sits between POS and
    POS+1. pysam fetch is 0-based half-open.
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


# ------------------------------------------------------------------ annotation
def load_gtf_genes(gtf_path):
    """gene (chrom, start, end, gene_name) tuples; GTF or sorted bgzipped."""
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


def load_exons(gtf_path):
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


def load_te_bed(bed_path):
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


def annotate_context(pos_rec, genes, exons, tes, binsearch):
    """exon / intron / intergenic (+ TE region flag)."""
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


def genes_sorted_check(genes, i, chrom):
    return i >= 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True,
                    help="SURVIVOR-merged SV VCF (pbsv source, INS with ALT seq)")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--te-lib", required=True,
                    help="TE consensus FASTA for RepeatMasker -lib (e.g. rhesus_te_lib.fa)")
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
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--max-records", type=int, default=0,
                    help="debug: cap number of INS records (0 = all)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    recs, samples = load_merged_vcf(args.input)
    if args.max_records:
        recs = recs[: args.max_records]
    log(f"sample columns: {samples}")

    # group labelling by SUPP_VEC: assumes control-first merge order
    vec = {"10": "ctrl_specific", "01": "exp_specific", "11": "shared"}
    for r in recs:
        r["group"] = vec.get(r["supp_vec"], f"other({r['supp_vec']})")

    # ---------------- classification (RepeatMasker) ----------------
    import bisect as _bs

    to_classify = [r for r in recs if r["ins_len"] >= args.min_len]
    log(f"classifying {len(to_classify)} INS >= {args.min_len} bp "
        f"(shorter kept as 'short_unclassified')")

    rm_dir = os.path.join(args.outdir, "repeatmasker")
    os.makedirs(rm_dir, exist_ok=True)
    ids = [f"ins{i}" for i in range(len(to_classify))]
    best = run_repeatmasker([r["seq"] for r in to_classify], ids, rm_dir,
                            args.te_lib, threads=args.threads)

    for r, sid in zip(to_classify, ids):
        h = best.get(sid)
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

    # family bucket: SINE/Alu, LINE/L1, SVA, LTR, DNA, simple_repeat, none
    def bucket(r):
        rep = (r["te_repeat"] or "").lower()
        cls = (r["te_class"] or "").lower()
        if not rep or rep == "none":
            return "unclassified"
        if "alu" in rep or "alu" in cls:
            return "Alu"
        if rep.startswith("l1") or "line" in cls or rep.startswith("lre"):
            return "L1"
        if "sva" in rep:
            return "SVA"
        if "ltr" in cls or rep.startswith("erv") or "maa" in rep:
            return "LTR"
        if "dna" in cls:
            return "DNA"
        if "(" in rep and ")" in rep:
            return "simple_repeat"
        return "other"
    for r in recs:
        r["family"] = bucket(r) if r.get("te_repeat") else (
            "short_unclassified" if r["ins_len"] < args.min_len else "unclassified")

    # ---------------- fingerprints ----------------
    fa_idx = None
    if args.fasta:
        try:
            import pysam
            fa_idx = pysam.FastaFile(args.fasta)
        except ImportError:
            log("pysam not available; TSD detection skipped")
    n_polya = 0
    for r in recs:
        r["polyA"] = bool(poly_a_tail(r["seq"])) if r["ins_len"] >= 200 else None
        if r["polyA"]:
            n_polya += 1
        if fa_idx is not None and r["family"] in ("L1", "Alu", "SVA"):
            r["tsd"] = microhomology_tsd(fa_idx, r["chrom"], r["pos"], r["seq"])
        else:
            r["tsd"] = None

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

    # full-length L1 candidates
    fl = [r for r in recs
          if r["family"] == "L1" and r["svlen"] >= args.full_length_l1]
    for r in fl:
        r["full_length_L1"] = True
    for r in recs:
        r.setdefault("full_length_L1", False)

    # ---------------- outputs ----------------
    cols = ["chrom", "pos", "id", "svlen", "ins_len", "supp_vec", "group",
            "family", "te_repeat", "te_class", "te_family", "te_cov",
            "te_strand", "polyA", "tsd", "context", "gene", "in_te_region",
            "full_length_L1"]
    df = pd.DataFrame([{c: r.get(c) for c in cols} for r in recs])
    df["tsd"] = df["tsd"].apply(lambda t: t["tsd"] if t else "")
    df.to_csv(os.path.join(args.outdir, "mei_events.tsv"), sep="\t", index=False)

    fam_group = df.groupby(["family", "group"]).size().unstack(fill_value=0)
    fam_group.to_csv(os.path.join(args.outdir, "mei_family_by_group.tsv"),
                     sep="\t")

    ctx_tbl = (df[df["family"].isin(["Alu", "L1", "SVA", "LTR", "DNA"])]
               .groupby(["family", "group", "context"]).size()
               .unstack(fill_value=0))
    ctx_tbl.to_csv(os.path.join(args.outdir, "mei_context.tsv"), sep="\t")

    fl_df = df[df["full_length_L1"]]
    fl_df.to_csv(os.path.join(args.outdir, "mei_full_length_L1.tsv"),
                 sep="\t", index=False)

    summary = {
        "input": args.input,
        "n_ins_total": len(recs),
        "n_classified": int((df["family"] != "unclassified").sum()
                            + (df["family"] == "short_unclassified").sum() * 0),
        "n_with_polyA": int(n_polya),
        "full_length_L1_candidates": len(fl_df),
        "family_counts_overall": df["family"].value_counts().to_dict(),
        "family_by_group": {fam: row.to_dict()
                            for fam, row in fam_group.iterrows()},
        "context_counts": df["context"].value_counts().to_dict(),
    }
    with open(os.path.join(args.outdir, "mei_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    log(f"done. outputs in {args.outdir}")
    log(f"families: {summary['family_counts_overall']}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Split a sawfish joint-call VCF (VCF 4.4, merged SV+CNV) into CNV and SV sets.

Record classification follows the sawfish VCF 4.4 ``SVCLAIM`` convention
(user guide "SV and CNV types"):

- CNV set: ``SVCLAIM`` contains ``D`` (depth-segmentation support, including
  merged ``DJ`` records) or ``SVTYPE=CNV`` (multi-sample gain/loss records).
- SV set (only with ``--keep-sv``): ``SVCLAIM`` contains ``J`` (breakpoint /
  junction support) or ``SVTYPE=BND``, plus records without ``SVCLAIM``.
  A ``DJ`` record is emitted in BOTH sets (it is both an SV and a CNV).

Outputs written under --out-dir (bgzipped + indexed):
  sawfish.cnv.vcf.gz(.csi)   - CNV records
  sawfish.sv.vcf.gz(.csi)    - SV records (only with --keep-sv)
  sawfish_cnv_summary.tsv    - record counts per category
"""

import argparse
import gzip
import os
import subprocess
import sys
import tempfile


def open_vcf(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info):
    d = {}
    for kv in info.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def classify(info):
    """Return (is_cnv, is_sv) for one record's INFO dict."""
    claim = str(info.get("SVCLAIM", "")).upper()
    svtype = str(info.get("SVTYPE", "")).upper()
    is_cnv = ("D" in claim) or svtype == "CNV"
    # 无 SVCLAIM 的记录按断点证据兜底(如 BND),但纯 CNV 记录除外
    is_sv = ("J" in claim) or svtype == "BND" or \
        ("SVCLAIM" not in info and svtype != "CNV")
    return is_cnv, is_sv


def bgzip_and_index(plain_path, gz_path, bgzip, bcftools):
    with open(gz_path, "wb") as out_fh:
        subprocess.run([bgzip, "-c", plain_path], stdout=out_fh, check=True)
    subprocess.run([bcftools, "index", "-f", gz_path], check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True,
                    help="sawfish joint-call output genotyped.sv.vcf.gz")
    ap.add_argument("--out-dir", required=True,
                    help="output directory for the split VCFs and summary")
    ap.add_argument("--keep-sv", action="store_true",
                    help="also emit the breakpoint-SV set (sawfish.sv.vcf.gz)")
    ap.add_argument("--bgzip", default="bgzip", help="bgzip executable")
    ap.add_argument("--bcftools", default="bcftools", help="bcftools executable")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    cnv_gz = os.path.join(args.out_dir, "sawfish.cnv.vcf.gz")
    sv_gz = os.path.join(args.out_dir, "sawfish.sv.vcf.gz")
    summary = os.path.join(args.out_dir, "sawfish_cnv_summary.tsv")

    counts = {"total": 0, "cnv": 0, "sv": 0, "both": 0,
              "DEL": 0, "DUP": 0, "CNV": 0, "BND": 0, "other": 0}

    tmpdir = tempfile.mkdtemp(prefix="sawfish_split_")
    cnv_plain = os.path.join(tmpdir, "cnv.vcf")
    sv_plain = os.path.join(tmpdir, "sv.vcf")

    with open_vcf(args.vcf) as fin, \
            open(cnv_plain, "w") as fcnv, \
            open(sv_plain, "w") as fsv:
        for line in fin:
            if line.startswith("#"):
                fcnv.write(line)
                if args.keep_sv:
                    fsv.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            counts["total"] += 1
            info = parse_info(fields[7])
            svtype = str(info.get("SVTYPE", "")).upper()
            counts[svtype if svtype in ("DEL", "DUP", "CNV", "BND") else "other"] += 1
            is_cnv, is_sv = classify(info)
            if is_cnv:
                counts["cnv"] += 1
                fcnv.write(line)
            if is_sv:
                counts["sv"] += 1
                if args.keep_sv:
                    fsv.write(line)
            if is_cnv and is_sv:
                counts["both"] += 1

    bgzip_and_index(cnv_plain, cnv_gz, args.bgzip, args.bcftools)
    if args.keep_sv:
        bgzip_and_index(sv_plain, sv_gz, args.bgzip, args.bcftools)

    with open(summary, "w") as fh:
        fh.write("category\tcount\n")
        for key in ("total", "cnv", "sv", "both",
                    "DEL", "DUP", "CNV", "BND", "other"):
            fh.write(f"{key}\t{counts[key]}\n")

    print(f"sawfish_split_vcf: {counts['total']} records -> "
          f"{counts['cnv']} CNV, {counts['sv']} SV"
          f"{'' if args.keep_sv else ' (not emitted)'}", file=sys.stderr)


if __name__ == "__main__":
    main()

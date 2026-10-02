#!/usr/bin/env python3
"""Convert a pbsv SV VCF into a Decoil-compatible (sniffles1-style) VCF.

Decoil 2.x ``reconstruct --sv-caller sniffles1`` requires fields that pbsv
does not emit:

- INFO/STRANDS   (read as ``STRANDS[0]`` by ``transform_record_sniffles1``)
- INFO/CHR2      (read unconditionally by ``pass_filter``)
- INFO/END       (missing on pbsv BND records; read by ``parsevcf``)
- FORMAT/DR, FORMAT/DV (read by ``pass_filter``/``parsevcf``; pbsv has AD)

BND orientation is derived from the VCF breakend ALT notation:

- ``t[p[`` -> ``+-``   (POS side keeps left, mate side keeps right)
- ``t]p]`` -> ``++``   (both sides keep left)
- ``]p]t`` -> ``-+``   (POS side keeps right, mate side keeps left)
- ``[p[t`` -> ``--``   (both sides keep right)

Non-BND records get a placeholder ``STRANDS=+-``: Decoil's ``addsv`` uses
fixed orientations for DEL/INS/DUP/INV and never reads STRANDS for them.

Input may be bgzip-compressed (.vcf.gz) or plain text; output is always
plain text because Decoil opens VCFs with a plain ``open()``.

Usage:
    pbsv2decoil.py <in.vcf[.gz]> <out.vcf>
"""

from __future__ import annotations

import gzip
import re
import sys

BND_ALT_RE = re.compile(
    r"^(?:(?P<pre>[A-Za-z.*]+)(?P<br1>[\[\]])(?P<mate1>[^\[\]:]+):(?P<pos1>\d+)(?P=br1)"
    r"|(?P<br2>[\[\]])(?P<mate2>[^\[\]:]+):(?P<pos2>\d+)(?P=br2)(?P<post>[A-Za-z.*]+))$"
)

# breakend ALT form -> Decoil strand code (POS side, mate side)
BND_STRANDS: dict[tuple[str, str], str] = {
    ("[", "right"): "+-",   # t[p[  — mate piece right of p joined after t
    ("]", "right"): "++",   # t]p]  — mate piece left of p joined after t
    ("]", "left"): "-+",    # ]p]t  — POS right side kept, mate left kept
    ("[", "left"): "--",    # [p[t  — POS right side kept, mate right kept
}

INFO_STRANDS_HDR = (
    '##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the SV">'
)
INFO_CHR2_HDR = (
    '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosomal position of the other breakpoint">'
)
FMT_DR_HDR = (
    '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Read depth supporting the reference allele">'
)
FMT_DV_HDR = (
    '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Read depth supporting the variant allele">'
)


def open_maybe_gzip(path: str):
    """Open plain or bgzip/gzip VCF for text reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_bnd_alt(alt: str) -> tuple[str, int, str, str]:
    """Extract (mate_chrom, mate_pos, orientation, side) from a breakend ALT.

    Parameters
    ----------
    alt : str
        VCF ALT allele for a BND record.

    Returns
    -------
    tuple
        (mate_chrom, mate_pos, bracket, side) where bracket is '[' or ']'
        and side is 'left' if the ref base is at the ALT start ('t[p['),
        'right' if at the ALT end (']p]t').

    Raises
    ------
    ValueError
        If ALT does not match breakend notation.
    """
    m = BND_ALT_RE.match(alt)
    if not m:
        raise ValueError(f"unparseable breakend ALT: {alt}")
    if m.group("pre") is not None:
        return m.group("mate1"), int(m.group("pos1")), m.group("br1"), "right"
    return m.group("mate2"), int(m.group("pos2")), m.group("br2"), "left"


def convert_info(info: str, svtype: str, chrom: str, alt: str) -> str:
    """Add CHR2 / STRANDS / (BND-only) END to a pbsv INFO column.

    Existing keys are preserved and never overwritten.
    """
    fields = [f for f in info.split(";") if f]
    keys = {f.split("=", 1)[0] for f in fields}
    add: list[str] = []

    if "CHR2" not in keys:
        if svtype == "BND":
            mate_chrom, _, _, _ = parse_bnd_alt(alt)
            add.append(f"CHR2={mate_chrom}")
        else:
            add.append(f"CHR2={chrom}")

    if "END" not in keys and svtype == "BND":
        _, mate_pos, _, _ = parse_bnd_alt(alt)
        add.append(f"END={mate_pos}")

    if "STRANDS" not in keys:
        if svtype == "BND":
            _, _, bracket, side = parse_bnd_alt(alt)
            add.append(f"STRANDS={BND_STRANDS[(bracket, side)]}")
        else:
            add.append("STRANDS=+-")

    return ";".join(fields + add) if add else info


def convert_format(fmt: str, sample: str) -> tuple[str, str]:
    """Replace AD with DR/DV in the FORMAT/sample columns.

    pbsv FORMAT is GT:AD:DP; Decoil reads GT/DR/DV (and DP for downsampling).
    DR = AD[ref], DV = AD[alt]. Records that already carry DR/DV pass through.
    """
    keys = fmt.split(":")
    if "DR" in keys and "DV" in keys:
        return fmt, sample
    if "AD" not in keys:
        raise ValueError(f"no AD in FORMAT {fmt}; cannot derive DR/DV")
    vals = sample.split(":")
    kv = dict(zip(keys, vals))
    ref_count, alt_count = kv["AD"].split(",")[:2]
    new_keys = ["GT", "DR", "DV"] + [k for k in keys if k not in ("GT", "AD")]
    kv["DR"], kv["DV"] = ref_count, alt_count
    out_keys = [k for k in new_keys if k in kv or k == "GT"]
    return ":".join(out_keys), ":".join(kv[k] for k in out_keys)


def main(in_path: str, out_path: str) -> None:
    """Convert ``in_path`` to Decoil format at ``out_path``.

    Parameters
    ----------
    in_path : str
        pbsv VCF, plain or bgzip-compressed.
    out_path : str
        Destination plain-text VCF.
    """
    n_records = 0
    with open_maybe_gzip(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("##"):
                fin_header_cache.add(line.split("ID=", 1)[1].split(",", 1)[0] if "ID=" in line else "")
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                for hdr in (INFO_STRANDS_HDR, INFO_CHR2_HDR, FMT_DR_HDR, FMT_DV_HDR):
                    hdr_id = hdr.split("ID=", 1)[1].split(",", 1)[0]
                    if hdr_id not in fin_header_cache:
                        fout.write(hdr + "\n")
                fout.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            chrom, alt, info, fmt, sample = cols[0], cols[4], cols[7], cols[8], cols[9]
            svtype = next(
                (f[7:] for f in info.split(";") if f.startswith("SVTYPE=")), None
            )
            if svtype is None:
                raise ValueError(f"record without SVTYPE: {line[:120]}")
            cols[7] = convert_info(info, svtype, chrom, alt)
            cols[8], cols[9] = convert_format(fmt, sample)
            fout.write("\t".join(cols) + "\n")
            n_records += 1

    print(f"[pbsv2decoil] {n_records} records -> {out_path}", file=sys.stderr)


fin_header_cache: set[str] = set()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])

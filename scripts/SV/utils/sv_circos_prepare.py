import re
from pathlib import Path
from typing import Dict, List, Tuple
import pysam
from pyfaidx import Fasta
import subprocess


class SVCircosPrepare():
    """
    SV Circos preparation utility class.
    """
    main_chrom = ["chr" + c for c in range(1, 20)] + ['X', 'Y', 'M']
    def __init__(self,vcf_path: str, fasta_path: str, outdir: str):
        self.vcf_path = vcf_path
        self.fasta_path = fasta_path
        self.outdir = outdir    

    def open_vcf_autotabix(self) -> pysam.VariantFile:
        """
        Open VCF with the following behavior:

        - .vcf:
            * read directly (no index, no compression)
        - .vcf.gz:
            * if index exists -> use
            * if index missing -> create tabix index
            * NEVER recompress the file
        """
        vcf_path = self.vcf_path
        vcf_path = Path(vcf_path)

        if not vcf_path.exists():
            raise FileNotFoundError(vcf_path)

        # ---------- case 1: plain VCF ----------
        if vcf_path.suffix == ".vcf":
            return pysam.VariantFile(str(vcf_path))

        # ---------- case 2: compressed VCF ----------
        if vcf_path.suffixes[-2:] == [".vcf", ".gz"]:
            tbi = vcf_path.with_suffix(".vcf.gz.tbi")
            csi = vcf_path.with_suffix(".vcf.gz.csi")

            if not tbi.exists() and not csi.exists():
                # create index only
                try:
                    subprocess.run(
                        ["tabix", "-p", "vcf", str(vcf_path)],
                        check=True
                    )
                except subprocess.CalledProcessError:
                    raise RuntimeError(
                        f"Failed to create tabix index for {vcf_path}\n"
                        f"Make sure this file is bgzip-compressed VCF."
                    )

            return pysam.VariantFile(str(vcf_path))

        raise ValueError(
            f"Unsupported VCF format: {vcf_path} "
            f"(supported: .vcf / .vcf.gz)"
        )

    def parse_fasta_chrom_sizes(self) -> Dict[str, int]:
        """
        Parse FASTA and return chromosome sizes.

        Returns
        -------
        dict: {chrom: length}
        """
        fasta_path = self.fasta_path
        fa = Fasta(fasta_path)
        chrom_sizes = {chrom: len(fa[chrom]) for chrom in fa.keys()}
        return chrom_sizes


    @staticmethod
    def write_circos_karyotype(
        chrom_sizes: Dict[str, int],
        out_file: str,
    ):
        """
        Write Circos karyotype file.
        """
        with open(out_file, "w") as f:
            for chrom, length in chrom_sizes.items():
                f.write(f"{chrom} 0 {length}\n")


    @staticmethod
    def parse_bnd_alt(alt: str) -> Tuple[str, int]:
        """
        Parse BND ALT field to extract partner chromosome and position.

        Example ALT:
            N]chr7:80000]
            [chr2:12345[N

        Returns
        -------
        (chrom, pos)
        """
        m = re.search(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]', alt)
        if not m:
            return None, None
        return m.group(1), int(m.group(2))


    def parse_sv_vcf(self) -> Dict[str, List[Tuple]]:
        """
        Parse SV VCF and extract coordinates by SV type.

        Returns
        -------
        dict:
            {
            "DEL": [(chr, start, end)],
            "DUP": [(chr, start, end)],
            "INV": [(chr, start, end)],
            "INS": [(chr, start, end)],
            "BND": [(chr1, pos1, chr2, pos2)]
            }
        """
        vcf = self.open_vcf_autotabix()

        sv_records = {
            "DEL": [],
            "DUP": [],
            "INV": [],
            "INS": [],
            "BND": [],
        }

        for rec in vcf:
            chrom1 = rec.chrom
            pos1 = rec.pos
            info = rec.info
            svtype = info.get("SVTYPE")

            if svtype in ("DEL", "DUP", "INV"):
                end = info.get("END")
                if end:
                    sv_records[svtype].append(
                        (chrom1, pos1, end)
                    )
            elif svtype == "INS":
                svlen = info.get("SVLEN")

                if svlen is None:
                    end = pos1 + 1
                else:
                    if isinstance(svlen, (tuple, list)):
                        svlen = svlen[0]
                    end = pos1 + abs(int(svlen))

                sv_records["INS"].append(
                    (chrom1, pos1, end)
                )

            elif svtype == "BND":
                alt = rec.alts[0]
                chrom2, pos2 = self.parse_bnd_alt(alt)
                if chrom2:
                    sv_records["BND"].append(
                        (chrom1, pos1, chrom2, pos2)
                    )

        return sv_records


    # -----------------------------
    # Circos 输入文件生成
    # -----------------------------
    def write_interval_file(
        sv_records: Dict[str, List[Tuple]],
        out_file: str,
        svtypes=("DEL", "DUP", "INV", "INS"),
    ):
        """
        Write interval-type SVs for Circos plots.
        """
        with open(out_file, "w") as f:
            for svtype in svtypes:
                for chrom, start, end in sv_records.get(svtype, []):
                    f.write(f"{chrom} {start} {end} {svtype}\n")


    def write_link_file(
        sv_records: Dict[str, List[Tuple]],
        out_file: str,
    ):
        """
        Write BND SVs as Circos link format.
        """
        with open(out_file, "w") as f:
            for chrom1, pos1, chrom2, pos2 in sv_records.get("BND", []):
                f.write(
                    f"{chrom1} {pos1} {pos1 + 1} "
                    f"{chrom2} {pos2} {pos2 + 1}\n"
                )

    def run_sv_to_circos(
        self
    ):
        """
        Full pipeline:
        - parse FASTA
        - parse SV VCF
        - generate Circos input files
        """
        outdir = self.outdir
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # FASTA
        chrom_sizes = self.parse_fasta_chrom_sizes()
        SVCircosPrepare.write_circos_karyotype(
            chrom_sizes,
            outdir / "karyotype.txt",
        )

        # VCF
        sv_records = self.parse_sv_vcf()

        self.write_interval_file(
            sv_records,
            outdir / "sv_intervals.txt",
        )

        SVCircosPrepare.write_link_file(
            sv_records,
            outdir / "sv_links.txt",
        )


# -----------------------------
# CLI 入口（可选）
# -----------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Parse SV VCF and FASTA for Circos input"
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--outdir", required=True)

    args = parser.parse_args()
    svCircorsPrepare = SVCircosPrepare(
        vcf_path=args.vcf,
        fasta_path=args.fasta,
        outdir=args.outdir,
    )
    svCircorsPrepare.run_sv_to_circos()

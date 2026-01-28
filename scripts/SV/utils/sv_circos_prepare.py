import re
from pathlib import Path
from typing import Dict, List, Tuple
import pysam
from pyfaidx import Fasta
import subprocess
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


class SVCircosPrepare():
    """
    SV Circos preparation utility class.
    """
    main_chrom = ["chr" + str(c) for c in range(1, 20)] + ['chrX', 'chrY', 'chrM']
    logger.info(f"Main chromosomes set: {main_chrom}")
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
            logger.info(f"Detected compressed VCF: {vcf_path}")
            tbi = vcf_path.with_suffix(".gz.tbi")
            csi = vcf_path.with_suffix(".gz.csi")
            logger.info(f"Checking for index files: {tbi} / {csi}")

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
        chrom_sizes = {chrom: len(fa[chrom]) for chrom in fa.keys() if chrom in self.main_chrom}
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

    def parse_sv_vcf_for_circle(self) -> Dict[str, List[Tuple]]:
        """
        Parse SV VCF and organize coordinates for Circle plot tracks.
        - Blocks: DEL (start, end)
        - Points: INS, DUP, INV (individual breakpoints)
        - Links: TRA (inter-chromosomal connections)
        """
        vcf = self.open_vcf_autotabix()
        
        sv_records = {
            "DEL": [],  # (chrom, start, end)
            "INS": [],  # (chrom, start, end)
            "DUP": [],  # (chrom, start, end)
            "INV": [],  # (chrom, start, end)
            "TRA": []    # (c1, p1, c2, p2)
        }

        for rec in vcf:
            chrom1 = rec.chrom
            pos1 = rec.pos
            info = rec.info
            svtype = info.get("SVTYPE")
            end = getattr(rec, 'stop', None)
            if end is None:
                end_val = info.get("END") # (end,)
                if end_val is not None:
                    end = end_val[0] if isinstance(end_val, (list, tuple)) else end_val

            if chrom1 not in self.main_chrom:
                continue

            if svtype in ("INS","DUP", "INV","DEL") and end is not None:
                if svtype == "INS":
                    sv_records[svtype].append((chrom1,pos1, pos1 + 1)) # INS 位置记录为单点，end = pos + 1,vcf中记录为插入位置+插入长度
                else:
                    sv_records[svtype].append((chrom1,pos1, int(end)))
                continue
            if svtype == "BND":
                alt = str(rec.alts[0])
                chrom2, pos2, _ = self.parse_complex_bnd(alt)
                
                if not chrom2 or chrom2 not in self.main_chrom:
                    continue

                if chrom1 != chrom2:
                    # 识别为 TRA Link (跨圆心连线)
                    link = tuple(sorted([(chrom1, pos1), (chrom2, pos2)]))
                    sv_records["TRA"].append((link[0][0], link[0][1], link[1][0], link[1][1]))
                else:
                    # 同染色体 BND 转为 DEL Block (外圈色块)
                    if self.is_deletion_link(alt):
                        start, end = min(pos1, pos2), max(pos1, pos2)
                        sv_records["DEL"].append((chrom1, start, end))

        # 去重处理
        for key in sv_records:
                    sv_records[key] = list(set(sv_records[key]))
        
        return sv_records

    def parse_complex_bnd(self, alt: str) -> Tuple:
        """Extract chrom2 and pos2 from BND ALT field."""
        match = re.search(r'([\[\]])(.+?):(\d+)([\[\]])', alt)
        if match:
            b1, chrom2, pos2, b2 = match.groups()
            return chrom2, int(pos2), b1
        return None, None, None

    def is_deletion_link(self, alt: str) -> bool:
        """Identify if BND ALT brackets indicate a deletion connection."""
        return ("[" in alt and not alt.startswith("[")) or ("]" in alt and alt.startswith("]"))


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
        sv_records = self.parse_sv_vcf_for_circle()

        # Write Circos input files
        with open(outdir / "blocks_del.bed", "w") as f_blocks:
            for chrom, start, end in sv_records["DEL"]:
                f_blocks.write(f"{chrom}\t{start}\t{end}\n")
        with open(outdir / "points_ins.bed", "w") as f_points_ins:
            for chrom, pos, end in sv_records["INS"]:
                f_points_ins.write(f"{chrom}\t{pos}\t{end}\n")
        with open(outdir / "points_dup.bed", "w") as f_points_dup:
            for chrom, pos, end in sv_records["DUP"]:
                f_points_dup.write(f"{chrom}\t{pos}\t{end}\n")
        with open(outdir / "points_inv.bed", "w") as f_points_inv:
            for chrom, pos, end in sv_records["INV"]:
                f_points_inv.write(f"{chrom}\t{pos}\t{end}\n")
        with open(outdir / "links_tra.bed", "w") as f_links:
            for c1, p1, c2, p2 in sv_records["TRA"]:
                f_links.write(f"{c1}\t{p1}\t{c2}\t{p2}\n")
        


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

from typing import Dict, List, Optional
from cyvcf2 import VCF


class SVVEPAnalyzer:
    """
    High-performance SV-focused VEP CSQ parser.
    Designed for large-scale structural variant analysis.
    """

    TARGET_FIELDS = {
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "VARIANT_CLASS",
        "CANONICAL"
    }

    DISRUPTIVE_TERMS = {
        "transcript_ablation",
        "exon_loss_variant",
        "feature_truncation",
        "feature_elongation"
    }

    def __init__(self, vcf_path: str):
        """
        Initialize the SV VEP analyzer for a given VCF file.

        Parses the CSQ header format and pre-computes field indices
        for high-performance lookups.

        Parameters
        ----------
        vcf_path : str
            Path to the VEP-annotated VCF file.
        """
        self.vcf_path = vcf_path
        self.vcf = VCF(vcf_path)
        self.csq_fields = self._extract_csq_fields()

        # 预计算字段索引（提高性能）
        self.field_index = {
            field: self.csq_fields.index(field)
            for field in self.TARGET_FIELDS
            if field in self.csq_fields
        }

    def _extract_csq_fields(self) -> List[str]:
        """
        Extract CSQ format field names from the VCF header.

        Scans the VCF header for the ``##INFO=<ID=CSQ`` line and
        parses the pipe-delimited field names.

        Returns
        -------
        List[str]
            Ordered list of CSQ field names.

        Raises
        ------
        ValueError
            If the CSQ header line is not found in the VCF.
        """
        for header_line in self.vcf.raw_header.split("\n"):
            if header_line.startswith("##INFO=<ID=CSQ"):
                prefix = "Format: "
                start = header_line.find(prefix)
                format_str = header_line[start + len(prefix):].rstrip('">')
                return format_str.split("|")

        raise ValueError("CSQ header not found in VCF.")

    def _classify_variant(self, consequences: List[str]) -> str:
        """
        Classify a structural variant into a functional impact category.

        Parameters
        ----------
        consequences : List[str]
            List of VEP consequence terms (e.g., "transcript_ablation").

        Returns
        -------
        str
            One of "gene_disrupting", "intergenic", or "gene_overlapping".
        """
        if any(term in self.DISRUPTIVE_TERMS for term in consequences):
            return "gene_disrupting"

        if "intergenic_variant" in consequences:
            return "intergenic"

        return "gene_overlapping"

    def analyze(self) -> Dict[str, int]:
        """
        Analyze all variants in the VCF and return summary statistics.

        Iterates over every variant, classifies it via the CSQ annotation,
        and tallies counts per category.

        Returns
        -------
        Dict[str, int]
            Summary counts with keys: total_sv, gene_disrupting,
            gene_overlapping, intergenic.
        """
        stats = {
            "total_sv": 0,
            "gene_disrupting": 0,
            "gene_overlapping": 0,
            "intergenic": 0,
        }

        for variant in self.vcf:
            stats["total_sv"] += 1

            csq = variant.INFO.get("CSQ")
            if not csq:
                continue

            # 只取第一个 transcript（通常你用了 --pick）
            first_entry = csq.split(",")[0]
            values = first_entry.split("|")

            consequences = []
            if "Consequence" in self.field_index:
                idx = self.field_index["Consequence"]
                if idx < len(values):
                    consequences = values[idx].split("&")

            category = self._classify_variant(consequences)
            stats[category] += 1

        return stats

if __name__ == "__main__":
    analyzer = SVVEPAnalyzer("/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB_only_annotated.vcf")
    summary = analyzer.analyze()

    print(summary)
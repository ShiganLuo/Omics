from typing import Dict,Union, Literal
try:
    from LogUtil import setup_logger
except ImportError:
    from .LogUtil import setup_logger
logger = setup_logger(__name__)
# ── Species alias → canonical genome version ──────────────────────────────
# Maps any common species alias (case-insensitive) to the unique genome
# version identifier used in config["genome"] keys.
# Example: "mouse" / "Mus musculus" / "mm10" → "GRCm39"
SPECIES_TO_GENOME: Dict[str, str] = {
    # mouse
    "mouse": "GRCm39",
    "mus musculus": "GRCm39",
    "mm": "GRCm39",
    "mm39": "GRCm39",
    "grcm39": "GRCm39",
    # human
    "human": "GRCh38",
    "homo sapiens": "GRCh38",
    "hg38": "GRCh38",
    "grch38": "GRCh38",
    # rhesus macaque
    "rhesus": "Mmul_10",
    "rhesus macaque": "Mmul_10",
    "macaca mulatta": "Mmul_10",
    "macaque": "Mmul_10",
    "mmul_10": "Mmul_10",
    "rhemac10": "Mmul_10",
    "mulatta": "Mmul_10"
}

GENOME_TO_SPECIES: Dict[str,str] = {
    "GRCm39": "mouse",
    "T2T_mhaESC_v1.5": "mouse",
    "GRCh38": "human",
    "Mmul_T2T": "macaque",
    "Mmul_10": "macaque"
}

def resolve_genome(organism: Union[str, int, float, None], direction: Literal["version", "alias"] = "version") -> str:
    """Resolve a species alias to its canonical genome version.

    If direction is "version", returns the canonical genome version.
    If direction is "alias", returns the first alias that maps to the given genome version.

    Raises ValueError if the alias is not recognised.
    """

    if isinstance(organism, (int, float)) or not organism:
        logger.info(f"Invalid organism: {organism}, return None")
        return "UNKNOWN"
    if direction == "version":
        key = organism.strip().lower()
        if key in SPECIES_TO_GENOME.keys():
            return SPECIES_TO_GENOME[key]
        if key in SPECIES_TO_GENOME.values():
            return key
    elif direction == "alias":
        key = organism.strip()
        if key in GENOME_TO_SPECIES.keys():
            return GENOME_TO_SPECIES[key]
        if key in GENOME_TO_SPECIES.values():
            return key
    raise ValueError(
        f"Unknown organism alias '{organism}'. "
        f"Known aliases: {sorted(set(SPECIES_TO_GENOME.values()))}"
    )

if __name__ == "__main__":
    pass
from typing import Dict,Union
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
}

def resolve_genome(organism: Union[str, float,None]) -> str:
    """Resolve a species alias to its canonical genome version.

    Raises ValueError if the alias is not recognised.
    """
    if isinstance(organism, float) or not organism:
        logger.info(f"Invalid organism: {organism}, return None")
        return None
    key = organism.strip().lower()
    if key in SPECIES_TO_GENOME.keys():
        return SPECIES_TO_GENOME[key]
    if key in SPECIES_TO_GENOME.values():
        return key
    raise ValueError(
        f"Unknown organism alias '{organism}'. "
        f"Known aliases: {sorted(set(SPECIES_TO_GENOME.values()))}"
    )

if __name__ == "__main__":
    pass
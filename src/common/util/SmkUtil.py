from typing import Dict,Union, Literal
try:
    from LogUtil import setup_logger
except ImportError:
    from .LogUtil import setup_logger
logger = setup_logger(__name__)
SPECIES_TO_GENOME: Dict[str, str] = {
    "mouse": "GRCm39",
    "mus musculus": "GRCm39",
    "mm": "GRCm39",
    "mm39": "GRCm39",
    "grcm39": "GRCm39",
    "human": "GRCh38",
    "homo sapiens": "GRCh38",
    "hg38": "GRCh38",
    "grch38": "GRCh38",
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

WORKFLOW_SPECIES_ALIAS_TO_GENOME: Dict[str, Dict[str, str]] = {
    "pacvar": {
        "mouse": "T2T_mhaESC_v1.5",
        "mus musculus": "T2T_mhaESC_v1.5",
        "mm": "T2T_mhaESC_v1.5",
        "macaque": "Mmul_T2T",
        "rhesus": "Mmul_T2T",
        "rhesus macaque": "Mmul_T2T",
        "macaca mulatta": "Mmul_T2T",
        "mulatta": "Mmul_T2T",
    },
}
WORKFLOW_SPECIES_ALIAS_TO_GENOME["fiberseq"] = WORKFLOW_SPECIES_ALIAS_TO_GENOME["pacvar"]

def resolve_genome(
        organism: Union[str, int, float, None],
        direction: Literal["version", "alias"] = "version",
        workflow: Union[str, None] = None,
    ) -> str:
    """Resolve a species alias to its canonical genome version.

    If direction is "version", returns the canonical genome version.
    If direction is "alias", returns the first alias that maps to the given genome version.

    ``workflow`` (e.g. "PacVar", "Fiberseq") reroutes SPECIES-NAME aliases to
    that workflow's preferred assembly (WORKFLOW_SPECIES_ALIAS_TO_GENOME):
    "mouse" -> T2T_mhaESC_v1.5 and "macaque" -> Mmul_T2T for the long-read
    assembly workflows instead of the standard GRCm39 / Mmul_10. Version
    aliases (mm39, grcm39, ...) and explicit genome keys are never rerouted —
    they always mean exactly what they say.

    Raises ValueError if the alias is not recognised.
    """

    if isinstance(organism, (int, float)) or not organism:
        logger.info(f"Invalid organism: {organism}, return None")
        return "UNKNOWN"
    if direction == "version":
        key = organism.strip().lower()
        if workflow:
            wf_map = WORKFLOW_SPECIES_ALIAS_TO_GENOME.get(str(workflow).strip().lower())
            if wf_map and key in wf_map:
                return wf_map[key]
        if key in SPECIES_TO_GENOME.keys():
            return SPECIES_TO_GENOME[key]
        if key in SPECIES_TO_GENOME.values():
            return key
        for genome in GENOME_TO_SPECIES:
            if genome.lower() == key:
                return genome
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
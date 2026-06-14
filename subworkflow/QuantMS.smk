shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os

indir = config.get("indir", "data/mzml")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
outfiles = config.get("outfiles", [])
samples = config.get("samples", [])
mzml_files = config.get("mzml_files", [])

rule all:
    input:
        outfiles

# Module config dicts
decoy_database_config = {
    "indir": indir,
    "outdir": f"{outdir}/decoy_database",
    "logdir": logdir,
    "Procedure": {
        "openms": config.get("Procedure", {}).get("openms")
    },
    "Params": {
        "decoy_database": config.get("Params", {}).get("decoy_database", {})
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}

search_engine_config = {
    "indir": indir,
    "outdir": f"{outdir}/search_engine",
    "logdir": logdir,
    "samples": samples,
    "mzml_files": mzml_files,
    "Procedure": {
        "comet": config.get("Procedure", {}).get("comet"),
        "msgf": config.get("Procedure", {}).get("msgf"),
        "sage": config.get("Procedure", {}).get("sage")
    },
    "Params": {
        "search_engine": config.get("Params", {}).get("search_engine", {})
    },
    "genome": {
        "decoy_fasta": f"{outdir}/decoy_database/{os.path.basename(config.get('genome', {}).get('fasta', 'protein.fasta'))}_decoy.fasta"
    }
}

psm_rescoring_config = {
    "indir": f"{outdir}/search_engine",
    "outdir": f"{outdir}/psm_rescoring",
    "logdir": logdir,
    "samples": samples,
    "Procedure": {
        "percolator": config.get("Procedure", {}).get("percolator")
    },
    "Params": {
        "psm_rescoring": config.get("Params", {}).get("psm_rescoring", {})
    }
}

psm_fdr_config = {
    "indir": f"{outdir}/psm_rescoring",
    "outdir": f"{outdir}/psm_fdr",
    "logdir": logdir,
    "samples": samples,
    "Procedure": {
        "openms": config.get("Procedure", {}).get("openms")
    },
    "Params": {
        "psm_fdr_control": config.get("Params", {}).get("psm_fdr_control", {})
    }
}

protein_inference_config = {
    "indir": f"{outdir}/psm_fdr",
    "outdir": f"{outdir}/protein_inference",
    "logdir": logdir,
    "samples": samples,
    "Procedure": {
        "epifany": config.get("Procedure", {}).get("epifany")
    },
    "Params": {
        "protein_inference": config.get("Params", {}).get("protein_inference", {})
    }
}

quantification_config = {
    "indir": f"{outdir}/protein_inference",
    "outdir": f"{outdir}/quantification",
    "logdir": logdir,
    "samples": samples,
    "quantification_method": config.get("quantification_method", "lfq"),
    "Procedure": {
        "proteomicslfq": config.get("Procedure", {}).get("proteomicslfq"),
        "proteinquantifier": config.get("Procedure", {}).get("proteinquantifier")
    },
    "Params": {
        "protein_quant": config.get("Params", {}).get("protein_quant", {}),
        "tmt": config.get("Params", {}).get("tmt", {}),
        "lfq": config.get("Params", {}).get("lfq", {}),
        "dia": config.get("Params", {}).get("dia", {})
    }
}

msstats_config = {
    "indir": f"{outdir}/quantification",
    "outdir": f"{outdir}/msstats",
    "logdir": logdir,
    "samples": samples,
    "quantification_method": config.get("quantification_method", "lfq"),
    "Procedure": {
        "msstats": config.get("Procedure", {}).get("msstats")
    },
    "Params": {
        "msstats": config.get("Params", {}).get("msstats", {}),
        "skip_post_msstats": config.get("Params", {}).get("skip_post_msstats", False)
    }
}

# Import modules
module decoy_database:
    snakefile: "../modules/openms/decoydatabase/decoydatabase.smk"
    config: decoy_database_config

module search_engine:
    snakefile: "../modules/openms/searchengine/searchengine.smk"
    config: search_engine_config

module psm_rescoring:
    snakefile: "../modules/openms/psmrescoring/psmrescoring.smk"
    config: psm_rescoring_config

module psm_fdr:
    snakefile: "../modules/openms/psmfdr/psmfdr.smk"
    config: psm_fdr_config

module protein_inference:
    snakefile: "../modules/openms/proteininference/proteininference.smk"
    config: protein_inference_config

module quantification:
    snakefile: "../modules/openms/quantification/quantification.smk"
    config: quantification_config

module msstats:
    snakefile: "../modules/openms/msstats/msstats.smk"
    config: msstats_config

# Use rules from modules
use rule decoy_database_result from decoy_database as QuantMS_decoy_database

use rule search_engine_result from search_engine as QuantMS_search_engine

use rule psm_rescoring_result from psm_rescoring as QuantMS_psm_rescoring

use rule psm_fdr_result from psm_fdr as QuantMS_psm_fdr

use rule protein_inference_result from protein_inference as QuantMS_protein_inference

use rule quantification_result from quantification as QuantMS_quantification

use rule msstats_result from msstats as QuantMS_msstats

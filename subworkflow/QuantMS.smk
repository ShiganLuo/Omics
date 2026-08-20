shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os

indir = config.get("indir")
ROOT_DIR = config.get("ROOT_DIR")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
outfiles = config.get("outfiles", [])
samples = config.get("samples", [])

rule all:
    input:
        outfiles


raw2mzml_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": indir,
    "outdir": f"{outdir}/common/2_raw2mzml",
    "logdir": logdir,
    "samples": samples,
    "raw_files": config.get("raw_files", []),
    "Procedure": {
        "msconvert": config.get("Procedure", {}).get("msconvert")
    },
    "Params": {
        "raw_to_mzml": config.get("Params", {}).get("raw_to_mzml")
    }
}

# Import modules
module raw2mzml:
    snakefile: "../modules/openms/raw2mzml/raw2mzml.smk"
    config: raw2mzml_config
use rule * from raw2mzml as QuantMS_*
logger.info(f"raw2mzml_config: {raw2mzml_config}")

decoy_database_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
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
        "fasta": config.get("genome", {}).get("fasta"),
        "decoy_fasta": config.get("genome", {}).get("decoy_fasta")
    }
}
module decoy_database:
    snakefile: "../modules/openms/decoydatabase/decoydatabase.smk"
    config: decoy_database_config
use rule decoy_database from decoy_database as QuantMS_decoy_database
logger.info(f"decoy_database_config: {decoy_database_config}")

search_engine_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "mzML_dir": raw2mzml_config["outdir"],
    "decoy_dir": decoy_database_config["outdir"],
    "outdir": f"{outdir}/common/3_search_engine",
    "logdir": logdir,
    "samples": samples,
    "Procedure": {
        "comet": config.get("Procedure", {}).get("comet"),
        "msgf": config.get("Procedure", {}).get("msgf"),
        "sage": config.get("Procedure", {}).get("sage")
    },
    "Params": {
        "search_engine": config.get("Params", {}).get("search_engine", {})
    },
    "genome": {
        "decoy_fasta": config.get("genome",{}).get("decoy_fasta")
    }
}

module search_engine:
    snakefile: "../modules/openms/searchengine/searchengine.smk"
    config: search_engine_config
use rule * from search_engine as QuantMS_*
logger.info(f"search_engine_config: {search_engine_config}")

psm_rescoring_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": search_engine_config["outdir"],
    "outdir": f"{outdir}/common/4_psm_rescoring",
    "logdir": logdir,
    "samples": samples,
    "Procedure": {
        "percolator": config.get("Procedure", {}).get("percolator")
    },
    "Params": {
        "psm_rescoring": config.get("Params", {}).get("psm_rescoring", {}),
        "search_engine": config.get("Params", {}).get("search_engine", {}).get("engine")
    }
}

module psm_rescoring:
    snakefile: "../modules/openms/psmrescoring/psmrescoring.smk"
    config: psm_rescoring_config
use rule * from psm_rescoring as QuantMS_*
logger.info(f"psm_rescoring_config: {psm_rescoring_config}")

psm_fdr_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": psm_rescoring_config["outdir"],
    "outdir": f"{outdir}/common/5_psm_fdr",
    "logdir": logdir,
    "samples": samples,
    "Procedure": {
        "openms": config.get("Procedure", {}).get("openms")
    },
    "Params": {
        "psm_fdr_control": config.get("Params", {}).get("psm_fdr_control", {})
    }
}

module psm_fdr:
    snakefile: "../modules/openms/psmfdr/psmfdr.smk"
    config: psm_fdr_config
use rule * from psm_fdr as QuantMS_*
logger.info(f"psm_fdr_config: {psm_fdr_config}")

protein_inference_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": psm_fdr_config["outdir"],
    "outdir": f"{outdir}/common/6_protein_inference",
    "logdir": logdir,
    "samples": samples,
    "Procedure": {
        "epifany": config.get("Procedure", {}).get("epifany")
    },
    "Params": {
        "protein_inference": config.get("Params", {}).get("protein_inference", {})
    }
}

module protein_inference:
    snakefile: "../modules/openms/proteininference/proteininference.smk"
    config: protein_inference_config
use rule * from protein_inference as QuantMS_*
logger.info(f"protein_inference_config: {protein_inference_config}")

quantification_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": protein_inference_config["outdir"],
    "outdir": f"{outdir}/common/7_quantification",
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

module quantification:
    snakefile: "../modules/openms/quantification/quantification.smk"
    config: quantification_config
use rule * from quantification as QuantMS_*
logger.info(f"quantification_config: {quantification_config}")

msstats_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": quantification_config["outdir"],
    "outdir": f"{outdir}/common/8_msstats",
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

module msstats:
    snakefile: "../modules/openms/msstats/msstats.smk"
    config: msstats_config
use rule * from msstats as QuantMS_*
logger.info(f"msstats_config: {msstats_config}")


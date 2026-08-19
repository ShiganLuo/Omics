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





msstats_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
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
    snakefile: "../modules/raw2mzml/raw2mzml.smk"
    config: raw2mzml_config


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
    snakefile: "../modules/decoydatabase/decoydatabase.smk"
    config: decoy_database_config

search_engine_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": decoy_database_config["outdir"],
    "outdir": f"{outdir}/search_engine",
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
        "decoy_fasta": f"{outdir}/decoy_database/{os.path.basename(config.get('genome', {}).get('fasta', 'protein.fasta'))}_decoy.fasta"
    }
}

module search_engine:
    snakefile: "../modules/searchengine/searchengine.smk"
    config: search_engine_config

psm_rescoring_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": search_engine_config["outdir"],
    "outdir": f"{outdir}/psm_rescoring",
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
    snakefile: "../modules/psmrescoring/psmrescoring.smk"
    config: psm_rescoring_config

psm_fdr_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": psm_rescoring_config["outdir"],
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

module psm_fdr:
    snakefile: "../modules/psmfdr/psmfdr.smk"
    config: psm_fdr_config

protein_inference_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": psm_fdr_config["outdir"],
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

module protein_inference:
    snakefile: "../modules/proteininference/proteininference.smk"
    config: protein_inference_config

quantification_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
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

module quantification:
    snakefile: "../modules/quantification/quantification.smk"
    config: quantification_config

module msstats:
    snakefile: "../modules/msstats/msstats.smk"
    config: msstats_config

# Use rules from modules
use rule raw2mzml_result from raw2mzml as QuantMS_raw2mzml

use rule decoy_database_result from decoy_database as QuantMS_decoy_database

use rule search_engine_result from search_engine as QuantMS_search_engine

use rule psm_rescoring_result from psm_rescoring as QuantMS_psm_rescoring

use rule psm_fdr_result from psm_fdr as QuantMS_psm_fdr

use rule protein_inference_result from protein_inference as QuantMS_protein_inference

use rule quantification_result from quantification as QuantMS_quantification

use rule msstats_result from msstats as QuantMS_msstats

shell.prefix("set -x; set -e;")
from snakemake.logging import logger
indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "logs")
outfiles = config.get("outfiles", [])
ROOT_DIR = config.get("ROOT_DIR","")
paired_samples = config.get("sample_params", {}).get("paired_samples", [])
single_samples = config.get("sample_params", {}).get("single_samples", [])
ip_samples = config.get("sample_params", {}).get("ip_samples", [])
input_samples = config.get("sample_params", {}).get("input_samples", [])
treated_ip_samples = config.get("sample_params", {}).get("treated_ip_samples", [])
treated_input_samples = config.get("sample_params", {}).get("treated_input_samples", [])

genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(genome, {}) if genome else {}

rule all:
    input:
        outfiles
trim_galore_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir":  f"{outdir}/cutadapt",
        "logdir": logdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        }
    }
module trim_galore:
    snakefile: "../modules/trim-galore/trim-galore.smk"
    config: trim_galore_config
logger.info(f"trim_galore_config: {trim_galore_config}")
use rule trimming_Paired from trim_galore as MERIP_trimming_Paired
use rule trimming_Single from trim_galore as MERIP_trimming_Single

hisat2_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": trim_galore_config["outdir"],
        "outdir":  f"{outdir}/hisat2",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "hisat2": config.get('Procedure',{}).get('hisat2')
        },
        "Params": {
            "hisat2": {
                "keep_intermediates": config.get('Params',{}).get('hisat2',{}).get('keep_intermediates', True)
            }
        },
        "genome": {
            "fasta": genome_ref.get("fasta"),
            "index_prefix": genome_ref.get("hisat2_index_prefix"),
        }
    }
module hisat2:
    snakefile: "../modules/hisat2/hisat2.smk"
    config: hisat2_config
logger.info(f"hisat2_config: {hisat2_config}")

use rule hisat2_align from hisat2 as MERIP_hisat2_align
use rule hisat2_index from hisat2 as MERIP_hisat2_index

igv_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": hisat2_config["outdir"],
        "outdir":  f"{outdir}/igv",
        "logdir": logdir,
        "Procedure": {
            "samtools": config.get('Procedure',{}).get('samtools'),
            "bamCoverage": config.get('Procedure',{}).get('bamCoverage')
        },
        "Params": {
            "bamCoverage": config.get('Params',{}).get('bamCoverage', {}),
        }
    }
module igv:
    snakefile: "../modules/igv/igv.smk"
    config: igv_config
logger.info(f"igv_config: {igv_config}")
use rule samtools_dedup from igv as MERIP_dedup
use rule wig from igv as MERIP_wig

sexcheck_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/sexcheck",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "sample_bams": {
        s: f"{igv_config['outdir']}/{s}/{s}.dedup.bam"
        for s in paired_samples + single_samples
    },
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "sex_chroms": (config.get("genome", {}) or {}).get("references", {})
            .get(config.get("genome", {}).get("default"), {})
            .get("sex_chroms"),
        "gtf": (config.get("genome", {}) or {}).get("references", {})
            .get(config.get("genome", {}).get("default"), {})
            .get("gtf")
    },
    "Params": {
        "sexcheck": config.get("Params", {}).get("sexcheck", {})
    },
}
module sexcheck:
    snakefile: "../modules/sexcheck/sexcheck.smk"
    config: sexcheck_config
logger.info(f"sexcheck_config: {sexcheck_config}")
use rule sexcheck_idxstats from sexcheck as MERIP_sexcheck_idxstats
use rule sexcheck from sexcheck as MERIP_sexcheck

qc_summary_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "qc_dirs": {
        "sexcheck": sexcheck_config["outdir"],
    },
    "qc_anchors": {
        "sexcheck": "sexcheck.tsv",
    },
}
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
logger.info(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as MERIP_qc_summary

regionplot_params = config.get("Params", {}).get("regionplot", {})
regionplot_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/regionplot",
    "logdir": f"{logdir}/group",
    "samples": single_samples + paired_samples,
    "tracks": [{"label": "bigwig", "dir": igv_config["outdir"], "suffix": ".bigwig"}],
    "Params": {
        "regionplot": regionplot_params
    },
    "genome": {
        "gtf": genome_ref.get("gtf"),
    },
    "resources": config.get("resources", {}),
}
module regionplot:
    snakefile: "../modules/regionplot/regionplot.smk"
    config: regionplot_config
logger.debug(f"regionplot_config: {regionplot_config}")
use rule regionplot from regionplot as MERIP_regionplot
use rule regionplot_result from regionplot as MERIP_regionplot_result

exomePeak_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": igv_config["outdir"],
        "outdir": f"{outdir}/exomePeak",
        "logdir": logdir,
        "gtf": genome_ref.get("gtf"),
        "ip_samples": ip_samples,
        "input_samples": input_samples,
        "treated_ip_samples": treated_ip_samples,
        "treated_input_samples": treated_input_samples,
    }
module exomePeak:
    snakefile: "../modules/exomePeak/exomePeak.smk"
    config: exomePeak_config
logger.info(f"exomePeak_config: {exomePeak_config}")
use rule diff_exomePeak from exomePeak as MERIP_diff_exomePeak
use rule call_exomePeak from exomePeak as MERIP_call_exomePeak

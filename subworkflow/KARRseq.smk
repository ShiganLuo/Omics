from snakemake.logging import logger
import os
shell.prefix("set -x; set -e;")
config.setdefault("OUTDIR_ROOT", config.get("Params", {}).get("workflow", {}).get("outdir", "output"))
indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "log")
outfiles = config.get("outfiles", [])
paired_samples = config.get("sample_params", {}).get("paired_samples", [])
single_samples = config.get("sample_params", {}).get("single_samples", [])
samples = paired_samples + single_samples
ROOT_DIR = config.get("ROOT_DIR", os.path.dirname(os.path.dirname(__file__)))

genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(genome, {}) if genome else {}

rule all:
    input:
        outfiles


if paired_samples:
    seqprep_config = {
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir": f"{outdir}/merge",
        "logdir": logdir,
        "samples": paired_samples,
        "resources": config.get("resources", {}),
        "Procedure": {
            "SeqPrep": config.get("Procedure", {}).get("SeqPrep") or "SeqPrep",
            "seqkit": config.get("Procedure", {}).get("seqkit") or "seqkit"
        }
    }
    module seqprep:
        snakefile: "../modules/seqprep/seqprep.smk"
        config: seqprep_config
    logger.info(f"seqprep_config: {seqprep_config}")
    use rule seqprep_merge from seqprep as KARRseq_seqprep_merge

star_index_dir = genome_ref.get("star_index_dir")

star_fastqs = {
    s: f"{outdir}/merge/{s}/{s}_merge.fastq.gz" for s in paired_samples
}
star_fastqs.update({s: f"{indir}/{s}.fastq.gz" for s in single_samples})

star_chimeric_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "outdir": f"{outdir}/bam",
    "logdir": logdir,
    "env": config.get("env", {}),
    "sample_fastqs": star_fastqs,
    "genome": {
        "star_index_dir": star_index_dir,
    },
    "Procedure": {
        "STAR": config.get("Procedure", {}).get("STAR") or "STAR",
        "samtools": config.get("Procedure", {}).get("samtools") or "samtools",
    },
    "Params": {
        "star": config.get("Params", {}).get("star", {}),
    },
    "resources": config.get("resources", {}),
}
module star_chimeric:
    snakefile: "../modules/star/star_chimeric/star_chimeric.smk"
    config: star_chimeric_config
logger.info(f"star_chimeric_config: {star_chimeric_config}")
use rule KARRseq_star_align from star_chimeric as KARRseq_star_align

sexcheck_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/sexcheck",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "sample_bams": {
        s: f"{outdir}/bam/{s}/{s}_Aligned.sortedByCoord.out.bam"
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
    "resources": config.get("resources", {}),
}
module sexcheck:
    snakefile: "../modules/sexcheck/sexcheck.smk"
    config: sexcheck_config
logger.info(f"sexcheck_config: {sexcheck_config}")
use rule sexcheck_idxstats from sexcheck as KARRseq_sexcheck_idxstats
use rule sexcheck from sexcheck as KARRseq_sexcheck

qc_summary_config = {
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
    "resources": config.get("resources", {}),
}
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
logger.info(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as KARRseq_qc_summary



karrseq_config = {
    "indir": f"{outdir}/bam",
    "outdir": f"{outdir}/chimeric",
    "logdir": logdir,
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools") or "samtools"
    },
    "Params": {
        "karrseq": config.get("Params", {}).get("karrseq", {})
    },
    "resources": config.get("resources", {}),
}
module karrseq:
    snakefile: "../modules/karrseq/karrseq.smk"
    config: karrseq_config
logger.info(f"karrseq_config: {karrseq_config}")
use rule karrseq_chimeric_to_pairs from karrseq as KARRseq_chimeric_to_pairs
use rule karrseq_remove_duplicates from karrseq as KARRseq_remove_duplicates
use rule karrseq_ligation from karrseq as KARRseq_ligation

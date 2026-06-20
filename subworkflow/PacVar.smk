shell.prefix("set -x; set -e;")
from snakemake.logging import logger
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
outfiles = config.get("outfiles", [])
samples = config.get("samples", [])
skip_snp = config.get("Params", {}).get("skip_snp", False)
skip_sv = config.get("Params", {}).get("skip_sv", False)
skip_phase = config.get("Params", {}).get("skip_phase", False)
skip_repeat = config.get("Params", {}).get("skip_repeat", False)
skip_telomere = config.get("Params", {}).get("skip_telomere", False)
snv_caller = config.get("Params", {}).get("snv_caller", "deepvariant")

rule all:
    input:
        outfiles


pbmm2_config = {
    "indir": indir,
    "outdir": f"{outdir}/bam/1_sorted_bam",
    "logdir": logdir,
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "Procedure": {
        "pbmm2": config.get("Procedure", {}).get("pbmm2")
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module pbmm2:
    snakefile: "../modules/pbmm2/pbmm2.smk"
    config: pbmm2_config
logger.info(f"pbmm2_config: {pbmm2_config}")
use rule pbmm2_align from pbmm2 as PacVar_pbmm2_align


gatk_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "indir": pbmm2_config["outdir"],
    "outdir": f"{outdir}/bam/2_markdup_bam",
    "logdir": logdir,
    "input_bam_substring": "sorted",
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools"),
    },
    "Params": {
        "gatk": {
            "addReadsGroup": {
                "RGLB": config.get("addReadsGroup", {}).get("RGLB"),
                "RGPL": config.get("addReadsGroup", {}).get("RGPL"),
                "RGPU": config.get("addReadsGroup", {}).get("RGPU")
            }
        }
    },
    "genome": {
        "fasta": config.get("genome",{}).get("fasta"),
        "fai_index": config.get("genome",{}).get("fai_index"),
        "dict_index": config.get("genome",{}).get("dict_index")
    }
}

module gatk_prepare:
    snakefile: "../modules/gatk/gatk_prepare.smk"
    config: gatk_prepare_config
logger.info(f"gatk_prepare parameters: {gatk_prepare_config}")
use rule gatk_index from gatk_prepare as PacVar_gatk_index
use rule addReadsGroup from gatk_prepare as PacVar_addReadsGroup
use rule MarkDuplicates from gatk_prepare as PacVar_MarkDuplicates

if not skip_snp:
    if snv_caller == "deepvariant":
        deepvariant_config = {
            "ROOT_DIR": ROOT_DIR,
            "indir": gatk_prepare_config["outdir"],
            "outdir": f"{outdir}/variation/germline_snv_indel",
            "logdir": logdir,
            "samples": samples,
            "bam_substring": "sorted_markdup",
            "Procedure": {
                "deepvariant": config.get("Procedure", {}).get("deepvariant")
            },
            "Params": {
                "deepvariant": config.get("Params", {}).get("deepvariant", {})
            },
            "genome": {
                "fasta": config.get("genome", {}).get("fasta"),
                "fai": config.get("genome", {}).get("fai")
            }
        }
        module deepvariant:
            snakefile: "../modules/deepvariant/deepvariant.smk"
            config: deepvariant_config
        logger.info(f"deepvariant_config: {deepvariant_config}")
        use rule deepvariant_run from deepvariant as PacVar_deepvariant_run
    elif snv_caller == "gatk4":
        gatk_germline_config = {
            "ROOT_DIR": ROOT_DIR,
            "indir": gatk_prepare_config["outdir"],
            "outdir": f"{outdir}/variation/germline_snv_indel",
            "logdir": logdir,
            "Procedure": {
                "gatk": config.get("Procedure", {}).get("gatk"),
                "samtools": config.get("Procedure", {}).get("samtools")
            },
            "genome": {
                "fasta": config.get("genome", {}).get("fasta"),
                "fai_index": config.get("genome", {}).get("fai"),
                "dict_index": config.get("genome", {}).get("dict")
            }
        }
        module gatk_germline:
            snakefile: "../modules/gatk/gatk_germline/gatk_germline.smk"
            config: gatk_germline_config
        logger.info(f"gatk_germline_config: {gatk_germline_config}")
        use rule HaplotypeCaller from gatk_germline as PacVar_HaplotypeCaller
        use rule filterHaplotypeCallerVcf from gatk_germline as PacVar_filterHaplotypeCallerVcf
    else:
        raise ValueError(f"Unsupported snv_caller: {snv_caller}")

# ============================================================
# Step 4: SV variant calling (pbsv)
# ============================================================
if not skip_sv:
    pbsv_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": gatk_prepare_config["outdir"],
        "outdir": f"{outdir}/variation/germline_sv",
        "logdir": logdir,
        "samples": samples,
        "bam_substring": "sorted_markdup",
        "Procedure": {
            "pbsv": config.get("Procedure", {}).get("pbsv")
        },
        "genome": {
            "fasta": config.get("genome", {}).get("fasta")
        }
    }
    module pbsv:
        snakefile: "../modules/pbsv/pbsv.smk"
        config: pbsv_config
    logger.info(f"pbsv_config: {pbsv_config}")
    use rule pbsv_discover from pbsv as PacVar_pbsv_discover
    use rule pbsv_call from pbsv as PacVar_pbsv_call
    
# ============================================================
# Step 5: Phasing (optional, requires both SNP + SV done)
# ============================================================
if not skip_phase and not skip_snp and not skip_sv:
    hiphase_snp_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": gatk_prepare_config["outdir"],
        "outdir": f"{outdir}/variation/germline_snv_indel",
        "logdir": logdir,
        "samples": samples,
        "bam_dir": gatk_prepare_config["outdir"],
        "input_bam_substring": "sorted_markdup",
        "input_vcf_substring": "filtered",
        "output_substring": "",
        "vcf_dir": f"{outdir}/variation/germline_snv_indel",
        "Procedure": {
            "hiphase": config.get("Procedure", {}).get("hiphase")
        },
        "genome": {
            "fasta": config.get("genome", {}).get("fasta")
        }
    }
    module hiphase_snp:
        snakefile: "../modules/hiphase/hiphase.smk"
        config: hiphase_snp_config
    logger.info(f"hiphase_snp_config: {hiphase_snp_config}")
    use rule hiphase_phase from hiphase_snp as PacVar_hiphase_snp

    hiphase_sv_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": gatk_prepare_config["outdir"],
        "outdir": f"{outdir}/variation/germline_sv",
        "logdir": logdir,
        "samples": samples,
        "bam_dir": gatk_prepare_config["outdir"],
        "vcf_dir": f"{outdir}/variation/germline_sv",
        "input_bam_substring": "sorted_markdup",
        "input_vcf_substring": "sv",
        "output_substring": "sv",
        "Procedure": {
            "hiphase": config.get("Procedure", {}).get("hiphase")
        },
        "genome": {
            "fasta": config.get("genome", {}).get("fasta")
        }
    }
    module hiphase_sv:
        snakefile: "../modules/hiphase/hiphase.smk"
        config: hiphase_sv_config
    logger.info(f"hiphase_sv_config: {hiphase_sv_config}")
    use rule hiphase_phase from hiphase_sv as PacVar_hiphase_sv

# ============================================================
# Step 6: Repeat characterization (optional)
# ============================================================
if not skip_repeat and config.get("genome", {}).get("repeat_bed") is not None:
    trgt_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": gatk_prepare_config["outdir"],
        "outdir": f"{outdir}/repeat/trgt",
        "logdir": logdir,
        "samples": samples,
        "bam_substring": "sorted_markdup",
        "Procedure": {
            "trgt": config.get("Procedure", {}).get("trgt")
        },
        "Params": {
            "trgt": config.get("Params", {}).get("trgt", {})
        },
        "genome": {
            "fasta": config.get("genome", {}).get("fasta"),
            "fai": config.get("genome", {}).get("fai"),
            "repeat_bed": config.get("reference", {}).get("repeat_bed")
        }
    }
    module trgt:
        snakefile: "../modules/trgt/trgt.smk"
        config: trgt_config
    logger.info(f"trgt_config: {trgt_config}")
    use rule trgt_genotype from trgt as PacVar_trgt_genotype
    use rule trgt_plot from trgt as PacVar_trgt_plot

# ============================================================
# Step 7: Telomere analysis (optional)
# ============================================================
if not skip_telomere:
    telomere_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": gatk_prepare_config["outdir"],
        "outdir": f"{outdir}/repeat/telomere",
        "logdir": logdir,
        "samples": samples,
        "bam_substring": "sorted_markdup",
        "Params": {
            "telogator2": config.get("Params", {}).get("telogator2", {})
        },
    }
    module telomere:
        snakefile: "../modules/telomere/telomere.smk"
        config: telomere_config
    logger.info(f"telomere_config: {telomere_config}")
    use rule telogator2_run from telomere as PacVar_telogator2_run

    centromere_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": indir,
        "outdir": f"{outdir}/repeat/centromere",
        "logdir": logdir,
        "samples": samples,
        "ROOT_DIR": config.get("ROOT_DIR", "."),
        "Params": {
            "hifiasm": config.get("Params", {}).get("hifiasm", {}),
            "RepeatMasker": config.get("Params", {}).get("RepeatMasker", {}),
        },
    }
    module centromere:
        snakefile: "../modules/centromere/centromere.smk"
        config: centromere_config
    logger.info(f"centromere_config: {centromere_config}")
    use rule * from centromere as PacVar_centromere_*

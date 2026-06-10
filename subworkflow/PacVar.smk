shell.prefix("set -x; set -e;")
from snakemake.logging import logger

indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
outfiles = config.get("outfiles", [])
samples = config.get("samples", [])
skip_snp = config.get("Params", {}).get("skip_snp", False)
skip_sv = config.get("Params", {}).get("skip_sv", False)
skip_phase = config.get("Params", {}).get("skip_phase", False)
skip_repeat = config.get("Params", {}).get("skip_repeat", False)
snv_caller = config.get("Params", {}).get("snv_caller", "deepvariant")

rule all:
    input:
        outfiles

# ============================================================
# Step 1: PBMM2 alignment
# ============================================================
pbmm2_config = {
    "indir": indir,
    "outdir": f"{outdir}/bam/1_raw_bam",
    "logdir": logdir,
    "samples": samples,
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

# ============================================================
# Step 2: SAMTOOLS sort + index
# ============================================================
samtools_sort_config = {
    "indir": pbmm2_config["outdir"],
    "outdir": f"{outdir}/bam/2_sorted_bam",
    "logdir": logdir,
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module samtools_sort:
    snakefile: "../modules/samtools/sort/samtools_sort.smk"
    config: samtools_sort_config
logger.info(f"samtools_sort_config: {samtools_sort_config}")
use rule bam_sort from samtools_sort as PacVar_bam_sort

# ============================================================
# Step 3: SNP variant calling
# ============================================================
if not skip_snp:
    if snv_caller == "deepvariant":
        deepvariant_config = {
            "indir": samtools_sort_config["outdir"],
            "outdir": f"{outdir}/snp/deepvariant",
            "logdir": logdir,
            "samples": samples,
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
            "indir": samtools_sort_config["outdir"],
            "outdir": f"{outdir}/snp/gatk4",
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
    else:
        raise ValueError(f"Unsupported snv_caller: {snv_caller}")

# ============================================================
# Step 4: SV variant calling (pbsv)
# ============================================================
if not skip_sv:
    pbsv_config = {
        "indir": samtools_sort_config["outdir"],
        "outdir": f"{outdir}/sv/pbsv",
        "logdir": logdir,
        "samples": samples,
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

    # bgzip + index for SV VCF
    tabix_sv_config = {
        "indir": f"{outdir}/sv/pbsv/call",
        "outdir": f"{outdir}/sv/pbsv/bgzip",
        "logdir": logdir,
        "Procedure": {
            "bgzip": config.get("Procedure", {}).get("bgzip"),
            "tabix": config.get("Procedure", {}).get("tabix")
        }
    }
    module tabix_sv:
        snakefile: "../modules/tabix/tabix.smk"
        config: tabix_sv_config
    logger.info(f"tabix_sv_config: {tabix_sv_config}")
    use rule tabix_bgzip from tabix_sv as PacVar_tabix_bgzip
    use rule tabix_index from tabix_sv as PacVar_tabix_index

# ============================================================
# Step 5: Phasing (optional, requires both SNP + SV done)
# ============================================================
if not skip_phase and not skip_snp and not skip_sv:
    hiphase_snp_config = {
        "indir": samtools_sort_config["outdir"],
        "outdir": f"{outdir}/phasing/snp",
        "logdir": logdir,
        "samples": samples,
        "bam_dir": samtools_sort_config["outdir"],
        "vcf_dir": f"{outdir}/snp/deepvariant" if snv_caller == "deepvariant" else f"{outdir}/snp/gatk4",
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
        "indir": samtools_sort_config["outdir"],
        "outdir": f"{outdir}/phasing/sv",
        "logdir": logdir,
        "samples": samples,
        "bam_dir": samtools_sort_config["outdir"],
        "vcf_dir": f"{outdir}/sv/pbsv/bgzip",
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
        "indir": samtools_sort_config["outdir"],
        "outdir": f"{outdir}/repeat/trgt",
        "logdir": logdir,
        "samples": samples,
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

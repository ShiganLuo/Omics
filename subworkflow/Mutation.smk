shell.prefix("set -x; set -e;")
from snakemake.logging import logger
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
outfiles = config.get("outfiles", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
rule all:
    input:
        outfiles

fastqc_raw_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": indir,
        "outdir":  f"{outdir}/quality/fastqc/raw",
        "logdir": logdir,
        "log_suffix": "raw.txt",
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "fastqc": config.get("Procedure", {}).get("fastqc") or "fastqc"
        }
    }
module fastqc_raw:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_raw_config
logger.info(f"fastqc_raw_config: {fastqc_raw_config}")
use rule fastqc from fastqc_raw as Mutation_fastqc_raw

cutadapt_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": indir,
        "outdir": f"{outdir}/fastq/cutadapt",
        "logdir": logdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        },
        "Params": {
            "trim_galore": {
                "quality": config.get('Params',{}).get("trim_galore", {}).get('quality')
            }
        },
    }

module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
logger.info(f"Cutadapt parameters: {cutadapt_config}")
use rule trimming_Paired from cutadapt as Mutation_trimming_Paired

fastqc_trimmed_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": cutadapt_config["outdir"],
        "outdir":  f"{outdir}/quality/fastqc/trimmed",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "log_suffix": "trimmed.txt",
        "Procedure": {
            "fastqc": config.get("Procedure", {}).get("fastqc")
        }
    }
module fastqc_trimmed:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_trimmed_config
logger.info(f"fastqc_trimmed_config: {fastqc_trimmed_config}")
use rule fastqc from fastqc_trimmed as Mutation_fastqc_trimmed

bwa_mem2_confg = {
    "indir": cutadapt_config["outdir"],
    "outdir":  f"{outdir}/bam/1_sorted_bam",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "bwaMem2": config.get("Procedure",{}).get("bwaMem2"),
        "samtools": config.get("Procedure",{}).get("samtools")
    },
    "genome": {
        "fasta": config.get("genome",{}).get("fasta"),
        "index_prefix": config.get("genome",{}).get("bwaMem2_index_prefix")
    }
}

module bwa_mem2:
    snakefile: "../modules/bwa-mem2/bwa-mem2.smk"
    config: bwa_mem2_confg
logger.info(f"BWA-MEM2 parameters: {bwa_mem2_confg}")
use rule bwaMem2_index from bwa_mem2 as Mutation_bwaMem2_index
use rule bwaMem2_alignment from bwa_mem2 as Mutation_bwaMem2_alignment

samtools_config = {
    "ROOT_DIR": ROOT_DIR,
    "indir": bwa_mem2_confg["outdir"],
    "outdir": bwa_mem2_confg["outdir"],
    "logdir": logdir,
    "Procedure": {
        "samtools": config.get("Procedure",{}).get("samtools")
    },
    "Params": {
        "samtools": {
            "onlykeep_properpair": config.get("Params",{}).get("samtools",{}).get("onlykeep_properpair")
        }
    }
}

module samtools:
    snakefile: "../modules/samtools/samtools.smk"
    config: samtools_config
logger.info(f"samtools parameters: {samtools_config}")
use rule bam_flagstat from samtools as Mutation_bam_flagstat


gatk_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "indir": bwa_mem2_confg["outdir"],
    "outdir": f"{outdir}/bam/2_markdup_bam",
    "logdir": logdir,
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
use rule gatk_index from gatk_prepare as Mutation_gatk_index
use rule addReadsGroup from gatk_prepare as Mutation_addReadsGroup
use rule MarkDuplicates from gatk_prepare as Mutation_MarkDuplicates

gatk_somatic_config = {
    "ROOT_DIR": ROOT_DIR,
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/somatic_snv_indel",
    "logdir": logdir,
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
    },
    "genome": {
        "fasta": config.get("genome",{}).get("fasta")
    }
}
module gatk_somatic:
    snakefile: "../modules/gatk/gatk_somatic/gatk_somatic.smk"
    config: gatk_somatic_config
logger.info(f"gatk_somatic parameters: {gatk_somatic_config}")
use rule somaticMutect2 from gatk_somatic as Mutation_somaticMutect2

gatk_germline_config = {
    "ROOT_DIR": ROOT_DIR,
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/germline_snv_indel",
    "logdir": logdir,
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
    },
    "genome": {
        "fasta": config.get("genome",{}).get("fasta"),
        "fai_index": config.get("genome",{}).get("fai_index"),
        "dict_index": config.get("genome",{}).get("dict_index")
    }
}
module gatk_germline:
    snakefile: "../modules/gatk/gatk_germline/gatk_germline.smk"
    config: gatk_germline_config
logger.info(f"gatk_germline parameters: {gatk_germline_config}")
use rule * from  gatk_germline as Mutation_*

if config.get("Params", {}).get("somatic_spectrum", {}).get("sample_somatic_vcf_dict", {}) and  config.get("Params", {}).get("somatic_spectrum", {}).get("sample_group_dict", {}):
    logger.info(f"sample_somatic_vcf_dict and sample_group_dict provided for somatic_spectrum. Will run somatic_spectrum module.")
    somatic_spectrum_config = {
        "indir": gatk_somatic_config["outdir"],
        "outdir": f"{outdir}/results/spectrum",
        "logdir": logdir,
        "ROOT_DIR":ROOT_DIR,
        "sample_somatic_vcf_dict": config.get("Params", {}).get("somatic_spectrum", {}).get("sample_somatic_vcf_dict", {}),
        "sample_group_dict": config.get("Params", {}).get("somatic_spectrum", {}).get("sample_group_dict", {}),
        "genome": {
            "fasta": config.get("genome",{}).get("fasta")
        }
    }
    module somatic_spectrum:
        snakefile: "../modules/spectrum/spectrum.smk"
        config: somatic_spectrum_config
    logger.info(f"somatic_spectrum parameters: {somatic_spectrum_config}")
    use rule somatic_spectrum from somatic_spectrum as Mutation_somatic_spectrum
else:
    logger.info(f"sample_somatic_vcf_dict or sample_group_dict not provided for somatic_spectrum. Skipping somatic_spectrum module.")

# ============================================
# Fragment Size Analysis (cfDNA)
# ============================================
skip_fragment_size = config.get("Params", {}).get("skip_fragment_size", False)
if not skip_fragment_size:
    fragment_size_config = {
        "indir": gatk_prepare_config["outdir"],
        "outdir": f"{outdir}/results/fragment_size",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
        "samples": paired_samples + single_samples,
        "Procedure": {
            "samtools": config.get("Procedure", {}).get("samtools")
        }
    }
    module fragment_size:
        snakefile: "../modules/fragment_size/fragment_size.smk"
        config: fragment_size_config
    logger.info(f"fragment_size_config: {fragment_size_config}")
    use rule samtools_stats from fragment_size as Mutation_samtools_stats
    use rule getFragmentSize from fragment_size as Mutation_getFragmentSize
    use rule plotFragmentSize from fragment_size as Mutation_plotFragmentSize
else:
    logger.info("Skipping fragment_size module (skip_fragment_size=True)")

# ============================================
# SV Detection with Manta (cfDNA)
# ============================================
skip_sv = config.get("Params", {}).get("skip_sv", False)
if not skip_sv:
    manta_config = {
        "indir": gatk_prepare_config["outdir"],
        "outdir": f"{outdir}/variation/somatic_sv",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
        "samples": paired_samples + single_samples,
        "Procedure": {
            "manta": config.get("Procedure", {}).get("manta"),
            "samtools": config.get("Procedure", {}).get("samtools")
        },
        "Params": {
            "manta": config.get("Params", {}).get("manta", {})
        },
        "genome": {
            "fasta": config.get("genome", {}).get("fasta")
        }
    }
    module manta:
        snakefile: "../modules/manta/manta.smk"
        config: manta_config
    logger.info(f"manta_config: {manta_config}")
    use rule manta_config from manta as Mutation_manta_config
    use rule manta_run from manta as Mutation_manta_run
else:
    logger.info("Skipping Manta SV detection (skip_sv=True)")

# ============================================
# CNV Detection with CNVkit (cfDNA)
# ============================================
skip_cnv = config.get("Params", {}).get("skip_cnv", False)
if not skip_cnv:
    control_samples = config.get("control_samples", [])
    cnvkit_config = {
        "indir": gatk_prepare_config["outdir"],
        "outdir": f"{outdir}/variation/germline_cnv",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
        "samples": paired_samples + single_samples,
        "control_samples": control_samples,
        "Procedure": {
            "cnvkit": config.get("Procedure", {}).get("cnvkit"),
            "samtools": config.get("Procedure", {}).get("samtools")
        },
        "Params": {
            "cnvkit": config.get("Params", {}).get("cnvkit", {})
        },
        "genome": {
            "fasta": config.get("genome", {}).get("fasta"),
            "access": config.get("genome", {}).get("access")
        }
    }
    module cnvkit:
        snakefile: "../modules/cnvkit/cnvkit.smk"
        config: cnvkit_config
    logger.info(f"cnvkit_config: {cnvkit_config}")
    if control_samples:
        use rule cnvkit_reference from cnvkit as Mutation_cnvkit_reference
    use rule cnvkit_batch from cnvkit as Mutation_cnvkit_batch
else:
    logger.info("Skipping CNVkit CNV detection (skip_cnv=True)")


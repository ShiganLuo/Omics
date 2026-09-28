shell.prefix("set -x; set -e;")
from snakemake.logging import logger
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("Params", {}).get("workflow").get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow").get("outdir", "output")
logdir = config.get("Params", {}).get("workflow").get("logdir", "logs")
samples = config.get("Params", {}).get("workflow").get("samples", [])
snv_caller = config.get("Params", {}).get("workflow").get("snv_caller", "deepvariant")
outfiles = config.get("outfiles", [])
# Genome reference resolution — nested `default` + `references` pattern
default_genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(default_genome, {}) if default_genome else {}

rule all:
    input:
        outfiles


genome_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "genome": {
        "fasta": genome_ref.get("fasta"),
    }
}
module genome:
    snakefile: "../modules/genome/genome.smk"
    config: genome_config
use rule genome_index from genome as PacVar_genome_index

pbmm2_config = {
    "indir": indir,
    "outdir": f"{outdir}/common/2_sorted_bam",
    "logdir": logdir,
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "Procedure": {
        "pbmm2": config.get("Procedure", {}).get("pbmm2")
    },
    "genome": {
        "fasta": genome_ref.get("fasta")
    }
}
module pbmm2:
    snakefile: "../modules/pbmm2/pbmm2.smk"
    config: pbmm2_config
logger.debug(f"pbmm2_config: {pbmm2_config}")
use rule pbmm2_align from pbmm2 as PacVar_pbmm2_align


gatk_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": pbmm2_config["outdir"],
    "outdir": f"{outdir}/common/3_markdup_bam",
    "logdir": logdir,
    "input_bam_substring": "sorted",
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools"),
    },
    "Params": {
        "gatk": {
            "addReadsGroup": {
                "RGLB": config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGLB"),
                "RGPL": config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGPL"),
                "RGPU": config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGPU")
            },
            "javaOptions": config.get("Params", {}).get("gatk", {}).get("javaOptions"),
            "tmp-dir": config.get("Params", {}).get("gatk", {}).get("tmp-dir")
        }
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "fai_index": genome_ref.get("fai_index"),
        "dict_index": genome_ref.get("dict_index")
    }
}

module gatk_prepare:
    snakefile: "../modules/gatk/gatk_prepare.smk"
    config: gatk_prepare_config
logger.debug(f"gatk_prepare parameters: {gatk_prepare_config}")
use rule gatk_index from gatk_prepare as PacVar_gatk_index
use rule addReadsGroup from gatk_prepare as PacVar_addReadsGroup
use rule MarkDuplicates from gatk_prepare as PacVar_MarkDuplicates


if snv_caller == "deepvariant":
    deepvariant_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
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
            "fasta": genome_ref.get("fasta"),
            "fai": genome_ref.get("fai")
        }
    }
    module deepvariant:
        snakefile: "../modules/deepvariant/deepvariant.smk"
        config: deepvariant_config
    logger.debug(f"deepvariant_config: {deepvariant_config}")
    use rule deepvariant_run from deepvariant as PacVar_deepvariant_run
elif snv_caller == "gatk4":
    # Future: GATK HaplotypeCaller path — currently nothing in node.py outfiles
    # references it, so no use rule here yet.
    pass

pbsv_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/germline_sv",
    "logdir": logdir,
    "samples": samples,
    "bam_substring": "sorted_markdup",
    "Procedure": {
        "pbsv": config.get("Procedure", {}).get("pbsv")
    },
    "genome": {
        "fasta": genome_ref.get("fasta")
    }
}
module pbsv:
    snakefile: "../modules/pbsv/pbsv.smk"
    config: pbsv_config
logger.debug(f"pbsv_config: {pbsv_config}")
use rule pbsv_discover from pbsv as PacVar_pbsv_discover
use rule pbsv_call from pbsv as PacVar_pbsv_call

hiphase_snp_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/germline_snv_indel",
    "logdir": logdir,
    "samples": samples,
    "bam_dir": gatk_prepare_config["outdir"],
    "input_bam_substring": "sorted_markdup",
    "input_vcf_substring": "",
    "output_substring": "",
    "vcf_dir": f"{outdir}/variation/germline_snv_indel",
    "Procedure": {
        "hiphase": config.get("Procedure", {}).get("hiphase")
    },
    "genome": {
        "fasta": genome_ref.get("fasta")
    }
}
module hiphase_snp:
    snakefile: "../modules/hiphase/hiphase.smk"
    config: hiphase_snp_config
logger.debug(f"hiphase_snp_config: {hiphase_snp_config}")
use rule hiphase_phase from hiphase_snp as PacVar_hiphase_snp

hiphase_sv_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
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
        "fasta": genome_ref.get("fasta")
    }
}
module hiphase_sv:
    snakefile: "../modules/hiphase/hiphase.smk"
    config: hiphase_sv_config
logger.debug(f"hiphase_sv_config: {hiphase_sv_config}")
use rule hiphase_phase from hiphase_sv as PacVar_hiphase_sv


trgt_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
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
        "fasta": genome_ref.get("fasta"),
        "fai": genome_ref.get("fai"),
        "repeat_bed": genome_ref.get("repeat_bed")
    }
}
module trgt:
    snakefile: "../modules/trgt/trgt.smk"
    config: trgt_config
logger.debug(f"trgt_config: {trgt_config}")
use rule trgt_genotype from trgt as PacVar_trgt_genotype
use rule trgt_plot from trgt as PacVar_trgt_plot


# ============================================================
# Step 7: Telomere analysis + Centromere
# ============================================================
centromere_outdir = f"{outdir}/repeat/centromere"

telomere_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/repeat/telomere",
    "logdir": logdir,
    "samples": samples,
    "bam_substring": "sorted_markdup",
    "assembly_dir": centromere_outdir,
    "Procedure": {
        "telogator2": config.get("Procedure", {}).get("telogator2"),
        "tidk": config.get("Procedure", {}).get("tidk")
    },
    "Params": {
        "telogator2": config.get("Params", {}).get("telogator2", {})
    },
}
module telomere:
    snakefile: "../modules/telomere/telomere.smk"
    config: telomere_config
logger.debug(f"telomere_config: {telomere_config}")
use rule telogator2_run from telomere as PacVar_telogator2_run
use rule assembly_telomere_scan from telomere as PacVar_assembly_telomere_scan
use rule read_density_telomere from telomere as PacVar_read_density_telomere
use rule tidk_scan from telomere as PacVar_tidk_scan

centromere_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": indir,
    "outdir": centromere_outdir,
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
logger.debug(f"centromere_config: {centromere_config}")
use rule * from centromere as PacVar_centromere_*


# ============================================================
# Step 8: eccDNA / ecDNA reconstruction (Decoil)
# ============================================================
# Reference:
#   Giurgiu et al., "Reconstructing extrachromosomal DNA structural heterogeneity
#   from long-read sequencing data using Decoil", Genome Research, 2024.
#   DOI: 10.1101/gr.279123.124.  https://github.com/madagiurgiu25/decoil-pre
decoil_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/eccdna",
    "logdir": logdir,
    "samples": samples,
    "bam_substring": "sorted_markdup",
    "Procedure": {
        "decoil": config.get("Procedure", {}).get("decoil")
    },
    "Params": {
        "decoil": config.get("Params", {}).get("decoil", {})
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "gtf": genome_ref.get("gtf")
    }
}
module decoil:
    snakefile: "../modules/decoil/decoil.smk"
    config: decoil_config
logger.debug(f"decoil_config: {decoil_config}")
use rule decoil_reconstruct from decoil as PacVar_decoil_reconstruct
use rule decoil_result from decoil as PacVar_decoil_result

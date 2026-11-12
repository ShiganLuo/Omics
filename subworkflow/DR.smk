shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "logs")
genome_paired_samples = config.get("Params", {}).get("workflow", {}).get("genome_paired_samples", {})
genome_single_samples = config.get("Params", {}).get("workflow", {}).get("genome_single_samples", {})
dr_pairs = config.get("Params", {}).get("workflow", {}).get("dr_pairs", {})          # {wes: rna}
drclass_samples = config.get("Params", {}).get("workflow", {}).get("drclass_samples", {})
wes_paired_samples = config.get("Params", {}).get("workflow", {}).get("wes_paired_samples", [])
wes_single_samples = config.get("Params", {}).get("workflow", {}).get("wes_single_samples", [])
wes_indir = config.get("Params", {}).get("workflow", {}).get("wes_indir", "")         # WES 原始 fastq 根
trimmer = config.get('Params', {}).get('workflow', {}).get('trimmer') or "cutadapt"
outfiles = config.get("outfiles", [])
default_genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(default_genome, {}) if default_genome else {}
genomes = set(genome_paired_samples.keys()).union(set(genome_single_samples.keys()))
logger.info(f"DR workflow: outdir={outdir}, rna_genomes={genomes}, wes_n={len(wes_paired_samples)+len(wes_single_samples)}, dr_pairs={dr_pairs}")

rule all:
    input:
        outfiles

# ############################################
# 链 0: WES 体细胞链 (复刻 Mutation.smk 必需段: fastqc->cutadapt->bwa-mem2
#       ->samtools->gatk_prepare(markdup)->targets->gatk_somatic->gatk_germline)
# 输出全部落在 {outdir}/WES/ 下, 与链 1 的 {outdir}/RNA/ 互不干扰.
# 输入 = 框架统一的 raw 整理目录 {outdir}/common/1_raw_data/ (WES 与 RNA
# 样本都在那里, 命名 <ID>_{1,2}.fq.gz).
# ############################################

wes_raw_dir = os.path.join(outdir, "common", "1_raw_data")

fastqc_raw_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": wes_raw_dir,
    "outdir": f"{outdir}/WES/quality/fastqc/raw",
    "logdir": f"{logdir}/sample",
    "log_suffix": "raw.txt",
    "paired_samples": wes_paired_samples,
    "single_samples": wes_single_samples,
    "Procedure": {"fastqc": config.get("Procedure", {}).get("fastqc") or "fastqc"},
}
module fastqc_raw:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_raw_config
use rule fastqc from fastqc_raw as DR_WES_fastqc_raw

wes_cutadapt_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": wes_raw_dir,
    "outdir": f"{outdir}/WES/common/2_trimmed_fastq",
    "logdir": f"{logdir}/sample",
    "Procedure": {"trim_galore": config.get('Procedure', {}).get('trim_galore')},
    "Params": {
        "trim_galore": {"quality": config.get('Params', {}).get("trim_galore", {}).get('quality')}
    },
}
module wes_cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: wes_cutadapt_config
use rule trimming_Paired from wes_cutadapt as DR_WES_trimming_Paired
use rule trimming_Single from wes_cutadapt as DR_WES_trimming_Single

fastqc_trimmed_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": wes_cutadapt_config["outdir"],
    "outdir": f"{outdir}/WES/quality/fastqc/trimmed",
    "logdir": f"{logdir}/sample",
    "paired_samples": wes_paired_samples,
    "single_samples": wes_single_samples,
    "log_suffix": "trimmed.txt",
    "Procedure": {"fastqc": config.get("Procedure", {}).get("fastqc")},
}
module fastqc_trimmed:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_trimmed_config
use rule fastqc from fastqc_trimmed as DR_WES_fastqc_trimmed

bwa_mem2_confg = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": wes_cutadapt_config["outdir"],
    "outdir": f"{outdir}/WES/common/3_sorted_bam",
    "logdir": f"{logdir}/sample",
    "paired_samples": wes_paired_samples,
    "single_samples": wes_single_samples,
    "Procedure": {
        "bwaMem2": config.get("Procedure", {}).get("bwaMem2"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "index_prefix": genome_ref.get("bwaMem2_index_prefix")
    }
}
module bwa_mem2:
    snakefile: "../modules/bwa-mem2/bwa-mem2.smk"
    config: bwa_mem2_confg
use rule bwaMem2_index from bwa_mem2 as DR_WES_bwaMem2_index
use rule bwaMem2_alignment from bwa_mem2 as DR_WES_bwaMem2_alignment

samtools_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": bwa_mem2_confg["outdir"],
    "outdir": bwa_mem2_confg["outdir"],
    "logdir": f"{logdir}/sample",
    "Procedure": {"samtools": config.get("Procedure", {}).get("samtools")},
    "Params": {
        "samtools": {
            "onlykeep_properpair": config.get("Params", {}).get("samtools", {}).get("onlykeep_properpair")
        }
    }
}
module samtools:
    snakefile: "../modules/samtools/samtools.smk"
    config: samtools_config
use rule bam_flagstat from samtools as DR_WES_bam_flagstat

_targets_bed = genome_ref.get("targets") or f"{outdir}/WES/common/5_targets_bed/targets.bed"
gatk_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": bwa_mem2_confg["outdir"],
    "outdir": f"{outdir}/WES/common/4_markdup_bam",
    "logdir": f"{logdir}/sample",
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
        "fasta": genome_ref.get("fasta"),
        "fai_index": genome_ref.get("fai_index"),
        "dict_index": genome_ref.get("dict_index"),
        "targets": _targets_bed,
    },
    "resources": config.get("resources", {}),
}
module gatk_prepare:
    snakefile: "../modules/gatk/gatk_prepare.smk"
    config: gatk_prepare_config
use rule gatk_index from gatk_prepare as DR_WES_gatk_index
use rule addReadsGroup from gatk_prepare as DR_WES_addReadsGroup
use rule MarkDuplicates from gatk_prepare as DR_WES_MarkDuplicates
use rule collectHsMetrics from gatk_prepare as DR_WES_collectHsMetrics

targets_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/WES/common/5_targets_bed",
    "logdir": f"{logdir}/group",
    "bam_suffix": ".sorted_markdup",
    "samples": wes_paired_samples + wes_single_samples,
    "Params": {"derive": config.get("Params", {}).get("derive_targets", {})},
}
module targets:
    snakefile: "../modules/targets/targets.smk"
    config: targets_config
use rule derive_targets from targets as DR_WES_derive_targets

gatk_somatic_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "tmp_dir": config.get("tmp_dir"),
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/WES/variation/somatic_snv_indel",
    "logdir": f"{logdir}/group",
    "Procedure": {"gatk": config.get("Procedure", {}).get("gatk")},
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "targets": _targets_bed,
        "germline_resource": genome_ref.get("germline_resource"),
        "strain_resource": genome_ref.get("strain_resource")
    },
    "pair_strains": config.get("Params", {}).get("workflow", {}).get("pair_strains", {}),
    "design_pairs": config.get("Params", {}).get("workflow", {}).get("paired_design", []),
    "Params": {
        "gatk": {
            "Mutect2": {
                "panel_of_normals": config.get("Params", {}).get("gatk", {}).get("Mutect2", {}).get("panel_of_normals")
            }
        }
    },
}
module gatk_somatic:
    snakefile: "../modules/gatk/gatk_somatic/gatk_somatic.smk"
    config: gatk_somatic_config
use rule somaticMutect2 from gatk_somatic as DR_WES_somaticMutect2
use rule germline_resource_strain from gatk_somatic as DR_WES_germline_resource_strain
use rule filterMutectCalls from gatk_somatic as DR_WES_filterMutectCalls

gatk_germline_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/WES/variation/germline_snv_indel",
    "logdir": f"{logdir}/sample",
    "Procedure": {"gatk": config.get("Procedure", {}).get("gatk")},
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "fai_index": genome_ref.get("fai_index"),
        "dict_index": genome_ref.get("dict_index"),
        "targets": _targets_bed
    }
}
module gatk_germline:
    snakefile: "../modules/gatk/gatk_germline/gatk_germline.smk"
    config: gatk_germline_config
use rule * from gatk_germline as DR_WES_germline_*

# ############################################
# 链 1: RNA 链 (cutadapt -> star 8_SNP_bam -> gatk_prepare markdup
#       -> gatk_RNAseq calling; 对齐 RNAseq.smk 接线, 无报告类产物)
# 注: run.py 框架已把 meta fastq 软链整理到 {outdir}/common/1_raw_data/
#     且命名为 <ID>_{1,2}.fq.gz, 直接消费 (indir 已被框架改写指向那里).
# star 索引单独一个模块实例, index 写到 {outdir}/RNA/genome/index/{genome},
# align 实例通过 index_dir 指向同一目录 (RNAseq.smk 的双实例模式).
# ############################################

star_config_for_index = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/RNA/genome",
    "logdir_index": os.path.join(logdir, "group"),
    "Procedure": {"star": config.get('Procedure', {}).get('star')},
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module star_for_index:
    snakefile: "../modules/star/polygenomes/star.smk"
    config: star_config_for_index
use rule star_index from star_for_index as DR_RNA_star_index

rna_cutadapt_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": indir,
    "outdir": f"{outdir}/RNA/common/2_trimmed_fastq",
    "logdir": os.path.join(logdir, "sample"),
    "genome_paired_samples": genome_paired_samples,
    "genome_single_samples": genome_single_samples,
    "Procedure": config.get("Procedure", {}),
    "Params": config.get("Params", {}),
    "resources": config.get("resources", {}),
}
module rna_cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: rna_cutadapt_config
use rule trimming_Paired from rna_cutadapt as DR_RNA_trimming_Paired
use rule trimming_Single from rna_cutadapt as DR_RNA_trimming_Single

star_config_for_SNP = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "indir": rna_cutadapt_config["outdir"],
    "outdir": f"{outdir}/RNA/common/8_SNP_bam",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_index": os.path.join(logdir, "group"),
    "isGenomeSubdir": False,
    "index_dir": f"{outdir}/RNA/genome",
    "genome_paired_samples": genome_paired_samples,
    "genome_single_samples": genome_single_samples,
    "env": config.get("env", {}),
    "Procedure": {"star": config.get('Procedure', {}).get('STAR_SNP')},
    "Params": config.get('Params', {}),
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module star_for_SNP:
    snakefile: "../modules/star/polygenomes/star.smk"
    config: star_config_for_SNP
use rule star_align from star_for_SNP as DR_RNA_star_align_for_SNP

gatk_prepare_SNP_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": star_config_for_SNP["outdir"],
    "outdir": f"{outdir}/RNA/common/10_markdup_deContaimation_bam",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module gatk_prepare_SNP:
    snakefile: "../modules/gatk/polygenomes/gatk_prepare.smk"
    config: gatk_prepare_SNP_config
use rule * from gatk_prepare_SNP as DR_RNA_SNP_*

gatk_RNAseq_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_SNP_config['outdir'],
    "outdir": f"{outdir}/RNA/variation/germline_snv_indel_RNAseq",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module gatk_RNAseq:
    snakefile: "../modules/gatk/polygenomes/gatk_RNAseq/gatk_RNAseq.smk"
    config: gatk_RNAseq_config
use rule * from gatk_RNAseq as DR_RNA_SNP_*

# ############################################
# 链 2: DRclass module (M3 pileup -> 状态机 -> O1/O2/O3)
# samples 表由 node.py runDR 显式生成(全部具体文件路径, WES/RNA 产物对齐到
# 链 0/链 1 的输出目录), 无目录扫描.
# ############################################

drclass_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "outdir": f"{outdir}/DRclass",
    "logdir": os.path.join(logdir, "sample"),
    "samples": drclass_samples,
    "genome": config.get('genome', {}),
    "Procedure": {
        "python3": config.get("Procedure", {}).get("python3") or "python3",
        "bcftools": config.get("Procedure", {}).get("bcftools") or "bcftools",
    },
    "Params": config.get("Params", {}),
    "resources": config.get("resources", {}),
}
module DRclass:
    snakefile: "../modules/DRclass/DRclass.smk"
    config: drclass_config
use rule * from DRclass as DR_drclass_*
logger.debug(f"drclass_config: {drclass_config}")

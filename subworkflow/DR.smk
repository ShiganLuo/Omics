shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "logs")
genome_paired_samples = config.get("Params", {}).get("workflow", {}).get("genome_paired_samples", {})
genome_single_samples = config.get("Params", {}).get("workflow", {}).get("genome_single_samples", {})
dr_pairs = config.get("Params", {}).get("workflow", {}).get("dr_pairs", {})          # {dna: rna}
drclass_samples = config.get("Params", {}).get("workflow", {}).get("drclass_samples", {})
dna_paired_samples = config.get("Params", {}).get("workflow", {}).get("dna_paired_samples", [])
dna_single_samples = config.get("Params", {}).get("workflow", {}).get("dna_single_samples", [])
dna_indir = config.get("Params", {}).get("workflow", {}).get("dna_indir", "")         # DNA 原始 fastq 根
trimmer = config.get('Params', {}).get('workflow', {}).get('trimmer') or "trim_galore"
aligner_TEtranscripts = config.get('Params', {}).get('workflow', {}).get('aligner_TEtranscripts') or "star"
outfiles = config.get("outfiles", [])
default_genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(default_genome, {}) if default_genome else {}
genomes = set(genome_paired_samples.keys()).union(set(genome_single_samples.keys()))
# per-genome all-samples map (TEtranscripts/StringTie/arriba/report) and the
# mixed maps for XenofilteR host/contaminant separation (RNAseq.smk 同款)
genome_samples = {}
for _g, _s in genome_paired_samples.items():
    genome_samples[_g] = list(_s)  # copy: never mutate config's own lists
for _g, _s in genome_single_samples.items():
    genome_samples.setdefault(_g, []).extend(_s)
_all_paired = {s for ss in genome_paired_samples.values() for s in ss}
_all_single = {s for ss in genome_single_samples.values() for s in ss}
genome_paired_samples_mixed = {g: _all_paired for g in genome_paired_samples}
genome_single_samples_mixed = {g: _all_single for g in genome_single_samples}
logger.info(f"DR workflow: outdir={outdir}, rna_genomes={genomes}, dna_n={len(dna_paired_samples)+len(dna_single_samples)}, dr_pairs={dr_pairs}, trimmer={trimmer}, aligner_TEtranscripts={aligner_TEtranscripts}")

rule all:
    input:
        outfiles

# ############################################
# 链 0: DNA 体细胞链(data_type=DNA 行) — 复刻 Mutation.smk 全部分析功能,
#       assay 列(WGS/Exome/Panel)驱动 CNVkit --seq-method 与捕获区域 targets,
#       不假定 DNA 测序类型:
#       fastqc->cutadapt->bwa-mem2->samtools->gatk_prepare(markdup)->targets
#       ->gatk_somatic->msi->gatk_germline->somatic_spectrum(+signature)
#       ->fragment_size->manta(germline+somatic SV)->cnvkit->burden
#       ->sexcheck->qc_summary->Mutation_report
# 输出全部落在 {outdir}/DNA/ 下, 与链 1 的 {outdir}/RNA/ 互不干扰.
# 输入 = node.py runDR 写入的 raw 整理根 dna_indir(meta data_type=DNA 行的
# 软链落点 {outdir}/DNA/common/1_raw_data/, 命名 <ID>_{1,2}.fq.gz).
# ############################################

dna_raw_dir = dna_indir

fastqc_raw_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": dna_raw_dir,
    "outdir": f"{outdir}/DNA/QC/fastqc/raw",
    "logdir": f"{logdir}/sample",
    "log_suffix": "raw.txt",
    "paired_samples": dna_paired_samples,
    "single_samples": dna_single_samples,
    "Procedure": {"fastqc": config.get("Procedure", {}).get("fastqc") or "fastqc"},
}
module fastqc_raw:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_raw_config
use rule fastqc from fastqc_raw as DR_DNA_fastqc_raw

# --- trimmer 三选一 (RNA 侧同款): cutadapt / trimmomatic / trim_galore ---
dna_trimmed_dir = f"{outdir}/DNA/common/2_trimmed_fastq"
if trimmer == "cutadapt":
    dna_cutadapt_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": dna_raw_dir,
        "outdir": dna_trimmed_dir,
        "logdir": os.path.join(logdir, "sample"),
        "Procedure": {
            "cutadapt": config.get('Procedure', {}).get('cutadapt')
        },
        "Params": {
            # 参数精确到 rule:Params.cutadapt.<rule>.<key>(模块读点)
            "cutadapt": config.get('Params', {}).get("cutadapt", {}),
        },
    }
    module dna_cutadapt:
        snakefile: "../modules/cutadapt/cutadapt.smk"
        config: dna_cutadapt_config
    logger.debug(f"dna_cutadapt_config: {dna_cutadapt_config}")
    use rule trimming_Paired from dna_cutadapt as DR_DNA_cutadapt_Paired
    use rule trimming_Single from dna_cutadapt as DR_DNA_cutadapt_Single
elif trimmer == "trimmomatic":
    dna_trimmomatic_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": dna_raw_dir,
        "outdir": dna_trimmed_dir,
        "logdir": os.path.join(logdir, "sample"),
        "Procedure": {
            "trimmomatic": config.get('Procedure', {}).get('trimmomatic')
        },
        "Params": {
            # 参数精确到 rule:整段透传(adapter_pe 挂 trimmomatic_Paired,
            # adapter_se 挂 trimmomatic_Single 各自 rule 子键)
            "trimmomatic": config.get('Params', {}).get("trimmomatic", {}),
        },
    }
    module dna_trimmomatic:
        snakefile: "../modules/trimmomatic/trimmomatic.smk"
        config: dna_trimmomatic_config
    logger.debug(f"dna_trimmomatic_config: {dna_trimmomatic_config}")
    use rule trimmomatic_Paired from dna_trimmomatic as DR_DNA_trimmomatic_Paired
    use rule trimmomatic_Single from dna_trimmomatic as DR_DNA_trimmomatic_Single
else:
    logger.info("use default trimmer: trim_galore")
    dna_trim_galore_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": dna_raw_dir,
        "outdir": dna_trimmed_dir,
        "logdir": os.path.join(logdir, "sample"),
        "Procedure": {
            "trim_galore": config.get('Procedure', {}).get('trim_galore')
        },
        "Params": {
            # 参数精确到 rule:Params.trim_galore.<rule>.<key>(与 RNAseq.json 一致)
            "trim_galore": config.get('Params', {}).get("trim_galore", {}),
        },
    }
    module dna_trim_galore:
        snakefile: "../modules/trim-galore/trim-galore.smk"
        config: dna_trim_galore_config
    logger.debug(f"dna_trim_galore_config: {dna_trim_galore_config}")
    use rule trimming_Paired from dna_trim_galore as DR_DNA_trim_galore_Paired
    use rule trimming_Single from dna_trim_galore as DR_DNA_trim_galore_Single

fastqc_trimmed_config = {

    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": dna_trimmed_dir,
    "outdir": f"{outdir}/DNA/quality/fastqc/trimmed",
    "logdir": f"{logdir}/sample",
    "paired_samples": dna_paired_samples,
    "single_samples": dna_single_samples,
    "log_suffix": "trimmed.txt",
    "Procedure": {"fastqc": config.get("Procedure", {}).get("fastqc")},
}
module fastqc_trimmed:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_trimmed_config
use rule fastqc from fastqc_trimmed as DR_DNA_fastqc_trimmed

bwa_mem2_confg = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": dna_trimmed_dir,
    "outdir": f"{outdir}/DNA/common/3_sorted_bam",
    "logdir": f"{logdir}/sample",
    "paired_samples": dna_paired_samples,
    "single_samples": dna_single_samples,
    "Procedure": {
        "bwaMem2": config.get("Procedure", {}).get("bwaMem2"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "index_prefix": genome_ref.get("bwaMem2_index_prefix")
    },
    "Params": {
        # temp 规范(modules.md):用户 Params.bwa-mem2.keep_intermediates 透传到
        # 模块读点 Params.bwa-mem2.keep_intermediates(缺省 True=兜底保留)
        "bwa-mem2": {
            "keep_intermediates": config.get("Params", {}).get("bwa-mem2", {}).get("keep_intermediates", True)
        }
    }
}
module bwa_mem2:
    snakefile: "../modules/bwa-mem2/bwa-mem2.smk"
    config: bwa_mem2_confg
use rule bwaMem2_index from bwa_mem2 as DR_DNA_bwaMem2_index
use rule bwaMem2_alignment from bwa_mem2 as DR_DNA_bwaMem2_alignment

samtools_config = {
    "resources": config.get("resources", {}),
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
use rule bam_flagstat from samtools as DR_DNA_bam_flagstat

_targets_bed = genome_ref.get("targets") or f"{outdir}/DNA/common/5_targets_bed/targets.bed"
gatk_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": bwa_mem2_confg["outdir"],
    "outdir": f"{outdir}/DNA/common/4_markdup_bam",
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
        },
        # temp 规范(modules.md):addReadsGroup 中间 bam/bai 受字段控制
        "gatk_prepare": {
            "keep_intermediates": config.get("Params", {}).get("gatk_prepare", {}).get("keep_intermediates", True)
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
use rule gatk_index from gatk_prepare as DR_DNA_gatk_index
use rule addReadsGroup from gatk_prepare as DR_DNA_addReadsGroup
use rule MarkDuplicates from gatk_prepare as DR_DNA_MarkDuplicates
use rule collectHsMetrics from gatk_prepare as DR_DNA_collectHsMetrics

targets_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/DNA/common/5_targets_bed",
    "logdir": f"{logdir}/group",
    "bam_suffix": ".sorted_markdup",
    "samples": dna_paired_samples + dna_single_samples,
    "Params": {"derive": config.get("Params", {}).get("derive_targets", {})},
}
module targets:
    snakefile: "../modules/targets/targets.smk"
    config: targets_config
use rule derive_targets from targets as DR_DNA_derive_targets

gatk_somatic_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "tmp_dir": config.get("tmp_dir"),
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/DNA/variation/somatic_snv_indel",
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
        },
        # temp 规范(modules.md):somaticMutect2 中间 vcf/stats 受字段控制
        "gatk_somatic": {
            "keep_intermediates": config.get("Params", {}).get("gatk_somatic", {}).get("keep_intermediates", True)
        }
    },
}
module gatk_somatic:
    snakefile: "../modules/gatk/gatk_somatic/gatk_somatic.smk"
    config: gatk_somatic_config
use rule somaticMutect2 from gatk_somatic as DR_DNA_somaticMutect2
use rule germline_resource_strain from gatk_somatic as DR_DNA_germline_resource_strain
use rule filterMutectCalls from gatk_somatic as DR_DNA_filterMutectCalls

# ============================================
# MSI (msisensorpro, Mutation.smk 同款)
# ============================================
msi_pairs = {
    rec["experimental"]: rec["control"]
    for rec in config.get("Params", {}).get("workflow", {}).get("paired_design", [])
    if rec.get("experimental") and rec.get("control")
}
msi_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/DNA/stability/msi",
    "logdir": f"{logdir}/group",
    "bam_substring": "sorted_markdup",
    "pairs": msi_pairs,
    "Procedure": {
        "msisensorpro": config.get("Procedure", {}).get("msisensorpro"),
    },
    "Params": {
        "msi": config.get("Params", {}).get("msi", {}),
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "fai": genome_ref.get("fai_index"),
    },
}
module msi:
    snakefile: "../modules/msi/msi.smk"
    config: msi_config
logger.debug(f"msi_config: {msi_config}")
use rule msi_scan from msi as DR_DNA_msi_scan
use rule msi_call from msi as DR_DNA_msi_call
use rule msi_plot from msi as DR_DNA_msi_plot
use rule msi_result from msi as DR_DNA_msi_result

gatk_germline_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/DNA/variation/germline_snv_indel",
    "logdir": f"{logdir}/sample",
    "Procedure": {"gatk": config.get("Procedure", {}).get("gatk")},
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "fai_index": genome_ref.get("fai_index"),
        "dict_index": genome_ref.get("dict_index"),
        "targets": _targets_bed
    },
    "Params": {
        # temp 规范(modules.md):HaplotypeCaller raw vcf 受字段控制;DR 的 raw.vcf.gz
        # 是 DRclass 交付物,模板不设字段(缺省 True=兜底保留),wrap 本身安全
        "gatk_germline": {
            "keep_intermediates": config.get("Params", {}).get("gatk_germline", {}).get("keep_intermediates", True)
        }
    }
}
module gatk_germline:
    snakefile: "../modules/gatk/gatk_germline/gatk_germline.smk"
    config: gatk_germline_config
use rule * from gatk_germline as DR_DNA_germline_*

# ============================================
# Somatic Mutation Spectrum + Signature (Mutation.smk 同款)
# ============================================
# Wired unconditionally; node.py registers the spectrum output only when
# design pairs exist (outfiles drive the DAG).
somatic_spectrum_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_somatic_config["outdir"],
    "outdir": f"{outdir}/DNA/stability/spectrum",
    "logdir": f"{logdir}/group",
    "sample_somatic_vcf_dict": config.get("Params", {}).get("somatic_spectrum", {}).get("sample_somatic_vcf_dict", {}),
    "sample_group_dict": config.get("Params", {}).get("somatic_spectrum", {}).get("sample_group_dict", {}),
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "signature_matrix": genome_ref.get("signature_matrix"),
    },
    "signature": config.get("Params", {}).get("somatic_spectrum", {}).get("signature", {}),
    # spectrum.smk reads config.plot.formats / config.signature.formats at the
    # module-config top level (see modules/spectrum/spectrum.json); forward the
    # workflow-side Params.somatic_spectrum.plot the same way as signature
    "plot": config.get("Params", {}).get("somatic_spectrum", {}).get("plot", {}),
}
module somatic_spectrum:
    snakefile: "../modules/spectrum/spectrum.smk"
    config: somatic_spectrum_config
logger.info(f"somatic_spectrum parameters: {somatic_spectrum_config}")
use rule somatic_spectrum from somatic_spectrum as DR_DNA_somatic_spectrum
use rule somatic_signature_refit from somatic_spectrum as DR_DNA_somatic_signature_refit
use rule somatic_signature_nmf from somatic_spectrum as DR_DNA_somatic_signature_nmf
use rule spectrum_result from somatic_spectrum as DR_DNA_spectrum_result

# ============================================
# Fragment Size Analysis (cfDNA, Mutation.smk 同款)
# ============================================
# Wired unconditionally; node.py controls whether it runs by registering
# (or omitting) the fragment-size outputs (Params.fragment_size.enabled).
fragment_size_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_suffix": ".sorted_markdup",
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/DNA/results/fragment_size",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": dna_paired_samples + dna_single_samples,
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    }
}
module fragment_size:
    snakefile: "../modules/fragment_size/fragment_size.smk"
    config: fragment_size_config
logger.info(f"fragment_size_config: {fragment_size_config}")
use rule samtools_stats from fragment_size as DR_DNA_samtools_stats
use rule getFragmentSize from fragment_size as DR_DNA_getFragmentSize
use rule plotFragmentSize from fragment_size as DR_DNA_plotFragmentSize

# ============================================
# SV Detection with Manta (cfDNA, Mutation.smk 同款)
# ============================================
# Wired unconditionally; node.py controls whether it runs by registering
# (or omitting) the Manta outputs (Params.manta.enabled).
manta_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_suffix": ".sorted_markdup",
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/DNA/variation/germline_sv",
    "logdir": f"{logdir}/sample",
    "samples": dna_paired_samples + dna_single_samples,
    "Procedure": {
        "manta": config.get("Procedure", {}).get("manta"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "Params": {
        "manta": config.get("Params", {}).get("manta", {})
    },
    "genome": {
        "fasta": genome_ref.get("fasta")
    }
}
module manta:
    snakefile: "../modules/manta/manta.smk"
    config: manta_config
logger.info(f"manta_config: {manta_config}")
use rule manta_config_germline from manta as DR_DNA_manta_config_germline
use rule manta_run_germline from manta as DR_DNA_manta_run_germline

# Somatic (tumor/normal) SV calling via the same manta module, instantiated
# with a separate outdir/logdir. Wired unconditionally; node.py registers the
# somaticSV output only when matched-normal design pairs exist (mirroring the
# somatic SNV / Mutect2 pair judgment), which drives the DAG.
manta_somatic_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_suffix": ".sorted_markdup",
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/DNA/variation/somatic_sv",
    "logdir": f"{logdir}/group",
    "samples": dna_paired_samples + dna_single_samples,
    "Procedure": {
        "manta": config.get("Procedure", {}).get("manta"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "Params": {
        "manta": config.get("Params", {}).get("manta", {})
    },
    "genome": {
        "fasta": genome_ref.get("fasta")
    }
}
module manta_somatic:
    snakefile: "../modules/manta/manta.smk"
    config: manta_somatic_config
logger.info(f"manta_somatic_config: {manta_somatic_config}")
use rule manta_config_somatic from manta_somatic as DR_DNA_manta_config_somatic
use rule manta_run_somatic from manta_somatic as DR_DNA_manta_run_somatic

# ============================================
# CNV Detection with CNVkit (cfDNA, Mutation.smk 同款)
# ============================================
# Wired unconditionally; node.py controls whether it runs by registering
# (or omitting) the CNVkit outputs (Params.cnvkit.enabled; control_samples
# additionally registers the reference panel output).
control_samples = config.get("Params", {}).get("workflow", {}).get("control_samples", [])
cnvkit_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_suffix": ".sorted_markdup",
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/DNA/variation/cnv",
    "logdir": f"{logdir}/group",
    "samples": dna_paired_samples + dna_single_samples,
    "control_samples": control_samples,
    "paired_design": config.get("Params", {}).get("workflow", {}).get("paired_design", []),
    "Procedure": {
        "cnvkit": config.get("Procedure", {}).get("cnvkit"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "Params": {
        "cnvkit": config.get("Params", {}).get("cnvkit", {})
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "access": genome_ref.get("access"),
        "targets": _targets_bed
    }
}
module cnvkit:
    snakefile: "../modules/cnvkit/cnvkit.smk"
    config: cnvkit_config
logger.info(f"cnvkit_config: {cnvkit_config}")
use rule cnvkit_targets from cnvkit as DR_DNA_cnvkit_targets
use rule cnvkit_reference from cnvkit as DR_DNA_cnvkit_reference
use rule cnvkit_batch from cnvkit as DR_DNA_cnvkit_batch
use rule cnvkit_paired from cnvkit as DR_DNA_cnvkit_paired

# ============================================
# Mutation Burden (SNV/InDel TMB + SV + CNV, Mutation.smk 同款)
# ============================================
# Wired unconditionally; node.py controls whether it runs by registering
# (or omitting) the burden outputs (Params.burden.enabled). The three
# sample -> file maps are filled by node.py from the somatic/cnv outputs
# above; hand-filled config entries win.
burden_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/DNA/stability/burden",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "sample_snv_vcf_dict": config.get("Params", {}).get("burden", {}).get("sample_snv_vcf_dict", {}),
    "sample_sv_vcf_dict": config.get("Params", {}).get("burden", {}).get("sample_sv_vcf_dict", {}),
    "sample_cnv_dict": config.get("Params", {}).get("burden", {}).get("sample_cnv_dict", {}),
    "Params": {
        "burden": config.get("Params", {}).get("burden", {})
    },
    "genome": {
        "targets": _targets_bed,
        "gtf": genome_ref.get("gtf")
    }
}
module burden:
    snakefile: "../modules/burden/burden.smk"
    config: burden_config
logger.info(f"burden_config: {burden_config}")
use rule burden_snv_indel from burden as DR_DNA_burden_snv_indel
use rule burden_sv from burden as DR_DNA_burden_sv
use rule burden_cnv from burden as DR_DNA_burden_cnv
use rule burden_aggregate from burden as DR_DNA_burden_aggregate

# ============================================
# Sex check (coverage on sex chromosomes, Mutation.smk 同款)
# ============================================
# Wired unconditionally over the markdup BAMs; node.py registers the
# sexcheck outputs (subworkflow.md QC wiring).
sexcheck_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/DNA/QC/sexcheck",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": dna_paired_samples + dna_single_samples,
    "sample_bams": {
        s: f"{gatk_prepare_config['outdir']}/{s}/{s}.sorted_markdup.bam"
        for s in dna_paired_samples + dna_single_samples
    },
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "sex_chroms": genome_ref.get("sex_chroms")
    },
    "Params": {
        "sexcheck": config.get("Params", {}).get("sexcheck", {})
    },
}
module sexcheck:
    snakefile: "../modules/sexcheck/sexcheck.smk"
    config: sexcheck_config
logger.info(f"sexcheck_config: {sexcheck_config}")
use rule sexcheck_idxstats from sexcheck as DR_DNA_sexcheck_idxstats
use rule sexcheck from sexcheck as DR_DNA_sexcheck

# ============================================
# QC Summary (cross-sample metric aggregation, Mutation.smk 同款)
# ============================================
# Wired unconditionally; node.py registers the qc_summary outputs so the
# aggregation always runs for the delivered output set (subworkflow.md
# 「QC 汇总与 report 模块规范」).
qc_summary_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/DNA/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": dna_paired_samples + dna_single_samples,
    "qc_dirs": {
        "fastqc_raw": fastqc_raw_config["outdir"],
        "fastqc_trimmed": fastqc_trimmed_config["outdir"],
        "flagstat": samtools_config["outdir"],
        "bwa_mem2": f"{logdir}/sample",
        "markdup": gatk_prepare_config["outdir"],
        "hs_metrics": gatk_prepare_config["outdir"],
        "msi": msi_config["outdir"],
        "sexcheck": sexcheck_config["outdir"],
    },
    "qc_anchors": {
        "flagstat": "{sample}/{sample}_flagstat.txt",
        "bwa_mem2": "{sample}/bwa_mem2-alignment.txt",
        "markdup": "{sample}/{sample}.Markdup-metrics.txt",
        "msi": "{sample}/{sample}_msi_pro.tsv",
        "sexcheck": "sexcheck.tsv",
    },
    # MSI anchors exist only for the tumor sample of each design pair
    "qc_anchor_samples": {
        "msi": [rec["experimental"] for rec in config.get("Params", {}).get("workflow", {}).get("paired_design", []) if rec.get("experimental")],
    },
}
# hs_metrics anchor + plot exist only for exome/hybrid runs (targets BED
# resolved by node.py); WGS runs keep the stage unanchored and best-effort.
if _targets_bed:
    qc_summary_config["qc_anchors"]["hs_metrics"] = "{sample}/{sample}.hs_metrics.txt"
# trimming QC: stage key follows the selected trimmer (三种去接头报告格式与
# 命名各异,各自解析器都已接入);锚点条件注入,PE/SE 命名不同的按样本类型拆
if trimmer == "cutadapt":
    qc_summary_config["qc_dirs"]["cutadapt"] = dna_trimmed_dir
    if dna_paired_samples:
        qc_summary_config["qc_anchors"]["cutadapt"] = "{sample}/{sample}.cutadapt_report.txt"
        qc_summary_config["qc_anchor_samples"]["cutadapt"] = dna_paired_samples
    if dna_single_samples:
        qc_summary_config["qc_dirs"]["cutadapt_single"] = dna_trimmed_dir
        qc_summary_config["qc_anchors"]["cutadapt_single"] = "{sample}/{sample}.single.cutadapt_report.txt"
        qc_summary_config["qc_anchor_samples"]["cutadapt_single"] = dna_single_samples
elif trimmer == "trimmomatic":
    qc_summary_config["qc_dirs"]["trimmomatic"] = dna_trimmed_dir
    qc_summary_config["qc_anchors"]["trimmomatic"] = "{sample}/trimmomatic_report.txt"
else:
    qc_summary_config["qc_dirs"]["trim_galore"] = dna_trimmed_dir
    if dna_paired_samples:
        qc_summary_config["qc_anchors"]["trim_galore"] = "{sample}/trimming_statistics_1.txt"
        qc_summary_config["qc_anchor_samples"]["trim_galore"] = dna_paired_samples
    if dna_single_samples:
        qc_summary_config["qc_dirs"]["trim_galore_single"] = dna_trimmed_dir
        qc_summary_config["qc_anchors"]["trim_galore_single"] = "{sample}/trimming_statistics.txt"
        qc_summary_config["qc_anchor_samples"]["trim_galore_single"] = dna_single_samples
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
logger.info(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as DR_DNA_qc_summary

# ============================================
# Final Report (PPT, Mutation.smk 同款)
# ============================================
# Wired unconditionally; node.py registers the pptx only when
# Params.report.enabled is set (downstream enabled pulls the upstream DAG).
Mutation_report_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/DNA",
    "logdir": f"{logdir}/group",
    "samples": dna_paired_samples + dna_single_samples,
    "paired_design": config.get("Params", {}).get("workflow", {}).get("paired_design", []),
    "qc_summary": qc_summary_config["outdir"],
    "somatic_dir": gatk_somatic_config["outdir"],
    "spectrum_dir": somatic_spectrum_config["outdir"],
    "burden_dir": burden_config["outdir"],
    "msi_dir": msi_config["outdir"],
    "fragment_dir": f"{fragment_size_config['outdir']}/fragment",
    # registered outputs -> DAG anchors inside the report module (module
    # filters them by its consumed dirs; only registered files enter, so
    # disabled stages stay best-effort 'not run' notes)
    "outfiles": outfiles,
    "Params": {
        "report": config.get("Params", {}).get("report", {}),
        # 生产者 formats 段直通(modules.md #9):report 脚本运行期对消费的
        # 上游图做 raster-first 变体解析(spectrum/signature 图归
        # somatic_spectrum、FragmentSize 归 fragment_size、msi_barplot 归 msi)
        "somatic_spectrum": config.get("Params", {}).get("somatic_spectrum", {}),
        "fragment_size": config.get("Params", {}).get("fragment_size", {}),
        "msi": config.get("Params", {}).get("msi", {}),
    }
}
module Mutation_report:
    snakefile: "../modules/Mutation_report/Mutation_report.smk"
    config: Mutation_report_config
logger.info(f"Mutation_report_config: {Mutation_report_config}")
use rule generate_report from Mutation_report as DR_DNA_generate_report

# ############################################
# 链 1: RNA 链 — 复刻 RNAseq.smk 全部分析功能:
#       trimmer(cutadapt/trimmomatic/trim_galore) -> TE 定量(hisat2|star
#       -> TEtranscripts) -> 定量(star -> featureCounts) -> 差异(DESeq2)
#       -> 富集(function GO/KEGG/GSEA/GSVA) -> 组装(hisat2 -> StringTie)
#       -> 融合(RmrRNA -> bowtie2 -> star chimera -> markdup -> arriba)
#       -> RNA SNP(star -> XenofilteR -> markdup -> gatk_RNAseq)
#       -> RNAseq_report + sexcheck + qc_summary
# 输出全部落在 {outdir}/RNA/ 下 (目录编号对齐 RNAseq.smk,带 RNA 前缀),
# 与链 0 的 {outdir}/DNA/ 互不干扰。
# 注: MetaUtil 按 meta data_type 列把 RNA 行软链整理到
#     {outdir}/RNA/common/1_raw_data/ 且命名为 <ID>_{1,2}.fq.gz,
#     直接消费 (indir 已被 node.py runDR 改写指向那里)。
# star 索引单独一个模块实例, index 写到 {outdir}/RNA/genome/index/{genome},
# 各 align 实例通过 index_dir 指向同一目录 (RNAseq.smk 的多实例模式)。
# RNA 差异对比 (DESeq2.group_pairs): node.py runDR 优先用 config 里手填的
# Params.DESeq2.group_pairs,否则按 RNA 行 group 列两两组合自动派生。
# ############################################

rna_star_config_for_index = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "outdir": f"{outdir}/RNA/genome",
    "logdir_index": os.path.join(logdir, "group"),
    "env": config.get("env", {}),
    "Procedure": {"star": config.get('Procedure', {}).get('star')},
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_star_for_index:
    snakefile: "../modules/star/polygenomes/star.smk"
    config: rna_star_config_for_index
use rule star_index from rna_star_for_index as DR_RNA_star_index

# --- trimmer 三选一 (RNAseq.smk 同款): cutadapt / trimmomatic / trim_galore ---
rna_trimmed_dir = f"{outdir}/RNA/common/2_trimmed_fastq"
if trimmer == "cutadapt":
    rna_cutadapt_config = {
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir": rna_trimmed_dir,
        "logdir": os.path.join(logdir, "sample"),
        "Procedure": {
            "cutadapt": config.get('Procedure', {}).get('cutadapt')
        },
        "Params": {
            # 参数精确到 rule:Params.cutadapt.<rule>.<key>(模块读点)
            "cutadapt": config.get('Params', {}).get("cutadapt", {}),
        },
        "resources": config.get("resources", {}),
    }
    module rna_cutadapt:
        snakefile: "../modules/cutadapt/cutadapt.smk"
        config: rna_cutadapt_config
    logger.debug(f"rna_cutadapt_config: {rna_cutadapt_config}")
    use rule trimming_Paired from rna_cutadapt as DR_RNA_cutadapt_Paired
    use rule trimming_Single from rna_cutadapt as DR_RNA_cutadapt_Single
elif trimmer == "trimmomatic":
    rna_trimmomatic_config = {
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir": rna_trimmed_dir,
        "logdir": os.path.join(logdir, "sample"),
        "Procedure": {
            "trimmomatic": config.get('Procedure', {}).get('trimmomatic')
        },
        "Params": {
            # 参数精确到 rule:整段透传(adapter_pe 挂 trimmomatic_Paired,
            # adapter_se 挂 trimmomatic_Single 各自 rule 子键)
            "trimmomatic": config.get('Params', {}).get("trimmomatic", {}),
        },
        "resources": config.get("resources", {}),
    }
    module rna_trimmomatic:
        snakefile: "../modules/trimmomatic/trimmomatic.smk"
        config: rna_trimmomatic_config
    logger.debug(f"rna_trimmomatic_config: {rna_trimmomatic_config}")
    use rule trimmomatic_Paired from rna_trimmomatic as DR_RNA_trimmomatic_Paired
    use rule trimmomatic_Single from rna_trimmomatic as DR_RNA_trimmomatic_Single
else:
    logger.info("use default trimmer: trim_galore")
    rna_trim_galore_config = {
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir": rna_trimmed_dir,
        "logdir": os.path.join(logdir, "sample"),
        "Procedure": {
            "trim_galore": config.get('Procedure', {}).get('trim_galore')
        },
        "Params": {
            # 参数精确到 rule:Params.trim_galore.<rule>.<key>(与 RNAseq.json 一致)
            "trim_galore": config.get('Params', {}).get("trim_galore", {}),
        },
        "resources": config.get("resources", {}),
    }
    module rna_trim_galore:
        snakefile: "../modules/trim-galore/trim-galore.smk"
        config: rna_trim_galore_config
    logger.debug(f"rna_trim_galore_config: {rna_trim_galore_config}")
    use rule trimming_Paired from rna_trim_galore as DR_RNA_trim_galore_Paired
    use rule trimming_Single from rna_trim_galore as DR_RNA_trim_galore_Single

# --- TE 定量: 比对器二选一 (hisat2 | star) -> TEtranscripts (TEcount/TElocal) ---
if aligner_TEtranscripts == 'hisat2':
    rna_hisat2_config_for_TEtranscripts = {
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "indir": rna_trimmed_dir,
        "outdir": f"{outdir}/RNA/common/3_TEtranscripts_bam",
        "env": config.get("env", {}),
        "logdir": os.path.join(logdir, "sample"),
        "logdir_combine": os.path.join(logdir, "group"),
        "genome_paired_samples": genome_paired_samples,
        "genome_single_samples": genome_single_samples,
        "Procedure": {
            "hisat2": config.get('Procedure', {}).get('hisat2'),
            "hisat2-build": config.get('Procedure', {}).get('hisat2_build')
        },
        "Params": {
            "hisat2": {
                "keep_intermediates": config.get('Params', {}).get('hisat2_TEtranscripts', {}).get('keep_intermediates', True),
                "score_min": config.get('Params', {}).get('hisat2_TEtranscripts', {}).get('score_min') or "L,0,-0.2",
                "flag_params": config.get('Params', {}).get('hisat2_TEtranscripts', {}).get('flag_params') or "--no-mixed --no-discordant",
                "k": config.get('Params', {}).get('hisat2_TEtranscripts', {}).get('k') or 100
            }
        },
        "genome": config.get('genome', {}),
        "resources": config.get("resources", {}),
    }
    module rna_hisat2_for_TEtranscripts:
        snakefile: "../modules/hisat2/polygenomes/hisat2.smk"
        config: rna_hisat2_config_for_TEtranscripts
    logger.debug(f"rna_hisat2_config_for_TEtranscripts: {rna_hisat2_config_for_TEtranscripts}")
    use rule hisat2_align from rna_hisat2_for_TEtranscripts as DR_RNA_hisat2_align_for_TEtranscripts
    use rule hisat2_index from rna_hisat2_for_TEtranscripts as DR_RNA_hisat2_index_for_TEtranscripts
    rna_te_align_outdir = rna_hisat2_config_for_TEtranscripts["outdir"]
elif aligner_TEtranscripts == 'star':
    rna_star_config_for_TEtranscripts = {
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "indir": rna_trimmed_dir,
        "outdir": f"{outdir}/RNA/common/3_TEtranscripts_bam",
        "env": config.get("env", {}),
        "logdir": os.path.join(logdir, "sample"),
        "logdir_index": os.path.join(logdir, "group"),
        "isGenomeSubdir": False,
        "index_dir": rna_star_config_for_index["outdir"],
        "genome_paired_samples": genome_paired_samples,
        "genome_single_samples": genome_single_samples,
        "Procedure": {
            "star": config.get('Procedure', {}).get('star')
        },
        "Params": {
            # 参数精确到 rule:Params.star.star_align.<key>(源读
            # STAR_TEtranscripts.star_align 新键,config 全量路径读取,
            # 无旧扁平键回退),缺省值保持 RNAseq.smk 同款默认
            "star": {
                "keep_intermediates": config.get('Params', {}).get('STAR_TEtranscripts', {}).get('keep_intermediates', True),
                "star_align": {
                    "alignEndsType": config.get('Params', {}).get('STAR_TEtranscripts', {}).get('star_align', {}).get('alignEndsType') or "Local",
                    "outFilterMultimapNmax": config.get('Params', {}).get('STAR_TEtranscripts', {}).get('star_align', {}).get('outFilterMultimapNmax') or 100,
                    "winAnchorMultimapNmax": config.get('Params', {}).get('STAR_TEtranscripts', {}).get('star_align', {}).get('winAnchorMultimapNmax') or 100,
                }
            }
        },
        "genome": config.get('genome', {}),
        "resources": config.get("resources", {}),
    }
    module rna_star_for_TEtranscripts:
        snakefile: "../modules/star/polygenomes/star.smk"
        config: rna_star_config_for_TEtranscripts
    logger.debug(f"rna_star_config_for_TEtranscripts: {rna_star_config_for_TEtranscripts}")
    use rule star_align from rna_star_for_TEtranscripts as DR_RNA_star_align_for_TEtranscripts
    rna_te_align_outdir = rna_star_config_for_TEtranscripts["outdir"]
else:
    raise ValueError(f"Unsupported aligner_TEtranscripts: {aligner_TEtranscripts}")

### gene and TE expression quantification using TEtranscripts
rna_TEtranscripts_config = {
    "indir": rna_te_align_outdir,
    "outdir": f"{outdir}/RNA/counts/TEtranscripts",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "genome_samples": genome_samples,
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "Params": {
        "TEtranscripts": {
            "keep_intermediates": config.get('Params',{}).get('TEtranscripts',{}).get('keep_intermediates', True)
        }
    },
    "Procedure": {
        "TEcount": config.get('Procedure', {}).get('TEcount') or 'TEcount',
        "TElocal": config.get('Procedure', {}).get('TElocal') or 'TElocal'
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_TEtranscripts:
    snakefile: "../modules/TEtranscripts/polygenomes/TEtranscripts.smk"
    config: rna_TEtranscripts_config
logger.debug(f"rna_TEtranscripts_config: {rna_TEtranscripts_config}")
use rule * from rna_TEtranscripts as DR_RNA_te_*

# --- 基因定量: star -> featureCounts (GSVA/定量表) ---
rna_star_config_for_featureCounts = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "indir": rna_trimmed_dir,
    "outdir": f"{outdir}/RNA/common/11_featureCounts_bam",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_index": os.path.join(logdir, "group"),
    "isGenomeSubdir": False,
    "index_dir": rna_star_config_for_index["outdir"],
    "genome_paired_samples": genome_paired_samples,
    "genome_single_samples": genome_single_samples,
    "env": config.get("env", {}),
    "Procedure": {
        "star": config.get('Procedure', {}).get('star')
    },
    "Params": {
        "star": config.get("Params", {}).get("STAR_featureCounts", {})
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_star_for_featureCounts:
    snakefile: "../modules/star/polygenomes/star.smk"
    config: rna_star_config_for_featureCounts
logger.debug(f"rna_star_config_for_featureCounts: {rna_star_config_for_featureCounts}")
use rule star_align from rna_star_for_featureCounts as DR_RNA_star_align_for_featureCounts

### gene expression quantification using featureCounts (for GSVA downstream)
rna_featureCounts_config = {
    "indir": rna_star_config_for_featureCounts["outdir"],
    "outdir": f"{outdir}/RNA/counts/featureCounts",
    "logdir": os.path.join(logdir, "group"),
    "genome_paired_samples": genome_paired_samples,
    "genome_single_samples": genome_single_samples,
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "Params": {
        "featureCounts": {
            "keep_intermediates": config.get('Params',{}).get('featureCounts',{}).get('keep_intermediates', True)
        }
    },
    "Procedure": {
        "featureCounts": config.get('Procedure', {}).get('featureCounts') or 'featureCounts'
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_featureCounts:
    snakefile: "../modules/featureCounts/polygenomes/featureCounts.smk"
    config: rna_featureCounts_config
logger.debug(f"rna_featureCounts_config: {rna_featureCounts_config}")
use rule featureCounts_paired_noMultiple from rna_featureCounts as DR_RNA_featureCounts_paired
use rule featureCounts_single_noMultiple from rna_featureCounts as DR_RNA_featureCounts_single
use rule featureCounts_merge from rna_featureCounts as DR_RNA_featureCounts_merge
use rule featureCounts_result from rna_featureCounts as DR_RNA_featureCounts_result

### differential expression analysis using DESeq2
# group_pairs 由 node.py runDR 写入 (config 手填优先, 否则按 group 列派生)
rna_DESeq2_config = {
    "indir": rna_TEtranscripts_config["outdir"],
    "outdir": f"{outdir}/RNA/diff_expression",
    "logdir": os.path.join(logdir, "group"),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "group_pairs": config.get("Params", {}).get("DESeq2", {}).get("group_pairs"),
    "Params": {
        "DESeq2": config.get("Params", {}).get("DESeq2", {})
    },
    "genome": config.get('genome', {}),
    "Procedure": {
        "DESeq2": config.get('Procedure', {}).get('DESeq2') or 'DESeq2'
    },
    "resources": config.get("resources", {}),
}
module rna_DESeq2:
    snakefile: "../modules/DESeq2/polygenomes/DESeq2.smk"
    config: rna_DESeq2_config
logger.debug(f"rna_DESeq2_config: {rna_DESeq2_config}")
use rule DESeq2_TEcount from rna_DESeq2 as DR_RNA_DESeq2_TEcount
use rule DESeq2_TElocal from rna_DESeq2 as DR_RNA_DESeq2_TElocal
use rule DESeq2_TElocal_manhattan from rna_DESeq2 as DR_RNA_DESeq2_TElocal_manhattan
use rule DESeq2_TElocal_locus_boxplot from rna_DESeq2 as DR_RNA_DESeq2_TElocal_locus_boxplot
use rule DESeq2_TElocal_subfamily_bar from rna_DESeq2 as DR_RNA_DESeq2_TElocal_subfamily_bar

### functional enrichment analysis using function module
rna_function_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": rna_DESeq2_config["outdir"],
    "counts_dir": rna_featureCounts_config["outdir"],
    "outdir": f"{outdir}/RNA/function",
    "logdir": os.path.join(logdir, "group"),
    "group_pairs": config.get("Params", {}).get("DESeq2", {}).get("group_pairs"),
    "genome": config.get('genome', {}),
    "Params": {
        "function": config.get("Params", {}).get("function", {})
    },
    "resources": config.get("resources", {}),
}
module rna_function:
    snakefile: "../modules/function/polygenomes/function.smk"
    config: rna_function_config
logger.debug(f"rna_function_config: {rna_function_config}")
use rule function_go_kegg from rna_function as DR_RNA_function_go_kegg
use rule function_gsea from rna_function as DR_RNA_function_gsea
use rule function_gsva from rna_function as DR_RNA_function_gsva

### transcript assembly using StringTie (hisat2 -> StringTie)
rna_hisat2_config_for_StringTie = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "indir": rna_trimmed_dir,
    "outdir": f"{outdir}/RNA/common/4_stringtie_bam",
    "env": config.get("env", {}),
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "genome_paired_samples": genome_paired_samples,
    "genome_single_samples": genome_single_samples,
    "Procedure": {
        "hisat2": config.get('Procedure', {}).get('hisat2')
    },
    "Params": {
        "hisat2": {
            "keep_intermediates": config.get('Params',{}).get('hisat2_stringtie', {}).get('keep_intermediates', True),
            "k": config.get('Params', {}).get('hisat2_stringtie', {}).get('k') or 1,
            "flag_params": config.get('Params', {}).get('hisat2_stringtie', {}).get('flag_params') or "-q --no-unal --dta",
        }
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_hisat2_for_StringTie:
    snakefile: "../modules/hisat2/polygenomes/hisat2.smk"
    config: rna_hisat2_config_for_StringTie
logger.debug(f"rna_hisat2_config_for_StringTie: {rna_hisat2_config_for_StringTie}")
use rule hisat2_align from rna_hisat2_for_StringTie as DR_RNA_hisat2_align_for_StringTie
use rule hisat2_index from rna_hisat2_for_StringTie as DR_RNA_hisat2_index_for_StringTie

rna_StringTie_config = {
    "indir": rna_hisat2_config_for_StringTie["outdir"],
    "outdir": f"{outdir}/RNA/transcripts",
    "env": config.get("env", {}),
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "genome_samples": genome_samples,
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "sample_groups": config.get("Params", {}).get("StringTie", {}).get('sample_groups'),
    "Params": {
        # wholesale section so formats (and future keys) reach the module;
        # StringTie.smk reads Params.StringTie.formats
        "StringTie": config.get("Params", {}).get("StringTie", {}),
    },
    "genome": config.get('genome', {}),
    "Procedure": {
        "stringtie": config.get("Procedure", {}).get("stringtie") or "stringtie"
    },
    "resources": config.get("resources", {}),
}
module rna_StringTie:
    snakefile: "../modules/StringTie/polygenomes/StringTie.smk"
    config: rna_StringTie_config
use rule * from rna_StringTie as DR_RNA_stringtie_*
logger.debug(f"rna_StringTie_config: {rna_StringTie_config}")

### fusion analysis: RmrRNA + bowtie2 -> star chimera -> markdup -> arriba
rna_rmrRNA_config = {
    "indir": rna_trimmed_dir,
    "outdir": f"{outdir}/RNA/genome",
    "logdir": os.path.join(logdir, "group"),
    "env": config.get("env", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "genome": config.get('genome', {}),
    "Params": {
        "RmrRNA": {
            "sam-append-comment": config.get('Params', {}).get('RmrRNA_rRNA', {}).get('sam-append-comment')
        }
    },
    "resources": config.get("resources", {}),
}
module rna_RmrRNA:
    snakefile: "../modules/RmrRNA/polygenomes/RmrRNA.smk"
    config: rna_rmrRNA_config
use rule * from rna_RmrRNA as DR_RNA_rmrna_*
logger.debug(f"rna_rmrRNA_config: {rna_rmrRNA_config}")

rna_bowtie2_rRNA_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": rna_trimmed_dir,
    "outdir": f"{outdir}/RNA/common/5_rRNA_fastq",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "genome_paired_samples": genome_paired_samples,
    "genome_single_samples": genome_single_samples,
    "Procedure": {
        "bowtie2": config.get('Procedure', {}).get('bowtie2'),
        "bowtie2-build": config.get('Procedure', {}).get('bowtie2-build')
    },
    "genome": config.get('genome', {}),
    "Params": {
        "bowtie2": {
            "keep_intermediates": config.get('Params',{}).get('bowtie2_rRNA', {}).get('keep_intermediates', True),
            "sam-append-comment": config.get('Params', {}).get('bowtie2_rRNA', {}).get('sam-append-comment')
        }
    },
    "resources": config.get("resources", {}),
}
module rna_bowtie2_for_rRNA:
    snakefile: "../modules/bowtie2/polygenomes/bowtie2.smk"
    config: rna_bowtie2_rRNA_config
use rule * from rna_bowtie2_for_rRNA as DR_RNA_bowtie2_*
logger.debug(f"rna_bowtie2_rRNA_config: {rna_bowtie2_rRNA_config}")

rna_star_config_for_fusion = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": rna_bowtie2_rRNA_config["outdir"],
    "outdir": f"{outdir}/RNA/common/6_fusion_bam",
    "logdir": os.path.join(logdir, "sample"),
    "isGenomeSubdir": True,
    "index_dir": rna_star_config_for_index["outdir"],
    "genome_paired_samples": genome_paired_samples,
    "genome_single_samples": genome_single_samples,
    "fastq_sample_suffix": "unmapped",
    "Procedure": {
        "STAR": config.get('Procedure', {}).get('STAR')
    },
    "Params": {
        # 参数精确到 rule:Params.star.star_align.<key>(源读 STAR_fusion.star_align
        # 新键,config 全量路径读取,无旧扁平键回退),缺省值保持 RNAseq.smk 同款默认
        "star": {
            "keep_intermediates": config.get('Params', {}).get('STAR_fusion', {}).get('keep_intermediates', True),
            "star_align": {
                "outFilterMultimapNmax": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('outFilterMultimapNmax') or 1,
                "outFilterMismatchNmax": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('outFilterMismatchNmax') or 3,
                "chimSegmentMin": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('chimSegmentMin') or 10,
                "chimOutType": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('chimOutType') or "Junctions WithinBAM SoftClip SeparateSAMold",
                "chimJunctionOverhangMin": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('chimJunctionOverhangMin') or 10,
                "outSAMstrandField": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('outSAMstrandField') or "intronMotif",
                "chimScoreMin": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('chimScoreMin') or 1,
                "chimScoreDropMax": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('chimScoreDropMax') or 30,
                "chimScoreJunctionNonGTAG": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('chimScoreJunctionNonGTAG') or 0,
                "chimScoreSeparation": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('chimScoreSeparation') or 1,
                "alignSJstitchMismatchNmax": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('alignSJstitchMismatchNmax') or "5 -1 5 5",
                "chimSegmentReadGapMax": config.get('Params', {}).get('STAR_fusion', {}).get('star_align', {}).get('chimSegmentReadGapMax') or 3,
            }
        }
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_star_for_fusion:
    snakefile: "../modules/star/polygenomes/star.smk"
    config: rna_star_config_for_fusion
use rule * from rna_star_for_fusion as DR_RNA_fusion_*
logger.debug(f"rna_star_config_for_fusion: {rna_star_config_for_fusion}")

rna_gatk_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": rna_star_config_for_fusion["outdir"],
    "outdir": f"{outdir}/RNA/common/7_markdup_bam",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_gatk_prepare_fusion:
    snakefile: "../modules/gatk/polygenomes/gatk_prepare.smk"
    config: rna_gatk_prepare_config
use rule * from rna_gatk_prepare_fusion as DR_RNA_fusion_*
logger.debug(f"rna_gatk_prepare_config: {rna_gatk_prepare_config}")

rna_arriba_config = {
    "indir": rna_gatk_prepare_config['outdir'],
    "outdir": f"{outdir}/RNA/fusion",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_substring": "sorted_markdup",
    "genome_samples": genome_samples,
    "genome": config.get('genome', {}),
    "Params": {
        "arriba": {
            "blacklist": config.get('Params', {}).get('arriba', {}).get('blacklist'),
            "known_fusions": config.get('Params', {}).get('arriba', {}).get('known_fusions'),
            "t": config.get('Params', {}).get('arriba', {}).get('t'),
            "d": config.get('Params', {}).get('arriba', {}).get('d'),
            "E": config.get('Params', {}).get('arriba', {}).get('E'),
            "p": config.get('Params', {}).get('arriba', {}).get('p'),
            # figure formats (arriba.smk reads Params.arriba.formats)
            "formats": config.get('Params', {}).get('arriba', {}).get('formats')
        }
    },
    "Procedure": {
        "arriba": config.get('Procedure', {}).get('arriba') or 'arriba'
    },
    "resources": config.get("resources", {}),
}
module rna_arriba:
    snakefile: "../modules/arriba/polygenomes/arriba.smk"
    config: rna_arriba_config
use rule * from rna_arriba as DR_RNA_fusion_*
logger.debug(f"rna_arriba_config: {rna_arriba_config}")

# --- RNA SNP calling: star -> XenofilteR -> markdup -> gatk_RNAseq ---
# star_for_SNP 用 mixed 样本表 (每个基因组下全部样本), XenofilteR 按
# sample_contamination (host/contaminant) 分离污染 reads (RNAseq.smk 同款)。
rna_star_config_for_SNP = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "indir": rna_trimmed_dir,
    "outdir": f"{outdir}/RNA/common/8_SNP_bam",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_index": os.path.join(logdir, "group"),
    "isGenomeSubdir": False,
    "index_dir": rna_star_config_for_index["outdir"],
    "genome_paired_samples": genome_paired_samples_mixed,
    "genome_single_samples": genome_single_samples_mixed,
    "env": config.get("env", {}),
    "Procedure": {
        "star": config.get('Procedure', {}).get('STAR_SNP')
    },
    # RNAseq.smk 同款:整段透传 + star.keep_intermediates 从 STAR_SNP 显式映射
    "Params": {**(config.get('Params') or {}),
               "star": {**((config.get('Params') or {}).get('star') or {}),
                        "keep_intermediates": (config.get('Params') or {}).get('STAR_SNP', {}).get('keep_intermediates', True)}},
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_star_for_SNP:
    snakefile: "../modules/star/polygenomes/star.smk"
    config: rna_star_config_for_SNP
logger.debug(f"rna_star_config_for_SNP: {rna_star_config_for_SNP}")
use rule star_align from rna_star_for_SNP as DR_RNA_star_align_for_SNP

rna_XenofilteR_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": rna_star_config_for_SNP["outdir"],
    "outdir": f"{outdir}/RNA/common/9_xenofilteR_bam",
    "logdir": os.path.join(logdir, "sample"),
    "Params": {
        "XenofilteR": config.get('Params', {}).get('XenofilteR', {})
    },
    "Procedure": {
        "XenofilteR": config.get('Procedure', {}).get('XenofilteR') or 'XenofilteR'
    },
    "resources": config.get("resources", {}),
}
module rna_XenofilteR:
    snakefile: "../modules/XenofilteR/XenofilteR.smk"
    config: rna_XenofilteR_config
use rule * from rna_XenofilteR as DR_RNA_xeno_*
logger.debug(f"rna_XenofilteR_config: {rna_XenofilteR_config}")

rna_gatk_prepare_SNP_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": rna_XenofilteR_config["outdir"],
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
module rna_gatk_prepare_SNP:
    snakefile: "../modules/gatk/polygenomes/gatk_prepare.smk"
    config: rna_gatk_prepare_SNP_config
use rule * from rna_gatk_prepare_SNP as DR_RNA_SNP_*
logger.debug(f"rna_gatk_prepare_SNP_config: {rna_gatk_prepare_SNP_config}")

rna_gatk_RNAseq_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": rna_gatk_prepare_SNP_config['outdir'],
    "outdir": f"{outdir}/RNA/variation/germline_snv_indel_RNAseq",
    "logdir": os.path.join(logdir, "sample"),
    "logdir_combine": os.path.join(logdir, "group"),
    "Params": {
        "gatk_RNAseq": {
            "keep_intermediates": config.get('Params',{}).get('gatk_RNAseq',{}).get('keep_intermediates', True)
        }
    },
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": config.get('genome', {}),
    "resources": config.get("resources", {}),
}
module rna_gatk_RNAseq:
    snakefile: "../modules/gatk/polygenomes/gatk_RNAseq/gatk_RNAseq.smk"
    config: rna_gatk_RNAseq_config
use rule * from rna_gatk_RNAseq as DR_RNA_SNP_*
logger.debug(f"rna_gatk_RNAseq_config: {rna_gatk_RNAseq_config}")

# --- report generation using RNAseq_report module ---
# Params.report 键名是模块读点; DR 里 DNA 链的 Mutation_report 共用
# Params.report 段,故 RNA 报告用独立的 Params.RNAseq_report 段映射进去,
# 两份报告的 title/genome/date 互不干扰。
rna_RNAseq_report_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": f"{outdir}/RNA",
    "outdir": os.path.join(outdir, "RNA", "results"),
    "logdir": logdir,
    "genome_samples": genome_samples,
    "genome_paired_samples": genome_paired_samples,
    "genome_single_samples": genome_single_samples,
    "Params": {
        "report": config.get("Params", {}).get("RNAseq_report", {}),
        "DESeq2": config.get("Params", {}).get("DESeq2", {}),
        # 生产者 formats 段直通(modules.md #9):report 模块 input 在 parse 期
        # 按上游 formats 做 raster-first 变体解析
        "StringTie": config.get("Params", {}).get("StringTie", {}),
        "arriba": config.get("Params", {}).get("arriba", {}),
        "function": config.get("Params", {}).get("function", {}),
    },
    "resources": config.get("resources", {}),
}
module rna_RNAseq_report:
    snakefile: "../modules/RNAseq_report/polygenomes/RNAseq_report.smk"
    config: rna_RNAseq_report_config
logger.debug(f"rna_RNAseq_report_config: {rna_RNAseq_report_config}")
use rule generate_report from rna_RNAseq_report as DR_RNA_generate_report

# ============================================
# Sex check (RNA marker-expression mode, RNAseq.smk 同款)
# ============================================
# 参考 RNAseq.smk:挂在 fusion 链 markdup BAM 上; Params.sexcheck_RNAseq
# 段映射到模块读点 Params.sexcheck (data_type=RNA 选择 marker 表达逻辑)。
rna_sexcheck_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/RNA/QC/sexcheck",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": [
        s
        for g in list(genome_paired_samples) + list(genome_single_samples)
        for s in genome_paired_samples.get(g, []) + genome_single_samples.get(g, [])
    ],
    "sample_bams": {
        s: f"{rna_gatk_prepare_config['outdir']}/{g}/{s}/{s}.sorted_markdup.bam"
        for g in list(genome_paired_samples) + list(genome_single_samples)
        for s in genome_paired_samples.get(g, []) + genome_single_samples.get(g, [])
    },
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        # multi-genome workflow: per-genome sex contigs {genome: [x, y]}
        "sex_chroms_map": {
            g: (refs.get("sex_chroms") or [None, None])[:2]
            for g, refs in (config.get("genome", {}) or {}).get("references", {}).items()
        },
        # RNA-mode marker lookup per genome {genome: gtf}
        "gtf_map": {
            g: refs.get("gtf")
            for g, refs in (config.get("genome", {}) or {}).get("references", {}).items()
            if refs.get("gtf")
        },
        # per-sample genome assignment (pure sample_id keys)
        "samples": {
            s: g
            for g in list(genome_paired_samples) + list(genome_single_samples)
            for s in genome_paired_samples.get(g, []) + genome_single_samples.get(g, [])
        },
    },
    "Params": {
        "sexcheck": config.get("Params", {}).get("sexcheck_RNAseq", {})
    },
}
module rna_sexcheck:
    snakefile: "../modules/sexcheck/sexcheck.smk"
    config: rna_sexcheck_config
logger.info(f"rna_sexcheck_config: {rna_sexcheck_config}")
use rule sexcheck_idxstats from rna_sexcheck as DR_RNA_sexcheck_idxstats
use rule sexcheck from rna_sexcheck as DR_RNA_sexcheck

# ============================================
# QC Summary (RNA cross-sample metric aggregation, RNAseq.smk 同款)
# ============================================
rna_qc_summary_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/RNA/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": [s for g in list(genome_paired_samples) + list(genome_single_samples) for s in genome_paired_samples.get(g, []) + genome_single_samples.get(g, [])],
    "qc_dirs": {
        "markdup": rna_gatk_prepare_config["outdir"],
        "sexcheck": rna_sexcheck_config["outdir"],
    },
    "qc_anchors": {
        "sexcheck": "sexcheck.tsv",
    },
}
# trimming QC: stage key follows the selected trimmer (RNA 侧 discovery-only)
if trimmer == "cutadapt":
    rna_qc_summary_config["qc_dirs"]["cutadapt"] = rna_trimmed_dir
elif trimmer == "trimmomatic":
    rna_qc_summary_config["qc_dirs"]["trimmomatic"] = rna_trimmed_dir
else:
    rna_qc_summary_config["qc_dirs"]["trim_galore"] = rna_trimmed_dir
module rna_qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: rna_qc_summary_config
logger.info(f"rna_qc_summary_config: {rna_qc_summary_config}")
use rule qc_summary from rna_qc_summary as DR_RNA_qc_summary

# ############################################
# 链 2: DRclass module (M3 pileup -> 状态机 -> O1/O2/O3)
# samples 表由 node.py runDR 显式生成(全部具体文件路径, DNA/RNA 产物对齐到
# 链 0/链 1 的输出目录), 无目录扫描.
# ############################################

drclass_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
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

# ############################################
# 链 3: TransEffect module (转录效应分析, 3.2.2.3 第四节):
#       DRclass O1 位点 + featureCounts 表达矩阵 + CNVkit 段
#       -> RNA 支持位点 VAF / RNA-DNA 比值 / 变异基因 Δ'(组内对照+随机 null)
#       / CNV 剂量一致性 / 稳定性程序基因反向检验 + 4 图
# 样本 = DRclass 已分类样本; sample_group_dict/gmt 由 node.py runDR 写入
# Params.TransEffect(node.py 仅在 Params.TransEffect.enabled 时注册输出,
# 输出决定流程走向)。
# ############################################

teffect_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/TransEffect",
    "logdir": f"{logdir}/group",
    "samples": sorted(drclass_samples),
    "sample_group_dict": config.get("Params", {}).get("TransEffect", {}).get("sample_group_dict", {}),
    "drclass_dir": f"{outdir}/DRclass",
    "counts_matrix": config.get("Params", {}).get("TransEffect", {}).get("counts_matrix")
        or (f"{outdir}/RNA/counts/featureCounts/{default_genome}/{default_genome}_featureCounts.tsv"
            if default_genome else ""),
    "cnv_dir": f"{outdir}/DNA/variation/cnv",
    "gmt": config.get("Params", {}).get("TransEffect", {}).get("gmt"),
    "gtf": genome_ref.get("gtf"),
    "Procedure": {
        "python3": config.get("Procedure", {}).get("python3") or "python3",
    },
    "Params": {
        "TransEffect": config.get("Params", {}).get("TransEffect", {}),
    },
}
module TransEffect:
    snakefile: "../modules/TransEffect/TransEffect.smk"
    config: teffect_config
use rule teffect_analyze from TransEffect as DR_teffect_analyze
use rule teffect_figures from TransEffect as DR_teffect_figures
logger.debug(f"teffect_config: {teffect_config}")

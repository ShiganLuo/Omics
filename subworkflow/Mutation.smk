shell.prefix("set -x; set -e;")
from snakemake.logging import logger
indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "logs")
outfiles = config.get("outfiles", [])
paired_samples = config.get("Params", {}).get("workflow", {}).get("paired_samples", [])
single_samples = config.get("Params", {}).get("workflow", {}).get("single_samples", [])
# Genome reference resolution - nested `default` + `references` pattern
default_genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(default_genome, {}) if default_genome else {}
ROOT_DIR = config.get("ROOT_DIR", ".")
rule all:
    input:
        outfiles

fastqc_raw_config = {
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir":  f"{outdir}/quality/fastqc/raw",
        "logdir": f"{logdir}/sample",
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
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir": f"{outdir}/common/2_trimmed_fastq",
        "logdir": f"{logdir}/sample",
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
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": cutadapt_config["outdir"],
        "outdir":  f"{outdir}/quality/fastqc/trimmed",
        "logdir": f"{logdir}/sample",
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
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": cutadapt_config["outdir"],
    "outdir":  f"{outdir}/common/3_sorted_bam",
    "logdir": f"{logdir}/sample",
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "bwaMem2": config.get("Procedure",{}).get("bwaMem2"),
        "samtools": config.get("Procedure",{}).get("samtools")
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "index_prefix": genome_ref.get("bwaMem2_index_prefix")
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
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": bwa_mem2_confg["outdir"],
    "outdir": bwa_mem2_confg["outdir"],
    "logdir": f"{logdir}/sample",
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
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": bwa_mem2_confg["outdir"],
    "outdir": f"{outdir}/common/4_markdup_bam",
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
        "dict_index": genome_ref.get("dict_index")
    }
}
module gatk_prepare:
    snakefile: "../modules/gatk/gatk_prepare.smk"
    config: gatk_prepare_config
logger.info(f"gatk_prepare parameters: {gatk_prepare_config}")
use rule gatk_index from gatk_prepare as Mutation_gatk_index
use rule addReadsGroup from gatk_prepare as Mutation_addReadsGroup
use rule MarkDuplicates from gatk_prepare as Mutation_MarkDuplicates

# Derive capture-region BED from markdup BAMs for exome/hybrid assays:
# consumed by GATK callers (-L) and CNVkit (--targets). Runs only when the
# assay is exome and no explicit genome.references.<g>.targets BED is set;
# node.py registers the output and fills the reference entry accordingly.
# Resolved by node.py: explicit BED path, the derived targets.bed, or None.
_targets_bed = genome_ref.get("targets") or None
targets_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/common/5_targets_bed",
    "logdir": f"{logdir}/group",
    "bam_suffix": ".sorted_markdup",
    "samples": paired_samples + single_samples,
    "Params": {
        "derive": config.get("Params", {}).get("derive_targets", {})
    },
}
module targets:
    snakefile: "../modules/targets/targets.smk"
    config: targets_config
logger.info(f"targets parameters: {targets_config}")
use rule derive_targets from targets as Mutation_derive_targets

gatk_somatic_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "tmp_dir": config.get("tmp_dir"),
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/somatic_snv_indel",
    "logdir": f"{logdir}/group",
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
    },
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
logger.info(f"gatk_somatic parameters: {gatk_somatic_config}")
use rule somaticMutect2 from gatk_somatic as Mutation_somaticMutect2
use rule germline_resource_strain from gatk_somatic as Mutation_germline_resource_strain
use rule filterMutectCalls from gatk_somatic as Mutation_filterMutectCalls

msi_pairs = {
    rec["experimental"]: rec["control"]
    for rec in config.get("Params", {}).get("workflow", {}).get("paired_design", [])
    if rec.get("experimental") and rec.get("control")
}
msi_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/stability/msi",
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
use rule msi_scan from msi as Mutation_msi_scan
use rule msi_call from msi as Mutation_msi_call
use rule msi_plot from msi as Mutation_msi_plot
use rule msi_result from msi as Mutation_msi_result

gatk_germline_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/germline_snv_indel",
    "logdir": f"{logdir}/sample",
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
    },
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
logger.info(f"gatk_germline parameters: {gatk_germline_config}")
use rule * from  gatk_germline as Mutation_*

# Wired unconditionally; node.py registers the spectrum output only when
# tumor/normal design pairs exist (outfiles drive the DAG).
somatic_spectrum_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": gatk_somatic_config["outdir"],
    "outdir": f"{outdir}/stability/spectrum",
    "logdir": f"{logdir}/group",
    "sample_somatic_vcf_dict": config.get("Params", {}).get("somatic_spectrum", {}).get("sample_somatic_vcf_dict", {}),
    "sample_group_dict": config.get("Params", {}).get("somatic_spectrum", {}).get("sample_group_dict", {}),
    "genome": {
        "fasta": genome_ref.get("fasta")
    }
}
module somatic_spectrum:
    snakefile: "../modules/spectrum/spectrum.smk"
    config: somatic_spectrum_config
logger.info(f"somatic_spectrum parameters: {somatic_spectrum_config}")
use rule somatic_spectrum from somatic_spectrum as Mutation_somatic_spectrum

# ============================================
# Fragment Size Analysis (cfDNA)
# ============================================
# Wired unconditionally; node.py controls whether it runs by registering
# (or omitting) the fragment-size outputs (Params.fragment_size.enabled).
fragment_size_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_suffix": ".sorted_markdup",
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/results/fragment_size",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
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

# ============================================
# SV Detection with Manta (cfDNA)
# ============================================
# Wired unconditionally; node.py controls whether it runs by registering
# (or omitting) the Manta outputs (Params.manta.enabled).
manta_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_suffix": ".sorted_markdup",
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/germline_sv",
    "logdir": f"{logdir}/sample",
    "samples": paired_samples + single_samples,
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
use rule manta_config_germline from manta as Mutation_manta_config_germline
use rule manta_run_germline from manta as Mutation_manta_run_germline

# Somatic (tumor/normal) SV calling via the same manta module, instantiated
# with a separate outdir/logdir. Wired unconditionally; node.py registers the
# somaticSV output only when matched-normal design pairs exist (mirroring the
# somatic SNV / Mutect2 pair judgment), which drives the DAG.
manta_somatic_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_suffix": ".sorted_markdup",
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/somatic_sv",
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
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
use rule manta_config_somatic from manta_somatic as Mutation_manta_config_somatic
use rule manta_run_somatic from manta_somatic as Mutation_manta_run_somatic

# ============================================
# CNV Detection with CNVkit (cfDNA)
# ============================================
# Wired unconditionally; node.py controls whether it runs by registering
# (or omitting) the CNVkit outputs (Params.cnvkit.enabled; control_samples
# additionally registers the reference panel output).
control_samples = config.get("Params", {}).get("workflow", {}).get("control_samples", [])
cnvkit_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "bam_suffix": ".sorted_markdup",
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/variation/cnv",
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
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
use rule cnvkit_targets from cnvkit as Mutation_cnvkit_targets
use rule cnvkit_reference from cnvkit as Mutation_cnvkit_reference
use rule cnvkit_batch from cnvkit as Mutation_cnvkit_batch
use rule cnvkit_paired from cnvkit as Mutation_cnvkit_paired

# ============================================
# Mutation Burden (SNV/InDel TMB + SV + CNV)
# ============================================
# Wired unconditionally; node.py controls whether it runs by registering
# (or omitting) the burden outputs (Params.burden.enabled). The three
# sample -> file maps are filled by node.py from the somatic/cnv outputs
# above; hand-filled config entries win.
burden_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/stability/burden",
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
use rule burden_snv_indel from burden as Mutation_burden_snv_indel
use rule burden_sv from burden as Mutation_burden_sv
use rule burden_cnv from burden as Mutation_burden_cnv
use rule burden_aggregate from burden as Mutation_burden_aggregate

# ============================================
# QC Summary (cross-sample metric aggregation)
# ============================================
# Wired unconditionally; node.py registers the qc_summary outputs so the
# aggregation always runs for the delivered output set (subworkflow.md
# 「QC 汇总与 report 模块规范」).
qc_summary_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "qc_dirs": {
        "fastqc_raw": fastqc_raw_config["outdir"],
        "fastqc_trimmed": fastqc_trimmed_config["outdir"],
        "cutadapt": cutadapt_config["outdir"],
        "flagstat": samtools_config["outdir"],
        "bwa_mem2": f"{logdir}/sample",
        "markdup": gatk_prepare_config["outdir"],
        "msi": msi_config["outdir"],
    },
}
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
logger.info(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as Mutation_qc_summary

# ============================================
# Final Report (PPT)
# ============================================
# Wired unconditionally; node.py registers the pptx only when
# Params.report.enabled is set (downstream enabled pulls the upstream DAG).
Mutation_report_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": outdir,
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "paired_design": config.get("Params", {}).get("workflow", {}).get("paired_design", []),
    "qc_summary": qc_summary_config["outdir"],
    "somatic_dir": gatk_somatic_config["outdir"],
    "spectrum_dir": somatic_spectrum_config["outdir"],
    "burden_dir": burden_config["outdir"],
    "msi_dir": msi_config["outdir"],
    "fragment_dir": f"{fragment_size_config['outdir']}/fragment",
    "Params": {
        "report": config.get("Params", {}).get("report", {})
    }
}
module Mutation_report:
    snakefile: "../modules/Mutation_report/Mutation_report.smk"
    config: Mutation_report_config
logger.info(f"Mutation_report_config: {Mutation_report_config}")
use rule generate_report from Mutation_report as Mutation_generate_report


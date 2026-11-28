shell.prefix("set -x; set -e;")
config.setdefault("OUTDIR_ROOT", config.get("Params", {}).get("workflow", {}).get("outdir", "output"))
include: "../modules/common/common.smk"
from snakemake.logging import logger
import os
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "logs")
genome_paired_samples = config.get("sample_params", {}).get("genome_paired_samples", {})
genome_single_samples = config.get("sample_params", {}).get("genome_single_samples", {})
dr_pairs = config.get("sample_params", {}).get("dr_pairs", {})
drclass_samples = config.get("sample_params", {}).get("drclass_samples", {})
dna_paired_samples = config.get("sample_params", {}).get("dna_paired_samples", [])
dna_single_samples = config.get("sample_params", {}).get("dna_single_samples", [])
dna_indir = config.get("Params", {}).get("workflow", {}).get("dna_indir", "")
trimmer = config.get('Params', {}).get('workflow', {}).get('trimmer') or "trim_galore"
aligner_TEtranscripts = config.get('Params', {}).get('workflow', {}).get('aligner_TEtranscripts') or "star"
outfiles = config.get("outfiles", [])
default_genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(default_genome, {}) if default_genome else {}
genomes = set(genome_paired_samples.keys()).union(set(genome_single_samples.keys()))
genome_samples = {}
for _g, _s in genome_paired_samples.items():
    genome_samples[_g] = list(_s)
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
    "tmp_dir": config.get("tmp_dir") or "/tmp",
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
    "pair_strains": config.get("sample_params", {}).get("pair_strains", {}),
    "design_pairs": config.get("sample_params", {}).get("paired_design", []),
    "Params": {
        "gatk": {
            "Mutect2": {
                "panel_of_normals": config.get("Params", {}).get("gatk", {}).get("Mutect2", {}).get("panel_of_normals")
            }
        },
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

msi_pairs = {
    rec["experimental"]: rec["control"]
    for rec in config.get("sample_params", {}).get("paired_design", [])
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
        "gatk_germline": {
            "keep_intermediates": config.get("Params", {}).get("gatk_germline", {}).get("keep_intermediates", True)
        }
    }
}
module gatk_germline:
    snakefile: "../modules/gatk/gatk_germline/gatk_germline.smk"
    config: gatk_germline_config
use rule * from gatk_germline as DR_DNA_germline_*

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

control_samples = config.get("sample_params", {}).get("control_samples", [])
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
    "paired_design": config.get("sample_params", {}).get("paired_design", []),
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
    "qc_anchor_samples": {
        "msi": [rec["experimental"] for rec in config.get("sample_params", {}).get("paired_design", []) if rec.get("experimental")],
    },
}
if _targets_bed:
    qc_summary_config["qc_anchors"]["hs_metrics"] = "{sample}/{sample}.hs_metrics.txt"
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

comparisons_view = design_comparisons(
    config.get("sample_params", {}).get("design") or {},
    config.get("sample_params", {}).get("organism_samples") or {})
_rna_only = {s for _ss in genome_paired_samples.values() for s in _ss} | \
            {s for _ss in genome_single_samples.values() for s in _ss}
comparisons_rna = {o: {n: c for n, c in _ps.items()
                       if set(c["control_samples"]) | set(c["experimental_samples"]) <= _rna_only}
                   for o, _ps in comparisons_view.items()}
rna_DESeq2_config = {
    "indir": rna_TEtranscripts_config["outdir"],
    "outdir": f"{outdir}/RNA/diff_expression",
    "logdir": os.path.join(logdir, "group"),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "comparisons": comparisons_rna,
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

rna_function_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": rna_DESeq2_config["outdir"],
    "counts_dir": rna_featureCounts_config["outdir"],
    "outdir": f"{outdir}/RNA/function",
    "logdir": os.path.join(logdir, "group"),
    "comparisons": comparisons_rna,
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
    "sample_groups": config.get("sample_params", {}).get("sample_groups"),
    "Params": {
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
        "sex_chroms_map": {
            g: (refs.get("sex_chroms") or [None, None])[:2]
            for g, refs in (config.get("genome", {}) or {}).get("references", {}).items()
        },
        "gtf_map": {
            g: refs.get("gtf")
            for g, refs in (config.get("genome", {}) or {}).get("references", {}).items()
            if refs.get("gtf")
        },
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


dr_report_contrasts = list(
    (comparisons_rna.get(default_genome) or {}).keys()
)
_dr_report_rna_samples = sorted(_all_paired | _all_single)
DR_report_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": outdir,
    "logdir": f"{logdir}/group",
    "samples": sorted(set(dna_paired_samples + dna_single_samples) | set(_dr_report_rna_samples)),
    "dna_samples": dna_paired_samples + dna_single_samples,
    "rna_samples": _dr_report_rna_samples,
    "paired_samples": sorted(_all_paired),
    "single_samples": sorted(_all_single),
    "paired_design": config.get("sample_params", {}).get("paired_design", []),
    "contrasts": dr_report_contrasts,
    "dna_dir": f"{outdir}/DNA",
    "rna_dir": f"{outdir}/RNA",
    "drclass_dir": drclass_config["outdir"],
    "teffect_dir": teffect_config["outdir"],
    "burden_dir": burden_config["outdir"],
    "genome_label": default_genome or "",
    "outfiles": outfiles,
    "Params": {
        "report": config.get("Params", {}).get("report", {}),
    },
}
module DR_report:
    snakefile: "../modules/DR_report/DR_report.smk"
    config: DR_report_config
logger.info(f"DR_report_config: {DR_report_config}")
use rule generate_report from DR_report as DR_generate_report
use rule report_result from DR_report as DR_report_result

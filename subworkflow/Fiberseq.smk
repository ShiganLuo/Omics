"""Fiber-seq subworkflow for single-molecule chromatin accessibility analysis.

Workflow:
  Step 1: fiberseq/1_prepare — predict-m6a (auto-skip if SPRQ) + add-nuc + fire
  Step 2: common/2_aligned  — pbmm2 align → aligned fiberseq BAM
  Step 3: fiberseq/3_extract — extract BED12 + call-peaks (FDR) + qc + fire BED
  Step 4: fiberseq/4_downstream — annotate, motif, footprint, nucpos, quant,
          tracks, coactuation, censat (Params.downstream.* gates in node.py)
  Step 5: fiberseq/5_haplotype — deepvariant + pbsv + hiphase → hap split
          extracts → haplotype diff (Params.downstream.haplotype.enabled)
  Step 6: enrichment — GO/KEGG of non-promoter peak genes
  Step 7: ldsc — partitioned heritability interface (Params.ldsc)

Reference:
  - Stergachis et al., 2020, Science (DOI: 10.1126/science.aaz1646)
  - Jha, Bohaczuk et al., 2024, Genome Research (DOI: 10.1101/gr.279095.124)
  - Single chromatin fiber profiling ..., Cell Rep Methods 2024 (DOI: 10.1016/j.crmeth.2024.100911)
  - DNA methylation influences human centromere positioning ..., Nat Genet 2025
    (DOI: 10.1038/s41588-025-02324-w)
  - https://fiberseq.github.io/
"""
shell.prefix("set -x; set -e;")
from snakemake.logging import logger

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
outfiles = config.get("outfiles") or []
samples = config.get("samples") or []

# ---- Resolve genome reference ----
genome_cfg = config.get("genome", {})
default_genome = genome_cfg.get("default", "")
genome_ref = genome_cfg.get("references", {}).get(default_genome, {})
fasta = genome_ref.get("fasta") or genome_cfg.get("fasta") or ""
logger.info(f"Fiber-seq genome: {default_genome}, fasta={fasta}")

rule all:
    input:
        outfiles


# ============================================================
# Step 1: Prepare fiberseq BAM → common/2_fiberseq_bam
#   predict-m6a (auto-skip if SPRQ) → add-nucleosomes → fire
# ============================================================
prepare_outdir = f"{outdir}/common/2_fiberseq_bam"

fibertools_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": indir,
    "outdir": prepare_outdir,
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "Params": {
        "fibertools": config.get("Params", {}).get("fibertools", {})
    },
    "genome": {
        "fasta": fasta
    }
}
module fibertools_prepare:
    snakefile: "../modules/fibertools/fibertools_prepare.smk"
    config: fibertools_prepare_config
logger.debug(f"Fiber-seq prepare config: {fibertools_prepare_config}")
use rule ft_predict_m6a from fibertools_prepare as Fiberseq_ft_predict_m6a
use rule ft_add_nucleosomes from fibertools_prepare as Fiberseq_ft_add_nucleosomes
use rule ft_fire from fibertools_prepare as Fiberseq_ft_fire


# ============================================================
# Step 2: Align to genome → common/3_align_bam
# ============================================================
aligned_dir = f"{outdir}/common/3_align_bam"

pbmm2_config = {
    "indir": prepare_outdir,
    "outdir": aligned_dir,
    "bam_suffix": ".fiberseq.fire.bam",
    "output_bam_suffix": ".fiberseq.fire.sorted",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "Procedure": {
        "pbmm2": config.get("Procedure", {}).get("pbmm2"),
        "samtools": config.get("Procedure", {}).get("samtools"),
    },
    "genome": {
        "fasta": fasta
    }
}
module pbmm2:
    snakefile: "../modules/pbmm2/pbmm2.smk"
    config: pbmm2_config
logger.debug(f"Fiber-seq pbmm2 config: {pbmm2_config}")
use rule pbmm2_align from pbmm2 as Fiberseq_pbmm2_align


# ============================================================
# Step 3: Extract BED + peaks + QC → results
# ============================================================
analysis_outdir = f"{outdir}/results"

fibertools_analysis_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": analysis_outdir,
    "bam_suffix": ".fiberseq.fire.sorted.bam",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "Params": {
        "fibertools": config.get("Params", {}).get("fibertools", {})
    },
    "genome": {
        "fasta": fasta,
        "fai": genome_ref.get("fai"),
    }
}
module fibertools_analysis:
    snakefile: "../modules/fibertools/fibertools_analysis.smk"
    config: fibertools_analysis_config
logger.debug(f"Fiber-seq analysis config: {fibertools_analysis_config}")
use rule ft_extract from fibertools_analysis as Fiberseq_ft_extract
use rule ft_call_peaks from fibertools_analysis as Fiberseq_ft_call_peaks
use rule ft_qc from fibertools_analysis as Fiberseq_ft_qc
use rule ft_extract_fire from fibertools_analysis as Fiberseq_ft_extract_fire


# ============================================================
# Step 4: Downstream interpretation (one module per analysis;
#   node.py gates outputs per Params.downstream.*.enabled)
# ============================================================
downstream = config.get("Params", {}).get("downstream", {})


def _fs_cfg(name, ds_key, genome_keys=()):
    """Config dict for one fiberseq_* downstream module."""
    return {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": analysis_outdir,
        "outdir": f"{analysis_outdir}/{name}",
        "logdir": f"{logdir}/sample",
        "samples": samples,
        "Params": {name: downstream.get(ds_key, {})},
        "genome": {k: genome_ref.get(k) for k in genome_keys},
    }


fiberseq_annotate_config = _fs_cfg("fiberseq_annotate", "annotate",
                                   ("gtf", "repeat_bed"))
module fiberseq_annotate:
    snakefile: "../modules/fiberseq/annotate/annotate.smk"
    config: fiberseq_annotate_config
logger.debug(f"fiberseq_annotate config: {fiberseq_annotate_config}")
use rule fiberseq_annotate_run from fiberseq_annotate as Fiberseq_annotate_run

fiberseq_motif_config = _fs_cfg("fiberseq_motif", "motif", ("fasta", "motif_file"))
module fiberseq_motif:
    snakefile: "../modules/fiberseq/motif/motif.smk"
    config: fiberseq_motif_config
logger.debug(f"fiberseq_motif config: {fiberseq_motif_config}")
use rule fiberseq_motif_scan from fiberseq_motif as Fiberseq_motif_scan
use rule fiberseq_motif_profile from fiberseq_motif as Fiberseq_motif_profile

fiberseq_footprint_config = _fs_cfg("fiberseq_footprint", "footprint")
module fiberseq_footprint:
    snakefile: "../modules/fiberseq/footprint/footprint.smk"
    config: fiberseq_footprint_config
logger.debug(f"fiberseq_footprint config: {fiberseq_footprint_config}")
use rule fiberseq_footprint_run from fiberseq_footprint as Fiberseq_footprint_run

fiberseq_nucpos_config = _fs_cfg("fiberseq_nucpos", "nucpos", ("gtf",))
module fiberseq_nucpos:
    snakefile: "../modules/fiberseq/nucpos/nucpos.smk"
    config: fiberseq_nucpos_config
logger.debug(f"fiberseq_nucpos config: {fiberseq_nucpos_config}")
use rule fiberseq_nucpos_run from fiberseq_nucpos as Fiberseq_nucpos_run

fiberseq_quant_config = _fs_cfg("fiberseq_quant", "quant")
module fiberseq_quant:
    snakefile: "../modules/fiberseq/quant/quant.smk"
    config: fiberseq_quant_config
logger.debug(f"fiberseq_quant config: {fiberseq_quant_config}")
use rule fiberseq_quant_run from fiberseq_quant as Fiberseq_quant_run
use rule fiberseq_quant_hap from fiberseq_quant as Fiberseq_quant_hap
use rule fiberseq_quant_merge from fiberseq_quant as Fiberseq_quant_merge

fiberseq_diff_config = _fs_cfg("fiberseq_diff", "diff")
module fiberseq_diff:
    snakefile: "../modules/fiberseq/diff/diff.smk"
    config: fiberseq_diff_config
logger.debug(f"fiberseq_diff config: {fiberseq_diff_config}")
use rule fiberseq_diff_run from fiberseq_diff as Fiberseq_diff_run

fiberseq_coactuation_config = _fs_cfg("fiberseq_coactuation", "coactuation")
module fiberseq_coactuation:
    snakefile: "../modules/fiberseq/coactuation/coactuation.smk"
    config: fiberseq_coactuation_config
logger.debug(f"fiberseq_coactuation config: {fiberseq_coactuation_config}")
use rule fiberseq_coactuation_run from fiberseq_coactuation as Fiberseq_coactuation_run

fiberseq_haplotype_config = _fs_cfg("fiberseq_haplotype", "haplotype")
module fiberseq_haplotype:
    snakefile: "../modules/fiberseq/haplotype/haplotype.smk"
    config: fiberseq_haplotype_config
logger.debug(f"fiberseq_haplotype config: {fiberseq_haplotype_config}")
use rule fiberseq_haplotype_diff from fiberseq_haplotype as Fiberseq_haplotype_diff

fiberseq_censat_config = _fs_cfg("fiberseq_censat", "censat",
                                 ("censat_bed", "mcpg_pileup"))
module fiberseq_censat:
    snakefile: "../modules/fiberseq/censat/censat.smk"
    config: fiberseq_censat_config
logger.debug(f"fiberseq_censat config: {fiberseq_censat_config}")
use rule fiberseq_censat_run from fiberseq_censat as Fiberseq_censat_run

fiberseq_track_config = _fs_cfg("fiberseq_track", "track", ("fai",))
module fiberseq_track:
    snakefile: "../modules/fiberseq/track/track.smk"
    config: fiberseq_track_config
logger.debug(f"fiberseq_track config: {fiberseq_track_config}")
use rule fiberseq_track_run from fiberseq_track as Fiberseq_track_run


# ============================================================
# Step 5: Haplotype chain (variant calling → phasing → per-hap analysis)
# ============================================================
variation_outdir = f"{outdir}/variation"

deepvariant_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": f"{variation_outdir}/germline_snv_indel",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_substring": "fiberseq.fire.sorted",
    "Procedure": {
        "deepvariant": config.get("Procedure", {}).get("deepvariant")
    },
    "Params": {
        "deepvariant": config.get("Params", {}).get("deepvariant", {})
    },
    "genome": {
        "fasta": fasta,
        "fai": genome_ref.get("fai"),
    }
}
module deepvariant:
    snakefile: "../modules/deepvariant/deepvariant.smk"
    config: deepvariant_config
logger.debug(f"Fiber-seq deepvariant config: {deepvariant_config}")
use rule deepvariant_run from deepvariant as Fiberseq_deepvariant_run

pbsv_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": f"{variation_outdir}/germline_sv",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_substring": "fiberseq.fire.sorted",
    "Procedure": {
        "pbsv": config.get("Procedure", {}).get("pbsv")
    },
    "genome": {
        "fasta": fasta
    }
}
module pbsv:
    snakefile: "../modules/pbsv/pbsv.smk"
    config: pbsv_config
logger.debug(f"Fiber-seq pbsv config: {pbsv_config}")
use rule pbsv_discover from pbsv as Fiberseq_pbsv_discover
use rule pbsv_call from pbsv as Fiberseq_pbsv_call

hiphase_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": f"{variation_outdir}/phased",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_dir": aligned_dir,
    "input_bam_substring": "fiberseq.fire.sorted",
    "vcf_dir": f"{variation_outdir}/germline_snv_indel",
    "input_vcf_substring": "",
    "output_substring": "",
    "Procedure": {
        "hiphase": config.get("Procedure", {}).get("hiphase")
    },
    "genome": {
        "fasta": fasta
    }
}
module hiphase:
    snakefile: "../modules/hiphase/hiphase.smk"
    config: hiphase_config
logger.debug(f"Fiber-seq hiphase config: {hiphase_config}")
use rule hiphase_phase from hiphase as Fiberseq_hiphase_phase

# Per-hap element extraction from the phased BAM
fibertools_hap_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": f"{variation_outdir}/phased",
    "outdir": analysis_outdir,
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_suffix": ".phased.bam",
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft",
        "samtools": config.get("Procedure", {}).get("samtools"),
    },
    "Params": {
        "fibertools": config.get("Params", {}).get("fibertools", {})
    },
    "genome": {
        "fasta": fasta
    }
}
module fibertools_hap:
    snakefile: "../modules/fibertools/fibertools_analysis.smk"
    config: fibertools_hap_config
logger.debug(f"Fiber-seq hap extract config: {fibertools_hap_config}")
use rule ft_split_hap from fibertools_hap as Fiberseq_ft_split_hap
use rule ft_extract_hap from fibertools_hap as Fiberseq_ft_extract_hap


# ============================================================
# Step 6: GO/KEGG + GSEA via the shared function module
#   (reuses function_go_kegg / function_gsea; input gene tables
#   are emitted by diff_accessibility in go-kegg.r format)
# ============================================================
function_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": f"{analysis_outdir}/fiberseq_diff",
    "outdir": f"{outdir}/function",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "group_pairs": config.get("group_pairs", {}),
    "genome": {
        "geneIDAnno": genome_ref.get("geneIDAnno"),
    },
    "Procedure": {
        "Rscript": config.get("Procedure", {}).get("Rscript")
    },
    "Params": {
        "function": config.get("Params", {}).get("function", {}),
    },
}
module function:
    snakefile: "../modules/function/function.smk"
    config: function_config
logger.debug(f"Fiber-seq function config: {function_config}")
use rule function_go_kegg from function as Fiberseq_function_go_kegg
use rule function_gsea from function as Fiberseq_function_gsea


# ============================================================
# Step 7: LDSC partitioned heritability (interface)
# ============================================================
ldsc_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/ldsc",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "ldsc": config.get("Procedure", {}).get("ldsc"),
        "munge_sumstats": config.get("Procedure", {}).get("munge_sumstats"),
        "make_annot": config.get("Procedure", {}).get("make_annot"),
    },
    "Params": {
        "ldsc": config.get("Params", {}).get("ldsc", {}),
    },
}
module ldsc:
    snakefile: "../modules/ldsc/ldsc.smk"
    config: ldsc_config
logger.debug(f"Fiber-seq ldsc config: {ldsc_config}")
use rule ldsc_annot from ldsc as Fiberseq_ldsc_annot
use rule ldsc_munge from ldsc as Fiberseq_ldsc_munge
use rule ldsc_h2_cts from ldsc as Fiberseq_ldsc_h2_cts
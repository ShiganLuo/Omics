shell.prefix("set -x; set -e;")
config.setdefault("OUTDIR_ROOT", config.get("Params", {}).get("workflow", {}).get("outdir", "output"))
include: "../modules/common/common.smk"
comparisons_view = design_comparisons(
    config.get("sample_params", {}).get("design") or {},
    config.get("sample_params", {}).get("organism_samples") or {})
comparisons_flat = {n: c for _oc in comparisons_view.values() for n, c in _oc.items()}
from snakemake.logging import logger

ROOT_DIR = config.get("ROOT_DIR", ".")
outfiles = config.get("outfiles") or []
indir = config.get("Params",{}).get("workflow",{}).get("indir") or "input"
outdir = config.get("Params",{}).get("workflow",{}).get("outdir") or "output"
logdir = config.get("Params",{}).get("workflow",{}).get("logdir") or "log"
samples = config.get("sample_params", {}).get("samples") or []

genome_cfg = config.get("genome", {})
default_genome = genome_cfg.get("default", "")
genome_ref = genome_cfg.get("references", {}).get(default_genome, {})
fasta = genome_ref.get("fasta") or genome_cfg.get("fasta") or ""
logger.debug(f"Fiber-seq genome: {default_genome}, fasta={fasta}")

rule all:
    input:
        outfiles

genome_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "logdir": f"{logdir}/group",
    "outdir": f"{outdir}/genome",
    "genome": {
        "fasta": fasta,
    },
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    }
}
module genome:
    snakefile: "../modules/genome/genome.smk"
    config: genome_config
use rule genome_index from genome as Fiberseq_genome_index

prepare_outdir = f"{outdir}/common/2_fiberseq_bam"

fibertools_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
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
    },
    "resources": config.get("resources", {}),
}
module fibertools_prepare:
    snakefile: "../modules/fibertools/fibertools_prepare.smk"
    config: fibertools_prepare_config
logger.debug(f"Fiber-seq prepare config: {fibertools_prepare_config}")
use rule ft_predict_m6a from fibertools_prepare as Fiberseq_ft_predict_m6a
use rule ft_add_nucleosomes from fibertools_prepare as Fiberseq_ft_add_nucleosomes
use rule ft_fire from fibertools_prepare as Fiberseq_ft_fire

aligned_dir = f"{outdir}/common/3_align_bam"

pbmm2_config = {
    "resources": config.get("resources", {}),
    "indir": prepare_outdir,
    "outdir": aligned_dir,
    "bam_suffix": ".fiberseq.fire.bam",
    "output_bam_suffix": ".fiberseq.fire.sorted",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
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

sexcheck_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/sexcheck",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": samples,
    "sample_bams": {
        s: f"{pbmm2_config['outdir']}/{s}/{s}.fiberseq.fire.sorted.bam"
        for s in samples
    },
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "sex_chroms": (config.get("genome", {}) or {}).get("references", {})
            .get(config.get("genome", {}).get("default"), {})
            .get("sex_chroms")
    },
    "Params": {
        "sexcheck": config.get("Params", {}).get("sexcheck", {})
    },
}
module sexcheck:
    snakefile: "../modules/sexcheck/sexcheck.smk"
    config: sexcheck_config
logger.debug(f"sexcheck_config: {sexcheck_config}")
use rule sexcheck_idxstats from sexcheck as Fiberseq_sexcheck_idxstats
use rule sexcheck from sexcheck as Fiberseq_sexcheck

qc_summary_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": samples,
    "qc_dirs": {
        "sexcheck": sexcheck_config["outdir"],
    },
    "qc_anchors": {
        "sexcheck": "sexcheck.tsv",
    },
}
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
logger.debug(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as Fiberseq_qc_summary

extract_outdir = f"{outdir}/common/4_elements"

fibertools_analysis_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": extract_outdir,
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
        "fai": genome_ref.get("fai_index"),
    },
    "tmp_dir": config.get("tmp_dir") or "/tmp",
}
module fibertools_analysis:
    snakefile: "../modules/fibertools/fibertools_analysis.smk"
    config: fibertools_analysis_config
logger.debug(f"Fiber-seq analysis config: {fibertools_analysis_config}")
use rule ft_extract from fibertools_analysis as Fiberseq_ft_extract
use rule ft_call_peaks from fibertools_analysis as Fiberseq_ft_call_peaks
use rule ft_qc from fibertools_analysis as Fiberseq_ft_qc
use rule ft_extract_fire from fibertools_analysis as Fiberseq_ft_extract_fire

analysis_outdir = f"{outdir}/fiberseq_analysis"
extract_config = fibertools_analysis_config
elements_dir = extract_config["outdir"]
phased_dir = f"{outdir}/common/5_phased_bam"

params = config.get("Params", {})

fiberseq_annotate_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "outdir": f"{analysis_outdir}/annotate",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_annotate": params.get("annotate", {})
    },
    "genome": {
        "gtf": genome_ref.get("gtf"),
        "repeat_bed": genome_ref.get("repeat_bed"),
    },
    "resources": config.get("resources", {}),
}
module fiberseq_annotate:
    snakefile: "../modules/fiberseq/annotate/annotate.smk"
    config: fiberseq_annotate_config
logger.debug(f"fiberseq_annotate config: {fiberseq_annotate_config}")
use rule fiberseq_annotate_run from fiberseq_annotate as Fiberseq_annotate_run

fiberseq_motif_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "outdir": f"{analysis_outdir}/motif",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_motif": params.get("motif", {})
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "motif_file": genome_ref.get("motif_file"),
    },
    "resources": config.get("resources", {}),
}
module fiberseq_motif:
    snakefile: "../modules/fiberseq/motif/motif.smk"
    config: fiberseq_motif_config
logger.debug(f"fiberseq_motif config: {fiberseq_motif_config}")
use rule fiberseq_motif_scan from fiberseq_motif as Fiberseq_motif_scan
use rule fiberseq_motif_profile from fiberseq_motif as Fiberseq_motif_profile

fiberseq_footprint_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "motif_dir": fiberseq_motif_config["outdir"],
    "outdir": f"{analysis_outdir}/footprint",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_footprint": params.get("footprint", {})
    },
    "genome": {},
    "resources": config.get("resources", {}),
}
module fiberseq_footprint:
    snakefile: "../modules/fiberseq/footprint/footprint.smk"
    config: fiberseq_footprint_config
logger.debug(f"fiberseq_footprint config: {fiberseq_footprint_config}")
use rule fiberseq_footprint_run from fiberseq_footprint as Fiberseq_footprint_run

fiberseq_nucpos_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "outdir": f"{analysis_outdir}/nucpos",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_nucpos": params.get("nucpos", {})
    },
    "genome": {
        "gtf": genome_ref.get("gtf"),
    },
    "resources": config.get("resources", {}),
}
module fiberseq_nucpos:
    snakefile: "../modules/fiberseq/nucpos/nucpos.smk"
    config: fiberseq_nucpos_config
logger.debug(f"fiberseq_nucpos config: {fiberseq_nucpos_config}")
use rule fiberseq_nucpos_run from fiberseq_nucpos as Fiberseq_nucpos_run

fiberseq_quant_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "hap_dir": phased_dir,
    "outdir": f"{analysis_outdir}/quant",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_quant": params.get("quant", {})
    },
    "genome": {},
    "resources": config.get("resources", {}),
}
module fiberseq_quant:
    snakefile: "../modules/fiberseq/quant/quant.smk"
    config: fiberseq_quant_config
logger.debug(f"fiberseq_quant config: {fiberseq_quant_config}")
use rule fiberseq_quant_run from fiberseq_quant as Fiberseq_quant_run
use rule fiberseq_quant_hap from fiberseq_quant as Fiberseq_quant_hap
use rule fiberseq_quant_merge from fiberseq_quant as Fiberseq_quant_merge

sigtest_dir = f"{analysis_outdir}/sigtest"

fiberseq_mods_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "outdir": f"{analysis_outdir}/mods",
    "logdir": f"{logdir}/group",
    "samples": samples,
    "te_bed": genome_ref.get("repeat_bed"),
    "tmp_dir": config.get("tmp_dir") or "/tmp",
    "sigtest_dir": sigtest_dir,
    "contrasts": (config.get("Params", {}).get("fiberseq_mods", {})
                  or {}).get("contrasts", []),
    "Params": {
        "fiberseq_mods": (params.get("fiberseq_mods") or {})
    },
    "resources": config.get("resources", {}),
}
module fiberseq_mods:
    snakefile: "../modules/fiberseq/mods/mods.smk"
    config: fiberseq_mods_config
logger.debug(f"fiberseq_mods config: {fiberseq_mods_config}")
use rule fiberseq_mods_longtable from fiberseq_mods as Fiberseq_mods_longtable
use rule fiberseq_mods_split from fiberseq_mods as Fiberseq_mods_split
use rule fiberseq_mods_top_loci from fiberseq_mods as Fiberseq_mods_top_loci

fiberseq_sigtest_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": sigtest_dir,
    "logdir": f"{logdir}/group",
    "sigtest_table": (config.get("Params", {}).get("sigtest", {})
                      or {}).get("sigtest_table", ""),
    "contrasts": (config.get("Params", {}).get("sigtest", {})
                  or {}).get("contrasts", []),
    "Params": {
        "sigtest": params.get("sigtest", {})
    },
    "resources": config.get("resources", {}),
}
module sigtest:
    snakefile: "../modules/sigtest/sigtest.smk"
    config: fiberseq_sigtest_config
logger.debug(f"fiberseq_sigtest config: {fiberseq_sigtest_config}")
use rule sigtest_run from sigtest as Fiberseq_sigtest_run
use rule sigtest_result from sigtest as Fiberseq_sigtest_result

fiberseq_diff_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "quant_dir": fiberseq_quant_config["outdir"],
    "annot_dir": fiberseq_annotate_config["outdir"],
    "outdir": f"{analysis_outdir}/diff",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_diff": params.get("diff", {})
    },
    "genome": {},
    "resources": config.get("resources", {}),
}
module fiberseq_diff:
    snakefile: "../modules/fiberseq/diff/diff.smk"
    config: fiberseq_diff_config
logger.debug(f"fiberseq_diff config: {fiberseq_diff_config}")
use rule fiberseq_diff_run from fiberseq_diff as Fiberseq_diff_run

fiberseq_coactuation_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "outdir": f"{analysis_outdir}/coactuation",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_coactuation": params.get("coactuation", {})
    },
    "genome": {},
    "resources": config.get("resources", {}),
}
module fiberseq_coactuation:
    snakefile: "../modules/fiberseq/coactuation/coactuation.smk"
    config: fiberseq_coactuation_config
logger.debug(f"fiberseq_coactuation config: {fiberseq_coactuation_config}")
use rule fiberseq_coactuation_run from fiberseq_coactuation as Fiberseq_coactuation_run

fiberseq_haplotype_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "quant_dir": fiberseq_quant_config["outdir"],
    "outdir": f"{analysis_outdir}/haplotype",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_haplotype": params.get("haplotype", {})
    },
    "genome": {},
    "resources": config.get("resources", {}),
}
module fiberseq_haplotype:
    snakefile: "../modules/fiberseq/haplotype/haplotype.smk"
    config: fiberseq_haplotype_config
logger.debug(f"fiberseq_haplotype config: {fiberseq_haplotype_config}")
use rule fiberseq_haplotype_diff from fiberseq_haplotype as Fiberseq_haplotype_diff

fiberseq_censat_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "outdir": f"{analysis_outdir}/censat",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_censat": params.get("censat", {})
    },
    "genome": {
        "censat_bed": genome_ref.get("censat_bed"),
        "mcpg_pileup": genome_ref.get("mcpg_pileup"),
    },
    "resources": config.get("resources", {}),
}
module fiberseq_censat:
    snakefile: "../modules/fiberseq/censat/censat.smk"
    config: fiberseq_censat_config
logger.debug(f"fiberseq_censat config: {fiberseq_censat_config}")
use rule fiberseq_censat_run from fiberseq_censat as Fiberseq_censat_run

fiberseq_telocen_norm_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": f"{analysis_outdir}/telocen_norm",
    "logdir": f"{logdir}/group",
    "samples": samples,
    "inputs": {
        "telocen_norm": {
            "bam": {s: f"{aligned_dir}/{s}/{s}.fiberseq.fire.sorted.bam"
                    for s in samples},
        },
    },
    "Params": {
        "fiberseq_telocen_norm": params.get("telocen_norm", {})
    },
    "genome": {
        "telocen_bed": genome_ref.get("telocen_bed"),
    },
    "resources": config.get("resources", {}),
}
module fiberseq_telocen_norm:
    snakefile: "../modules/fiberseq/telocen_norm/telocen_norm.smk"
    config: fiberseq_telocen_norm_config
logger.debug(f"fiberseq_telocen_norm config: {fiberseq_telocen_norm_config}")
use rule fiberseq_telocen_norm_run from fiberseq_telocen_norm as Fiberseq_telocen_norm_run
use rule fiberseq_telocen_norm_fig from fiberseq_telocen_norm as Fiberseq_telocen_norm_fig

fiberseq_epi_drift_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": f"{analysis_outdir}/epi_drift",
    "logdir": f"{logdir}/group",
    "samples": samples,
    "inputs": {
        "epi_drift": {
            "bam": {s: f"{aligned_dir}/{s}/{s}.fiberseq.fire.sorted.bam"
                    for s in samples},
        },
    },
    "Params": {
        "fiberseq_epi_drift": params.get("epi_drift", {})
    },
    "genome": {
        "fai_index": genome_ref.get("fai_index"),
    },
    "resources": config.get("resources", {}),
}
module fiberseq_epi_drift:
    snakefile: "../modules/fiberseq/epi_drift/epi_drift.smk"
    config: fiberseq_epi_drift_config
logger.debug(f"fiberseq_epi_drift config: {fiberseq_epi_drift_config}")
use rule fiberseq_epi_drift_run from fiberseq_epi_drift as Fiberseq_epi_drift_run
use rule fiberseq_epi_drift_fig from fiberseq_epi_drift as Fiberseq_epi_drift_fig

fiberseq_track_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": elements_dir,
    "outdir": f"{analysis_outdir}/track",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Params": {
        "fiberseq_track": params.get("track", {})
    },
    "genome": {
        "fasta": fasta,
        "fai": genome_ref.get("fai_index"),
    },
    "resources": config.get("resources", {}),
}
module fiberseq_track:
    snakefile: "../modules/fiberseq/track/track.smk"
    config: fiberseq_track_config
logger.debug(f"fiberseq_track config: {fiberseq_track_config}")
use rule fiberseq_track_run from fiberseq_track as Fiberseq_track_run

assembly_outdir = f"{outdir}/assembly"
hifiasm_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": "",
    "outdir": assembly_outdir,
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_substring": "",
    "inputs": {
        "hifiasm": {
            "bam": {s: f"{outdir}/common/1_raw_data/{s}/{s}.bam" for s in samples},
            "pbi": {s: f"{outdir}/common/1_raw_data/{s}/{s}.bam.pbi" for s in samples},
        },
    },
    "Params": {
        "hifiasm": {
            "keep_intermediates": config.get("Params", {}).get("hifiasm", {}).get("keep_intermediates", True),
            "qc": config.get("Params", {}).get("hifiasm", {}).get("qc", {}),
        }
    },
}
module hifiasm:
    snakefile: "../modules/hifiasm/hifiasm.smk"
    config: hifiasm_config
logger.debug(f"hifiasm_config: {hifiasm_config}")
use rule hifiasm_assemble from hifiasm as Fiberseq_hifiasm_assemble
use rule hifiasm_result from hifiasm as Fiberseq_hifiasm_result
use rule hifiasm_qc from hifiasm as Fiberseq_hifiasm_qc
use rule hifiasm_qc_result from hifiasm as Fiberseq_hifiasm_qc_result

centromere_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": assembly_outdir,
    "assembly_dir": hifiasm_config["outdir"],
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": samples,
    "Params": {
        "RepeatMasker": config.get("Params", {}).get("RepeatMasker", {}),
    },
}
module centromere:
    snakefile: "../modules/centromere/centromere.smk"
    config: centromere_config
logger.debug(f"centromere_config: {centromere_config}")
use rule * from centromere as Fiberseq_centromere_*

telomere_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": f"{outdir}/telomere/length",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": samples,
    "comparisons": comparisons_flat,
    "bam_substring": "fiberseq.fire.sorted",
    "assembly_dir": config.get("assembly_dir") or assembly_outdir,
    "genome": {"fasta": fasta, "subtel_ref": genome_ref.get("subtel_ref")},
    "inputs": {
        "telogator2": {
            "bam": {s: f"{aligned_dir}/{s}/{s}.fiberseq.fire.sorted.bam" for s in samples},
            "bai": {s: f"{aligned_dir}/{s}/{s}.fiberseq.fire.sorted.bai" for s in samples},
        },
        "read_density": {
            "bam": {s: f"{outdir}/common/1_raw_data/{s}/{s}.bam" for s in samples},
        },
    },
    "Procedure": {
        "telogator2": config.get("Procedure", {}).get("telogator2"),
        "tidk": config.get("Procedure", {}).get("tidk")
    },
    "Params": {
        "telomere": config.get("Params", {}).get("telomere", {})
    },
}
module telomere:
    snakefile: "../modules/telomere/length/length.smk"
    config: telomere_config
logger.debug(f"telomere_config: {telomere_config}")
use rule telogator2_run from telomere as Fiberseq_telogator2_run
use rule assembly_telomere_scan from telomere as Fiberseq_assembly_telomere_scan
use rule read_density_telomere from telomere as Fiberseq_read_density_telomere
use rule tidk_init from telomere as Fiberseq_tidk_init
use rule tidk_scan from telomere as Fiberseq_tidk_scan
use rule read_density_combined from telomere as Fiberseq_read_density_combined
use rule telomere_method_overview from telomere as Fiberseq_telomere_method_overview
use rule telomere_depth_qc from telomere as Fiberseq_telomere_depth_qc

telomere_composition_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/telomere/composition",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": samples,
    "inputs": {
        "composition": {
            "bam": {s: f"{aligned_dir}/{s}/{s}.fiberseq.fire.sorted.bam" for s in samples},
        },
    },
    "Procedure": {},
    "Params": {"telomere": config.get("Params", {}).get("telomere", {})},
}
module telomere_composition:
    snakefile: "../modules/telomere/composition/composition.smk"
    config: telomere_composition_config
use rule composition_run from telomere_composition as Fiberseq_composition_run

telomere_regions_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/telomere/regions",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "inputs": {
        "regions": {
            "scan_tsv": {s: f"{outdir}/telomere/length/assembly_scan/{s}/{s}_assembly_telomere.tsv"
                         for s in samples},
        },
    },
    "Procedure": {},
    "Params": {"telomere": config.get("Params", {}).get("telomere", {})},
}
module telomere_regions:
    snakefile: "../modules/telomere/regions/regions.smk"
    config: telomere_regions_config
use rule regions_telomere from telomere_regions as Fiberseq_regions_telomere

variation_outdir = f"{outdir}/variation"

deepvariant_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
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
        "fai": genome_ref.get("fai_index"),
    },
    "tmp_dir": config.get("tmp_dir") or "/tmp",
}
module deepvariant:
    snakefile: "../modules/deepvariant/deepvariant.smk"
    config: deepvariant_config
logger.debug(f"Fiber-seq deepvariant config: {deepvariant_config}")
use rule deepvariant_run from deepvariant as Fiberseq_deepvariant_run

pbsv_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": f"{variation_outdir}/germline_sv",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_substring": "fiberseq.fire.sorted",
    "Procedure": {
        "pbsv": config.get("Procedure", {}).get("pbsv")
    },
    "Params": {
        "pbsv": config.get("Params", {}).get("pbsv", {})
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

sawfish_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": f"{variation_outdir}/germline_cnv",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_substring": "fiberseq.fire.sorted",
    "Procedure": {
        "sawfish": config.get("Procedure", {}).get("sawfish"),
        "bgzip": config.get("Procedure", {}).get("bgzip"),
        "bcftools": config.get("Procedure", {}).get("bcftools"),
        "python": config.get("Procedure", {}).get("python")
    },
    "Params": {
        "sawfish": config.get("Params", {}).get("sawfish", {})
    },
    "genome": {
        "fasta": fasta,
        "fai_index": genome_ref.get("fai_index")
    }
}
module sawfish:
    snakefile: "../modules/sawfish/sawfish.smk"
    config: sawfish_config
logger.debug(f"Fiber-seq sawfish config: {sawfish_config}")
use rule sawfish_discover from sawfish as Fiberseq_sawfish_discover
use rule sawfish_joint_call from sawfish as Fiberseq_sawfish_joint_call
use rule sawfish_cnv_filter from sawfish as Fiberseq_sawfish_cnv_filter

hiphase_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": phased_dir,
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_dir": aligned_dir,
    "input_bam_substring": "fiberseq.fire.sorted",
    "vcf_dir": deepvariant_config["outdir"],
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

fibertools_hap_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": phased_dir,
    "outdir": phased_dir,
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

# function 模块对比视图:group 名排序定角色(匹配 diff_accessibility 约定)
def _fn_sorted_comparisons(view_src):
    view = {}
    for _oc in (view_src or {}).values():
        for _n, _c in _oc.items():
            _g1, _g2 = sorted((_c["control_group_name"], _c["experimental_group_name"]))
            _s1 = _c["control_samples"] if _c["control_group_name"] == _g1 else _c["experimental_samples"]
            _s2 = _c["experimental_samples"] if _s1 == _c["control_samples"] else _c["control_samples"]
            view[_n] = {"control_group_name": _g1, "experimental_group_name": _g2,
                        "control_samples": _s1, "experimental_samples": _s2}
    return view

function_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": fiberseq_diff_config["outdir"],
    "outdir": f"{outdir}/function",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "comparisons": _fn_sorted_comparisons(comparisons_view),
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

ldsc_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
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

bamCoverage_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": pbmm2_config["outdir"],
    "outdir": f"{outdir}/common/6_coverage_bw",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_substring": "fiberseq.fire.sorted",
    "bw_suffix": ".coverage.bw",
    "Procedure": {
        "bamCoverage": config.get("Procedure", {}).get("bamCoverage")
    },
    "Params": {
        "bamCoverage": config.get("Params", {}).get("bamCoverage", {})
    }
}
module bamCoverage:
    snakefile: "../modules/deeptools/bamCoverage.smk"
    config: bamCoverage_config
logger.debug(f"bamCoverage_config: {bamCoverage_config}")
use rule bamCoverage from bamCoverage as Fiberseq_bamCoverage

regionplot_params = config.get("Params", {}).get("regionplot", {})
regionplot_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/regionplot",
    "logdir": f"{logdir}/group",
    "samples": samples,
    "tracks": [{"label": "accessibility", "dir": f"{analysis_outdir}/track", "suffix": ".accessibility.bw"},
              {"label": "m6a_density", "dir": f"{analysis_outdir}/track", "suffix": ".m6a_density.bw"},
              {"label": "coverage", "dir": bamCoverage_config["outdir"], "suffix": ".coverage.bw"}],
    "Params": {
        "regionplot": regionplot_params
    },
    "genome": {
        "gtf": genome_ref.get("gtf"),
    },
    "resources": config.get("resources", {}),
}
module regionplot:
    snakefile: "../modules/regionplot/regionplot.smk"
    config: regionplot_config
logger.debug(f"regionplot_config: {regionplot_config}")
use rule regionplot from regionplot as Fiberseq_regionplot
use rule regionplot_result from regionplot as Fiberseq_regionplot_result

fiberbrowser_params = config.get("Params", {}).get("fiberbrowser", {})
fiberbrowser_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/fiberbrowser",
    "logdir": f"{logdir}/group",
    "samples": samples,
    "sample_sams": {
        s: f"{pbmm2_config['outdir']}/{s}/{s}.fiberseq.fire.sorted.bam"
        for s in samples
    },
    "Params": {
        "fiberbrowser": fiberbrowser_params
    },
    "resources": config.get("resources", {}),
}
module fiberbrowser:
    snakefile: "../modules/fiberbrowser/fiberbrowser.smk"
    config: fiberbrowser_config
logger.debug(f"fiberbrowser_config: {fiberbrowser_config}")
use rule fiberbrowser from fiberbrowser as Fiberseq_fiberbrowser
use rule fiberbrowser_result from fiberbrowser as Fiberseq_fiberbrowser_result

fibercluster_params = config.get("Params", {}).get("fibercluster", {})
fibercluster_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/fibercluster",
    "logdir": f"{logdir}/group",
    "samples": samples,
    "indir": extract_outdir,
    "regions_bed": fibercluster_params.get("regions_bed") or "",
    "settings": fibercluster_params.get("settings") or [],
    "Params": {
        "fibercluster": fibercluster_params
    },
    "resources": config.get("resources", {}),
}
module fibercluster:
    snakefile: "../modules/fibercluster/fibercluster.smk"
    config: fibercluster_config
logger.debug(f"fibercluster_config: {fibercluster_config}")
use rule fiber_cluster from fibercluster as Fiberseq_fiber_cluster
use rule fiber_report from fibercluster as Fiberseq_fiber_report
use rule fibercluster_result from fibercluster as Fiberseq_fibercluster_result

decoil_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": pbmm2_config["outdir"],
    "outdir": f"{outdir}/eccdna",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_substring": "fiberseq.fire.sorted",
    "vcf_dir": pbsv_config["outdir"],
    "vcf_suffix": ".sv.vcf.gz",
    "bw_dir": bamCoverage_config["outdir"],
    "bw_suffix": ".coverage.bw",
    "Procedure": {
        "decoil": config.get("Procedure", {}).get("decoil")
    },
    "Params": {
        "decoil": config.get("Params", {}).get("decoil", {}),
        "igv": config.get("Params", {}).get("igv", {})
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
use rule decoil_reconstruct from decoil as Fiberseq_decoil_reconstruct
use rule decoil_result from decoil as Fiberseq_decoil_result
use rule decoil_plot from decoil as Fiberseq_decoil_plot
use rule decoil_structure from decoil as Fiberseq_decoil_structure
use rule decoil_compare from decoil as Fiberseq_decoil_compare
use rule decoil_igv from decoil as Fiberseq_decoil_igv
use rule decoil_ucsc from decoil as Fiberseq_decoil_ucsc

sv_params = config.get("Params", {}).get("sv", {})
_mei = sv_params.get("mei", {})
if isinstance(_mei.get("te_buckets"), dict):
    sv_params = {**sv_params, "mei": {**_mei,
                 "te_buckets": _mei["te_buckets"].get(default_genome, [])}}

sv_outdir = f"{variation_outdir}/germline_sv_analysis"

sv_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": pbsv_config["outdir"],
    "outdir": sv_outdir,
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "comparisons": comparisons_flat,
    "vcf_dir": pbsv_config["outdir"],
    "vcf_suffix": ".sv.vcf.gz",
    "Params": {
        "sv": sv_params
    },
    "genome": {
        "fasta": genome_ref.get("fasta"),
        "gtf": genome_ref.get("gtf"),
        "fai": genome_ref.get("fai_index"),
        "custom_gtf": genome_ref.get("custom_gtf"),
        "custom_fasta": genome_ref.get("custom_fasta"),
        "te_lib": genome_ref.get("te_lib"),
        "repeat_bed": genome_ref.get("repeat_bed")
    }
}
module sv:
    snakefile: "../modules/sv/sv.smk"
    config: sv_config
logger.debug(f"Fiber-seq sv config: {sv_config}")
use rule sv_exp_specific from sv as Fiberseq_sv_exp_specific
use rule sv_mei from sv as Fiberseq_sv_mei
use rule sv_tra_te from sv as Fiberseq_sv_tra_te
use rule sv_exp_circos from sv as Fiberseq_sv_exp_circos
use rule sv_gene_model from sv as Fiberseq_sv_gene_model
use rule sv_diff_analysis from sv as Fiberseq_sv_diff_analysis
use rule sv_oncoprint from sv as Fiberseq_sv_oncoprint
use rule sv_result from sv as Fiberseq_sv_result

svgokegg_params = config.get("Params", {}).get("function", {}).get("sv_go_kegg", {})
sv_go_kegg_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": sv_outdir,
    "outdir": sv_outdir,
    "logdir": f"{logdir}/group",
    "comparisons": sv_config.get("comparisons", {}),
    "Params": {
        "function": {
            "sv_go_kegg": svgokegg_params,
        }
    },
}
module sv_go_kegg:
    snakefile: "../modules/function/sv_go_kegg/sv_go_kegg.smk"
    config: sv_go_kegg_config
logger.debug(f"Fiber-seq sv_go config: {sv_go_kegg_config}")
use rule sv_go_kegg from sv_go_kegg as Fiberseq_sv_go_kegg

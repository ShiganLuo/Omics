shell.prefix("set -x; set -e;")
from snakemake.logging import logger
ROOT_DIR = config.get("ROOT_DIR", ".")
_wf = config.get("Params", {}).get("workflow", {})
indir = _wf.get("indir") or "data/fastq"
outdir = _wf.get("outdir") or "output"
logdir = _wf.get("logdir") or "logs"
paired_samples = _wf.get("paired_samples") or []
single_samples = _wf.get("single_samples") or []
genome_paired_samples = _wf.get("genome_paired_samples") or {}
genome_single_samples = _wf.get("genome_single_samples") or {}
genome_pairs = _wf.get("genome_pairs") or []
genome = config.get("genome", {})
outfiles = config.get("outfiles", [])
rule all:
    input:
        outfiles

aligner = config.get('Procedure',{}).get('aligner')
SOAPnuke_cofig = {
    "resources": config.get("resources", {}),
        "indir": indir,
        "outdir":  f"{outdir}/common/2_trimmed_fastq",
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "logdir": logdir
}
module SOAPnuke:
    snakefile: "../modules/SOAPnuke/SOAPnuke.smk"
    config: SOAPnuke_cofig
logger.info(f"SOAPnuke config: {SOAPnuke_cofig}")
use rule soapnuke_filter_paired from SOAPnuke as CoCulture_soapnuke_filter_paired
use rule soapnuke_filter_single from SOAPnuke as CoCulture_soapnuke_filter_single

    
hisat2_config = {

    
    "resources": config.get("resources", {}),
        "indir": SOAPnuke_cofig["outdir"],
        "outdir":  f"{outdir}/common/3_raw_bam",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "genome_paired_samples": genome_paired_samples,
        "genome_single_samples": genome_single_samples,
        "Procedure": {
            "hisat2": config.get('Procedure',{}).get('hisat2')
        },
        # polygenomes modules resolve refs via genome.references.<genome>.fasta
        "genome": {"references": genome.get("references", {})}
    }

module hisat2:
    snakefile: "../modules/hisat2/polygenomes/hisat2.smk"
    config: hisat2_config
logger.info(f"hisat2_config: {hisat2_config}")

use rule hisat2_align from hisat2 as CoCulture_hisat2_align
use rule hisat2_index from hisat2 as CoCulture_hisat2_index

disambiguate_config = {

    "resources": config.get("resources", {}),
        "indir": hisat2_config["outdir"],
        "outdir": f"{outdir}/disambiguate",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "genome_pairs": genome_pairs,
        "single_samples": single_samples,
        "paired_samples": paired_samples,
        "bam": {
            "aligner": config.get('Procedure',{}).get('aligner') or 'hisat2'
        },
        "Params": {
            "disambiguate": {
                "keep_intermediates": config.get('Params',{}).get('disambiguate',{}).get('keep_intermediates', True)
            }
        },
        "Procedure": {
            "ngs_disambiguate": config.get('Procedure',{}).get('ngs_disambiguate') or 'ngs_disambiguate'
        }
    }
module disambiguate:
    snakefile: "../modules/disambiguate/disambiguate.smk"
    config: disambiguate_config
logger.info(f"disambiguate config: {disambiguate_config}")
use rule ngs_disambiguate from disambiguate as CoCulture_ngs_disambiguate
use rule disambiguate_sort_rename from disambiguate as CoCulture_disambiguate_sort_rename
use rule disambiguate_report from disambiguate as CoCulture_disambiguate_report

# ============================================
# Sex check (coverage on sex chromosomes)
# ============================================
# Wired unconditionally over the alignment BAMs; node.py registers the
# sexcheck outputs (subworkflow.md QC wiring).
# genome/sample tuples ("GRCh38/TestSample1") expanded from the per-genome dicts
_sexcheck_samples = [
    f"{g}/{s}"
    for g in list(genome_paired_samples) + list(genome_single_samples)
    for s in genome_paired_samples.get(g, []) + genome_single_samples.get(g, [])
]
sexcheck_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/sexcheck",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": _sexcheck_samples,
    "sample_bams": {
        s: f"{hisat2_config['outdir']}/{s.split('/')[0]}/{s.split('/')[1]}/{s.split('/')[1]}.bam"
        for s in _sexcheck_samples
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
        "sexcheck": config.get("Params", {}).get("sexcheck", {})
    },
}
module sexcheck:
    snakefile: "../modules/sexcheck/sexcheck.smk"
    config: sexcheck_config
logger.info(f"sexcheck_config: {sexcheck_config}")
use rule sexcheck_idxstats from sexcheck as CoCulture_sexcheck_idxstats
use rule sexcheck from sexcheck as CoCulture_sexcheck


# ============================================
# QC Summary (cross-sample metric aggregation)
# ============================================
# Wired unconditionally; node.py registers the qc_summary outputs (see
# subworkflow.md QC wiring spec). Only stages with qc_summary parsers are
# declared; unregistered stages stay best-effort scans.
qc_summary_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
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
logger.info(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as CoCulture_qc_summary



TEtranscripts_config = {



    "resources": config.get("resources", {}),
        "indir": disambiguate_config["outdir"],
        "outdir": f"{outdir}/TEtranscripts",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "genome_pairs": disambiguate_config["genome_pairs"],
        "samples": single_samples + paired_samples,
        # disambiguate TEtranscripts resolves refs via genome.<genome>.<field>
        # (flat map), while the workflow config keeps them under
        # genome.references.<genome>
        "genome": genome.get("references", {}),
        "Procedure": {
            "TEcount": config.get('Procedure',{}).get('TEcount') or 'TEcount',
            "TElocal": config.get('Procedure',{}).get('TElocal') or 'TElocal'
        },
        "Params": {
            "TEtranscripts": {
                "keep_intermediates": config.get('Params',{}).get('TEtranscripts',{}).get('keep_intermediates', True)
            }
        }
    }
module TEtranscripts:
    snakefile: "../modules/TEtranscripts/disambiguate/TEtranscripts.smk"
    config: TEtranscripts_config
logger.info(f"TEtranscripts config: {TEtranscripts_config}")
use rule * from TEtranscripts as CoCulture_*

shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
aligner_TEtranscripts = config.get('Params',{}).get('workflow', {}).get('aligner_TEtranscripts') or "star"
trimmer = config.get('Params',{}).get('workflow', {}).get('trimmer') or "cutadapt"
outfiles = config.get("outfiles", [])
rule all:
    input:
        outfiles
if trimmer == "cutadapt":
    cutadapt_config = {
            "indir": indir,
            "outdir":  f"{outdir}/fastq/cutadapt",
            "logdir": logdir,
            "Procedure": {
                "trim_galore": config.get('Procedure',{}).get('trim_galore')
            }
        }
    module cutadapt:
        snakefile: "../modules/cutadapt/cutadapt.smk"
        config: cutadapt_config
    logger.info(f"cutadapt_config: {cutadapt_config}")
    use rule trimming_Paired from cutadapt as RNAseq_trimming_Paireds
    use rule trimming_Single from cutadapt as RNAseq_trimming_Single
elif trimmer == "trimmomatic":
    trimmomatic_config = {
            "indir": indir,
            "outdir":  f"{outdir}/fastq/trimmomatic",
            "logdir": logdir,
            "Procedure": {
                "trimmomatic": config.get('Procedure',{}).get('trimmomatic')
            },
            "Params": {
                "trimmomatic": {
                    "adapter_pe": config.get('Params',{}).get("trimmomatic", {}).get('adapter_pe'),
                    "adapter_se": config.get('Params',{}).get("trimmomatic", {}).get('adapter_se')
                }
            }
        }
    module trimmomatic:
        snakefile: "../modules/trimmomatic/trimmomatic.smk"
        config: trimmomatic_config
    logger.info(f"trimmomatic_config: {trimmomatic_config}")
    use rule trimmomatic_Paired from trimmomatic as RNAseq_trimmomatic_Paireds
    use rule trimmomatic_Single from trimmomatic as RNAseq_trimmomatic_Singles
else:
    raise ValueError(f"Unsupported trimmer: {trimmer}")

if aligner_TEtranscripts == 'hisat2':
    hisat2_config_for_TEtranscripts = {
            "indir": cutadapt_config["outdir"] if trimmer == "cutadapt" else trimmomatic_config["outdir"],
            "outdir":  f"{outdir}/TEtranscripts/bam",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "hisat2": config.get('Procedure',{}).get('hisat2')
            },
            "Params": {
                "hisat2": {
                    "score_min": config.get('Params',{}).get('hisat2_TEtranscripts', {}).get('score_min') or "L,0,-0.2",
                    "flag_params": config.get('Params',{}).get('hisat2_TEtranscripts', {}).get('flag_params') or "--no-mixed --no-discordant",
                    "k": config.get('Params',{}).get('hisat2_TEtranscripts', {}).get('k') or 100
                }
            },
            "genome": {
                "fasta": config.get('genome',{}).get('fasta'),
                "index_prefix": config.get('genome',{}).get('hisat2_index_prefix')
            }
        }
    module hisat2_for_TEtranscripts:
        snakefile: "../modules/hisat2/hisat2.smk"
        config: hisat2_config_for_TEtranscripts
    logger.info(f"hisat2_config: {hisat2_config_for_TEtranscripts}")
    use rule hisat2_align from hisat2_for_TEtranscripts as RNAseq_hisat2_align_for_TEtranscripts
    use rule hisat2_index from hisat2_for_TEtranscripts as RNAseq_hisat2_index_for_TEtranscripts
elif aligner_TEtranscripts == 'star':
    star_config_for_TEtranscripts = {
            "indir": cutadapt_config["outdir"] if trimmer == "cutadapt" else trimmomatic_config["outdir"],
            "outdir":  f"{outdir}/TEtranscripts/bam",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "star": config.get('Procedure',{}).get('star')
            },
            "Params": {
                "STAR": {
                    "alignEndsType": config.get('Params',{}).get('STAR_TEtranscripts', {}).get('alignEndsType') or "Local",
                    "outFilterMultimapNmax": config.get('Params',{}).get('STAR_TEtranscripts', {}).get('outFilterMultimapNmax') or 100,
                    "winAnchorMultimapNmax": config.get('Params',{}).get('STAR_TEtranscripts', {}).get('winAnchorMultimapNmax') or 100
                }
            },
            "genome": {
                "fasta": config.get('genome',{}).get('fasta'),
                "gtf": config.get('genome',{}).get('gtf'),
                "index_dir": config.get('genome',{}).get('star_index_dir')
            }
        }
    module star_for_TEtranscripts:
        snakefile: "../modules/star/star.smk"
        config: star_config_for_TEtranscripts
    logger.info(f"star_config: {star_config_for_TEtranscripts}")
    use rule star_align from star_for_TEtranscripts as RNAseq_star_align_for_TEtranscripts
    use rule star_index from star_for_TEtranscripts as RNAseq_star_index_for_TEtranscripts
else:
    raise ValueError(f"Unsupported aligner_TEtranscripts: {aligner_TEtranscripts}")


TEtranscripts_config = {
        "indir": star_config_for_TEtranscripts["outdir"] if aligner_TEtranscripts == 'star' else hisat2_config_for_TEtranscripts["outdir"],
        "outdir":  f"{outdir}/TEtranscripts",
        "logdir": logdir,
        "samples": single_samples + paired_samples,
        "ROOT_DIR": ROOT_DIR,
        "Procedure": {
            "TEcount": config.get('Procedure',{}).get('TEcount') or 'TEcount',
            "TElocal": config.get('Procedure',{}).get('TElocal') or 'TElocal'
        },
        "genome": {
            "gtf": config.get('genome',{}).get('gtf'),
            "TEind": config.get('genome',{}).get('TEind'),
            "TE_gtf": config.get('genome',{}).get('TE_gtf')
        }
    }

module TEtranscripts:
    snakefile: "../modules/TEtranscripts/TEtranscripts.smk"
    config: TEtranscripts_config
logger.info(f"TEtranscripts_config: {TEtranscripts_config}")
use rule * from TEtranscripts as RNAseq_*

DESeq2_config = {
        "indir": TEtranscripts_config["outdir"],
        "outdir":  f"{outdir}/DESeq2",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
        "control_samples": config.get("control_samples", []),
        "control_group_name": config.get("control_group_name", "control"),
        "treatment_samples": config.get("treatment_samples", []),
        "experimental_group_name": config.get("experimental_group_name", "treatment"),
        "genome": {
            "geneIDAnno": config.get('genome',{}).get('geneIDAnno')
        },
        "Procedure": {
            "DESeq2": config.get('Procedure',{}).get('DESeq2') or 'DESeq2'
        }
    }
module DESeq2:
    snakefile: "../modules/DESeq2/DESeq2.smk"
    config: DESeq2_config
logger.info(f"DESeq2_config: {DESeq2_config}")
use rule DESeq2_TEcount from DESeq2 as RNAseq_DESeq2_TEcount


hisat2_config_for_StringTie = {
    "indir": cutadapt_config["outdir"] if trimmer == "cutadapt" else trimmomatic_config["outdir"],
    "outdir":  f"{outdir}/stringtie/bam",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "hisat2": config.get('Procedure',{}).get('hisat2')
    },
    "Params": {
        "hisat2": {
            "k": config.get('Params',{}).get('hisat2_stringtie', {}).get('k') or 1,
            "flag_params": config.get('Params',{}).get('hisat2_stringtie', {}).get('flag_params') or "-q --no-unal --dta",
        }
    },
    "genome": {
        "fasta": config.get('genome',{}).get('fasta'),
        "index_prefix": config.get('genome',{}).get('hisat2_index_prefix')
    }
}
module hisat2_for_StringTie:
    snakefile: "../modules/hisat2/hisat2.smk"
    config: hisat2_config_for_StringTie
logger.info(f"hisat2_config_for_StringTie: {hisat2_config_for_StringTie}")
use rule hisat2_align from hisat2_for_StringTie as RNAseq_hisat2_align_for_StringTie
use rule hisat2_index from hisat2_for_StringTie as RNAseq_hisat2_index_for_StringTie
 
StringTie_config = {
        "indir": hisat2_config_for_StringTie["outdir"],
        "outdir":  f"{outdir}/stringtie",
        "logdir": logdir,
        "samples": single_samples + paired_samples,
        "ROOT_DIR": ROOT_DIR,
        "sample_groups": config.get('sample_groups'),
        "genome": {
            "gtf": config.get('genome',{}).get('gtf'),
            "TE_gtf": config.get('genome',{}).get('TE_gtf')
        },
        "Procedure": {
            "stringtie": config.get("Procedure", {}).get("stringtie") or "stringtie"
        }
    }
module StringTie:
    snakefile: "../modules/StringTie/StringTie.smk"
    config: StringTie_config
use rule * from StringTie as RNAseq_*
logger.info(f"StringTie_config: {StringTie_config}")


rmrRNA_config = {
    "indir": cutadapt_config["outdir"] if trimmer == "cutadapt" else trimmomatic_config["outdir"],
    "outdir":  f"{outdir}/rRNA",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "ROOT_DIR": ROOT_DIR,
    "genome": {
        "fasta": config.get('genome',{}).get('fasta'),
        "gtf": config.get('genome',{}).get('gtf')
    },
    "Params": {
        "RmrRNA": {
            "sam-append-comment": config.get('Params',{}).get('RmrRNA_rRNA', {}).get('sam-append-comment')
        }
    }
}

module RmrRNA:
    snakefile: "../modules/RmrRNA/RmrRNA.smk"
    config: rmrRNA_config
use rule * from RmrRNA as RNAseq_rRNA_*

logger.info(f"rmrRNA_config: {rmrRNA_config}")
bowtie2_rRNA_config = {
    "indir": cutadapt_config["outdir"] if trimmer == "cutadapt" else trimmomatic_config["outdir"],
    "outdir":  f"{outdir}/rRNA",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "bowtie2": config.get('Procedure',{}).get('bowtie2'),
        "bowtie2-build": config.get('Procedure',{}).get('bowtie2-build')
    },
    "genome": {
        "fasta": config.get('genome',{}).get('rRNA_fasta') if config.get('genome',{}).get('rRNA_fasta') else f"{rmrRNA_config['outdir']}/rRNA.fasta",
        "index_prefix": config.get('genome',{}).get('bowtie2_index_prefix_for_rRNA') if config.get('genome',{}).get('rRNA_fasta') else None
    },
    "Params": {
        "bowtie2": {
            "sam-append-comment": config.get('Params',{}).get('bowtie2_rRNA', {}).get('sam-append-comment')
        }
    }
}
module bowtie2_for_rRNA:
    snakefile: "../modules/bowtie2/bowtie2.smk"
    config: bowtie2_rRNA_config
use rule * from bowtie2_for_rRNA as RNAseq_rRNA_*

star_config_for_fusion = {
    "indir": bowtie2_rRNA_config["outdir"],
    "outdir":  f"{outdir}/fusion/star",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "fastq_sample_suffix": "unmapped",
    "Procedure": {
        "STAR": config.get('Procedure',{}).get('STAR')
    },
    "Params": {
        "STAR": {
            "outFilterMultimapNmax": config.get('Params',{}).get('STAR_fusion', {}).get('outFilterMultimapNmax') or 1,
            "outFilterMismatchNmax": config.get('Params',{}).get('STAR_fusion', {}).get('outFilterMismatchNmax') or 3,
            "chimSegmentMin": config.get('Params',{}).get('STAR_fusion', {}).get('chimSegmentMin') or 10,
            "chimOutType": config.get('Params',{}).get('STAR_fusion', {}).get('chimOutType') or "Junctions WithinBAM SoftClip SeparateSAMold",
            "chimJunctionOverhangMin": config.get('Params',{}).get('STAR_fusion', {}).get('chimJunctionOverhangMin') or 10,
            "outSAMstrandField": config.get('Params',{}).get('STAR_fusion', {}).get('outSAMstrandField') or "intronMotif",
            "chimScoreMin": config.get('Params',{}).get('STAR_fusion', {}).get('chimScoreMin') or 1,
            "chimScoreDropMax": config.get('Params',{}).get('STAR_fusion', {}).get('chimScoreDropMax') or 30,
            "chimScoreJunctionNonGTAG": config.get('Params',{}).get('STAR_fusion', {}).get('chimScoreJunctionNonGTAG') or 0,
            "chimScoreSeparation": config.get('Params',{}).get('STAR_fusion', {}).get('chimScoreSeparation') or 1,
            "alignSJstitchMismatchNmax": config.get('Params',{}).get('STAR_fusion', {}).get('alignSJstitchMismatchNmax') or "5 -1 5 5",
            "chimSegmentReadGapMax": config.get('Params',{}).get('STAR_fusion', {}).get('chimSegmentReadGapMax') or 3,
        }
    },
    "genome": {
        "fasta": config.get('genome',{}).get('fasta'),
        "gtf": config.get('genome',{}).get('gtf'),
        "index_dir": config.get('genome',{}).get('star_index_dir')
    }
}
module star_for_fusion:
    snakefile: "../modules/star/star.smk"
    config: star_config_for_fusion
use rule * from star_for_fusion as RNAseq_fusion_*
logger.info(f"star_config_for_fusion: {star_config_for_fusion}")

gatk_prepare_config = {
    "indir": star_config_for_fusion["outdir"],
    "outdir":  f"{outdir}/fusion/bam",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "fasta": config.get('genome',{}).get('fasta'),
        "gtf": config.get('genome',{}).get('gtf')
    }
}

module gatk_prepare:
    snakefile: "../modules/gatk/gatk_prepare.smk"
    config: gatk_prepare_config
use rule * from gatk_prepare as RNAseq_fusion_*
logger.info(f"gatk_prepare_config: {gatk_prepare_config}")

arriba_config = {
        "indir": f"{gatk_prepare_config['outdir']}/bam-sorted-Markdup",
        "outdir":  f"{outdir}/fusion/arriba",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
        "samples": single_samples + paired_samples,
        "genome": {
            "fasta": config.get('genome',{}).get('fasta'),
            "gtf": config.get('genome',{}).get('gtf')
        },
        "Params": {
            "arriba": {
                "blacklist": config.get('Params',{}).get('arriba',{}).get('blacklist'),
                "known_fusions": config.get('Params',{}).get('arriba',{}).get('known_fusions'),
                "t": config.get('Params',{}).get('arriba',{}).get('t'),
                "d": config.get('Params',{}).get('arriba',{}).get('d'),
                "E": config.get('Params',{}).get('arriba',{}).get('E'),
                "p": config.get('Params',{}).get('arriba',{}).get('p')
            }
        },
        "Procedure": {
            "arriba": config.get('Procedure',{}).get('arriba') or 'arriba'
        }
    }
module arriba:
    snakefile: "../modules/arriba/arriba.smk"
    config: arriba_config
use rule * from arriba as RNAseq_fusion_*
logger.info(f"arriba_config: {arriba_config}")

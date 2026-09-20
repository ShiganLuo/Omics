shell.prefix("set -x; set -e;")
from snakemake.logging import logger

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
outfiles = config.get("outfiles", [])
samples = config.get("samples", [])
skip_quant = config.get("Params", {}).get("skip_quant", False)

genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(genome, {})

rule all:
    input:
        outfiles

# ============================================================
# Step 1: Align with minimap2
# ============================================================
minimap2_config = {
    "indir": indir,
    "outdir": f"{outdir}/alignment",
    "logdir": logdir,
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "Procedure": {
        "minimap2": config.get("Procedure", {}).get("minimap2"),
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "Params": {
        "minimap2": config.get("Params", {}).get("minimap2", {})
    },
    "genome": {
        "fasta": genome_ref.get("fasta")
    }
}
module minimap2:
    snakefile: "../modules/minimap2/minimap2.smk"
    config: minimap2_config
logger.info(f"minimap2_config: {minimap2_config}")
use rule minimap2_align from minimap2 as LRtranscriptome_minimap2_align

# ============================================================
# Step 2-4: StringTie long-read pipeline (assemble + merge + quant)
# ============================================================
longread_config = {
    "indir": minimap2_config["outdir"],
    "upstream_indir": minimap2_config["outdir"],
    "outdir": f"{outdir}/stringtie",
    "logdir": logdir,
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "merged_gtf": f"{outdir}/stringtie/stringtie_merged.gtf",
    "Procedure": {
        "stringtie": config.get("Procedure", {}).get("stringtie") or "stringtie"
    },
    "genome": {
        "gtf": genome_ref.get("gtf")
    }
}
module longread:
    snakefile: "../modules/StringTie/longread/longread.smk"
    config: longread_config
logger.info(f"longread_config: {longread_config}")
use rule stringtie_lr_assemble from longread as LRtranscriptome_stringtie_assemble
use rule stringtie_lr_merge from longread as LRtranscriptome_stringtie_merge
use rule stringtie_lr_result from longread as LRtranscriptome_stringtie_result

if not skip_quant:
    longread_quant_config = {
        "indir": minimap2_config["outdir"],
        "upstream_indir": minimap2_config["outdir"],
        "outdir": f"{outdir}/quantification",
        "logdir": logdir,
        "samples": samples,
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "merged_gtf": f"{outdir}/stringtie/stringtie_merged.gtf",
        "Procedure": {
            "stringtie": config.get("Procedure", {}).get("stringtie") or "stringtie"
        },
        "genome": {
            "gtf": genome_ref.get("gtf")
        }
    }
    module longread_quant:
        snakefile: "../modules/StringTie/longread/longread.smk"
        config: longread_quant_config
    logger.info(f"longread_quant_config: {longread_quant_config}")
    use rule stringtie_lr_quant from longread_quant as LRtranscriptome_stringtie_quant
    use rule stringtie_lr_result from longread_quant as LRtranscriptome_quant_result
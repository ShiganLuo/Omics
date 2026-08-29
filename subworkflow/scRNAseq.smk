from snakemake.logging import logger
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
outfiles = config.get("outfiles") or []
paired_samples = config.get("paired_samples") or []
single_samples = config.get("single_samples") or []
aligner = config.get("aligner") or "cellranger"
counters = config.get("counters") or ["scTE", "cellranger"]
rule all:
    input:
        outfiles
genome = config.get("genome", {}).get("default")
raw_bam_outdir = f"{outdir}/common/2_raw_bam"
h5ad_outdir = f"{outdir}/common/3_raw_h5ad"
if aligner == "star":
    logger.info("Using STAR aligner for scRNA-seq workflow. only scTE counter is supported for STAR aligner.")
    star_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir": raw_bam_outdir,
        "logdir": f"{logdir}/sample",
        "logdir_index": f"{logdir}/group",
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "omics_type": "scRNAseq",
        "Params": {
            "star": {
                "outSAMattributes": config.get('Params', {}).get('star', {}).get('outSAMattributes'),
                "outFilterMultimapNmax": config.get('Params', {}).get('star', {}).get('outFilterMultimapNmax'),
                "winAnchorMultimapNmax": config.get('Params', {}).get('star', {}).get('winAnchorMultimapNmax'),
                "outMultimapperOrder": config.get('Params', {}).get('star', {}).get('outMultimapperOrder'),
                "runRNGseed": config.get('Params', {}).get('star', {}).get('runRNGseed'),
                "outSAMmultNmax": config.get('Params', {}).get('star', {}).get('outSAMmultNmax'),
                "soloType": config.get('Params', {}).get('star', {}).get('soloType'),
                "soloCBwhitelist": config.get('Params', {}).get('star', {}).get('soloCBwhitelist'),
                "soloBarcodeReadLength": config.get('Params', {}).get('star', {}).get('soloBarcodeReadLength'),
                "limitSjdbInsertNsj": config.get('Params', {}).get('star', {}).get('limitSjdbInsertNsj')
            }
        },
        "genome": {
            "index_dir": config.get("genome", {}).get("references",{}).get(genome, {}).get("index_dir"),
            "gtf": config.get("genome", {}).get("references",{}).get(genome, {}).get("gtf"),
            "fasta": config.get("genome",{}).get("references",{}).get(genome, {}).get("fasta")
        }
    }
    module star:
        snakefile: "../modules/star/star.smk"
        config: star_config
    logger.info(f"Using STAR aligner with scTE counter for scRNA-seq workflow. STAR config: {star_config}")
    use rule * from star as scRNAseq_*
    scTE_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": star_config["outdir"],
        "outdir": h5ad_outdir,
        "logdir": f"{logdir}/sample",
        "logdir_index": f"{logdir}/group",
        "Params": {
            "scTE": config.get("Params", {}).get("scTE", {}),
        },
        "genome": {
            "scTE_index": config.get("genome", {}).get("references", {}).get(genome, {}).get("scTE_index"),
            "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf"),
            "te_bed": config.get("genome", {}).get("references", {}).get(genome, {}).get("te_bed")
        },
        "Procedure": {
            "scTE": config.get("Procedure", {}).get("scTE") or "scTE"
        }
    }
    module scTE:
        snakefile: "../modules/scTE/scTE.smk"
        config: scTE_config
    logger.info(f"Using scTE counter for scRNA-seq workflow. scTE config: {scTE_config}")
    use rule * from scTE as scRNAseq_*

elif aligner == "cellranger":
    cellranger_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir": raw_bam_outdir,
        "logdir": f"{logdir}/sample",
        "logdir_ref": f"{logdir}/group",
        "h5ad_outdir": h5ad_outdir,
        "cellranger_input_dict": config.get("cellranger_input_dict", {}),
        "Params": {
            "cellranger": config.get("Params", {}).get("cellranger", {}),
        },
        "genome": {
            "cellranger_ref_dir": config.get("genome", {}).get("references",{}).get(genome,{}).get("cellranger_ref_dir"),
            "fasta": config.get("genome", {}).get("references",{}).get(genome,{}).get("fasta"),
            "gtf": config.get("genome", {}).get("references",{}).get(genome,{}).get("gtf")
        },
        "Procedure": {
            "cellranger": config.get("Procedure", {}).get("cellranger") or "cellranger"
        }
    }
    module cellranger:
        snakefile: "../modules/cellranger/cellranger.smk"
        config: cellranger_config
    logger.info(f"Using Cell Ranger aligner for scRNA-seq workflow. Cell Ranger config: {cellranger_config}")
    use rule * from cellranger as scRNAseq_*

    if "scTE" in counters:
        scTE_config = {
            "ROOT_DIR": ROOT_DIR,
            "env": config.get("env", {}),
            "indir": cellranger_config["outdir"],
            "outdir": h5ad_outdir,
            "logdir": f"{logdir}/sample",
            "logdir_index": f"{logdir}/group",
            "Params": {
                "scTE": config.get("Params", {}).get("scTE", {}),
            },
            "genome": {
                "scTE_index": config.get("genome", {}).get("references", {}).get(genome, {}).get("scTE_index"),
                "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf"),
                "te_bed": config.get("genome", {}).get("references", {}).get(genome, {}).get("te_bed")
            },
            "Procedure": {
                "scTE": config.get("Procedure", {}).get("scTE") or "scTE"
            }
        }
        module scTE:
            snakefile: "../modules/scTE/scTE.smk"
            config: scTE_config
        logger.info(f"Using scTE counter for scRNA-seq workflow. scTE config: {scTE_config}")
        use rule * from scTE as scRNAseq_*
    elif "cellranger" in counters:
        logger.info("Using Cell Ranger counter for scRNA-seq workflow. No additional configuration needed.")
    else:
        raise ValueError(f"Unsupported counter: {counters}. Please choose either 'scTE' or 'cellranger' for counter.")

else:
    raise ValueError(f"Unsupported aligner or counter: {aligner}, {counters}. Please choose either 'star' or 'cellranger' for aligner and 'scTE' or 'cellranger' for counter.")


scanpy_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": h5ad_outdir,
    "outdir": f"{outdir}/common/4_qc_h5ad",
    "outdir_combine": f"{outdir}/common/5_combine_h5ad",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "Params": {
        "scanpy": config.get("Params", {}).get("scanpy", {}),
    },
    "Procedure": {
        "python": config.get("Procedure", {}).get("python") or "python",
    },
}
module scanpy:
    snakefile: "../modules/scanpy/scanpy.smk"
    config: scanpy_config
logger.info(f"Using scanpy for scRNA-seq downstream analysis. scanpy config: {scanpy_config}")
use rule * from scanpy as scRNAseq_*
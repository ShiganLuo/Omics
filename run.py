import argparse
from copy import deepcopy
import json
import os
from src.common.MetaUtil import MetadataUtils
from src.common.LogUtil import setup_logger
from src.common.CmdUtil import _run_cmd
import logging
from typing import Dict, Any, Optional
logger = logging.getLogger(__name__)

def _load_model_json(model_json_file: str) -> Dict[str, Any]:
    """Load model JSON template from disk."""
    with open(model_json_file, 'r', encoding='utf-8') as f:
        return json.load(f)
    
def runCoCulture(
    input_json: str,
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,
):
    datajsonTemplate = _load_model_json(input_json)
    datajson = deepcopy(datajsonTemplate)
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    outfiles = []
    paired_samples = []
    single_samples = []
    single_sample_genome_pairs = []
    paired_sample_genome_pairs = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            paired_sample_genome_pairs.append((sample_id, sample_info.organism))
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_2.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            single_sample_genome_pairs.append((sample_id, sample_info.organism))
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}.single.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")

    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["single_sample_genome_pairs"] = single_sample_genome_pairs
    datajson["paired_sample_genome_pairs"] = paired_sample_genome_pairs
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runMERIP(
    input_json: str,
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,
):
    """
    Function: Prepare input JSON for MERIP workflow based on the provided model JSON template and sample information.
    Parameters:
    - input_json: Path to the model JSON template file.
    - samples_info_dict: A dictionary containing sample information, where keys are sample IDs and values are objects with attributes 'layout' and 'design'.
    - indir: Input directory containing raw data (e.g., FASTQ files).
    - outdir: Output directory where results will be stored.
    Returns:
    - instance_json: Path to the generated input JSON file that will be used for the MERIP workflow.
    """
    datajsonTemplate = _load_model_json(input_json)
    datajson = deepcopy(datajsonTemplate)
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir

    paired_samples = []
    single_samples = []
    ip_samples = []
    input_samples = []
    treated_ip_samples = []
    treated_input_samples = []
    outfiles = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}_1.fq.gz")
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}_2.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/dedup/{sample_id}.dedup.bam")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}.single.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/dedup/{sample_id}.dedup.bam")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
        
        if sample_info.design == "ip":
            ip_samples.append(sample_id)
        elif sample_info.design == "input":
            input_samples.append(sample_id)
        elif sample_info.design == "treated_ip":
            treated_ip_samples.append(sample_id)
        elif sample_info.design == "treated_input":
            treated_input_samples.append(sample_id)
        else:
            logger.error(f"Unknown design type for sample {sample_id}: {sample_info.design}")
    outfiles.append(f"{outdir}/exomePeak/sig_diff_peak_gene_names.xls")
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["ip_samples"] = ip_samples
    datajson["input_samples"] = input_samples
    datajson["treated_ip_samples"] = treated_ip_samples
    datajson["treated_input_samples"] = treated_input_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runRNAseq(
    input_json: str,
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,
):
    datajsonTemplate = _load_model_json(input_json)
    datajson = deepcopy(datajsonTemplate)
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    outfiles = []
    paired_samples = []
    single_samples = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    outfiles.append(f"{outdir}/TEtranscripts/TEcount/all_TEcount.tsv")
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def parse_args():
    parser = argparse.ArgumentParser(description="workflow")
    parser.add_argument('-m','--meta', type=str, required=True, help='meta input file or data dir which condatain fastq file')
    parser.add_argument('-w','--workflow_name', type=str, choices=["CoCulture", "MERIP", "RNAseq"],default='CoCulture' ,help='workflow name')
    parser.add_argument('-o','--output_dir', type=str, required=True, help='output dir')
    parser.add_argument('-t','--threads', type=int, default=10, help='threads')
    parser.add_argument('--dry-run', action='store_true', help='dry run')
    parser.add_argument('--log', type=str, default='workflow.log', help='log file')
    parser.add_argument('--conda-prefix', type=str, default='/data/pub/zhousha/env/mutation_0.1', help='conda prefix for snakemake')
    parser.add_argument('--rerun-trigger', type=str, default="input", choices=["code", "input", "mtime", "params", "software-env"],help='snakemake rerun-triggers, e.g.  code, input, mtime, params, software-env')
    args = parser.parse_args()
    return args
if __name__ == "__main__":
    args = parse_args()
    logger = setup_logger("root",args.log)
    ROOT_DIR = os.path.dirname(__file__)
    outdir = os.path.join(args.output_dir, args.workflow_name)
    abs_outdir = os.path.abspath(outdir)
    if os.path.isfile(args.meta):
        metadataUtil = MetadataUtils(
            outdir = abs_outdir,
            meta = args.meta,
        )
    else:
        metadataUtil = MetadataUtils(
            outdir = abs_outdir,
            fastq_dir = args.meta,
        )
    samples_info_dict, pairs, raw_fastq_dir = metadataUtil.run()
    os.makedirs(abs_outdir, exist_ok=True)
    if args.workflow_name == "CoCulture":
        model_json = os.path.join(ROOT_DIR, "config/CoCulture.json")
        input_json = runCoCulture(model_json, samples_info_dict, raw_fastq_dir, abs_outdir)
        if args.dry_run:
            logger.info(f"Dry run mode, generated input json: {input_json}")
            cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/CoCulture.smk", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda","--dry-run"]
        else:
            cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/CoCulture.smk", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda"]
        _run_cmd(cmds)
    elif args.workflow_name == "MERIP":
        model_json = os.path.join(ROOT_DIR, "config/MERIP.json")
        input_json = runMERIP(model_json, samples_info_dict, raw_fastq_dir, abs_outdir)
        if args.dry_run:
            logger.info(f"Dry run mode, generated input json: {input_json}")
            cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/MERIP.smk", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda","--dry-run"]
        else:
            cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/MERIP.smk", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda"]
        _run_cmd(cmds)
    elif args.workflow_name == "RNAseq":
        model_json = os.path.join(ROOT_DIR, "config/RNAseq.json")
        input_json = runRNAseq(model_json, samples_info_dict, raw_fastq_dir, abs_outdir)
        if args.dry_run:
            logger.info(f"Dry run mode, generated input json: {input_json}")
            cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/RNAseq.smk", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda","--dry-run"]
        else:
            cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/RNAseq.smk", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda"]
        _run_cmd(cmds)
    else:
        logger.error(f"Unknown workflow name: {args.workflow_name}")
    pass
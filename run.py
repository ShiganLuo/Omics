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
    **kwargs
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
    **kwargs
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
    **kwargs
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

def runCLIP(
    input_json: str,
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,
    **kwargs
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
            outfiles.append(f"{outdir}/cutadapt/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/cutadapt/{sample_id}_2.fq.gz")
            outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/cutadapt/{sample_id}.single.fq.gz")
            outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json
    pass
def parse_args():
    parser = argparse.ArgumentParser(description="workflow")
    parser.add_argument('-m','--meta', type=str, required=True, help='meta input file or data dir which condatain fastq file')
    parser.add_argument('-w','--workflow_name', type=str, choices=["CoCulture", "MERIP", "RNAseq", "CLIP"],default='CoCulture' ,help='workflow name')
    parser.add_argument('-o','--output_dir', type=str, required=True, help='output dir')
    parser.add_argument('-t','--threads', type=int, default=10, help='threads')
    parser.add_argument('--dry-run', action='store_true', help='dry run')
    parser.add_argument('--log', type=str, default='workflow.log', help='log file')
    parser.add_argument('--conda-prefix', type=str, default='/data/pub/zhousha/env/mutation_0.1', help='conda prefix for snakemake')
    parser.add_argument('--rerun-trigger', type=str, default="input", choices=["code", "input", "mtime", "params", "software-env"],help='snakemake rerun-triggers, e.g.  code, input, mtime, params, software-env')
    parser.add_argument('--conda-frontend', type=str, choices=["conda", "mamba"], default="mamba", help='conda frontend for snakemake')
    parser.add_argument('--config', type=str, help='extra config json file to override workflow defaults')
    # 支持 --key=value 和 --key value 两种形式的额外参数
    args, unknown = parser.parse_known_args()
    extra_args = {}
    i = 0
    while i < len(unknown):
        arg = unknown[i]
        if arg.startswith('--'):
            key = arg[2:]
            if '=' in key:
                k, v = key.split('=', 1)
                extra_args[k] = v
            elif i + 1 < len(unknown) and not unknown[i + 1].startswith('--'):
                extra_args[key] = unknown[i + 1]
                i += 1
            else:
                extra_args[key] = True
        i += 1
    args.extra_args = extra_args
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
    
    # 根据 workflow_name 选择 json 模板
    model_json = os.path.join(ROOT_DIR, f"config/{args.workflow_name}.json")
    # 读取默认参数
    with open(model_json, 'r', encoding='utf-8') as f:
        workflow_config = json.load(f)
    # 如果有 --config 指定的 json，合并覆盖
    if args.config:
        with open(args.config, 'r', encoding='utf-8') as f:
            user_config = json.load(f)
        workflow_config.update(user_config)
    # 命令行额外参数覆盖
    workflow_config.update(args.extra_args)
    # 移除关键路径参数，防止与函数签名冲突
    for k in ["indir", "outdir", "logdir"]:
        workflow_config.pop(k, None)
    # 传递给下游函数
    if args.workflow_name == "CoCulture":
        input_json = runCoCulture(model_json, samples_info_dict, raw_fastq_dir, abs_outdir, **workflow_config)
        smk = "CoCulture.smk"
    elif args.workflow_name == "MERIP":
        input_json = runMERIP(model_json, samples_info_dict, raw_fastq_dir, abs_outdir, **workflow_config)
        smk = "MERIP.smk"
    elif args.workflow_name == "RNAseq":
        input_json = runRNAseq(model_json, samples_info_dict, raw_fastq_dir, abs_outdir, **workflow_config)
        smk = "RNAseq.smk"
    elif args.workflow_name == "CLIP":
        input_json = runCLIP(model_json, samples_info_dict, raw_fastq_dir, abs_outdir, **workflow_config)
        smk = "CLIP.smk"
    else:
        logger.error(f"Unknown workflow name: {args.workflow_name}")
        exit(1)
    if args.dry_run:
        logger.info(f"Dry run mode, generated input json: {input_json}")
        cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/{smk}", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda","--dry-run", args.conda_frontend]
    else:
        cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/{smk}", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda", args.conda_frontend]
    _run_cmd(cmds)
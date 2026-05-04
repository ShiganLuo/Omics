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

def smart_cast(val):
    """尝试将字符串转换为 int/float/bool，否则原样返回"""
    if isinstance(val, str):
        if val.lower() in {"true", "false"}:
            return val.lower() == "true"
        try:
            if val.startswith("0") and len(val) > 1 and not val.startswith("0."):
                return val  # 避免八进制等
            return int(val)
        except Exception:
            pass
        try:
            return float(val)
        except Exception:
            pass
    return val

def dict_set_by_path(d, keys, value):
    """递归设置嵌套字典的值，keys为key列表，自动类型转换"""
    for k in keys[:-1]:
        if k not in d or not isinstance(d[k], dict):
            d[k] = {}
        d = d[k]
    d[keys[-1]] = smart_cast(value)

def parse_dot_args(extra_args):
    """从extra_args中提取点号语法参数，返回{(k1,k2,...):v}"""
    dot_args = {}
    for k, v in list(extra_args.items()):
        if '.' in k:
            dot_args[tuple(k.split('.'))] = v
    return dot_args


def _load_model_json(model_json_file: str) -> Dict[str, Any]:
    """Load model JSON template from disk."""
    with open(model_json_file, 'r', encoding='utf-8') as f:
        return json.load(f)
    
def runCoCulture(
    datajson: Dict[str,Any],
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,

):
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
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_2.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount_name.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount_name.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}.single.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount_name.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount_name.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    outfiles.append(f"{outdir}/disambiguate/disambiguate_qc.tsv")
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runMERIP(
    datajson: Dict[str, Any],
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
    datajson: Dict[str, Any],
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

def runCLIP(
    datajson: Dict[str, Any],
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,
):
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
            if datajson["aligner"] == "star":
                outfiles.append(f"{outdir}/star/{sample_id}/{sample_id}.bam")
            elif datajson["aligner"] == "hisat2":
                outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            outfiles.append(f"{outdir}/fastqc/raw/{sample_id}/fastqc.raw.txt")
            outfiles.append(f"{outdir}/fastqc/trimmed/{sample_id}/fastqc.trimmed.txt")
            outfiles.append(f"{outdir}/PureCLIP/{sample_id}.pureclip.sites.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.plus.bw")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.minus.bw")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/cutadapt/{sample_id}.single.fq.gz")
            if datajson["aligner"] == "star":
                outfiles.append(f"{outdir}/star/{sample_id}/{sample_id}.bam")
            elif datajson["aligner"] == "hisat2":
                outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            outfiles.append(f"{outdir}/fastqc/raw/{sample_id}/fastqc.raw.txt")
            outfiles.append(f"{outdir}/fastqc/trimmed/{sample_id}/fastqc.trimmed.txt")
            outfiles.append(f"{outdir}/PureCLIP/{sample_id}.pureclip.sites.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.plus.bw")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.minus.bw")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    outfiles.append(f"{outdir}/track/igv_track_iclip.html")
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    # parameters suggest by https://doi.org/10.1016/j.ymeth.2019.11.008
    datajson["Params"]["bamCoverage"]["offset"] = "-1"
    datajson["Params"]["bamCoverage"]["binSize"] = 1
    datajson["Params"]["bamCoverage"]["normalizeUsing"] = "CPM"
    datajson["Params"]["bamCoverage"]["extendReads"] = 1
    datajson["Params"]["STAR"]["alignEndsType"] = "Extend5pOfRead1"
    datajson["Params"]["STAR"]["outFilterMismatchNoverReadLmax"] = 0.04
    datajson["Params"]["STAR"]["outFilterMismatchNmax"] = 999
    datajson["Params"]["STAR"]["outFilterMultimapNmax"] = 999
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

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
    workflow_config = _load_model_json(model_json)
    # 命令行额外参数覆盖（先处理非点号参数）
    flat_args = {k: v for k, v in args.extra_args.items() if '.' not in k}
    workflow_config.update(flat_args)
    # 递归合并点号语法参数
    dot_args = parse_dot_args(args.extra_args)
    for key_tuple, v in dot_args.items():
        logger.info(f"Setting config parameter {'.'.join(key_tuple)} to {v} from command line")
        dict_set_by_path(workflow_config, list(key_tuple), v)
    
    # 传递给 downstream 函数，直接传递 workflow_config（已合并参数）
    if args.workflow_name == "CoCulture":
        input_json = runCoCulture(deepcopy(workflow_config), samples_info_dict, raw_fastq_dir, abs_outdir)
        smk = "CoCulture.smk"
    elif args.workflow_name == "MERIP":
        input_json = runMERIP(deepcopy(workflow_config), samples_info_dict, raw_fastq_dir, abs_outdir)
        smk = "MERIP.smk"
    elif args.workflow_name == "RNAseq":
        input_json = runRNAseq(deepcopy(workflow_config), samples_info_dict, raw_fastq_dir, abs_outdir)
        smk = "RNAseq.smk"
    elif args.workflow_name == "CLIP":
        input_json = runCLIP(deepcopy(workflow_config), samples_info_dict, raw_fastq_dir, abs_outdir)
        smk = "CLIP.smk"
    else:
        logger.error(f"Unknown workflow name: {args.workflow_name}")
        exit(1)
    if args.dry_run:
        logger.info(f"Dry run mode, generated input json: {input_json}")
        cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/{smk}", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda","--dry-run", "--conda-frontend",args.conda_frontend]
    else:
        cmds = ["snakemake","-s", f"{ROOT_DIR}/subworkflow/{smk}", "--configfile", input_json, "--cores", str(args.threads), "--conda-prefix", args.conda_prefix, "--rerun-triggers", args.rerun_trigger, "--use-conda", "--conda-frontend",args.conda_frontend]
    _run_cmd(cmds)
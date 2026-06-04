import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.LogUtil import setup_logger
from common.CmdUtil import _run_cmd
import logging
import re
from pathlib import Path
from typing import Literal
from typing import Union, List
logger = setup_logger("SVCircosPipeline", level=logging.INFO)
ImageFormat = Literal["png", "pdf", "svg"]
def run_cricos_pipeline(
        vcf:str,
        fasta: str,
        outdir: str,
        image_formats: Union[ImageFormat, List[ImageFormat]] = "png",
        outfile_name_mode: Literal["parent","sample"] = "sample",  # parent / sample
        prepare_script: str = "sv_circos_prepare.py",
        circos_script: str = "circos.r",
        ins_bin_size: int = 100000,
        ins_plot_type: Literal["points","bar"] = "bar",
        genome: str = "mm39",
        cytoband_file: str = "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/mm39.cytoBand.txt",
):
    """
    参数化的 SV Circos 绘图流程
    :param indir: 预处理结果输入目录
    :param fasta: 参考基因组 FASTA 文件路径
    :param outdir: Circos 绘图结果输出目录
    :param prepare_script: SV Circos 准备脚本路径
    :param circos_script: SV Circos 绘图脚本路径
    """
    os.makedirs(outdir, exist_ok=True)
    prepare_script = os.path.join(os.path.dirname(__file__),"utils", prepare_script)
    circos_script = os.path.join(os.path.dirname(__file__), "utils",circos_script)
    logger.info(f"Processing VCF: {vcf}")
    if outfile_name_mode == "parent":
        name = Path(vcf).parent.name
    else:
        name = re.sub(r'(\.sv)?\.vcf(\.gz)?$', '', Path(vcf).name)
    logger.info(f"Sample name: {name}")
    logger.info(">>> Starting SV Circos Preparation")
    prepare_cmd = [
        "python", prepare_script,
        "--vcf", vcf,
        "--fasta", fasta,
        "--outdir", outdir
    ]
    _run_cmd(prepare_cmd)
    logger.info(">>> Starting Circos Plotting")
    for image_format in image_formats:
        outImage = os.path.join(outdir, f"{name}_sv_circos.{image_format}")
        circos_cmd = [
            "Rscript", circos_script,
            "--input_dir", str(outdir),
            "--output", str(outImage),
            "--genome", genome,
            "--ins_bin_size", ins_bin_size,
            "--ins_plot_type", ins_plot_type,
            "--cytoband", cytoband_file
        ]
        try:
            _run_cmd(circos_cmd)
        except Exception as e:
            logger.error(f"Error occurred while running Circos for {name}: {e}")

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Run SV Circos Pipeline")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA file")
    parser.add_argument("--outdir", required=True, help="Output directory for Circos results")
    parser.add_argument("--vcf_pattern", default="**/unphased/*.vcf*", help="Glob pattern to find VCF files")
    parser.add_argument("--outfile_name_mode", choices=["parent", "sample"], default="sample", help="Naming mode for output directories")
    parser.add_argument("--image_formats", nargs="+", default=["png"], help="Image formats for Circos output (e.g. png pdf svg)")
    parser.add_argument("--vcf", required=True, help="Specific VCF file to process")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    fasta = args.fasta
    vcf = args.vcf
    outdir = args.outdir
    outfile_name_mode = args.outfile_name_mode
    image_formats = args.image_formats
    vcf = args.vcf
    run_cricos_pipeline(vcf=vcf, fasta=fasta, outdir=outdir, image_formats=image_formats, outfile_name_mode=outfile_name_mode)


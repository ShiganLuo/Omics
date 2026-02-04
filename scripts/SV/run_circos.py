
from utils.run_cmd import run_cmd_list
import logging
import re
from pathlib import Path
from typing import Literal
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)
current_dir = Path(__file__).parent.absolute()
def run_cricos_pipelien(
        indir: str,
        fasta: str,
        outdir: str,
        vcf_pattern: str = "**/unphased/*.vcf*",
        outfile_name_mode: Literal["parent","sample"] = "sample",  # parent / sample
        prepare_script: str = "utils/sv_circos_prepare.py",
        circos_script: str = "utils/circos.r",
        ins_bin_size: int = 100000,
        ins_plot_type: Literal["points","bar"] = "bar"
):
    """
    参数化的 SV Circos 绘图流程
    :param indir: 预处理结果输入目录
    :param fasta: 参考基因组 FASTA 文件路径
    :param outdir: Circos 绘图结果输出目录
    :param prepare_script: SV Circos 准备脚本路径
    :param circos_script: SV Circos 绘图脚本路径
    """
    prepare_script = str(Path(current_dir) / prepare_script)
    circos_script = str(Path(current_dir) / circos_script)
    vcfs = [
        str(p) for p in Path(indir).glob(vcf_pattern) 
        if p.suffix in [".vcf", ".gz"] and p.stat().st_size > 0
    ]
    logger.info(f"Found {len(vcfs)} VCF files for Circos preparation.")
    for vcf in vcfs:
        logger.info(f"Processing VCF: {vcf}")
        if outfile_name_mode == "parent":
            name = Path(vcf).parent.name
        else:
            name = re.sub(r'(\.sv)?\.vcf(\.gz)?$', '', Path(vcf).name)
        logger.info(f"Sample name: {name}")
        sample_outdir = Path(outdir) / name
        sample_outdir.mkdir(parents=True, exist_ok=True)
        logger.info(">>> Starting SV Circos Preparation")
        prepare_cmd = [
            "python", prepare_script,
            "--vcf", vcf,
            "--fasta", fasta,
            "--outdir", str(sample_outdir)
        ]
        run_cmd_list(prepare_cmd)
        logger.info(">>> Starting Circos Plotting")
        outImage = sample_outdir / f"{name}_sv_circos.png"
        circos_cmd = [
            "Rscript", circos_script,
            "--input_dir", str(sample_outdir),
            "--output", str(outImage),
            "--genome", "mm39",
            "--ins_bin_size", ins_bin_size,
            "--ins_plot_type", ins_plot_type
        ]
        run_cmd_list(circos_cmd)


if __name__ == "__main__":
    indir = "/data/pub/zhousha/Totipotent20251031/data/Pacbio"
    fasta = "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa"
    outdir = "/data/pub/zhousha/Totipotent20251031/PacBio/circos"
    run_cricos_pipelien(indir, fasta, outdir,vcf_pattern="**/unphased/*.vcf*", outfile_name_mode="sample")
    indir = "/data/pub/zhousha/Totipotent20251031/PacBio/SV"
    outdir = "/data/pub/zhousha/Totipotent20251031/PacBio/circos"
    run_cricos_pipelien(indir, fasta, outdir,vcf_pattern="**/PlaB_only.vcf*", outfile_name_mode="parent")

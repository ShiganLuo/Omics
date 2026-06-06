import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.LogUtil import setup_logger
from common.CmdUtil import _run_cmd
import logging
import re
from pathlib import Path
from typing import List, Literal, Optional
logger = setup_logger("SVCircosPipeline", level=logging.INFO)
PlotFormat = Literal["png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg", "pgf", "raw", "rgba"]
def run_cricos_pipeline(
        vcf: str,
        fasta: str,
        outdir: str,
        image_formats: Optional[List[PlotFormat]] = None,
        outfile_name_mode: Literal["parent","sample"] = "sample",  # parent / sample
        prepare_script: str = "sv_circos_prepare.py",
        circos_script: str = "circos.r",
        ins_bin_size: int = 100000,
        ins_plot_type: Literal["points","bar"] = "bar",
        genome: str = "mm39",
        cytoband_file: str = "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/mm39.cytoBand.txt",
):
    """Run the SV Circos plotting pipeline.

    Prepare SV data from a VCF file and render a Circos plot using
    R scripts in the specified output formats.

    Parameters
    ----------
    vcf : str
        Path to the input VCF file.
    fasta : str
        Path to the reference genome FASTA file.
    outdir : str
        Output directory for Circos results.
    image_formats : list of str, optional
        Image formats for output plots. Defaults to ``["png"]``.
    outfile_name_mode : {"parent", "sample"}, optional
        Naming mode for the output file. ``"parent"`` uses the parent
        directory name; ``"sample"`` uses the VCF filename. Default is
        ``"sample"``.
    prepare_script : str, optional
        Filename of the SV Circos preparation script. Default is
        ``"sv_circos_prepare.py"``.
    circos_script : str, optional
        Filename of the Circos R plotting script. Default is
        ``"circos.r"``.
    ins_bin_size : int, optional
        Bin size (bp) for insertion density plotting. Default is 100000.
    ins_plot_type : {"points", "bar"}, optional
        Plot type for insertion density. Default is ``"bar"``.
    genome : str, optional
        Genome assembly identifier. Default is ``"mm39"``.
    cytoband_file : str, optional
        Path to the cytoband annotation file. Default points to the
        mm39 cytoband file.

    Returns
    -------
    None
    """
    if image_formats is None:
        image_formats = ["png"]
    os.makedirs(outdir, exist_ok=True)
    prepare_script = os.path.join(os.path.dirname(__file__),"utils", prepare_script)
    circos_script = os.path.join(os.path.dirname(__file__), "utils",circos_script)
    logger.info(f"Processing VCF: {vcf}")
    if outfile_name_mode == "parent":
        name = Path(vcf).parent.name
    else:
        name = re.sub(r'(\.sv)?\.vcf(\.gz)?$', '', Path(vcf).name)
    logger.info(f"Sample name: {name}")
    logger.info("Starting SV Circos Preparation")
    prepare_cmd = [
        "python", prepare_script,
        "--vcf", vcf,
        "--fasta", fasta,
        "--outdir", outdir
    ]
    _run_cmd(prepare_cmd)
    logger.info("Starting Circos Plotting")
    for fmt in image_formats:
        outImage = os.path.join(outdir, f"{name}_sv_circos.{fmt}")
        circos_cmd = [
            "Rscript", circos_script,
            "--input_dir", str(outdir),
            "--output", str(outImage),
            "--genome", genome,
            "--ins_bin_size", str(ins_bin_size),
            "--ins_plot_type", ins_plot_type,
            "--cytoband", cytoband_file
        ]
        try:
            _run_cmd(circos_cmd)
            logger.info(f"Circos plot saved to: {outImage}")
        except Exception as e:
            logger.error(f"Error occurred while running Circos for {name}: {e}")

def parse_args():
    """Parse command-line arguments for the SV Circos pipeline.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes: ``fasta``, ``outdir``,
        ``vcf_pattern``, ``outfile_name_mode``, ``formats``, and ``vcf``.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Run SV Circos Pipeline")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA file")
    parser.add_argument("--outdir", required=True, help="Output directory for Circos results")
    parser.add_argument("--vcf_pattern", default="**/unphased/*.vcf*", help="Glob pattern to find VCF files")
    parser.add_argument("--outfile_name_mode", choices=["parent", "sample"], default="sample", help="Naming mode for output directories")
    parser.add_argument("-f", "--format", action="append", dest="formats", metavar="FMT", help="Image output format (png, pdf, svg, ...). Can be specified multiple times. Default: png.")
    parser.add_argument("--vcf", required=True, help="Specific VCF file to process")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    fasta = args.fasta
    vcf = args.vcf
    outdir = args.outdir
    outfile_name_mode = args.outfile_name_mode
    image_formats = args.formats
    vcf = args.vcf
    run_cricos_pipeline(vcf=vcf, fasta=fasta, outdir=outdir, image_formats=image_formats, outfile_name_mode=outfile_name_mode)


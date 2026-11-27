import pandas as pd
import pybedtools
import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from typing import Union
import glob
import re
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def query_gene_for_vcf(
        vcf:str,
        gtf:str,
        outfile:str
):
    vcf = pd.read_csv(vcf, comment="#", sep="\t",
                    usecols=[0,1,3,4], names=["chrom","pos","ref","alt"])

    vcf["start"] = vcf["pos"] - 1
    vcf["end"] = vcf["pos"]

    vcf_bed = pybedtools.BedTool.from_dataframe(vcf[["chrom","start","end","ref","alt"]])
    gtf_bed = pybedtools.BedTool(gtf)

    res = vcf_bed.intersect(gtf_bed, wa=True, wb=True).to_dataframe()
    res.to_csv(outfile,sep="\t",index=False)

def gene_pred(
        gtf:Union[str,Path], 
        db_dir:Union[str,Path],
        buildver:Union[str,Path]
        ):
    """
    Python 等价函数：
    gtfToGenePred -genePredExt <gtf> <refGene>

    Args:
        gtf: 输入 GTF 文件 (str | Path)
        db_dir: 输出目录 (str | Path)
        buildver: 自定义版本
    return:
        db_dir / f"{buildver}_refGene.txt"
    """
    gtf = Path(gtf)
    db_dir = Path(db_dir)
    ref_gene = db_dir / f"{buildver}_refGene.txt"

    ref_gene.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "gtfToGenePred",
        "-genePredExt",
        str(gtf),
        str(ref_gene)
    ]

    subprocess.run(cmd, check=True)
    return ref_gene


def retrieve(
        retrieve_pl:Union[str,Path], 
        genome:Union[str,Path], 
        ref_gene:Union[str,Path], 
        db_dir:Union[str,Path],
        buildver:Union[str,Path]
        ):
    """
    Python 等价函数：
    perl retrieve_seq_from_fasta.pl --format refGene --seqfile <genome> <refGene> --out <refGeneMrna>

    Args:
        retrieve_pl: perl 脚本路径 (str | Path)
        genome: 基因组 fasta (str | Path)
        ref_gene: refGene.txt (str | Path)
        db_dir: 输出目录 (str | Path)
        buildver: 自定义版本
    return:
        db_dir / f"{buildver}_refGeneMrna.fa"
    """
    ref_gene = Path(ref_gene)
    db_dir = Path(db_dir)
    ref_gene_mrna = db_dir / f"{buildver}_refGeneMrna.fa"

    ref_gene_mrna.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "perl",
        str(retrieve_pl),
        "--format", "refGene",
        "--seqfile", str(genome),
        str(ref_gene),
        "--out", str(ref_gene_mrna)
    ]

    subprocess.run(cmd, check=True)
    return ref_gene_mrna
def convert_to_avinput(
        convert_pl:Union[str,Path], 
        input_file:Union[str,Path], 
        out_dir:Union[str,Path]
        ):
    """
    将 VCF/BED 转换为 ANNOVAR .avinput 文件
    """
    input_file = Path(input_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / (input_file.stem + ".avinput")
    logger.info(f"Converting {input_file.name} -> {out_file.name}")
    with open(out_file, "w") as out_f:
        subprocess.run([
            "perl",
            str(convert_pl),
            "-format", "vcf4",
            "-withfreq",
            str(input_file)
        ], stdout=out_f, check=True)
    return out_file

def annotate_avinput(
        table_pl:Union[str,Path], 
        avinput_file:Union[str,Path], 
        db_dir:Union[str,Path], 
        buildver:str, 
        out_dir:Union[str,Path]):
    """
    调用 table_annovar.pl 注释 .avinput 文件
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / avinput_file.stem
    logger.info(f"Annotating {avinput_file.name} -> {out_prefix}_*.csv")
    subprocess.run([
        "perl",
        str(table_pl),
        str(avinput_file),
        str(db_dir),
        "-buildver", buildver,
        "-out", str(out_prefix),
        "-remove",
        "-protocol", "refGene",
        "-operation", "g",
        "-nastring", ".",
        "-csvout"
    ], check=True)
    return out_prefix

def process_file(file_path:Union[str,Path], 
                 convert_pl:Union[str,Path], 
                 table_pl:Union[str,Path], 
                 db_dir:Union[str,Path], 
                 buildver:str, 
                 avinput_outdir:Union[str,Path], 
                 table_outdir:Union[str,Path]):
    file_path = Path(file_path)
    if db_dir.exists():
        pass
    else:
        raise ValueError(f"db_dir is not exists:{str(db_dir)}")
    try:
        avinput_file = convert_to_avinput(convert_pl, file_path, avinput_outdir)
        annotate_avinput(table_pl, avinput_file, db_dir, buildver, table_outdir)
        logger.info(f"Finished {file_path.name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error processing {file_path.name}: {e}")

def prepare_db(
        gtf:Union[str,Path],
        retrieve_pl:Union[str,Path],
        genome:Union[str,Path],
        db_dir:Union[str,Path], 
        buildver:str, 
        ):
    gtf = Path(gtf)
    try:
        ref_gene = gene_pred(gtf,db_dir,buildver)
        retrieve(retrieve_pl,genome,ref_gene,db_dir,buildver)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error processing {gtf.name}: {e}")

def batch_process(input_dir:str, 
                  convert_pl:str, 
                  table_pl:str, 
                  db_dir:str, 
                  buildver:str, 
                  avinput_outdir:str, 
                  table_outdir:str, 
                  threads=4):
    """
    批量处理整个文件夹中的 VCF/BED 文件，支持多线程
    """
    input_dir = Path(input_dir)
    files = [f for f in input_dir.glob("*") if f.suffix.lower() in [".vcf", ".bed"]]
    logger.info(f"Found {len(files)} files to process in {input_dir}")

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(process_file, f, convert_pl, table_pl, db_dir, buildver, avinput_outdir, table_outdir)
            for f in files
        ]
        for future in as_completed(futures):
            future.result()

if __name__ == "__main__":
    vcf_dict = {
        "TLSC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/GSE166216/toti_only.vcf.gz",
        "ciTotiSC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/GSE185005/toti_only.vcf.gz",
        "ci8CLC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/GSE204801/toti_only.vcf.gz",
        "hTBLC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/GSE224794/toti_only.vcf.gz"
    }

    mouse = ["TLSC","ciTotiSC"]
    human = ["ci8CLC","hTBLC"]

    gtf_GRCh38 = "/disk5/luosg/Reference/GENCODE/human/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf"
    gtf_GRCm39 = "/disk5/luosg/Reference/GENCODE/mouse/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf"
    TEgtf_GRCh38 = "/disk5/luosg/Reference/GENCODE/human/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf"
    TEgtf_GRCm39 = "/disk5/luosg/Reference/GENCODE/mouse/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf"
    outdir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotation"
    query_gene_for_vcf("/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/ci8CLC/toti_only.vcf.gz",
                       TEgtf_GRCh38,"a.bed")
    
    convert_pl = "/opt/annovar/convert2annovar.pl"
    table_pl = "/opt/annovar/table_annovar.pl"
    
    db_dir_GRCh38 = "/disk5/luosg/Reference/GENCODE/human/GRCh38/annovar/GRCh38"
    buildver_GRCh38 = "GRCh38"
    db_dir_GRCm39 = "/disk5/luosg/Reference/GENCODE/mouse/GRCm39/annovar/GRCm39"
    buildver_GRCm39 = "GRCm39"

    threads = 4
    

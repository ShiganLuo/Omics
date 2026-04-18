import pandas as pd
import sys
import os
import glob
from typing import Optional
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.SepUtil import detect_delimiter
import logging
logger = logging.getLogger(__name__)
def export_md5(
        input_file:str, 
        output_file:str,
        col_number:int=3
    ):
    """
    Function: Export MD5 checksums from a TSV file to a text file in the format "MD5  filename".
    Parameters:
    - input_file: Path to the input TSV file containing file names and their corresponding MD5 checksums.
    - output_file: Path to the output text file where the MD5 checksums and file names will be written.
    Note:
    - filex
    """
    df = pd.read_csv(input_file, sep='\t')
    md5_lines = []
    relation_cols = [(f"file{i}_name", f"file{i}_md5") for i in range(1, col_number + 1)]
    for i in range(len(df)):
        for file_col, md5_col in relation_cols:
            file = str(df.at[i, file_col]) if file_col in df.columns else ''
            md5 = str(df.at[i, md5_col]) if md5_col in df.columns else ''
            if file and file != 'nan' and md5 and md5 != 'nan':
                md5_lines.append(f"{md5}  {file}")
    with open(output_file, 'w') as f:
        f.write('\n'.join(md5_lines) + '\n')

def meta_generate(
    meta_experiment_file:str,
    fq_dir:str,
    outfile:str,
    condition_col:str = "experiment_title",
    condition_value:Optional[str] = "Cell line co-culture (mouse and human ips)",
    sample_id_col:str = "sample_name",
    data_id_col:str = "sample_name",
):
    """
    Function: Generate pipeline metadata from the cngb experiment metadata file.
    Parameters:
    - meta_experiment_file: Path to the metadata file containing experiment information.
    - fq_dir: Path to the directory containing the fastq files.
    - condition_col: Column name in the metadata file to filter by condition (default: "Experiment Title").
    - condition_value: Value in the condition column to filter the metadata (default: "Cell line co-culture (mouse and human ips)").
    Returns:
    - A DataFrame containing the filtered metadata with added file paths for the fastq files.
    """
    meta_cols = ["sample_id", "data_id","design", "fastq_1", "fastq_2"]
    sep = detect_delimiter(meta_experiment_file)
    df = pd.read_csv(meta_experiment_file, sep=sep)
    if condition_value:
        df = df[df[condition_col] == condition_value]
    meta_dict = {col: [] for col in meta_cols}
    for _,row in df.iterrows():
        meta_dict["sample_id"].append(row[sample_id_col])
        meta_dict["data_id"].append(row[data_id_col])
        fastq_files = glob.glob(os.path.join(fq_dir, f"{row[data_id_col]}*.fq.gz"))
        fastq_1 = fastq_files[0] if len(fastq_files) > 0 else None
        fastq_2 = fastq_files[1] if len(fastq_files) > 1 else None
        meta_dict["fastq_1"].append(fastq_1)
        meta_dict["fastq_2"].append(fastq_2)
        meta_dict["design"] = None
    df_out = pd.DataFrame(meta_dict)
    df_out.to_csv(outfile,sep="\t",index=False)
    
    pass
if __name__ == "__main__":
    # meta = "/data/pub/zhousha/20260411_RNAseq/data/meta/metadata_CNP0003135_experiment.tsv"
    # outfile = "/data/pub/zhousha/20260411_RNAseq/data/fastq/md5sum.txt"
    # export_md5(meta, outfile, col_number=3)
    meta = "/data/pub/zhousha/20260411_RNAseq/data/meta.tsv"
    outfile = "/data/pub/zhousha/20260411_RNAseq/data/meta_input.tsv"
    fqDir = "/data/pub/zhousha/20260411_RNAseq/data/fastq"
    meta_generate(
        meta_experiment_file=meta,
        outfile=outfile,
        fq_dir=fqDir,
        condition_value=None
    )
    pass
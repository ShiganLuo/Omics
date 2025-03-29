import pandas as pd
from intervaltree import IntervalTree
from typing import Callable
import re
import argparse

def parse_gtf(gtf_file:str, gene_set:set):
    """
    parse_gtf function: parase gtf file.The gtf file must have the gene_id attribute. sep with tab.
    and it should only contain one featrue type in the 3rd column in case of the bed overlap between different feature types.
    if fact it's no need to worry about that some rows's region is overlapped as long as it's real
    """
    print("parse_gtf")
    bed_records = []
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue  # 跳过注释行
            fields = line.strip().split("\t")
            feature_type = fields[2]  # 可能是 exon 或 TE
            
            # 解析 gene_id
            match = re.search(r'gene_id "([^"]+)"', fields[8])
            if match:
                gene_id = match.group(1)
                # print(gene_id)
                if gene_id in gene_set:
                    # print(gene_id)
                    chrom = fields[0]
                    start = int(fields[3]) - 1  # GTF 是 1-based，BED 需要 0-based
                    end = int(fields[4])
                    strand = fields[6]
                    bed_records.append((chrom, start, end, gene_id, feature_type, strand))
    return bed_records
def parase_gtfForTEcount(gtf_file:str, gene_set:set):
    """
    parase_gtfForTEcount function: parase gtf file.The gtf file must have the gene_id,family_id,class_id(because the first column of
    of TEcount is gene_id:family_id:class_id) in the ninth column. sep with tab.
    and it should only contain one featrue type in the 3rd column in case of the bed overlap between different feature types.
    """
    print("parase_gtfForTEcount")
    bed_records = []
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue  # 跳过注释行
            fields = line.strip().split("\t")
            feature_type = fields[2]  # 可能是 exon 或 TE
            
            # 解析 gene_id
            match = re.search(r'gene_id "([^"]+)";.*?family_id "([^"]+)";.*?class_id "([^"]+)";', fields[8])
            if match:
                attr_dict = {
                    "gene_id": match.group(1),
                    "family_id": match.group(2) or "NA",
                    "class_id": match.group(3) or "NA"
                }
                combined_id = f"{attr_dict['gene_id']}:{attr_dict['family_id']}:{attr_dict['class_id']}"

                if combined_id in gene_set:
                    # print(combined_id)
                    chrom = fields[0]
                    start = int(fields[3]) - 1  # GTF 是 1-based，BED 需要 0-based
                    end = int(fields[4])
                    strand = fields[6]
                    bed_records.append((chrom, start, end, combined_id, feature_type, strand))
    return bed_records

def getBed(infile:str,outfile:str,gtf:str,writeMode:str,func:Callable[[str,set],list]=parse_gtf):
    """
    getBed function: get the bed file of gene_list from the gtf file by calling the parse_gtf or parase_gtfForTEcount function.
    the first columan of infile must contain the gene_id.
    the mode parameter is the mode of open the outfile.
    """
    if func is None:
        raise ValueError("A parsing function must be provided")
    with open(infile, "r") as f:
        gene_ids = set(line.strip().split("\t")[0].strip('"') for line in f)
    print(gene_ids)
    bed = func(gtf, gene_ids)
    bed.sort()
    with open(outfile, writeMode) as f:
        for record in bed:
            f.write("\t".join(map(str, record)) + "\n")

    print(f"BED file: {outfile}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extract the gene region for TE and exon")
    parser.add_argument('--input', type=str, required=True, help='Path to input file')
    parser.add_argument('--output', type=str, required=True,action='append', help='Path to output file')
    parser.add_argument('--exonGtf', type=str, required=True, help='Path to exon gtf,must only contain exon annotation')
    parser.add_argument('--TEGtf', type=str, required=True, help='Path to TE gtf')
    parser.add_argument('--mode', type=str,default="SNP", choices=["SNP", "StringTie"],help='how to write file,combine or separate')
    args = parser.parse_args()
    # exon_gtf = "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/exon.gtf"
    # TE_gtf = "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf"
    # countfile = "/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/RNASNP202503TEcount_common.cntTable"
    # outfile = "/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/RNASNP202503TEcount_common.bed"
    infile = args.input
    outfileList = args.output
    exon_gtf = args.exonGtf
    TE_gtf = args.TEGtf
    mode = args.mode
    if mode == "SNP":
        outfile = outfileList[0]
        getBed(infile,outfile,exon_gtf,writeMode="w",func=parse_gtf)
        getBed(infile,outfile,TE_gtf,writeMode="a+",func=parase_gtfForTEcount)
    elif mode == "StringTie":
        if len(outfileList) == 2:
            getBed(infile,outfileList[0],exon_gtf,writeMode="w",func=parse_gtf)
            getBed(infile,outfileList[1],TE_gtf,writeMode="w",func=parase_gtfForTEcount)
        else:
            print("Mode StringTie must provide two input addresses, the first is the regular gtf, the second is the TE gtf")
            exit(1)
    else:
        print("please provide corrct mode,SNP or StringTie")
        exit(1)
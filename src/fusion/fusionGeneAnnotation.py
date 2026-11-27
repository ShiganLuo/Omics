import pandas as pd
from intervaltree import IntervalTree
from typing import Dict
import argparse
def annTree(annFile:str):
    """
    annTree function: create a IntervalTree for the annotation file(gtf format).
    """
    print("annotationTree function")
    ann_df = pd.read_csv(annFile, sep="\t", header=None,comment='#')
    tree = {}
    for _, row in ann_df.iterrows():
        chr = row[0]
        feature = row[2]
        if chr not in tree:
            tree[chr] = IntervalTree()
            gene_id = row[8].split('"')[1]
            start = row[3]
            end = row[4]
            if start < end:
                tree[chr].addi(start, end, gene_id)
            else:
                print(f"{gene_id} some row annotate false,start >= end")
    return tree
# Read the input bed file and annotate it
def query_annotation(targetFile:str, tree: Dict[str, IntervalTree],outfile:str):
    """
    query_annotation function: query the annotation file with the bed file
    col: 1:chr,2:start,3:end
    """
    print("query_annotation function")
    df_target = pd.read_csv(targetFile, sep="\t", header=None)
    df_unique = df_target.drop_duplicates(subset=df_target.columns[:3].tolist())
    print(len(df_target), len(df_unique))
    results = []
    for _, row in df_unique.iterrows():
        chr, start, end = row[0], row[1], row[2]
        start = start + 1
        if chr in tree:
            print(chr)
            intervals = tree[chr].overlap(start, end)
            for interval in intervals:
                annotation_data = interval.data 
                result_row = {
                    'chr': chr,
                    'start': start,
                    'end': end,
                    'gene_id': annotation_data
                }
                results.append(result_row)
    df_result = pd.DataFrame(results)
    df_result.to_csv(outfile, sep="\t", index=False, header=False)
    return df_result
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extract the gene region for TE and exon")
    parser.add_argument('--input', type=str, required=True, help='Path to input file')
    parser.add_argument('--output', type=str, required=True, help='Path to output file')
    parser.add_argument('--gtf', type=str, required=True, help='Path to gtf file')
    args = parser.parse_args()
    infile = args.input
    outfile = args.output
    annfile = args.gtf
    print("------begin annotation------")
    Tree = annTree(annfile)
    query_annotation(infile, Tree, outfile)
    print("------annotation done------")
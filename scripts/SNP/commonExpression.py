import pandas as pd
import argparse

def commonExpression(infile:str,outfile:str):
    """
    commonExpression function
    This script is used to filter out the rows with zeros in the columns except the first column.
    use numpy to speed up the process.
    """
    df = pd.read_csv(infile,sep="\t",header=0)
    mask = ~(df.iloc[:, 1:].to_numpy() == 0).any(axis=1)
    df_common = df.loc[mask]
    df_common.to_csv(outfile,sep="\t",index=False)
    return outfile
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extract the gene row which expressed in all samples")
    parser.add_argument('--input', type=str, required=True, help='Path to input file')
    parser.add_argument('--output', type=str, required=True, help='Path to output file')
    args = parser.parse_args()
    # infile = "/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/RNASNP202503TEcount.cntTable"
    # outfile = "/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/RNASNP202503TEcount_common.cntTable"
    infile = args.input
    outfile = args.output
    # print(infile)
    # print(outfile)
    commonExpression(infile,outfile)
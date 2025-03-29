import pandas as pd
import os
from functools import reduce
import argparse


def combineTE(infiles: list[str],outfile: str,split:str):
    """
    combineTEcount function
    This script is used to combine the {sample_id}TEcount.cntTable files from different samples into one file and
    rename the second column to {sample_id}.
    """
    dfList = []
    for infile in infiles:
        newColName = os.path.basename(infile).split(f"{split}")[0]
        print(newColName)
        #默认会去除引号，可以使用 quotechar 参数来指定如何处理引号，quotechar='"'保留
        df = pd.read_csv(infile,sep="\t",header=0)
        df.columns.values[1] = newColName
        dfList.append(df)
    combineDf = reduce(lambda left, right: pd.merge(left, right, on='gene/TE', how='inner'), dfList)
    combineDf.to_csv(outfile,sep="\t",index=False)
    return combineDf

def fileList(path:str,pattern:str):
    """
    fileList function
    This function is used to get all the {sample_id}TEcount.cntTable files from the given path.
    """
    from pathlib import Path
    matching_files = Path(path).rglob(f"*{pattern}")
    return [str(f) for f in matching_files]
def getTE(df:pd.DataFrame,outfile:str):
    pattern = r'^[^:]+:[^:]+:[^:]+:[^:]+$'
    df_TE = df[df['gene/TE'].str.match(pattern,na=False)]
    df_TE.to_csv(outfile,sep="\t",index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to combine the output of TEcount and TElocal.")
    parser.add_argument('-p','--procedure',type=str,required=True,choices=['TEcount', 'TElocal','TEcountStringTie'],
                        help="name of your procedur,TEcount or TElocal")
    parser.add_argument('-i','--inputDir',type=str,required=True,help="Directory path of your input file")
    parser.add_argument('-o','--output',type=str,required=True,help="Path to your output file")
    args = parser.parse_args()
    if args.procedure == "TEcount":
        ## TEcount
        infiles = fileList(args.inputDir,"TEcount.cntTable")
        print(f"TEcount.cntTable file path : {infiles}")
        combineTE(infiles,args.output,"TEcount")
        print(f"write file into {args.output}")
    elif args.procedure == "TElocal":
        ### TElocal
        infiles = fileList(args.inputDir,"TElocal.cntTable")
        print(f"TEcount.cntTable file path : {infiles}")
        combineTE(infiles,args.output,"TElocal")
        print(f"write file into {args.output}")
        # df = pd.read_csv(outfile,sep="\t",header=0)
        # outfile = "output/RNASNP202503TElocal_TE.cntTable"
        # getTE(df,outfile)
    elif args.procedure == "TEcountStringTie":
        ## TEStringtie
        infiles = fileList(args.inputDir,"TEcountStringTie.cntTable")
        print(f"TEcountStringTie.cntTable file path : {infiles}")
        combineTE(infiles,args.output,"TEcountStringTie")
        print(f"write file into {args.output}")
    else:
        print(f"Please check your procedure name {args.procedure},muste be TEcount or TElocal or TEcountStringTie")
        exit(1)

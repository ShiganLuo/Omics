import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt

def caculateSnpForEachGene(controlFiles:str,experimentFiles:str,geneId:str="ENSMUSG"):
    """
    caculateSnpForEachGene
    This script is used to caculate the number of total,TE and exon records in each gene and group them into a DataFrame 
    fron the Intersectbed files. 
    infile:
     - sep with tab
     - the thirteen column of the bed file is used to distinguish TE and exon records.
     - TE records have 2 colons in the thirteen column, 
     - while exon records start with "ENSMUSG" on behalf of mouse."ENSG":human
    """
    control = []
    for infile in controlFiles:
        df = pd.read_csv(infile,sep="\t",header=None)
        df_Gene = df[df[13].str.startswith(geneId, na=False)]
        control.append(df_Gene)
    df_control = pd.concat(control,axis=0)
    ans_control = df_control[13].value_counts()

    experiment = []
    for infile in experimentFiles:
        df = pd.read_csv(infile,sep="\t",header=None)
        df_Gene = df[df[13].str.startswith(geneId, na=False)]
        experiment.append(df_Gene)
    df_experiment = pd.concat(experiment,axis=0)
    ans_experiment = df_experiment[13].value_counts()
    ans = pd.merge(ans_control.reset_index(), ans_experiment.reset_index(), how='outer', on=13)
    ans = ans.fillna(0)
    ans.columns = ['gene_id', 'control', 'experiment']
    return ans
def plotSnpForEachGene(df_SNP:pd.DataFrame,annFile:str,outImage:str,outfile:str=None):
    """
    plotSnpForEachGene function:
    plot the result of caculateSnpForEachGene
    """
    df_ann = pd.read_csv(annFile,sep="\t")
    df = pd.merge(df_SNP,df_ann, how='inner',on='gene_id')
    df = df[['gene_id','gene_name','gene_type','control','experiment']]
    df['difference'] = df['experiment'] - df['control']
    if outfile is not None:
        df.to_csv(outfile,sep="\t",index=False)
    df_up = df[df['difference'] > 0]
    df_down = df[df['difference'] < 0]

    thresholdUp = df_up['difference'].quantile(0.95)
    thresholdDown = df_down['difference'].quantile(0.05)

    df_up = df_up[df_up['difference'] > thresholdUp]
    df_down = df_down[df_down['difference'] < thresholdDown]

    df_up = df_up.sort_values('difference',ascending=False)
    df_down = df_down.sort_values('difference',ascending=True)

    n_up = len(df_up)
    n_down = len(df_down)
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(14, 18),dpi=300)
    ax1.bar(df_up['gene_name'], df_up['difference'],color='red',label=f"Increased mutations gene number: {n_up}")
    ax1.set_xticks(range(0, len(df_up['gene_name']), 25))
    ax1.set_xticklabels(range(0, len(df_up['gene_name']), 25), fontsize=22)
    ax1.tick_params(axis='y', labelsize=22) 
    ax1.legend(loc='upper right',fontsize=25)
    ax1.set_ylabel("count",fontsize=25)
    ax1.set_xlabel("Increased mutations gene",fontsize=25)
    ax2.bar(df_down['gene_name'],df_down['difference'],color='blue',label=f"Reduced mutations gene number: {n_down}")
    
    ax2.set_xticks(range(0, len(df_down['gene_name']), 25))
    ax2.set_xticklabels(range(0, len(df_down['gene_name']),25), fontsize=22)
    ax2.tick_params(axis='y', labelsize=22) 
    ax2.legend(loc='lower right',fontsize=25)
    ax2.set_ylabel("count",fontsize=25)
    ax2.set_xlabel("Reduced mutations gene",fontsize=25)
    plt.tight_layout()
    fig.savefig(outImage)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="plot common vcf and annovar file")
    parser.add_argument('--control', type=str, required=True,nargs='+', help='Paths to control files')
    parser.add_argument('--experiment', type=str, required=True,nargs='+', help='Path to experimet files')
    parser.add_argument('--outdir', type=str,required=True,help='Path to output directory')
    parser.add_argument('--species', type=str,default="mouse",help='species name mouse or human')
    args = parser.parse_args()
    controlFiles = args.control
    experimentFiles = args.experiment
    species = args.species
    outdir = args.outdir
    outfile = outdir + "alteredGene.csv"
    outImage = outdir + "alteredGene.jpeg"

    if species == "mouse":
        ans = caculateSnpForEachGene(controlFiles,experimentFiles,geneId="ENSMUSG")
        annFile = "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv"
        plotSnpForEachGene(ans,annFile,outImage,outfile)
    elif species == "human":
        ans = caculateSnpForEachGene(controlFiles,experimentFiles,geneId="ENSG")
        annFile = "/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/geneIDAnnotation.csv"
        plotSnpForEachGene(ans,annFile,outImage,outfile)
    else:
        raise ValueError("请输入正确的物种,mouse or human")


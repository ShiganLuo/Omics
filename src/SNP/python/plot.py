import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import re
import argparse
def CaculateNumber(controlFiles:list[str],experimentFiles:list[str],controlName:str,experiemtName:str,outfile:str=None,geneId:str="ENSMUSG"):
    """
    CaculateNumber function"
    This script is used to caculate the number of total,TE and exon records in each sample and group them into a DataFrame 
    fron the Intersectbed files. 
    infile:
     - sep with tab
     - the thirteen column of the bed file is used to distinguish TE and exon records.
     - TE records have 2 colons in the thirteen column, 
     - while exon records start with "ENSMUSG" on behalf of mouse.ENSG:human
    """
    rows = []
    
    for infile in controlFiles:
        df = pd.read_csv(infile,sep="\t",header=None)
        countTE = df[df[13].str.count(":") == 2].shape[0]
        countExon = df[df[13].str.startswith(geneId, na=False)].shape[0]
        countTotal = df.shape[0]
        if countTE + countExon != countTotal:
            raise ValueError(f"Error: {infile} has {countTotal} records, but {countTE} TE records and {countExon} Gene records")
        sample = os.path.basename(infile).split("Common")[0]
        rows.append({"Sample":sample,"Total":countTotal,"TE":countTE,"Gene":countExon,"group":f"{controlName}"})
        
    for infile in experimentFiles:
        df = pd.read_csv(infile,sep="\t",header=None)
        countTE = df[df[13].str.count(":") == 2].shape[0]
        countExon = df[df[13].str.startswith(geneId, na=False)].shape[0]
        countTotal = df.shape[0]
        if countTE + countExon != countTotal:
            raise ValueError(f"Error: {infile} has {countTotal} records, but {countTE} TE records and {countExon} Gene records")
        sample = os.path.basename(infile).split("Common")[0]
        rows.append({"Sample":sample,"Total":countTotal,"TE":countTE,"Gene":countExon,"group":f"{experiemtName}"})
    df = pd.DataFrame(rows)
    if outfile != None:
        df.to_csv(outfile,index=False)
    return df

def plotNumber(df:pd.DataFrame,outfile:str,controlName:str,experimentName:str):
    """
    plotNumber function: plot the output of CaculateNumber function
    """
    df_melted = df.melt(id_vars=["Sample", "group"], value_vars=["TE", "Gene"], var_name="Category", value_name="Count")
    group_colors = {
    controlName: "blue",
    experimentName: "red"
    }

    plt.figure(figsize=(14, 10),dpi=300)
    sns.set_style("whitegrid")
    ax = sns.barplot(x="Sample", y="Count", hue="Category", data=df_melted, palette="Set2")
    handles_cat, labels_cat = ax.get_legend_handles_labels()
    for xtick, sample in zip(ax.get_xticklabels(), df["Sample"]):
        group = df.loc[df["Sample"] == sample, "group"].values[0]
        xtick.set_color(group_colors[group])
    
    handles_group = [mpatches.Patch(color=group_colors[controlName], label=controlName),
           mpatches.Patch(color=group_colors[experimentName], label=experimentName)]
    legend1 = plt.legend(handles=handles_cat, title="feature", loc="upper right",bbox_to_anchor=(1.1, 1),fonstsize = 12)
    legend2 = plt.legend(handles=handles_group, title="Group", loc="lower right",bbox_to_anchor=(1.1, 1),fonstsize = 12)
    plt.gca().add_artist(legend1)

    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Sample", fontsize=15)
    plt.ylabel("Count", fontsize=15)
    plt.savefig(outfile)

def TopAlterFeature(controlFiles:list[str],experimentFiles:list[str],controlName:str,experiemtName:str,number:int,outfile:str=None):
    """
    GeneStatistic function: get The top alterned feature from controlFiles and experimentFiles.
    infile: The thirteenth is feature,file sep with tab.
    parameter: number decide get top number
    """
    def getTop(df:pd.DataFrame,number:int,group:str):
        pattern = r'^[^:]+:[^:]+:[^:]+$'
        geneSNPCount = df[13].value_counts()
        seriesTE = geneSNPCount[geneSNPCount.index.str.match(pattern, na=False)]
        seriesGene = geneSNPCount[~geneSNPCount.index.str.match(pattern, na=False)]
        seriesTE = seriesTE.sort_values(ascending = False)
        seriesGene = seriesGene.sort_values(ascending = False)
        topTE = seriesTE.head(number)
        topTE_df = pd.DataFrame({
            'feature' : topTE.index,
            'type' : ['TE'] * len(topTE),
            'number' : topTE.values
        })
        topGene = seriesGene.head(number)
        topGene_df = pd.DataFrame({
            'feature' : topGene.index,
            'type' : ['gene_id']*len(topGene),
            'number': topGene.values
        })
        result = pd.concat([topTE_df,topGene_df],axis=0)
        result['group'] = group
        return result

    control_list = []
    for infile in controlFiles:
        df = pd.read_csv(infile,sep="\t",header=None)
        control_list.append(df)
    df_control = pd.concat(control_list,axis=0)
    ans1 = getTop(df_control,number,controlName)
    expreiment_list = []
    for infile in experimentFiles:
        df = pd.read_csv(infile,sep="\t",header=None)
        expreiment_list.append(df)
    df_experiment = pd.concat(expreiment_list,axis=0)
    ans2 = getTop(df_experiment,number,experiemtName)
    ans = pd.concat([ans1,ans2],axis=0,ignore_index=True)
    if outfile != None:
        ans.to_csv(outfile,index=False)
    return ans

def plotTopfeature(df:pd.DataFrame,outfile:str):
    """
    plotTopfeatrue function:Comparison of top altered TE and gene between Groups

    """ 
    df_gene = df[df["type"] == "gene_id"].groupby("feature")["number"].sum().reset_index()
    df_te = df[df["type"] == "TE"].groupby("feature")["number"].sum().reset_index()

    df_gene = df_gene.sort_values(by="number", ascending=False)
    df_te = df_te.sort_values(by="number", ascending=False)

    sorted_features = list(df_gene["feature"]) + list(df_te["feature"])

    df["feature"] = pd.Categorical(df["feature"], categories=sorted_features, ordered=True)
   
    
    plt.figure(figsize=(18, 12),dpi=300)
    ax = sns.barplot(data=df, x="feature", y="number", hue="group", dodge=True)
    group_colors = {
    "Gene": "blue",
    "TE": "red"
    }
    
    pattern = r'^[^:]+:[^:]+:[^:]+$'
    for xtick in ax.get_xticklabels():
        feature = xtick.get_text()
        if re.match(pattern, feature):
            print(f"{group_colors['TE']}:",feature)
            xtick.set_color(group_colors['TE'])
        else:
            print(f"{group_colors['Gene']}:",feature)
            xtick.set_color(group_colors['Gene'])  
    handles_cat, labels_cat = ax.get_legend_handles_labels()
    handles_group = [mpatches.Patch(color=group_colors["Gene"], label="Gene"),
                     mpatches.Patch(color=group_colors["TE"], label="TE")]
    legend1 = plt.legend(handles=handles_cat, title="Group", loc="upper right",bbox_to_anchor=(1, 1))
    legend2 = plt.legend(handles=handles_group, title="feature", loc="lower right",bbox_to_anchor=(1, 0.2))
    plt.gca().add_artist(legend1)
    plt.xticks(rotation=25)
    plt.xlabel("Feature")
    plt.ylabel("Number")
    plt.title("Comparison of top altered TE and gene between Groups")
    plt.savefig(outfile)
   
def caculateIntervalDistribution(controlFiles:list[str],experimentFiles:list[str],controlNmae:str,experimentName:str,outfile:str=None):
    """
    caculateIntervalDistribution function: make a statistics of the snp distributition,such as ncRNA_exonic,UTR5
    infile:standard multianno.csv produced by the table annotation function of annovar。sep with comma(,)
    no need to pay attentition to TE,because there is no TE in function annotation by annovar
    """
    def csv_read(infile:str):
        data = []
        with open(infile, "r") as file:
            next(file)
            for line in file:
                columns = line.strip().split(",")
                if len(columns) >= 6:
                    data.append(columns[5])
        return pd.Series(data).value_counts()
    controlList = []
    for infile in controlFiles:
        df = pd.read_csv(infile,sep=",")
        funcIntervalSNP = df['Func.refGene'].value_counts()
        sample = os.path.basename(infile).split(".GRCm39_multianno")[0]
        
        df_count = pd.DataFrame({
            'Func' : funcIntervalSNP.index,
            'Number': funcIntervalSNP.values,
            'Sample' : [sample]*len(funcIntervalSNP)
        })
        controlList.append(df_count)
    df_control = pd.concat(controlList,axis=0,ignore_index=True)
    df_control['group'] = controlNmae
    experimentList = []
    for infile in experimentFiles:
        df = pd.read_csv(infile,sep=",")
        funcIntervalSNP = df['Func.refGene'].value_counts()
        sample = os.path.basename(infile).split(".GRCm39_multianno")[0]
        
        df_count = pd.DataFrame({
            'Func' : funcIntervalSNP.index,
            'Number': funcIntervalSNP.values,
            'Sample' : [sample]*len(funcIntervalSNP)
        })
        experimentList.append(df_count)
    df_experiment = pd.concat(experimentList,axis=0,ignore_index=True)
    df_experiment['group'] = experimentName
    newdf = pd.concat([df_control,df_experiment],axis=0,ignore_index=True)
    if outfile != None:
        newdf.to_csv(outfile,index=False)
        
    return newdf

def plotIntervalDistribution(df:pd.DataFrame,outfile:str):

    groups = df["group"].unique()
    num_groups = len(groups)
    fig, axes = plt.subplots(1, num_groups, figsize=(8 * num_groups, 6),dpi=300)

    if num_groups == 1:
        axes = [axes]
    cmap = plt.colormaps.get_cmap("tab20") 
    unique_funcs = df["Func"].unique()
    color_map = {func: cmap(i / len(unique_funcs)) for i, func in enumerate(unique_funcs)}
    def autopct_format(pct):
        return f'{pct:.1f}%' if pct >= 5 else ''  
    for ax, group in zip(axes, groups):
        subset = df[df["group"] == group]
        func_counts = subset.groupby("Func")["Number"].sum()

       
        colors = [color_map[func] for func in func_counts.index]
        wedges, texts, autotexts = ax.pie(
            func_counts, 
            autopct=autopct_format,
            startangle=140, 
            colors=colors, 
            wedgeprops={'edgecolor': 'black'},
            labeldistance=1.2 
        )
        ax.set_title(f"Functional Distribution\n{group}")
        ax.set_ylabel("")
    fig.legend(
        wedges, func_counts.index, 
        title="Func", 
        loc="center",
        bbox_to_anchor=(0.5, 0.5),
        ncol=1,
    )
    
    plt.tight_layout()
    fig.savefig(outfile)




        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="plot common vcf and annovar file")
    parser.add_argument('--mode', type=str,required=True,default="vcf",choices=["vcf","annovar"],help='experiment group name')
    parser.add_argument('--control', type=str, required=True,nargs='+', help='Paths to control files')
    parser.add_argument('--experiment', type=str, required=True,nargs='+', help='Path to experimet files')
    parser.add_argument('--outdir', type=str,required=True,help='Path to output directory')
    parser.add_argument('--controlName', type=str,required=True,help='control group name')
    parser.add_argument('--experimentName', type=str,required=True,help='experiment group name')
    parser.add_argument('--species', type=str,default="mouse",help='species name mouse or human')
    args = parser.parse_args()
    control = args.control
    experiment = args.experiment
    outdir = args.outdir
    controlName = args.controlName
    experimentName = args.experimentName
    species = args.species
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if args.mode == "vcf":
        outfile = os.path.join(outdir,"NumberBar.csv")
        if species == "mouse":
            df = CaculateNumber(control,experiment,controlName,experimentName,outfile)
        elif species == "human":
            df = CaculateNumber(control,experiment,controlName,experimentName,outfile,geneId="ENSG")
        else:
            raise ValueError("species must be mouse or human")
        outfile = os.path.join(outdir,"NumberBar.png")
        plotNumber(df,outfile,controlName,experimentName)

        outfile = os.path.join(outdir,"TopAlterFeature.csv")
        df = TopAlterFeature(control,experiment,controlName,experimentName,5,outfile)
        outfile = os.path.join(outdir,"TopAlterFeature.png")
        plotTopfeature(df,outfile)
    elif args.mode == "annovar":
        outfile = os.path.join(outdir,"IntervalDistribution.csv")
        df = caculateIntervalDistribution(control,experiment,controlName,experimentName,outfile)
        outfile = os.path.join(outdir,"IntervalDistribution.png")
        plotIntervalDistribution(df,outfile)
    else:
        raise ValueError("mode must be vcf or annovar")



    
import pandas as pd
import glob
from pathlib import Path
def count_genes(multiano_path: str,output_path: str,col:str = "GeneName.symbol"):
    """
    统计 ANNOVAR multiano.csv 中 col 对应的基因出现次数
    """

    # 读取注释文件
    df = pd.read_csv(multiano_path)

    if col not in df.columns:
        raise ValueError(f"文件中未找到 {col} 列")

    gene_count = {}

    for genes in df[col].fillna(""):
        # 多基因用逗号分隔，例如 "TP53,RP11-34P13.7"
        for g in str(genes).split(","):
            g = g.strip()
            if g == "" or g == ".":
                continue
            gene_count[g] = gene_count.get(g, 0) + 1

    # 转成 DataFrame 并排序
    df_out = (
        pd.DataFrame.from_dict(gene_count, orient="index", columns=["count"])
        .sort_values("count", ascending=False)
    )

    df_out.to_csv(output_path)

    print(f"统计完成，共 {len(df_out)} 个基因，输出到: {output_path}")

if __name__ == "__main__":

    human_cell = ["ci8CLC","hTBLC"]
    mouse_cell = ["TLSC","ciTotiSC"]
    ### toti only
    # files = {
    #     "ci8CLC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/table/toti_only.vcf.GRCh38_multianno_with_names.csv",
    #     "hTBLC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/hTBLC/table/toti_only.vcf.GRCh38_multianno_with_names.csv",
    #     "TLSC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/TLSC/table/toti_only.vcf.GRCm39_multianno_with_names.csv",
    #     "ciTotiSC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ciTotiSC/table/toti_only.vcf.GRCm39_multianno_with_names.csv"
    # }
    # outdir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate"
    # for cell,infile in files.items():
    #     outfile = f"{outdir}/{cell}/{cell}_gene_counts.csv"
    #     count_genes(infile, outfile)
    
    ### all except toti only
    indir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect"
    infiles = glob.glob(f"{indir}/**/*_multianno_with_names.csv",recursive=True)
    target_names = ["pluripotency","totipotency"]
    for infile in infiles:
        infile = Path(infile)
        if infile.parent.parent.name in target_names:
            outfile = infile.parent.parent / f"{infile.parent.parent.parent.name}_{infile.parent.parent.name}_gene_counts.csv"
            print(outfile)
            count_genes(infile, outfile)
        else:
            print("skip:\t" + str(infile))
            continue
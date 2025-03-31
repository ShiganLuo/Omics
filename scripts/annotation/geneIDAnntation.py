import pandas as pd

def geneIDAnnotation(gtf:str,outfile:str):
    # 解析 GTF 文件（跳过注释行）
    gtf = pd.read_csv(gtf, sep="\t", comment="#", header=None, names=[
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
    ])

    # 解析 gene_id 和 gene_name
    def extract_gene_info(attr):
        gene_id, gene_name = None, None
        attributes = attr.split(";")
        for item in attributes:
            item = item.strip()
            if item.startswith("gene_id"):
                gene_id = item.split('"')[1]
            elif item.startswith("gene_name"):
                gene_name = item.split('"')[1]
            elif item.startswith("gene_type"):
                gene_type = item.split('"')[1]
        return gene_id, gene_name, gene_type

    # 应用解析函数
    gtf[["gene_id", "gene_name","gene_type"]] = gtf["attribute"].apply(lambda x: pd.Series(extract_gene_info(x)))

    # 删除无效基因（可能有些行没有 gene_id）
    gtf_filtered = gtf.dropna(subset=["gene_id"]).drop_duplicates(subset=["gene_id", "gene_name","gene_type"], keep="first")

    # 仅保留 gene_id 和 gene_name，并按原顺序输出
    gtf_filtered[["gene_id", "gene_name","gene_type"]].to_csv(outfile, sep="\t", index=False, header=True)

    print(f"去重后的基因信息已保存至 {outfile}")

if __name__ == '__main__':
    # gtf =  "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/gencode.vM36.primary_assembly.annotation.gtf"
    # outfile = "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv"
    gtf = "/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/gencode.v47.primary_assembly.annotation.gtf"
    outfile = "/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/geneIDAnnotation.csv"
    geneIDAnnotation(gtf,outfile)

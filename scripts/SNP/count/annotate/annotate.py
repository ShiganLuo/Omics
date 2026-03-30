import pandas as pd
import pickle
import os
from typing import Dict
import glob
from pathlib import Path


def get_file_signature(file_path: str) -> Dict:
    """返回 GTF 文件签名，用于判断缓存是否失效"""
    stat = os.stat(file_path)
    return {
        "path": os.path.abspath(file_path),
        "size": stat.st_size,
        "mtime": stat.st_mtime,
    }


def save_cache(cache_path: str, gene_map: dict, signature: dict):
    """保存缓存：包括基因映射和 GTF 签名"""
    with open(cache_path, "wb") as f:
        pickle.dump({"signature": signature, "gene_map": gene_map}, f)


def load_cache_if_valid(cache_path: str, gtf_signature: dict):
    """如果缓存存在且与当前 GTF 一致，则加载，否则返回 None"""
    if not os.path.exists(cache_path):
        return None

    try:
        with open(cache_path, "rb") as f:
            cache = pickle.load(f)
    except:
        return None  # 缓存损坏

    cached_sig = cache.get("signature", {})
    if cached_sig == gtf_signature:
        print(f"缓存匹配，直接加载：{cache_path}")
        return cache["gene_map"]

    print("检测到缓存与当前 GTF 不一致 → 忽略缓存并重建")
    return None


def parse_gtf_gene_map(gtf_path: str) -> dict:
    """从 GTF 文件解析 gene_id → gene_name"""
    gene_map = {}

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue

            info = fields[8]
            attrs = {}

            for kv in info.split(';'):
                kv = kv.strip()
                if not kv:
                    continue
                parts = kv.replace('"', '').split(' ')
                if len(parts) == 2:
                    key, val = parts
                    attrs[key] = val

            gid = attrs.get("gene_id")
            gname = attrs.get("gene_name")

            if gid and gname:
                gene_map[gid] = gname

    print(f"GTF 解析完成，共 {len(gene_map)} 个基因")
    return gene_map


def load_gtf_gene_map(gtf_path: str, cache_path="gene_map.pkl") -> dict:
    """
    智能加载 gene_map：如果缓存对应同一个 GTF，则直接使用；
    如果缓存无效，则重新解析并保存。
    """
    gtf_signature = get_file_signature(gtf_path)

    # 尝试加载缓存
    gene_map = load_cache_if_valid(cache_path, gtf_signature)

    if gene_map is not None:
        return gene_map

    # 解析 GTF
    print("开始解析 GTF（可能耗时较长）...")
    gene_map = parse_gtf_gene_map(gtf_path)

    # 保存缓存
    save_cache(cache_path, gene_map, gtf_signature)
    print(f"缓存已写入：{cache_path}")

    return gene_map


def translate_gene_ids(df:pd.DataFrame, gene_map: dict,col:str):
    """将 Gene.refGene 列从 gene_id 转成 gene_name"""

    def convert(gene_ids):
        if pd.isna(gene_ids):
            return gene_ids
        ids = gene_ids.split(',')
        names = [gene_map.get(gid, gid) for gid in ids]
        return ",".join(names)

    return df[col].apply(convert)


def convert_annovar_gene_ids(multiano_path, gtf_path,
                             cache_path="gene_map.pkl",
                             save_path=None):
    df = pd.read_csv(multiano_path)
    print("读取 multiano.csv 行数:", len(df))

    gene_map = load_gtf_gene_map(gtf_path, cache_path)

    df["GeneName.symbol"] = translate_gene_ids(df, gene_map,"Gene.refGene")

    if save_path:
        df.to_csv(save_path, index=False)
        print("结果已保存：", save_path)

    return df

def extract_gene_name_or_keep(df: pd.DataFrame,pattern:str) -> pd.DataFrame:
    """
    根据特定模式检查 'gene_name' 列。
    如果匹配，提取第一个冒号之前的内容；如果不匹配，保留原样。
    新列替换 'gene_name' 列。
    """
    df['gene_name'] = df['gene_name'].str.replace(
        pat=pattern,
        repl=r'\1',   # 替换模式：r'\1' 代表正则表达式中第一个捕获组的内容
        regex=True
    )
    
    return df
def convert_TEtranscripts_gene_ids(
        TEtranscripts_path:str,
        gtf_path:str,
        cache_path:str="gene_map.pkl",
        save_path:str=None):
    df = pd.read_csv(TEtranscripts_path,sep="\t")
    df["gene/TE"] = df["gene/TE"].astype(str).str.strip()
    df["gene/TE"] = df["gene/TE"].str.replace(
        pat=r'^\"(.*)\"$',
        repl=r'\1',
        regex=True
    )
    print("读取 TEtranscripts file 行数:", len(df))
    gene_map = load_gtf_gene_map(gtf_path, cache_path)
    df['gene_name'] = [gene_map.get(gid, gid) for gid in df["gene/TE"]]
    if "TEcount" in TEtranscripts_path:
        pattern = r'^([^:]+):[^:]+:[^:]+$'
        df = extract_gene_name_or_keep(df,pattern)
    elif "TElocal" in TEtranscripts_path:
        pattern = r'^([^:]+):[^:]+:[^:]+:[^:]+$'
        df = extract_gene_name_or_keep(df,pattern)
    else:
        raise ValueError("目前只支持TEtranscripts和TElocal的定量输出文件")
    cols = df.columns.tolist()
    cols.remove('gene/TE')
    name_index = cols.index('gene_name')
    cols.pop(name_index)
    cols.insert(0, 'gene_name')
    df_reordered = df[cols]
    if save_path is not None:
        df_reordered.to_csv(save_path,sep="\t",index=False)
    return df_reordered


def convert_matrix_gene_ids(
    matrx_path:str,
    gtf_path:str,
    cache_path:str="gene_map.pkl",
    save_path:str=None
):
    """
    将基因表达矩阵的基因id转换为gene name
    """
    df = pd.read_csv(matrx_path,sep="\t",index_col=0)
    gene_map = load_gtf_gene_map(gtf_path, cache_path)
    df.index = [gene_map.get(gid, gid) for gid in df.index]
    if save_path is not None:
        df.to_csv(save_path,sep="\t")

if __name__ == "__main__":
    human_gtf = "/disk5/luosg/Reference/GENCODE/human/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf"
    mouse_gtf = "/disk5/luosg/Reference/GENCODE/mouse/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf"
    
    ### TEtranscripts and TElocal
    # convert_TEtranscripts_gene_ids("/disk5/luosg/Totipotent20251031/output/counts/TEcount/human/all_TEcount.cntTable",
    #                                human_gtf,
    #                                "GRCh38_gene_map.pkl",
    #                                "/disk5/luosg/Totipotent20251031/output/counts/TEcount/human/all_TEcount_name.cntTable")
    # convert_TEtranscripts_gene_ids("/disk5/luosg/Totipotent20251031/output/counts/TElocal/human/all_TElocal.cntTable",
    #                             human_gtf,
    #                             "GRCh38_gene_map.pkl",
    #                             "/disk5/luosg/Totipotent20251031/output/counts/TElocal/human/all_TElocal_name.cntTable")
    # convert_TEtranscripts_gene_ids("/disk5/luosg/Totipotent20251031/output/counts/TEcount/mouse/all_TEcount.cntTable",
    #                         mouse_gtf,
    #                         "GRCm39_gene_map.pkl",
    #                         "/disk5/luosg/Totipotent20251031/output/counts/TEcount/mouse/all_TEcount_name.cntTable")
    # convert_TEtranscripts_gene_ids("/disk5/luosg/Totipotent20251031/output/counts/TElocal/mouse/all_TElocal.cntTable",
    #                     mouse_gtf,
    #                     "GRCm39_gene_map.pkl",
    #                     "/disk5/luosg/Totipotent20251031/output/counts/TElocal/mouse/all_TElocal_name.cntTable")
    convert_matrix_gene_ids("/disk5/luosg/Totipotent20251031/output/result/ci8CLC/DESeq2/TEcount_Gene.tsv",
                            gtf_path=human_gtf,
                            cache_path="GRCh38_gene_map.pkl",
                            save_path="/disk5/luosg/Totipotent20251031/output/result/ci8CLC/DESeq2/TEcount_Gene_name.tsv")
    
    convert_matrix_gene_ids("/disk5/luosg/Totipotent20251031/output/result/hTBLC/DESeq2/TEcount_Gene.tsv",
                        gtf_path=human_gtf,
                        cache_path="GRCh38_gene_map.pkl",
                        save_path="/disk5/luosg/Totipotent20251031/output/result/hTBLC/DESeq2/TEcount_Gene_name.tsv")
    
    convert_matrix_gene_ids("/disk5/luosg/Totipotent20251031/output/result/TLSC/DESeq2/TEcount_Gene.tsv",
                        gtf_path=mouse_gtf,
                        cache_path="GRCm39_gene_map.pkl",
                        save_path="/disk5/luosg/Totipotent20251031/output/result/TLSC/DESeq2/TEcount_Gene_name.tsv")
    
    convert_matrix_gene_ids("/disk5/luosg/Totipotent20251031/output/result/ciTotiSC/DESeq2/TEcount_Gene.tsv",
                    gtf_path=mouse_gtf,
                    cache_path="GRCm39_gene_map.pkl",
                    save_path="/disk5/luosg/Totipotent20251031/output/result/ciTotiSC/DESeq2/TEcount_Gene_name.tsv")
    #### totionly
    # files = {
    #     "ci8CLC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/table/toti_only.vcf.GRCh38_multianno.csv",
    #     "hTBLC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/hTBLC/table/toti_only.vcf.GRCh38_multianno.csv",
    #     "TLSC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/TLSC/table/toti_only.vcf.GRCm39_multianno.csv",
    #     "ciTotiSC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ciTotiSC/table/toti_only.vcf.GRCm39_multianno.csv"
    # }

    # human_cell = ["ci8CLC","hTBLC"]
    # mouse_cell = ["TLSC","ciTotiSC"]
    # outdir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate"
    # for cell,infile in files.items():
    #     if cell in human_cell:
    #         outfile = f"{outdir}/{cell}/toti_only.vcf.GRCh38_multianno_with_names.csv"
    #         convert_annovar_gene_ids(
    #             multiano_path=infile,
    #             gtf_path=human_gtf,
    #             cache_path="GRCh38_gene_map.pkl",
    #             save_path=outfile
    #         )
    #     elif cell in mouse_cell:
    #         outfile = f"{outdir}/{cell}/toti_only.vcf.GRCm39_multianno_with_names.csv"
    #         convert_annovar_gene_ids(
    #             multiano_path=infile,
    #             gtf_path=mouse_gtf,
    #             cache_path="GRCm39_gene_map.pkl",
    #             save_path=outfile
    #         )
    #     else:
    #         raise ValueError("not support cell")
    ### 
    # mouse = ["TLSC","ciTotiSC"]
    # human = ["ci8CLC","hTBLC"]
    # indir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect"
    # infiles = glob.glob(f"{indir}/**/*_multianno.csv",recursive=True)
    # target_names = ["pluripotency","totipotency"]
    # for infile in infiles:
    #     infile = Path(infile)
    #     outdir = infile.parent
    #     if infile.parent.parent.name in target_names:
    #         basePrefix = infile.with_suffix('')
    #         outfile = f"{str(basePrefix)}_with_names.csv" 
    #         print(infile)
    #         print(outfile)
    #         if infile.parent.parent.parent.name in human:
    #             print("human\n")
    #             convert_annovar_gene_ids(
    #                 multiano_path=infile,
    #                 gtf_path=human_gtf,
    #                 cache_path="GRCh38_gene_map.pkl",
    #                 save_path=outfile
    #             )
    #         elif infile.parent.parent.parent.name in mouse:
    #             print("mouse\n")
    #             convert_annovar_gene_ids(
    #                 multiano_path=infile,
    #                 gtf_path=mouse_gtf,
    #                 cache_path="GRCm39_gene_map.pkl",
    #                 save_path=outfile
    #             )
    #         else:
    #             raise ValueError("something wrong!")
    #     else:
    #         print("skip: " + str(infile))


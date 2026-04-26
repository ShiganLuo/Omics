import argparse
import pandas as pd
import pickle
import os
from typing import Dict
import logging
import sys
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,  # 指定输出到 stdout 而不是 stderr
	datefmt='%Y-%m-%d %H:%M:%S'
)


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
        logging.info(f"缓存匹配，直接加载：{cache_path}")
        return cache["gene_map"]

    logging.info("检测到缓存与当前 GTF 不一致 → 忽略缓存并重建")
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

    logging.info(f"GTF 解析完成，共 {len(gene_map)} 个基因")
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
    logging.info("开始解析 GTF（可能耗时较长）...")
    gene_map = parse_gtf_gene_map(gtf_path)

    # 保存缓存
    save_cache(cache_path, gene_map, gtf_signature)
    logging.info(f"缓存已写入：{cache_path}")

    return gene_map


def translate_gene_ids(df:pd.DataFrame, gene_map: dict,col:str):
    """
    Function: convert gene IDs in a DataFrame column to gene names using a provided mapping.
    Parameters:
        - df: Input DataFrame containing gene IDs.
        - gene_map: A dictionary mapping gene IDs to gene names.
        - col: The name of the column in df that contains the gene IDs to be converted.
    Returns:A Series with gene names corresponding to the gene IDs in the specified column.
    """

    def convert(gene_ids):
        if pd.isna(gene_ids):
            return gene_ids
        ids = gene_ids.split(',')
        names = [gene_map.get(gid, gid) for gid in ids]
        return ",".join(names)

    return df[col].apply(convert)



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
        geneId_col:str = "gene/TE",
        cache_path:str="gene_map.pkl",
        save_path:str=None
    ) -> pd.DataFrame:
    df = pd.read_csv(TEtranscripts_path,sep="\t")
    df[geneId_col] = df[geneId_col].astype(str).str.strip()
    df[geneId_col] = df[geneId_col].str.replace(
        pat=r'^\"(.*)\"$',
        repl=r'\1',
        regex=True
    )
    logging.info(f"读取 TEtranscripts file 行数:{len(df)}")
    gene_map = load_gtf_gene_map(gtf_path, cache_path)
    df['gene_name'] = translate_gene_ids(df,gene_map,geneId_col)
    if "TEcount" in TEtranscripts_path:
        pattern = r'^([^:]+):[^:]+:[^:]+$'
        df = extract_gene_name_or_keep(df,pattern)
    elif "TElocal" in TEtranscripts_path:
        pattern = r'^([^:]+):[^:]+:[^:]+:[^:]+$'
        df = extract_gene_name_or_keep(df,pattern)
    else:
        raise ValueError("目前只支持TEtranscripts和TElocal的定量输出文件")
    cols = df.columns.tolist()
    cols.remove(geneId_col)
    name_index = cols.index('gene_name')
    cols.pop(name_index)
    cols.insert(0, 'gene_name')
    df_reordered = df[cols]
    if save_path is not None:
        df_reordered.to_csv(save_path,sep="\t",index=False)
    return df_reordered



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Gene IDs to gene names and aggregate expression values")
    parser.add_argument("-c", "--count", help="featureCounts output file path")
    parser.add_argument("-g", "--gtf", help="Path to the GTF annotation file")
    parser.add_argument("-o", "--output", help="Output file path")
    args = parser.parse_args()
    cache_path = f"{os.path.dirname(args.gtf)}/gene_map.pkl"
    gene_id2name_df = convert_TEtranscripts_gene_ids(
                                                    TEtranscripts_path=args.count, 
                                                     gtf_path=args.gtf,
                                                     cache_path=cache_path,
                                                     save_path=args.output)



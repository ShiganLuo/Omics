import pandas as pd
import os
import pickle
import logging
from typing import Dict
logger = logging.getLogger(__name__)
def get_file_signature(file_path: str) -> Dict:
    """返回 GTF 文件签名，用于判断缓存是否失效"""
    stat = os.stat(file_path)
    return {
        "path": os.path.abspath(file_path),
        "size": stat.st_size,
        "mtime": stat.st_mtime,
    }


def save_cache(cache_path: str, gene_map: Dict[str, str], signature: Dict):
    """保存缓存：包括基因映射和 GTF 签名"""
    with open(cache_path, "wb") as f:
        pickle.dump({"signature": signature, "gene_map": gene_map}, f)


def load_cache_if_valid(cache_path: str, gtf_signature: Dict):
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
        logger.info(f"cache matched, loading directly: {cache_path}")
        return cache["gene_map"]

    logger.info("cache mismatch detected → ignoring cache and rebuilding")
    return None


def parse_gtf_gene_map(gtf_path: str) -> Dict:
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

    logger.info(f"GTF 解析完成，共 {len(gene_map)} 个基因")
    return gene_map


def load_gtf_gene_map(gtf_path: str, cache_path="gene_map.pkl") -> Dict[str, str]:
    """
    Function to load gene_id to gene_name mapping from GTF file, with caching mechanism.
     - gtf_path: path to the GTF file
     - cache_path: path to the cache file (default: gene_map.pkl in the same directory as GTF)
     Returns a dictionary mapping gene_id to gene_name.
     The function first checks if a valid cache exists (matching GTF signature). If so, it loads the gene map from cache. Otherwise, it parses the GTF file to build the gene map, saves it to cache, and returns it.
    """
    gtf_signature = get_file_signature(gtf_path)

    gene_map = load_cache_if_valid(cache_path, gtf_signature)

    if gene_map is not None:
        return gene_map

    logger.info("begin parsing GTF to build gene map...")
    gene_map = parse_gtf_gene_map(gtf_path)

    cache_path = os.path.join(os.path.dirname(gtf_path), "gene_map.pkl")
    save_cache(cache_path, gene_map, gtf_signature)
    logger.info(f"cache written: {cache_path}")

    return gene_map


def translate_gene_ids(df:pd.DataFrame, gene_map: Dict[str,str],col:str):
    """将 Gene.refGene 列从 gene_id 转成 gene_name"""

    def convert(gene_ids):
        if pd.isna(gene_ids):
            return gene_ids
        ids = gene_ids.split(',')
        names = [gene_map.get(gid, gid) for gid in ids]
        return ",".join(names)

    return df[col].apply(convert)



def convert_gene_id_to_name(
    infile:str,
    gtf_path:str,
    gene_id_col:str="name",
    ):
    if not os.path.exists(infile):
        logger.error(f"输入文件不存在：{infile}")
        return None
    if not os.path.exists(gtf_path):
        logger.error(f"GTF 文件不存在：{gtf_path}")
        return None
    gtf_cache_path = os.path.join(os.path.dirname(infile), "gene_map.pkl")
    gene_map = load_gtf_gene_map(gtf_path, gtf_cache_path)
    df = pd.read_csv(infile, sep="\t")
    df["gene_name"] = translate_gene_ids(df, gene_map, gene_id_col)
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="translate gene_id to gene_name in a TSV file based on GTF annotation")
    parser.add_argument("--infile", required=True, help="input TSV file path (must contain a column with gene IDs)")
    parser.add_argument("--gtf", required=True, help="GTF file path")
    parser.add_argument("--outfile", required=True, help="output file path (TSV format)")
    parser.add_argument("--gene_id_col", default="name", help="column name for gene IDs (default: 'name')")
    args = parser.parse_args()

    result_df = convert_gene_id_to_name(args.infile, args.gtf, args.gene_id_col)
    if result_df is not None:
        result_df.to_csv(args.outfile, sep="\t", index=False)
        logger.info(f"转换完成，结果已保存：{args.outfile}")

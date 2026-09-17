"""Gene ID to gene name translation utilities.

Provides functions for mapping Ensembl gene IDs (e.g. ``ENSMUSG00000051951``)
to human-readable gene symbols (e.g. ``Xkr4``) using GENCODE/Ensembl GTF
annotations.  A pickle-based cache avoids redundant GTF parsing.

Typical usage::

    from src.annotation.gene_id2name import convert_featurecounts_gene_ids

    df = convert_featurecounts_gene_ids(
        count="featureCounts_output.txt",
        gtf_path="/path/to/gencode.gtf",
        save_path="counts_by_gene_name.tsv",
    )

Supported input formats:
- featureCounts output (Ensembl gene IDs in ``Geneid`` column)
- DESeq2 result tables (Ensembl gene IDs as row index)
- TEtranscripts/TElocal quantification (mixed gene and TE identifiers)
- ANNOVAR multianno CSV
"""

import logging
import os
import pickle
import sys
from typing import Dict, Optional, Union

import pandas as pd

logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,  # Output to stdout instead of stderr
	datefmt='%Y-%m-%d %H:%M:%S'
)


def get_file_signature(file_path: str) -> Dict:
    """Compute a file signature for cache invalidation.

    The signature captures the absolute path, file size, and last modification
    time. Two signatures can be compared with ``==`` to detect whether the
    underlying file has changed since the cache was created.

    Parameters
    ----------
    file_path : str
        Path to the file whose signature is to be computed.

    Returns
    -------
    dict
        File signature with keys ``"path"`` (str), ``"size"`` (int), and
        ``"mtime"`` (float).
    """
    stat = os.stat(file_path)
    return {
        "path": os.path.abspath(file_path),
        "size": stat.st_size,
        "mtime": stat.st_mtime,
    }


def save_cache(cache_path: str, gene_map: dict, signature: dict):
    """Persist a gene-id-to-name mapping alongside its GTF signature.

    The output pickle contains a single dict with ``"signature"`` and
    ``"gene_map"`` keys, allowing :func:`load_cache_if_valid` to verify
    freshness on the next load.

    Parameters
    ----------
    cache_path : str
        Destination path for the pickle file.
    gene_map : dict
        Mapping of gene_id strings to gene_name strings.
    signature : dict
        File signature dict produced by :func:`get_file_signature`.
    """
    with open(cache_path, "wb") as f:
        pickle.dump({"signature": signature, "gene_map": gene_map}, f)


def load_cache_if_valid(cache_path: str, gtf_signature: dict):
    """Load a cached gene mapping if the pickle exists and its GTF signature matches.

    Parameters
    ----------
    cache_path : str
        Path to the pickle cache file.
    gtf_signature : dict
        Current GTF file signature produced by :func:`get_file_signature`.

    Returns
    -------
    dict or None
        The cached gene_id -> gene_name mapping, or ``None`` when the cache
        is missing, corrupted, or stale (GTF changed).
    """
    if not os.path.exists(cache_path):
        return None

    try:
        with open(cache_path, "rb") as f:
            cache = pickle.load(f)
    except:
        return None  # Corrupted cache

    cached_sig = cache.get("signature", {})
    if cached_sig == gtf_signature:
        logging.info(f"Cache matched, loading directly: {cache_path}")
        return cache["gene_map"]

    logging.info("Cache does not match current GTF, ignoring cache and rebuilding")
    return None


def parse_gtf_gene_map(gtf_path: str) -> dict:
    """Build a gene_id -> gene_name mapping by parsing a GENCODE/Ensembl GTF.

    Only lines with feature type ``gene`` are examined.  For each such line the
    ``gene_id`` and ``gene_name`` attributes are extracted from the ninth column
    (GFF3-style ``key "value"`` pairs separated by semicolons).

    Parameters
    ----------
    gtf_path : str
        Path to a GTF annotation file (gzipped files are **not** supported;
        decompress first).

    Returns
    -------
    dict
        Mapping from Ensembl gene_id (e.g. ``"ENSMUSG00000051951"``) to
        gene_name/symbol (e.g. ``"Xkr4"``).

    Notes
    -----
    This is a pure-Python line-by-line parser and does not depend on pandas or
    gffutils.  Parsing a full GENCODE GTF typically takes 5-15 s.
    """
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

    logging.info(f"GTF parsing complete, {len(gene_map)} genes found")
    return gene_map


def load_gtf_gene_map(gtf_path: str, cache_path="gene_map.pkl") -> dict:
    """Load a gene_id -> gene_name mapping, using a pickle cache when valid.

    On first call (or when the GTF has changed), the full GTF is parsed via
    :func:`parse_gtf_gene_map` and the result is persisted to *cache_path*.
    Subsequent calls with the same GTF return the cached mapping immediately.

    Parameters
    ----------
    gtf_path : str
        Path to the GTF annotation file.
    cache_path : str, optional
        Path for the pickle cache file.  Defaults to ``"gene_map.pkl"`` in the
        current working directory.

    Returns
    -------
    dict
        gene_id -> gene_name mapping.
    """
    gtf_signature = get_file_signature(gtf_path)

    # Try loading from cache
    gene_map = load_cache_if_valid(cache_path, gtf_signature)

    if gene_map is not None:
        return gene_map

    # Parse GTF
    logging.info("Parsing GTF (this may take a while)...")
    gene_map = parse_gtf_gene_map(gtf_path)

    # Save cache
    save_cache(cache_path, gene_map, gtf_signature)
    logging.info(f"Cache written: {cache_path}")

    return gene_map


def translate_gene_ids(df:pd.DataFrame, gene_map: dict,col:str):
    """Translate gene IDs in a DataFrame column to gene names.

    Comma-separated ID lists are supported: each token is mapped independently
    and the results are joined back with commas.  IDs not present in
    *gene_map* are left unchanged.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a column of gene IDs.
    gene_map : dict
        Mapping of gene_id -> gene_name, typically from
        :func:`load_gtf_gene_map`.
    col : str
        Name of the column in *df* that holds the gene IDs.

    Returns
    -------
    pd.Series
        Series of translated gene names, aligned with *df*'s index.
    """

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
    """Convert ANNOVAR multianno gene IDs to gene names.

    Reads a ``multianno.csv`` produced by ANNOVAR's ``table_annovar.pl``,
    maps the ``Gene.refGene`` column from Ensembl gene IDs to gene symbols,
    and appends a ``GeneName.symbol`` column.

    Parameters
    ----------
    multiano_path : str
        Path to the ANNOVAR multianno CSV file.
    gtf_path : str
        Path to the GTF annotation file for gene_id -> gene_name mapping.
    cache_path : str, optional
        Pickle cache path (default ``"gene_map.pkl"``).
    save_path : str or None, optional
        If provided, write the augmented DataFrame to this CSV path.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with an additional ``GeneName.symbol`` column.
    """
    df = pd.read_csv(multiano_path)
    logging.info(f"Read multiano.csv rows: {len(df)}")

    gene_map = load_gtf_gene_map(gtf_path, cache_path)

    df["GeneName.symbol"] = translate_gene_ids(df, gene_map,"Gene.refGene")

    if save_path:
        df.to_csv(save_path, index=False)
        logging.info(f"Results saved: {save_path}")

    return df

def extract_gene_name_or_keep(df: pd.DataFrame,pattern:str) -> pd.DataFrame:
    """Apply a regex substitution to the ``gene_name`` column in-place.

    When *pattern* contains a capture group, the replacement ``r'\\1'`` extracts
    only the captured portion.  Rows that do not match are left unchanged.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame that must contain a ``gene_name`` column.  Modified in-place.
    pattern : str
        Regex pattern with at least one capture group, e.g.
        ``r'^([^:]+):[^:]+:[^:]+$'`` to extract the subfamily name from a
        TEcount identifier.

    Returns
    -------
    pd.DataFrame
        The same *df* with ``gene_name`` updated.
    """
    df['gene_name'] = df['gene_name'].str.replace(
        pat=pattern,
        repl=r'\1',   # Replacement pattern: r'\1' refers to the first capture group
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
    """Convert TEtranscripts/TElocal gene IDs to gene names.

    Reads a TEcount or TElocal quantification TSV, translates the gene/TE
    column via a GTF-derived mapping, and handles TE identifier formats:

    - **TEcount**: ``Class:Family:Subfamily`` → extract ``Subfamily``
    - **TElocal**: ``Class:Family:Subfamily:Locus`` → extract ``Locus``

    Gene rows (without colons) pass through unchanged.

    Parameters
    ----------
    TEtranscripts_path : str
        Path to the TEcount or TElocal output TSV file.
    gtf_path : str
        Path to the GTF annotation for gene_id -> gene_name mapping.
    geneId_col : str, optional
        Column name containing gene/TE identifiers (default ``"gene/TE"``).
    cache_path : str, optional
        Pickle cache path (default ``"gene_map.pkl"``).
    save_path : str or None, optional
        If provided, write the reordered result to this TSV path.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``gene_name`` as the first column, original ``gene/TE``
        column removed, and all other columns preserved.

    Raises
    ------
    ValueError
        If the filename does not contain ``"TEcount"`` or ``"TElocal"``.
    """
    df = pd.read_csv(TEtranscripts_path,sep="\t")
    df[geneId_col] = df[geneId_col].astype(str).str.strip()
    df[geneId_col] = df[geneId_col].str.replace(
        pat=r'^\"(.*)\"$',
        repl=r'\1',
        regex=True
    )
    logging.info(f"Read TEtranscripts file rows: {len(df)}")
    gene_map = load_gtf_gene_map(gtf_path, cache_path)
    df['gene_name'] = translate_gene_ids(df,gene_map,geneId_col)
    if "TEcount" in TEtranscripts_path:
        pattern = r'^([^:]+):[^:]+:[^:]+$'
        df = extract_gene_name_or_keep(df,pattern)
    elif "TElocal" in TEtranscripts_path:
        pattern = r'^([^:]+):[^:]+:[^:]+:[^:]+$'
        df = extract_gene_name_or_keep(df,pattern)
    else:
        raise ValueError("Only TEtranscripts and TElocal quantification output files are supported")
    cols = df.columns.tolist()
    cols.remove(geneId_col)
    name_index = cols.index('gene_name')
    cols.pop(name_index)
    cols.insert(0, 'gene_name')
    df_reordered = df[cols]
    if save_path is not None:
        df_reordered.to_csv(save_path,sep="\t",index=False)
    return df_reordered


def convert_featurecounts_gene_ids(
    count: Union[str,pd.DataFrame],
    gtf_path: str,
    cache_path: str = "gene_map.pkl",
    save_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Convert featureCounts Geneid to gene_name and aggregate expression values
    at gene level.

    This function maps Ensembl gene IDs (Geneid) produced by featureCounts
    to gene symbols (gene_name) using a GTF annotation file, and then collapses
    multiple rows belonging to the same gene_name by summing expression values.

    Notes
    -----
    Why summation instead of averaging?

    In featureCounts output, a single gene_name may correspond to multiple
    Geneid entries due to:
        - multiple Ensembl gene versions (e.g. ENSMUSGxxxx.x)
        - multiple transcript-derived features aggregated at gene level
        - duplicated or split annotations in the GTF

    Expression values (counts, CPM, TPM, etc.) represent *abundance* or
    *read support* for genomic features. When multiple Geneid entries map
    to the same gene_name, their values should be **summed** to obtain the
    total gene-level expression.

    Averaging would artificially reduce expression levels and break the
    biological interpretation, because:
        - expression is additive, not an intensity per feature
        - downstream analyses (DESeq2, edgeR, limma, etc.) assume summed counts
        - TPM/CPM values are proportional to total transcript abundance

    Therefore, summation is the correct and standard approach for collapsing
    transcript-/ID-level data into gene-level expression matrices.

    Parameters
    ----------
    count : str or pd.DataFrame
        Path to a featureCounts output TSV file, or an already-loaded
        DataFrame.  Must contain a ``"Geneid"`` column.
    gtf_path : str
        Path to the GTF annotation file used for Geneid -> gene_name mapping.
    cache_path : str, optional
        Pickle cache path for the parsed gene ID mapping
        (default ``"gene_map.pkl"``).
    save_path : str or None, optional
        If provided, save the converted gene-level table to this path as TSV.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``gene_name`` as the first column, ``Geneid`` removed,
        and all count/sample columns preserved.

    Raises
    ------
    ValueError
        If the input DataFrame does not contain a ``"Geneid"`` column.

    See Also
    --------
    convert_DESeq2_gene_ids : Similar conversion for DESeq2 result tables.
    convert_TEtranscripts_gene_ids : Similar conversion for TEtranscripts output.
    """
    if isinstance(count, pd.DataFrame):
        df = count.copy()
    else:
        df = pd.read_csv(count, sep="\t")

    if "Geneid" not in df.columns:
        raise ValueError("Input count data must contain 'Geneid' column")

    df["Geneid"] = df["Geneid"].astype(str).str.strip()
    logging.info(f"Read count rows: {len(df)}")

    gene_map = load_gtf_gene_map(gtf_path, cache_path)
    df["gene_name"] = translate_gene_ids(df, gene_map, "Geneid")

    cols = df.columns.tolist()
    cols.remove("gene_name")
    cols.remove("Geneid")
    df = df[["gene_name"] + cols]
    # df = (
    #     df
    #     .groupby('gene_name', as_index=False)
    #     .sum(numeric_only=True)
    # )
    if save_path is not None:
        df.to_csv(save_path, sep="\t", index=False)

    return df

def convert_DESeq2_gene_ids(
    deseq2_path: str,
    save_path: str,
    gtf_path: str,
    cache_path: str = "gene_map.pkl",
):
    """
    Convert DESeq2 result table gene IDs to gene names.

    Reads a DESeq2 output TSV (with row names as Geneid), translates Ensembl
    gene IDs to gene symbols, and reorders columns so ``gene_name`` comes first.

    Parameters
    ----------
    deseq2_path : str
        Path to the DESeq2 results TSV file.  The first column (index) is
        treated as the Geneid.
    save_path : str
        Output path.  Format is inferred from the extension: ``.tsv`` writes
        tab-separated, ``.xlsx`` writes Excel, anything else writes
        comma-separated CSV.
    gtf_path : str
        Path to the GTF annotation file used for Geneid -> gene_name mapping.
    cache_path : str, optional
        Pickle cache path for the parsed gene ID mapping
        (default ``"gene_map.pkl"``).

    Returns
    -------
    pd.DataFrame
        DataFrame with ``gene_name`` as the first column, ``Geneid`` removed,
        and all DESeq2 result columns (baseMean, log2FC, pvalue, padj, etc.)
        preserved.

    See Also
    --------
    convert_featurecounts_gene_ids : Similar conversion for featureCounts output.
    """
    df = pd.read_csv(deseq2_path, index_col=0,sep="\t")
    df.reset_index(inplace=True,names="Geneid")
    logging.info(f"Read DESeq2 rows: {len(df)}")
    gene_map = load_gtf_gene_map(gtf_path, cache_path)
    df["gene_name"] = translate_gene_ids(df, gene_map, "Geneid")

    cols = df.columns.tolist()
    cols.remove("gene_name")
    cols.remove("Geneid")
    df = df[["gene_name"] + cols]

    if save_path is not None:
        if save_path.endswith(".tsv"):
            df.to_csv(save_path, sep="\t", index=False)
        elif save_path.endswith(".xlsx"):
            df.to_excel(save_path, index=False)
        else:        
            df.to_csv(save_path, index=False)

    return df

if __name__ == "__main__":
    human_gtf = "/home/luosg/Database/Reference/human/GENCODE/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf"
    human_map = "/home/luosg/Database/Reference/human/GENCODE/GRCh38/gene_map.pkl"
    mouse_gtf = "/home/luosg/Database/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf"
    mouse_map = "/home/luosg/Database/Reference/mouse/GENCODE/GRCm39/gene_map.pkl"
    convert_DESeq2_gene_ids(
        deseq2_path="/home/luosg/Data/genomeStability/output/RNAseq/diff_expression/Scramble_vs_Rn7sk/TEcount_Gene.tsv",
        save_path="/home/luosg/Data/genomeStability/output/RNAseq/diff_expression/Scramble_vs_Rn7sk/TEcount_Gene_name.xlsx",
        gtf_path=mouse_gtf,
        cache_path=mouse_map
   )

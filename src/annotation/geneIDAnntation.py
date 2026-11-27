import pandas as pd
import logging
import sys
from pathlib import Path
from typing import Union
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,
	datefmt='%Y-%m-%d %H:%M:%S'
)


def geneIDAnnotation(gtf_path: Union[str, Path]) -> pd.DataFrame:
	gtf = pd.read_csv(
		gtf_path,
		sep="	",
		comment="#",
		header=None,
		names=["seqname", "_b", "_c", "_d", "_e", "_f", "_g", "_h", "attribute"],
		usecols=[8],
	)

	gtf['gene_id'] = gtf['attribute'].str.extract(r'gene_id\s*"(.*?)"')
	gtf['gene_name'] = gtf['attribute'].str.extract(r'gene_name\s*"(.*?)"')
	gtf['gene_type'] = gtf['attribute'].str.extract(r'gene_type\s*"(.*?)"')
	if gtf['gene_type'].isnull().all():
		gtf['gene_type'] = gtf['attribute'].str.extract(r'gene_biotype\s*"(.*?)"')
	gtf_filtered = gtf.dropna(subset=["gene_id"])
	gtf_filtered = gtf_filtered.drop_duplicates(
		subset=["gene_id", "gene_name", "gene_type"],
		keep="first"
	)
	na_mask = gtf_filtered['gene_name'].isnull() | gtf_filtered['gene_name'].str.strip().eq('')
	gtf_filtered.loc[na_mask, 'gene_name'] = gtf_filtered.loc[na_mask, 'gene_id']
	return gtf_filtered[["gene_id", "gene_name", "gene_type"]]


def teBed(gtf_path: Union[str, Path]) -> pd.DataFrame:
	"""Extract per-locus coordinates and TE classification from a RepeatMasker-derived TE GTF.

	The TE GTF (e.g. rheMac10_rmsk_TE.gtf) contains one row per exon segment with attributes:
		gene_id "<subfamily_id>"; transcript_id "<instance_id>";
		family_id "<family>"; class_id "<class>";

	Returns a DataFrame with one row per segment line (1-based inclusive coords):
		chrom, start, end, strand, gene_id (= subfamily_id),
		family_id, class_id

	Rows lacking gene_id (e.g. comment/region lines, non-TE features) are removed.
	"""
	gtf = pd.read_csv(
		gtf_path,
		sep="	",
		comment="#",
		header=None,
		names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
	)

	gtf["gene_id"] = gtf["attribute"].str.extract(r'gene_id\s*"(.*?)"')
	gtf = gtf.dropna(subset=["gene_id"])
	gtf["family_id"] = gtf["attribute"].str.extract(r'family_id\s*"(.*?)"')
	gtf["class_id"] = gtf["attribute"].str.extract(r'class_id\s*"(.*?)"')

	gtf_filtered = gtf.copy()
	gtf_filtered["start"] = gtf_filtered["start"] - 1
	return gtf_filtered[["chrom", "start", "end", "strand", "gene_id", "family_id", "class_id"]]



if __name__ == '__main__':
	te_gtf = "/home/luosg/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf"
	te_outfile = "/home/luosg/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.bed"
	te_df = teBed(te_gtf)
	te_df.to_csv(te_outfile, sep="\t", index=False, header=True)
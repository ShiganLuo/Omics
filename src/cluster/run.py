import re

import pandas as pd

from expression_cluster import run_clustering


# ── Species config ──────────────────────────────────────────────────────────
SPECIES = {
    "GRCh38": {
        "EED": "/data/pub/zhousha/Totipotent20251031/data/EED/EED_GRCh38_featureCounts.tsv",
        "p2t": "/data/pub/zhousha/Totipotent20251031/output/pluripotency2totipotency/RNAseq/counts/featureCounts/GRCh38/GRCh38_featureCounts.tsv",
    },
    "GRCm39": {
        "EED": "/data/pub/zhousha/Totipotent20251031/data/EED/EED_GRCm39_featureCounts.tsv",
        "p2t": "/data/pub/zhousha/Totipotent20251031/output/pluripotency2totipotency/RNAseq/counts/featureCounts/GRCm39/GRCm39_featureCounts.tsv",
    },
}
OUTPUT_DIR = "/data/pub/zhousha/Totipotent20251031/output/pluripotency2totipotency/RNAseq/cluster/featureCounts"

# Early embryo patterns to exclude (keep stem cell lines)
EMBRYO_PATTERNS = re.compile(
    r"Oocyte|Zygote|\d+-cell|Morula|blastocyst", re.IGNORECASE
)


def infer_cell_type(sample_name: str) -> str:
    """Strip lab prefix and replicate suffix to get cell type.

    Examples: WIBR3_E8_hESC1 -> hESC, ci8CLC_r1 -> ci8CLC, mESC1 -> mESC
    """
    s = sample_name
    # strip trailing replicate suffix
    s = re.sub(r"_r\d+$", "", s)
    s = re.sub(r"-\d+$", "", s)
    s = re.sub(r"_\d+$", "", s)
    s = re.sub(r"\d+$", "", s)
    # strip known lab prefixes (WIBR3_E8_, etc.)
    s = re.sub(r"^[A-Z]+\d*_[A-Z]\d*_", "", s)
    return s


def load_featurecounts(path: str) -> pd.DataFrame:
    """Load featureCounts output, skip comment/metadata columns, clean sample names."""
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    df = df.drop(columns=["Chr", "Start", "End", "Strand", "Length"], errors="ignore")
    df.columns = [re.sub(r".*/([^/]+)\.[^.]+$", r"\1", c) if "/" in c else c for c in df.columns]
    return df


def load_species_counts(genome: str) -> pd.DataFrame:
    """Load and merge EED + p2t featureCounts for one species, keep only embryo samples."""
    paths = SPECIES[genome]
    frames = []
    for label, path in paths.items():
        df = load_featurecounts(path)
        print(f"  {label} ({genome}): {df.shape[0]} genes x {df.shape[1]} samples")
        frames.append(df)
    merged = pd.concat(frames, axis=1)
    merged = merged.groupby(level=0).sum()

    # Filter: keep stem cell lines, exclude early embryo stages
    stem_cols = [c for c in merged.columns if not EMBRYO_PATTERNS.search(c)]
    excluded = [c for c in merged.columns if EMBRYO_PATTERNS.search(c)]
    merged = merged[stem_cols]
    print(f"  Combined ({genome}): {merged.shape[0]} genes x {merged.shape[1]} stem cell samples")
    if excluded:
        print(f"  Excluded embryo: {excluded}")
    return merged


if __name__ == "__main__":
    for genome in SPECIES:
        print(f"Processing {genome}")

        df = load_species_counts(genome)
        prefix = f"{OUTPUT_DIR}/{genome}_stem"
        matrix_path = f"{prefix}_counts.tsv"
        df.to_csv(matrix_path, sep="\t")
        print(f"Saved to {matrix_path}")

        cell_types = {col: infer_cell_type(col) for col in df.columns}
        print(f"Cell types: {sorted(set(cell_types.values()))}")

        run_clustering(
            matrix=matrix_path,
            output_prefix=prefix,
            method="hierarchical",
            n_clusters=5,
            min_mean=1.0,
            log_transform=True,
            zscore=True,
            transpose=True,
            metric="cosine",
            linkage_method="complete",
            heatmap=True,
            dendrogram=True,
            pca=True,
            cell_types=cell_types,
        )

"""Convert Cell Ranger output to h5ad format for Scanpy downstream analysis.

Reads the filtered_feature_bc_matrix from a Cell Ranger count run and
produces a standard AnnData .h5ad file with:
  - obs: cell barcodes
  - var: gene symbols / feature IDs
  - X: count matrix (sparse)
"""
import argparse
import anndata as ad
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(
        description="Convert Cell Ranger filtered matrix to h5ad"
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to filtered_feature_bc_matrix directory (barcodes/features/matrix tsv.gz)"
    )
    parser.add_argument("--output", required=True, help="Output .h5ad path")
    parser.add_argument("--sample-id", default="", help="Sample ID for .obs annotation")
    args = parser.parse_args()

    adata = sc.read_10x_mtx(args.input, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()
    if args.sample_id:
        adata.obs["sample_id"] = args.sample_id
    adata.write_h5ad(args.output)
    print(f"Wrote {adata.n_obs} cells x {adata.n_vars} genes to {args.output}")


if __name__ == "__main__":
    main()

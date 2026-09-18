import pandas as pd

def combine_tpm(

):
    """
    Combine TPM (Transcripts Per Million) values from multiple samples into a single DataFrame.

    Returns:
        pd.DataFrame: A DataFrame containing combined TPM values with genes as rows and samples as columns.
    """
    EED_GRCh38 = "/data/pub/zhousha/Totipotent20251031/data/EED/EED_GRCh38_tpm.tsv"
    EED_GRCm39 = "/data/pub/zhousha/Totipotent20251031/data/EED/EED_GRCm39_tpm.tsv"
    p2t_GRCh38 = "/data/pub/zhousha/Totipotent20251031/output/pluripotency2totipotency/RNAseq/function/GRCh38/gsva/expression_tpm.tsv"
    p2t_GRCm39 = "/data/pub/zhousha/Totipotent20251031/output/pluripotency2totipotency/RNAseq/function/GRCm39/gsva/expression_tpm.tsv"
    df_EED_GRCh38 = pd.read_csv(EED_GRCh38, sep="\t", index_col=0).groupby(level=0).sum()
    df_EED_GRCm39 = pd.read_csv(EED_GRCm39, sep="\t", index_col=0).groupby(level=0).sum()
    df_p2t_GRCh38 = pd.read_csv(p2t_GRCh38, sep="\t", index_col=0).groupby(level=0).sum()
    df_p2t_GRCm39 = pd.read_csv(p2t_GRCm39, sep="\t", index_col=0).groupby(level=0).sum()
    df_combine = pd.concat([df_EED_GRCh38, df_EED_GRCm39, df_p2t_GRCh38, df_p2t_GRCm39], axis=1)
    print(df_combine.columns)
    return df_combine

if __name__ == "__main__":
    combine_tpm()
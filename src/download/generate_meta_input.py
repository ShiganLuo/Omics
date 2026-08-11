import os
import pandas as pd
from typing import List, Optional

# gsm_metadata 中的标识列，不参与 design 构建
_ID_COLS = {"GSM", "GSE", "SRA", "BioSample"}


def generate_meta_input(
        gsm_meta_path: str,
        sra_meta_path: str,
        fastq_dir: str,
        output_path: str,
        sample_id_col: str = "Sample_id",
        data_id_col: str = "Data_id",
        layout_col: str = "Layout",
        gsm_key: str = "GSM",
        gsm_sra_col: str = "SRA",
        sra_srx_col: str = "SRX",
        design_sep: str = "_",
    ) -> None:
    """
    Function:
        合并 gsm_metadata.csv 和 sra_metadata.csv，生成 meta_input.tsv。
        design 列由 gsm_metadata 中所有非标识列（排除 GSM/GSE/SRA/BioSample）的非空值
        用 "_" 拼接作为组名，同组内按 GSM 出序分配 .repN。
    Parameters:
        gsm_meta_path: gsm_metadata.csv 路径
        sra_meta_path: sra_metadata.csv 路径
        fastq_dir: fastq 文件目录
        output_path: 输出 TSV 路径
    """
    gsm = pd.read_csv(gsm_meta_path)
    sra = pd.read_csv(sra_meta_path)

    # merge: gsm.SRA == sra.SRX
    merged = gsm.merge(sra, left_on=gsm_sra_col, right_on=sra_srx_col, how="inner")

    # --- build design group from all non-id columns in gsm_metadata ---
    design_cols = [c for c in gsm.columns if c not in _ID_COLS]

    def build_group(row):
        parts = []
        for col in design_cols:
            val = str(row[col]).strip()
            if val and val.lower() != "nan":
                parts.append(val)
        return design_sep.join(parts) if parts else "unclassified"

    merged["__group"] = merged.apply(build_group, axis=1)

    # assign rep number within each group (by GSM, technical reps share rep number)
    rep_map = {}
    for group, sub in merged.groupby("__group"):
        seen = {}
        idx = 0
        for _, row in sub.iterrows():
            gsm_val = row[gsm_key]
            if gsm_val not in seen:
                idx += 1
                seen[gsm_val] = idx
            rep_map[(group, gsm_val)] = idx

    merged["design"] = merged.apply(
        lambda r: f'{r["__group"]}.rep{rep_map[(r["__group"], r[gsm_key])]}', axis=1
    )

    # --- build fastq paths + output ---
    records = []
    for _, row in merged.iterrows():
        data_id = row[data_id_col]
        sample_id = row[sample_id_col]
        layout = str(row[layout_col]).upper() if layout_col in row else "PAIRED"
        fq1 = os.path.join(fastq_dir, f"{data_id}_1.fastq.gz")
        fq2 = os.path.join(fastq_dir, f"{data_id}_2.fastq.gz")
        fq_single = os.path.join(fastq_dir, f"{data_id}.fastq.gz")
        design = row["design"]
        if layout == "PAIRED":
            if os.path.exists(fq1) and os.path.exists(fq2):
                records.append([sample_id, data_id, fq1, fq2, design])
            elif os.path.exists(fq_single):
                records.append([sample_id, data_id, fq_single, "", design])
        else:
            if os.path.exists(fq_single):
                records.append([sample_id, data_id, fq_single, "", design])
            elif os.path.exists(fq1):
                records.append([sample_id, data_id, fq1, "", design])

    out_df = pd.DataFrame(records, columns=["sample_id", "data_id", "fastq_1", "fastq_2", "design"])
    out_df.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Generate meta_input.tsv from gsm_metadata + sra_metadata + fastq.")
    p.add_argument("-g", "--gsm_meta_path", required=True, help="Path to gsm_metadata.csv")
    p.add_argument("-i", "--sra_meta_path", required=True, help="Path to sra_metadata.csv")
    p.add_argument("-d", "--fastq_dir", required=True, help="Directory containing FASTQ files")
    p.add_argument("-o", "--output_path", required=True, help="Output TSV path")
    p.add_argument("--design_sep", default="_", help="Separator for design group name (default: _)")
    args = p.parse_args()
    generate_meta_input(
        gsm_meta_path=args.gsm_meta_path,
        sra_meta_path=args.sra_meta_path,
        fastq_dir=args.fastq_dir,
        output_path=args.output_path,
        design_sep=args.design_sep,
    )

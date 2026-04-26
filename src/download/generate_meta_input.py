import os
import pandas as pd
from pathlib import Path

def generate_meta_input(
        sra_meta_path:str, 
        fastq_dir:str, 
        output_path:str,
        sample_id_col:str = 'Sample_id',
        data_id_col:str = 'Data_id',
        Layout_col:str = 'Layout'
    ) -> None:
    """
    Function: 
        Generate a meta input file for downstream analysis based on SRA metadata and available FASTQ files.
    Parameters:
        sra_meta_path (str): Path to the SRA metadata CSV file.
        fastq_dir (str): Directory containing the FASTQ files.
        output_path (str): Path to save the generated meta input TSV file.
    """
    df = pd.read_csv(sra_meta_path)
    records = []
    for _, row in df.iterrows():
        data_id = row[data_id_col] if data_id_col in row else row.get('SRR', '')
        sample_id = row[sample_id_col] if sample_id_col in row else ''
        layout = row[Layout_col] if Layout_col in row else row.get('library_type', 'PAIRED')
        fq1 = os.path.join(fastq_dir, f"{data_id}_1.fastq.gz")
        fq2 = os.path.join(fastq_dir, f"{data_id}_2.fastq.gz")
        fq_single = os.path.join(fastq_dir, f"{data_id}.fastq.gz")
        design = ""
        if str(layout).upper() == 'PAIRED':
            if os.path.exists(fq1) and os.path.exists(fq2):
                records.append([sample_id, data_id, fq1, fq2, design])
            else:
                # fallback to single file if paired not found
                if os.path.exists(fq_single):
                    records.append([sample_id, data_id, fq_single, '', design])
        else:
            if os.path.exists(fq_single):
                records.append([sample_id, data_id, fq_single, '', design])
            else:
                # fallback to paired if single not found
                if os.path.exists(fq1):
                    records.append([sample_id, data_id, fq1, '', design])
    out_df = pd.DataFrame(records, columns=['sample_id', 'data_id', 'fastq_1', 'fastq_2', 'design'])
    out_df.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate meta input file for downstream analysis.")
    parser.add_argument("-i","--sra_meta_path", required=True, help="Path to the SRA metadata CSV file.")
    parser.add_argument("-d","--fastq_dir", required=True, help="Directory containing the FASTQ files.")
    parser.add_argument("-o","--output_path", required=True, help="Path to save the generated meta input TSV file.")
    args = parser.parse_args()
    generate_meta_input(
        sra_meta_path=args.sra_meta_path,
        fastq_dir=args.fastq_dir,
        output_path=args.output_path
    )

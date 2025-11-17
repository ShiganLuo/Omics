"""
not filter common expression gene to analysis SNP
the number of snp in different vcf file is varied, can't use to analyze
"""
import pandas as pd
import argparse
def totalSNP(gzvcf:str,outImage:str):
    df = pd.read_csv(gzvcf,sep="\t",header=None,comment="#")
    
    return df

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate total SNP from gzvcf file")
    parser.add_argument("--control", type=str, help="Input gzvcf file")
    parser.add_argument("--experiment", type=str, help="Input gzvcf file")
    parser.add_argument("--outImage", type=str, default=None, help="Output image file (optional)")
    
    control = ["/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633379.vcf.gz",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633380.vcf.gz"]
    experiment = [
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633381.vcf.gz",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633382.vcf.gz",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633383.vcf.gz",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633384.vcf.gz",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633385.vcf.gz",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633386.vcf.gz"]
    args = parser.parse_args()
    
    for i in control + experiment:
        df_snp = totalSNP(i, "a.png")
        print(f"{i} : {len(df_snp)}")
    # df_snp = totalSNP(args.gzvcf, args.outImage)
    # print(df_snp.head())  # Display the first few rows of the DataFrame

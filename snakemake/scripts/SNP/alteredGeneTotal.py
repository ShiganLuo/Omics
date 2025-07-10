import pandas as pd
def alteredGeneTotal(df:pd.DataFrame):
    altered = df['difference'].sum()
    print(altered)

if __name__ == '__main__':
    infiles = ["/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/SNP/alteredGene.csv",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE183522/SNP/alteredGene.csv",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE185005/SNP/alteredGene.csv",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE204801/SNP/alteredGene.csv",
               "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/SNP/alteredGene.csv"]
    for infile in infiles:
        df = pd.read_csv(infile,sep="\t")
        alteredGeneTotal(df)
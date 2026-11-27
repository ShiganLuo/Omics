#!/bin/bash
function mutationRate {
    vcf=$1
   
    sample=$(basename $vcf .vcf.gz)
    SNPNumber=$(bcftools view -v snps $vcf | grep -v "^#" | wc -l)
    echo "${sample} SNPs: $SNPNumber"
}
vcf=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633379.vcf.gz
commonBed=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/counts/mouseTEcountCommon.bed
for i in /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/*.vcf.gz;do
    mutationRate $i $commonBed
done

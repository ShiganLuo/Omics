#!/bin/bash
function mutationRate {
    vcf=$1
    # TotalBairs=$2
   
    sample=$(basename $vcf .vcf.gz)
    SNPNumber=$(bcftools view -v snps $vcf | grep -v "^#" | wc -l)
    echo "${sample} SNPs: $SNPNumber"
    # commonBed=$2
    # # sort -k1,1V -k2,2n $commonBed > /tmp/sorted_common.bed
    # # TotalBairs=$(bedtools merge -i /tmp/sorted_common.bed | awk '{sum += $3 - $2} END {print sum}')
    # # echo "SNPs: $SNPNumber"
    # # echo "Total Base Pairs: $TotalBairs"
    # mutationRate=$(echo "scale=6; ${SNPNumber} / ${TotalBairs} * 1000000" | bc)
    # echo "${sample} Mutation rate: ${mutationRate} SNPs per million base pairs"
}
vcf=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/SRR13633379.vcf.gz
commonBed=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/counts/mouseTEcountCommon.bed
for i in /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/filter/vcf/mouse/*.vcf.gz;do
    # echo "Processing $i"
    mutationRate $i $commonBed
done
# mutationRate $vcf $commonBed
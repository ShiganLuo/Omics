#!/bin/bash
#conda:RNA-SNP
function vcfIntersectBed(){
    vcf=$1
    bed=$2
    out=$3
    if [[ ! -f ${vcf}.tbi ]];then
        tabix -p vcf ${vcf}
    fi
    bedtools intersect -a ${vcf} -b ${bed} -wa -wb > ${out}
}
export -f vcfIntersectBed
bed=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/RNASNP202503TEcount_common.bed
# vcf=/ChIP_seq_2/StemCells/RNASNP202503/output/filter/vcf/mouse/SRR17762738.vcf.gz
# outfile="output/vcf/$(basename -s .vcf.gz ${vcf}).bed"
# vcfIntersectBed ${vcf} ${bed} a.bed
# echo ${outfile}

infiles=$(find /ChIP_seq_2/StemCells/RNASNP202503/output/filter/vcf/mouse -type f -name "*.vcf.gz")
for i in ${infiles[@]};do
    outfile="output/vcf/$(basename -s .vcf.gz ${i}).bed"
    vcfIntersectBed ${i} ${bed} ${outfile}
done
# vcfIntersectBed ${vcf} ${bed} a.bed
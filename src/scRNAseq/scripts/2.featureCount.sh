#!/bin/bash
featureCounts=/home/lsg/tools/subread-2.0.7/bin/featureCounts
gtf=/home/lsg/Data/glioblastoma/data/star/genes.gtf
star=/home/lsg/Data/glioblastoma/output/star/star.txt
sample=('GBM27' 'GBM28' 'GBM29')
output=/home/lsg/Data/glioblastoma/output/featureCounts
mapfile -t arr < ${star}
for i in {0..2};do

    echo "$featureCounts -t exon -g gene_id -a ${gtf} -o ${output}/${sample[${i}]} ${arr[${i}]}"
done | parallel -j 3

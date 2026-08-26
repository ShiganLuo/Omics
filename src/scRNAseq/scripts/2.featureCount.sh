#!/bin/bash
featureCounts=/home/lsg/tools/subread-2.0.7/bin/featureCounts
gtf=/home/lsg/Data/glioblastoma/data/star/genes.gtf
star=/home/lsg/Data/glioblastoma/output/star/star.txt
sample=('GBM27' 'GBM28' 'GBM29')
output=/home/lsg/Data/glioblastoma/output/featureCounts
mapfile -t arr < ${star}
for i in {0..2};do

    # echo "$scTE -i ${arr[${i}]} -o  ${sample[${i}]} -p 20 -x ${hg38} --min_counts 1 --min_genes 1 -CB CB -UMI UB"
    echo "$featureCounts -t exon -g gene_id -a ${gtf} -o ${output}/${sample[${i}]} ${arr[${i}]}"
done | parallel -j 3

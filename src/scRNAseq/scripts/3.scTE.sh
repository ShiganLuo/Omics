#!/bin/bash
scTE_build=/home/lsg/tools/scTE/bin/scTE_build
scTE=/home/lsg/tools/scTE/bin/scTE

hg38=/home/lsg/Data/glioblastoma/output/index/hg38.exclusive.idx





cellranger=/home/lsg/Data/glioblastoma/output/bam/bam_clean.txt
mapfile -t arr < ${cellranger}
sample=('GBM27' 'GBM28' 'GBM29')
for i in {0..2};do
    echo "$scTE -i ${arr[${i}]} -o ${sample[${i}]} -x ${hg38} --min_counts 1 --min_genes 1 -CB CB -UMI UB"
done | parallel -j 3
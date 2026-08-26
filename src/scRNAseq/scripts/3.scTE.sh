#!/bin/bash
scTE_build=/home/lsg/tools/scTE/bin/scTE_build
scTE=/home/lsg/tools/scTE/bin/scTE
# index
# scTE_build -g hg38 -o /home/lsg/Data/glioblastoma/output/index/hg38 # Human
#########################################star#####################################

hg38=/home/lsg/Data/glioblastoma/output/index/hg38.exclusive.idx
# # analysis

# star=/home/lsg/Data/glioblastoma/output/star/star.txt
# output=/home/lsg/Data/glioblastoma/output/scTE

# mapfile -t arr < ${star}
# sample=('GBM27' 'GBM28' 'GBM29')
# for i in {0..2};do
#     # echo ${arr[${i}]}
#     # echo ${sample[${i}]}

#     # echo "$scTE -i ${arr[${i}]} -o  ${sample[${i}]} -p 20 -x ${hg38} --min_counts 1 --min_genes 1 -CB CB -UMI UB"
#     echo "$scTE -i ${arr[${i}]} -o  ${sample[${i}]} -p 20 -x ${hg38} --min_counts 1 --min_genes 1 -CB CR -UMI UR"
# done | parallel -j 3 


#########################################or cellranger#####################################
cellranger=/home/lsg/Data/glioblastoma/output/bam/bam_clean.txt
mapfile -t arr < ${cellranger}
sample=('GBM27' 'GBM28' 'GBM29')
for i in {0..2};do
    # echo ${arr[${i}]}
    # echo ${sample[${i}]}
    # echo "$scTE -i ${arr[${i}]} -o  ${sample[${i}]} -p 20 -x ${hg38} --min_counts 1 --min_genes 1 -CB CB -UMI UB"
    echo "$scTE -i ${arr[${i}]} -o ${sample[${i}]} -x ${hg38} --min_counts 1 --min_genes 1 -CB CB -UMI UB"
done | parallel -j 3
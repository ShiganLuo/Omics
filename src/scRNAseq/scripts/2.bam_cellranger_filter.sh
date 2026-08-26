#!/bin/bash
in_dir=/home/lsg/Data/glioblastoma/output/bam
# echo ${in_dir}/GBM27/outs/possorted_genome_bam.bam
# echo ${in_dir}/GBM28/outs/possorted_genome_bam.bam
# echo ${in_dir}/GBM29/outs/possorted_genome_bam.bam
sample=('GBM27' 'GBM28' 'GBM29')
for i in {0..2};do
    # echo "samtools view ${in_dir}/${sample[${i}]}/outs/possorted_genome_bam.bam -h | awk '/^@/ || /CB:/' | samtools view -h -b > ${in_dir}/${sample[${i}]}/outs/possorted_genome_bam.clean.bam"
    echo "samtools view ${in_dir}/${sample[${i}]}/outs/possorted_genome_bam.clean.bam -h | awk '/^@/ || /UB:/' | samtools view -h -b > ${in_dir}/${sample[${i}]}/outs/possorted_genome_bam.finalclean.bam"
done | parallel -j 3
# samtools view ${in_dir}/GBM27/outs/possorted_genome_bam.bam -h | awk '/^@/ || (/UB:/ && /CB:/)' | samtools view -h -b > ${in_dir}/GBM27/outs/possorted_genome_bam.clean.bam
# samtools view ${in_dir}/GBM28/outs/possorted_genome_bam.bam -h | awk '/^@/ || (/UB:/ && /CB:/)' | samtools view -h -b > ${in_dir}/GBM28/outs/possorted_genome_bam.clean.bam
# samtools view ${in_dir}/GBM29/outs/possorted_genome_bam.bam -h | awk '/^@/ || (/UB:/ && /CB:/)' | samtools view -h -b > ${in_dir}/GBM29/outs/possorted_genome_bam.clean.bam
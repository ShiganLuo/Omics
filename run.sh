#!/bin/bash
source /home/luosg/Data/genomeStability/workflow/RNA-SNP/scripts/map/star.sh
genome_index=/home/luosg/genome/Human/GENCODE/GRCh38/STAR_index
genome_fa=/home/luosg/genome/Human/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa
genome_gtf=/home/luosg/genome/Human/GENCODE/GRCh38/gencode.v49.primary_assembly.annotation.gtf
STAR=/opt/STAR-2.7.11b/bin/Linux_x86_64/STAR
star_index ${genome_index} ${genome_fa} ${genome_gtf} 149 30 ${STAR}
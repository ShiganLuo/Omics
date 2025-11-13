#!/bin/bash
source /disk5/luosg/Totipotent20251031/workflow/RNA-SNP/scripts/map/star.sh
genome_index=/disk5/luosg/Reference/GENCODE/human/GRCh38/index_star2_7
genome_fa=/disk5/luosg/Reference/GENCODE/human/GRCh38/GRCh38.primary_assembly.genome.fa
genome_gtf=/disk5/luosg/Reference/GENCODE/human/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf
read_length=149
index_core=25
STAR=/opt/STAR-2.7.11b/bin/Linux_x86_64/STAR
star_index ${genome_index} ${genome_fa} ${genome_gtf} ${read_length} ${index_core} ${STAR}
genome_index=/disk5/luosg/Reference/GENCODE/mouse/GRCm39/index_star2_7
genome_fa=/disk5/luosg/Reference/GENCODE/mouse/GRCm39/GRCm39.primary_assembly.genome.fa
genome_gtf=/disk5/luosg/Reference/GENCODE/mouse/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf
star_index ${genome_index} ${genome_fa} ${genome_gtf} ${read_length} ${index_core} ${STAR}

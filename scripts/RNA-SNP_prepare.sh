#!/bin/bash
source /home/luosg/Data/genomeStability/workflow/RNA-SNP/scripts/map/star.sh
genome_index=/home/luosg/genome/Human/GENCODE/GRCh38/STAR_index
genome_fa=/home/luosg/genome/Human/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa
genome_gtf=/home/luosg/genome/Human/GENCODE/GRCh38/gencode.v49.primary_assembly.annotation.gtf
read_length=149
index_core=25
STAR=/opt/STAR-2.7.11b/bin/Linux_x86_64/STAR
# star_index ${genome_index} ${genome_fa} ${genome_gtf} ${read_length} ${index_core} ${STAR}
genome_index=/home/luosg/genome/Mouse/GENCODE/GRCm39/STAR_index
genome_fa=/home/luosg/genome/Mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa
genome_gtf=/home/luosg/genome/Mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.annotation.gtf
# star_index ${genome_index} ${genome_fa} ${genome_gtf} ${read_length} ${index_core} ${STAR}

# source /home/luosg/Data/genomeStability/workflow/RNA-SNP/scripts/gatk/gatkPrepare.sh
# gatk=/opt/gatk-4.6.0.0/gatk
# samtools=/usr/bin/samtools
# genome=/home/luosg/genome/Human/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa
# gatkIndex ${gatk} ${samtools} ${genome}
# genome=/home/luosg/genome/Mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa
# gatkIndex ${gatk} ${samtools} ${genome}

source /home/luosg/Data/genomeStability/workflow/RNA-SNP/scripts/SNP/bash/annovar.sh
retrieve_pl=/home/luosg/tools/annovar/retrieve_seq_from_fasta.pl
genome=/home/luosg/genome/Human/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa
gtf=/home/luosg/genome/Human/GENCODE/GRCh38/gencode.v49.primary_assembly.annotation.gtf
refGene=/home/luosg/genome/Human/GENCODE/GRCh38/annovar/GRCh38/GRCh38_refGene.txt
refGeneMrna=/home/luosg/genome/Human/GENCODE/GRCh38/annovar/GRCh38/GRCh38_refGeneMrna.fa
# GenePred ${gtf} ${refGene}
# retrieve ${retrieve_pl} ${genome} ${refGene} ${refGeneMrna}

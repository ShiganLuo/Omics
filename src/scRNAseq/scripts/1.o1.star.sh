hg38fasta=/home/lsg/Data/glioblastoma/data/ref/refdata-gex-GRCh38-2024-A/fasta/genome.fa
output_index=/home/lsg/Data/glioblastoma/output/index/star
gtf=/home/lsg/Data/glioblastoma/data/star/genes.gtf
#index
# STAR --runMode genomeGenerate \
# --runThreadN 50 \
# --genomeDir ${output_index} \
# --genomeFastaFiles ${hg38fasta} \
# --sjdbGTFfile ${gtf} \
# --sjdbOverhang 99
#sjdbOverhang: max(ReadLength) -1;most comman sense 99 is proper
#comparsion
fq=/home/lsg/Data/glioblastoma/output/index/fq.txt
outdir=/home/lsg/Data/glioblastoma/output/star
cat ${fq} | parallel -j 3 --colsep ' ' \
    STAR --outSAMtype BAM SortedByCoordinate \
    --runThreadN 30 \
    --genomeDir ${output_index} \
    --readFilesIn {1} {2} \
    --outFileNamePrefix ${outdir}/{3} \
    --outSAMattributes NH HI AS nM CR CY UR UY \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 100 \
    --winAnchorMultimapNmax 100 \
    --outMultimapperOrder Random \
    --runRNGseed 777 \
    --outSAMmultNmax 1 \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist None \
    --soloBarcodeReadLength 98



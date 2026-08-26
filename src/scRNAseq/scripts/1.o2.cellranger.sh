cellranger=/home/lsg/tools/cellranger-8.0.1/bin/cellranger
ref=/home/lsg/Data/glioblastoma/data/ref/refdata-gex-GRCh38-2024-A
fqdir=/home/lsg/Data/glioblastoma/data/fq
# wget -c -P /home/lsg/Data/glioblastoma/data/ref/ https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
# refMd5=$(md5sum /home/lsg/Data/glioblastoma/data/ref/refdata-gex-GRCh38-2024-A.tar.gz)
# echo -e "hg38基因组校验值\n${refMd5}"
sample=('GBM29' 'GBM28' 'GBM27')
parallel -j 3 --dry-run \
    $cellranger count \
    --fastqs=${fqdir} \
    --sample={} \
    --transcriptome=${ref} \
    --id={} \
    --output-dir=/home/lsg/Data/glioblastoma/output/bam/{} \
    --create-bam=true \
    --localcores=24 \
    --nosecondary ::: "${sample[@]}"

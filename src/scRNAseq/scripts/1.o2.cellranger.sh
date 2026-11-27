cellranger=/home/lsg/tools/cellranger-8.0.1/bin/cellranger
ref=/home/lsg/Data/glioblastoma/data/ref/refdata-gex-GRCh38-2024-A
fqdir=/home/lsg/Data/glioblastoma/data/fq
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

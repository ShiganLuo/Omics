#!/bin/bash
# python /disk5/luosg/Totipotent20251031/workflow/vcf/sta/TE/TE_SNP_diff.py \
#     --vcf_dir /disk5/luosg/Totipotent20251031/output/SNP/vcf/filter/Homo_sapiens \
#     --group_file /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE224794.tsv \
#     --te_gtf /disk5/luosg/Reference/GENCODE/human/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf \
#     --outdir /disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE224794 \
#     --sample_frac 1

# python /disk5/luosg/Totipotent20251031/workflow/vcf/sta/TE/TE_SNP_diff.py \
#     --vcf_dir /disk5/luosg/Totipotent20251031/output/SNP/vcf/filter/Homo_sapiens \
#     --group_file /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE204801.tsv \
#     --te_gtf /disk5/luosg/Reference/GENCODE/human/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf \
#     --outdir /disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE204801 \
#     --sample_frac 1

# python /disk5/luosg/Totipotent20251031/workflow/vcf/sta/TE/TE_SNP_diff.py \
#     --vcf_dir /disk5/luosg/Totipotent20251031/output/SNP/vcf/filter/Mus_musculus \
#     --group_file /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE166216.tsv \
#     --te_gtf /disk5/luosg/Reference/GENCODE/mouse/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf \
#     --outdir /disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE166216 \
#     --sample_frac 1

# python /disk5/luosg/Totipotent20251031/workflow/vcf/sta/TE/TE_SNP_diff.py \
#     --vcf_dir /disk5/luosg/Totipotent20251031/output/SNP/vcf/filter/Mus_musculus \
#     --group_file /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE185005.tsv \
#     --te_gtf /disk5/luosg/Reference/GENCODE/mouse/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf \
#     --outdir /disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE185005 \
#     --sample_frac 1


python /disk5/luosg/Totipotent20251031/workflow/vcf/sta/TE/TE_SNP_diff.py \
    --vcf_dir /disk5/luosg/Totipotent20251031/PRJNA663159/SNV \
    --group_file /disk5/luosg/Totipotent20251031/data/PRJNA663159/group.tsv \
    --te_gtf /disk5/luosg/Reference/GENCODE/mouse/GRCm38/GRCm38_GENCODE_rmsk_TE.gtf \
    --outdir /disk5/luosg/Totipotent20251031/PRJNA663159/SNV/TE \
    --sample_frac 1
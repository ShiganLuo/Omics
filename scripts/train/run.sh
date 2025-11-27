#!/bin/bash
latent_factor_analysis=/home/luosg/Data/genomeStability/workflow/RNA-SNP/scripts/train/latent_factor_analysis.py
resid=/home/luosg/Data/genomeStability/output/result/matrix/log2tpm229_514_resid.tsv
python ${latent_factor_analysis} \
    --resid ${resid} \
    --method pca \
    --k 5 \
    --prefix /home/luosg/Data/genomeStability/output/result/pca/pca_out


python ${latent_factor_analysis} \
    --resid ${resid} \
    --method fa \
    --k 5 \
    --prefix /home/luosg/Data/genomeStability/output/result/factor/fa 

python ${latent_factor_analysis} \
    --resid ${resid} \
    --method nmf \
    --k 5 \
    --prefix /home/luosg/Data/genomeStability/output/result/nmf/nmf_out

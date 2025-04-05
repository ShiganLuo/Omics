#!/bin/bash
#RNA-SNP
# after RNASNP-noRemov.smk
# you should modify workflow/assets/group.csv to match your actual group
# also --pattern in the following command,the first is control, the second is experiment
set -e
function downstream_mouse(){
    pattern=$1
    outdir=$2
    group=$3
    title=$4 # gene kd/ko
    echo ${pattern}
    echo ${outdir}
    echo ${group}
    echo ${title}
    # ###
    # # Rscript workflow/scripts/DESeq2/DESeq2.r --mode TEcount \
    # #     --matrix output/counts/mouseTEcount.cntTable \
    # #     --group workflow/assets/group.csv \
    # #     --pattern "siCtr" "siMat2a" \
    # #     --outdir output/ \
    # #     --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    # #     --figure pca heatmap volcano \
    # #     --TEcountMode all
    # Rscript workflow/scripts/DESeq2/DESeq2.r --mode TEcount \
    #     --matrix output/${outdir}/counts/mouseTEcount.cntTable \
    #     --group ${group} \
    #     --pattern ${pattern} \
    #     --outdir output/${outdir}/ \
    #     --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    #     --figure pca heatmap volcano \
    #     --TEcountMode all

    # # Rscript workflow/scripts/DESeq2/DESeq2.r --mode TElocal \
    # #     --matrix output/counts/mouseTElocal.cntTable \
    # #     --group workflow/assets/group.csv \
    # #     --pattern "siCtr" "siMat2a" \
    # #     --outdir output/ \
    # #     --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    # #     --figure heatmap volcano
    # Rscript workflow/scripts/DESeq2/DESeq2.r --mode TElocal \
    #     --matrix output/${outdir}/counts/mouseTElocal.cntTable \
    #     --group ${group} \
    #     --pattern ${pattern} \
    #     --outdir output/${outdir}/ \
    #     --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    #     --figure heatmap volcano
    
    # # Rscript workflow/scripts/DESeq2/go-kegg.r \
    # #     --mode go kegg \
    # #     --up output/DESeq2/upDown/TEcount_Gene_up.csv \
    # #     --down output/DESeq2/upDown/TEcount_Gene_down.csv \
    # #     --outdir output/ \
    # #     --title "Mat2a knock-down"
    # Rscript workflow/scripts/DESeq2/go-kegg.r \
    #     --mode go kegg \
    #     --up output/${outdir}/DESeq2/upDown/TEcount_Gene_up.csv \
    #     --down output/${outdir}/DESeq2/upDown/TEcount_Gene_down.csv \
    #     --outdir output/${outdir}/

    # # Rscript workflow/scripts/DESeq2/gsea.r \
    # #     -i output/DESeq2/TEcount_Gene.csv \
    # #     -g data/gsea/mh.all.v2024.1.Mm.symbols.gmt \
    # #     -o output/ \
    # #     -a /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    # #     -t "Hallmark pathways altered by Mat2a in mouse"
    # # mainTitle="Hallmark pathways altered after ${title} in mouse"
    # Rscript workflow/scripts/DESeq2/gsea.r \
    #     -m "Gene" \
    #     -i output/${outdir}/DESeq2/TEcount_Gene.csv \
    #     -g data/gsea/mh.all.v2024.1.Mm.symbols.gmt \
    #     -o output/${outdir}/ \
    #     -a /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv 

    # # gseapy prerank \
    # #     -r output/DESeq2/GSEA/TEcount_Gene_GSEA.rnk \
    # #     -g data/gsea/TwoGenes.gmt \
    # #     -o output/DESeq2/GSEA/TwoGenes \
    # #     -f jpeg
    # gseapy prerank \
    #     -r output/${outdir}/DESeq2/GSEA/TEcount_Gene_GSEA.rnk \
    #     -g data/gsea/TwoGenes.gmt \
    #     -o output/${outdir}/DESeq2/GSEA/TwoGenes \
    #     -f jpeg
    
    Rscript workflow/scripts/DESeq2/TEsite_subfamily.r \
        --up output/${outdir}/DESeq2/upDown/TElocal_TE_up.csv \
        --down output/${outdir}/DESeq2/upDown/TElocal_TE_down.csv \
        --outdir output/${outdir}/ \
        --graphTitle ""

}
export -f downstream_mouse
## GSE166216
pattern="mESC 2C"
outdir=GSE166216
group=workflow/assets/GSE166216.csv
title=GSE166216
# downstream_mouse "${pattern}" ${outdir} ${group} ${title}
## GSE183522
pattern="EPS TPS"
outdir=GSE183522
group=workflow/assets/GSE183522.csv
title=GSE183522
# downstream_mouse "${pattern}" ${outdir} ${group} ${title}
## GSE185005
pattern="mESC ciTotiSC"
outdir=GSE185005
group=workflow/assets/GSE185005.csv
title=GSE185005
# downstream_mouse "${pattern}" ${outdir} ${group} ${title}

function downstream_human(){
    pattern=$1
    outdir=$2
    group=$3
    title=$4 # gene kd/ko
    echo ${pattern}
    echo ${outdir}
    echo ${group}
    echo ${title}
    ###
    Rscript workflow/scripts/DESeq2/DESeq2.r --mode TEcount \
        --matrix output/${outdir}/counts/humanTEcount.cntTable \
        --group ${group} \
        --pattern ${pattern} \
        --outdir output/${outdir}/ \
        --annotation /ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/geneIDAnnotation.csv\
        --figure pca heatmap volcano \
        --TEcountMode all

    Rscript workflow/scripts/DESeq2/DESeq2.r --mode TElocal \
        --matrix output/${outdir}/counts/humanTElocal.cntTable \
        --group ${group} \
        --pattern ${pattern} \
        --outdir output/${outdir}/ \
        --annotation /ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/geneIDAnnotation.csv \
        --figure heatmap volcano
    
    Rscript workflow/scripts/DESeq2/go-kegg.r \
        --mode go kegg \
        --up output/${outdir}/DESeq2/upDown/TEcount_Gene_up.csv \
        --down output/${outdir}/DESeq2/upDown/TEcount_Gene_down.csv \
        --outdir output/${outdir}/ \
        --species human

    Rscript workflow/scripts/DESeq2/gsea.r \
        -m "Gene" \
        -i output/${outdir}/DESeq2/TEcount_Gene.csv \
        -g data/gsea/h.all.v2024.1.Hs.symbols.gmt \
        -o output/${outdir}/ \
        -a /ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/geneIDAnnotation.csv

    gseapy prerank \
        -r output/${outdir}/DESeq2/GSEA/TEcount_Gene_GSEA.rnk \
        -g data/gsea/EightGenes3.gmt \
        -o output/${outdir}/DESeq2/GSEA/EightGenes \
        -f jpeg
    
    Rscript workflow/scripts/DESeq2/TEsite_subfamily.r \
        --up output/${outdir}/DESeq2/upDown/TElocal_TE_up.csv \
        --down output/${outdir}/DESeq2/upDown/TElocal_TE_down.csv \
        --outdir output/${outdir}/ \
        --graphTitle ""
}
export -f downstream_human
## GSE204801
pattern="prEpiSC ci8CLC"
outdir=GSE204801
group=workflow/assets/GSE204801.csv
title=GSE204801
# downstream_human "${pattern}" ${outdir} ${group} ${title}
## GSE224794
pattern="hESC hTBLC"
outdir=GSE224794
group=workflow/assets/GSE224794.csv
title=GSE224794
# downstream_human "${pattern}" ${outdir} ${group} ${title}

function mouse_SNP(){
    ##
    group=$1
    controlStr=$2
    experimentStr=$3
    outdir=$4
    item=$5
    ###  vcf
    mapfile -t control < <(awk -v pat=${controlStr} 'NR > 1 && $2 == pat {print $1}' ${group})
    mapfile -t experiment < <(awk -v pat=${experimentStr} 'NR > 1 && $2 == pat {print $1}' ${group})
    control_vcf=("${control[@]/#/output/${item}/filter/vcf/mouse/}")
    control_vcf=("${control_vcf[@]/%/Common.vcf}")
    experiment_vcf=("${experiment[@]/#/output/${item}/filter/vcf/mouse/}")
    experiment_vcf=("${experiment_vcf[@]/%/Common.vcf}")
    echo ${control_vcf[@]}
    echo ${experiment_vcf[@]}
    python workflow/scripts/SNP/plot.py \
        --mode vcf \
        --control ${control_vcf[@]} \
        --experiment ${experiment_vcf[@]}\
        --outdir ${outdir} \
        --controlName ${controlStr}\
        --experimentName ${experimentStr}
    ### annovar
    declare -a control_annovar
    declare -a experiment_annovar
    for sample in "${control[@]}"; do
        control_annovar+=("output/${item}/annovar/mouse/${sample}/${sample}.GRCm39_multianno.csv")
    done
    for sample in "${experiment[@]}"; do
        experiment_annovar+=("output/${item}/annovar/mouse/${sample}/${sample}.GRCm39_multianno.csv")
    done
    echo ${control_annovar[@]}
    echo ${experiment_annovar[@]}
    python workflow/scripts/SNP/plot.py \
        --mode annovar \
        --control ${control_annovar[@]} \
        --experiment ${experiment_annovar[@]}\
        --outdir ${outdir} \
        --controlName ${controlStr}\
        --experimentName ${experimentStr}
}
export -f mouse_SNP
## GSE166216
group=workflow/assets/GSE166216.csv
controlStr=mESC
experimentStr=2C
outdir=output/GSE166216/SNP/
item=GSE166216
# mouse_SNP ${group} ${controlStr} ${experimentStr} ${outdir} ${item}

## GSE183522
group=workflow/assets/GSE183522.csv
controlStr=EPS
experimentStr=TPS
outdir=output/GSE183522/SNP/
item=GSE183522
# mouse_SNP ${group} ${controlStr} ${experimentStr} ${outdir} ${item}
## GSE185005
group=workflow/assets/GSE185005.csv
controlStr=mESC
experimentStr=ciTotiSC
outdir=output/GSE185005/SNP/
item=GSE185005
# mouse_SNP ${group} ${controlStr} ${experimentStr} ${outdir} ${item}

function human_SNP(){
    ##
    group=$1
    controlStr=$2
    experimentStr=$3
    outdir=$4
    item=$5
    ###  vcf
    mapfile -t control < <(awk -v pat=${controlStr} 'NR > 1 && $2 == pat {print $1}' ${group})
    mapfile -t experiment < <(awk -v pat=${experimentStr} 'NR > 1 && $2 == pat {print $1}' ${group})
    control_vcf=("${control[@]/#/output/${item}/filter/vcf/human/}")
    control_vcf=("${control_vcf[@]/%/Common.vcf}")
    experiment_vcf=("${experiment[@]/#/output/${item}/filter/vcf/human/}")
    experiment_vcf=("${experiment_vcf[@]/%/Common.vcf}")
    echo ${control_vcf[@]}
    echo ${experiment_vcf[@]}
    python workflow/scripts/SNP/plot.py \
        --mode vcf \
        --control ${control_vcf[@]} \
        --experiment ${experiment_vcf[@]} \
        --outdir ${outdir} \
        --controlName ${controlStr} \
        --experimentName ${experimentStr} \
        --species human
    ### annovar
    declare -a control_annovar
    declare -a experiment_annovar
    for sample in "${control[@]}"; do
        control_annovar+=("output/${item}/annovar/human/${sample}/${sample}.GRCh38_multianno.csv")
    done
    for sample in "${experiment[@]}"; do
        experiment_annovar+=("output/${item}/annovar/human/${sample}/${sample}.GRCh38_multianno.csv")
    done
    echo ${control_annovar[@]}
    echo ${experiment_annovar[@]}
    python workflow/scripts/SNP/plot.py \
        --mode annovar \
        --control ${control_annovar[@]} \
        --experiment ${experiment_annovar[@]}\
        --outdir ${outdir} \
        --controlName ${controlStr} \
        --experimentName ${experimentStr} \
        --species human
}
export -f human_SNP
## GSE204801
group=workflow/assets/GSE204801.csv
controlStr=prEpiSC
experimentStr=ci8CLC
outdir=output/GSE204801/SNP/
item=GSE204801
# human_SNP ${group} ${controlStr} ${experimentStr} ${outdir} ${item}
## GSE224794
group=workflow/assets/GSE224794.csv
controlStr=hESC
experimentStr=hTBLC
outdir=output/GSE224794/SNP/
item=GSE224794
human_SNP ${group} ${controlStr} ${experimentStr} ${outdir} ${item}
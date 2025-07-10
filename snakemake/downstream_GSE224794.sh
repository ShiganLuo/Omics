#!/bin/bash
## conda: RNA-SNP
function combine_TEcount() {
    combineTE="/ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/scripts/combineTE.py"
    indir=$1
    outfile=$2
    log=$3
    python ${combineTE} -p TEcount -i ${indir} -o ${outfile} > ${log} 2>&1
    
}
export -f combine_TEcount
indir=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts
outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcount.cntTable
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/combine_TEcount.log
# combine_TEcount ${indir} ${outfile} ${log}

function combine_TElocal() {
    combineTE="/ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/scripts/combineTE.py"
    indir=$1
    outfile=$2
    log=$3
    python ${combineTE} -p TElocal -i ${indir} -o ${outfile} > ${log} 2>&1
}
export -f combine_TElocal
indir=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts
outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTElocal.cntTable
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/combine_TElocal.log
# combine_TElocal ${indir} ${outfile} ${log}

function stringTieMerge() {
    outfile=$1
    refGtf=$2
    log=$3
    gtf=(/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/2pass/*/human/*gtf)

    stringtie --merge ${gtf[@]} -o ${outfile} -G ${refGtf} > ${log} 2>&1
}
export -f stringTieMerge
outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/2pass/human.gtf
refGtf="/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/gencode.v47.primary_assembly.annotation.gtf"
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/stringTieMerge.log
# stringTieMerge ${outfile} ${refGtf} ${log}

# rule TEcountStringTie:
#     input:
#         bam = outdir + "/counts/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam",
#         gtf = outdir + "/2pass/human.gtf"
#     output:
#         project = outdir + "/counts/{sample_id}/human/{sample_id}TEcountStringTie.cntTable"
#     params:
#         project = "{sample_id}TEcountStringTie",
#         outdir = outdir + "/counts/{sample_id}/human",
#         TE_gtf = config['TEtranscripts']['human']['TE_gtf'],
#     threads:5 #防止过多并行运行爆内存
#     log:
#         log=outdir+"/log/human/{sample_id}/TEtranscriptsStringTie.log"
#     conda:
#         config['conda']['TE']
#     shell:
#         """
#         TEcount --sortByPos --format BAM --mode multi \
#         -b {input.bam} --GTF {input.gtf} --TE {params.TE_gtf} \
#         --project {params.project} --outdir {params.outdir} \
#         > {log.log} 2>&1
#         """
function TEcountStringTie() {
    # conda : TE
    bam=$1
    gtf=$2
    project=$3
    outdir=$4
    TE_gtf=$5
    log=$6

    TEcount --sortByPos --format BAM --mode multi \
        -b ${bam} \
        --GTF ${gtf} --TE ${TE_gtf} \
        --project ${project} --outdir ${outdir} > ${log} 2>&1
}
export -f TEcountStringTie

# for bam in /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/*/human/*Aligned.sortedByCoord.out.bam; do
#     sample_id=$(basename -s Aligned.sortedByCoord.out.bam ${bam})
#     # echo ${bam}
#     # echo ${sample_id}
#     gtf=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/2pass/human.gtf
#     project=${sample_id}TEcountStringTie
#     outdir=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/${sample_id}/human
#     TE_gtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf
#     log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/${sample_id}/TEtranscriptsStringTie.log
#     echo "TEcountStringTie ${bam} ${gtf} ${project} ${outdir} ${TE_gtf} ${log}"
# done | parallel -j 5


function combine_TEcountStringTie() {
    combineTE="/ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/scripts/combineTE.py"
    indir=$1
    outfile=$2
    log=$3
    python ${combineTE} -p TEcountStringTie -i ${indir} -o ${outfile} > ${log} 2>&1
}
indir=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts
outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcountStringTie.cntTable
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/combine_TEcountStringTie.log
# combine_TEcountStringTie ${indir} ${outfile} ${log}
function getStringtieBed() {
    script="/ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/scripts/SNP/getBed.py"
    infile=$1
    gtf=$2
    genefile=$3
    TEfile=$4
    TE_gtf=$5
    log=$6
    python ${script} \
        --mode StringTie \
        --input ${infile} \
        --output ${genefile} \
        --output ${TEfile} \
        --Gtf ${gtf} \
        --TEGtf ${TE_gtf} > ${log} 2>&1
}
export -f getStringtieBed
infile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcountStringTie.cntTable
gtf=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/2pass/human.gtf
genefile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/2pass/human_STG.bed
TEfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/2pass/human_TE.bed
TE_gtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/getStringtieBed.log
# getStringtieBed ${infile} ${gtf} ${genefile} ${TEfile} ${TE_gtf} ${log}

outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/2pass/human_StgTEOverlap.bed
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human_StgTEOverlap.log
# bedtools intersect -a ${genefile} -b ${TEfile} -wa -wb > ${outfile} 2>${log}

script=/ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/scripts/SNP/run-NormCountMat.R
infile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcount.cntTable
outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcountCPM.cntTable
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/TEcoutCPM.log
# /usr/bin/Rscript ${script} -i ${infile} -o ${outfile} > ${log} 2>&1

script=/ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/scripts/SNP/commonExpression.py
infile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcountCPM.cntTable
outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcountCommon.cntTable
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/commonExpression.log
# python ${script} --input ${infile} --output ${outfile} --threshold 5 > ${log} 2>&1

script=/ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/scripts/SNP/getBed.py
infile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcountCommon.cntTable
outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcountCommon.bed
gtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/gencode.v47.primary_assembly.annotation.gtf
TE_gtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf
log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/getBed.log
# python ${script} \
#     --input ${infile} \
#     --output ${outfile} \
#     --Gtf ${gtf} \
#     --TEGtf ${TE_gtf} > ${log} 2>&1

# for vcf in /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/filter/vcf/human/*.vcf.gz; do
#     sample_id=$(basename -s .vcf.gz ${vcf})
#     bed=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcountCommon.bed
#     outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/filter/vcf/human/${sample_id}Common.vcf
#     log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/${sample_id}/vcfIntersectBed.log
#     echo "bedtools intersect -a ${vcf} -b ${bed} -wa -wb > ${outfile} 2>${log}"
# done | parallel -j 5


# for vcf in /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/filter/vcf/human/*Common.vcf; do
#     sample_id=$(basename -s Common.vcf ${vcf})
#     avinput=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/annovar/human/${sample_id}/${sample_id}.avinput
#     log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/${sample_id}/annovar_convert.log
#     convert="/opt/annovar/convert2annovar.pl"
#     echo "/usr/bin/perl ${convert} \
#         -format vcf4 \
#         -withfreq ${vcf} > ${avinput} 2>${log}"
# done | parallel -j 5


# for avinput in /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/annovar/human/*/*.avinput; do
#     sample_id=$(basename -s .avinput ${avinput})
#     outfile=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/annovar/human/${sample_id}/${sample_id}.GRCh38_multianno.csv
#     db=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/annovar/GRCh38/
#     buildver=GRCh38
#     table=/opt/annovar/table_annovar.pl
#     out=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/annovar/human/${sample_id}/${sample_id}
#     log=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/log/human/${sample_id}/annovar_table.log
#     echo "/usr/bin/perl ${table} \
#         ${avinput} ${db} \
#         -buildver ${buildver} \
#         -out ${out} \
#         -remove -protocol refGene \
#         -operation g \
#         -nastring . \
#         -csvout > ${log} 2>&1"
# done | parallel -j 5



Rscript /ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/scripts/DESeq2/DESeq2.r --mode TEcount \
    --matrix /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/counts/humanTEcount.cntTable \
    --group /ChIP_seq_2/StemCells/RNASNP_PT/workflow/snakemake/assets/ZLC.csv \
    --pattern hESC ZLC \
    --outdir /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/ZLC/ \
    --annotation /ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/geneIDAnnotation.csv \
    --figure pca heatmap volcano \
    --TEcountMode all


# gseapy prerank \
# -r output/${outdir}/DESeq2/GSEA/TEcount_Gene_GSEA.rnk \
# -g /ChIP_seq_2/Data/index/geneSets/ZLC.gmt \
# -o output/${outdir}/DESeq2/GSEA/CRGs \
# -f jpeg
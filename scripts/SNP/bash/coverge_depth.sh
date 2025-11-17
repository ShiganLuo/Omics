#!/bin/bash
bam=$1
outfile=$2
genomeSize=$3
function depth(){
    bam=$1
    sample_id=$(basename -s .bam ${bam})
    outfile=$2
    genomeSize=$3
        # 只在输出文件不存在时写入标题行
    if [ ! -f ${outfile} ]; then
        echo -e "sample\tgenome_depth\tcoverage_depth\tcoverage\tcoverage_ratio" > ${outfile}
    fi
    # the bam must be sorted
    samtools depth ${bam} | \
        awk -v sample=${sample_id} -v genomeSize=${genomeSize} '{sum+=$3} END { print sample  "\t" sum/genomeSize "\t" sum/NR "\t" NR "\t" NR/genomeSize } ' \
        >> ${outfile}
}
export -f depth
# depth ${bam} ${outfile} ${genomeSize}
# depth /ChIP_seq_2/Data/serum/cfDNA/output/RG/V350318972.bam /ChIP_seq_2/Data/serum/cfDNA/output/RG/stats/depth_coverge.csv

# #### human
# ### GSE204801
# for i in /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE204801/xenofilterR/bam/Filtered_bams/*.bam; do
#     # echo "Processing $i"
#     sample_id=$(basename -s Aligned.sortedByCoord.out_Filtered.bam ${i})
#     # echo ${sample_id}
#     echo "depth ${i} /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE204801/xenofilterR/stats/${sample_id}_depth_coverge.csv 2913022398"
# done | parallel -j 5
# ### GSE224794
# for i in /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/xenofilterR/bam/Filtered_bams/*.bam; do
#     # echo "Processing $i"
#     sample_id=$(basename -s Aligned.sortedByCoord.out_Filtered.bam ${i})
#     # echo ${sample_id}
#     echo "depth ${i} /ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/xenofilterR/stats/${sample_id}_depth_coverge.csv 2913022398"
# done | parallel -j 5
# ### ZLC
# for i in /ChIP_seq_2/StemCells/RNASNP_PT/output/ZLC/xenofilterR/bam/Filtered_bams/*.bam; do
#     # echo "Processing $i"
#     sample_id=$(basename -s Aligned.sortedByCoord.out_Filtered.bam ${i})
#     # echo ${sample_id}
#     echo "depth ${i} /ChIP_seq_2/StemCells/RNASNP_PT/output/ZLC/xenofilterR/stats/${sample_id}_depth_coverge.csv 2913022398"
# done | parallel -j 2
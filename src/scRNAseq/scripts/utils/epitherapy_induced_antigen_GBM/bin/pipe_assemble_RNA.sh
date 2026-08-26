ibam=${1}


samtools view -q 255 $ibam | stringtie - -o ${ibam%.bam}.gtf -e -b ${ibam%.bam}_stats -p 2 -m 100 -c 1 -G reference_merged_candidates.gtf







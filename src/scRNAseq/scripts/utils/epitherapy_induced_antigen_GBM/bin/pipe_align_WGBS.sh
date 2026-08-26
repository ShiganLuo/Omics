fastq1=${1}
fastq2=${2}
odir=${3}

#1. Trim adaptors
cutadapt -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC -o "Trimmed_"$fastq1 $fastq1
cutadapt -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC -o "Trimmed_"$fastq2 $fastq2

#2. Align
bismark -q --score-min L,0,-0.2 --ignore-quals --no-mixed --no-discordant --dovetail --minins 0 --maxins 2000 --bowtie2 -o $odir -1 "Trimmed_"$fastq1 -2 "Trimmed_"$fastq2

#3. Deduplicate
deduplicate_bismark -p --bam "Trimmed_"${fastq1%fastq*} --output_dir $odir


#4. Extract methylation call
bismark_methylation_extractor -paired-end --bedgraph --no_overlap --comprehensive --merge_non_CpG --report -o $odir --gzip "Trimmed_"${fastq1%fastq*}"_bismark_bt2_pe.bam"

fastq1=${1}
fastq2=${2}


#1. trim adaptors
cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --quality-cutoff=15,10 --minimum-length=36 -o "Trimmed_"$fastq1 -p "Trimmed_"$fastq2 $fastq1 $fastq2 > $fastq1"_cutadapt.log"

#2. align
bowtie2 --very-sensitive -p 4 -X 2000 -x ~/reference/bowtie/hg38/hg38_bw -1 "Trimmed_"$fastq1 -2 "Trimmed_"$fastq2 | samtools view -u - | samtools sort - > "Trimmed_"$fastq1".bam"
samtools index "Trimmed_"$fastq1".bam"

#3. Remove chrM
samtools view -b "Trimmed_"$fastq1".bam" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | samtools sort > "Trimmed_"$fastq1"_nochrM.bam"
samtools index "Trimmed_"$fastq1"_nochrM.bam"

#4. Deduplicate
java -jar /opt/apps/picard/2.8.1/picard.jar MarkDuplicates I="Trimmed_"$fastq1"_nochrM.bam" O="Trimmed_"$fastq1"_nochrM_nodup.bam" M="Trimmed_"$fastq1"_nochrM.bam_dups.txt" REMOVE_DUPLICATES=true


#5. Extract properly paired reads
samtools view -h -b -q 10 -f 2 "Trimmed_"$fastq1"_nochrM_nodup.bam" > "Trimmed_"$fastq1"_nochrM_nodup_qffilter.bam"
samtools index "Trimmed_"$fastq1"_nochrM_nodup_qffilter.bam"

#6. Convert bam to tagAlign
samtools sort -n "Trimmed_"$fastq1"_nochrM_nodup_qffilter.bam" | bedtools bamtobed -bedpe -i stdin | awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){print $1,$2-5,$6-5}}' > "Trimmed_"$fastq1"_nochrM_nodup_qffilter.tagAlign"

#7. call peaks
macs2 callpeak -t "Trimmed_"$fastq1"_nochrM_nodup_qffilter.tagAlign" -f BED -n "Trimmed_"$fastq1"_nochrM_nodup_qffilter" -q .05 --nomodel --shift 37 --extsize 73 -B --keep-dup all --call-summits

#8. Remove black list
bedtools intersect -v -a "Trimmed_"$fastq1"_nochrM_nodup_qffilter_peaks.narrowPeak" -b ~/reference/blacklist/hg38.blacklist.bed > "Trimmed_"$fastq1"_nochrM_nodup_qffilter_peaks_noBL.narrowPeak"


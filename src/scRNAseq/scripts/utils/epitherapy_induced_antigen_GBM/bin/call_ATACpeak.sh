ibam=$1
blacklist=$2
bedtools bamtobed -i $ibam | awk -v OFS="\t" '{if($6=="+"){ $2=$2+4; }else if($6=="-"){ $3=$3-5;} print $0; }' > ${ibam%.bam}".PE2SE_tn5.tagAlign"
macs2 callpeak -t ${ibam%.bam}".PE2SE_tn5.tagAlign" -f BED -n ${ibam%.bam}".PE2SE_tn5" -q .05 --nomodel --shift -37 --extsize 73 -B --keep-dup all --call-summits
bedtools intersect -v -a ${ibam%.bam}".PE2SE_tn5_peaks.narrowPeak" -b $blacklist  > ${ibam%.bam}".PE2SE_tn5_peaks_noBL.narrowPeak"



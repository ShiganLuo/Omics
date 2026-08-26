fastq1=${1}
fastq2=${2}


cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50  -o "Trimmed_"$fastq1 -p "Trimmed_"$fastq2 $fastq1 $fastq2 > $fastq1"_cutadapt.log"

STAR  --runMode alignReads  --runThreadN 4  --genomeDir ~/reference/STAR_index_hg38_gencodeV22/ --readFilesIn "Trimmed_"$fastq1 "Trimmed_"$fastq2 --readFilesCommand zcat --outFileNamePrefix "Trimmed_"$fastq1 --outSAMtype BAM   SortedByCoordinate   --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMheaderHD @HD VN:1.4 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33  --alignIntronMax 500000  --alignMatesGapMax 1000000 --twopassMode Basic






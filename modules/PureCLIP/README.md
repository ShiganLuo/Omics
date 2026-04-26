[文档](https://pureclip.readthedocs.io/en/latest/GettingStarted/preprocessing.html)

We remove PCR duplicates based on the mapping position and random barcodes using UMI, which requires the read ID in the format _@HISEQ:87:00000000_BARCODE read1_. Therefore we append the barcode to the read ID prior the mapping, using the following command:

```sh
zcat rep1/reads.R1.trimmed2.fastq.gz > rep1/reads.R1.trimmed2.fastq
zcat rep1/reads.R2.trimmed2.fastq.gz > rep1/reads.R2.trimmed2.fastq
awk -v l=10 'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }; substr($1, 1, 1) != "@" {print}; ' rep1/reads.R1.trimmed2.fastq  | gzip > rep1/reads.R1.trimmed2.bc.fastq.gz
awk -v l=10 'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }; substr($1, 1, 1) != "@" {print}; ' rep1/reads.R2.trimmed2.fastq  | gzip > rep1/reads.R2.trimmed2.bc.fastq.gz
```
where l=10 denotes the used barcode length.
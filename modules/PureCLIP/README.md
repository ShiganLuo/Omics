[文档](https://pureclip.readthedocs.io/en/latest/GettingStarted/preprocessing.html)

## 去除重复
We remove PCR duplicates based on the mapping position and random barcodes using UMI, which requires the read ID in the format _@HISEQ:87:00000000_BARCODE read1_. Therefore we append the barcode to the read ID prior the mapping, using the following command:

```sh
zcat rep1/reads.R1.trimmed2.fastq.gz > rep1/reads.R1.trimmed2.fastq
zcat rep1/reads.R2.trimmed2.fastq.gz > rep1/reads.R2.trimmed2.fastq
awk -v l=10 'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }; substr($1, 1, 1) != "@" {print}; ' rep1/reads.R1.trimmed2.fastq  | gzip > rep1/reads.R1.trimmed2.bc.fastq.gz
awk -v l=10 'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }; substr($1, 1, 1) != "@" {print}; ' rep1/reads.R2.trimmed2.fastq  | gzip > rep1/reads.R2.trimmed2.bc.fastq.gz
```
where l=10 denotes the used barcode length.

## 主要输出结果

chr: Name of the chromosome or scaffold.
start: Position of crosslink site.
end: Position behind crosslink site (start+1).
state: ‘3’ (0-based, corresponds to state ‘4’ in the PureCLIP publication)
score: log posterior probability ratio of the first and second likely state.
strand: + or -

| 文件中的 state | 实际含义             |
| ---------- | ---------------- |
| 0          | 背景               |
| 1          | 覆盖               |
| 2          | 过渡               |
| 3          | crosslink site |

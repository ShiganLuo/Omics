#!/bin/bash
function starPrepare(){
    fa=$1
    readLength=$2
    genome_length=$(awk '!/^>/ {total += length($0)} END {print total}' ${fa})
    number_of_refs=$(grep -c "^>" ${fa})

    value=$(echo "${genome_length} ${number_of_refs} ${read_length}" | awk '{m=($1/$2 > $3 ? $1/$2 : $3); printf "%.0f", log(m)/log(2)}')
    result=$(echo "$value" | awk '{print ($1>18)?18:$1}')
    echo "--genomeChrBinNbits: ${result}"

    value=$(echo "$genome_length" | awk '{printf "%.2f", (log($1)/log(2))/2 - 1}')
    result=$(echo "$value" | awk '{print ($1>14)?14:$1}')
    echo "--genomeSAindexNbases: ${result}"
}
export -f starPrepare
# fa=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/transcripts.fa
# read_length=150
# starPrepare ${fa} ${read_length}
fa=$1
read_length=$2
# starPrepare ${fa} ${read_length}







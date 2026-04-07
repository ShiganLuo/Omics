#!/bin/bash
#!/bin/bash
function gatkIndex() {
    local gatk=$1
    local samtools=$2
    local genome=$3
    
    ${gatk} CreateSequenceDictionary -R ${genome}
    ${samtools} faidx ${genome}
}
export -f gatkIndex
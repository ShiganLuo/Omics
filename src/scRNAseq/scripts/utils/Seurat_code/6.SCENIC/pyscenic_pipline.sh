#!/bin/bash
export PATH=/opt/anaconda3/envs/pyscenic/bin/:$PATH

seurat_path=$1
prefix=$2
group=$3
output_path=$4

input=$seurat_path/${prefix}.rds
output=$output_path/${prefix}_${group}

echo $input
echo $output

Rscript seurat2count.R -i $input -o ${output} -g ${group}
python pyscenic_pipline.py -i $output -s human -o $output

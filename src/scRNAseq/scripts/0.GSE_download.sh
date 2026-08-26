#!/bin/bash
key=/home/lsg/anaconda3/envs/scRNA-seq/etc/asperaweb_id_dsa.openssh
src=anonftp@ftp.ncbi.nlm.nih.gov:/geo/series/GSE139nnn/GSE139448/suppl/GSE139448_RAW.tar
destination=/home/lsg/Data/glioblastoma/data/matrix/GSE139448_RAW.tar
ascp -v -k 1 -T -l 200m -i $key $src $destination
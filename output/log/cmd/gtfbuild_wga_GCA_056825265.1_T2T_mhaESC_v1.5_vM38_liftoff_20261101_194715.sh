#!/bin/bash
set -euo pipefail
minimap2 -cx asm5 -t 48 --secondary=no /data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5_genomic.fa /data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa > /data/16T/luosg/output/gtfbuild/Mmul_T2T_rmsk/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.wga.paf
echo "gtfbuild_wga for GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff at 20261101_194715 completed successfully"

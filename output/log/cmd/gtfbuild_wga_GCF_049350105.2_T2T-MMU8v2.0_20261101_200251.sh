#!/bin/bash
set -euo pipefail
minimap2 -cx asm5 -t 48 --secondary=no /data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/GCF_049350105.2_T2T-MMU8v2.0_genomic.fa /data/16T/Database/Reference/mulatta/ENSEMBL/Mmul_10/Macaca_mulatta.Mmul_10.dna.toplevel.fa > /data/16T/luosg/output/gtfbuild/Mmul_T2T_rmsk/GCF_049350105.2_T2T-MMU8v2.0.wga.paf
echo "gtfbuild_wga for GCF_049350105.2_T2T-MMU8v2.0 at 20261101_200251 completed successfully"

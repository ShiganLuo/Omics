#!/bin/bash
set -euo pipefail
paftools.js liftover -q 5 -l 50000 /data/16T/luosg/output/gtfbuild/Mmul_T2T_rmsk/GCF_049350105.2_T2T-MMU8v2.0.wga.paf /data/16T/Database/Reference/mulatta/ENSEMBL/Mmul_10/rheMac10_rmsk_TE.bed > /data/16T/luosg/output/gtfbuild/Mmul_T2T_rmsk/rheMac10_rmsk_TE.lifted.paftools.bed
python /data/16T/luosg/workflow/Omics/modules/gtfbuild/bin/gtfbuild_restore_bed.py --lifted /data/16T/luosg/output/gtfbuild/Mmul_T2T_rmsk/rheMac10_rmsk_TE.lifted.paftools.bed --source /data/16T/Database/Reference/mulatta/ENSEMBL/Mmul_10/rheMac10_rmsk_TE.bed --out /data/16T/luosg/output/gtfbuild/Mmul_T2T_rmsk/rheMac10_rmsk_TE.lifted.bed --genome GCF_049350105.2_T2T-MMU8v2.0
in_n=$(grep -vc "^#" /data/16T/Database/Reference/mulatta/ENSEMBL/Mmul_10/rheMac10_rmsk_TE.bed || true)
out_n=$(grep -vc "^#" /data/16T/luosg/output/gtfbuild/Mmul_T2T_rmsk/rheMac10_rmsk_TE.lifted.bed || true)
awk -v s="rheMac10_rmsk_TE" -v i="$in_n" -v o="$out_n" 'BEGIN{printf "metric\tvalue\nsource_bed\t%s\nin_intervals\t%s\nout_intervals\t%s\nlift_rate\t%.4f\n", s, i, o, (i>0? o/i : 0)}' > /data/16T/luosg/output/gtfbuild/Mmul_T2T_rmsk/rheMac10_rmsk_TE.lifted.stats.tsv
echo "gtfbuild_liftover_bed for rheMac10_rmsk_TE at 20261101_201925 completed successfully"

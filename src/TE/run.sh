#!/usr/bin/env bash
# TE locus Manhattan plots from TElocal DESeq2 results.
# Inputs per contrast: DESeq2 *_TE.tsv (counts + padj), group.tsv, and the
# species TElocal locus index (.locInd). TPM / delta TPM / mean TPM are
# computed inside te_manhattan.py (--deseq mode).
set -euo pipefail

DNA_PY=/data/pub/zhousha/miniforge3/envs/RNA/bin/python
TE_DIR=/data/pub/zhousha/Totipotent20251031/workflow/Omics/src/TE
DIFF=/data/pub/zhousha/Totipotent20251031/output/pluripotent2totipotent/RNAseq/diff_expression
DB=/data/pub/zhousha/Database/Reference
LOCI_GRCh38=$DB/human/GENCODE/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf.locInd
LOCI_GRCm39=$DB/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf.locInd

run_one () {
    local deseq="$1" group="$2" ctrl="$3" exp="$4" loci="$5" outdir="$6" prefix="$7" ylabel="$8"; shift 8
    local hl_args=()
    local fam
    for fam in "$@"; do hl_args+=(--highlight-family "$fam"); done
    "$DNA_PY" "$TE_DIR/te_manhattan.py" \
        --deseq "$deseq" --group "$group" --control "$ctrl" --experiment "$exp" \
        --loci-tsv "$loci" \
        --signif-col padj --signif-thresh 0.05 \
        --value-thresh 0 --burst-sign up \
        --y-transform log2 \
        --y-break-auto \
        --label-top 10 \
        --ylabel "$ylabel" --hline 0 \
        "${hl_args[@]}" \
        -o "$outdir" --out-prefix "${prefix}.manhattan" -f png
}

run_one "$DIFF/GRCh38/hESC_vs_hTBLC/TElocal/hESC_vs_hTBLC.TElocal_TE.tsv" \
    "$DIFF/GRCh38/hESC_vs_hTBLC/TElocal/group.tsv" hESC hTBLC "$LOCI_GRCh38" \
    "$DIFF/GRCh38/hESC_vs_hTBLC/TElocal/manhattan" hESC_vs_hTBLC \
    "log2(1+|delta TPM|)" \
    L1PA2 L1PA3 L1PA5 HERVH

run_one "$DIFF/GRCh38/prEpiSC_vs_ci8CLC/TElocal/prEpiSC_vs_ci8CLC.TElocal_TE.tsv" \
    "$DIFF/GRCh38/prEpiSC_vs_ci8CLC/TElocal/group.tsv" prEpiSC ci8CLC "$LOCI_GRCh38" \
    "$DIFF/GRCh38/prEpiSC_vs_ci8CLC/TElocal/manhattan" prEpiSC_vs_ci8CLC \
    "log2(1+|delta TPM|)"

run_one "$DIFF/GRCm39/mESC_vs_ciTotiSC/TElocal/mESC_vs_ciTotiSC.TElocal_TE.tsv" \
    "$DIFF/GRCm39/mESC_vs_ciTotiSC/TElocal/group.tsv" mESC ciTotiSC "$LOCI_GRCm39" \
    "$DIFF/GRCm39/mESC_vs_ciTotiSC/TElocal/manhattan" mESC_vs_ciTotiSC \
    "log2(1+|delta TPM|)" \
    MERVL MT2

run_one "$DIFF/GRCm39/mESC_vs_TLSC/TElocal/mESC_vs_TLSC.TElocal_TE.tsv" \
    "$DIFF/GRCm39/mESC_vs_TLSC/TElocal/group.tsv" mESC TLSC "$LOCI_GRCm39" \
    "$DIFF/GRCm39/mESC_vs_TLSC/TElocal/manhattan" mESC_vs_TLSC \
    "log2(1+|delta TPM|)" \
    MERVL MT2

echo ALL_FOUR_DONE

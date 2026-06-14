#!/bin/bash
set -euo pipefail

ref_mouse=/data/pub/zhousha/env/mutation_0.1/f120c8f185498867119c167234a292e3_/lib/python3.10/site-packages/source/resources/non-human/telogator-ref-mouse.fa.gz
max_parallel=2  # 同时运行的样本数

function telogator2_run() {
    local fastq="$1"
    local output="$2"
    local threads="$3"
    mkdir -p "${output}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting telogator2 for ${fastq}"
    telogator2 \
        -i "${fastq}" \
        -o "${output}" \
        -r hifi \
        -c "${threads}" \
        -t "${ref_mouse}" \
        > "${output}/telogator2.log" 2>&1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished telogator2 for ${fastq}"
}
export -f telogator2_run
export ref_mouse

# 样本列表: fastq_path output_dir threads_per_sample
samples=(
    "/data/pub/zhousha/20260207_Exome/output/PacVar/fastq/raw_fastq/DMSO_P20.bam /data/pub/zhousha/20260207_Exome/output/PacVar/repeat/telomere/DMSO_P20 16"
    "/data/pub/zhousha/20260207_Exome/output/PacVar/fastq/raw_fastq/PlaB_P20.bam /data/pub/zhousha/20260207_Exome/output/PacVar/repeat/telomere/PlaB_P20 16"
)

# 并行执行
pids=()
for sample in "${samples[@]}"; do
    read -r fastq outdir threads <<< "${sample}"
    telogator2_run "${fastq}" "${outdir}" "${threads}" &
    pids+=($!)

    # 控制并行数
    while (( $(jobs -rp | wc -l) >= max_parallel )); do
        wait -n 2>/dev/null || true
    done
done

# 等待所有任务完成
failed=0
for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        echo "ERROR: Job ${pid} failed"
        failed=1
    fi
done

if (( failed )); then
    echo "ERROR: Some telogator2 jobs failed"
    exit 1
fi

echo "All telogator2 jobs completed successfully"

#!/bin/bash

function vep_annotation_install() {
    local species=$1 # Example: mus_musculus
    local assembly=$2 # Example: GRCm39
    vep_install -a cf -s ${species} -y ${assembly}
}
export -f vep_annotation_install
# vep_annotation_install mus_musculus GRCm39

# ==========================================
# 1. 合并 SV 文件 (使用 SURVIVOR)
# 参数: <输出文件名> <距离阈值,默认500> <输入VCF列表...>
# ==========================================
function merge_sv() {
    local out_vcf=$1
    local dist=${2:-500}
    shift 2
    local vcf_files=("$@")
    local tmp_dir=$(mktemp -d) # 创建临时目录存放解压后的文件

    echo "[INFO] Preparing files..."
    local decompressed_list=()
    
    # 遍历输入文件，如果是 .gz 则解压
    for vcf in "${vcf_files[@]}"; do
        local base_name=$(basename "$vcf" .gz)
        local tmp_vcf="$tmp_dir/$base_name"
        
        if [[ "$vcf" == *.gz ]]; then
            zcat "$vcf" > "$tmp_vcf"
        else
            cp "$vcf" "$tmp_vcf"
        fi
        decompressed_list+=("$tmp_vcf")
    done

    # 创建临时列表文件
    local tmp_list="$tmp_dir/tmp_vcf_list.txt"
    printf "%s\n" "${decompressed_list[@]}" > "$tmp_list"

    echo "[INFO] Merging ${#vcf_files[@]} VCFs with distance ${dist}bp..."
    SURVIVOR merge "$tmp_list" "$dist" 1 1 1 0 50 "$out_vcf"
    
    # 清理临时文件
    rm -rf "$tmp_dir"
    echo "[SUCCESS] Merged VCF saved to: $out_vcf"
}
export -f merge_sv
out_vcf=/disk5/luosg/Totipotent20251031/PacBio/SV/merged_sv.vcf
vcfs=("/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/DMSO.sv.vcf.gz" "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/PlaB.sv.vcf.gz")
merge_sv "$out_vcf" 500 "${vcfs[@]}"

# ==========================================
# 1. 合并 SV 文件 (使用 SURVIVOR)
# 参数: <输出文件名> <距离阈值,默认500> <输入VCF列表...>
# ==========================================
function merge_sv_jasmine() {
    local out_vcf=$1
    shift
    local vcf_files=("$@")

    # 创建文件列表
    local tmp_list="vcf_list_jasmine.txt"
    printf "%s\n" "${vcf_files[@]}" > "$tmp_list"

    echo "[INFO] Merging with Jasmine..."
    # --output_genotypes: 保留各个样本的基因型（对于比较组间差异必选）
    # --closeness: 相当于 SURVIVOR 的距离阈值（默认约 1000bp，可调）
    jasmine file_list="$tmp_list" out_file="$out_vcf" --output_genotypes
    
    rm "$tmp_list"
    echo "[SUCCESS] Merged VCF saved to: $out_vcf"
}

# ==========================================
# 2. 提取组特有变异 (基于 SUPP_VEC)
# 参数: <合并后的VCF> <期望的SUPP_VEC值,如10> <输出文件名>
# ==========================================
function extract_specific_sv() {
    local in_vcf=$1
    local vec=$2
    local out_vcf=$3

    echo "[INFO] Extracting SVs with SUPP_VEC=${vec}..."
    # 使用 bcftools 提取包含特定 SUPP_VEC 的行
    bcftools view -i "INFO/SUPP_VEC='${vec}'" "$in_vcf" > "$out_vcf"
    echo "[SUCCESS] Extracted VCF saved to: $out_vcf"
}

export -f extract_specific_sv

# ==========================================
# 3. 运行 VEP 注释 (针对 GRCm39)
# 参数: <输入VCF> <输出VCF> <Cache目录路径>
# ==========================================
function annotate_sv_vep() {
    local in_vcf=$1
    local out_vcf=$2
    local cache_dir=$3

    echo "[INFO] Running VEP annotation for GRCm39..."
    vep \
        -i "$in_vcf" \
        -o "$out_vcf" \
        --cache \
        --dir_cache "$cache_dir" \
        --species mus_musculus \
        --assembly GRCm39 \
        --format vcf \
        --vcf \
        --force_overwrite \
        --everything \
        --pick \
        --per_gene \
        --offline \
        --buffer_size 5000
    
    echo "[SUCCESS] Annotation complete: $out_vcf"
}
export -f annotate_sv_vep




import os
import pysam
import pandas as pd
import bisect

def compute_exonic_mutation_burden_from_vcf_dir(
        vcf_dir: str,
        gtf_file: str,
        output_csv: str = "exonic_mutation_burden.csv",
        qual_filter: float = 30,
        main_chroms: list = None,
        chunksize: int = 100000
    ) -> pd.DataFrame:

    if main_chroms is None:
        main_chroms = ["chr" + str(i) for i in range(1,23)] + ["chrX","chrY","chrM"]

    # ---------------------------------------------------------------
    # 1) 构建外显子区间，并预处理二分查找所需结构
    # ---------------------------------------------------------------
    def build_exon_intervals(gtf_file, main_chroms):
        exon_dict = {}
        with open(gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] != "exon":
                    continue
                chrom = fields[0]
                if chrom not in main_chroms:
                    continue
                start = int(fields[3])
                end = int(fields[4])
                exon_dict.setdefault(chrom, []).append((start, end))

        # 合并区间
        def merge(intervals):
            if not intervals:
                return []
            intervals.sort(key=lambda x: x[0])
            merged = [intervals[0]]
            for st, ed in intervals[1:]:
                pst, ped = merged[-1]
                if st <= ped:
                    merged[-1] = (pst, max(ped, ed))
                else:
                    merged.append((st, ed))
            return merged

        for c in exon_dict:
            merged = merge(exon_dict[c])
            # 生成便于二分查找的起点和终点列表
            starts = [s for s, e in merged]
            ends = [e for s, e in merged]
            exon_dict[c] = {"intervals": merged, "starts": starts, "ends": ends}

        return exon_dict

    exon_dict = build_exon_intervals(gtf_file, main_chroms)
    total_exon_length = sum(e - s + 1 for chrom in exon_dict for s, e in exon_dict[chrom]["intervals"])

    # ---------------------------------------------------------------
    # 2) 判断 POS 是否落在外显子区间（用二分法）
    # ---------------------------------------------------------------
    def is_pos_in_exon(chrom, pos):
        if chrom not in exon_dict:
            return False
        starts = exon_dict[chrom]["starts"]
        ends = exon_dict[chrom]["ends"]
        idx = bisect.bisect_right(starts, pos) - 1
        if idx >= 0 and pos <= ends[idx]:
            return True
        return False

    # ---------------------------------------------------------------
    # 3) 计算外显子突变数
    # ---------------------------------------------------------------
    def count_exonic_mutations(vcf_file, qual_filter, main_chroms, chunksize=100000):
        count = 0
        try:
            # 标准VCF
            vcf = pysam.VariantFile(vcf_file)
            vcf_chroms = set(vcf.header.contigs)
            valid_chroms = set(main_chroms) & vcf_chroms
            for chrom in valid_chroms:
                for rec in vcf.fetch(chrom):
                    if rec.qual is not None and rec.qual < qual_filter:
                        continue
                    if is_pos_in_exon(chrom, rec.pos):
                        count += 1
        except:
            # 非标准VCF文本，分块读取
            if vcf_file.endswith(".gz"):
                file_like = vcf_file
                compression = "gzip"
            else:
                file_like = vcf_file
                compression = None

            reader = pd.read_csv(
                file_like,
                sep=None,
                engine="python",
                compression=compression,
                chunksize=chunksize
            )

            required_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]

            for chunk in reader:
                if not all(c in chunk.columns for c in required_cols):
                    raise ValueError(f"{vcf_file} 缺少必要列: {required_cols}")
                chunk = chunk[chunk["CHROM"].isin(main_chroms)]
                # chunk = chunk[chunk["QUAL"].notna() & (chunk["QUAL"] >= qual_filter)]

                for chrom, group in chunk.groupby("CHROM"):
                    for pos in group["POS"]:
                        if is_pos_in_exon(chrom, pos):
                            count += 1

        return count

    # ---------------------------------------------------------------
    # 4) 批量处理目录
    # ---------------------------------------------------------------
    results = []
    vcf_files = [f for f in os.listdir(vcf_dir)
                 if f.endswith(".vcf") or f.endswith(".vcf.gz") or f.endswith("_SNV.txt.gz")]

    for vcf_name in vcf_files:
        vcf_path = os.path.join(vcf_dir, vcf_name)
        mut_count = count_exonic_mutations(
            vcf_path, qual_filter, main_chroms, chunksize
        )
        mut_burden = mut_count / (total_exon_length / 1e6)
        results.append({
            "sample": vcf_name.replace(".vcf.gz","").replace(".vcf","").replace("_SNV.txt.gz",""),
            "mutation_count": mut_count,
            "mutation_burden_per_Mb": mut_burden
        })

    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    return df



if __name__ == "__main__":
    df_result = compute_exonic_mutation_burden_from_vcf_dir(
        vcf_dir="/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN/GSE224794/vcf_nopon",
        gtf_file="/disk5/luosg/Reference/GENCODE/human/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf",
        output_csv="/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN/GSE224794/exonic_mutation_burden.csv"
    )
    df_result = compute_exonic_mutation_burden_from_vcf_dir(
        vcf_dir="/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN/GSE204801/vcf_nopon",
        gtf_file="/disk5/luosg/Reference/GENCODE/human/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf",
        output_csv="/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN/GSE204801/exonic_mutation_burden.csv"
    )
    # df_result = compute_exonic_mutation_burden_from_vcf_dir(
    #     vcf_dir="/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN/GSE185005/vcf_nopon",
    #     gtf_file="/disk5/luosg/Reference/GENCODE/mouse/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf",
    #     output_csv="/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN/GSE185005/exonic_mutation_burden.csv"
    # )
    # df_result = compute_exonic_mutation_burden_from_vcf_dir(
    #     vcf_dir="/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN/GSE166216/vcf_nopon",
    #     gtf_file="/disk5/luosg/Reference/GENCODE/mouse/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf",
    #     output_csv="/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN/GSE166216/exonic_mutation_burden.csv"
    # )
    # df_result = compute_exonic_mutation_burden_from_vcf_dir(
    #     vcf_dir="/disk5/luosg/Totipotent20251031/PRJNA663159/SNV",
    #     gtf_file="/disk5/luosg/Reference/GENCODE/mouse/GRCm38/gencode.vM23.primary_assembly.annotation.gtf",
    #     output_csv="/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/PRJNA663159_exonic_mutation_burden.csv"
    # )

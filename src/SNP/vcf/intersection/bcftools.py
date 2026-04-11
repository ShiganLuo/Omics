import pandas as pd
import os
import subprocess

def ensure_sorted_indexed(vcf_path):
    """
    确保 VCF 文件已排序并建立索引，返回压缩索引后的文件路径（BGZF）
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"{vcf_path} 不存在")
    
    if vcf_path.endswith(".vcf.gz"):
        tbi_path = vcf_path + ".tbi"
        if not os.path.exists(tbi_path):
            print(f"索引不存在，创建索引: {vcf_path}")
            subprocess.run(["tabix", "-p", "vcf", vcf_path], check=True)
        return vcf_path
    else:
        sorted_vcf = vcf_path.replace(".vcf", ".sorted.vcf.gz")
        if not os.path.exists(sorted_vcf):
            print(f"排序并压缩 VCF: {vcf_path}")
            subprocess.run(["bcftools", "sort", "-O", "z", "-o", sorted_vcf, vcf_path], check=True)
            subprocess.run(["tabix", "-p", "vcf", sorted_vcf], check=True)
        return sorted_vcf

def bcftools_supplement(meta: str,
                        indir: str = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/filter",
                        outdir: str = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect",
                        group_col: str = "study_id"):
    df_meta = pd.read_csv(meta, sep="\t")
    df_group = df_meta.groupby(group_col)
    
    for study, df in df_group:
        study_outdir = os.path.join(outdir, study)
        os.makedirs(study_outdir, exist_ok=True)
        
        df_pluri = df[df["condition"] == "pluripotency"][["sample_id", "organism"]]
        df_toti = df[df["condition"] == "totipotency"][["sample_id", "organism"]]
        
        def get_condition_intersect(df_cond, cond_name):
            vcf_list = []
            for _, row in df_cond.iterrows():
                organism = row["organism"].replace(" ", "_")
                sample_id = row["sample_id"]
                vcf_path = os.path.join(indir, organism, f"{sample_id}.vcf.gz")
                if os.path.exists(vcf_path):
                    vcf_path = ensure_sorted_indexed(vcf_path)
                    vcf_list.append(vcf_path)
                else:
                    print(f"Warning: {vcf_path} not found!")
            if not vcf_list:
                return None
            
            intersect_vcf = os.path.join(study_outdir, f"{cond_name}_intersect.vcf.gz")
            
            if len(vcf_list) == 1:
                # 单样本，直接复制并索引
                subprocess.run(["cp", vcf_list[0], intersect_vcf], check=True)
                subprocess.run(["tabix", "-p", "vcf", intersect_vcf], check=True)
            else:
                # 直接用 bcftools isec 输出 BGZF，并保留 INFO（第一文件）
                cmd = ["bcftools", "isec", f"-n=+{len(vcf_list)}", "-w1"] + vcf_list + ["-O", "z", "-o", intersect_vcf]
                subprocess.run(cmd, check=True)
                subprocess.run(["tabix", "-p", "vcf", intersect_vcf], check=True)
            return intersect_vcf
        
        pluri_vcf = get_condition_intersect(df_pluri, "pluripotency")
        toti_vcf = get_condition_intersect(df_toti, "totipotency")
        
        if pluri_vcf and toti_vcf:
            diff_vcf = os.path.join(study_outdir, "toti_only.vcf.gz")
            # 补集：totipotency 中有,pluripotency 中没有
            full_command = f"bcftools isec -C {toti_vcf} {pluri_vcf} -O v -w1 | bgzip -c > {diff_vcf}"
            subprocess.run(full_command, shell=True, check=True)
            subprocess.run(["tabix", "-p", "vcf", diff_vcf], check=True)

if __name__ == "__main__":
    meta = "/disk5/luosg/Totipotent20251031/data/Totipotent.tsv"
    bcftools_supplement(meta)

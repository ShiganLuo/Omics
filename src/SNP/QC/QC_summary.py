import os
from glob import glob
import pandas as pd

def parse_star_log(log_file):
    """
    解析 STAR Log.final.out 提取 QC 指标
    返回字典
    """
    qc = {}
    with open(log_file, "r") as f:
        for line in f:
            if "|" not in line:
                continue
            key, value = [x.strip() for x in line.split("|")]
            # 去掉多余空格，转数字
            value_clean = value.replace("%","").replace(",","")
            try:
                if "." in value_clean:
                    value_num = float(value_clean)
                else:
                    value_num = int(value_clean)
            except:
                value_num = value_clean
            qc[key] = value_num
    return qc


def batch_star_qc(base_dir):
    """
    扫描 base_dir 下所有样本 STAR log 文件，输出 QC 表格
    """
    samples = sorted([d for d in os.listdir(base_dir) 
                      if os.path.isdir(os.path.join(base_dir,d))])
    all_qc = []
    for s in samples:
        log_file = glob(os.path.join(base_dir, s, "Log.final.out"))
        if not log_file:
            continue
        qc = parse_star_log(log_file[0])
        qc["sample"] = s
        all_qc.append(qc)

    df = pd.DataFrame(all_qc)
    # 可按 total reads 排序
    df = df.sort_values("Number of input reads", ascending=False)
    return df


if __name__ == "__main__":
    qc_df = batch_star_qc("RNAseq")
    print(qc_df.head())
    qc_df.to_csv("STAR_QC_summary.csv", index=False)

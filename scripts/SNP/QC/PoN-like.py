import subprocess
from pathlib import Path
from collections import Counter
import logging
from typing import List
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

# -------------------------
# 工具函数
# -------------------------
def run_cmd(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(result.stderr)
        raise RuntimeError(result.stderr)
    return result.stdout


def extract_sites(vcf: Path) -> list[str]:
    """从 VCF 提取 Chr:Pos:Ref:Alt"""
    logging.info(f"提取位点: {vcf}")
    out = run_cmd([
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
        str(vcf)
    ])
    sites = []
    for line in out.strip().splitlines():
        chrom, pos, ref, alt = line.split("\t")
        sites.append(f"{chrom}:{pos}:{ref}:{alt}")
    return sites


# -------------------------
# 1) 主函数：生成 PON
# -------------------------
def generate_pon(
        vcfs: list[Path],
        outdir: Path,
        freq_threshold: int = 2
):
    """
    根据指定的 norm.vcf.gz 文件生成 PON：

    Parameters
    ----------
    vcfs : list[Path]
        指定的 *.norm.vcf.gz 文件列表
    outdir : Path
        输出目录
    freq_threshold : int
        出现次数阈值（默认 >=2）
    """
    outdir.mkdir(parents=True, exist_ok=True)

    if not vcfs:
        raise ValueError("输入文件列表为空！")

    logging.info(f"接收到 {len(vcfs)} 个 VCF 用于生成 PON")

    all_sites = []

    for vcf in vcfs:
        sites = extract_sites(vcf)
        all_sites.extend(sites)

    # 合并计数
    counter = Counter(all_sites)

    # 保存 all_sites.counts.txt
    count_path = outdir / "all_sites.counts.txt"
    with open(count_path, "w") as f:
        for site, cnt in counter.most_common():
            f.write(f"{cnt} {site}\n")
    logging.info(f"写入: {count_path}")

    # 选择 PON
    pon_sites = [s for s, cnt in counter.items() if cnt >= freq_threshold]

    # 写入 pon_sites.txt
    pon_txt = outdir / "pon_sites.txt"
    with open(pon_txt, "w") as f:
        for s in pon_sites:
            f.write(s + "\n")
    logging.info(f"写入: {pon_txt}")

    # 写入 BED-like
    pon_bed = outdir / "pon_sites.bed"
    with open(pon_bed, "w") as f:
        for s in pon_sites:
            chrom, pos, ref, alt = s.split(":")
            pos = int(pos)
            f.write(f"{chrom}\t{pos-1}\t{pos}\t{ref}\t{alt}\n")
    logging.info(f"写入: {pon_bed}")

    logging.info("PON 生成完成！")
    return pon_bed


# -------------------------
# 2) 搜索 norm.vcf.gz（按分组过滤）
# -------------------------
def collect_norm_vcfs(
        indir: Path,
        samples: List[str],
        recursive: bool = False
) -> list[Path]:
    """
    根据分组文件和输入目录，搜寻对应的 *.norm.vcf.gz 文件

    group_file 内容示例：
        A01
        A02
        B01

    匹配文件：文件名包含任意一个 sampleID 且以 .norm.vcf.gz 结尾
    """
    indir = Path(indir)
    if recursive:
        all_files = list(indir.rglob("*.norm.vcf.gz"))
    else:
        all_files = list(indir.glob("*.norm.vcf.gz"))

    logging.info(f"在目录中找到 {len(all_files)} 个 norm.vcf.gz")

    matched = []

    for v in all_files:
        name = v.name
        if any(s in name for s in samples):
            matched.append(v)

    logging.info(f"匹配到 {len(matched)} 个样本对应 VCF 文件")
    return matched

def remove_pon_from_vcfs(vcfs: List[Path], pon_bed: Path, outdir: Path):
    """
    对每个 VCF 去除 PoN 位点，并建立索引

    Parameters
    ----------
    vcfs : List[Path]
        要处理的 VCF 文件列表
    pon_bed : Path
        PoN 位点的 BED 文件
    outdir : Path
        输出目录
    """
    outdir.mkdir(parents=True, exist_ok=True)

    for vcf in vcfs:
        out_vcf = outdir / (vcf.stem + ".noPON.vcf.gz")
        logging.info(f"从 {vcf.name} 去除 PoN 位点 -> {out_vcf.name}")

        # bcftools view -T ^pon_sites.bed sample.vcf.gz -Oz -o sample.noPON.vcf.gz
        run_cmd([
            "bcftools", "view",
            "-T", f"^{pon_bed}",   # ^表示排除 BED 位点
            str(vcf),
            "-Oz",
            "-o", str(out_vcf)
        ])

        # tabix 索引
        run_cmd([
            "tabix", "-p", "vcf", str(out_vcf)
        ])


if __name__ == "__main__":
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/data/Totipotency/Totipotent.tsv", sep="\t")

    for group, newdf in df.groupby("study_id"):
        samples = newdf["sample_id"].to_list()
        orgs = newdf["organism"].unique()
        print(newdf)

        if len(orgs) > 1:
            raise ValueError(f"metadata is wrong, study {group} has multiple organisms: {orgs}")

        if orgs[0] == "Mus musculus":
            vcfs = collect_norm_vcfs("/disk5/luosg/Totipotent20251031/output/SNP/vcf/norm/mouse", samples)
        elif orgs[0] == "Homo sapiens":
            vcfs = collect_norm_vcfs("/disk5/luosg/Totipotent20251031/output/SNP/vcf/norm/human", samples)
        else:
            raise ValueError(f"Unknown organism {orgs[0]} in study {group}")

        outdir = Path("/disk5/luosg/Totipotent20251031/output/SNP/vcf/PoN") / group
        pon_bed = generate_pon(vcfs, outdir,freq_threshold=len(vcfs) - 1)
        outdir = outdir / "vcf_nopon"
        remove_pon_from_vcfs(vcfs,pon_bed,outdir)


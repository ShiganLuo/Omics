import subprocess
from pathlib import Path
import logging
from typing import Union,List
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] [%(threadName)s] %(message)s"
)
from concurrent.futures import ThreadPoolExecutor, as_completed
def run_cmd(cmd):
    """运行系统命令并在失败时抛错"""
    logging.info("CMD: " + " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)

def ensure_bgzip(vcf: Path) -> Path:
    """
    原地 bgzip 压缩：
    - 输入 .vcf    -> 原地生成 .vcf.gz（覆盖）
    - 输入 .vcf.gz -> 不做任何处理
    返回压缩后的 .vcf.gz 路径
    """
    # 已经压缩，直接返回
    if vcf.suffix == ".gz":
        return vcf

    # 未压缩，原地压缩
    logging.info(f"{vcf.name} 未压缩，原地 bgzip...")
    run_cmd(["bgzip", str(vcf)])

    # bgzip 会生成 vcf.gz 并删除原文件
    gz_path = vcf.with_suffix(vcf.suffix + ".gz")
    return gz_path

def ensure_tabix(gz: Path):
    """确保 tabix 索引存在"""
    tbi = gz.with_suffix(gz.suffix + ".tbi")
    if tbi.exists():
        logging.info(f"已存在索引：{tbi.name}")
    else:
        logging.info(f"未发现索引，正在建立 tabix 索引...")
        run_cmd(["tabix", "-p", "vcf", str(gz)])

def process_vcf(vcf_path: Path, ref_fa: Path, outdir: Path):
    """
    对 VCF 执行：
    1）必要时 bgzip 压缩到输出目录
    2）必要时 tabix 建索引
    3）bcftools norm 左对齐 + 多等位拆分
    输出均放入 outdir
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Step 1: 确保压缩
    gz = ensure_bgzip(vcf_path)

    # Step 2: 确保索引
    ensure_tabix(gz)

    # Step 3: 判断是否已经 norm
    if ".norm." in gz.name:
        logging.info("文件名包含 .norm.，认为已经标准化，直接返回。")
        return gz

    norm_out = outdir / (gz.stem + ".norm.vcf.gz")

    logging.info("执行 bcftools norm...")

    run_cmd([
        "bcftools", "norm",
        "-m", "-both",
        "-f", str(ref_fa),
        "-Oz",
        "-o", str(norm_out),
        str(gz)
    ])

    # Step 4: 索引新的 norm 文件
    ensure_tabix(norm_out)

    return norm_out
def search_vcf(
        indir: Union[str, Path],
        suffix: str,
        recursive: bool = False
) -> List[Path]:
    """
    搜索指定目录下所有以 suffix 结尾的文件。

    Parameters
    ----------
    indir : Union[str, Path]
        搜索目录
    suffix : str
        文件后缀，如 '.vcf' 或 '.vcf.gz'
    recursive : bool
        是否递归子目录（默认 False）

    Returns
    -------
    List[Path] : 所有匹配文件的列表
    """
    indir = Path(indir)

    if not indir.is_dir():
        raise NotADirectoryError(f"{indir} is not a directory")

    if recursive:
        files = [p for p in indir.rglob("*") if p.name.endswith(suffix)]
    else:
        files = [p for p in indir.glob("*") if p.name.endswith(suffix)]

    return files



def process_all_vcf(
        indir: Union[str, Path],
        outdir: Union[str, Path],
        ref_fa: Union[str, Path],
        suffix: str = ".vcf.gz",
        threads: int = 4,
        recursive: bool = False
):
    indir = Path(indir)
    outdir = Path(outdir)
    ref_fa = Path(ref_fa)

    # 搜索 VCF 文件
    vcfs = search_vcf(indir, suffix=suffix, recursive=recursive)
    if not vcfs:
        logging.warning(f"目录中没有找到后缀 '{suffix}' 的文件: {indir}")
        return

    logging.info(f"共找到 {len(vcfs)} 个 VCF 文件，开始多线程处理...")

    # 多线程执行
    futures = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for vcf in vcfs:
            futures.append(
                executor.submit(process_vcf, vcf, ref_fa, outdir)
            )

        for fut in as_completed(futures):
            try:
                result = fut.result()
                logging.info(f"完成：{result}")
            except Exception as e:
                logging.error(f"处理某个文件时出错: {e}")
if __name__ == "__main__":
    # process_all_vcf(
    #     indir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/filter/Homo_sapiens",
    #     outdir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/norm/human",
    #     ref_fa = "/disk5/luosg/Reference/GENCODE/human/GRCh38/GRCh38.primary_assembly.genome.fa",
    # )
    process_all_vcf(
        indir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/filter/Mus_musculus",
        outdir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/norm/mouse",
        ref_fa = "/disk5/luosg/Reference/GENCODE/mouse/GRCm39/GRCm39.primary_assembly.genome.fa",
    )
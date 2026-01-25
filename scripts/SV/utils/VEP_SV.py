import subprocess
import os
import shutil
import tempfile
import gzip
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


class VEP_SV:
    def __init__(self, vep_cache_dir, species="mus_musculus", assembly="GRCm39"):
        """
        初始化分析类
        :param vep_cache_dir: VEP 缓存的根目录
        :param species: 物种名称 (如 mus_musculus)
        :param assembly: 基因组版本 (如 GRCm39)
        """
        self.vep_cache_dir = str(Path(vep_cache_dir).expanduser().resolve()) # 支持 ~ 和绝对路径
        self.species = species
        self.assembly = assembly
        
        # 自动创建缓存根目录
        if not os.path.exists(self.vep_cache_dir):
            os.makedirs(self.vep_cache_dir, exist_ok=True)
            logger.info(f"Created VEP cache directory: {self.vep_cache_dir}")


    def _run_cmd(self, cmd:list):
        """
        执行外部命令，返回 stdout
        - 命令不存在：给出清晰提示
        - 命令执行失败：打印 stdout / stderr
        """
        cmd_str = " ".join(cmd)
        cmd_bin = cmd[0]

        logger.info(f"Running: {cmd_str}")

        # 1️⃣ 预检查：命令是否存在（比 FileNotFoundError 更友好）
        if shutil.which(cmd_bin) is None:
            logger.error(f"Command not found: '{cmd_bin}'")
            logger.error("Please make sure it is installed and in $PATH")
            raise RuntimeError(f"Command not found: {cmd_bin}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )

            if result.stdout:
                logger.info(f"Command Output:\n{result.stdout}")

            return result.stdout

        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with return code {e.returncode}")
            logger.error(f"STDOUT:\n{e.stdout or '[empty]'}")
            logger.error(f"STDERR:\n{e.stderr or '[empty]'}")
            raise RuntimeError(
                f"Command execution failed: {cmd_str}"
            ) from e


    def vep_annotation_install(self):
        """安装类指定的 VEP 缓存"""
        cmd = [
            "vep_install", "-a", "cf", 
            "-s", self.species, 
            "-y", self.assembly, 
            "-c", self.vep_cache_dir
        ]
        self._run_cmd(cmd)

    def merge_sv_survivor(self, vcf_files:list, out_vcf:str, dist:int = 500, min_support:int = 1):
            """
            强制使用本地绝对路径，解决 SURVIVOR 读取失败问题(SURVIVOR在读取临时文件时有时有bug,多个样本时尤其明显,会被认为是同一个样本)
            """
            # 在输出目录旁创建一个实体的临时文件夹
            work_dir = os.path.dirname(os.path.abspath(out_vcf))
            tmp_folder = os.path.join(work_dir, "survivor_tmp")
            os.makedirs(tmp_folder, exist_ok=True)
            
            decompressed_list = []
            try:
                for i, vcf in enumerate(vcf_files):
                    # 显式获取样本名并生成绝对路径
                    s_name = f"Sample_{i}" # 强制用不同名字避免任何冲突
                    tmp_vcf = os.path.abspath(os.path.join(tmp_folder, f"S{i}.vcf"))
                    
                    logger.info(f"Decompressing {vcf} to {tmp_vcf}")
                    if vcf.endswith(".gz"):
                        with gzip.open(vcf, 'rb') as f_in, open(tmp_vcf, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    else:
                        shutil.copy2(vcf, tmp_vcf)
                    decompressed_list.append(tmp_vcf)

                # 写入列表文件，强制使用 \n 换行
                tmp_list = os.path.abspath(os.path.join(tmp_folder, "vcf_list.txt"))
                with open(tmp_list, "w", newline='\n') as f:
                    for path in decompressed_list:
                        f.write(f"{path}\n")

                logger.info(f"Final VCF list path: {tmp_list}")
                logger.info(f"VCF list content:\n{open(tmp_list).read()}")

                # 执行 SURVIVOR
                cmd = ["SURVIVOR", "merge", tmp_list, str(dist), str(min_support), "1", "1", "0", "50", out_vcf]
                self._run_cmd(cmd)

                # 验证结果
                with open(out_vcf, 'r') as f:
                    for line in f:
                        if line.startswith("#CHROM"):
                            cols = line.strip().split('\t')
                            logger.info(f"Merge Check - Columns: {len(cols)}, Samples: {cols[9:]}")
                            break
            finally:
                # 如果成功了就清理，没成功你可以注释掉这行进去看文件在不在
                if os.path.exists(out_vcf) and os.path.getsize(out_vcf) > 0:
                    shutil.rmtree(tmp_folder)
                    logger.info("Cleaned up temporary files.")

    def extract_specific_sv(self, in_vcf:str, out_vcf:str, vec:str = "10"):
        """
        提取特定的 SUPP_VEC 变异
        SUPP_VEC 是一个二进制字符串，表示每个样本中变异的支持情况，由 SURVIVOR 生成
        """
        cmd = ["bcftools", "view", "-i", f"INFO/SUPP_VEC='{vec}'", in_vcf, "-o", out_vcf]
        self._run_cmd(cmd)
        logger.info(f"Extracted SVs (VEC={vec}) to: {out_vcf}")

    def annotate_sv_vep(self, in_vcf:str, out_vcf:str):
        """运行 VEP 注释，使用类初始化的参数"""
        # 检查特定物种的缓存子目录是否存在
        species_cache = os.path.join(self.vep_cache_dir, self.species)
        logger.info(species_cache)
        if not os.path.exists(species_cache):
            logger.warning(f"Cache for {self.species} not found. Attempting install...")
            self.vep_annotation_install()

        cmd = [
            "vep", "-i", in_vcf, "-o", out_vcf,
            "--cache", "--dir_cache", self.vep_cache_dir,
            "--species", self.species,
            "--assembly", self.assembly,
            "--format", "vcf", "--vcf", "--force_overwrite",
            "--everything", "--pick", "--per_gene", "--offline"
        ]
        self._run_cmd(cmd)
        logger.info(f"VEP annotation finished: {out_vcf}")

if __name__ == "__main__":
    analysis = VEP_SV(
        vep_cache_dir="/home/luosg/.vep",
        species="mus_musculus",
        assembly="GRCm39"
    )
    
    raw_vcfs = [
        "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/DMSO.sv.vcf.gz",
        "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/PlaB.sv.vcf.gz"
    ]
    work_dir = "/disk5/luosg/Totipotent20251031/PacBio/SV"
    os.makedirs(work_dir, exist_ok=True)

    merged = os.path.join(work_dir, "merged_sv.vcf")
    specific = os.path.join(work_dir, "PlaB_only.vcf")
    annotated = os.path.join(work_dir, "PlaB_only_annotated.vcf")

    # 流程运行
    analysis.merge_sv_survivor(raw_vcfs, merged)
    analysis.extract_specific_sv(merged, specific, vec="01")
    analysis.annotate_sv_vep(specific, annotated)
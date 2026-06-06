import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from common.LogUtil import setup_logger
import subprocess
import shutil
import gzip
import logging
from pathlib import Path
import pandas as pd

logger = setup_logger("VEP_SV", level=logging.INFO)


class VEP_SV:
    def __init__(self, vep_cache_dir, species="mus_musculus", assembly="GRCm39"):
        """
        Initialize the VEP_SV analysis class.

        Parameters
        ----------
        vep_cache_dir : str
            Root directory for the VEP cache. Created automatically if it
            does not exist.
        species : str, optional
            Species name used by VEP. Default is ``"mus_musculus"``.
        assembly : str, optional
            Genome assembly version. Default is ``"GRCm39"``.

        Returns
        -------
        None
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
        Execute an external command and return its stdout.

        Checks that the command binary exists before execution and raises
        informative errors on failure.

        Parameters
        ----------
        cmd : list of str
            Command and arguments as a list (e.g. ``["vep", "-i", "in.vcf"]``).

        Returns
        -------
        str
            Standard output of the executed command.

        Raises
        ------
        RuntimeError
            If the command binary is not found or execution fails.
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
        """
        Install VEP cache for the configured species and assembly.

        Calls ``vep_install`` with the instance's species, assembly, and
        cache directory settings.

        Returns
        -------
        None
        """
        cmd = [
            "vep_install", "-a", "cf", 
            "-s", self.species, 
            "-y", self.assembly, 
            "-c", self.vep_cache_dir
        ]
        self._run_cmd(cmd)

    def merge_sv_survivor(self, vcf_files:list, out_vcf:str, dist:int = 500, min_support:int = 1):
        """
        Merge SVs from multiple VCF files using SURVIVOR.

        Handles decompression of gzipped VCFs, writes a temporary file list,
        and cleans up intermediate files on success. Uses absolute paths for
        all temporary files to avoid SURVIVOR reading issues.

        Parameters
        ----------
        vcf_files : list of str
            Input VCF file paths (plain or gzip-compressed).
        out_vcf : str
            Output merged VCF file path.
        dist : int, optional
            Maximum distance between SVs to merge (bp). Default is ``500``.
        min_support : int, optional
            Minimum number of supporting samples. Default is ``1``.

        Returns
        -------
        None
            Writes the merged VCF to *out_vcf*.
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
        Extract SVs matching a specific SUPP_VEC value using bcftools.

        SUPP_VEC is a binary string indicating per-sample support, generated
        by SURVIVOR during merging.

        Parameters
        ----------
        in_vcf : str
            Input VCF file path.
        out_vcf : str
            Output VCF file path for extracted variants.
        vec : str, optional
            SUPP_VEC value to filter by. Default is ``"10"``.

        Returns
        -------
        None
            Writes filtered VCF to *out_vcf*.
        """
        cmd = ["bcftools", "view", "-i", f"INFO/SUPP_VEC='{vec}'", in_vcf, "-o", out_vcf]
        self._run_cmd(cmd)
        logger.info(f"Extracted SVs (VEC={vec}) to: {out_vcf}")

    def annotate_sv_vep(
            self, 
            in_vcf:str, 
            outfile:str,
            result_format:str = "vcf"
        ):
        """
        Annotate structural variants using Ensembl VEP.

        Automatically installs the VEP cache for the configured species if
        not already present. Runs VEP in offline mode with ``--everything``
        and ``--pick`` flags.

        Parameters
        ----------
        in_vcf : str
            Input VCF file containing SVs.
        outfile : str
            Output file path for annotated results.
        result_format : str, optional
            Output format passed to VEP via ``--{result_format}``.
            Default is ``"vcf"``.

        Returns
        -------
        None
            Writes annotated output to *outfile*.
        """
        # 检查特定物种的缓存子目录是否存在
        species_cache = os.path.join(self.vep_cache_dir, self.species)
        logger.info(species_cache)
        if not os.path.exists(species_cache):
            logger.warning(f"Cache for {self.species} not found. Attempting install...")
            self.vep_annotation_install()

        cmd = [
            "vep", "-i", in_vcf, "-o", outfile,
            "--cache", "--dir_cache", self.vep_cache_dir,
            "--species", self.species,
            "--assembly", self.assembly,
            "--format", "vcf", f"--{result_format}", "--force_overwrite",
            "--everything", "--pick", "--per_gene", "--offline"
        ]
        self._run_cmd(cmd)
        logger.info(f"VEP annotation finished: {outfile}")


def read_vep_tab(
        table_file:str,
        col_line_prefiex: str = "#Uploaded_variation"
):
    """
    Read a VEP tab-delimited annotation result file into a DataFrame.

    Locates the header line by searching for a line starting with
    *col_line_prefiex*, then reads the remaining rows as a DataFrame.

    Parameters
    ----------
    table_file : str
        Path to the VEP tab-delimited output file.
    col_line_prefiex : str, optional
        Prefix used to identify the header line. Default is
        ``"#Uploaded_variation"``.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the VEP annotation results.

    Raises
    ------
    ValueError
        If no header line matching *col_line_prefiex* is found.
    """
    header_line = None
    with open(table_file) as f:
        for i, line in enumerate(f):
            if line.startswith(col_line_prefiex):
                header_line = i
                break
    if header_line is None:
        raise ValueError(f"Header line starting with '{col_line_prefiex}' not found in {table_file}")
    df = pd.read_csv(table_file, sep='\t', skiprows=header_line)
    return df

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
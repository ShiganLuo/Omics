import os
import shutil
import re
import logging
import pandas as pd
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
from enum import Enum, unique
import argparse
import math
logger = logging.getLogger(__name__)
@unique
class FastqMode(str, Enum):
    FASTQ_META = "FASTQ_META"
    FASTQ_DIR = "FASTQ_DIR"
@unique
class Layout(str, Enum):
    SE = "SE"
    PE = "PE"
    UNKNOWN = "UNKNOWN"

@unique
class MERIPDesign(str, Enum):
    IP = "ip"
    INPUT = "input"
    TREATED_IP = "treated_ip"
    TREATED_INPUT = "treated_input"

@dataclass
class SampleInfo:
    sample_id: str = ""
    organism: str = ""
    layout: Layout = Layout.UNKNOWN
    fastq_1: Optional[Path] = None
    fastq_2: Optional[Path] = None

DESIGN_PATTERN = re.compile(r"^(ctr|exp)_(.+)$")
class MetadataUtils:
    """
    Utilities for variant-analysis metadata parsing and FASTQ preparation.

    note:
    Each data_id corresponds to a single FASTQ file, 
    while the relationship between sample_id and data_id can be either one-to-one or one-to-many.

    Features:
    - Supports meta with explicit fastq paths or only sample_id + design.
    - Validates fastq existence.
    - Determines SE/PE.
    - Handles sample_id + data_id read merging.
    - Establishes standardized symlinks in work directory.
    """

    def __init__(
        self,
        outdir: str,
        meta: Optional[str] = None,
        fastq_dir: Optional[str] = None,
        required_cols: set = {"sample_id", "fastq_1", "fastq_2"},
        data_id_col: str = "data_id",
        design_col: str = "design"
    ):
        """
        Function: Initialize MetadataUtils.
        Parameters:
            - outdir: Output directory for processed FASTQ and logs.
            - meta: Path to metadata file (CSV/TSV) containing sample information and optionally FASTQ paths.
            - fastq_dir: Directory containing FASTQ files (if not specified in meta).
            - required_cols: Set of required columns in the metadata file. Default includes 'sample_id', 'fastq_1', 'fastq_2'.
            - data_id_col: Column name in metadata that represents unique FASTQ identifiers (default: 'data_id').
            - design_col: sample compare mode
        Note:
            - fq_pattern: Glob pattern to identify FASTQ files in fastq_dir (default: '*fq.gz').

        """
        if not meta and not fastq_dir:
            raise ValueError("Either meta or fastq_dir must be provided.")
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.meta = Path(meta) if meta else None
        self.fastq_dir = Path(fastq_dir) if fastq_dir else None
        self.required_cols = required_cols        
        self.data_id_col = data_id_col
        self.design_col = design_col
        self.samples_dict = defaultdict(SampleInfo)

    def load_meta(self) -> pd.DataFrame:
        """
        function: load metadata from meta file, sep can be \t or ,
        """
        with open(self.meta, "r", encoding="utf-8") as f:
            head = f.read(2048)

        sep = "\t" if head.count("\t") >= head.count(",") else ","
        df = pd.read_csv(self.meta, sep=sep)

        return df



    def build_design_pairs(
            self
        ) -> List[Tuple[str, str, str]]:
        """
        Determine ctr/exp pairs based on the design stored in self.samples_dict.
        ctr_x vs exp_x1, exp_x2 ...
        only get the first ctr for each exp, if there are multiple ctr with the same tag, will log a warning and only use the first one.
        return a list of tuple: (organism, ctr_sample_id, exp_sample_id)
        """
        groups: Dict[str, Dict[str, List[SampleInfo]]] = defaultdict(lambda: defaultdict(list))
        design_col = self.design_col
        for sample_id, info in self.samples_dict.items():
            design_val = getattr(info, design_col, "")
            if design_val is None:
                logger.info(f"{sample_id} design value is None, skipping it")
                continue
            if isinstance(design_val, bytes):
                design_val = design_val.decode("utf-8")

            if isinstance(design_val, float) and math.isnan(design_val):
                logger.info(f"{sample_id} design value is None, skipping it")
                continue
            m  = DESIGN_PATTERN.match(design_val)
            if not m:
                logger.warning(f"Invalid design format for {sample_id}: {design_val}")
                continue
            role, tag = m.groups()
            groups[tag][role].append(info)
        
        pairs = []
        for tag, g in groups.items():
            if "ctr" not in g or "exp" not in g:
                logger.warning(f"Incomplete design group for tag '{tag}': missing ctr or exp")
                continue
            for exp_sample_info in g["exp"]:
                if len(g["ctr"]) > 1:
                    logger.warning(f"Multiple ctr samples for tag '{tag}': {g['ctr']}. Only using the first one: {g['ctr'][0].sample_id}")
                pairs.append((exp_sample_info.organism,g["ctr"][0].sample_id, exp_sample_info.sample_id))  # 每个 exp 对应 ctr
        return pairs


    def prepare_fastq_meta(
            self, 
            df: pd.DataFrame,
            outdir:Path,
            sample_id_col:str = 'sample_id', 
            data_id_col:str = 'data_id',
            design_col:str = 'design',
            fastq_r1_col:str = 'fastq_1',
            fastq_r2_col:str = "fastq_2",
            organism_col:str = "organism"
            ) -> None:
        """
            data_id represents a unique FASTQ file.
            If the relationship between sample_id and fastq is one-to-one, a symbolic link is created with the filename prefixed by sample_id.
            If the relationship is one-to-many, FASTQ files corresponding to different data_ids are merged and renamed using the sample_id prefix.

            supplement smaple_id,layout,fastq_1 or fastq_2 information
        """

        if data_id_col not in df.columns:
            df[data_id_col] = df[sample_id_col]

        if not self.required_cols.issubset(df.columns):
            raise ValueError(f"Metadata must contain columns: {self.required_cols}")


        raw_fq_dir = outdir / "raw_fastq"
        raw_fq_dir.mkdir(parents=True, exist_ok=True)

        newdf = df.groupby(sample_id_col)

        for sample_id, df_sample in newdf:
            data_ids = df_sample[data_id_col].values
            if len(data_ids) < 1:
                raise ValueError(f"something wrong: {sample_id} have no {data_id_col} meta")
            
            self.samples_dict[sample_id].sample_id = sample_id
            self.samples_dict[sample_id].design = df_sample[design_col].values[0] if design_col in df_sample.columns else ""
            self.samples_dict[sample_id].organism = df_sample[organism_col].values[0] if organism_col in df_sample.columns else "UNKNOWN"

            if len(data_ids) == 1:
                logger.info(f"Detect the relationship between {sample_id} and {data_ids[0]} is one-to-one")
                origin_r1 = df_sample[fastq_r1_col].values[0]
                origin_r2 = df_sample[fastq_r2_col].values[0] if fastq_r2_col in df_sample.columns else None
                origin_r1 = Path(origin_r1) if origin_r1 else None
                origin_r2 = Path(origin_r2) if origin_r2 else None
                
                if pd.notna(origin_r1) and pd.notna(origin_r2):
                    logger.info(f"Detect {data_ids[0]} is Paired END")
                    self.samples_dict[sample_id].layout = Layout.PE
                    rename_r1 = raw_fq_dir / f"{sample_id}_1.fq.gz"
                    rename_r2 = raw_fq_dir / f"{sample_id}_2.fq.gz"
                    self._link_file(origin_r1,rename_r1)
                    self._link_file(origin_r2,rename_r2)
                    self.samples_dict[sample_id].fastq_1 = rename_r1
                    self.samples_dict[sample_id].fastq_2 = rename_r2
                elif pd.notna(origin_r1):
                    logger.info(f"Detect {data_ids[0]} is Single End")
                    self.samples_dict[sample_id].layout = Layout.SE
                    rename_r1 = raw_fq_dir / f"{sample_id}.single.fq.gz"
                    self._link_file(origin_r1,rename_r1)
                    self.samples_dict[sample_id].fastq_1 = rename_r1
                else:
                    logger.warning(f"{sample_id} have no fastqs, skip it")
                    continue
            else:
                logger.info(f"Detect the relationship between {sample_id} and {data_ids[0]} is one-to-many")
                origin_r1_list = sorted([r for r in df_sample[fastq_r1_col].values if pd.notna(r)])
                origin_r2_list = sorted([r for r in df_sample[fastq_r2_col].values if pd.notna(r)]) if fastq_r2_col in df_sample.columns else []
                
                origin_r1_list_path = [Path(r1) for r1 in origin_r1_list]
                origin_r2_list_path = [Path(r2) for r2 in origin_r2_list]

                if len(origin_r1_list_path) > 0 and len(origin_r2_list_path) > 0:
                    logger.info(f"Detect the fastq of {sample_id} is Paired END")
                    self.samples_dict[sample_id].layout = Layout.PE
                    merge_rename_r1 = raw_fq_dir / f"{sample_id}_1.fq.gz"
                    merge_rename_r2 = raw_fq_dir / f"{sample_id}_2.fq.gz"
                    self._merge_files(origin_r1_list_path, merge_rename_r1)
                    self._merge_files(origin_r2_list_path, merge_rename_r2)
                    self.samples_dict[sample_id].fastq_1 = merge_rename_r1
                    self.samples_dict[sample_id].fastq_2 = merge_rename_r2
                elif len(origin_r1_list_path) > 0:
                    logger.info(f"Detect the fastq of {sample_id} is Single END")
                    self.samples_dict[sample_id].layout = Layout.SE
                    merge_rename_r1 = raw_fq_dir / f"{sample_id}.single.fq.gz"
                    self._merge_files(origin_r1_list_path, merge_rename_r1)
                    self.samples_dict[sample_id].fastq_1 = merge_rename_r1
                else:
                    logger.warning(f"{sample_id} have no fastqs, skip it")
                    continue


    def prepare_fastq_dir(
        self,
        fq_dir: Path,
        outdir: Path,
        fq_pattern: str = r"\.f(ast)?q.gz$"
    ) -> None:
        """
        自动检测 FASTQ，处理多 Lane 合并或单文件软连，并填充 self.samples_dict。
        修复：
        1. 单端测序文件 sample_id 不应带 _1/_2 后缀。
        2. 单端文件命名为 {sample_id}.fq.gz，双端为 {sample_id}_1.fq.gz/{sample_id}_2.fq.gz。
        3. 能正确识别单端和双端。
        """
        temp_files = defaultdict(lambda: {"fastq_1": [], "fastq_2": []})

        logger.info(f"Scanning directory: {fq_dir} with pattern: {fq_pattern}")

        for fq_file in fq_dir.rglob("*"):
            fq_name = fq_file.name
            if not re.search(fq_pattern, fq_name):
                logger.debug(f"Skipping non-FASTQ file: {fq_name}")
                continue

            # 优先识别 _R1/_R2 或 _1/_2，sample_id 不带 lane/read后缀
            m = re.match(r"(.+?)(?:_R?([12]))[^/]*\.f(ast)?q(?:\.gz)?$", fq_name)
            if m:
                sample_id, read_num = m.group(1), m.group(2)
                if read_num == "1":
                    temp_files[sample_id]["fastq_1"].append(fq_file)
                elif read_num == "2":
                    temp_files[sample_id]["fastq_2"].append(fq_file)
            else:
                # 单端：去掉扩展名
                sample_id = re.sub(r"\.(f(ast)?q)(\.gz)?$", "", fq_name)
                temp_files[sample_id]["fastq_1"].append(fq_file)
                logger.warning(f"File {fq_name} did not match R1 or R2 patterns, treat as SE: sample_id={sample_id}")

        raw_fq_dir = outdir / "raw_fastq"
        raw_fq_dir.mkdir(parents=True, exist_ok=True)

        for sample_id, reads in temp_files.items():
            sample_info = self.samples_dict[sample_id]
            sample_info.sample_id = sample_id

            files_r1 = sorted(reads["fastq_1"])
            files_r2 = sorted(reads["fastq_2"])

            if files_r1 and files_r2:
                # PE
                target_r1 = raw_fq_dir / f"{sample_id}_1.fq.gz"
                target_r2 = raw_fq_dir / f"{sample_id}_2.fq.gz"
                if len(files_r1) > 1:
                    logger.info(f"[{sample_id}] Merging {len(files_r1)} R1 files into {target_r1.name}")
                    self._merge_files(files_r1, target_r1)
                else:
                    logger.info(f"[{sample_id}] Creating symlink for {target_r1.name}")
                    self._link_file(files_r1[0], target_r1)
                if len(files_r2) > 1:
                    logger.info(f"[{sample_id}] Merging {len(files_r2)} R2 files into {target_r2.name}")
                    self._merge_files(files_r2, target_r2)
                else:
                    logger.info(f"[{sample_id}] Creating symlink for {target_r2.name}")
                    self._link_file(files_r2[0], target_r2)
                sample_info.fastq_1 = target_r1
                sample_info.fastq_2 = target_r2
                sample_info.layout = Layout.PE
            elif files_r1:
                # SE
                target_se = raw_fq_dir / f"{sample_id}.single.fq.gz"
                if len(files_r1) > 1:
                    logger.info(f"[{sample_id}] Merging {len(files_r1)} SE files into {target_se.name}")
                    self._merge_files(files_r1, target_se)
                else:
                    logger.info(f"[{sample_id}] Creating symlink for {target_se.name}")
                    self._link_file(files_r1[0], target_se)
                sample_info.fastq_1 = target_se
                sample_info.layout = Layout.SE
            else:
                logger.warning(f"Sample {sample_id} has no FASTQ files, skipping.")
                continue

            logger.info(f"Sample {sample_id} layout inferred as: {sample_info.layout}")

        logger.info(f"Successfully processed {len(self.samples_dict)} samples.")


    def _merge_files(self, files: List[Path], out: Path):
        if out.exists():
            logger.info(f"[SKIP] Merged file already exists: {out}")
            return
        logger.info(f"[MERGE] Creating {out} from {len(files)} files")
        with open(out, "wb") as w:
            for f in sorted(files):
                logger.info(f"  -> Merging file: {f}")
                with open(f, "rb") as r:
                    shutil.copyfileobj(r, w) # stream copy to handle large files efficiently

    def _link_file(self, src: Path, dst: Path):
        if dst.is_symlink():
            if dst.resolve() == src.resolve():
                logger.info(f"[SKIP] Link already correct: {dst}")
                return
            dst.unlink()

        elif dst.exists():
            raise RuntimeError(f"Destination exists and is not symlink: {dst}")

        os.symlink(src.resolve(), dst)
        logger.info(f"[LINK] {dst} -> {src}")

    def group_pairs_by_organism(
        self, pairs: List[Tuple[str, str]], samples: Dict[str, SampleInfo]
    ) -> Dict[str, List[Tuple[str, str]]]:
        out = defaultdict(list)
        for ctr, exp in pairs:
            org = samples.get(ctr, SampleInfo()).organism or "UNKNOWN"
            out[org].append((ctr, exp))
        return out

    def run(self):
        if self.meta:
            df = self.load_meta()
            self.prepare_fastq_meta(df = df,
                                    outdir = self.outdir,
                                    data_id_col = self.data_id_col,
                                )
            if df[self.design_col].isnull().all():
                logger.info(f"meta {self.design_col} is all none, skip build_design_pairs")
                pairs = []
            else:
                pairs = self.build_design_pairs()
            return self.samples_dict, pairs, str(self.outdir / "raw_fastq")
        elif self.fastq_dir:
            self.prepare_fastq_dir(self.fastq_dir,self.outdir)
            return self.samples_dict, [], str(self.outdir / "raw_fastq")
        else:
            raise ValueError("Either meta or fastq_dir must be provided.")

    
def main():
    parser = argparse.ArgumentParser(description="Metadata Variants Utils")
    parser.add_argument("--meta", help="Path to metadata file (CSV/TSV)")
    parser.add_argument("--outdir", required=True, help="Output directory for processed FASTQ and logs")
    parser.add_argument("--fastq_dir", help="Directory containing FASTQ files (if not specified in meta)")
    parser.add_argument("--log", help="Path to log file (default: stdout)")

    args = parser.parse_args()

    metadataUtils = MetadataUtils(
        meta=args.meta,
        outdir=args.outdir,
        fastq_dir=args.fastq_dir,
        log=args.log
    )
    res = metadataUtils.run()
    return res    
if __name__ == "__main__":
    main()

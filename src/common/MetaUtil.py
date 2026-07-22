import os
import shutil
import re
import logging
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, Union
import argparse
import math
try:
    from type import FastqMode, Layout,MERIPDesign, SampleInfo, DesignPair
except Exception:
    from .type import FastqMode, Layout,MERIPDesign, SampleInfo, DesignPair

logger = logging.getLogger(__name__)

DESIGN_PATTERN = re.compile(r"^(ctr|ctrl|exp)_(.+)$")
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
        fastq_required_cols: set = {"sample_id", "fastq_1", "fastq_2"},
        pacbio_required_cols: set = {"sample_id", "bam", "pbi"},
        data_id_col: str = "data_id",
        design_col: str = "design",
        group_col: str = "group",
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
            - group_col: Column name in metadata that represents sample groups (default: 'group').
        Note:
            - fq_pattern: Glob pattern to identify FASTQ files in fastq_dir (default: '*fq.gz').

        """
        if not meta and not fastq_dir:
            raise ValueError("Either meta or fastq_dir must be provided.")
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.meta = Path(meta) if meta else None
        self.fastq_dir = Path(fastq_dir) if fastq_dir else None
        self.fastq_required_cols = fastq_required_cols
        self.pacbio_required_cols = pacbio_required_cols
        self.data_id_col = data_id_col
        self.design_col = design_col
        self.group_col = group_col
        self.samples_dict = defaultdict(SampleInfo)
        self.raw_fq_dir = self.outdir / "common" / "1_raw_fastq"
        self.raw_fq_dir.mkdir(parents=True, exist_ok=True)

    def load_meta(self, meta:Union[Path,str]) -> pd.DataFrame:
        """
        function: load metadata from meta file, sep can be \t or ,
        """
        with open(meta, "r", encoding="utf-8") as f:
            head = f.read(2048)

        sep = "\t" if head.count("\t") >= head.count(",") else ","
        df = pd.read_csv(meta, sep=sep)

        return df



    def build_design_pairs(
            self
        ) -> List[DesignPair]:
        """
        Determine ctr/exp pairs based on the design stored in self.samples_dict.

        Design format:  ctr_TAG  or  ctrl_TAG  for control,  exp_TAG  for experiment.
        Tags are underscore-delimited token sets.  A control matches an experiment
        when their token sets intersect (i.e. they share at least one token).

        Examples (all produce a pair):
            ctrl_WT       + exp_WT          -> match (token "WT" shared)
            ctrl_WT_KO    + exp_WT          -> match (token "WT" shared)
            ctrl_WT_KO    + exp_KO          -> match (token "KO" shared)
            ctrl_WT_KO    + exp_WT_KO       -> match ("WT" and "KO" shared)

        No match:
            ctrl_WT       + exp_ABC         -> no common token
            ctrl_ABC      + exp_WT          -> no common token

        If multiple control samples share the same tag, only the first is used
        (a warning is logged).

        return a list of DesignPair objects
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
            design_val = str(design_val).strip()
            m  = DESIGN_PATTERN.match(design_val)
            if not m:
                logger.warning(f"Invalid design format for {sample_id}: {design_val}")
                continue
            role, tag = m.groups()
            # normalise role: "ctrl" -> "ctr" for uniform key
            role = "ctr" if role in ("ctr", "ctrl") else "exp"
            groups[tag][role].append(info)

        # Pre-compute token sets for each tag
        ctr_tags = {tag: set(tag.split("_")) for tag, g in groups.items() if "ctr" in g}
        exp_tags = {tag: set(tag.split("_")) for tag, g in groups.items() if "exp" in g}

        pairs = []
        seen = set()  # deduplicate (ctr_sample_id, exp_sample_id)
        for exp_tag, exp_token_set in exp_tags.items():
            best_ctr_sample = None
            for ctr_tag, ctr_token_set in ctr_tags.items():
                if exp_token_set & ctr_token_set:  # non-empty intersection
                    ctr_samples = groups[ctr_tag]["ctr"]
                    if len(ctr_samples) > 1:
                        logger.warning(
                            f"Multiple ctr samples for tag '{ctr_tag}': "
                            f"{[s.sample_id for s in ctr_samples]}. Only using the first one."
                        )
                    best_ctr_sample = ctr_samples[0]
                    break  # take the first matching control
            if best_ctr_sample is None:
                logger.warning(f"No matching control found for exp tag '{exp_tag}'")
                continue
            for exp_sample_info in groups[exp_tag]["exp"]:
                pair_key = (best_ctr_sample.sample_id, exp_sample_info.sample_id)
                if pair_key in seen:
                    continue
                seen.add(pair_key)
                designPair = DesignPair(
                    organism=exp_sample_info.organism,
                    ctr_sample_id=best_ctr_sample.sample_id,
                    exp_sample_id=exp_sample_info.sample_id,
                    exp_group=exp_sample_info.group
                )
                pairs.append(designPair)
        return pairs


    def prepare_fastq_meta(
            self, 
            df: pd.DataFrame,
            sample_id_col:str = 'sample_id', 
            data_id_col:str = 'data_id',
            design_col:str = 'design',
            fastq_r1_col:str = 'fastq_1',
            fastq_r2_col:str = "fastq_2",
            organism_col:str = "organism",
            workflow_col:str = "workflow",
            group_col:str = "group"
            ) -> None:
        """
            data_id represents a unique FASTQ file.
            If the relationship between sample_id and fastq is one-to-one, a symbolic link is created with the filename prefixed by sample_id.
            If the relationship is one-to-many, FASTQ files corresponding to different data_ids are merged and renamed using the sample_id prefix.

            supplement smaple_id,layout,fastq_1 or fastq_2 information
        """

        if data_id_col not in df.columns:
            df[data_id_col] = df[sample_id_col]

        if not self.fastq_required_cols.issubset(df.columns):
            raise ValueError(f"Metadata must contain columns: {self.fastq_required_cols}")


        raw_fq_dir = self.raw_fq_dir

        df_group = df.groupby(sample_id_col)

        for sample_id, df_sample in df_group:
            sample_id = str(sample_id)
            data_ids = df_sample[data_id_col].values
            if len(data_ids) < 1:
                raise ValueError(f"something wrong: {sample_id} have no {data_id_col} meta")
            
            self.samples_dict[sample_id].sample_id = sample_id
            self.samples_dict[sample_id].design = df_sample[design_col].values[0] if design_col in df_sample.columns else ""
            self.samples_dict[sample_id].organism = df_sample[organism_col].values[0] if organism_col in df_sample.columns else "UNKNOWN"
            self.samples_dict[sample_id].workflow = df_sample[workflow_col].values[0] if workflow_col in df_sample.columns else None
            self.samples_dict[sample_id].group = df_sample[group_col].values[0] if group_col in df_sample.columns else None
            if len(data_ids) == 1:
                logger.info(f"Detect the relationship between {sample_id} and {data_ids[0]} is one-to-one")
                origin_r1 = df_sample[fastq_r1_col].values[0]
                origin_r2 = df_sample[fastq_r2_col].values[0] if fastq_r2_col in df_sample.columns else None
                origin_r1 = Path(origin_r1) if os.path.exists(origin_r1) else None
                origin_r2 = Path(origin_r2) if origin_r2 and os.path.exists(origin_r2) else None

                if origin_r1 and origin_r2:
                    logger.info(f"Detect {data_ids[0]} is Paired END")
                    self.samples_dict[sample_id].layout = Layout.PE
                    rename_r1 = raw_fq_dir / sample_id / f"{sample_id}_1.fq.gz"
                    rename_r2 = raw_fq_dir / sample_id /  f"{sample_id}_2.fq.gz"
                    self._link_file(origin_r1,rename_r1)
                    self._link_file(origin_r2,rename_r2)
                    self.samples_dict[sample_id].fastq_1 = rename_r1
                    self.samples_dict[sample_id].fastq_2 = rename_r2
                elif origin_r1:
                    logger.info(f"Detect {data_ids[0]} is Single End")
                    self.samples_dict[sample_id].layout = Layout.SE
                    rename_r1 = raw_fq_dir / sample_id /  f"{sample_id}.single.fq.gz"
                    self._link_file(origin_r1,rename_r1)
                    self.samples_dict[sample_id].fastq_1 = rename_r1
                else:
                    logger.warning(f"{sample_id} have no fastqs, skip it")
                    continue
            else:
                logger.info(f"Detect the relationship between {sample_id} and {data_ids[0]} is one-to-many")
                origin_r1_list = sorted([r for r in df_sample[fastq_r1_col].values if r])
                origin_r2_list = sorted([r for r in df_sample[fastq_r2_col].values if r]) if fastq_r2_col in df_sample.columns else []
                
                origin_r1_list_path = [Path(r1) for r1 in origin_r1_list]
                origin_r2_list_path = [Path(r2) for r2 in origin_r2_list]

                if len(origin_r1_list_path) > 0 and len(origin_r2_list_path) > 0:
                    logger.info(f"Detect the fastq of {sample_id} is Paired END")
                    self.samples_dict[sample_id].layout = Layout.PE
                    merge_rename_r1 = raw_fq_dir / sample_id /  f"{sample_id}_1.fq.gz"
                    merge_rename_r2 = raw_fq_dir / sample_id /  f"{sample_id}_2.fq.gz"
                    self._merge_files(origin_r1_list_path, merge_rename_r1)
                    self._merge_files(origin_r2_list_path, merge_rename_r2)
                    self.samples_dict[sample_id].fastq_1 = merge_rename_r1
                    self.samples_dict[sample_id].fastq_2 = merge_rename_r2
                elif len(origin_r1_list_path) > 0:
                    logger.info(f"Detect the fastq of {sample_id} is Single END")
                    self.samples_dict[sample_id].layout = Layout.SE
                    merge_rename_r1 = raw_fq_dir / sample_id /  f"{sample_id}.single.fq.gz"
                    self._merge_files(origin_r1_list_path, merge_rename_r1)
                    self.samples_dict[sample_id].fastq_1 = merge_rename_r1
                else:
                    logger.warning(f"{sample_id} have no fastqs, skip it")
                    continue

    def prepare_pacbio_meta(self, df: pd.DataFrame, sample_id_col:str = 'sample_id', bam_col:str = 'bam', pbi_col:str = 'pbi') -> None:
        """
        Prepare PacBio BAM metadata. For each sample_id, create a symlink to the BAM file and its PBI index in the raw_fq_dir.
        """
        if not self.pacbio_required_cols.issubset(df.columns):
            raise ValueError(f"Metadata must contain columns: {self.pacbio_required_cols}")
        for sample_id, df_sample in df.groupby(sample_id_col):
            sample_id = str(sample_id)
            bam_path = df_sample[bam_col].values[0]
            pbi_path = df_sample[pbi_col].values[0]

            if not bam_path or not pbi_path:
                logger.warning(f"{sample_id} is missing BAM or PBI path, skipping.")
                continue

            bam_path = Path(bam_path)
            pbi_path = Path(pbi_path)

            if not bam_path.exists() or not pbi_path.exists():
                logger.warning(f"BAM or PBI file for {sample_id} does not exist, skipping.")
                continue

            target_bam = self.raw_fq_dir / sample_id /  f"{sample_id}.bam"
            target_pbi = self.raw_fq_dir / sample_id /  f"{sample_id}.bam.pbi"

            self._link_file(bam_path, target_bam)
            self._link_file(pbi_path, target_pbi)

            self.samples_dict[sample_id].sample_id = sample_id
            self.samples_dict[sample_id].pacbio_bam = target_bam
            self.samples_dict[sample_id].pacbio_pbi = target_pbi
            self.samples_dict[sample_id].layout = Layout.SE  # Treat BAM as SE for downstream processing
            logger.info(f"Prepared PacBio metadata for sample {sample_id}")

    def prepare_fastq_dir(
        self,
        fq_dir: Path,
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

        raw_fq_dir = self.raw_fq_dir

        for sample_id, reads in temp_files.items():
            sample_info = self.samples_dict[sample_id]
            sample_info.sample_id = sample_id

            files_r1 = sorted(reads["fastq_1"])
            files_r2 = sorted(reads["fastq_2"])

            if files_r1 and files_r2:
                # PE
                target_r1 = raw_fq_dir / sample_id / f"{sample_id}_1.fq.gz"
                target_r2 = raw_fq_dir / sample_id / f"{sample_id}_2.fq.gz"
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
                target_se = raw_fq_dir / sample_id /  f"{sample_id}.single.fq.gz"
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
        out.parent.mkdir(exist_ok=True,parents=True)
        with open(out, "wb") as w:
            for f in sorted(files):
                logger.info(f"  -> Merging file: {f}")
                with open(f, "rb") as r:
                    shutil.copyfileobj(r, w) # stream copy to handle large files efficiently

    def _link_file(self, src: Path, dst: Path):
        """
        Function: soft link file
        """
        dst.parent.mkdir(parents=True,exist_ok=True)
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
            df = self.load_meta(self.meta)
            if "bam" in df.columns and "pbi" in df.columns:
                logger.info("Detected BAM/PBI columns in metadata, preparing PacBio metadata")
                self.prepare_pacbio_meta(df = df, sample_id_col = "sample_id", bam_col = "bam", pbi_col = "pbi")
            else:
                self.prepare_fastq_meta(df = df, data_id_col = self.data_id_col)
            if self.design_col not in df.columns or df[self.design_col].isnull().all():
                logger.info(f"meta {self.design_col} is all none, skip build_design_pairs")
                pairs = []
            else:
                pairs = self.build_design_pairs()
            return self.samples_dict, pairs, str(self.raw_fq_dir)
        elif self.fastq_dir:
            self.prepare_fastq_dir(self.fastq_dir)
            return self.samples_dict, [], str(self.raw_fq_dir)
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
        fastq_dir=args.fastq_dir
    )
    res = metadataUtils.run()
    return res    
if __name__ == "__main__":
    main()

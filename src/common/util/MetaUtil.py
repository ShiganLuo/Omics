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
    from type import FastqMode, Layout, SampleInfo, DesignPair, CompareGroupPair, CellrangerInput
    from LogUtil import setup_logger
except Exception:
    from .type import FastqMode, Layout, SampleInfo, DesignPair, CompareGroupPair, CellrangerInput
    from .LogUtil import setup_logger

logger = setup_logger(__name__, level=logging.DEBUG)

DESIGN_PATTERN = re.compile(r"^(ctr|ctrl|exp)_(.+)$")
REP_PATTERN = re.compile(r"^rep\d+$", re.IGNORECASE)
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
        self.cellranger_input_dict: Dict[str, CellrangerInput] = {}
        self.raw_fq_dir = self.outdir / "common" / "1_raw_data"
        self.raw_fq_dir.mkdir(parents=True, exist_ok=True)

    def load_meta(self, meta:Union[Path,str]) -> pd.DataFrame:
        """Load metadata from a TSV or CSV file.

        Auto-detects the delimiter by comparing the counts of tabs and commas
        in the first 2048 bytes of the file.

        Args:
            meta: Path to the metadata file (TSV or CSV).

        Returns:
            pd.DataFrame: Parsed metadata table.
        """
        with open(meta, "r", encoding="utf-8") as f:
            head = f.read(2048)

        sep = "\t" if head.count("\t") >= head.count(",") else ","
        df = pd.read_csv(meta, sep=sep)

        return df

    def build_design_pairs(
            self
        ) -> Tuple[List[DesignPair], List[CompareGroupPair]]:
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

        return a list of DesignPair objects and a list of CompareGroupPair objects
        """
        group_dict: Dict[str, Dict[str, List[SampleInfo]]] = defaultdict(lambda: defaultdict(list))
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
            role, contrast = m.groups()
            # normalise role: "ctrl" -> "ctr" for uniform key
            role = "ctr" if role in ("ctr", "ctrl") else "exp"
            group_dict[contrast][role].append(info)
        # Pre-compute token sets for each tag
        ctr_contrast_tokens = {contrast: set(contrast.split("_")) for contrast, role in group_dict.items() if "ctr" in role}
        exp_contrast_tokens = {contrast: set(contrast.split("_")) for contrast, role in group_dict.items() if "exp" in role}
        logger.debug("group_dict: %s", group_dict)
        sample_pairs = []
        group_pairs = []
        seen = set()  # deduplicate (ctr_sample_id, exp_sample_id)
        for exp_contrast, exp_token_set in exp_contrast_tokens.items():
            best_ctr_sample = None
            best_ctr_contrast = None
            for ctr_contrast, ctr_token_set in ctr_contrast_tokens.items():
                shared_tokens = exp_token_set & ctr_token_set
                if shared_tokens:  # non-empty intersection
                    non_rep_shared = {t for t in shared_tokens if not REP_PATTERN.match(t)}
                    if not non_rep_shared:
                        logger.warning(
                            f"Suspicious match: exp '{exp_contrast}' and ctr '{ctr_contrast}' "
                            f"share only rep-like tokens {shared_tokens}. "
                            "Consider providing a 'group' column or adjusting design format."
                        )
                    ctr_samples = group_dict[ctr_contrast]["ctr"]
                    if len(ctr_samples) > 1:
                        logger.warning(
                            f"Multiple ctr samples for tag '{ctr_contrast}': "
                            f"{[s.sample_id for s in ctr_samples]}. Only using the first one."
                        )
                    best_ctr_sample = ctr_samples[0]
                    best_ctr_contrast = ctr_contrast
                    break  # take the first matching control
            if best_ctr_sample is None or best_ctr_contrast is None:
                logger.warning(f"No matching control found for exp group '{exp_contrast}'")
                continue
            for exp_sample_info in group_dict[exp_contrast]["exp"]:
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
                sample_pairs.append(designPair)
            
            ctr_group_name = group_dict.get(best_ctr_contrast, {}).get("ctr", [None])[0].group if best_ctr_contrast and group_dict.get(best_ctr_contrast, {}).get("ctr", [None])[0] else f"{best_ctr_contrast}_ctr"
            exp_group_name = group_dict.get(exp_contrast, {}).get("exp", [None])[0].group if exp_contrast and group_dict.get(exp_contrast, {}).get("exp", [None])[0] else f"{exp_contrast}_exp"

            ctr_samples = group_dict.get(best_ctr_contrast, {}).get("ctr", [])
            exp_samples = group_dict.get(exp_contrast, {}).get("exp", [])
            compare_group_pair = CompareGroupPair(
                ctr_group_token=best_ctr_contrast,
                exp_group_token=exp_contrast,
                ctr_group_name=ctr_group_name,
                exp_group_name=exp_group_name,
                ctr_sample_ids=[s.sample_id for s in ctr_samples],
                exp_sample_ids=[s.sample_id for s in exp_samples]
            )
            group_pairs.append(compare_group_pair)

        # --- Validation: warn if final results look wrong ---
        if exp_contrast_tokens and not group_pairs:
            logger.warning(
                "Design column has experimental samples but no group_pairs were generated. "
                "Check that control samples (ctr_/ctrl_) exist with matching tokens."
            )
        for gp in group_pairs:
            if gp.ctr_group_name == gp.exp_group_name:
                logger.warning(
                    f"Control and experimental group have the same name: '{gp.ctr_group_name}'. "
                    "group_pairs will collapse into a single key."
                )

        logger.debug("group_pairs: %s", group_pairs)
        logger.debug("sample_pairs: %s", sample_pairs)
        return sample_pairs, group_pairs


    def _collect_raw_files(self) -> List[str]:
        """Collect resolved (non-symlink) input file paths from raw_fq_dir.

        Scans raw_fq_dir/{sample_id}/ for input files and returns
        os.path.realpath() resolved paths for container bind-mounting.
        """
        result: List[str] = []
        if not self.raw_fq_dir.is_dir():
            return result
        for sample_dir in sorted(self.raw_fq_dir.iterdir()):
            if not sample_dir.is_dir():
                continue
            for f in sorted(sample_dir.iterdir()):
                if f.is_file():
                    result.append(str(f.resolve()))
        return result


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
        """Prepare FASTQ metadata from a dataframe with explicit file paths.

        Handles the sample_id ↔ data_id relationship:
          - One-to-one: symlink the FASTQ to raw_fq_dir/{sample_id}/{sample_id}_{1,2}.fq.gz
          - One-to-many: cat-merge multiple data_id FASTQs into a single file per read

        Populates samples_dict with sample_id, layout (PE/SE), fastq_1, fastq_2,
        design, organism, workflow, and group.

        Args:
            df: Metadata dataframe with at least sample_id, fastq_1, fastq_2 columns.
            sample_id_col: Column name for sample identifier.
            data_id_col: Column name for per-FASTQ identifier. Defaults to sample_id if absent.
            design_col: Column name for experimental design (e.g. ctr_age, exp_age).
            fastq_r1_col: Column name for R1 FASTQ path.
            fastq_r2_col: Column name for R2 FASTQ path.
            organism_col: Column name for organism.
            workflow_col: Column name for workflow identifier.
            group_col: Column name for sample group.
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
            if group_col in df_sample.columns:
                self.samples_dict[sample_id].group = df_sample[group_col].values[0]
            elif design_col in df_sample.columns:
                _design_val = str(df_sample[design_col].values[0]).strip()
                _m = DESIGN_PATTERN.match(_design_val)
                if _m:
                    self.samples_dict[sample_id].group = _m.group(2).replace("_", " ")
                else:
                    self.samples_dict[sample_id].group = None
            else:
                self.samples_dict[sample_id].group = None
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
            elif len(data_ids) > 1:
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
            else:
                logger.warning(f"{sample_id} have no fastqs, skip it")
                continue

    def prepare_pacbio_meta(
        self, df: pd.DataFrame,
        sample_id_col: str = 'sample_id',
        bam_col: str = 'bam',
        pbi_col: str = 'pbi',
    ) -> None:
        """Prepare PacBio BAM metadata: symlink BAM + PBI index per sample.

        Args:
            df: Metadata dataframe with sample_id, bam, pbi columns.
            sample_id_col: Column name for sample identifier.
            bam_col: Column name for BAM file path.
            pbi_col: Column name for PBI index file path.

        Raises:
            ValueError: If required columns (sample_id, bam, pbi) are missing.
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

    def prepare_ms_meta(
            self,
            df: pd.DataFrame,
            sample_id_col: str = 'sample_id',
            ms_file_col: str = 'ms_file',
        ) -> None:
        """Prepare Mass Spectrometry metadata: symlink MS files per sample.

        Supports .raw, .mzML, .mgf and other MS file formats.

        Args:
            df: Metadata dataframe with sample_id and ms_file columns.
            sample_id_col: Column name for sample identifier.
            ms_file_col: Column name for MS file path.

        Raises:
            ValueError: If ms_file_col is missing from df.
        """
        if ms_file_col not in df.columns:
            raise ValueError(f"Metadata must contain column: {ms_file_col}")
        for sample_id, df_sample in df.groupby(sample_id_col):
            sample_id = str(sample_id)
            organism = df_sample.get('organism', pd.Series(["UNKNOWN"])).values[0]
            self.samples_dict[sample_id].organism = organism
            ms_file_path = df_sample[ms_file_col].values[0]

            if not ms_file_path:
                logger.warning(f"{sample_id} is missing MS file path, skipping.")
                continue

            ms_file_path = Path(ms_file_path)

            if not ms_file_path.exists():
                logger.warning(f"MS file for {sample_id} does not exist, skipping.")
                continue

            target_ms_file = self.raw_fq_dir / sample_id /  f"{sample_id}{ms_file_path.suffix}"

            self._link_file(ms_file_path, target_ms_file)

            self.samples_dict[sample_id].sample_id = sample_id
            self.samples_dict[sample_id].ms_file = target_ms_file
            logger.info(f"Prepared MS metadata for sample {sample_id}")

    def prepare_scRNAseq_meta(
        self,
        df: pd.DataFrame,
        sample_id_col: str = 'sample_id',
        fastq_dir_col: str = 'fastq_dir',
        sample_prefix_col: str = 'sample_prefix',
        fq_pattern: str = r"_R?([12])_\d+\.f(ast)?q\.gz$",
    ):
        """Prepare scRNA-seq metadata: scan fastq_dir, merge multi-lane R1/R2.

        For each sample:
          1. Find all FASTQ files matching sample_prefix + fq_pattern in fastq_dir
          2. Single-lane: symlink to raw_fq_dir/{sample_id}/{sample_id}_{1,2}.fq.gz
          3. Multi-lane: cat-merge all R1 into one file, all R2 into one file
          4. Populate samples_dict with merged paths, layout=PE, cellranger fastq_dir

        Output structure:
          raw_fq_dir/{sample_id}/{sample_id}_1.fq.gz  (merged R1)
          raw_fq_dir/{sample_id}/{sample_id}_2.fq.gz  (merged R2)

        Args:
            fq_pattern: regex to match FASTQ filenames and capture read number (1 or 2).
                        Default matches Illumina naming: _R1_001.fastq.gz, _R2_002.fastq.gz
        """
        required = {sample_id_col, fastq_dir_col, sample_prefix_col}
        if not required.issubset(df.columns):
            raise ValueError(f"Metadata must contain columns: {required}")

        for _, row in df.iterrows():
            sample_id = str(row[sample_id_col])
            fastq_dir = Path(row[fastq_dir_col])
            sample_prefix = str(row[sample_prefix_col])
            organism = str(row.get("organism", "UNKNOWN"))
            self.samples_dict[sample_id].organism = organism
            if not fastq_dir.exists():
                logger.warning(f"FASTQ directory for {sample_id} does not exist: {fastq_dir}")
                continue

            # Scan for matching FASTQ files
            reads = {"1": [], "2": []}
            for fq_file in sorted(fastq_dir.iterdir()):
                if not fq_file.is_file():
                    continue
                if not fq_file.name.startswith(sample_prefix):
                    continue
                m = re.search(fq_pattern, fq_file.name)
                if m:
                    read_num = m.group(1)
                    reads[read_num].append(fq_file)

            r1_files = reads["1"]
            r2_files = reads["2"]

            if not r1_files:
                logger.warning(f"No R1 FASTQ found for {sample_id} (prefix={sample_prefix}) in {fastq_dir}")
                continue

            sample_dir = self.raw_fq_dir / sample_id
            target_r1 = sample_dir / f"{sample_id}_1.fq.gz"
            target_r2 = sample_dir / f"{sample_id}_2.fq.gz"

            # R1
            if len(r1_files) == 1:
                logger.info(f"[{sample_id}] Linking R1: {r1_files[0].name}")
                self._link_file(r1_files[0], target_r1)
            else:
                logger.info(f"[{sample_id}] Merging {len(r1_files)} R1 files")
                self._merge_files(r1_files, target_r1)

            # R2
            if r2_files:
                if len(r2_files) == 1:
                    logger.info(f"[{sample_id}] Linking R2: {r2_files[0].name}")
                    self._link_file(r2_files[0], target_r2)
                else:
                    logger.info(f"[{sample_id}] Merging {len(r2_files)} R2 files")
                    self._merge_files(r2_files, target_r2)
                self.samples_dict[sample_id].layout = Layout.PE
            else:
                logger.info(f"[{sample_id}] Single-end (no R2 found)")
                self.samples_dict[sample_id].layout = Layout.SE

            self.samples_dict[sample_id].sample_id = sample_id
            self.samples_dict[sample_id].fastq_1 = target_r1
            if r2_files:
                self.samples_dict[sample_id].fastq_2 = target_r2
            # Keep original dir/prefix for Cell Ranger (--fastqs / --sample)
            self.samples_dict[sample_id].fastq_dir = fastq_dir
            self.samples_dict[sample_id].sample_prefix = sample_prefix

            # Also store design/group/organism/tissue if present
            for col in ("design", "group", "organism", "tissue"):
                if col in df.columns:
                    setattr(self.samples_dict[sample_id], col, row[col])

            # Build cellranger_input for this sample
            self.cellranger_input_dict[sample_id] = CellrangerInput(
                fastq_dir=str(fastq_dir),
                sample_prefix=sample_prefix,
            )

            logger.info(f"[{sample_id}] R1={len(r1_files)} file(s), R2={len(r2_files)} file(s), "
                        f"layout={self.samples_dict[sample_id].layout}")

    def prepare_fastq_dir(
        self,
        fq_dir: Path,
        fq_pattern: str = r"\.f(ast)?q.gz$"
    ) -> None:
        """Auto-detect FASTQ files in a directory, merge multi-lane, populate samples_dict.

        Scans fq_dir recursively for FASTQ files matching fq_pattern.
        Groups files by sample_id (extracted from filename) and read number
        (R1/R2 detected via _R1/_R2 or _1/_2 suffixes).

        Output naming:
          - PE: raw_fq_dir/{sample_id}/{sample_id}_1.fq.gz, {sample_id}_2.fq.gz
          - SE: raw_fq_dir/{sample_id}/{sample_id}.single.fq.gz

        Multi-lane files for the same sample are cat-merged automatically.

        Args:
            fq_dir: Directory containing FASTQ files (searched recursively).
            fq_pattern: Regex pattern to identify FASTQ files (default: *.fastq.gz / *.fq.gz).
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
        """Concatenate multiple files into a single output file (binary stream copy).

        Skips if the output file already exists. Files are sorted before merging
        to ensure deterministic output.

        Args:
            files: List of input file paths to concatenate.
            out: Output file path.
        """
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
        """Create a symbolic link from dst to src.

        If dst is already a symlink pointing to src, this is a no-op.
        If dst is a symlink pointing elsewhere, it is replaced.
        If dst exists as a regular file, a RuntimeError is raised.

        Args:
            src: Source file path (resolved to absolute before linking).
            dst: Destination symlink path.

        Raises:
            RuntimeError: If dst exists and is not a symlink.
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
        """Group (ctr_sample, exp_sample) pairs by organism.

        Args:
            pairs: List of (control_sample_id, experiment_sample_id) tuples.
            samples: Dict mapping sample_id to SampleInfo.

        Returns:
            Dict mapping organism name to list of (ctr, exp) pairs.
        """
        out = defaultdict(list)
        for ctr, exp in pairs:
            org = samples.get(ctr, SampleInfo()).organism or "UNKNOWN"
            out[org].append((ctr, exp))
        return out

    def run(self):
        """Main entry point: auto-detect metadata type and prepare samples.

        Dispatches to the appropriate prepare_* method based on detected columns:
          - bam + pbi columns → prepare_pacbio_meta
          - ms_file column → prepare_ms_meta
          - fastq_dir + sample_prefix → prepare_scRNAseq_meta
          - otherwise → prepare_fastq_meta

        After FASTQ preparation, builds design pairs if a design column exists.

        Returns:
            Tuple of (samples_dict, sample_pairs, group_pairs, raw_fq_dir, raw_files).
        """
        if self.meta:
            df = self.load_meta(self.meta)
            if "bam" in df.columns and "pbi" in df.columns:
                logger.info("Detected BAM/PBI columns in metadata, preparing PacBio metadata")
                self.prepare_pacbio_meta(df = df, sample_id_col = "sample_id", bam_col = "bam", pbi_col = "pbi")
            elif "ms_file" in df.columns:
                logger.info("Detected ms_file columns in metadata, preparing MS metadata")
                self.prepare_ms_meta(df = df, sample_id_col = 'sample_id', ms_file_col = 'ms_file')
            elif "fastq_dir" in df.columns and "sample_prefix" in df.columns:
                logger.info("Detected scRNA-seq columns in metadata, preparing scRNA-seq metadata")
                self.prepare_scRNAseq_meta(df = df, sample_id_col = 'sample_id', fastq_dir_col = 'fastq_dir', sample_prefix_col = 'sample_prefix')
            else:
                self.prepare_fastq_meta(df = df, data_id_col = self.data_id_col)
            if self.design_col not in df.columns or df[self.design_col].isnull().all():
                logger.info(f"meta {self.design_col} is all none, skip build_design_pairs")
                sample_pairs, group_pairs = [], []
            else:
                sample_pairs, group_pairs = self.build_design_pairs()
            return self.samples_dict, sample_pairs, group_pairs, str(self.raw_fq_dir), self._collect_raw_files(), self.cellranger_input_dict
        elif self.fastq_dir:
            self.prepare_fastq_dir(self.fastq_dir)
            return self.samples_dict, [], [], str(self.raw_fq_dir), self._collect_raw_files(), {}
        else:
            raise ValueError("Either meta or fastq_dir must be provided.")

    
def main():
    """CLI entry point for MetadataUtils.

    Usage:
        python MetaUtil.py --meta meta.tsv --outdir output/
        python MetaUtil.py --fastq_dir /path/to/fq --outdir output/
    """
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

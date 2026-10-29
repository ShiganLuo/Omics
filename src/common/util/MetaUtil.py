import os
import shutil
import re
import logging
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import Any, List, Tuple, Dict, Optional, Union
import argparse
import math
try:
    from type import FastqMode, Layout, SampleInfo, DesignPair, CompareGroupPair, CellrangerInput
    from LogUtil import setup_logger
    from SmkUtil import resolve_genome
except Exception:
    from .type import FastqMode, Layout, SampleInfo, DesignPair, CompareGroupPair, CellrangerInput
    from .LogUtil import setup_logger
    from .SmkUtil import resolve_genome

logger = setup_logger(__name__, level=logging.INFO)

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
        bam_required_cols: set = {"sample_id", "bam", "bai"},
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
            - pacbio_required_cols: Set of required columns for PacBio files. Default includes 'sample_id', 'bam', 'pbi'.
            - bam_required_cols: Set of required columns for BAM files. Default includes 'sample_id', 'bam'.
            - data_id_col: Column name in metadata that represents unique FASTQ identifiers (default: 'data_id').
            - design_col: sample compare mode
            - group_col: Column name in metadata that represents sample groups (default: 'group').
        Note:
            - fq_pattern: Glob pattern to identify FASTQ files in fastq_dir (default: '*fq.gz').

        """
        if not meta and not fastq_dir:
            raise ValueError("Either meta or fastq_dir must be provided.")
        self.outdir = outdir
        self.meta = Path(meta) if meta else None
        self.fastq_dir = Path(fastq_dir) if fastq_dir else None
        self.fastq_required_cols = fastq_required_cols
        self.pacbio_required_cols = pacbio_required_cols
        self.bam_required_cols = bam_required_cols
        self.data_id_col = data_id_col
        self.design_col = design_col
        self.group_col = group_col
        self.samples_dict = defaultdict(SampleInfo)
        self.cellranger_input_dict: Dict[str, CellrangerInput] = {}
        self.raw_fq_dir = os.path.join(self.outdir, "common", "1_raw_data")
        self.raw_bam_dir = os.path.join(self.outdir, "common", "2_raw_bam")


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
            # Every control whose tokens intersect the experiment's tokens
            # yields a comparison (exp_E14_KSR matches both ctr_E14 and ctr_KSR).
            matched_ctr_contrasts: List[str] = []
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
                    matched_ctr_contrasts.append(ctr_contrast)
            if not matched_ctr_contrasts:
                logger.warning(f"No matching control found for exp group '{exp_contrast}'")
                continue
            for ctr_contrast in matched_ctr_contrasts:
                ctr_samples = group_dict[ctr_contrast]["ctr"]
                for exp_sample_info in group_dict[exp_contrast]["exp"]:
                    pair_key = (ctr_samples[0].sample_id, exp_sample_info.sample_id)
                    if pair_key in seen:
                        continue
                    seen.add(pair_key)
                    designPair = DesignPair(
                        organism=exp_sample_info.organism,
                        ctr_sample_id=ctr_samples[0].sample_id,
                        exp_sample_id=exp_sample_info.sample_id,
                        exp_group=exp_sample_info.group
                    )
                    sample_pairs.append(designPair)

                ctr_samples = group_dict[ctr_contrast]["ctr"]
                exp_samples = group_dict[exp_contrast]["exp"]
                ctr_group_name = ctr_samples[0].group if ctr_samples else f"{ctr_contrast}_ctr"
                exp_group_name = exp_samples[0].group if exp_samples else f"{exp_contrast}_exp"

                compare_group_pair = CompareGroupPair(
                    organism=exp_samples[0].organism if exp_samples else "UNKNOWN",
                    ctr_group_token=ctr_contrast,
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
        """Collect resolved input file paths from raw_fq_dir.

        Scans raw_fq_dir/{sample_id}/ for input files and returns
        resolved paths for container bind-mounting.
        """
        result: List[str] = []
        raw_fq_dir = os.fspath(self.raw_fq_dir)

        if not os.path.isdir(raw_fq_dir):
            return result

        for sample_id in sorted(os.listdir(raw_fq_dir)):
            sample_dir = os.path.join(raw_fq_dir, sample_id)
            if not os.path.isdir(sample_dir):
                continue

            for filename in sorted(os.listdir(sample_dir)):
                filepath = os.path.join(sample_dir, filename)
                if os.path.isfile(filepath):
                    result.append(os.path.realpath(filepath))
        return result


    @staticmethod
    def _parse_design(
            design_raw: Any = None,
            group_raw: Any = None,
        ) -> Tuple[str, Optional[str]]:
        """Resolve one sample's design string and group from raw meta cells.

        Args:
            design_raw: Raw design cell value (e.g. ``"ctr_WT"``, ``"exp_KO"``).
                Missing/NaN values resolve to ``""``.
            group_raw: Raw group cell value. When absent, empty, or NaN, the
                group falls back to the design tag (the text after the
                ctr_/ctrl_/exp_ prefix, underscores replaced by spaces) when
                the design matches DESIGN_PATTERN.

        Returns:
            Tuple of (design, group): the cleaned design string and the group
            name, or ``None`` when no group can be determined.
        """
        design = ""
        if design_raw is not None and not pd.isna(design_raw):
            if isinstance(design_raw, bytes):
                design_raw = design_raw.decode("utf-8")
            design = str(design_raw)
        group: Optional[str] = None
        if group_raw is not None and not pd.isna(group_raw):
            group_str = str(group_raw).strip()
            if group_str:
                group = group_str
        if group is None:
            m = DESIGN_PATTERN.match(design.strip())
            if m:
                group = m.group(2).replace("_", " ")
        return design, group

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
            group_col:str = "group",
            contaminated_organism_col:str = "contaminated_organism"
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
            sample_id = str(sample_id).strip()
            if ' ' in sample_id or '\t' in sample_id:
                raise ValueError(f"sample_id '{sample_id}' contains whitespace, please remove it from meta file")
            data_ids = df_sample[data_id_col].values
            if len(data_ids) < 1:
                raise ValueError(f"something wrong: {sample_id} have no {data_id_col} meta")
            
            self.samples_dict[sample_id].sample_id = sample_id
            design_value, group_value = self._parse_design(
                df_sample[design_col].values[0] if design_col in df_sample.columns else None,
                df_sample[group_col].values[0] if group_col in df_sample.columns else None,
            )
            self.samples_dict[sample_id].design = design_value
            self.samples_dict[sample_id].group = group_value
            host_organism_value = "UNKNOWN"
            try:
                host_organism_value = str(df_sample[organism_col].values[0]) if organism_col in df_sample.columns else "UNKNOWN"
                contaminated_organism_raw = df_sample[contaminated_organism_col].values[0] if contaminated_organism_col in df_sample.columns else None
                self.samples_dict[sample_id].organism = resolve_genome(host_organism_value)
                if contaminated_organism_raw is not None and not pd.isna(contaminated_organism_raw) and str(contaminated_organism_raw).strip():
                    self.samples_dict[sample_id].contaminated_organism = resolve_genome(str(contaminated_organism_raw))
                else:
                    self.samples_dict[sample_id].contaminated_organism = None 
            except ValueError as e:
                logger.error(f"Failed to resolve genome for organism '{host_organism_value}' in sample '{sample_id}': {e}")
                # some pipeline may not rely on organism, so we can set it to UNKNOWN and continue
                self.samples_dict[sample_id].organism = "UNKNOWN"
                self.samples_dict[sample_id].contaminated_organism = None
            workflow_value = str(df_sample[workflow_col].values[0]) if workflow_col in df_sample.columns else None
            self.samples_dict[sample_id].workflow = workflow_value
            if len(data_ids) == 1:
                logger.info(f"Detect the relationship between {sample_id} and {data_ids[0]} is one-to-one")
                origin_r1 = str(df_sample[fastq_r1_col].values[0]) if fastq_r1_col in df_sample.columns else None
                origin_r2 = str(df_sample[fastq_r2_col].values[0]) if fastq_r2_col in df_sample.columns else None
                origin_r1 = origin_r1 if origin_r1 and os.path.exists(origin_r1) else None
                origin_r2 = origin_r2 if origin_r2 and os.path.exists(origin_r2) else None

                if origin_r1 and origin_r2:
                    logger.info(f"Detect {data_ids[0]} is Paired END")
                    self.samples_dict[sample_id].layout = Layout.PE
                    rename_r1 = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}_1.fq.gz")
                    rename_r2 = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}_2.fq.gz")
                    self._link_file(origin_r1,rename_r1)
                    self._link_file(origin_r2,rename_r2)
                    self.samples_dict[sample_id].fastq_1 = rename_r1
                    self.samples_dict[sample_id].fastq_2 = rename_r2
                elif origin_r1:
                    logger.info(f"Detect {data_ids[0]} is Single End")
                    self.samples_dict[sample_id].layout = Layout.SE
                    rename_r1 = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}.single.fq.gz")
                    self._link_file(origin_r1,rename_r1)
                    self.samples_dict[sample_id].fastq_1 = rename_r1
                else:
                    logger.warning(f"{sample_id} have no fastqs, skip it")
                    continue
            elif len(data_ids) > 1:
                logger.info(f"Detect the relationship between {sample_id} and {data_ids[0]} is one-to-many")
                origin_r1_list = sorted([r for r in df_sample[fastq_r1_col].values if pd.notna(r)])
                origin_r2_list = sorted([r for r in df_sample[fastq_r2_col].values if pd.notna(r)]) if fastq_r2_col in df_sample.columns else []
                
                origin_r1_list_path = [r1 for r1 in origin_r1_list]
                origin_r2_list_path = [r2 for r2 in origin_r2_list]

                if len(origin_r1_list_path) > 0 and len(origin_r2_list_path) > 0:
                    logger.info(f"Detect the fastq of {sample_id} is Paired END")
                    self.samples_dict[sample_id].layout = Layout.PE
                    merge_rename_r1 = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}_1.fq.gz")
                    merge_rename_r2 = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}_2.fq.gz")
                    self._merge_files(origin_r1_list_path, merge_rename_r1)
                    self._merge_files(origin_r2_list_path, merge_rename_r2)
                    self.samples_dict[sample_id].fastq_1 = merge_rename_r1
                    self.samples_dict[sample_id].fastq_2 = merge_rename_r2
                elif len(origin_r1_list_path) > 0:
                    logger.info(f"Detect the fastq of {sample_id} is Single END")
                    self.samples_dict[sample_id].layout = Layout.SE
                    merge_rename_r1 = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}.single.fq.gz")
                    self._merge_files(origin_r1_list_path, merge_rename_r1)
                    self.samples_dict[sample_id].fastq_1 = merge_rename_r1
                else:
                    logger.warning(f"{sample_id} have no fastqs, skip it")
                    continue
            else:
                logger.warning(f"{sample_id} have no fastqs, skip it")
                continue

    def prepare_bam_meta(
        self, df: pd.DataFrame,
        sample_id_col: str = 'sample_id',
        bam_col: str = 'bam',
        pbi_col: str = 'pbi',
        bai_col: str = 'bai',
        design_col: str = 'design',
        group_col: str = 'group'
    ) -> None:
        """Prepare BAM metadata: symlink BAM files per sample. may be used for PacBio CCS data or other BAM-based workflows.
        Args:
            df: Metadata dataframe with sample_id, bam, pbi columns.
            sample_id_col: Column name for sample identifier.
            bam_col: Column name for BAM file path.
            pbi_col: Column name for PBI index file path.

            bai_col: Column name for BAI index file path.
            design_col: Column name for experimental design (e.g. ctrl_WT, exp_KO).
            group_col: Column name for sample group.

        Raises:
            ValueError: If required columns (sample_id, bam, pbi) are missing.
        """
        if self.pacbio_required_cols.issubset(df.columns):
            for sample_id, df_sample in df.groupby(sample_id_col):
                sample_id = str(sample_id).strip()
                if ' ' in sample_id or '\t' in sample_id:
                    raise ValueError(f"sample_id '{sample_id}' contains whitespace, please remove it from meta file")
                bam_path = str(df_sample[bam_col].values[0]) if bam_col in df_sample.columns else None
                pbi_path = str(df_sample[pbi_col].values[0]) if pbi_col in df_sample.columns else None
                organism = str(df_sample.get('organism', pd.Series(["UNKNOWN"])).values[0])
                tissue = str(df_sample.get('tissue', pd.Series(["UNKNOWN"])).values[0])
                self.samples_dict[sample_id].tissue = tissue
                design_value, group_value = self._parse_design(
                    df_sample[design_col].values[0] if design_col in df_sample.columns else None,
                    df_sample[group_col].values[0] if group_col in df_sample.columns else None,
                )
                self.samples_dict[sample_id].design = design_value
                self.samples_dict[sample_id].group = group_value
                try:
                    self.samples_dict[sample_id].organism = resolve_genome(organism)
                except ValueError as e:
                    logger.error(f"Failed to resolve genome for organism '{organism}' in sample '{sample_id}': {e}")
                    # some pipeline may not rely on organism, so we can set it to UNKNOWN and continue
                    self.samples_dict[sample_id].organism = "UNKNOWN"
                if not bam_path or not os.path.exists(bam_path) or not pbi_path or not os.path.exists(pbi_path):
                    logger.warning(f"{sample_id} is missing BAM path or PBI path, skipping.")
                    continue
                else:
                    logger.debug(f"Preparing PacBio metadata for sample {sample_id}: BAM={bam_path}, PBI={pbi_path}")
  
                    target_bam = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}.bam")
                    target_pbi = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}.bam.pbi")

                    self._link_file(bam_path, target_bam)
                    self._link_file(pbi_path, target_pbi)

                    self.samples_dict[sample_id].sample_id = sample_id
                    self.samples_dict[sample_id].pacbio_bam = target_bam
                    self.samples_dict[sample_id].pacbio_pbi = target_pbi
                    self.samples_dict[sample_id].layout = Layout.SE  # Treat BAM as SE for downstream processing
        elif self.bam_required_cols.issubset(df.columns):
            for sample_id, df_sample in df.groupby(sample_id_col):
                sample_id = str(sample_id).strip()
                if ' ' in sample_id or '\t' in sample_id:
                    raise ValueError(f"sample_id '{sample_id}' contains whitespace, please remove it from meta file")
                bam_path = str(df_sample[bam_col].values[0]) if bam_col in df_sample.columns else None
                bai_path = str(df_sample[bai_col].values[0]) if bai_col in df_sample.columns else None
                organism = str(df_sample.get('organism', pd.Series(["UNKNOWN"])).values[0])
                tissue = str(df_sample.get('tissue', pd.Series(["UNKNOWN"])).values[0])
                cellranger_dir = str(df_sample.get('cellranger_dir', pd.Series([None])).values[0])
                self.samples_dict[sample_id].tissue = tissue
                self.samples_dict[sample_id].cellranger_dir = cellranger_dir
                design_value, group_value = self._parse_design(
                    df_sample[design_col].values[0] if design_col in df_sample.columns else None,
                    df_sample[group_col].values[0] if group_col in df_sample.columns else None,
                )
                self.samples_dict[sample_id].design = design_value
                self.samples_dict[sample_id].group = group_value
                try:
                    self.samples_dict[sample_id].organism = resolve_genome(organism)
                except ValueError as e:
                    logger.error(f"Failed to resolve genome for organism '{organism}' in sample '{sample_id}': {e}")
                    # some pipeline may not rely on organism, so we can set it to UNKNOWN and continue
                    self.samples_dict[sample_id].organism = "UNKNOWN"
                if not bam_path or not os.path.exists(bam_path) or not bai_path or not os.path.exists(bai_path):
                    logger.warning(f"{sample_id} is missing BAM path or BAI path, skipping.")
                    continue
                else:
                    logger.debug(f"Preparing BAM metadata for sample {sample_id}: BAM={bam_path}")
                    self.samples_dict[sample_id].sample_id = sample_id
                    self.samples_dict[sample_id].layout = Layout.SE  # Treat BAM as SE for downstream processing

                    if cellranger_dir and os.path.exists(cellranger_dir):
                        logger.debug(f"Sample {sample_id} has a cellranger_dir specified: {cellranger_dir}. skipping symlinking of BAM/BAI files for this sample.")
                        target_dir = os.path.join(self.raw_bam_dir, sample_id)
                        self._link_file(cellranger_dir, target_dir)
                    else:
                        logger.debug(f"Sample {sample_id} does not have a valid cellranger_dir. Proceeding to symlink BAM/BAI files.")
                        target_bam = os.path.join(self.raw_bam_dir, sample_id, f"{sample_id}.bam")
                        target_bai = os.path.join(self.raw_bam_dir, sample_id, f"{sample_id}.bam.bai")
                        self._link_file(bam_path, target_bam)
                        self._link_file(bai_path, target_bai)
                        self.samples_dict[sample_id].bam = target_bam
                        self.samples_dict[sample_id].bai = target_bai
        else:
            raise ValueError(f"Metadata must contain either columns: {self.pacbio_required_cols} or {self.bam_required_cols}")

    def prepare_ms_meta(
            self,
            df: pd.DataFrame,
            sample_id_col: str = 'sample_id',
            ms_file_col: str = 'ms_file',
            design_col: str = 'design',
            group_col: str = 'group',
        ) -> None:
        """Prepare Mass Spectrometry metadata: symlink MS files per sample.

        Supports .raw, .mzML, .mgf and other MS file formats.

        Args:
            df: Metadata dataframe with sample_id and ms_file columns.
            sample_id_col: Column name for sample identifier.
            ms_file_col: Column name for MS file path.
            design_col: Column name for experimental design (e.g. ctrl_WT, exp_KO).
            group_col: Column name for sample group.

        Raises:
            ValueError: If ms_file_col is missing from df.
        """
        if ms_file_col not in df.columns:
            raise ValueError(f"Metadata must contain column: {ms_file_col}")
        for sample_id, df_sample in df.groupby(sample_id_col):
            sample_id = str(sample_id).strip()
            if ' ' in sample_id or '\t' in sample_id:
                raise ValueError(f"sample_id '{sample_id}' contains whitespace, please remove it from meta file")
            host_organism_value = str(df_sample.get('organism', pd.Series(["UNKNOWN"])).values[0])
            try:
                organism = resolve_genome(host_organism_value)
            except ValueError as e:
                logger.error(f"Failed to resolve genome for organism '{host_organism_value}' in sample '{sample_id}': {e}")
                # some pipeline may not rely on organism, so we can set it to UNKNOWN and continue
                organism = "UNKNOWN"
            self.samples_dict[sample_id].organism = organism
            design_value, group_value = self._parse_design(
                df_sample[design_col].values[0] if design_col in df_sample.columns else None,
                df_sample[group_col].values[0] if group_col in df_sample.columns else None,
            )
            self.samples_dict[sample_id].design = design_value
            self.samples_dict[sample_id].group = group_value
            ms_file_path = str(df_sample[ms_file_col].values[0]) if ms_file_col in df_sample.columns else None

            if not ms_file_path or not os.path.exists(ms_file_path):
                logger.warning(f"{sample_id} is missing MS file path or file does not exist, skipping.")
                continue
            suffix = os.path.splitext(ms_file_path)[1]
            target_ms_file = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}{suffix}")

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
            sample_id = str(row[sample_id_col]).strip()
            if ' ' in sample_id or '\t' in sample_id:
                raise ValueError(f"sample_id '{sample_id}' contains whitespace, please remove it from meta file")
            fastq_dir = Path(row[fastq_dir_col])
            sample_prefix = str(row[sample_prefix_col])
            try:
                organism = resolve_genome(str(row.get("organism", "UNKNOWN")))
            except ValueError as e:
                logger.error(f"Failed to resolve genome for organism '{str(row.get('organism', 'UNKNOWN'))}' in sample '{sample_id}': {e}")
                organism = "UNKNOWN"
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

            sample_dir = os.path.join(self.raw_fq_dir, sample_id)
            target_r1 = os.path.join(sample_dir, f"{sample_id}_1.fq.gz")
            target_r2 = os.path.join(sample_dir, f"{sample_id}_2.fq.gz")

            # R1
            if len(r1_files) == 1:
                logger.debug(f"[{sample_id}] Linking R1: {r1_files[0].name}")
                self._link_file(r1_files[0], target_r1)
            else:
                logger.debug(f"[{sample_id}] Merging {len(r1_files)} R1 files")
                self._merge_files(r1_files, target_r1)

            # R2
            if r2_files:
                if len(r2_files) == 1:
                    logger.debug(f"[{sample_id}] Linking R2: {r2_files[0].name}")
                    self._link_file(r2_files[0], target_r2)
                else:
                    logger.info(f"[{sample_id}] Merging {len(r2_files)} R2 files")
                    self._merge_files(r2_files, target_r2)
                self.samples_dict[sample_id].layout = Layout.PE
            else:
                logger.debug(f"[{sample_id}] Single-end (no R2 found)")
                self.samples_dict[sample_id].layout = Layout.SE

            self.samples_dict[sample_id].sample_id = sample_id
            self.samples_dict[sample_id].fastq_1 = target_r1
            if r2_files:
                self.samples_dict[sample_id].fastq_2 = target_r2
            # Keep original dir/prefix for Cell Ranger (--fastqs / --sample)
            self.samples_dict[sample_id].fastq_dir = str(fastq_dir)
            self.samples_dict[sample_id].sample_prefix = sample_prefix

            # Also store organism/tissue if present
            for col in ("organism", "tissue"):
                if col in df.columns:
                    setattr(self.samples_dict[sample_id], col, row[col])
            self.samples_dict[sample_id].design, self.samples_dict[sample_id].group = self._parse_design(
                row.get(self.design_col), row.get("group"))

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
                target_r1 = os.path.join(raw_fq_dir, sample_id, f"{sample_id}_1.fq.gz")
                target_r2 = os.path.join(raw_fq_dir, sample_id, f"{sample_id}_2.fq.gz")
                if len(files_r1) > 1:
                    logger.info(f"[{sample_id}] Merging {len(files_r1)} R1 files into {os.path.basename(target_r1)}")
                    self._merge_files(files_r1, target_r1)
                else:
                    logger.info(f"[{sample_id}] Creating symlink for {os.path.basename(target_r1)}")
                    self._link_file(files_r1[0], target_r1)
                if len(files_r2) > 1:
                    logger.info(f"[{sample_id}] Merging {len(files_r2)} R2 files into {os.path.basename(target_r2)}")
                    self._merge_files(files_r2, target_r2)
                else:
                    logger.info(f"[{sample_id}] Creating symlink for {os.path.basename(target_r2)}")
                    self._link_file(files_r2[0], target_r2)
                sample_info.fastq_1 = target_r1
                sample_info.fastq_2 = target_r2
                sample_info.layout = Layout.PE
            elif files_r1:
                # SE
                target_se = os.path.join(self.raw_fq_dir, sample_id, f"{sample_id}.single.fq.gz")
                if len(files_r1) > 1:
                    logger.info(f"[{sample_id}] Merging {len(files_r1)} SE files into {os.path.basename(target_se)}")
                    self._merge_files(files_r1, target_se)
                else:
                    logger.info(f"[{sample_id}] Creating symlink for {os.path.basename(target_se)}")
                    self._link_file(files_r1[0], target_se)
                sample_info.fastq_1 = target_se
                sample_info.layout = Layout.SE
            else:
                logger.warning(f"Sample {sample_id} has no FASTQ files, skipping.")
                continue

            logger.info(f"Sample {sample_id} layout inferred as: {sample_info.layout}")

        logger.info(f"Successfully processed {len(self.samples_dict)} samples.")


    def _merge_files(self, files: List[str], out: str) -> None:
        """Concatenate multiple files into a single output file.

        Skips if the output file already exists. Files are sorted before merging
        to ensure deterministic output.

        Args:
            files: List of input file paths to concatenate.
            out: Output file path.
        """
        if os.path.exists(out):
            logger.debug(f"[SKIP] Merged file already exists: {out}")
            return

        logger.info(f"[MERGE] Creating {out} from {len(files)} files")
        os.makedirs(os.path.dirname(out), exist_ok=True)

        with open(out, "wb") as w:
            for f in sorted(files):
                logger.info(f"  -> Merging file: {f}")
                with open(f, "rb") as r:
                    shutil.copyfileobj(r, w)

    def _link_file(self, src: str, dst: str) -> None:
        """Create a symbolic link to a file or directory.

        If dst is already a symlink pointing to src, this is a no-op.
        If dst is a symlink pointing elsewhere, it is replaced.
        If dst exists and is not a symlink, a RuntimeError is raised.
        """
        src = os.path.realpath(src)
        dst = os.fspath(dst)

        os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)

        if os.path.islink(dst):
            if os.path.realpath(dst) == src:
                logger.debug(f"[SKIP] Link already correct: {dst}")
                return
            os.unlink(dst)
        elif os.path.exists(dst):
            raise RuntimeError(f"Destination exists and is not symlink: {dst}")

        os.symlink(src, dst, target_is_directory=os.path.isdir(src))
        logger.info(f"LINK: {dst} -> {src}")

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
          - bam + pbi column → prepare_bam_meta
          - bam + bai column → prepare_bam_meta
          - ms_file column → prepare_ms_meta
          - fastq_dir + sample_prefix → prepare_scRNAseq_meta
          - otherwise → prepare_fastq_meta

        After FASTQ preparation, builds design pairs if a design column exists.

        Returns:
            Tuple of (samples_dict, sample_pairs, group_pairs, raw_fq_dir, raw_files).
        """
        if self.meta:
            df = self.load_meta(self.meta)
            if "fastq_1" in df.columns or "fastq_2" in df.columns:
                logger.info("Detected fastq_1/fastq_2 columns in metadata, preparing FASTQ metadata")
                self.prepare_fastq_meta(df = df, sample_id_col = "sample_id", data_id_col = self.data_id_col)
            elif "fastq_dir" in df.columns and "sample_prefix" in df.columns:
                logger.info("Detected fastq_dir/sample_prefix columns in metadata, preparing scRNA-seq metadata")
                self.prepare_scRNAseq_meta(df = df, sample_id_col = "sample_id", fastq_dir_col = "fastq_dir", sample_prefix_col = "sample_prefix")
            elif "bam" in df.columns:
                logger.info("Detected BAM column in metadata, preparing BAM metadata")
                self.prepare_bam_meta(df = df, sample_id_col = "sample_id", bam_col = "bam", pbi_col = "pbi")
            elif "ms_file" in df.columns:
                logger.info("Detected ms_file columns in metadata, preparing MS metadata")
                self.prepare_ms_meta(df = df, sample_id_col = 'sample_id', ms_file_col = 'ms_file')
            else:
                logger.info("No recognized FASTQ, BAM/PBI, or MS columns found in metadata, preparing FASTQ metadata by default")
                raise ValueError("Metadata must contain either fastq_1/fastq_2, fastq_dir/sample_prefix, bam/pbi, or ms_file columns.")

            if self.design_col not in df.columns or df[self.design_col].isnull().all():
                logger.info(f"meta {self.design_col} is all none, skip build_design_pairs")
                sample_pairs, group_pairs = [], []
            else:
                sample_pairs, group_pairs = self.build_design_pairs()
            return self.samples_dict, sample_pairs, group_pairs, self.raw_fq_dir, self._collect_raw_files(), self.cellranger_input_dict
        elif self.fastq_dir:
            self.prepare_fastq_dir(self.fastq_dir)
            return self.samples_dict, [], [], self.raw_fq_dir, self._collect_raw_files(), {}
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

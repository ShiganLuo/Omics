from dataclasses import dataclass
from enum import Enum, unique
from typing import Optional, List, Dict
from pathlib import Path

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
    organism: str = "UNKNOWN"
    contaminated_organism: Optional[str] = None
    layout: Layout = Layout.UNKNOWN
    fastq_1: Optional[str] = None # Path to the first FASTQ file (for SE or PE)
    fastq_2: Optional[str] = None # Path to the second FASTQ file (for PE)
    workflow: Optional[str] = None
    group: Optional[str] = None
    design: Optional[str] = None
    pacbio_bam: Optional[str] = None # PacBio BAM file
    pacbio_pbi: Optional[str] = None # PacBio index file
    bam: Optional[str] = None # BAM file
    bai: Optional[str] = None # BAM index file
    ms_file: Optional[str] = None # Mass Spectrometry file (.raw, .mzML, .mgf, etc.)
    fastq_dir: Optional[str] = None # Directory containing FASTQ files for scRNA-seq
    cellranger_dir: Optional[str] = None # Directory containing Cell Ranger output for scRNA-seq
    sample_prefix: Optional[str] = None # Prefix for FASTQ files in the directory for scRNA-seq
    tissue: Optional[str] = None # Tissue type for scRNA-seq grouping (e.g. PBMC, brain)

@dataclass
class DesignPair:
    organism: str
    ctr_sample_id: str
    exp_sample_id: str
    exp_group: Optional[str] = None

@dataclass
class CompareGroupPair:
    organism: str
    ctr_group_token: str
    exp_group_token: str
    ctr_group_name: str
    exp_group_name: str
    ctr_sample_ids: List[str]
    exp_sample_ids: List[str]

@dataclass
class CellrangerInput:
    fastq_dir: str
    sample_prefix: str


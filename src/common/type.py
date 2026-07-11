from dataclasses import dataclass
from enum import Enum, unique
from typing import Optional
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
    organism: str = ""
    layout: Layout = Layout.UNKNOWN
    fastq_1: Optional[Path] = None
    fastq_2: Optional[Path] = None
    workflow: Optional[str] = None
    group: Optional[str] = None
    design: Optional[str] = None
    pacbio_bam: Optional[Path] = None
    pacbio_pbi: Optional[Path] = None

@dataclass
class DesignPair:
    organism: str
    ctr_sample_id: str
    exp_sample_id: str
    exp_group: Optional[str] = None
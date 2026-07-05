#!/usr/bin/env python3
"""Detect whether a FASTQ file already has UMI extracted into the read header.

UMI extraction tools (umi_tools extract, fumi_tools copy_umi) append the UMI
sequence to the read name with an underscore separator:

    @LH00326:331:22THVWLT3:7:1101:1407:1064_ANGGCTTCCCTG 1:N:0:CAAGCTAG+CGCTATGT
                                             ^^^^^^^^^^^^
                                             UMI (12bp ACGT)

This module provides a function to detect this pattern so that downstream
rules can skip re-extraction and symlink the input directly.
"""

from __future__ import annotations

import gzip
import re
from typing import Optional

# Pattern: underscore followed by one or more ACGT/N characters at the end of
# the read name (before the first whitespace).  The ACGT/N constraint avoids
# false positives on instrument IDs that may contain underscores.
_UMI_SUFFIX_RE = re.compile(r"_([ACGTNacgtn]+)$")


def has_umi_in_header(
    fastq_path: str,
    umi_length: Optional[int] = None,
) -> bool:
    """Check whether the first read in *fastq_path* already carries a UMI.

    Parameters
    ----------
    fastq_path : str
        Path to a (possibly gzipped) FASTQ file.
    umi_length : int or None
        If given, the UMI segment must be exactly this many bases long.
        If ``None``, any ``_<ACGTN>+`` suffix is accepted.

    Returns
    -------
    bool
        ``True`` when a UMI is detected in the read name.
    """
    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as fh:
        first_line = fh.readline().rstrip("\n")

    if not first_line.startswith("@"):
        return False

    # Read name is everything before the first whitespace
    read_name = first_line.split()[0]  # includes leading '@'
    name_body = read_name[1:]  # strip '@'

    m = _UMI_SUFFIX_RE.search(name_body)
    if m is None:
        return False

    umi_seq = m.group(1)
    if umi_length is not None and len(umi_seq) != umi_length:
        return False

    return True

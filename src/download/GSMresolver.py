from typing import Dict, List, Optional, Set
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import csv
from common.SepUtil import detect_delimiter
import logging

logger = logging.getLogger(__name__)

class GSMResolver:
    """
    Resolve GSM accession IDs from multiple input sources with clear priority.

    This helper class unifies GSM collection logic so downstream download/extract
    workflows can consume a single normalized GSM list regardless of where IDs
    come from.

    Supported sources:
    1. Explicit GSM list from CLI (`gsm`)
    2. Structured or plain-text file (`gsm_file`)
    3. Existing HTML filenames in a directory (`html_dir`)
    4. Piped stdin input (only when no IDs were found above)

    Notes
    -----
    - IDs are de-duplicated via a set.
    - Final output is always sorted for deterministic behavior.
    - Column matching in tabular files is case-insensitive.
    """
    def __init__(
        self,
        gsm: Optional[List[str]] = None,
        gsm_file: Optional[str] = None,
        html_dir: Optional[str] = None,
        id_column: str = "GSM"
    ):
        self.gsm = gsm
        self.gsm_file = gsm_file
        self.html_dir = html_dir
        self.id_column = id_column.lower()

    def resolve(self) -> List[str]:
        """
        Resolve GSM IDs from all configured sources.

        Resolution order:
        1. Add IDs from `self.gsm` when provided;
        2. Merge IDs parsed from `self.gsm_file`;
        3. Merge IDs inferred from `self.html_dir`;
        4. If still empty and stdin is piped, read from stdin.

        Returns
        -------
        List[str]
            Sorted GSM accession list with duplicates removed.

        Raises
        ------
        ValueError
            If no GSM can be resolved from any available source.
        """
        gsms: Set[str] = set()

        if self.gsm:
            gsms.update(self.gsm)

        if self.gsm_file:
            gsms.update(self._from_file(self.gsm_file))

        if self.html_dir:
            gsms.update(self._from_html_dir(self.html_dir))

        if not gsms and not sys.stdin.isatty():
            gsms.update(self._from_stdin())

        if not gsms:
            raise ValueError("No GSM found from any source")

        return sorted(gsms)

    def _from_file(self, path: str) -> Set[str]:
        """
        Parse GSM IDs from a file path.

        Parsing behavior depends on file extension:
        - `.csv` / `.tsv`: use `csv.DictReader`, auto-detect delimiter,
          and read values from the configured ID column.
        - other extensions: treat as plain text, one ID per non-empty line,
          skipping comment lines beginning with `#`.

        Parameters
        ----------
        path : str
            Input file path containing GSM identifiers.

        Returns
        -------
        Set[str]
            Parsed GSM IDs with surrounding whitespace stripped.

        Raises
        ------
        ValueError
            If tabular input does not contain the configured ID column.
        """
        gsms: Set[str] = set()
        ext = os.path.splitext(path)[1].lower()

        if ext in (".csv", ".tsv"):
            delimiter = detect_delimiter(path)
            with open(path, encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter=delimiter)
                logger.info(f"Reading {path} with delimiter {delimiter}")
                col = self._find_id_column(reader.fieldnames)
                if not col:
                    raise ValueError(f"No column '{self.id_column}' in {path}")
                for row in reader:
                    val = row.get(col)
                    if val:
                        gsms.add(val.strip())
        else:
            with open(path, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        gsms.add(line)

        return gsms

    def _from_html_dir(self, html_dir: str) -> Set[str]:
        """
        Infer GSM IDs from cached HTML filenames.

        The method scans the directory for files ending with `.html` and
        returns basename values with the extension removed.

        Parameters
        ----------
        html_dir : str
            Directory containing downloaded HTML files.

        Returns
        -------
        Set[str]
            Identifier set extracted from HTML filenames.
        """
        return {
            fn.replace(".html", "")
            for fn in os.listdir(html_dir)
            if fn.endswith(".html")
        }

    def _from_stdin(self) -> Set[str]:
        """
        Read GSM IDs from standard input stream.

        Empty lines are ignored and trailing whitespace is stripped.

        Returns
        -------
        Set[str]
            GSM identifiers collected from piped stdin input.
        """
        return {line.strip() for line in sys.stdin if line.strip()}

    def _find_id_column(self, columns):
        """
        Find the matching ID column name from table headers.

        Matching is case-insensitive against `self.id_column` while returning
        the original header name for correct dictionary lookup.

        Parameters
        ----------
        columns : Iterable[str] | None
            Header names reported by `csv.DictReader`.

        Returns
        -------
        str | None
            The matched original column name, or None if not found.
        """
        for c in columns or []:
            if c.lower() == self.id_column:
                return c
        return None

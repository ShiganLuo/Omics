from typing import Dict, List, Optional, Set
import sys
import os
import csv
class GSMResolver:
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
        gsms: Set[str] = set()
        ext = os.path.splitext(path)[1].lower()

        if ext in (".csv", ".tsv"):
            delimiter = "," if ext == ".csv" else "\t"
            with open(path, encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter=delimiter)
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
        return {
            fn.replace(".html", "")
            for fn in os.listdir(html_dir)
            if fn.endswith(".html")
        }

    def _from_stdin(self) -> Set[str]:
        return {line.strip() for line in sys.stdin if line.strip()}

    def _find_id_column(self, columns):
        for c in columns or []:
            if c.lower() == self.id_column:
                return c
        return None

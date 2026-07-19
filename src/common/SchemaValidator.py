"""Schema validator for Omics workflow config JSONs.

Each workflow has its own schema file: config/<Workflow>.schema.json
The schema mirrors the config structure (top-level fields + genome section).
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List


class SchemaValidator:
    """Validate workflow configs against per-workflow schema files."""

    def __init__(self):
        self._schema: Dict[str, Any] = {}
        self._schema_dir: str = ""

    def load(self, schema_path: str) -> None:
        """Load a schema JSON file."""
        with open(schema_path, "r", encoding="utf-8") as f:
            self._schema = json.load(f)
        self._schema_dir = str(Path(schema_path).parent)

    def load_workflow(self, workflow_name: str) -> None:
        """Load schema for a specific workflow from config/<wf>.schema.json."""
        schema_path = os.path.join(self._schema_dir, f"{workflow_name}.schema.json")
        if not os.path.exists(schema_path):
            raise FileNotFoundError(f"Schema not found: {schema_path}")
        self.load(schema_path)

    @property
    def schema(self) -> Dict[str, Any]:
        if not self._schema:
            raise RuntimeError("Schema not loaded. Call load() first.")
        return self._schema

    def get_path_fields(self) -> Dict[str, dict]:
        """Return all path-type fields from the loaded schema.

        Returns:
            Dict mapping dotted key (e.g. "genome.fasta") to schema entry.
        """
        result = {}
        for key, defn in self.schema.items():
            if key == "genome":
                continue
            if isinstance(defn, dict) and defn.get("path"):
                result[key] = defn

        genome = self.schema.get("genome", {})
        for field_name, defn in genome.items():
            if isinstance(defn, dict) and defn.get("path"):
                result[f"genome.{field_name}"] = defn

        return result

    def validate(self, config: Dict[str, Any]) -> List[str]:
        """Validate config against the loaded schema.

        Returns:
            List of error messages. Empty = valid.
        """
        errors = []

        # Check top-level fields
        for key, defn in self.schema.items():
            if key == "genome":
                continue
            if not isinstance(defn, dict):
                continue
            value = config.get(key)
            if defn.get("required") and (value is None or value == ""):
                errors.append(f"Missing or empty required field: '{key}'")
            if not defn.get("nullable", True) and value is None:
                errors.append(f"Non-nullable field is null: '{key}'")

        # Check genome fields
        genome_schema = self.schema.get("genome", {})
        genome_cfg = config.get("genome", {})
        for field_name, defn in genome_schema.items():
            if not isinstance(defn, dict):
                continue
            value = genome_cfg.get(field_name)
            dotted = f"genome.{field_name}"
            if defn.get("required") and (value is None or value == ""):
                errors.append(f"Missing or empty required field: '{dotted}'")
            if not defn.get("nullable", True) and value is None:
                errors.append(f"Non-nullable field is null: '{dotted}'")

        return errors

    def generate_test_paths(self, test_data: str, genome: str) -> Dict[str, str]:
        """Generate test file paths for all path-type fields.

        Scans all config/*.schema.json files, collects path-type fields,
        creates placeholder files, returns {config_key: abs_path} mapping.
        """
        test_data_path = Path(test_data).resolve()
        ref_dir = test_data_path / "ref"
        index_dir = test_data_path / "index"
        ref_dir.mkdir(parents=True, exist_ok=True)
        index_dir.mkdir(parents=True, exist_ok=True)

        INDEX_MAP = {
            "hisat2": {"exts": [f".{i}.ht2" for i in range(1, 9)]},
            "bwaMem2": {"exts": [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]},
            "bowtie2": {"exts": [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]},
            "bowtie2_for_rRNA": {"exts": [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]},
        }

        mapping: Dict[str, str] = {}
        config_dir = Path(self._schema_dir)

        for schema_file in sorted(config_dir.glob("*.schema.json")):
            with open(schema_file, "r", encoding="utf-8") as f:
                schema = json.load(f)

            # Top-level path fields
            for key, defn in schema.items():
                if key == "genome" or not isinstance(defn, dict):
                    continue
                if defn.get("path") and key not in mapping:
                    path_type = defn["path"]
                    if path_type == "dir":
                        p = test_data_path / key
                        p.mkdir(parents=True, exist_ok=True)
                    else:
                        p = ref_dir / f"{genome}.{key}"
                        p.parent.mkdir(parents=True, exist_ok=True)
                        p.touch(exist_ok=True)
                    mapping[key] = str(p)

            # Genome path fields
            genome_schema = schema.get("genome", {})
            for field_name, defn in genome_schema.items():
                if not isinstance(defn, dict) or not defn.get("path"):
                    continue
                config_key = f"genome.{field_name}"
                if config_key in mapping:
                    continue

                path_type = defn["path"]
                if path_type == "file":
                    # Determine filename
                    if "smallrna" in field_name:
                        name = "smallrna"
                    elif "rRNA" in field_name:
                        name = "rRNA"
                    elif "access" in field_name:
                        name = "access"
                    elif "repeat" in field_name:
                        name = "repeat"
                    elif "decoy" in field_name:
                        name = "decoy"
                    elif "TE" in field_name:
                        name = f"GRCm39.{field_name}"
                    else:
                        name = genome
                    p = ref_dir / name
                    p.parent.mkdir(parents=True, exist_ok=True)
                    p.touch(exist_ok=True)
                    mapping[config_key] = str(p)

                elif path_type == "dir":
                    if "smallrna" in field_name:
                        d = index_dir / "star" / "smallrna"
                    elif "star" in field_name:
                        d = index_dir / "star" / genome
                    else:
                        d = index_dir / field_name.replace("_dir", "")
                    d.mkdir(parents=True, exist_ok=True)
                    for fname in ["Genome", "SA", "SAindex"]:
                        (d / fname).touch(exist_ok=True)
                    mapping[config_key] = str(d)

                elif path_type == "prefix":
                    # Determine tool from field name
                    tool = field_name.replace("_index_prefix", "").replace("_for_rRNA", "")
                    if tool in INDEX_MAP:
                        exts = INDEX_MAP[tool]["exts"]
                    else:
                        exts = [".1.ht2"]
                    pfx_dir = index_dir / tool / genome
                    pfx_dir.mkdir(parents=True, exist_ok=True)
                    pfx = pfx_dir / genome
                    for ext in exts:
                        (pfx_dir / f"{genome}{ext}").touch(exist_ok=True)
                    mapping[config_key] = str(pfx)

        # Write chrom.sizes
        chrom_sizes = ref_dir / "chrom.sizes"
        if not chrom_sizes.exists():
            with open(chrom_sizes, "w") as f:
                f.write("chr1\t195471971\nchr2\t182113224\nchr3\t159970021\n")

        return mapping

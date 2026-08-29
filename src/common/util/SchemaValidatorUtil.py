"""Schema validator for Omics workflow config JSONs.

Each workflow has its own schema file: config/<Workflow>.schema.json
The schema mirrors the config structure (top-level fields + genome section).
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def smart_cast(val):
    """Convert a string value to int, float, or bool if possible.

    Recursively converts list elements. Avoids octal misinterpretation
    for strings starting with '0' (except '0.').

    Args:
        val: Input value (str, list, or other).

    Returns:
        Converted value: bool for 'true'/'false', int, float, or original string.
    """
    if isinstance(val, list):
        return [smart_cast(v) for v in val]
    if isinstance(val, str):
        if val.lower() in {"true", "false"}:
            return val.lower() == "true"
        try:
            if val.startswith("0") and len(val) > 1 and not val.startswith("0."):
                return val  # 避免八进制等
            return int(val)
        except Exception:
            pass
        try:
            return float(val)
        except Exception:
            pass
    return val


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

    # ------------------------------------------------------------------
    # Type resolution & cast
    # ------------------------------------------------------------------

    @staticmethod
    def _is_field_def(obj: Any) -> bool:
        """Check if a dict looks like a field definition (has 'type' key)."""
        return isinstance(obj, dict) and "type" in obj

    @staticmethod
    def _resolve_field(node: dict, key: str) -> Optional[Any]:
        """Resolve *key* inside a schema node.

        Lookup order:
          1. node["properties"][key]   — Style A (RNAseq/scRNAseq)
          2. node["additionalProperties"]["properties"][key] — dynamic keys
          3. node[key]                 — Style B/C (direct fields)

        Returns the child node (dict with 'type' or nested 'properties'),
        or None if not found.
        """
        if not isinstance(node, dict):
            return None

        props = node.get("properties")
        if isinstance(props, dict) and key in props:
            return props[key]

        add_props = node.get("additionalProperties")
        if isinstance(add_props, dict):
            inner_props = add_props.get("properties")
            if isinstance(inner_props, dict) and key in inner_props:
                return inner_props[key]

        if key in node:
            return node[key]

        return None

    def get_field_type(self, dotted_key: str) -> Optional[Dict[str, Any]]:
        """Return the schema definition for *dotted_key* (e.g. ``"fasta"`` or
        ``"Params.star.outFilterMultimapNmax"``).

        Returns the definition dict (with at least a ``"type"`` entry), or
        ``None`` if the key is not covered by the current schema.
        """
        parts = dotted_key.split(".")
        if not parts:
            return None

        # ---- traverse to the parent node --------------------------------
        node = self.schema
        for part in parts[:-1]:
            resolved = self._resolve_field(node, part)
            if resolved is None:
                return None
            # If resolved is a field def with properties, step into it.
            if isinstance(resolved, dict) and "properties" in resolved:
                node = resolved
            elif isinstance(resolved, dict) and "type" in resolved:
                # Typed field with child keys → still a container (e.g. "type":"dict")
                child_keys = {k for k in resolved.keys() if k not in (
                    "type", "required", "nullable", "path", "items", "description",
                    "properties", "additionalProperties"
                )}
                if child_keys:
                    node = resolved
                else:
                    return None  # Pure leaf; can't go deeper.
            else:
                node = resolved

        # ---- resolve the final segment ----------------------------------
        leaf = self._resolve_field(node, parts[-1])
        if leaf is None:
            return None

        # If leaf is itself a field def, return it directly.
        if self._is_field_def(leaf):
            return leaf

        # Leaf is a container (e.g. genome subtree without 'type').
        # We can still extract structural info.
        if isinstance(leaf, dict):
            if "properties" in leaf:
                return {"type": "object", **leaf}
            # Pure subtree — not a typed field.
            return None

        return None

    # ------------------------------------------------------------------
    # Index file extensions (for prefix path validation)
    # ------------------------------------------------------------------

    INDEX_MAP = {
        "hisat2": [f".{i}.ht2" for i in range(1, 9)],
        "bwaMem2": [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"],
        "bowtie2": [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"],
        "bowtie2_for_rRNA": [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"],
    }

    # ------------------------------------------------------------------
    # Cast & validate extra_args
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_field_def(fd: Dict[str, Any]) -> Tuple[str, bool, Optional[Dict], Optional[str]]:
        """Extract (type_str, nullable, items_def, path_type) from a field definition."""
        type_str = fd.get("type", "str")
        nullable = fd.get("nullable", True)
        items = fd.get("items") if isinstance(fd.get("items"), dict) else None
        path_type = fd.get("path") if isinstance(fd.get("path"), str) else None
        return type_str, nullable, items, path_type

    def cast_extra_args(self, extra_args: Dict[str, Any]) -> Tuple[Dict[str, Any], List[str]]:
        """Cast *extra_args* values to their schema-declared types.

        Returns ``(casted_dict, errors)``.  An empty *errors* list means all
        values are valid.  Keys not present in the schema are passed through
        ``smart_cast`` unchanged (no error raised).
        """
        result: Dict[str, Any] = {}
        errors: List[str] = []

        for key, value in extra_args.items():
            field_def = self.get_field_type(key)

            # ---- key not in schema → fall back to smart_cast -------------
            if field_def is None:
                result[key] = smart_cast(value)
                continue

            type_str, nullable, items_def, path_type = self._parse_field_def(field_def)

            # ---- nullable guard -----------------------------------------
            if (value is None or value == "") and not nullable:
                errors.append(f"字段 '{key}' 不允许为空")
                result[key] = value
                continue

            # ---- cast by target type ------------------------------------
            if type_str in ("list", "array"):
                if not isinstance(value, list):
                    value = [value]
                if items_def:
                    item_type = items_def.get("type", "str")
                    casted_items: List[Any] = []
                    for i, item in enumerate(value):
                        cv, cerr = self._cast_scalar(
                            item, item_type, f"{key}[{i}]"
                        )
                        casted_items.append(cv)
                        errors.extend(cerr)
                    value = casted_items
                result[key] = value

            elif type_str == "int":
                casted, cerr = self._cast_scalar(value, "int", key)
                result[key] = casted
                errors.extend(cerr)

            elif type_str == "float":
                casted, cerr = self._cast_scalar(value, "float", key)
                result[key] = casted
                errors.extend(cerr)

            elif type_str == "bool":
                casted, cerr = self._cast_scalar(value, "bool", key)
                result[key] = casted
                errors.extend(cerr)

            elif type_str == "str":
                result[key] = value  # CLI strings are always valid

            elif type_str in ("dict", "object"):
                if isinstance(value, str):
                    try:
                        result[key] = json.loads(value)
                    except (json.JSONDecodeError, ValueError):
                        errors.append(
                            f"字段 '{key}' 期望 dict/object，"
                            f"无法解析 JSON: '{value}'"
                        )
                        result[key] = value
                else:
                    result[key] = value

            elif type_str == "null":
                result[key] = value

            else:
                # Unknown type string — pass through
                result[key] = smart_cast(value)

            # ---- path validation ----------------------------------------
            if path_type and result[key] is not None and result[key] != "":
                path_errors = self._validate_path(key, result[key], path_type)
                errors.extend(path_errors)

            # ---- enum validation ----------------------------------------
            allowed = field_def.get("enum")
            if allowed is not None and isinstance(allowed, list):
                enum_errors = self._validate_enum(key, result[key], allowed)
                errors.extend(enum_errors)

            # ---- range validation ---------------------------------------
            # 1. Explicit from schema
            lo = field_def.get("minimum")
            hi = field_def.get("maximum")
            # 2. Inferred from field name if schema doesn't declare
            if lo is None and hi is None and isinstance(result[key], (int, float)):
                inferred = self._infer_range(key)
                if inferred:
                    lo, hi = inferred
            if (lo is not None or hi is not None) and isinstance(result[key], (int, float)):
                range_errors = self._validate_range(key, result[key], lo, hi)
                errors.extend(range_errors)

        return result, errors

    @staticmethod
    def _cast_scalar(
        value: Any, target_type: str, field_name: str
    ) -> Tuple[Any, List[str]]:
        """Cast a single scalar value.  Returns ``(casted_value, errors)``."""
        errors: List[str] = []

        if target_type == "int":
            try:
                return int(value), errors
            except (ValueError, TypeError):
                errors.append(
                    f"字段 '{field_name}' 期望 int，收到 '{value}'"
                )
                return value, errors

        if target_type == "float":
            try:
                return float(value), errors
            except (ValueError, TypeError):
                errors.append(
                    f"字段 '{field_name}' 期望 float，收到 '{value}'"
                )
                return value, errors

        if target_type == "bool":
            _BOOL_MAP = {
                "true": True, "false": False,
                "yes": True,  "no": False,
                "1": True,    "0": False,
            }
            if isinstance(value, str) and value.lower() in _BOOL_MAP:
                return _BOOL_MAP[value.lower()], errors
            if isinstance(value, bool):
                return value, errors
            if isinstance(value, int) and value in (0, 1):
                return bool(value), errors
            errors.append(
                f"字段 '{field_name}' 期望 bool，收到 '{value}'，"
                f"有效值: true/false/yes/no/1/0"
            )
            return value, errors

        # Fallback — shouldn't be reached for known types
        return smart_cast(value), errors

    @classmethod
    def _validate_path(
        cls, field_name: str, value: Any, path_type: str
    ) -> List[str]:
        """Validate that a path value points to an existing file/dir/prefix.

        Args:
            field_name: Dotted key for error messages.
            value: The path value (str or list of str).
            path_type: One of "file", "dir", "prefix".

        Returns:
            List of error messages (empty if valid).
        """
        errors: List[str] = []

        # Handle list of paths
        if isinstance(value, list):
            for i, item in enumerate(value):
                if isinstance(item, str):
                    errors.extend(
                        cls._validate_single_path(f"{field_name}[{i}]", item, path_type)
                    )
            return errors

        if isinstance(value, str):
            return cls._validate_single_path(field_name, value, path_type)

        return errors

    @classmethod
    def _validate_single_path(
        cls, field_name: str, value: str, path_type: str
    ) -> List[str]:
        """Validate a single path value."""
        errors: List[str] = []
        value = value.strip()

        if not value:
            return errors

        if path_type == "file":
            if not os.path.isfile(value):
                errors.append(
                    f"字段 '{field_name}' path=file，文件不存在: '{value}'"
                )

        elif path_type == "dir":
            if not os.path.isdir(value):
                errors.append(
                    f"字段 '{field_name}' path=dir，目录不存在: '{value}'"
                )

        elif path_type == "prefix":
            # Check if at least one index file exists for this prefix.
            # Infer tool name from field name (e.g. "hisat2_index_prefix" → "hisat2")
            # Prefer longer matches (bowtie2_for_rRNA over bowtie2)
            tool = None
            for name in sorted(cls.INDEX_MAP, key=len, reverse=True):
                if name in field_name:
                    tool = name
                    break
            if tool:
                exts = cls.INDEX_MAP[tool]
                has_any = any(os.path.exists(value + ext) for ext in exts)
                if not has_any:
                    errors.append(
                        f"字段 '{field_name}' path=prefix，索引文件不存在: "
                        f"'{value}' (检查了 {tool} 扩展名: {', '.join(exts[:3])}...)"
                    )
            else:
                # Unknown tool — just check if prefix path's parent dir exists
                parent = os.path.dirname(value)
                if parent and not os.path.isdir(parent):
                    errors.append(
                        f"字段 '{field_name}' path=prefix，父目录不存在: '{parent}'"
                    )

        return errors

    # ------------------------------------------------------------------
    # Enum constraint
    # ------------------------------------------------------------------

    @staticmethod
    def _validate_enum(
        field_name: str, value: Any, allowed: List[Any]
    ) -> List[str]:
        """Check that *value* is one of *allowed*.

        Args:
            field_name: Dotted key for error messages.
            value: The value to check.
            allowed: List of permitted values.

        Returns:
            List of error messages (empty if valid).
        """
        if value in allowed:
            return []
        display = ", ".join(repr(v) for v in allowed[:10])
        suffix = "..." if len(allowed) > 10 else ""
        return [
            f"字段 '{field_name}' 不在允许值范围内: "
            f"收到 '{value}'，有效值: {display}{suffix}"
        ]

    # ------------------------------------------------------------------
    # Range constraint
    # ------------------------------------------------------------------

    # Patterns for inferring numeric ranges from field names.
    # Each entry: (substring, min, max, description)
    # None means unbounded.
    _RANGE_PATTERNS: List[tuple] = [
        # probability / ratio fields → [0, 1]
        ("resolution",  0.0, 1.0, "probability"),
        ("fdr",         0.0, 1.0, "FDR"),
        ("p_cut",       0.0, 1.0, "p-value cutoff"),
        ("pvalue",      0.0, 1.0, "p-value"),
        ("p_adj",       0.0, 1.0, "adjusted p-value"),
        ("doublet_rate", 0.0, 1.0, "doublet rate"),
        ("max_pct_mt",  0.0, 100.0, "percentage"),
        ("identity",    0.0, 1.0, "identity"),
        # positive integers
        ("min_genes",   1, None, "positive integer"),
        ("max_genes",   1, None, "positive integer"),
        ("n_pcs",       1, None, "positive integer"),
        ("n_neighbors", 1, None, "positive integer"),
        ("n_top_genes", 1, None, "positive integer"),
        ("threads",     1, None, "positive integer"),
        ("top",         1, None, "positive integer"),
        ("threshold",   1, None, "positive integer"),
        ("quality",     1, None, "positive integer"),
    ]

    @classmethod
    def _infer_range(
        cls, field_name: str
    ) -> Optional[Tuple[Optional[float], Optional[float]]]:
        """Infer (min, max) bounds from a field name.

        Returns None if no pattern matches.
        """
        lower = field_name.lower()
        for pattern, lo, hi, _ in cls._RANGE_PATTERNS:
            if pattern in lower:
                return (lo, hi)
        return None

    @staticmethod
    def _validate_range(
        field_name: str, value: Any,
        minimum: Optional[float], maximum: Optional[float],
    ) -> List[str]:
        """Check that *value* falls within [minimum, maximum].

        Unbounded sides use None.  Non-numeric values are skipped.

        Args:
            field_name: Dotted key for error messages.
            value: The value to check.
            minimum: Lower bound (inclusive), or None.
            maximum: Upper bound (inclusive), or None.

        Returns:
            List of error messages (empty if valid).
        """
        if not isinstance(value, (int, float)):
            return []

        errors: List[str] = []
        if minimum is not None and value < minimum:
            errors.append(
                f"字段 '{field_name}' 值 {value} 小于下限 {minimum}"
            )
        if maximum is not None and value > maximum:
            errors.append(
                f"字段 '{field_name}' 值 {value} 大于上限 {maximum}"
            )
        return errors

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
                    # Determine tool from field name (prefer longer match)
                    tool = None
                    for name in sorted(self.INDEX_MAP, key=len, reverse=True):
                        if name in field_name:
                            tool = name
                            break
                    if tool:
                        exts = self.INDEX_MAP[tool]
                    else:
                        exts = [".1.ht2"]
                        tool = field_name.split("_index_prefix")[0] if "_index_prefix" in field_name else field_name
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

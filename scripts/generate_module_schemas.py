#!/usr/bin/env python3
"""Generate <stem>.schema.json next to every module config JSON.

Schema nodes mirror the config structure; leaf types are inferred from the
JSON values (bool/int/float/str/list/dict). String keys ending in _dir (or
named env_dir) get path="dir"; keys containing file-ish hints (fasta, gtf,
bed, tsv, ...) get path="file". outdir/logdir/indir are marked required.

Re-run after editing any module config JSON to refresh its schema.

Usage:
    python scripts/generate_module_schemas.py
"""
from __future__ import annotations

import json
import os
import sys

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "modules")

FILE_HINTS = ("fasta", "gtf", "bed", "tsv", "txt", "csv", "fa", "fq", "fastq", "bam", "vcf")
DIR_HINTS = ("_dir",)


def leaf_def(key: str, value):
    # ROOT_DIR is filled at runtime (Omics root dir) — always a directory string
    if key == "ROOT_DIR":
        return {"type": "str", "path": "dir", "nullable": True, "required": False}
    if value is None:
        return {"type": "null", "nullable": True, "required": False}
    if isinstance(value, bool):
        return {"type": "bool", "nullable": True, "required": False}
    if isinstance(value, int):
        return {"type": "int", "nullable": True, "required": False}
    if isinstance(value, float):
        return {"type": "float", "nullable": True, "required": False}
    if isinstance(value, str):
        kl = key.lower()
        if any(h in kl for h in DIR_HINTS) or kl == "env_dir":
            return {"type": "str", "path": "dir", "nullable": True, "required": False}
        if any(h in kl for h in FILE_HINTS):
            return {"type": "str", "path": "file", "nullable": True, "required": False}
        return {"type": "str", "nullable": True, "required": False}
    if isinstance(value, list):
        return {"type": "list", "nullable": True, "required": False}
    if isinstance(value, dict):
        return {"type": "dict", "nullable": True, "required": False,
                "properties": {k: leaf_def(k, v) for k, v in value.items()}}
    return {"type": "str", "nullable": True, "required": False}


def main() -> int:
    base = os.path.normpath(BASE)
    if not os.path.isdir(base):
        print(f"modules dir not found: {base}")
        return 1
    count = 0
    for mod in sorted(os.listdir(base)):
        mdir = os.path.join(base, mod)
        if not os.path.isdir(mdir) or mod.startswith("."):
            continue
        for root, dirs, files in os.walk(mdir):
            dirs[:] = [d for d in dirs if not d.startswith(".")]
            for f in sorted(files):
                if not f.endswith(".json") or "schema" in f:
                    continue
                jpath = os.path.join(root, f)
                try:
                    cfg = json.load(open(jpath, encoding="utf-8"))
                except (json.JSONDecodeError, UnicodeDecodeError):
                    print(f"skip unparseable: {jpath}")
                    continue
                if not isinstance(cfg, dict):
                    continue
                stem = f[:-5]
                schema = {k: leaf_def(k, v) for k, v in cfg.items()}
                for k in ("outdir", "logdir"):
                    if isinstance(schema.get(k), dict):
                        schema[k]["required"] = True
                spath = os.path.join(root, f"{stem}.schema.json")
                json.dump(schema, open(spath, "w", encoding="utf-8"),
                          indent=2, ensure_ascii=False)
                count += 1
    print(f"generated {count} module schemas")
    return 0


if __name__ == "__main__":
    sys.exit(main())

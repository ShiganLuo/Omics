#!/usr/bin/env python3
"""Standalone invocation of a single Omics module.

A module (modules/<name>/) is an independent unit with its own config JSON
(modules/<name>/**/<name>.json), own Snakefile, and its own inputs — it has
nothing to do with a workflow's meta.tsv. This script runs one module
directly:

    python run_module.py --module cutadapt \
        --config modules/cutadapt/cutadapt.json \
        -o /path/to/run_dir -t 8 --sdm apptainer

What it borrows from run.py's design (per request):
    - server-aware path switching: module JSON paths referencing a foreign
      server's Database dir (SERVER_DATABASE_DIR in node.py) are rewritten to
      the current server's;
    - container bind-path auto-collection (run.py::_collect_bind_paths) and
      --singularity-args merging (run.py::_merge_singularity_args);
    - mem_mb auto-detection and conda/apptainer backend handling identical to
      run.py's build_snakemake_cmd().

What it does NOT do: no meta.tsv parsing, no workflow dispatch, no rule-all
DAG of a whole workflow. The module's own Snakefile is the entrypoint and its
config JSON is the sole configuration (edit it, or override keys with dotted
--key.value args).

Target selection: most module Snakefiles have no `rule all` (their first rule
carries wildcards, which snakemake rejects as a default target). Without
--until/--target, this script expands each rule's output paths (evaluated with
outdir bound, wildcards kept as placeholders) and substitutes --samples values
for {sample_id}, passing the concrete files as targets. Use --rules to limit
which rules contribute targets, --until for a rule-level cutoff, or --target
for a single file.

Input convention: module JSONs describe inputs relative to `indir`
(e.g. cutadapt expects indir/<sample>/<sample>_1.fq.gz). Point `indir` at
your data (dotted override or edit the JSON) and `outdir` at the run dir.
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from typing import Any, Dict, List, Optional

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT_DIR)

from node import SERVER_DATABASE_DIR, _detect_current_server_db, _replace_server_paths  # noqa: E402
from run import (  # noqa: E402
    _collect_bind_paths,
    _detect_singularity,
    _detect_total_mem_mb,
    _merge_singularity_args,
)
from src.common.util.CmdUtil import _run_cmd  # noqa: E402
from src.common.util.SchemaValidatorUtil import SchemaValidator  # noqa: E402
from src.common.util.LogUtil import setup_logger  # noqa: E402
from src.common.util.SchemaValidatorUtil import smart_cast  # noqa: E402

logger = setup_logger("run_module", level=logging.DEBUG)




def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(__doc__ or "").splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-m","--module", default=None,
                        help="module name; smk/json auto-discovered under "
                             "modules/<name>/. Optional when --snakefile is "
                             "given (module inferred from the path).")
    parser.add_argument("--config", default=None,
                        help="module config JSON (default: auto-discovered "
                             "next to the snakefile)")
    parser.add_argument("--snakefile", default=None,
                        help="module Snakefile override (e.g. "
                             "modules/star/star.smk for the single-genome "
                             "variant instead of the default polygenomes/ one). "
                             "Implies --module (inferred from the path under "
                             "modules/); its sibling <stem>.json is used as "
                             "the config unless --config is given.")
    parser.add_argument("-o", "--outdir", default=None,
                        help="run/output directory (overrides config outdir; "
                             "default: ./<module>_run)")
    parser.add_argument("-i", "--indir", default=None,
                        help="input directory (overrides config indir)")
    parser.add_argument("-t", "--threads", type=int, default=8)
    parser.add_argument("--until", default=None,
                        help="stop after this rule (snakemake --until; its "
                             "upstream is built automatically)")
    parser.add_argument("--target", default=None,
                        help="build this output file/rule instead of `all`")
    parser.add_argument("--samples", nargs="+", default=None,
                        help="sample_id wildcard values; expands the module's "
                             "per-sample rules into concrete targets so a "
                             "module without `rule all` is runnable directly")
    parser.add_argument("--rules", nargs="+", default=None,
                        help="restrict to these rule(s) (module-local names); "
                             "default: every wildcarded per-sample rule")
    parser.add_argument("--genome", default=None,
                        help="restrict to this genome wildcard value")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--list-rules", action="store_true",
                        help="list the module's rules and exit")
    parser.add_argument("--sdm", action="store_const", const="apptainer", default=None,
                        help="apptainer SIF backend (recommended)")
    parser.add_argument("--conda-prefix", default=None)
    parser.add_argument("--conda-frontend", choices=["conda", "mamba"], default="mamba")
    parser.add_argument("--singularity-args", default=None)
    parser.add_argument("--rerun-trigger", nargs="+",
                        default=["code", "input", "mtime", "params", "software-env"])
    parser.add_argument("--mem-mb", type=int, default=None, dest="mem_mb")
    parser.add_argument("--no-schema-validate", dest="schema_validate",
                        action="store_false", default=True,
                        help="skip module config schema validation")
    parser.add_argument("--snakemake-args", nargs=argparse.REMAINDER, default=[])
    parser.add_argument("--log", default=None)
    args, unknown = parser.parse_known_args(argv)
    # dotted overrides: --Params.cutadapt.quality 30 / --key=value
    extra: Dict[str, Any] = {}
    i = 0
    while i < len(unknown):
        arg = unknown[i]
        if arg.startswith("--"):
            key = arg[2:]
            if "=" in key:
                k, v = key.split("=", 1)
                extra[k] = v
            else:
                # dotted/unknown flags take exactly one value in --key value
                # form; multi-value input uses --key=v1,v2 (comma-split), so a
                # following bare token can never be swallowed as a value.
                if i + 1 < len(unknown) and not unknown[i + 1].startswith("--"):
                    extra[key] = unknown[i + 1]
                    i += 1
                else:
                    extra[key] = True
        i += 1
    args.extra_args = extra
    return args


def discover_module_files(module: str,
                          snakefile_override: Optional[str] = None) -> tuple[str, str]:
    """Locate (snakefile, config_json) for a module.

    Default snakefile: modules/<m>/**/<m>.smk, preferring polygenomes/ (the
    variant every multi-genome subworkflow imports); the module-root copy is
    the single-genome legacy variant. With --snakefile, that exact file is
    used instead.

    The config JSON travels WITH the snakefile: sibling <stem>.json in the
    same directory (polygenomes/star.smk -> polygenomes/star.json,
    star.smk -> star.json). Only --config overrides this binding.
    """
    module_dir = os.path.join(ROOT_DIR, "modules", module)
    if not os.path.isdir(module_dir):
        raise SystemExit(f"Module not found: {module_dir}")

    if snakefile_override:
        snakefile = os.path.normpath(os.path.join(ROOT_DIR, snakefile_override)) \
            if not os.path.isabs(snakefile_override) else snakefile_override
        if not os.path.isfile(snakefile):
            raise SystemExit(f"--snakefile not found: {snakefile}")
        stem = os.path.splitext(os.path.basename(snakefile))[0]
        sibling_json = os.path.join(os.path.dirname(snakefile), f"{stem}.json")
        config = sibling_json if os.path.isfile(sibling_json) else None
        if config is None:
            raise SystemExit(
                f"No config JSON next to {snakefile} "
                f"(expected {sibling_json}); pass --config explicitly.")
        return snakefile, config

    def _includes_resolve(smk: str) -> bool:
        """True when every include: in the smk resolves to an existing file."""
        base = os.path.dirname(smk)
        for inc in re.findall(r'^include:\s*["\']([^"\']+)', open(smk).read(), re.M):
            if not os.path.isfile(os.path.normpath(os.path.join(base, inc))):
                return False
        return True

    def _all_smks() -> List[str]:
        hits: List[str] = []
        for root, _dirs, files in os.walk(module_dir):
            for f in files:
                if f.endswith(".smk"):
                    hits.append(os.path.join(root, f))
        return sorted(hits, key=len)

    def _pick(candidates: List[str]) -> Optional[str]:
        # polygenomes/ first (the variant subworkflows import), then shortest
        for h in candidates:
            if "polygenomes" in h and _includes_resolve(h):
                return h
        for h in candidates:
            if _includes_resolve(h):
                return h
        return None

    exact = [h for h in _all_smks() if os.path.basename(h) == f"{module}.smk"]
    chosen = _pick(exact)
    if chosen is None:
        # no exact-name smk resolves (e.g. fibertools: only polygenomes/fibertools.smk
        # exists but its include path is broken) — fall back to the module's other
        # entry-point smks (fibertools_prepare.smk / fibertools_analysis.smk).
        def _same_module(basename: str) -> bool:
            stem = os.path.splitext(basename)[0]
            norm = lambda s: s.replace("-", "_")
            return norm(stem).startswith(norm(module) + "_") or norm(stem) == norm(module)

        others = [h for h in _all_smks() if _same_module(os.path.basename(h))]
        chosen = _pick(others)
        if chosen is None:
            raise SystemExit(
                f"No runnable {module}*.smk under {module_dir}: "
                f"all candidates have unresolvable include paths.")
    snakefile = chosen
    stem = os.path.splitext(os.path.basename(snakefile))[0]
    config = os.path.join(os.path.dirname(snakefile), f"{stem}.json")
    if not os.path.isfile(config):
        raise SystemExit(f"Config not found next to snakefile: {config}")
    return snakefile, config


def load_module_config(config_path: str, module: str) -> Dict[str, Any]:
    """Load the module JSON with run.py-parity fixes.

    - server database path migration (node.py helpers)
    - ROOT_DIR + env: guarantee `env.env_dir` so sif() resolves SIFs even when
      the module JSON lacks an env section
    """
    import json
    with open(config_path, "r", encoding="utf-8") as f:
        cfg: Dict[str, Any] = json.load(f)

    cfg["ROOT_DIR"] = ROOT_DIR

    # env section: modules commonly omit it; the workflow configs carry it.
    env = cfg.setdefault("env", {})
    if not env.get("env_dir"):
        for db_dir in SERVER_DATABASE_DIR.values():
            candidate = os.path.join(db_dir, "env")
            if os.path.isdir(candidate):
                env["env_dir"] = candidate
                break
    # sif() fallback uses workflow.basedir (= wrapper dir when running via the
    # generated wrapper), producing <env_dir>/.run_module/<stem>.sif — wrong.
    # Pre-resolve each conda yaml stem of this module to an explicit SIF path
    # so sif() branch 1 (explicit mapping) always hits.
    env_dir = env.get("env_dir")
    if env_dir:
        module_root = os.path.join(ROOT_DIR, "modules", module)
        for root, _dirs, files in os.walk(module_root):
            for fname in files:
                if fname.endswith(".yaml"):
                    stem = os.path.splitext(fname)[0]
                    env.setdefault(stem, os.path.join(env_dir, module, f"{stem}.sif"))

    # server path migration (same as _init_raw_json)
    current_db = _detect_current_server_db()
    if current_db is not None:
        foreign_dirs = sorted(
            (d for d in SERVER_DATABASE_DIR.values() if d != current_db),
            key=len, reverse=True,
        )
        if foreign_dirs:
            _replace_server_paths(cfg, foreign_dirs, current_db)
    return cfg


def apply_overrides(cfg: Dict[str, Any], extra: Dict[str, Any]) -> None:
    """Apply dotted --key.value overrides (smart_cast, like run.py).

    Comma-separated values become lists: --genome_paired_samples.GRCm39=S1,S2
    """
    for key, value in extra.items():
        if isinstance(value, str) and "," in value:
            value = value.split(",")
        if "." in key:
            d = cfg
            keys = key.split(".")
            for k in keys[:-1]:
                d = d.setdefault(k, {})
            leaf = keys[-1]
            # *_samples keys (and their per-genome children) hold lists —
            # downstream code does membership tests against them.
            samples_context = any("sample" in k for k in keys[:-1]) or \
                re.search(r"samples?(_\w+)?$", leaf)
            if isinstance(value, str) and samples_context:
                value = [value]
            d[leaf] = smart_cast(value)
        else:
            cfg[key] = smart_cast(value)


RULE_RE = re.compile(r"rule\s+(\w+)\s*:(.*?)(?=\nrule\s|\Z)", re.S)


class _Wildcards:
    """Duck-typed wildcards for eval(): attribute access returns '{name}'."""

    def __getattr__(self, name: str) -> str:
        return "{" + name + "}"


def _eval_output_exprs(body: str, outdir: str) -> List[str]:
    """Evaluate a rule's output: entries with outdir bound; wildcards stay as
    '{name}' placeholders via a duck-typed wildcards object."""
    exprs: List[str] = []
    # strip the 'output:' prefix; each line like `name = expr,` or bare `expr,`
    payload = re.sub(r"^\s*output:\s*", "", body)
    # join continuation lines into one expression list
    for stmt in _split_output_stmts(payload):
        stmt = stmt.strip().rstrip(",")
        if not stmt:
            continue
        # named output: name = expr -> take RHS
        if re.match(r"^\w+\s*=\s*", stmt):
            stmt = re.sub(r"^\w+\s*=\s*", "", stmt)
        snakemake_markers = {
            # output markers in module smks; the path inside is what matters
            "directory": lambda x: x,
            "protected": lambda x: x,
            "ancient": lambda x: x,
            "temp": lambda x: x,
            "touch": lambda x: x,
            "report": lambda x: x,
        }
        try:
            value = eval(stmt, {"outdir": outdir, "wildcards": _Wildcards(),
                                "os": os, "config": {}, **snakemake_markers})
            if isinstance(value, (list, tuple)):
                exprs.extend(str(v) for v in value)
            else:
                exprs.append(str(value))
        except Exception:  # noqa: BLE001 — complex expressions are skipped
            continue
    return exprs


def _split_output_stmts(payload: str) -> List[str]:
    """Split output entries on top-level commas (respecting quotes/brackets)."""
    stmts: List[str] = []
    depth = 0
    quote = None
    cur = []
    for ch in payload:
        if quote:
            cur.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "\"'":
            quote = ch
            cur.append(ch)
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            stmts.append("".join(cur))
            cur = []
            continue
        cur.append(ch)
    if cur:
        stmts.append("".join(cur))
    return stmts


WRAPPER_TEMPLATE = """# Auto-generated by run_module.py — include the module as-is and add the
# missing wildcard-free default target. The module smk is NOT modified.
include: "{snakefile}"

rule all:
    input:
{inputs}
"""


def default_targets(snakefile: str, config_path: str, outdir: str,
                    samples: List[str], rules: Optional[List[str]],
                    sample_genomes: Optional[Dict[str, List[str]]] = None,
                    args_genome: Optional[str] = None) -> List[str]:
    """Concrete output patterns for the module's terminal rules.

    Collects each rule's output expressions (evaluated with outdir bound,
    wildcards kept as '{name}' placeholders); per-sample patterns expand to one
    pattern per sample. These become the wrapper's `rule all` inputs — snakemake
    resolves them at DAG-build time with real wildcard values, so input
    functions see concrete wildcards (unlike positional-file targets).
    """
    text = open(snakefile).read()
    genome_values = _genome_values_from_config(config_path)
    patterns: List[str] = []
    for m in RULE_RE.finditer(text):
        rule, body = m.group(1), m.group(2)
        if rules and rule not in rules:
            continue
        out_m = re.search(r"\boutput:\s*(.*?)(?=\n\s{4}\w+:|\Z)", body, re.S)
        if not out_m:
            continue
        for path in _eval_output_exprs(out_m.group(1), outdir):
            wcs = re.findall(r"\{(\w+)(?:,[^}]*)?\}", path)
            if "sample_id" in wcs:
                if not samples:
                    logger.warning("rule %s output needs --samples; skipped", rule)
                    continue
                for s in samples:
                    genomes = sample_genomes.get(s, []) if sample_genomes else genome_values
                    genomes = [g for g in genomes if not args_genome or g == args_genome]
                    if "genome" in wcs:
                        for g in genomes:
                            patterns.append(path
                                            .replace("{sample_id}", s)
                                            .replace("{genome}", g))
                    else:
                        patterns.append(path.replace("{sample_id}", s))
            elif "genome" in wcs:
                for g in genome_values:
                    patterns.append(path.replace("{genome}", g))
            elif not wcs:
                patterns.append(path)
    return sorted(set(patterns))


def _genome_values_from_config(config_path: str) -> List[str]:
    """Genome wildcard values: union of genome_paired/single_samples keys,
    restricted to genomes with a real fasta reference when any is configured
    (template placeholders like 'path/to/genome.fa' don't count)."""
    import json
    with open(config_path) as f:
        cfg = json.load(f)
    values = set()
    for key in ("genome_paired_samples", "genome_single_samples"):
        values.update((cfg.get(key) or {}).keys())

    refs = (cfg.get("genome") or {}).get("references") or {}
    if refs:
        usable = {g for g, ref in refs.items()
                  if isinstance(ref, dict) and ref.get("fasta")
                  and os.path.isfile(ref["fasta"])}
        if usable:
            values &= usable
    return sorted(values)


def _sample_genomes_from_config(config_path: str) -> Dict[str, List[str]]:
    """{sample_id -> [genome, ...]} from genome_paired/single_samples."""
    import json
    with open(config_path) as f:
        cfg = json.load(f)
    mapping: Dict[str, List[str]] = {}
    for key in ("genome_paired_samples", "genome_single_samples"):
        for genome, samples in (cfg.get(key) or {}).items():
            for s in samples or []:
                mapping.setdefault(s, []).append(genome)
    return mapping


def write_wrapper(snakefile: str, patterns: List[str], outdir: str) -> str:
    """Write the wrapper Snakefile with rule all: input: <patterns>.

    Each pattern is emitted as a Python string whose {wildcards} are expanded
    by snakemake's expand() inside rule all — using str.format semantics via
    an f-string so runtime wildcards are already concrete per DAG job.
    """
    wrapper_dir = os.path.join(outdir, ".run_module")
    os.makedirs(wrapper_dir, exist_ok=True)
    wrapper = os.path.join(wrapper_dir, "all.smk")
    lines = []
    for p in patterns:
        escaped = p.replace('"', '\\"')
        # plain string inside expand(): snakemake enumerates wildcard values
        # from the included module's rules at DAG-build time.
        lines.append(f'        expand("{escaped}"),')
    body = WRAPPER_TEMPLATE.format(snakefile=snakefile, inputs="\n".join(lines))
    with open(wrapper, "w", encoding="utf-8") as f:
        f.write(body)
    return wrapper


def build_cmd(wrapper: str, config_path: str, args: argparse.Namespace,
              mem_mb: Optional[int]) -> List[str]:
    """Assemble the snakemake command (mirrors run.py::build_snakemake_cmd).

    *wrapper* is the generated Snakefile that includes the module and defines
    `rule all`; the module smk itself is never modified.
    """
    use_singularity = args.sdm is not None or _detect_singularity(args.snakemake_args)
    cmd = [
        "snakemake", "-s", wrapper,
        "--configfile", config_path,
        "--cores", str(args.threads),
        "--rerun-triggers", *args.rerun_trigger,
    ]
    if mem_mb:
        cmd += ["--resources", f"mem_mb={mem_mb}"]
    if not use_singularity:
        if not args.conda_prefix:
            raise SystemExit("--conda-prefix is required when not using --sdm apptainer")
        cmd += ["--conda-prefix", args.conda_prefix, "--use-conda",
                "--conda-frontend", args.conda_frontend]
    else:
        cmd += ["--sdm", "apptainer"]
        singularity_args = _merge_singularity_args(
            _collect_bind_paths(config_path), args.singularity_args,
        )
        if singularity_args:
            cmd += ["--singularity-args", singularity_args]
    if args.dry_run:
        cmd.append("--dry-run")
    if args.list_rules:
        cmd.append("--list-rules")
    if args.until:
        cmd += ["--until", args.until]
    if args.target:
        cmd.append(args.target)
    cmd += args.snakemake_args
    return cmd


def infer_module_from_snakefile(snakefile: str) -> str:
    """Module name = the path component right after modules/ in the snakefile
    path. modules/star/polygenomes/star.smk -> star;
    modules/gatk/polygenomes/gatk_prepare.smk -> gatk.
    """
    parts = os.path.normpath(snakefile).split(os.sep)
    if "modules" in parts:
        i = parts.index("modules")
        if i + 1 < len(parts):
            return parts[i + 1]
    raise SystemExit(f"Cannot infer module name from snakefile: {snakefile}")


def main(argv: List[str] | None = None) -> None:
    args = parse_args(argv)
    if not args.module and not args.snakefile:
        raise SystemExit("--module or --snakefile is required")
    if not args.module:
        sf_abs = args.snakefile if os.path.isabs(args.snakefile) \
            else os.path.normpath(os.path.join(ROOT_DIR, args.snakefile))
        args.module = infer_module_from_snakefile(sf_abs)
        logger.info("module inferred from --snakefile: %s", args.module)
    module = args.module

    snakefile, config_path = discover_module_files(module, args.snakefile)
    logger.info("module=%s snakefile=%s config=%s", module, snakefile, config_path)

    cfg = load_module_config(config_path, module)
    apply_overrides(cfg, args.extra_args)

    # schema validation: <config_stem>.schema.json next to the module JSON
    schema_path = os.path.splitext(config_path)[0] + ".schema.json"
    schema_errors: List[str] = []
    if args.schema_validate and os.path.isfile(schema_path):
        import json as _json
        schema = _json.load(open(schema_path))
        schema_validator = SchemaValidator()
        schema_validator._schema = schema
        schema_errors = schema_validator.validate(cfg)
        if schema_errors:
            for e in schema_errors:
                logger.error(f"[schema] {e}")
            raise SystemExit(
                f"Module config failed schema validation ({schema_path}): "
                f"{len(schema_errors)} error(s). Fix the config or pass "
                "--no-schema-validate to skip.")
    else:
        logger.info("schema validation skipped for %s (no schema file or flag)", schema_path)

    # outdir/indir: CLI wins over config; outdir default ./<module>_run
    outdir = os.path.abspath(args.outdir or cfg.get("outdir") or f"{module}_run")
    if args.outdir or not cfg.get("outdir"):
        cfg["outdir"] = outdir
    if args.indir:
        cfg["indir"] = os.path.abspath(args.indir)
    # Log directory defaults to <outdir>/log. Module JSONs may carry relative
    # logdir values ("logs"); anchor them to outdir (snakemake resolves
    # relative paths against its cwd, which is NOT the outdir).
    logdir = cfg.get("logdir") or os.path.join(outdir, "log")
    if not os.path.isabs(logdir):
        logdir = os.path.join(outdir, logdir)
    cfg["logdir"] = logdir
    for key in list(cfg):
        if "logdir" in key and key != "logdir":
            v = cfg.get(key)
            if isinstance(v, str) and v and not os.path.isabs(v):
                cfg[key] = os.path.join(outdir, v)
    os.makedirs(cfg["outdir"], exist_ok=True)
    os.makedirs(cfg["logdir"], exist_ok=True)

    # write the effective config next to the module JSON's own copy? No — write
    # a run copy into outdir so the original module JSON stays pristine.
    import json
    effective_config = os.path.join(outdir, f"{module}.config.json")
    with open(effective_config, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False)
    logger.info("effective config: %s", effective_config)

    mem_mb = args.mem_mb or _detect_total_mem_mb()

    # Default target: a generated wrapper that includes the module smk and adds
    # `rule all` over the module's per-sample outputs (expand() fixes wildcards
    # at DAG-build time, so input functions see concrete values). The module
    # smk is never modified.
    wrapper = None
    if not args.target and not args.until and not args.list_rules:
        patterns = default_targets(snakefile, effective_config, cfg["outdir"],
                                   args.samples or [], args.rules,
                                   _sample_genomes_from_config(effective_config),
                                   args.genome)
        if not patterns:
            raise SystemExit(
                "No output patterns generated. Pass --samples (for per-sample "
                "rules) and/or --rules to select outputs, or use "
                "--until/--target.")
        wrapper = write_wrapper(snakefile, patterns, cfg["outdir"])
        logger.info("wrapper: %s (patterns: %d)", wrapper, len(patterns))

    cmd = build_cmd(wrapper or snakefile, effective_config, args, mem_mb)
    logger.info("%s", cmd)
    _run_cmd(cmd, cwd=ROOT_DIR)


if __name__ == "__main__":
    main()

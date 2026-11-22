import copy
import glob
from itertools import combinations
import os
import re
import time
from typing import Any, Dict, List, Optional, Sequence
import logging
import json
from src.common.util.type import DesignPair, CompareGroupPair, SampleInfo, CellrangerInput
from src.common.util.LogUtil import setup_logger
from src.common.util.SmkUtil import resolve_genome
from src.common.util.EnvUtil import is_path_like
logger = setup_logger(__name__, level=logging.DEBUG)
SERVER_DATABASE_DIR = {
    "zhang_c2": "/home/luosg/Database",
    "math_zhou": "/data/pub/zhousha/Database",
    "lv": "/disk5/luosg/Database",
    "lv_zheda": "/data/16T/Database"
}

# Genome reference key -> Dfam species alias for <alias>_te_lib.fa lookups
# (libraries under <server_db>/Reference/famdb/, built with
# modules/sv/bin/utils/famdb.py from dfam39_full.7.h5).
GENOME_TE_SPECIES = {
    "GRCh38": "human",
    "GRCm39": "mouse",
    "T2T_mhaESC_v1.5": "mouse",
    "Mmul_10": "rhesus",
    "Mmul_T2T": "rhesus",
}

# Genome reference key -> Ensembl VEP species/assembly names for the
# burden module's TMB annotation (modules/burden/bin/vep_ann.py). Genomes
# without an entry must set Params.burden.vep_species/vep_assembly
# explicitly — never guess a species.
GENOME_VEP_SPECIES = {
    "GRCh38": ("homo_sapiens", "GRCh38"),
    "GRCm39": ("mus_musculus", "GRCm39"),
    "Mmul_10": ("macaca_mulatta", "Mmul_10"),
}

# Genome reference key -> sex-chromosome contig names (X, Y) for
# sex-inference and any sex-aware analysis (modules/sexcheck first
# consumer). Naming is assembly-dependent: UCSC-style chrX/chrY
# (GENCODE GRCm39/GRCh38), Ensembl bare X/Y (Mmul_10), NCBI accessions
# (T2T assemblies: CM166267.1/CM166268.1, NC_133426.1/NC_133427.1 —
# verified against the local FASTA headers 2026-11). Filled into
# genome.references.<g>.sex_chroms when unset; genomes without an entry
# fall back to the module's name sniffing (chrX/X, chrY/Y) — add an
# entry here when a new assembly's sex contigs use neither convention.
GENOME_SEX_CHROMS = {
    "GRCh38": ("chrX", "chrY"),
    "GRCm39": ("chrX", "chrY"),
    "GRCm38": ("chrX", "chrY"),
    "Mmul_10": ("X", "Y"),
    "T2T_mhaESC_v1.5": ("CM166267.1", "CM166268.1"),
    "Mmul_T2T": ("NC_133426.1", "NC_133427.1"),
}

# Genome reference key -> RepeatMasker/Dfam species term for the sv.mei
# Dfam classification path (Params.sv.mei.dfam_h5; the famdb.py export
# term and RM -species search term, e.g. 'Mus musculus'). Filled into
# Params.sv.mei.rm_species when unset — genomes without an entry must
# set it explicitly, never guess a species (same policy as VEP).
GENOME_RM_SPECIES = {
    "GRCh38": "Homo sapiens",
    "GRCm39": "Mus musculus",
    "T2T_mhaESC_v1.5": "Mus musculus",
    "Mmul_10": "Macaca mulatta",
    "Mmul_T2T": "Macaca mulatta",
}

# Genome reference key -> telogator2 non-human reference alias. telogator2
# resolves `telogator-ref-<alias>.fa.gz` in its install; the alias MUST match
# that filename. Filled into Params.telomere.telogator2.species when unset —
# genomes without an entry must set it explicitly, never guess a species.
# (verify aliases against the telogator2 install's non-human/ refs.)
GENOME_TELOGATOR2_SPECIES = {
    "GRCh38": "human",
    "GRCm39": "mouse",
    "Mmul_10": "macaca_mulatta",
    "Mmul_T2T": "macaca_mulatta",
}


# MEI TE-bucket vocabularies live in the config templates
# (Params.sv.mei.te_buckets keyed by genome reference key, like
# genome.references) — node.py only selects the active genome's list.
# Bucket schema and matching semantics are documented in
# modules/sv/bin/mei_analysis.py above parse_te_buckets.

# Genome reference key -> TE-bucket preset name (for the TE-lib sentinel
# check below only).
GENOME_TE_BUCKET_PRESET = {
    "GRCh38": "primate",
    "GRCm39": "mouse",
    "T2T_mhaESC_v1.5": "mouse",
    "Mmul_10": "primate",
    "Mmul_T2T": "primate",
}

# Preset -> sentinel family-name markers that a sane TE consensus library
# MUST contain (see utils/famdb.py SENTINEL_SUBFAMILIES; mirrored here so
# the check runs at config time without importing famdb).
TE_LIB_SENTINELS = {
    "mouse": (["l1m", "l1mus"], ["b1"], ["b2"], ["iap"], ["etn", "mervl"]),
    "primate": (["l1p", "l1hs", "l1_"], ["alu"], ["sva"], ["erv", "ltr"]),
}

def _validate_te_lib(datajson: Dict[str, Any], workflow: str) -> None:
    """Warn when the genome's TE library lacks species sentinel families.

    Deterministic for the wrong-branch failure mode (a library with only
    ancestor-level families); cannot detect sequence-level corruption.
    """
    genome_key = str(datajson.get("genome", {}).get("default"))
    preset = GENOME_TE_BUCKET_PRESET.get(genome_key, "primate")
    te_lib = ((datajson.get("genome", {}).get("references", {})
               .get(genome_key) or {}).get("te_lib"))
    groups = TE_LIB_SENTINELS.get(preset)
    if not groups or not te_lib or not os.path.exists(te_lib):
        return
    names = set()
    with open(te_lib) as fh:
        for line in fh:
            if line.startswith(">"):
                names.add(line[1:].split()[0].lower())
                if len(names) > 5000:
                    break
    missing = ["/".join(g) for g in groups
               if not any(any(m in n for n in names) for m in g)]
    if missing:
        logger.warning(f"{workflow}: TE library {te_lib} lacks sentinel "
            f"families for preset '{preset}': {', '.join(missing)}. It was "
            f"likely curated on the wrong lineage branch — every insertion "
            f"may classify as noise. Rebuild with: python "
            f"modules/sv/bin/utils/famdb.py -i <dfam.h5> -s "
            f"{'mouse' if preset == 'mouse' else 'human'} -o {te_lib}")


def _fill_mei_rm_species(datajson: Dict[str, Any], workflow: str) -> None:
    """Auto-fill Params.sv.mei.rm_species from genome.default.

    The sv.mei Dfam path (Params.sv.mei.dfam_h5) classifies insertions
    with RepeatMasker against a Dfam FamDB database and needs the
    species term for the curated famdb export ('Mus musculus' ...).
    An explicit config value always wins; for genomes without a
    GENOME_RM_SPECIES entry the value must be set explicitly — never
    guess a species (same policy as the burden VEP settings).
    """
    mei = (datajson.get("Params", {}).get("sv", {}).get("mei") or {})
    if mei.get("rm_species") or not mei.get("dfam_h5"):
        return
    genome_key = str(datajson.get("genome", {}).get("default"))
    species = GENOME_RM_SPECIES.get(genome_key)
    if not species:
        raise ValueError(
            f"{workflow}: genome '{genome_key}' has no known RepeatMasker "
            f"species for the sv.mei Dfam classification; set "
            f"'Params.sv.mei.rm_species' in the config JSON, or drop "
            f"Params.sv.mei.dfam_h5 to keep the te_lib path."
        )
    datajson.setdefault("Params", {}).setdefault("sv", {}).setdefault(
        "mei", {})["rm_species"] = species
    logger.info(
        f"{workflow}: Params.sv.mei.rm_species := '{species}' "
        f"(genome '{genome_key}')")


def _fill_telomere_species(datajson: Dict[str, Any], workflow: str) -> None:
    """Auto-fill Params.telomere.telogator2.species from genome.default.

    telogator2 picks its subtelomere reference by species
    (`telogator-ref-<alias>.fa.gz`). Explicit config wins; genomes without a
    GENOME_TELOGATOR2_SPECIES entry must set it explicitly — never guess.
    """
    tel = (datajson.get("Params", {}).get("telomere", {}) or {})
    if (tel.get("telogator2", {}) or {}).get("species"):
        return
    genome_key = str(datajson.get("genome", {}).get("default"))
    species = GENOME_TELOGATOR2_SPECIES.get(genome_key)
    if not species:
        raise ValueError(
            f"{workflow}: genome '{genome_key}' has no known telogator2 species "
            f"alias; set 'Params.telomere.telogator2.species' in the config JSON "
            f"(must match telogator2's telogator-ref-<alias>.fa.gz).")
    datajson.setdefault("Params", {}).setdefault("telomere", {}).setdefault(
        "telogator2", {})["species"] = species
    logger.info(
        f"{workflow}: Params.telomere.telogator2.species := '{species}' "
        f"(genome '{genome_key}')")


# Genome reference key -> telogator2-format subtelomere tip reference
# (contigs '<build>_<chr><arm>', ~500 kb per arm tip — telogator2's -t
# contract, see modules/telomere/length/length.smk get_telogator2_ref).
# A whole-genome FASTA is NOT a valid -t; each genome needs its own tip
# reference. Filled into Params.telomere.telogator2.subtel_ref when unset —
# genomes without an entry must set it explicitly, never guess.
GENOME_TELOGATOR2_SUBTEL_REF = {
    "Mmul_T2T": "/data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/telogator-ref-macaca_mulatta.fa.gz",
}


def _fill_telomere_subtel_ref(datajson: Dict[str, Any], workflow: str) -> None:
    """Auto-fill Params.telomere.telogator2.subtel_ref from genome.default.

    telogator2 -t takes a per-arm subtelomere tip reference, not the genome
    FASTA (its anchor parsing hard-codes '<build>_<chr><arm>' contig names).
    Explicit config wins; genomes without a GENOME_TELOGATOR2_SUBTEL_REF
    entry fall through to telogator2's built-in refs (human default,
    mouse/maize aliases) or fail fast in the module — never guess.
    """
    tel = (datajson.get("Params", {}).get("telomere", {}) or {})
    if (tel.get("telogator2", {}) or {}).get("subtel_ref"):
        return
    genome_key = str(datajson.get("genome", {}).get("default"))
    subtel_ref = GENOME_TELOGATOR2_SUBTEL_REF.get(genome_key)
    if not subtel_ref:
        return
    datajson.setdefault("Params", {}).setdefault("telomere", {}).setdefault(
        "telogator2", {})["subtel_ref"] = subtel_ref
    logger.info(
        f"{workflow}: Params.telomere.telogator2.subtel_ref := '{subtel_ref}' "
        f"(genome '{genome_key}')")


# Genome reference key -> diploid chromosome count, the denominator of
# read_density's "Avg telomere per arm" (= total_tel_bp / n_chrom_arms;
# modules/telomere/length/README.md §3.3). Module default 40 is the mouse
# value (20 chromosomes, diploid count); rhesus is 44 per the same README.
# Filled into Params.telomere.read_density.n_chrom_arms when unset — genomes
# without an entry must set it explicitly, never guess (same policy as
# telogator2).
GENOME_N_CHROM_ARMS = {
    "GRCh38": 46,
    "GRCm39": 40,
    "T2T_mhaESC_v1.5": 40,
    "Mmul_10": 44,
    "Mmul_T2T": 44,
}


def _fill_telomere_n_chrom_arms(datajson: Dict[str, Any], workflow: str) -> None:
    """Auto-fill Params.telomere.read_density.n_chrom_arms from genome.default.

    Only the coverage-dependent "Avg telomere per arm" metric divides by it
    (module README §3.3; the robust Telomeric fraction is unaffected).
    Explicit config wins; genomes without a GENOME_N_CHROM_ARMS entry must set
    it explicitly — never guess a chromosome count.
    """
    tel = (datajson.get("Params", {}).get("telomere", {}) or {})
    if (tel.get("read_density", {}) or {}).get("n_chrom_arms") is not None:
        return
    genome_key = str(datajson.get("genome", {}).get("default"))
    n_arms = GENOME_N_CHROM_ARMS.get(genome_key)
    if not n_arms:
        raise ValueError(
            f"{workflow}: genome '{genome_key}' has no known chromosome-arm "
            f"count; set 'Params.telomere.read_density.n_chrom_arms' in the "
            f"config JSON (diploid chromosome count, telomere length README "
            f"§3.3).")
    datajson.setdefault("Params", {}).setdefault("telomere", {}).setdefault(
        "read_density", {})["n_chrom_arms"] = n_arms
    logger.info(
        f"{workflow}: Params.telomere.read_density.n_chrom_arms := '{n_arms}' "
        f"(genome '{genome_key}')")


def _fill_repeatmasker_species(datajson: Dict[str, Any], workflow: str) -> None:
    """Auto-fill Params.RepeatMasker species (RepeatMasker -species) from genome.

    Reuses GENOME_RM_SPECIES (same term as the sv.mei Dfam path). Explicit
    config wins; genomes without an entry must set it explicitly — never guess.
    """
    rm = (datajson.get("Params", {}).get("RepeatMasker", {}) or {})
    if rm.get("species") or (rm.get("repeatmasker_init", {}) or {}).get("species"):
        return
    genome_key = str(datajson.get("genome", {}).get("default"))
    species = GENOME_RM_SPECIES.get(genome_key)
    if not species:
        raise ValueError(
            f"{workflow}: genome '{genome_key}' has no known RepeatMasker species; "
            f"set 'Params.RepeatMasker.species' in the config JSON.")
    datajson.setdefault("Params", {}).setdefault("RepeatMasker", {})["species"] = species
    logger.info(
        f"{workflow}: Params.RepeatMasker.species := '{species}' "
        f"(genome '{genome_key}')")



def _detect_current_server_db() -> str | None:
    """Detect current server's database dir by checking which paths in SERVER_DATABASE_DIR exist."""
    for db_dir in SERVER_DATABASE_DIR.values():
        if os.path.isdir(db_dir):
            return db_dir
    return None

# [NAME=]chrom:start-end region row (thousands separators tolerated)
_REGION_RE = re.compile(r"([^:]+):([\d,]+)-([\d,]+)")


def _enabled(datajson: Dict[str, Any], step: str) -> bool:
    """Read Params.<step>.enabled (default True).

    Each downstream step has its own gate under `Params.<step>` so users
    can switch off individual analyses (e.g. eccDNA) without affecting
    the rest of the pipeline.
    """
    block = datajson.get("Params", {}).get(step)
    if isinstance(block, dict) and "enabled" in block:
        return bool(block["enabled"])
    return True



def _figs(base: str, fmts) -> List[str]:
    """Per-format figure files for one extension-less figure base."""
    return [f"{base}.{fmt}" for fmt in fmts]


def _register_regionplot(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        samples: List[str],
        outdir: str,
        workflow: str,
    ) -> List[str]:
    """Register regionplot panel outputs under the Params.regionplot gate.

    Shared by the six bw-producing workflows (Fiberseq/PacVar/CLIP/MERIP/
    ncRNAseq/PeakCalling). Requires ``Params.regionplot.regions`` when
    enabled, auto-fills ``sample_group_dict`` from the meta group column
    (a hand-filled config value wins), and returns the outfile list to
    append (per panel: every requested format + tracks.ini +
    track_sources.tsv + regions.bed).

    Parameters
    ----------
    datajson : the workflow config being materialized (mutated: the
        ``Params.regionplot`` block gets ``sample_group_dict`` written back).
    samples_info_dict : meta sample info (source of the group column).
    samples : sample ids that own bigwigs in this workflow.
    outdir : workflow output root.
    workflow : workflow name (for error messages).

    Returns
    -------
    List[str]
        Outfile paths to register; empty when the gate is off.
    """
    if not _enabled(datajson, "regionplot"):
        return []
    # 参数精确到 rule:Params.regionplot.regionplot.<key>(消费 rule regionplot;
    # config 全量路径读取,无 module 级共享回退)。
    # regions: BED rows or a [NAME=]chrom:start-end list
    regions_file = (datajson.get("Params", {}).get("regionplot", {})
                    .get("regionplot", {}).get("regions"))
    if not regions_file:
        raise ValueError(
            f"{workflow}: Params.regionplot.enabled=true requires "
            "Params.regionplot.regionplot.regions (a BED of chrom/start/end[/name] "
            "rows or a list of [NAME=]chrom:start-end lines); set it in "
            "the workflow config JSON, then rerun.")
    # group map: hand-filled config wins, else derive from meta group column
    group_dict = dict(datajson.get("Params", {}).get("regionplot", {})
                      .get("regionplot", {}).get("sample_group_dict") or {})
    for sid in samples:
        if sid not in group_dict:
            group = getattr(samples_info_dict.get(sid), "group", None)
            if group:
                group_dict[sid] = group
    # 写回精确到 rule:Params.regionplot.regionplot.sample_group_dict
    datajson.setdefault("Params", {}).setdefault("regionplot", {}).setdefault(
        "regionplot", {})["sample_group_dict"] = group_dict
    # formats normalize (same expression as the module smk)
    formats = (datajson.get("Params", {}).get("regionplot", {})
               .get("regionplot", {}).get("formats") or ["png", "svg"])
    if isinstance(formats, str):
        formats = [formats]
    if not isinstance(formats, list) or not all(isinstance(f, str) for f in formats):
        raise ValueError(f"{workflow}: Params.regionplot.regionplot.formats must be a list of strings")
    # panel names mirror the module's regions parsing (BED col4 or
    # chrom_start_end; list NAME= prefix or the sanitized chr:start-end)
    panels: List[str] = []
    seen: Dict[str, int] = {}
    with open(regions_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            line = line.rstrip("\n")
            if "\t" in line:
                f = line.split("\t")
                name = (f[3].strip() if len(f) >= 4 and f[3].strip()
                        else f"{f[0]}_{f[1]}_{f[2]}")
            else:
                s = line.strip()
                if "=" in s:
                    name = s.split("=", 1)[0]
                else:
                    m = _REGION_RE.fullmatch(s)
                    if m is None:
                        raise ValueError(
                            f"{workflow}: malformed region row {s!r} in "
                            f"{regions_file} (expected [NAME=]chrom:start-end)")
                    name = (f"{m.group(1)}_"
                            f"{int(m.group(2).replace(',', ''))}_"
                            f"{int(m.group(3).replace(',', ''))}")
            for ch in ("/", " ", ":", "="):
                name = name.replace(ch, "_")
            if name in seen:
                seen[name] += 1
                name = f"{name}.{seen[name]}"
            else:
                seen[name] = 0
            panels.append(name)
    if not panels:
        raise ValueError(
            f"{workflow}: Params.regionplot.regionplot.regions has no rows: {regions_file}")
    out = []
    for panel in panels:
        for fmt in formats:
            out.append(f"{outdir}/regionplot/{panel}/{panel}.{fmt}")
        out.append(f"{outdir}/regionplot/{panel}/tracks.ini")
        out.append(f"{outdir}/regionplot/{panel}/track_sources.tsv")
        out.append(f"{outdir}/regionplot/{panel}/regions.bed")
    return out


def _register_fiberbrowser(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        samples: List[str],
        outdir: str,
        workflow: str,
    ) -> List[str]:
    """Register fiberbrowser figure outputs under the Params.fiberbrowser gate.

    Figure delivery is opt-in (default OFF; unlike _enabled's default-on —
    a presentation deliverable that is not requested must not materialize).
    When enabled: requires ``Params.fiberbrowser.loci`` (BED rows or
    [NAME=]chrom:start-end lines); parses it into the shell-safe
    ``loci_map`` ({panel_name: "chrom:start-end"}) written back to the
    config — the module consumes it for output naming and per-locus engine
    calls; auto-fills ``sample_group_dict`` from the meta group column (a
    hand-filled config value wins).  Returns the figure outfiles to append
    (locus x variant x format).

    Parameters
    ----------
    datajson : the workflow config being materialized (mutated: the
        ``Params.fiberbrowser`` block gets ``loci_map`` and
        ``sample_group_dict`` written back).
    samples_info_dict : meta sample info (source of the group column).
    samples : sample ids whose alignments feed the figures.
    outdir : workflow output root.
    workflow : workflow name (for error messages).

    Returns
    -------
    List[str]
        Outfile paths to register; empty when the gate is off.
    """
    if not (datajson.get("Params", {}).get("fiberbrowser", {}) or {}).get("enabled"):
        return []
    loci_file = (datajson.get("Params", {}).get("fiberbrowser", {}) or {}).get("loci")
    if not loci_file:
        raise ValueError(
            f"{workflow}: Params.fiberbrowser.enabled=true requires "
            "Params.fiberbrowser.loci (a BED of chrom/start/end[/name] rows "
            "or a list of [NAME=]chrom:start-end lines); set it in the "
            "workflow config JSON, then rerun.")
    # group map: hand-filled config wins, else derive from meta group column
    # (参数精确到 rule:Params.fiberbrowser.fiberbrowser.<key>;config 全量路径
    # 读取,无 module 级共享回退)
    group_dict = dict(datajson.get("Params", {}).get("fiberbrowser", {})
                      .get("fiberbrowser", {}).get("sample_group_dict") or {})
    for sid in samples:
        if sid not in group_dict:
            group = getattr(samples_info_dict.get(sid), "group", None)
            if group:
                group_dict[sid] = group
    formats = (datajson.get("Params", {}).get("fiberbrowser", {})
               .get("fiberbrowser", {}).get("formats") or ["png"])
    if isinstance(formats, str):
        formats = [formats]
    if not isinstance(formats, list) or not all(isinstance(f, str) for f in formats):
        raise ValueError(f"{workflow}: Params.fiberbrowser.fiberbrowser.formats must be a list of strings")
    variants = (datajson.get("Params", {}).get("fiberbrowser", {})
                .get("fiberbrowser", {}).get("variants") or ["plain", "zoom"])
    if isinstance(variants, str):
        variants = [variants]
    variant_suffix = {"plain": "", "zoom": "_zoom"}
    bad = [v for v in variants if v not in variant_suffix]
    if bad:
        raise ValueError(
            f"{workflow}: Params.fiberbrowser.fiberbrowser.variants entries must be "
            f"'plain' or 'zoom', got {bad}")
    # loci parse; names become output path components AND receipt arguments,
    # so sanitize to a shell-safe charset and dedupe
    loci_map: Dict[str, str] = {}
    seen: Dict[str, int] = {}
    with open(loci_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            line = line.rstrip("\n")
            if "\t" in line:
                f = line.split("\t")
                name = (f[3].strip() if len(f) >= 4 and f[3].strip()
                        else f"{f[0]}_{f[1]}_{f[2]}")
                region = f"{f[0]}:{f[1]}-{f[2]}"
                if _REGION_RE.fullmatch(region) is None:
                    raise ValueError(
                        f"{workflow}: malformed BED row {line!r} in "
                        f"{loci_file} (expected chrom<TAB>start<TAB>end"
                        "[<TAB>name])")
            else:
                s = line.strip()
                if "=" in s:
                    name, region = s.split("=", 1)
                else:
                    region = s
                    m = _REGION_RE.fullmatch(s)
                    if m is None:
                        raise ValueError(
                            f"{workflow}: malformed locus row {s!r} in "
                            f"{loci_file} (expected [NAME=]chrom:start-end)")
                    name = (f"{m.group(1)}_"
                            f"{int(m.group(2).replace(',', ''))}_"
                            f"{int(m.group(3).replace(',', ''))}")
            name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
            if name in seen:
                seen[name] += 1
                name = f"{name}.{seen[name]}"
            else:
                seen[name] = 0
            loci_map[name] = region.replace(",", "")
    if not loci_map:
        raise ValueError(
            f"{workflow}: Params.fiberbrowser.loci has no rows: {loci_file}")
    # 写回精确到 rule:Params.fiberbrowser.fiberbrowser.<key>(模块 rule fiberbrowser
    # 经 rule_param 读取,无 module 级共享回退)
    fiberbrowser = (datajson.setdefault("Params", {}).setdefault("fiberbrowser", {})
                    .setdefault("fiberbrowser", {}))
    fiberbrowser["loci_map"] = loci_map
    fiberbrowser["sample_group_dict"] = group_dict
    out = []
    for name in loci_map:
        for v in variants:
            for fmt in formats:
                out.append(
                    f"{outdir}/fiberbrowser/{name}/{name}{variant_suffix[v]}.{fmt}")
    return out


def _register_fibercluster(
        datajson: Dict[str, Any],
        outdir: str,
        workflow: str,
    ) -> List[str]:
    """Register fibercluster outputs under the Params.fibercluster gate.

    Analysis deliverable, opt-in like fiberbrowser (default OFF — an
    analysis that is not requested must not materialize). When enabled:
    requires ``Params.fibercluster.regions_bed`` (a BED of
    chrom/start/end/name rows; region names drive the per-region outputs),
    a non-empty ``settings`` label list, a
    ``channel_by_setting``/``k_by_setting`` entry for EVERY label (the
    per-run -C/-k; a missing entry fails HERE naming the key instead of
    mid-run), and ``fiber_report.canonical``. Returns the outfile paths:
    one done marker per setting (the module tracks per-region outputs by
    marker — region names live in the BED, not in wildcards) plus the
    report-package files.

    Parameters
    ----------
    datajson : the workflow config being materialized (read only).
    outdir : workflow output root.
    workflow : workflow name (for error messages).

    Returns
    -------
    List[str]
        Outfile paths to register; empty when the gate is off.
    """
    if not (datajson.get("Params", {}).get("fibercluster", {}) or {}).get("enabled"):
        return []
    p = datajson.get("Params", {}).get("fibercluster", {}) or {}
    regions_bed = p.get("regions_bed")
    if not regions_bed:
        raise ValueError(
            f"{workflow}: Params.fibercluster.enabled=true requires "
            "Params.fibercluster.regions_bed (a BED of chrom/start/end/name "
            "rows; region names drive the per-region outputs); set it in the "
            "workflow config JSON, then rerun.")
    settings = p.get("settings") or []
    if isinstance(settings, str):
        settings = [settings]
    if not settings:
        raise ValueError(
            f"{workflow}: Params.fibercluster.enabled=true requires a non-empty "
            "Params.fibercluster.settings label list (one engine run per label).")
    fc = p.get("fiber_cluster", {}) or {}
    for label in settings:
        for key in ("channel_by_setting", "k_by_setting"):
            if label not in (fc.get(key) or {}):
                raise ValueError(
                    f"{workflow}: Params.fibercluster.fiber_cluster.{key} has no "
                    f"entry for settings label '{label}' — add it "
                    + ("(msp|fire)." if key == "channel_by_setting" else "(int)."))
    canonical = (p.get("fiber_report", {}) or {}).get("canonical")
    if not canonical:
        raise ValueError(
            f"{workflow}: Params.fibercluster.enabled=true requires "
            "Params.fibercluster.fiber_report.canonical (the settings label "
            "whose figures/summary become the deliverable).")
    if canonical not in settings:
        raise ValueError(
            f"{workflow}: Params.fibercluster.fiber_report.canonical "
            f"'{canonical}' is not among the settings labels {settings}.")
    base = f"{outdir}/fibercluster"
    out = [f"{base}/{label}/fiber_cluster.done" for label in settings]
    out += [
        f"{base}/summary/loci_summary.tsv",
        f"{base}/summary/cluster_proportions.png",
        f"{base}/sensitivity/sensitivity_table.tsv",
        f"{base}/figures.done",
    ]
    return out


def _resolve_genome_ref(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        workflow: str,
    ) -> Dict[str, Any]:
    """Resolve and normalize ``genome.default``, return its reference block.

    ``genome.default`` (which may be null) is user input and may be any
    species alias; it is normalized through ``resolve_genome()`` exactly like
    ``sample_info.organism`` is during meta parsing (e.g. ``"mouse"`` /
    ``"mm39"`` -> ``"GRCm39"``) so it matches the ``genome.references`` keys,
    and unrecognised names fall back to the raw string so custom reference
    keys (e.g. ``"T2T_mhaESC_v1.5"``) keep working. When ``genome.default``
    is null, the normalized sample organism recorded in the meta is used as
    the fallback. A mismatch between the two is reported but the config value
    is kept.

    The normalized key is written back to ``datajson["genome"]["default"]``
    so the subworkflow resolves the same reference entry.

    Raises
    ------
    ValueError
        When neither source resolves to a ``genome.references`` entry.
    """
    genome_cfg = datajson.get("genome") or {}
    references = genome_cfg.get("references") or {}
    raw_default = genome_cfg.get("default")
    sample_genomes = sorted({
        info.organism for info in samples_info_dict.values()
        if info.organism and info.organism != "UNKNOWN"
    })

    normalized: Optional[str] = None
    if raw_default:
        try:
            normalized = resolve_genome(str(raw_default))
        except ValueError:
            # Custom reference key (e.g. T2T_mhaESC_v1.5): use verbatim.
            normalized = str(raw_default)
        if normalized != str(raw_default):
            logger.info(
                f"{workflow}: genome.default '{raw_default}' normalized to "
                f"'{normalized}' (sample_info.organism convention)"
            )
        if (sample_genomes and normalized != str(raw_default)
                and normalized not in sample_genomes):
            logger.warning(
                f"{workflow}: genome.default '{normalized}' does not match the "
                f"sample organism(s) recorded in the meta {sample_genomes}; "
                f"keeping the config value. Check the meta organism column or "
                f"genome.default before trusting variant calls."
            )
    elif sample_genomes:
        normalized = sample_genomes[0]
        if len(sample_genomes) > 1:
            logger.warning(
                f"{workflow}: samples carry several organisms {sample_genomes}; "
                f"using '{normalized}' for genome.default"
            )
        logger.info(
            f"{workflow}: genome.default is null; using the sample organism "
            f"'{normalized}' recorded in the meta"
        )
    else:
        raise ValueError(
            f"{workflow}: genome.default is null and the meta has no "
            f"'organism' column; set 'genome.default' in the config JSON or "
            f"record the organism in the meta, then rerun."
        )

    ref = references.get(normalized)
    if not ref:
        raise ValueError(
            f"{workflow}: 'genome.references' has no entry '{normalized}' "
            f"(from genome.default '{raw_default}'); "
            f"known entries: {sorted(references)}"
        )

    # Auto-fill derived reference resources: a null fai_index next to an existing
    # FASTA is filled with <fasta>.fai so the subworkflow's genome_index rule
    # (samtools faidx) is scheduled when the index is missing. An explicit
    # template/config value always wins; the fill only happens when the key
    # is absent or null. (Legacy key 'fai' is no longer read: configs and
    # subworkflows are unified on fai_index/dict_index.)
    fasta_path = ref.get("fasta")
    if fasta_path and not ref.get("fai_index"):
        ref["fai_index"] = f"{fasta_path}.fai"
        logger.info(
            f"{workflow}: genome.fai_index not set; auto-filled '{ref['fai_index']}' "
            f"(genome_index rule builds it if missing)"
        )

    # Auto-fill sex-chromosome contig names from the genome key (modules/
    # sexcheck + future sex-aware analyses); an explicit value wins.
    if not ref.get("sex_chroms"):
        mapped = GENOME_SEX_CHROMS.get(normalized)
        if mapped:
            ref["sex_chroms"] = list(mapped)
            logger.info(
                f"{workflow}: genome.references.{normalized}.sex_chroms "
                f"auto-filled: {list(mapped)}"
            )

    genome_cfg["default"] = normalized
    datajson["genome"] = genome_cfg
    return ref

def _replace_server_paths(obj: Any, foreign_dirs: List[str], current_db: str) -> None:
    """Recursively replace foreign server database paths with current server's.

    For every string value that is path-like and starts with a foreign server's
    database directory, replace the prefix with the current server's database
    directory.  Operates in-place on dicts and lists.

    Args:
        obj: The dict/list/str to walk.
        foreign_dirs: Sorted (longest-first) foreign database directory prefixes.
        current_db: Current server's database directory.
    """
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, str) and is_path_like(v):
                for foreign in foreign_dirs:
                    if v.startswith(foreign):
                        obj[k] = current_db + v[len(foreign):]
                        logger.info(f"[server_path] {v} -> {obj[k]}")
                        break
            else:
                _replace_server_paths(v, foreign_dirs, current_db)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            if isinstance(v, str) and is_path_like(v):
                for foreign in foreign_dirs:
                    if v.startswith(foreign):
                        obj[i] = current_db + v[len(foreign):]
                        logger.info(f"[server_path] {v} -> {obj[i]}")
                        break
            else:
                _replace_server_paths(v, foreign_dirs, current_db)

def _init_raw_json(
    datajson: Dict[str, Any],
    indir: str,
    outdir: str,
    raw_files: List[str],
    workflow: Optional[str] = None
) -> str:
    """Set common fields for raw.json: ROOT_DIR, indir, outdir, logdir, raw_files.

    Also detects the current server and rewrites any paths that reference a
    foreign server's database directory to use the current server's directory.
    """
    root_dir = os.path.dirname(__file__)
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    tmp_dir = datajson.get("tmp_dir") or os.path.join(outdir, "tmp")
    if not datajson.get("tmp_dir"):
        os.makedirs(tmp_dir, exist_ok=True)
    # Common fields for raw.json: ROOT_DIR, Params.workflow.{indir,outdir,
    # logdir,raw_files}, tmp_dir. Workflow-level fields live under
    # Params.workflow; top level carries only outfiles/ROOT_DIR/tmp_dir/
    # resources/env/Procedure/genome/Params.
    datajson["ROOT_DIR"] = root_dir
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    wf["indir"] = indir
    wf["outdir"] = outdir
    wf["logdir"] = logdir
    wf["raw_files"] = raw_files
    datajson["tmp_dir"] = tmp_dir

    # --- server database path migration ---
    current_db = _detect_current_server_db()
    if current_db is not None:
        # Build foreign dirs list, longest prefix first to avoid partial matches
        foreign_dirs = sorted(
            (d for d in SERVER_DATABASE_DIR.values() if d != current_db),
            key=len, reverse=True,
        )
        if foreign_dirs:
            _replace_server_paths(datajson, foreign_dirs, current_db)

    return root_dir

def _classify_layouts(samples_info_dict: Dict[str, SampleInfo]) -> tuple[List[str], List[str]]:
    """Classify samples into paired/single by layout."""
    paired: List[str] = []
    single: List[str] = []
    for sid, info in samples_info_dict.items():
        if info.layout == "PE":
            paired.append(sid)
        elif info.layout == "SE":
            single.append(sid)
        else:
            logger.error(f"Unknown layout type for sample {sid}: {info.layout}")
    return paired, single

def _write_raw_json(datajson: Dict[str, Any], outdir: str) -> str:
    """Write raw.json and return its path."""
    path = os.path.join(outdir, "raw.json")
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(datajson, f, indent=2, ensure_ascii=False)
    return path

def runCoCulture(
        datajson: Dict[str,Any],
        samples_info_dict:Dict[str, SampleInfo],
        indir:str,
        outdir: str,
        raw_files: List[str],
    ):
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)
    outfiles = []
    paired_samples = []
    single_samples = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            outfiles.append(f"{outdir}/common/2_trimmed_fastq/{sample_id}/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/common/2_trimmed_fastq/{sample_id}/{sample_id}_2.fq.gz")
            outfiles.append(f"{outdir}/common/3_raw_bam/GRCh38/{sample_id}/{sample_id}.bam")
            outfiles.append(f"{outdir}/common/3_raw_bam/GRCm39/{sample_id}/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount_name.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount_name.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/common/2_trimmed_fastq/{sample_id}/{sample_id}.single.fq.gz")
            outfiles.append(f"{outdir}/common/3_raw_bam/GRCm39/{sample_id}/{sample_id}.bam")
            outfiles.append(f"{outdir}/common/3_raw_bam/GRCh38/{sample_id}/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount_name.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount_name.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    outfiles.append(f"{outdir}/disambiguate/disambiguate_qc.tsv")
    # dual-genome co-culture: every sample is aligned against BOTH references
    # (disambiguate later separates organism-specific reads), so the per-genome
    # sample maps used by the polygenomes aligner cover all samples per genome
    _co_genomes = ("GRCh38", "GRCm39")
    datajson["Params"].setdefault("workflow", {})["genome_paired_samples"] = {g: list(paired_samples) for g in _co_genomes}
    datajson["Params"].setdefault("workflow", {})["genome_single_samples"] = {g: list(single_samples) for g in _co_genomes}
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    wf["paired_samples"] = paired_samples
    wf["single_samples"] = single_samples
    return _write_raw_json(datajson, outdir)

def runMERIP(
        datajson: Dict[str, Any],
        samples_info_dict:Dict[str, SampleInfo],
        indir:str,
        outdir: str,
        raw_files: List[str],
    ):
    """
    Function: Prepare input JSON for MERIP workflow based on the provided model JSON template and sample information.
    Parameters:
    - input_json: Path to the model JSON template file.
    - samples_info_dict: A dictionary containing sample information, where keys are sample IDs and values are objects with attributes 'layout' and 'design'.
    - indir: Input directory containing raw data (e.g., FASTQ files).
    - outdir: Output directory where results will be stored.
    Returns:
    - instance_json: Path to the generated input JSON file that will be used for the MERIP workflow.
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)

    paired_samples = []
    single_samples = []
    ip_samples = []
    input_samples = []
    treated_ip_samples = []
    treated_input_samples = []
    outfiles = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}_1.fq.gz")
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}_2.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/{sample_id}/{sample_id}.dedup.bam")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}.single.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/{sample_id}/{sample_id}.dedup.bam")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
        
        if sample_info.design == "ip":
            ip_samples.append(sample_id)
        elif sample_info.design == "input":
            input_samples.append(sample_id)
        elif sample_info.design == "treated_ip":
            treated_ip_samples.append(sample_id)
        elif sample_info.design == "treated_input":
            treated_input_samples.append(sample_id)
        else:
            logger.error(f"Unknown design type for sample {sample_id}: {sample_info.design}")
    outfiles.append(f"{outdir}/exomePeak/sig_diff_peak_gene_names.xls")
    outfiles += _register_regionplot(datajson, samples_info_dict,
                                     ip_samples + input_samples + treated_ip_samples + treated_input_samples,
                                     outdir, "MERIP")
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    wf["paired_samples"] = paired_samples
    wf["single_samples"] = single_samples
    wf["ip_samples"] = ip_samples
    wf["input_samples"] = input_samples
    wf["treated_ip_samples"] = treated_ip_samples
    wf["treated_input_samples"] = treated_input_samples
    return _write_raw_json(datajson, outdir)


def runCLIP(
        datajson: Dict[str, Any],
        samples_info_dict:Dict[str, SampleInfo],
        indir:str,
        outdir: str,
        raw_files: List[str],
    ):
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)
    outfiles = []
    paired_samples = []
    single_samples = []
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    # CLI --extra-args aligner=... lands at config top level and wins over
    # Params.workflow.aligner; normalize into Params.workflow for raw.json.
    aligner = datajson.pop("aligner", None) or wf.get("aligner") or "star"
    wf["aligner"] = aligner
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}_2.fq.gz")
            if aligner == "star":
                outfiles.append(f"{outdir}/star/{sample_id}/{sample_id}.bam")
            elif aligner == "hisat2":
                outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            outfiles.append(f"{outdir}/fastqc/raw/{sample_id}/{sample_id}.fastqc.raw.txt")
            outfiles.append(f"{outdir}/fastqc/trimmed/{sample_id}/{sample_id}.fastqc.trimmed.txt")
            outfiles.append(f"{outdir}/PureCLIP/{sample_id}.pureclip.sites.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.plus.bw")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.minus.bw")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}.single.fq.gz")
            if aligner == "star":
                outfiles.append(f"{outdir}/star/{sample_id}/{sample_id}.bam")
            elif aligner == "hisat2":
                outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            outfiles.append(f"{outdir}/fastqc/raw/{sample_id}/{sample_id}.fastqc.raw.txt")
            outfiles.append(f"{outdir}/fastqc/trimmed/{sample_id}/{sample_id}.fastqc.trimmed.txt")
            outfiles.append(f"{outdir}/PureCLIP/{sample_id}.pureclip.sites.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.plus.bw")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.minus.bw")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    outfiles.append(f"{outdir}/track/igv_track_iclip.html")
    outfiles.append(f"{outdir}/track/ucsc_track_iclip.txt")
    outfiles += _register_regionplot(datajson, samples_info_dict,
                                     paired_samples + single_samples,
                                     outdir, "CLIP")
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles
    wf["paired_samples"] = paired_samples
    wf["single_samples"] = single_samples
    # parameters suggest by https://doi.org/10.1016/j.ymeth.2019.11.008
    datajson["Params"]["bamCoverage"]["offset"] = "-1"
    datajson["Params"]["bamCoverage"]["binSize"] = 1
    # normalizeUsing 精确到 rule:bamCoverage/wig 各自子键
    # (Params.bamCoverage.<rule>.normalizeUsing;未迁移的 binSize/offset/
    # extendReads 仍是模块级扁平键)
    _bc = datajson["Params"].setdefault("bamCoverage", {})
    _bc.setdefault("bamCoverage", {})["normalizeUsing"] = "CPM"
    _bc.setdefault("wig", {})["normalizeUsing"] = "CPM"
    datajson["Params"]["bamCoverage"]["extendReads"] = 1
    # star 参数精确到 rule:写入 star_align 子段(Params.star.star_align.<key>,
    # 消费 rule star_align 经 rule_param 读取,无 module 级共享回退)
    _star_align = datajson["Params"].setdefault("star", {}).setdefault("star_align", {})
    _star_align["alignEndsType"] = "Extend5pOfRead1"
    _star_align["outFilterMismatchNoverReadLmax"] = 0.04
    _star_align["outFilterMismatchNmax"] = 999
    _star_align["outFilterMultimapNmax"] = 999
    datajson["Params"]["igv"]["js"] = "/data/pub/zhousha/Reference/igv.min.js"
    datajson["Params"]["igv"]["publicPathMap"] = {
        "/data/pub/zhousha/": "/data/",
        "/data/pub/zhousha/Reference/": "/ref/"
    }
    return _write_raw_json(datajson, outdir)

def runPacVar(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        group_pairs: List[CompareGroupPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ):
    """Prepare input JSON for PacVar (PacBio variant calling) workflow.

    Step gating is driven entirely by `Params.<step>.enabled` flags
    (snv, sv, phase, repeat, telomere, eccdna) and `Params.snv_caller`.
    The subworkflow declares every module unconditionally; this function
    decides which outfiles (and therefore which DAG leaves) are requested.

    `Params.sv.enabled` additionally gates the SV comparison analysis
    (modules/sv). `Params.sv.comparisons` follows the DESeq2 group_pairs
    design (name -> group names + sample name lists); it is filled from the
    meta group pairs when not hand-authored in the config.
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files, "PacVar")
    outfiles: List[str] = []
    samples: List[str] = []
    params = datajson.get("Params", {})

    snv_caller = params.get("snv_caller", "deepvariant")
    enable_snp = _enabled(datajson, "snv")
    enable_sv = _enabled(datajson, "sv")
    enable_phase = _enabled(datajson, "phase")
    enable_repeat = _enabled(datajson, "repeat")
    enable_telomere = _enabled(datajson, "telomere")
    enable_eccdna = _enabled(datajson, "decoil")  # Decoil is the eccDNA reconstruction step

    _genome_ref = _resolve_genome_ref(datajson, samples_info_dict, "PacVar")
    _genome_name = datajson["genome"]["default"]
    datajson["Params"].setdefault("sv", {}).setdefault("exp_specific", {})
    datajson["Params"]["sv"]["exp_specific"]["species"] = _genome_name
    datajson["Params"]["sv"]["exp_specific"]["assembly"] = _genome_name
    for sample_id, sample_info in samples_info_dict.items():
        samples.append(sample_id)
        # SNP calling
        if enable_snp:
            if snv_caller == "deepvariant":
                outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.vcf.gz")
                outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.vcf.gz.csi")
            elif snv_caller == "gatk4":
                outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.filtered.vcf.gz")
                outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.filtered.vcf.gz.csi")
        # SV calling
        if enable_sv:
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.vcf.gz")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.vcf.gz.csi")
        # phasing
        if enable_phase and enable_snp and enable_sv:
            outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.phased.vcf.gz")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.phased.vcf.gz")
        # repeat characterization
        if enable_repeat and _genome_ref.get("repeat_bed"):
            outfiles.append(f"{outdir}/repeat/trgt/genotype/{sample_id}/{sample_id}.trgt.vcf.gz")
            _trgt_fmts = (datajson.get("Params", {}).get("trgt", {})
                          .get("trgt_plot", {}).get("formats") or ["png"])
            if isinstance(_trgt_fmts, str):
                _trgt_fmts = [_trgt_fmts]
            outfiles += _figs(
                f"{outdir}/repeat/trgt/plot/{sample_id}/{sample_id}.trgt.repeat",
                _trgt_fmts)

    # eccDNA / ecDNA reconstruction via Decoil
    # Decoil is reconstruction-only: SV VCF comes from pbsv upstream and the
    # coverage bigwig from the subworkflow's `decoil_coverage` rule (see
    # subworkflow/PacVar.smk). Only declare decoil's own reconstruction
    # outputs here.
    genome_gtf = _genome_ref.get("gtf")
    if enable_eccdna and genome_gtf:
        for sample_id in samples:
            outfiles.append(f"{outdir}/eccdna/{sample_id}/decoil/reconstruct.bed")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/decoil/reconstruct.ecDNA.bed")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/decoil/summary.txt")
            # Visualization: per-sample stats + interactive IGV HTML
            _decoil_plot_fmts = (datajson.get("Params", {}).get("decoil", {})
                                 .get("plot", {}).get("formats") or ["png"])
            if isinstance(_decoil_plot_fmts, str):
                _decoil_plot_fmts = [_decoil_plot_fmts]
            outfiles += _figs(
                f"{outdir}/eccdna/{sample_id}/visualization/decoil_stats",
                _decoil_plot_fmts)
            _decoil_structure_fmts = (datajson.get("Params", {}).get("decoil", {})
                                      .get("structure", {}).get("formats") or ["png"])
            if isinstance(_decoil_structure_fmts, str):
                _decoil_structure_fmts = [_decoil_structure_fmts]
            outfiles += _figs(
                f"{outdir}/eccdna/{sample_id}/visualization/decoil_structure",
                _decoil_structure_fmts)
            outfiles.append(f"{outdir}/eccdna/{sample_id}/visualization/decoil_structure_genes.tsv")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/visualization/igv_ecdna.html")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/visualization/ucsc_track.txt")
        # Cross-sample comparison (aggregated): decoil module writes these
        # under {outdir}/eccdna/compare/ (rule decoil_compare), not
        # {outdir}/eccdna/visualization/.
        _decoil_compare_fmts = (datajson.get("Params", {}).get("decoil", {})
                                .get("compare", {}).get("formats") or ["png"])
        if isinstance(_decoil_compare_fmts, str):
            _decoil_compare_fmts = [_decoil_compare_fmts]
        outfiles += _figs(f"{outdir}/eccdna/compare/decoil_compare",
                          _decoil_compare_fmts)
        outfiles.append(f"{outdir}/eccdna/compare/decoil_compare_shared_unique.tsv")
    elif enable_eccdna:
        logger.warning(
            "PacVar: Params.eccdna.enabled=true but genome.gtf is missing; "
            "Decoil requires a gene annotation GTF. Skipping eccDNA outfiles."
        )

    # telomere & centromere analysis
    if enable_telomere:
        # Which telomere-length methods to run (Params.telomere.method, list).
        # The telomere module declares every rule unconditionally; this list
        # decides which outputs (DAG leaves) are requested, so Snakemake runs
        # only the selected methods (tool choice is a workflow concern, not the
        # module's). Empty/absent => all four.
        _TEL_METHODS = ("telogator2", "assembly_scan", "read_density", "tidk")
        tel_methods = params.get("telomere", {}).get("method") or list(_TEL_METHODS)
        if isinstance(tel_methods, str):
            tel_methods = [tel_methods]
        _bad = [m for m in tel_methods if m not in _TEL_METHODS]
        if _bad:
            raise ValueError(
                f"PacVar: Params.telomere.method has unknown method(s) {_bad}; "
                f"valid methods: {list(_TEL_METHODS)}."
            )
        # figure formats (Params.telomere.formats), for multi-format plot outputs
        _tel_fmts = params.get("telomere", {}).get("formats") or ["png"]
        if isinstance(_tel_fmts, str):
            _tel_fmts = [_tel_fmts]
        for sample_id in samples:
            # Telogator2 (per-chromosome-arm)
            if "telogator2" in tel_methods:
                outfiles.append(f"{outdir}/telomere/length/telogator2/{sample_id}/tlens_by_allele.tsv")
                outfiles.append(f"{outdir}/telomere/length/telogator2/{sample_id}/all_final_alleles.png")
                outfiles.append(f"{outdir}/telomere/length/telogator2/{sample_id}/violin_atl.png")
            # Approach A: Assembly contig end scanning
            if "assembly_scan" in tel_methods:
                outfiles.append(f"{outdir}/telomere/length/assembly_scan/{sample_id}/{sample_id}_assembly_telomere_stats.txt")
                for _f in _tel_fmts:
                    outfiles.append(f"{outdir}/telomere/length/assembly_scan/{sample_id}/{sample_id}_assembly_telomere_length.{_f}")
            # Approach B: Read-level k-mer density
            if "read_density" in tel_methods:
                outfiles.append(f"{outdir}/telomere/length/read_density/{sample_id}/{sample_id}_read_telomere_stats.txt")
                for _f in _tel_fmts:
                    outfiles.append(f"{outdir}/telomere/length/read_density/{sample_id}/{sample_id}_read_telomere_content.{_f}")
            # Approach C: tidk
            if "tidk" in tel_methods:
                outfiles.append(f"{outdir}/telomere/length/tidk/{sample_id}/{sample_id}_tidk_telomeres.tsv")
                for _f in _tel_fmts:
                    outfiles.append(f"{outdir}/telomere/length/tidk/{sample_id}/{sample_id}_tidk_telomeres.{_f}")
            if "read_density" in tel_methods:
                for _f in _tel_fmts:
                    outfiles.append(f"{outdir}/telomere/length/read_density/combined.{_f}")
            if len(tel_methods) >= 4:
                for _f in _tel_fmts:
                    outfiles.append(f"{outdir}/telomere/length/combined_methods.{_f}")
            _dq_bed = ((params.get("telomere", {}).get("depth_qc", {})
                        or {}).get("telomere_bed"))
            if _dq_bed:
                outfiles.append(f"{outdir}/telomere/length/depth_qc/{sample_id}/{sample_id}_telomere_region_depth_stats.txt")
            # composition + regions 子模块输出
            outfiles.append(f"{outdir}/telomere/composition/telomere_summary.tsv")
            outfiles.append(f"{outdir}/telomere/composition/telomere_composition.tsv")
            for _f in _tel_fmts:
                outfiles.append(f"{outdir}/telomere/composition/telomere_length_violin.{_f}")
            if "assembly_scan" in tel_methods:
                outfiles.append(f"{outdir}/telomere/regions/{sample_id}/{sample_id}.telomere.bed")
            # Centromere (separate module; not gated by telomere method)
            outfiles.append(f"{outdir}/repeat/centromere/{sample_id}/{sample_id}.centromere_stats.txt")

    # SV comparison analysis (modules/sv): DESeq2-style design —
    # Params.sv.comparisons maps "{control_group_name}_vs_{experimental_group_name}"
    # to {control_group_name, experimental_group_name, control_samples,
    # experimental_samples}. Sample names only: the sv rules resolve each
    # sample's SV VCF from vcf_dir/vcf_suffix. A hand-filled config block
    # wins; otherwise it is derived from the meta group pairs, like
    # runRNAseq fills Params.DESeq2.group_pairs.
    sv_params = datajson.get("Params", {}).get("sv", {})
    _validate_te_lib(datajson, "PacVar")
    _fill_mei_rm_species(datajson, "PacVar")
    _fill_telomere_species(datajson, "PacVar")
    _fill_telomere_subtel_ref(datajson, "PacVar")
    _fill_telomere_n_chrom_arms(datajson, "PacVar")
    _fill_repeatmasker_species(datajson, "PacVar")
    comparisons = dict(sv_params.get("comparisons") or {})
    if enable_sv and not comparisons:
        for gp in group_pairs:
            name = f"{gp.ctr_group_name}_vs_{gp.exp_group_name}"
            comparisons[name] = {
                "control_group_name": gp.ctr_group_name,
                "experimental_group_name": gp.exp_group_name,
                "control_samples": gp.ctr_sample_ids,
                "experimental_samples": gp.exp_sample_ids,
            }
    datajson["Params"].setdefault("sv", {})["comparisons"] = comparisons
    if enable_sv and not comparisons:
        logger.warning(
            "PacVar: Params.sv.enabled=true but no SV comparisons are defined; "
            "fill Params.sv.comparisons in the config JSON or the meta design "
            "column (ctrl_TAG / exp_TAG). Skipping SV comparison outputs."
        )
    if enable_sv and comparisons:
        sv_outdir = f"{outdir}/variation/germline_sv_analysis"
        sv_genes = (datajson.get("Params", {}).get("sv", {})
                    .get("sv_gene_model", {}).get("genes") or [])
        sv_oncoprint = sv_params.get("oncoprint") or {}
        known_samples = set(samples_info_dict)
        for name, c in comparisons.items():
            for key in ("control_samples", "experimental_samples"):
                unknown = [s for s in (c.get(key) or []) if s not in known_samples]
                if unknown:
                    logger.warning(
                        f"PacVar: Params.sv.comparisons['{name}'].{key} lists "
                        f"samples not in the meta: {unknown}"
                    )
            base = f"{sv_outdir}/{name}"
            outfiles.append(f"{base}/{name}_merged.vcf")
            outfiles.append(f"{base}/{name}_only.vcf")
            outfiles.append(f"{base}/{name}_annotated.tab")
            # TE enrichment (RepeatMasker subfamily test) + per-SV stats from
            # sv_exp_specific — registered so the rule's declared directory
            # outputs are requested by the DAG
            outfiles.append(f"{base}/repeatmasker")
            outfiles.append(f"{base}/statistics")
            outfiles.append(f"{base}/circos")
            # MEI suite (sv_mei rule): classified insertions + 6-figure
            # suite; gate via Params.sv.mei.enabled (default on)
            if (datajson.get("Params", {}).get("sv", {}).get("mei", {}) or {}).get("enabled", True):
                outfiles.append(f"{base}/mei/mei.done")
            if sv_genes:
                outfiles.append(f"{base}/gene_model")
        outfiles.append(f"{sv_outdir}/sv_diff_analysis")
        if sv_oncoprint.get("cosmic_file"):
            outfiles.append(f"{sv_outdir}/sv_oncoprint/oncoprint_matrix.csv")
        # SV-affected genes GO/KEGG (modules/function/sv_go; macaque-aware
        # clusterProfiler route). Registered only when the gate is on.
        svgokegg_params = datajson.get("Params", {}).get("function", {}).get("sv_go_kegg", {}) or {}
        if svgokegg_params.get("enabled"):
            svgokegg_species = resolve_genome(_genome_name,"alias")
            if svgokegg_species not in ("human", "mouse", "macaque"):
                raise ValueError(
                    "PacVar: Params.function.sv_go_kegg.enabled=true requires "
                    "Params.function.sv_go_kegg.species to be one of "
                    f"human/mouse/macaque (got: {svgokegg_species!r})."
                )
            datajson["Params"]["function"]["sv_go_kegg"]["species"] = svgokegg_species
            for name in comparisons:
                outfiles.append(f"{sv_outdir}/{name}/sv_go_kegg/sv_go_kegg.done")

    gatk_tmp_dir = os.path.join(outdir, "tmp")
    os.makedirs(gatk_tmp_dir, exist_ok=True)
    # tmp-dir 精确到 rule:镜像写入全部消费 rule 的子段
    # (Params.gatk.<rule>.tmp-dir;无 module 级共享回退)
    _gatk = datajson["Params"].setdefault("gatk", {})
    for _grule in ("gatk_index", "addReadsGroup", "MarkDuplicates",
                   "collectHsMetrics", "BaseRecalibrator", "ApplyBQSR",
                   "HaplotypeCaller", "filterHaplotypeCallerVcf",
                   "SplitNCigarReads", "VarientCalling", "vcf_filter",
                   "somaticMutect2", "filterMutectCalls",
                   "gatk_population_haplotype_caller",
                   "gatk_population_joint_genotyping"):
        _gatk.setdefault(_grule, {})["tmp-dir"] = gatk_tmp_dir
    datajson["Params"]["workflow"]["samples"] = samples
    outfiles += _register_regionplot(datajson, samples_info_dict, samples, outdir, "PacVar")
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles
    return _write_raw_json(datajson, outdir)

def runMutation(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        designPairs: List[DesignPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ):
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files, workflow="Mutation")
    _genome_ref = _resolve_genome_ref(datajson, samples_info_dict, "Mutation")
    # Sequencing assay type -> CNVkit --seq-method (meta `assay` column wins over
    # the config value; all samples must share one assay for a single cnvkit run).
    _ASSAY_TO_CNVKIT_METHOD = {
        "wgs": "wgs", "exome": "hybrid", "wes": "hybrid",
        "panel": "amplicon", "amplicon": "amplicon",
    }
    assays = sorted({(info.assay or "").strip().lower() for info in samples_info_dict.values()} - {""})
    if assays:
        if len(assays) > 1:
            raise ValueError(
                f"Mutation: meta assay column holds mixed assay types {assays}; "
                "CNVkit needs a single sequencing method per run."
            )
        method = _ASSAY_TO_CNVKIT_METHOD.get(assays[0])
        if method is None:
            raise ValueError(
                f"Mutation: unknown meta assay '{assays[0]}'; expected one of "
                f"{sorted(_ASSAY_TO_CNVKIT_METHOD)}"
            )
        configured = datajson["Params"].setdefault("cnvkit", {}).get("method")
        if configured != method:
            logger.info(
                f"Mutation: meta assay '{assays[0]}' -> Params.cnvkit.method "
                f"'{configured}' overridden to '{method}'"
            )
        datajson["Params"]["cnvkit"]["method"] = method
    outfiles = []
    paired_samples = []
    single_samples = []
    mutect2_samples = []
    sample_somatic_vcf_dict = {}
    sample_group_dict = {}
    paired_design = []
    # Exome assays get a capture-region BED for GATK -L and CNVkit --targets:
    # an explicit genome.references.<g>.targets BED wins; otherwise derive it
    # from markdup BAM coverage (upstream rule) and register the output.
    _ref_entry = datajson["genome"]["references"][datajson["genome"]["default"]]
    if _ref_entry.get("targets"):
        _targets = _ref_entry["targets"]
    elif assays and datajson["Params"]["cnvkit"].get("method") == "hybrid":
        _targets = f"{outdir}/common/5_targets_bed/targets.bed"
        outfiles.append(_targets)
        _ref_entry["targets"] = _targets
    else:
        _targets = None
        _ref_entry["targets"] = None
    for designPair in designPairs:
        somatic_file = f"{outdir}/variation/somatic_snv_indel/{designPair.ctr_sample_id}_vs_{designPair.exp_sample_id}/{designPair.ctr_sample_id}_vs_{designPair.exp_sample_id}.filtered.vcf.gz"
        sample_somatic_vcf_dict[designPair.exp_sample_id] = somatic_file
        sample_group_dict[designPair.exp_sample_id] = designPair.exp_group
        outfiles.append(somatic_file)
        mutect2_samples.append(designPair.ctr_sample_id)
        mutect2_samples.append(designPair.exp_sample_id)
        paired_design.append({"control": designPair.ctr_sample_id, "experimental": designPair.exp_sample_id})

    datajson["Params"]["workflow"]["paired_design"] = paired_design
    # Strain-aware germline filtering (meta ``strain`` column -> per-pair strain
    # map for gatk_somatic; empty when the column is absent/left blank).
    pair_strains = {}
    for designPair in designPairs:
        strain = getattr(samples_info_dict.get(designPair.exp_sample_id), "strain", None)
        if strain:
            pair_strains[designPair.exp_sample_id] = strain
    datajson["Params"]["workflow"]["pair_strains"] = pair_strains
    if pair_strains:
        logger.info(f"Mutation: strain-aware germline filtering enabled for pairs: {pair_strains}")
    # CNVkit calling mode resolution (Params.cnvkit.mode; "auto" decides here
    # from the design): matched normals in the design pairs -> paired; else a
    # normal panel (Params.workflow.control_samples) -> pooled; else flat.
    # An explicit mode wins over "auto"; impossible combos hard-error below.
    _cnv_cfg = datajson["Params"].setdefault("cnvkit", {})
    _paired_map = {
        rec["experimental"]: rec["control"]
        for rec in paired_design
        if rec.get("experimental") and rec.get("control")
    }
    _cnv_mode = _cnv_cfg.get("mode") or "auto"
    if _cnv_mode == "auto":
        _cnv_mode = "paired" if _paired_map else (
            "pooled" if (datajson["Params"]["workflow"].get("control_samples") or []) else "flat"
        )
        _cnv_cfg["mode"] = _cnv_mode
        logger.info(f"Mutation: cnvkit mode auto -> {_cnv_mode} "
                    f"(design pairs: {len(_paired_map)}, "
                    f"control_samples: {len(datajson['Params']['workflow'].get('control_samples') or [])})")
    if _cnv_mode == "paired" and not _paired_map:
        raise ValueError(
            "Mutation: Params.cnvkit.mode 'paired' requires matched normal "
            "samples in the meta design; none found. Use mode 'pooled' (with "
            "control_samples) or 'flat'."
        )
    if _cnv_mode == "pooled" and not (datajson["Params"]["workflow"].get("control_samples") or []):
        raise ValueError(
            "Mutation: Params.cnvkit.mode 'pooled' requires "
            "Params.workflow.control_samples; none found. Use mode 'flat' "
            "or provide a normal panel."
        )
    if _cnv_mode == "paired":
        # Paired mode must not also build a pooled reference alongside the
        # per-pair references (conflicting reference semantics).
        datajson["Params"]["workflow"]["control_samples"] = []
    # 参数精确到 rule:resolve 结果(以及 meta assay 推导的 method)镜像写入
    # rule 子段——modules/cnvkit 只读 Params.cnvkit.<rule>.mode/.method
    # (rule_param,无 module 级回退),扁平键仅留作 node.py 自身门控/兼容。
    for _cnv_rule in ("cnvkit_reference", "cnvkit_batch", "cnvkit_paired"):
        _cnv_sub = _cnv_cfg.setdefault(_cnv_rule, {})
        _cnv_sub["mode"] = _cnv_mode
        if _cnv_cfg.get("method"):
            _cnv_sub["method"] = _cnv_cfg["method"]
    for sample_id, sample_info in samples_info_dict.items():

        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            outfiles.append(f"{outdir}/common/3_sorted_bam/{sample_id}/{sample_id}_flagstat.txt")
            if sample_id in mutect2_samples:
                logger.info(f"Sample {sample_id} is involved in mutect2 analysis, skipping germline workflow for this sample.")
                continue
            outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.filtered.vcf.gz")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/common/3_sorted_bam/{sample_id}/{sample_id}_flagstat.txt")
            if sample_id in mutect2_samples:
                logger.info(f"Sample {sample_id} is involved in mutect2 analysis, skipping germline workflow for this sample.")
                continue
            outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.filtered.vcf.gz")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    if designPairs:
        _spectrum_fmts = (datajson.get("Params", {}).get("somatic_spectrum", {})
                          .get("plot", {}).get("formats") or ["png"])
        if isinstance(_spectrum_fmts, str):
            _spectrum_fmts = [_spectrum_fmts]
        outfiles += _figs(f"{outdir}/stability/spectrum/somatic_spectrum_stacked_bar",
                          _spectrum_fmts)
        outfiles += _figs(f"{outdir}/stability/spectrum/somatic_spectrum_heatmap",
                          _spectrum_fmts)
        outfiles.append(f"{outdir}/stability/spectrum/somatic_spectrum.csv")
        outfiles.append(f"{outdir}/stability/spectrum/somatic_spectrum_per_sample.csv")
        outfiles.append(f"{outdir}/stability/spectrum/somatic_spectrum_sample_stats.csv")
        # Signature analysis (modules/spectrum): refit against COSMIC and/or
        # de novo NMF. Gated by Params.somatic_spectrum.signature.enabled
        # (default off); mode selects which outputs are registered.
        _sig_cfg = datajson["Params"].get("somatic_spectrum", {}).get("signature", {}) or {}
        if _sig_cfg.get("enabled"):
            _sig_mode = str(_sig_cfg.get("mode", "both")).lower()
            if _sig_mode in ("refit", "both"):
                _sig_fmts = (datajson.get("Params", {}).get("somatic_spectrum", {})
                             .get("signature", {}).get("formats") or ["png"])
                if isinstance(_sig_fmts, str):
                    _sig_fmts = [_sig_fmts]
                for _f in ("signature_refit_exposures.csv",
                           "signature_refit_exposure_fractions.csv",
                           "signature_refit_family_exposures.csv",
                           "signature_refit_summary.csv"):
                    outfiles.append(f"{outdir}/stability/spectrum/{_f}")
                for _b in ("signature_refit_exposure_stacked_bar",
                           "signature_refit_exposure_heatmap",
                           "signature_refit_signature_cosine"):
                    outfiles += _figs(f"{outdir}/stability/spectrum/{_b}",
                                      _sig_fmts)
            if _sig_mode in ("nmf", "both"):
                _sig_fmts = (datajson.get("Params", {}).get("somatic_spectrum", {})
                             .get("signature", {}).get("formats") or ["png"])
                if isinstance(_sig_fmts, str):
                    _sig_fmts = [_sig_fmts]
                for _f in ("signature_nmf_signatures.csv",
                           "signature_nmf_exposures.csv",
                           "signature_nmf_rank_diagnostics.csv",
                           "signature_nmf_cosmic_similarity.csv",
                           "signature_nmf_cosmic_best_match.csv"):
                    outfiles.append(f"{outdir}/stability/spectrum/{_f}")
                for _b in ("signature_nmf_rank_diagnostics",
                           "signature_nmf_exposures_stacked_bar",
                           "signature_nmf_cosmic_match"):
                    outfiles += _figs(f"{outdir}/stability/spectrum/{_b}",
                                      _sig_fmts)
    
    all_samples = paired_samples + single_samples
    
    # Hybrid-selection QC (modules/gatk collectHsMetrics): per-sample
    # CollectHsMetrics on the markdup BAM over the capture region. Requires
    # a targets BED/interval_list (exome/hybrid runs; node.py fills the
    # reference entry above); WGS runs resolve _targets to None and skip.
    if _targets:
        for sample_id in all_samples:
            outfiles.append(f"{outdir}/common/4_markdup_bam/{sample_id}/{sample_id}.hs_metrics.txt")
    
    # Per-module gates: Params.<module>.enabled (default true)
    if _enabled(datajson, "fragment_size"):
        outfiles.append(f"{outdir}/results/fragment_size/fragment/FragmentSize.txt")
        _frag_fmts = (datajson.get("Params", {}).get("fragment_size", {})
                      .get("plotFragmentSize", {}).get("formats") or ["png"])
        if isinstance(_frag_fmts, str):
            _frag_fmts = [_frag_fmts]
        outfiles += _figs(f"{outdir}/results/fragment_size/fragment/FragmentSize",
                          _frag_fmts)

    if _enabled(datajson, "manta"):
        # Matched-normal design pairs -> Manta somatic mode (somaticSV.vcf.gz),
        # mirroring the somatic SNV (Mutect2) pair judgment above; paired
        # samples are skipped from germline SV exactly like mutect2_samples.
        for designPair in designPairs:
            outfiles.append(f"{outdir}/variation/somatic_sv/{designPair.ctr_sample_id}_vs_{designPair.exp_sample_id}/results/variants/somaticSV.vcf.gz")
        for sample_id in all_samples:
            if sample_id in mutect2_samples:
                logger.info(f"Sample {sample_id} is involved in manta somatic analysis, skipping germline SV workflow for this sample.")
                continue
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/results/variants/candidateSV.vcf.gz")

    if _enabled(datajson, "msi") and designPairs:
        for designPair in designPairs:
            outfiles.append(
                f"{outdir}/stability/msi/{designPair.exp_sample_id}/"
                f"{designPair.exp_sample_id}_msi_pro.tsv")
        outfiles.append(f"{outdir}/stability/msi/msi_input.txt")
        # formats 精确到 rule:Params.msi.msi_plot.formats(config 全量路径读取,
        # 无 module 级共享回退)。
        _msi_fmts = (datajson.get("Params", {}).get("msi", {})
                     .get("msi_plot", {}).get("formats") or ["jpg"])
        if isinstance(_msi_fmts, str):
            _msi_fmts = [_msi_fmts]
        outfiles += _figs(f"{outdir}/stability/msi/msi_barplot", _msi_fmts)

    if _enabled(datajson, "cnvkit"):
        # Mode-aware registration: paired -> per-pair tumors only (no pooled
        # reference); pooled -> all samples + shared reference.cnn; flat ->
        # all samples, no reference.
        if _cnv_mode == "paired":
            for sample_id in _paired_map:
                outfiles.append(f"{outdir}/variation/cnv/{sample_id}/{sample_id}.cnr")
                outfiles.append(f"{outdir}/variation/cnv/{sample_id}/{sample_id}.cns")
        else:
            for sample_id in all_samples:
                outfiles.append(f"{outdir}/variation/cnv/{sample_id}/{sample_id}.cnr")
                outfiles.append(f"{outdir}/variation/cnv/{sample_id}/{sample_id}.cns")
            if _cnv_mode == "pooled":
                outfiles.append(f"{outdir}/variation/cnv/reference/reference.cnn")
    
    # Mutation burden (modules/burden): SNV/InDel TMB + SV + CNV burden for
    # each tumor sample. The three sample -> file maps derive from the somatic
    # and CNV outputs above (hand-filled Params.burden entries win); the CNV
    # sample set follows the resolved CNVkit mode (paired -> tumors only).
    if _enabled(datajson, "burden"):
        burden_cfg = datajson["Params"].setdefault("burden", {})
        snv_vcf_dict = dict(sample_somatic_vcf_dict)  # {tumor: somatic SNV VCF}
        sv_vcf_dict = {
            designPair.exp_sample_id:
                f"{outdir}/variation/somatic_sv/{designPair.ctr_sample_id}_vs_{designPair.exp_sample_id}/results/variants/somaticSV.vcf.gz"
            for designPair in designPairs
        }
        if _cnv_mode == "paired":
            burden_cnv_samples = list(_paired_map)
        else:
            burden_cnv_samples = list(all_samples)
        cnv_dict = {
            sample_id: f"{outdir}/variation/cnv/{sample_id}/{sample_id}.cns"
            for sample_id in burden_cnv_samples
        }
        for key, derived in (("sample_snv_vcf_dict", snv_vcf_dict),
                             ("sample_sv_vcf_dict", sv_vcf_dict),
                             ("sample_cnv_dict", cnv_dict)):
            if not burden_cfg.get(key):
                burden_cfg[key] = derived
        # VEP annotation settings (Params.burden.burden_snv_indel.annotate,
        # default on;config 全量路径读取,无 module 级共享回退):
        # species/assembly are derived from the genome; an explicit config
        # value wins and an unmapped genome hard-errors instead of guessing.
        _burden_snvi = burden_cfg.setdefault("burden_snv_indel", {})
        if _burden_snvi.get("annotate", True):
            if not _burden_snvi.get("vep_species") \
                    or not _burden_snvi.get("vep_assembly"):
                _genome_key = str(datajson["genome"]["default"])
                mapped = GENOME_VEP_SPECIES.get(_genome_key)
                if not mapped:
                    raise ValueError(
                        f"Mutation: genome '{_genome_key}' has no known VEP "
                        f"species/assembly for the burden annotation; set "
                        f"'Params.burden.burden_snv_indel.vep_species' and "
                        f"'Params.burden.burden_snv_indel.vep_assembly' in the config JSON, or "
                        f"disable annotation with Params.burden.burden_snv_indel.annotate=false."
                    )
                _burden_snvi.setdefault("vep_species", mapped[0])
                _burden_snvi.setdefault("vep_assembly", mapped[1])
                logger.info(
                    f"Mutation: burden VEP settings auto-filled for genome "
                    f"'{_genome_key}': species={mapped[0]} assembly={mapped[1]}"
                )
        for sample_id in burden_cfg["sample_snv_vcf_dict"]:
            outfiles.append(f"{outdir}/stability/burden/{sample_id}/{sample_id}.snv_indel_burden.tsv")
        for sample_id in burden_cfg["sample_sv_vcf_dict"]:
            outfiles.append(f"{outdir}/stability/burden/{sample_id}/{sample_id}.sv_burden.tsv")
        for sample_id in burden_cfg["sample_cnv_dict"]:
            outfiles.append(f"{outdir}/stability/burden/{sample_id}/{sample_id}.cnv_burden.tsv")
        if (burden_cfg["sample_snv_vcf_dict"] or burden_cfg["sample_sv_vcf_dict"]
                or burden_cfg["sample_cnv_dict"]):
            outfiles.append(f"{outdir}/stability/burden/burden_summary.tsv")

    datajson["Params"]["somatic_spectrum"]["sample_somatic_vcf_dict"] = sample_somatic_vcf_dict
    datajson["Params"]["somatic_spectrum"]["sample_group_dict"] = sample_group_dict

    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")

    # Final PPT report + results workbook (modules/Mutation_report), gated
    # by Params.report.enabled like the other report modules.
    if _enabled(datajson, "report"):
        # 写回精确到 rule:Params.report.generate_report.date(与
        # runPeakCalling/runRIPseq/runRNAseq 的写回同一键形)
        datajson["Params"].setdefault("report", {}).setdefault(
            "generate_report", {}).setdefault(
            "date", time.strftime("%Y-%m-%d", time.localtime()))
        outfiles.append(f"{outdir}/Mutation_report.pptx")
        outfiles.append(f"{outdir}/Mutation_report.xlsx")
        # 报告自有图导出清单(报告模块规范:figures/<name>.<fmt> 条件产出,
        # 以 manifest.tsv 锚定 DAG 与磁盘一致)
        outfiles.append(f"{outdir}/figures/manifest.tsv")

    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    datajson["outfiles"] = outfiles
    datajson["Params"]["workflow"]["paired_samples"] = paired_samples
    datajson["Params"]["workflow"]["single_samples"] = single_samples
    datajson["Params"]["workflow"]["control_samples"] = datajson["Params"]["workflow"].get("control_samples", []) or []
    return _write_raw_json(datajson, outdir)

def runKARRseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ):
    """Prepare input JSON for KARRseq (Kethoxal-Assisted RNA-RNA interaction sequencing) workflow."""
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)

    paired_samples = []
    single_samples = []
    outfiles = []

    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            # Final ligation pairs
            outfiles.append(f"{outdir}/chimeric/{sample_id}/{sample_id}.dedup.ligation.pairs.gz")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/chimeric/{sample_id}/{sample_id}.dedup.ligation.pairs.gz")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")

    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    wf["paired_samples"] = paired_samples
    wf["single_samples"] = single_samples
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles

    return _write_raw_json(datajson, outdir)

def runPeakCalling(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        design_pairs:List[DesignPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ):
    """Prepare input JSON for PeakCalling (ChIP-seq/DIP-seq peak calling) workflow.
    
    Workflow steps:
    1. FastQC (raw)
    2. Trimming (trim_galore)
    3. FastQC (trimmed)
    4. Bowtie2 alignment
    5. AddReadsGroup + MarkDuplicates (GATK4)
    6. BigWig tracks (bamCoverage)
    7. MACS3 peak calling
    8. FRiP score
    9. deeptools enrichment heatmap
    10. HOMER peak annotation
    
    Supports both ChIP-seq and DIP-seq experiments.
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)

    ip_samples = []
    input_samples = []
    sample_ip_input_map = {}
    outfiles = []
    for design_pair in design_pairs:
        sample_ip_input_map[design_pair.exp_sample_id] = design_pair.ctr_sample_id
        ip_samples.append(design_pair.exp_sample_id)
        input_samples.append(design_pair.ctr_sample_id)

    paired_samples, single_samples = _classify_layouts(samples_info_dict)

    # Auto-detect organism and set genome default (RNAseq pattern)
    organisms = set()
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.organism:
            organisms.add(sample_info.organism)
    if len(organisms) == 1:
        try:
            organism = resolve_genome(next(iter(organisms)))
        except Exception as e:
            logger.error(f"Failed to resolve genome for organism {next(iter(organisms))}: {e}")
            raise ValueError(f"Failed to resolve genome for organism {next(iter(organisms))}: {e}")
        datajson["genome"]["default"] = organism
        datajson["Params"].setdefault("report", {}).setdefault("generate_report", {})["genome"] = organism
    elif len(organisms) > 1:
        logger.warning(f"Multiple organisms detected: {organisms}, using config default genome")
    datajson["Params"].setdefault("report", {}).setdefault("generate_report", {})["date"] = time.strftime("%Y-%m-%d", time.localtime())

    # Heatmap samples: from computeMatrix.samples, fallback to all ip_samples
    cm_params = datajson.get("Params", {}).get("computeMatrix", {})
    cm_heatmap_samples = cm_params.get("samples", None) or ip_samples
    # Heatmap output suffix based on regions mode
    cm_regions = cm_params.get("regions", "tss")
    if cm_regions == "tss":
        hm_suffix = "tss"
    elif cm_regions == "peaks":
        hm_suffix = "peaks"
    elif isinstance(cm_regions, dict) and "genes" in cm_regions:
        hm_suffix = "per_gene" if cm_regions.get("per_gene") else "genes"
    else:
        hm_suffix = "heatmap"
    # Gene names for per_gene mode
    cm_gene_names = cm_regions.get("genes", []) if isinstance(cm_regions, dict) else []

    if input_samples:
        for ip_sample in ip_samples:
            outfiles.append(f"{outdir}/QC/1_raw_fastqc/{ip_sample}/{ip_sample}.fastqc.raw.txt")
            outfiles.append(f"{outdir}/QC/2_trimmed_fastqc/{ip_sample}/{ip_sample}.fastqc.trimmed.txt")
            # Step 5: AddReadsGroup + MarkDuplicates (GATK4)
            outfiles.append(f"{outdir}/common/4_markdup_bam/{ip_sample}/{ip_sample}.sorted_markdup.bam")
            # Step 6: BigWig tracks
            outfiles.append(f"{outdir}/results/tracks/{ip_sample}/{ip_sample}.bigwig")
            # Step 7: MACS3 peak calling
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_broad_peaks.broadPeak")
            # Step 8: FRiP score
            outfiles.append(f"{outdir}/QC/3_frip_score/{ip_sample}/{ip_sample}.FRiP.txt")
            # Step 9: deeptools enrichment heatmap
            if ip_sample in cm_heatmap_samples:
                if hm_suffix == "per_gene" and cm_gene_names:
                    for gene_name in cm_gene_names:
                        _hm_gene_fmts = (datajson.get("Params", {}).get("plotHeatmap", {})
                                         .get("heatmap_gene", {}).get("formats") or ["png"])
                        if isinstance(_hm_gene_fmts, str):
                            _hm_gene_fmts = [_hm_gene_fmts]
                        outfiles += _figs(
                            f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{gene_name}_heatmap",
                            _hm_gene_fmts)
                else:
                    _hm_fmts = (datajson.get("Params", {}).get("plotHeatmap", {})
                                .get("heatmap", {}).get("formats") or ["png"])
                    if isinstance(_hm_fmts, str):
                        _hm_fmts = [_hm_fmts]
                    outfiles += _figs(
                        f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{hm_suffix}_heatmap",
                        _hm_fmts)
            # Step 10: HOMER annotation
            outfiles.append(f"{outdir}/results/annotation/{ip_sample}/{ip_sample}_peaks.annotatePeaks.txt")
            # Step 11: Peak-TE overlap
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_te_subfamily_overlap.tsv")
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_peak_centric_te.tsv")
        _macs3_cutoff_fmts = (datajson.get("Params", {}).get("macs3", {})
                              .get("macs3_cutoff_plot", {}).get("formats") or ["png"])
        if isinstance(_macs3_cutoff_fmts, str):
            _macs3_cutoff_fmts = [_macs3_cutoff_fmts]
        outfiles += _figs(f"{outdir}/results/peaks/cutoff_analysis",
                          _macs3_cutoff_fmts)
        _pte_fmts = (datajson.get("Params", {}).get("peak_te_overlap", {})
                     .get("peak_te_overlap_fig", {}).get("formats") or ["png"])
        if isinstance(_pte_fmts, str):
            _pte_fmts = [_pte_fmts]
        outfiles += _figs(f"{outdir}/results/te_overlap/te_subfamily_overlap_combined",
                          _pte_fmts)
        for input_sample in input_samples:
            outfiles.append(f"{outdir}/common/4_markdup_bam/{input_sample}/{input_sample}.sorted_markdup.bam")
            outfiles.append(f"{outdir}/results/tracks/{input_sample}/{input_sample}.bigwig")
    else:
        logger.warning("No Input samples found. MACS3 will run without control.")
        for ip_sample in ip_samples:
            sample_ip_input_map[ip_sample] = None
            outfiles.append(f"{outdir}/common/4_markdup_bam/{ip_sample}/{ip_sample}.sorted_markdup.bam")
            outfiles.append(f"{outdir}/results/tracks/{ip_sample}/{ip_sample}.bigwig")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_broad_peaks.broadPeak")
            outfiles.append(f"{outdir}/QC/3_frip_score/{ip_sample}/{ip_sample}.FRiP.txt")
            # Step 9: deeptools enrichment heatmap
            if ip_sample in cm_heatmap_samples:
                if hm_suffix == "per_gene" and cm_gene_names:
                    for gene_name in cm_gene_names:
                        _hm_gene_fmts = (datajson.get("Params", {}).get("plotHeatmap", {})
                                         .get("heatmap_gene", {}).get("formats") or ["png"])
                        if isinstance(_hm_gene_fmts, str):
                            _hm_gene_fmts = [_hm_gene_fmts]
                        outfiles += _figs(
                            f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{gene_name}_heatmap",
                            _hm_gene_fmts)
                else:
                    _hm_fmts = (datajson.get("Params", {}).get("plotHeatmap", {})
                                .get("heatmap", {}).get("formats") or ["png"])
                    if isinstance(_hm_fmts, str):
                        _hm_fmts = [_hm_fmts]
                    outfiles += _figs(
                        f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{hm_suffix}_heatmap",
                        _hm_fmts)
            # Step 10: HOMER annotation
            outfiles.append(f"{outdir}/results/annotation/{ip_sample}/{ip_sample}_peaks.annotatePeaks.txt")
            # Step 11: Peak-TE overlap
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_te_subfamily_overlap.tsv")
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_peak_centric_te.tsv")
    _pte_fmts = (datajson.get("Params", {}).get("peak_te_overlap", {})
                 .get("peak_te_overlap_fig", {}).get("formats") or ["png"])
    if isinstance(_pte_fmts, str):
        _pte_fmts = [_pte_fmts]
    outfiles += _figs(f"{outdir}/results/te_overlap/te_subfamily_overlap_combined",
                      _pte_fmts)
    outfiles.append(f"{outdir}/results/tracks/ucsc_track.txt")
    outfiles.append(f"{outdir}/results/tracks/igv_track.html")
    outfiles.append(f"{outdir}/PeakCalling_report.pptx")
    outfiles.append(f"{outdir}/PeakCalling_report.xlsx")
    outfiles.append(f"{outdir}/figures/manifest.tsv")
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    wf["paired_samples"] = paired_samples
    wf["single_samples"] = single_samples
    wf["samples"] = paired_samples + single_samples
    wf["ip_samples"] = ip_samples
    wf["input_samples"] = input_samples
    wf["sample_ip_input_map"] = sample_ip_input_map
    outfiles += _register_regionplot(datajson, samples_info_dict,
                                     ip_samples + input_samples,
                                     outdir, "PeakCalling")
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles

    return _write_raw_json(datajson, outdir)

def runRIPseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        design_pairs:List[DesignPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ):
    """Prepare input JSON for the RIPseq (RNA Immunoprecipitation sequencing) workflow.
    
    Workflow steps:
    1. FastQC (raw)
    2. Trimming (trim_galore)
    3. FastQC (trimmed)
    4. HISAT2 alignment
    5. AddReadsGroup + MarkDuplicates (GATK4)
    6. BigWig tracks (bamCoverage)
    7. MACS3 peak calling
    8. FRiP score
    9. deeptools enrichment heatmap
    10. HOMER peak annotation
    
    Supports IP vs Input peak analysis of RNA-protein binding (RIP-seq).
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)

    ip_samples = []
    input_samples = []
    sample_ip_input_map = {}
    outfiles = []
    for design_pair in design_pairs:
        sample_ip_input_map[design_pair.exp_sample_id] = design_pair.ctr_sample_id
        ip_samples.append(design_pair.exp_sample_id)
        input_samples.append(design_pair.ctr_sample_id)

    paired_samples, single_samples = _classify_layouts(samples_info_dict)

    # Auto-detect organism and set genome default (RNAseq pattern)
    organisms = set()
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.organism:
            organisms.add(sample_info.organism)
    if len(organisms) == 1:
        try:
            organism = resolve_genome(next(iter(organisms)))
        except Exception as e:
            logger.error(f"Failed to resolve genome for organism {next(iter(organisms))}: {e}")
            raise ValueError(f"Failed to resolve genome for organism {next(iter(organisms))}: {e}")
        datajson["genome"]["default"] = organism
        datajson["Params"].setdefault("report", {}).setdefault("generate_report", {})["genome"] = organism
    elif len(organisms) > 1:
        logger.warning(f"Multiple organisms detected: {organisms}, using config default genome")
    datajson["Params"].setdefault("report", {}).setdefault("generate_report", {})["date"] = time.strftime("%Y-%m-%d", time.localtime())

    # Heatmap samples: from computeMatrix.samples, fallback to all ip_samples
    cm_params = datajson.get("Params", {}).get("computeMatrix", {})
    cm_heatmap_samples = cm_params.get("samples", None) or ip_samples
    # Heatmap output suffix based on regions mode
    cm_regions = cm_params.get("regions", "tss")
    # Gene names for per_gene mode
    cm_gene_names = cm_regions.get("genes", []) if isinstance(cm_regions, dict) else []
    if cm_regions == "tss":
        hm_suffix = "tss"
    elif cm_regions == "peaks":
        hm_suffix = "peaks"
    elif isinstance(cm_regions, dict) and "genes" in cm_regions:
        # Per-gene heatmap files are only registered when the gene list is
        # non-empty (mirrors the deeptools module's HEATMAP_GENE_MODE gating);
        # otherwise the module emits the combined {sample}_genes_heatmap.png.
        hm_suffix = "per_gene" if (cm_regions.get("per_gene") and cm_gene_names) else "genes"
    else:
        hm_suffix = "heatmap"

    if input_samples:
        for ip_sample in ip_samples:
            outfiles.append(f"{outdir}/QC/1_raw_fastqc/{ip_sample}/{ip_sample}.fastqc.raw.txt")
            outfiles.append(f"{outdir}/QC/2_trimmed_fastqc/{ip_sample}/{ip_sample}.fastqc.trimmed.txt")
            # Step 5: AddReadsGroup + MarkDuplicates (GATK4)
            outfiles.append(f"{outdir}/common/4_markdup_bam/{ip_sample}/{ip_sample}.sorted_markdup.bam")
            # Step 6: BigWig tracks
            outfiles.append(f"{outdir}/results/tracks/{ip_sample}/{ip_sample}.bigwig")
            # Step 7: MACS3 peak calling
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_broad_peaks.broadPeak")
            # Step 8: FRiP score
            outfiles.append(f"{outdir}/QC/3_frip_score/{ip_sample}/{ip_sample}.FRiP.txt")
            # Step 9: deeptools enrichment heatmap
            if ip_sample in cm_heatmap_samples:
                if hm_suffix == "per_gene" and cm_gene_names:
                    for gene_name in cm_gene_names:
                        _hm_gene_fmts = (datajson.get("Params", {}).get("plotHeatmap", {})
                                         .get("heatmap_gene", {}).get("formats") or ["png"])
                        if isinstance(_hm_gene_fmts, str):
                            _hm_gene_fmts = [_hm_gene_fmts]
                        outfiles += _figs(
                            f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{gene_name}_heatmap",
                            _hm_gene_fmts)
                else:
                    _hm_fmts = (datajson.get("Params", {}).get("plotHeatmap", {})
                                .get("heatmap", {}).get("formats") or ["png"])
                    if isinstance(_hm_fmts, str):
                        _hm_fmts = [_hm_fmts]
                    outfiles += _figs(
                        f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{hm_suffix}_heatmap",
                        _hm_fmts)
            # Step 10: HOMER annotation
            outfiles.append(f"{outdir}/results/annotation/{ip_sample}/{ip_sample}_peaks.annotatePeaks.txt")
            # Step 11: Peak-TE overlap
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_te_subfamily_overlap.tsv")
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_peak_centric_te.tsv")
        _macs3_cutoff_fmts = (datajson.get("Params", {}).get("macs3", {})
                              .get("macs3_cutoff_plot", {}).get("formats") or ["png"])
        if isinstance(_macs3_cutoff_fmts, str):
            _macs3_cutoff_fmts = [_macs3_cutoff_fmts]
        outfiles += _figs(f"{outdir}/results/peaks/cutoff_analysis",
                          _macs3_cutoff_fmts)
        _pte_fmts = (datajson.get("Params", {}).get("peak_te_overlap", {})
                     .get("peak_te_overlap_fig", {}).get("formats") or ["png"])
        if isinstance(_pte_fmts, str):
            _pte_fmts = [_pte_fmts]
        outfiles += _figs(f"{outdir}/results/te_overlap/te_subfamily_overlap_combined",
                          _pte_fmts)
        for input_sample in input_samples:
            outfiles.append(f"{outdir}/common/4_markdup_bam/{input_sample}/{input_sample}.sorted_markdup.bam")
            outfiles.append(f"{outdir}/results/tracks/{input_sample}/{input_sample}.bigwig")
    else:
        logger.warning("No Input samples found. MACS3 will run without control.")
        for ip_sample in ip_samples:
            sample_ip_input_map[ip_sample] = None
            outfiles.append(f"{outdir}/common/4_markdup_bam/{ip_sample}/{ip_sample}.sorted_markdup.bam")
            outfiles.append(f"{outdir}/results/tracks/{ip_sample}/{ip_sample}.bigwig")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_broad_peaks.broadPeak")
            outfiles.append(f"{outdir}/QC/3_frip_score/{ip_sample}/{ip_sample}.FRiP.txt")
            # Step 9: deeptools enrichment heatmap
            if ip_sample in cm_heatmap_samples:
                if hm_suffix == "per_gene" and cm_gene_names:
                    for gene_name in cm_gene_names:
                        _hm_gene_fmts = (datajson.get("Params", {}).get("plotHeatmap", {})
                                         .get("heatmap_gene", {}).get("formats") or ["png"])
                        if isinstance(_hm_gene_fmts, str):
                            _hm_gene_fmts = [_hm_gene_fmts]
                        outfiles += _figs(
                            f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{gene_name}_heatmap",
                            _hm_gene_fmts)
                else:
                    _hm_fmts = (datajson.get("Params", {}).get("plotHeatmap", {})
                                .get("heatmap", {}).get("formats") or ["png"])
                    if isinstance(_hm_fmts, str):
                        _hm_fmts = [_hm_fmts]
                    outfiles += _figs(
                        f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{hm_suffix}_heatmap",
                        _hm_fmts)
            # Step 10: HOMER annotation
            outfiles.append(f"{outdir}/results/annotation/{ip_sample}/{ip_sample}_peaks.annotatePeaks.txt")
            # Step 11: Peak-TE overlap
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_te_subfamily_overlap.tsv")
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_peak_centric_te.tsv")
    _pte_fmts = (datajson.get("Params", {}).get("peak_te_overlap", {})
                 .get("peak_te_overlap_fig", {}).get("formats") or ["png"])
    if isinstance(_pte_fmts, str):
        _pte_fmts = [_pte_fmts]
    outfiles += _figs(f"{outdir}/results/te_overlap/te_subfamily_overlap_combined",
                      _pte_fmts)
    outfiles.append(f"{outdir}/results/tracks/ucsc_track.txt")
    outfiles.append(f"{outdir}/results/tracks/igv_track.html")
    # Report (subworkflow.md: pptx + xlsx pair, gated by Params.report.enabled)
    if _enabled(datajson, "report"):
        outfiles.append(f"{outdir}/RIPseq_report.pptx")
        outfiles.append(f"{outdir}/RIPseq_report.xlsx")
        outfiles.append(f"{outdir}/figures/manifest.tsv")
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    wf["paired_samples"] = paired_samples
    wf["single_samples"] = single_samples
    wf["samples"] = paired_samples + single_samples
    wf["ip_samples"] = ip_samples
    wf["input_samples"] = input_samples
    wf["sample_ip_input_map"] = sample_ip_input_map
    outfiles += _register_regionplot(datajson, samples_info_dict,
                                     ip_samples + input_samples,
                                     outdir, "RIPseq")
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles

    return _write_raw_json(datajson, outdir)

def runQuantMS(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ):
    """Prepare input JSON for QuantMS (quantitative proteomics) workflow.
    
    Workflow steps:
    1. Raw conversion (optional, when raw_manifest/raw inputs are provided)
    2. Decoy database generation
    3. Database search engines (Comet, MSGF+, Sage)
    4. PSM rescoring (Percolator)
    5. PSM FDR control
    6. Protein inference (EpiFany)
    7. Protein quantification (ProteomicsLFQ or ProteinQuantifier)
    8. Statistical analysis (MSstats)
    
    Supports TMT, LFQ, DIA, and raw-to-mzML entry points.
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files, "QuantMS")

    samples: List[str] = []
    outfiles: List[str] = []

    quantification_method = datajson["Params"]["workflow"]["quantification_method"]
    # Search engine follows the data type implied by the quantification method:
    # DIA spectra need a DIA-capable engine (sage), DDA (tmt/lfq) defaults to
    # comet. An explicit config engine wins, except comet-on-DIA which is
    # corrected (a DDA engine on wide-window spectra yields garbage).
    search_engine = datajson["Params"]["search_engine"].get("engine")
    _method_default = "sage" if quantification_method == "dia" else "comet"
    if not search_engine:
        search_engine = _method_default
    elif quantification_method == "dia" and search_engine == "comet":
        print("[runQuantMS] quantification_method=dia cannot use engine=comet "
              "(DDA engine on DIA spectra); switching to sage")
        search_engine = "sage"
    datajson["Params"]["search_engine"]["engine"] = search_engine

    # DIA with DIA-NN replaces the whole DDA identification chain
    # (search engine -> rescoring -> FDR -> protein inference -> quant).
    use_diann = (
        quantification_method == "dia"
        and (datajson.get("Params", {}).get("diann", {}) or {}).get("enabled", False)
    )
    if use_diann:
        print("[runQuantMS] quantification_method=dia with Params.diann.enabled: "
              "using DIA-NN, skipping the DDA identification chain")

    for sample_id in samples_info_dict:
        samples.append(sample_id)
        if use_diann:
            continue
        if search_engine== "comet":
            outfiles.append(f"{outdir}/common/3_search_engine/{sample_id}/{sample_id}_comet.idXML")
        elif search_engine == "msgf":
            outfiles.append(f"{outdir}/common/3_search_engine/{sample_id}/{sample_id}_msgf.idXML")
        elif search_engine == "sage":
            outfiles.append(f"{outdir}/common/3_search_engine/{sample_id}/{sample_id}_sage.idXML")
        elif search_engine == "msfragger":
            outfiles.append(f"{outdir}/common/3_search_engine/{sample_id}/{sample_id}.pepXML")
            outfiles.append(f"{outdir}/common/3_search_engine/{sample_id}/{sample_id}_msfragger.idXML")
        else:
            raise ValueError(f"Unknown search engine: {search_engine}")
        outfiles.append(f"{outdir}/common/4_psm_rescoring/{sample_id}/{sample_id}_scored.idXML")
        outfiles.append(f"{outdir}/common/5_psm_fdr/{sample_id}/{sample_id}_filtered.idXML")
        outfiles.append(f"{outdir}/common/6_protein_inference/{sample_id}/{sample_id}_protein.idXML")
    if use_diann:
        outfiles.append(f"{outdir}/common/3_diann/diann.report.tsv")
        outfiles.append(f"{outdir}/common/3_diann/diann.report.stats.tsv")
        outfiles.append(f"{outdir}/common/3_diann/diann_qc.tsv")
        outfiles.append(f"{outdir}/common/3_diann/diann_protein_matrix.tsv")
    elif quantification_method == "tmt":
        # Per-plex protein tables + long abundance table. The merged
        # consensusXML/mzTab triplet is only produced when an experimental
        # design is configured (FileMerger would collapse per-plex sample
        # identity otherwise).
        if datajson.get("Params", {}).get("workflow", {}).get("expdesign", ""):
            outfiles.append(f"{outdir}/common/7_quantification/tmt_quantification.mzTab")
            outfiles.append(f"{outdir}/common/7_quantification/tmt_quantification.consensusXML")
            outfiles.append(f"{outdir}/common/7_quantification/tmt_quantification.proteins.csv")
        else:
            outfiles.append(f"{outdir}/common/7_quantification/tmt_protein_abundances.csv")
    elif quantification_method == "lfq":
        outfiles.append(f"{outdir}/common/7_quantification/lfq_quantification.mzTab")
        outfiles.append(f"{outdir}/common/7_quantification/lfq_msstats.csv")
        if datajson.get("Params", {}).get("add_triqler_output", False):
            outfiles.append(f"{outdir}/common/7_quantification/lfq_quantification.triqler.tsv")
    elif quantification_method == "dia":
        outfiles.append(f"{outdir}/common/7_quantification/dia_quantification.mzTab")
        outfiles.append(f"{outdir}/common/7_quantification/dia_msstats.csv")
        if datajson.get("Params", {}).get("add_triqler_output", False):
            outfiles.append(f"{outdir}/common/7_quantification/dia_quantification.triqler.tsv")

    if not datajson.get("Params", {}).get("skip_post_msstats", False) and not use_diann:
        # TMT MSstats needs a plex->channel design; without one there is
        # nothing to convert (the long abundance table is the deliverable).
        if quantification_method == "tmt" and not datajson.get("Params", {}).get("workflow", {}).get("expdesign", ""):
            pass
        else:
            outfiles.append(f"{outdir}/common/8_msstats/msstats_input.csv")

    datajson["Params"]["workflow"]["samples"] = samples

    datajson["outfiles"] = outfiles

    return _write_raw_json(datajson, outdir)

def runRNAseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        group_pairs: List[CompareGroupPair],
        indir:str,
        outdir: str,
        raw_files: List[str],
    ):
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files, "RNAseq")
    outfiles = []
    genome_paired_samples = {}
    genome_single_samples = {}
    sample_groups = {}
    sample_contamination = {}
    datajson["Params"]["DESeq2"]["group_pairs"] = {}
    # read enabled flags from each module's own Params (default True)
    enable_te = _enabled(datajson, "TEtranscripts")
    enable_fusion = _enabled(datajson, "arriba")
    enable_variation = _enabled(datajson, "gatk_RNAseq")
    enable_deseq2 = _enabled(datajson, "DESeq2")
    enable_TEcount = _enabled(datajson, "TEcount")
    enable_TElocal = _enabled(datajson, "TElocal")
    enable_report = _enabled(datajson, "report")
    enable_StringTie = _enabled(datajson, "StringTie")
    for group_pair in group_pairs:
        datajson["Params"]["DESeq2"]["group_pairs"].setdefault(group_pair.organism, {})
        datajson["Params"]["DESeq2"]["group_pairs"][group_pair.organism].setdefault(f"{group_pair.ctr_group_name}_vs_{group_pair.exp_group_name}", {
            "control_group_name": group_pair.ctr_group_name,
            "experimental_group_name": group_pair.exp_group_name,
            "control_samples": group_pair.ctr_sample_ids,
            "experimental_samples": group_pair.exp_sample_ids
        })
        if enable_deseq2 and enable_TEcount:
            outfiles.append(f"{outdir}/diff_expression/{group_pair.organism}/{group_pair.ctr_group_name}_vs_{group_pair.exp_group_name}/TEcount/DESeq2.done")
        if enable_deseq2 and enable_TElocal:
            te_figures_dir = f"{outdir}/diff_expression/{group_pair.organism}/{group_pair.ctr_group_name}_vs_{group_pair.exp_group_name}/TElocal"
            outfiles.append(f"{te_figures_dir}/DESeq2.done")
            _manhattan_fmts = (datajson.get("Params", {}).get("DESeq2", {})
                               .get("DESeq2_TElocal_manhattan", {}).get("formats") or ["png"])
            if isinstance(_manhattan_fmts, str):
                _manhattan_fmts = [_manhattan_fmts]
            outfiles += _figs(
                f"{te_figures_dir}/manhattan/{group_pair.ctr_group_name}_vs_{group_pair.exp_group_name}.manhattan",
                _manhattan_fmts)
            # locus boxplot / subfamily bar only run when their per-contrast
            # family patterns are configured (Params.DESeq2.<rule>.families;
            # rule 键 = 真实规则名,与 modules/DESeq2 读点一致)。
            for te_fig_rule, te_fig_dir in (("DESeq2_TElocal_locus_boxplot", "locus_box"),
                                            ("DESeq2_TElocal_subfamily_bar", "subfamily_bar")):
                te_fig_families = (datajson.get("Params", {}).get("DESeq2", {}).get(te_fig_rule, {}) or {}).get("families", {}) or {}
                if te_fig_families.get(group_pair.organism, {}).get(f"{group_pair.ctr_group_name}_vs_{group_pair.exp_group_name}"):
                    te_fig_fmts = (datajson.get("Params", {}).get("DESeq2", {})
                                   .get(te_fig_rule, {}).get("formats") or ["png"])
                    if isinstance(te_fig_fmts, str):
                        te_fig_fmts = [te_fig_fmts]
                    outfiles += _figs(
                        f"{te_figures_dir}/{te_fig_dir}/{group_pair.ctr_group_name}_vs_{group_pair.exp_group_name}.{te_fig_dir}",
                        te_fig_fmts)
        sample_groups.setdefault(group_pair.ctr_group_name, []).extend(group_pair.ctr_sample_ids)
        sample_groups.setdefault(group_pair.exp_group_name, []).extend(group_pair.exp_sample_ids)

        if datajson.get("Params", {}).get("function", {}).get(group_pair.organism, {}).get("enabled", False):
            pair_dir = f"{outdir}/function/{group_pair.organism}/{group_pair.ctr_group_name}_vs_{group_pair.exp_group_name}"
            _fn_go_fmts = (datajson.get("Params", {}).get("function", {})
                           .get(group_pair.organism, {})
                           .get("function_go_kegg", {}).get("formats") or ["png"])
            if isinstance(_fn_go_fmts, str):
                _fn_go_fmts = [_fn_go_fmts]
            outfiles += _figs(f"{pair_dir}/go_back_to_back", _fn_go_fmts)
            outfiles += _figs(f"{pair_dir}/kegg_back_to_back", _fn_go_fmts)
            _fn_gsea_fmts = (datajson.get("Params", {}).get("function", {})
                             .get(group_pair.organism, {})
                             .get("function_gsea", {}).get("formats") or ["jpeg"])
            if isinstance(_fn_gsea_fmts, str):
                _fn_gsea_fmts = [_fn_gsea_fmts]
            outfiles += _figs(f"{pair_dir}/GSEA/TEcount_Gene_GSEA", _fn_gsea_fmts)
            outfiles.append(f"{pair_dir}/go_up.csv")
            outfiles.append(f"{pair_dir}/go_down.csv")
            outfiles.append(f"{pair_dir}/kegg_up.csv")
            outfiles.append(f"{pair_dir}/kegg_down.csv")
            outfiles.append(f"{pair_dir}/up_genes.txt")
            outfiles.append(f"{pair_dir}/down_genes.txt")
            outfiles.append(f"{pair_dir}/GSEA/TEcount_Gene_GSEA.csv")
    datajson["Params"]["StringTie"]["sample_groups"] = sample_groups
    Organisms = set()
    for sample_id, sample_info in samples_info_dict.items():
        Organisms.add(sample_info.organism)
        sample_contamination.setdefault(sample_id,{}).setdefault("host", sample_info.organism)
        if sample_info.contaminated_organism:
            sample_contamination.setdefault(sample_id,{}).setdefault("contaminant", sample_info.contaminated_organism)
        else:
            sample_contamination.setdefault(sample_id,{}).setdefault("contaminant", sample_info.organism)
        if enable_variation:
            outfiles.append(f"{outdir}/variation/germline_snv_indel_RNAseq/{sample_info.organism}/{sample_id}/{sample_id}.filtered.vcf.gz")
        if sample_info.layout == "PE":
            genome_paired_samples.setdefault(sample_info.organism, []).append(sample_id)
            if enable_StringTie:
                outfiles.append(f"{outdir}/transcripts/{sample_info.organism}/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt")
            if enable_fusion:
                outfiles.append(f"{outdir}/fusion/{sample_info.organism}/{sample_id}/{sample_id}_passed_fusions.tsv")
        elif sample_info.layout == "SE":
            genome_single_samples.setdefault(sample_info.organism, []).append(sample_id)
            if enable_StringTie:
                outfiles.append(f"{outdir}/transcripts/{sample_info.organism}/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt")
            if enable_fusion:
                outfiles.append(f"{outdir}/fusion/{sample_info.organism}/{sample_id}/{sample_id}_passed_fusions.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    datajson["Params"]["XenofilteR"]["sample_contamination"] = sample_contamination
    datajson["genome"]["default"] = ""
    datajson["Params"].setdefault("report", {}).setdefault("generate_report", {})["genome"] = ""
    for organism in Organisms:
        try:
            organism = resolve_genome(organism)
        except Exception as e:
            logger.error(f"Failed to resolve organism {organism}: {e}")
            raise ValueError(f"Failed to resolve organism {organism}: {e}")
        datajson["genome"]["default"] += organism
        datajson["Params"].setdefault("report", {}).setdefault("generate_report", {})["genome"] += organism
        # species 精确到 rule:镜像写入 function_go_kegg/function_gsva 子键
        if organism in ["Homo sapiens", "human", "GRCh38"]:
            _fn_species = "human"
        elif organism in ["Mus musculus", "mouse", "GRCm39"]:
            _fn_species = "mouse"
        else:
            raise ValueError(f"function module don't support {organism}, only support human or mouse(Homo sapiens or Mus musculus)")
        _fn_org = datajson["Params"].setdefault("function", {}).setdefault(organism, {})
        for _fn_rule in ("function_go_kegg", "function_gsva"):
            _fn_org.setdefault(_fn_rule, {})["species"] = _fn_species
    datajson["Params"].setdefault("report", {}).setdefault("generate_report", {})["date"] = time.strftime("%Y-%m-%d", time.localtime())
    for organism in Organisms:
        if enable_fusion:
            outfiles.append(f"{outdir}/fusion/{organism}/arriba_report/arriba_fusion_report.html")
        if enable_StringTie:
            outfiles.append(f"{outdir}/transcripts/{organism}/stringtie_merged.gtf")
            if len(group_pairs) > 0:
                for _st_fig in ("TE_chimeric_group_stacked",
                                "TE_chimeric_te_type_top",
                                "TE_chimeric_te_type_by_group"):
                    _st_fmts = (datajson.get("Params", {}).get("StringTie", {})
                                .get("formats") or ["png"])
                    if isinstance(_st_fmts, str):
                        _st_fmts = [_st_fmts]
                    outfiles += _figs(
                        f"{outdir}/transcripts/{organism}/TE_chimeric/{_st_fig}",
                        _st_fmts)
                outfiles.append(f"{outdir}/transcripts/{organism}/TE_chimeric/TE_chimeric_sample_summary.tsv")
                outfiles.append(f"{outdir}/transcripts/{organism}/TE_chimeric/TE_chimeric_group_summary.tsv")
                outfiles.append(f"{outdir}/transcripts/{organism}/TE_chimeric/TE_chimeric_te_type_counts.tsv")
        if enable_te or enable_deseq2:
            outfiles.append(f"{outdir}/counts/TEtranscripts/{organism}/{organism}_TEcount.tsv")
            outfiles.append(f"{outdir}/counts/TEtranscripts/{organism}/{organism}_TEcount_name.tsv")
        if enable_TElocal:
            outfiles.append(f"{outdir}/counts/TEtranscripts/{organism}/{organism}_TElocal.tsv")
            outfiles.append(f"{outdir}/counts/TEtranscripts/{organism}/{organism}_TElocal_name.tsv")

            # outfiles.append(f"{outdir}/counts/TEtranscripts/{organism}/{organism}_TElocal.tsv")
        # featureCounts quantification
        outfiles.append(f"{outdir}/counts/featureCounts/{organism}/{organism}_featureCounts.tsv")
        # GSVA enrichment analysis
        if datajson.get("Params", {}).get("function", {}).get(organism, {}).get("enabled", False):
            outfiles.append(f"{outdir}/function/{organism}/gsva/gsva_scores.tsv")
            _fn_gsva_fmts = (datajson.get("Params", {}).get("function", {})
                             .get(organism, {})
                             .get("function_gsva", {}).get("formats") or ["png"])
            if isinstance(_fn_gsva_fmts, str):
                _fn_gsva_fmts = [_fn_gsva_fmts]
            outfiles += _figs(f"{outdir}/function/{organism}/gsva/gsva_heatmap",
                              _fn_gsva_fmts)
        if enable_report:
            outfiles.append(f"{outdir}/results/{organism}/RNAseq_report.pptx")
            outfiles.append(f"{outdir}/results/{organism}/figures/manifest.tsv")
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles
    datajson["Params"]["workflow"]["genome_paired_samples"] = genome_paired_samples
    datajson["Params"]["workflow"]["genome_single_samples"] = genome_single_samples
    return _write_raw_json(datajson, outdir)

def runtRNAseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        meta: str,
        raw_files: List[str],
    ):
    """Prepare input JSON for tRNAseq (mim-tRNAseq) workflow.

    mim-tRNAseq is an all-in-one pipeline for tRNA sequencing analysis:
    tRNA clustering, GSNAP alignment, coverage analysis, modification
    quantification, CCA analysis, and DESeq2 differential expression.

    The pipeline processes all samples together via a sample data sheet.
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)
    # data_dir 精确到 rule:Params.mimseq.mimseq_tRNAtools.data_dir(config
    # 全量路径读取与写回,无 module 级共享回退)
    if not (datajson.get("Params", {}).get("mimseq", {})
            .get("mimseq_tRNAtools", {}) or {}).get("data_dir"):
        datajson.setdefault("Params", {}).setdefault("mimseq", {}).setdefault(
            "mimseq_tRNAtools", {})["data_dir"] = os.path.join(
            root_dir, "modules", "mimseq", "mimseq", "data")
    datajson.setdefault("Params", {}).setdefault("workflow", {})["meta"] = meta
    samples = []
    for sample_id in samples_info_dict:
        samples.append(sample_id)

    outfiles = [f"{outdir}/mimseq/mimseq.done"]

    datajson["Params"]["workflow"]["samples"] = samples
    datajson["outfiles"] = outfiles

    return _write_raw_json(datajson, outdir)

def runncRNAseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        design_pairs:List[DesignPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ):
    """Prepare input JSON for ncRNAseq (small/non-coding RNA-seq) workflow.

    Pipeline: jla-demultiplexer -> trim_galore -> subsample -> STAR (star / star_3pass / star_3pass_gene) -> featureCounts + Tailer.
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files, "ncRNAseq")
    sample_ip_input_map = {}
    ip_samples = []
    input_samples = []
    outfiles = []
    paired_samples = []
    single_samples = []
    layouts = set()
    organisms = set()
    for design_pair in design_pairs:
        sample_ip_input_map[design_pair.exp_sample_id] = design_pair.ctr_sample_id
        ip_samples.append(design_pair.exp_sample_id)
        input_samples.append(design_pair.ctr_sample_id)
    aligner = datajson.get("Params", {}).get("workflow",{}).get("aligner") or "star"
    if aligner == "star_3pass":
        bam_subdir = "common/3_raw_bam"
    elif aligner == "star_3pass_gene":
        bam_subdir = "common/4_per_gene_bam"
    else:
        bam_subdir = "common/3_raw_bam"
    for sample_id, sample_info in samples_info_dict.items():
        organisms.add(sample_info.organism)
        outfiles.append(f"{outdir}/{bam_subdir}/{sample_id}/{sample_id}.bam")
        if aligner in ["star_3pass", "star_3pass_gene"]:
            outfiles.append(f"{outdir}/{bam_subdir}/{sample_id}/{sample_id}_tail.csv")
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            layouts.add("PE")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            layouts.add("SE")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    if len(organisms) != 1:
        raise ValueError(f"meta don't support multiple organsim temporarily, please check your meta file, found: {organisms}")
    for organism in organisms:
        try:
            datajson["Params"]["track"]["igv"]["default"] = resolve_genome(organism)
            datajson["genome"]["default"] = resolve_genome(organism)
        except Exception as e:
            logger.error(f"Failed to resolve genome for organism {organism}: {e}")
            raise ValueError(f"Failed to resolve genome for organism {organism}: {e}")
    for layout in layouts:
        if layout == "PE":
            outfiles.append(f"{outdir}/counts/all_paired_featureCounts.tsv")
        elif layout == "SE":
            outfiles.append(f"{outdir}/counts/all_single_featureCounts.tsv")
    all_samples = paired_samples + single_samples
    outfiles.append(f"{outdir}/tracks/igv_track.html")
    outfiles.append(f"{outdir}/tracks/ucsc_track.txt")
    outfiles.append(f"{outdir}/ncRNAseq_report.pptx")
    outfiles.append(f"{outdir}/ncRNAseq_report_files.xlsx")
    outfiles.append(f"{outdir}/figures/manifest.tsv")
    outfiles += _register_regionplot(datajson, samples_info_dict,
                                     all_samples, outdir, "ncRNAseq")
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles
    datajson["Params"]["workflow"]["samples"] = all_samples
    datajson["Params"]["workflow"]["raw_files"] = raw_files
    datajson["Params"]["workflow"]["paired_samples"] = paired_samples
    datajson["Params"]["workflow"]["single_samples"] = single_samples
    datajson["Params"]["workflow"]["ip_samples"] = ip_samples
    datajson["Params"]["workflow"]["input_samples"] = input_samples
    datajson["Params"]["workflow"]["sample_ip_input_map"] = sample_ip_input_map
    return _write_raw_json(datajson, outdir)

def runscRNAseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
        cellranger_input_dict: Dict[str, CellrangerInput]
    ):
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files, "scRNAseq")
    datajson["Params"]["workflow"]["cellranger_input_dict"] = {k: v.__dict__ for k, v in cellranger_input_dict.items()}
    outfiles = []
    organisms = set()
    # CLI --extra-args counters=/aligner= land at config top level and win
    # over Params.workflow; normalize into Params.workflow for raw.json.
    _wf = datajson["Params"]["workflow"]
    counters = datajson.pop("counters", None) or _wf.get("counters") or ["scTE", "cellranger"]
    _wf["counters"] = counters
    paired_samples = []
    single_samples = []
    cellranger_dirs = []
    first_sample_id = next(iter(samples_info_dict))
    first_sample_info = samples_info_dict[first_sample_id]
    if first_sample_info.cellranger_dir and os.path.exists(first_sample_info.cellranger_dir):
        datajson["Params"]["workflow"]["input_mode"] = "cellranger_dir"
    elif first_sample_info.bam and os.path.exists(first_sample_info.bam):
        datajson["Params"]["workflow"]["input_mode"] = "bam"
    else:
        datajson["Params"]["workflow"]["input_mode"] = "fastq"
    for sample_id, sample_info in samples_info_dict.items():
        organisms.add(sample_info.organism)
        if sample_info.cellranger_dir and os.path.exists(sample_info.cellranger_dir):
            cellranger_dirs.append(sample_info.cellranger_dir)
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            if "scTE" in counters:
                outfiles.append(f"{outdir}/common/3_raw_h5ad/{sample_id}/{sample_id}_scTE.h5ad")
                outfiles.append(f"{outdir}/common/4_qc_h5ad/{sample_id}/{sample_id}_scTE_qc.h5ad")
            if "cellranger" in counters:
                outfiles.append(f"{outdir}/common/3_raw_h5ad/{sample_id}/{sample_id}_cellranger.h5ad")
                outfiles.append(f"{outdir}/common/4_qc_h5ad/{sample_id}/{sample_id}_cellranger_qc.h5ad")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            if "scTE" in counters:
                outfiles.append(f"{outdir}/common/3_raw_h5ad/{sample_id}/{sample_id}_scTE.h5ad")
                outfiles.append(f"{outdir}/common/4_qc_h5ad/{sample_id}/{sample_id}_scTE_qc.h5ad")
            if "cellranger" in counters:
                outfiles.append(f"{outdir}/common/3_raw_h5ad/{sample_id}/{sample_id}_cellranger.h5ad")
                outfiles.append(f"{outdir}/common/4_qc_h5ad/{sample_id}/{sample_id}_cellranger_qc.h5ad")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    datajson["Params"]["workflow"]["paired_samples"] = paired_samples
    datajson["Params"]["workflow"]["single_samples"] = single_samples
    datajson["Params"]["workflow"]["cellranger_dirs"] = cellranger_dirs

    # scCODA compare sharing: Params.sccoda.<counter>.sccoda_run.compare
    # defaults to the pseudobulk comparison design (identical
    # tissue:ref:treat[=label] entry shape) — one
    # --Params.scanpy.<counter>.scanpy_pseudobulk.compare serves both
    # analyses. A hand-filled sccoda compare wins over the shared one.
    for counter in counters:
        sccoda_run = (datajson["Params"].setdefault("sccoda", {})
                      .setdefault(counter, {}).setdefault("sccoda_run", {}))
        if not sccoda_run.get("compare"):
            pb = ((datajson["Params"].get("scanpy", {}).get(counter, {}) or {})
                  .get("scanpy_pseudobulk", {}) or {})
            shared = pb.get("compare") or []
            if shared:
                sccoda_run["compare"] = list(shared)
                logger.info(
                    f"[sccoda] Params.sccoda.{counter}.sccoda_run.compare is "
                    f"empty — sharing the pseudobulk design ({len(shared)} entries)")

    if len(organisms) != 1:
        raise ValueError(f"meta don't support multiple organsim temporarily, please check your meta file, found: {organisms}")
    try:
        organism = resolve_genome(next(iter(organisms)))
    except Exception as e:
        logger.error(f"Failed to resolve genome for organism {next(iter(organisms))}: {e}")
        raise ValueError(f"Failed to resolve genome for organism {next(iter(organisms))}: {e}")

    datajson["genome"]["default"] = organism
    datajson["Params"]["scTE"]["genome"] = organism
    datajson["Params"]["cellranger"]["mkref"]["genome_name"] = organism
    datajson["Params"]["cellranger"]["mkref"]["version"] = time.strftime("%Y-%m-%d", time.localtime())
    counters = _wf["counters"]
    aligner = datajson.pop("aligner", None) or _wf.get("aligner") or "cellranger"
    _wf["aligner"] = aligner
    if "scTE" in counters:
        datajson["Params"]["scanpy"]["scTE"]["scanpy_qc"]["min_genes"] = 200
        datajson["Params"]["scanpy"]["scTE"]["scanpy_qc"]["max_genes"] = 6000
        datajson["Params"]["scanpy"]["scTE"]["scanpy_qc"]["scrublet"] = True
        if aligner == "star":
            datajson["Params"]["scTE"]["cb_tag"] = "CR"
            datajson["Params"]["scTE"]["umi_tag"] = "UR"
        elif aligner == "cellranger":
            datajson["Params"]["scTE"]["cb_tag"] = "CB"
            datajson["Params"]["scTE"]["umi_tag"] = "UB"
        else:
            raise ValueError(f"Unsupported aligner for scTE counter: {aligner}, must be star or cellranger")
    elif "cellranger" in counters:
        # cellranger: keep defaults (min_genes=200, max_genes=6000, scrublet=True)
        datajson["Params"]["scanpy"]["cellranger"]["scanpy_qc"]["min_genes"] = 200
        datajson["Params"]["scanpy"]["cellranger"]["scanpy_qc"]["max_genes"] = 6000
        datajson["Params"]["scanpy"]["cellranger"]["scanpy_qc"]["scrublet"] = True
    else:
        raise ValueError(f"Unsupported counter type: {counters}, must be scTE or cellranger")
    # Inject LLM config from environment variables (config takes priority over env)
    env_llm_method = os.environ.get("LLM_METHOD", None)
    env_llm_model = os.environ.get("LLM_MODEL", None)
    env_llm_api_key = os.environ.get("LLM_API_KEY", None)
    env_llm_base_url = os.environ.get("LLM_BASE_URL", None)
    # Build tissue_samples for scanpy downstream
    tissue_samples = {}
    for sid in paired_samples + single_samples:
        tissue = getattr(samples_info_dict.get(sid), "tissue", None) or "unknown"
        tissue_samples.setdefault(tissue, []).append(sid)
    for tissue in tissue_samples.keys():
        for counter in counters:
            if env_llm_method and env_llm_model and env_llm_api_key and env_llm_base_url:
                outfiles.append(f"{outdir}/common/5_combine_h5ad/{tissue}/{tissue}_{counter}_advanced.h5ad")
            else:
                outfiles.append(f"{outdir}/common/5_combine_h5ad/{tissue}/{tissue}_{counter}_merged.h5ad")
            # Pseudo-bulk DEG (PyDESeq2) — gated per counter by
            # Params.scanpy.<counter>.scanpy_pseudobulk.enabled; consumes the
            # scanpy_auto h5ad, so registering the report target pulls the
            # qc→merge→auto chain into the DAG
            pb_cfg = ((datajson["Params"].get("scanpy", {}).get(counter, {}) or {})
                      .get("scanpy_pseudobulk", {}) or {})
            # Register the report only for tissues owning at least one
            # compare entry: entries are tissue-scoped
            # "tissue:ref:treat[=label]" strings and the per-tissue rule
            # filters by prefix, so a tissue with no entries cannot run.
            pb_compare = pb_cfg.get("compare") or []
            if pb_cfg.get("enabled") and any(
                    str(c).split(":", 1)[0] == tissue for c in pb_compare):
                outfiles.append(f"{outdir}/diff/pseudobulk/{tissue}/{counter}/pseudobulk_report.html")
            # scCODA differential abundance (modules/sccoda) — gated per
            # counter by Params.sccoda.<counter>.sccoda_run.enabled; consumes
            # the scanpy_auto h5ad (obs: sample_id + cell_type), so
            # registering the report target pulls the qc->merge->auto chain
            # into the DAG. Same tissue-prefix gate expression as the rule's
            # compare filter.
            sc_cfg = ((datajson["Params"].get("sccoda", {}) or {}).get(counter, {})
                      or {}).get("sccoda_run", {}) or {}
            sc_compare = sc_cfg.get("compare") or []
            if sc_cfg.get("enabled") and any(
                    str(c).split(":", 1)[0] == tissue for c in sc_compare):
                outfiles.append(f"{outdir}/diff/abundance/{tissue}/{counter}/sccoda_report.html")
            # Manual cluster-annotation overrides (scanpy_reannotate): an
            # in-place label fix over the auto h5ad (auto.h5ad.bak keeps
            # one generation). Register the marker file only for tissues
            # with an override JSON configured (producer rule
            # scanpy_reannotate's own cluster_annotations key).
            ca = ((datajson["Params"].get("scanpy", {}).get(counter, {}) or {})
                  .get("scanpy_reannotate", {}) or {}).get("cluster_annotations", {}) or {}
            if tissue in ca:
                outfiles.append(f"{outdir}/common/5_combine_h5ad/{tissue}/{tissue}_{counter}_auto.h5ad.reannotated")
    datajson["Params"]["scanpy"]["tissue_samples"] = tissue_samples
    datajson["Params"]["sccoda"]["tissue_samples"] = tissue_samples
    # GO/KEGG enrichment over the pseudobulk tree (modules/function/
    # pseudobulk_go): one sentinel target pulls the per-unit enrichment jobs
    # behind the runtime discovery checkpoint. Gated by
    # Params.function.pseudobulk_go.enabled.
    pbgo_cfg = (datajson["Params"].get("function", {}) or {}).get("pseudobulk_go", {}) or {}
    if pbgo_cfg.get("enabled"):
        outfiles.append(f"{outdir}/function/go_kegg/enrichment.done")
    # TE-dominant cluster identification (modules/tedominant): one
    # group-level run over all tissues, gated by
    # Params.tedominant.enabled; consumes the combine-stage h5ad pair
    # (registering the table + done marker pulls the qc->merge->auto
    # chain for every tissue into the DAG).
    tedom_cfg = datajson["Params"].get("tedominant", {}) or {}
    if tedom_cfg.get("enabled"):
        outfiles.append(f"{outdir}/diff/figures_te/te_dominant_clusters.csv")
        outfiles.append(f"{outdir}/diff/figures_te/markers_dotplot.png")
        outfiles.append(f"{outdir}/diff/figures_te/te_dominant.done")
        for tissue in tissue_samples.keys():
            outfiles.append(f"{outdir}/diff/figures_te/cluster_burden_{tissue}.png")
            outfiles.append(f"{outdir}/diff/figures_te/umap_te_{tissue}.png")
    # TE profile (modules/teprofile): three-layer TE expression
    # characterization over the same h5ad pairs; gated by
    # Params.teprofile.enabled.
    tep_cfg = datajson["Params"].get("teprofile", {}) or {}
    if tep_cfg.get("enabled"):
        outfiles.append(f"{outdir}/diff/te_profile/te_profile_baseline.csv")
        outfiles.append(f"{outdir}/diff/te_profile/te_profile_clusters.csv")
        outfiles.append(f"{outdir}/diff/te_profile/te_profile.done")
        for tissue in tissue_samples.keys():
            outfiles.append(f"{outdir}/diff/te_profile/baseline_te_fraction_{tissue}.png")
            outfiles.append(f"{outdir}/diff/te_profile/cluster_deviation_{tissue}.png")
    # outfiles.append(f"{outdir}/scRNAseq_report.pptx")
    datajson["outfiles"] = outfiles
    return _write_raw_json(datajson, outdir)



def _contrasts_from_group_pairs(
        group_pairs: List[CompareGroupPair],
) -> Dict[str, Dict[str, List[str]]]:
    """Meta-derived contrast map: name -> {group -> [samples]}.

    Shared by every group-level Fiber-seq analysis (mods, diff, ...):
    one source of grouping (the meta design/group columns that drive
    group_pairs), no per-module design files.
    """
    contrasts: Dict[str, Dict[str, List[str]]] = {}
    for gp in group_pairs:
        contrasts[f"{gp.ctr_group_name}_vs_{gp.exp_group_name}"] = {
            gp.ctr_group_name: list(gp.ctr_sample_ids),
            gp.exp_group_name: list(gp.exp_sample_ids),
        }
    return contrasts


def runFiberseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        group_pairs: List[CompareGroupPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ) -> str:
    """Prepare input JSON for Fiber-seq workflow.

    Pipeline: ft predict-m6a → ft fire → pbmm2 align → ft extract/call-peaks/qc
    Reference: Stergachis et al., 2020, Science (DOI: 10.1126/science.aaz1646).
    Guide: https://fiberseq.github.io/

    `Params.sv.enabled` gates the SV comparison analysis (modules/sv).
    `Params.sv.comparisons` follows the DESeq2 group_pairs design
    (name -> group names + sample name lists); it is filled from the meta
    group pairs when not hand-authored in the config.
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files, "Fiberseq")

    outfiles = []
    samples = []
    genome_versions = set()
    default_ref = _resolve_genome_ref(datajson, samples_info_dict, "Fiberseq")
    ### switch 
    enable_annotate = _enabled(datajson, "annotate")
    enable_motif = _enabled(datajson, "motif")
    enable_footprint = _enabled(datajson, "footprint")
    enable_nucpos = _enabled(datajson, "nucpos")
    enable_track = _enabled(datajson, "track")
    enable_coactuation = _enabled(datajson, "coactuation")
    enable_censat = _enabled(datajson, "censat")
    enable_quant = _enabled(datajson, "quant")
    enable_diff = _enabled(datajson, "diff")
    enable_haplotype = _enabled(datajson, "haplotype")
    enable_lsdc = _enabled(datajson, "ldsc")
    enable_function = _enabled(datajson, "function")
    enable_sv = _enabled(datajson, "sv")
    enable_telomere = _enabled(datajson, "telomere")
    enable_telocen_norm = _enabled(datajson, "telocen_norm")
    enable_epi_drift = _enabled(datajson, "epi_drift")
    want_motif = enable_motif or enable_footprint
    want_quant = enable_quant or enable_diff or enable_haplotype
    want_eccdna = _enabled(datajson, "decoil") and default_ref.get("gtf")

    ### collect data
    for sample_id, sample_info in samples_info_dict.items():
        genome_versions.add(sample_info.organism)
        samples.append(sample_id)
        # Step 1: prepare fiberseq BAM (predict-m6a + fire)
        outfiles.append(f"{outdir}/common/2_fiberseq_bam/{sample_id}/{sample_id}.fiberseq.fire.bam")
        # Step 2: aligned BAM (pbmm2)
        outfiles.append(f"{outdir}/common/3_align_bam/{sample_id}/{sample_id}.fiberseq.fire.sorted.bam")
        outfiles.append(f"{outdir}/common/3_align_bam/{sample_id}/{sample_id}.fiberseq.fire.sorted.bai")
        # Step 3: raw results
        outfiles.append(f"{outdir}/common/4_elements/{sample_id}/{sample_id}.m6a.bed.gz")
        outfiles.append(f"{outdir}/common/4_elements/{sample_id}/{sample_id}.nuc.bed.gz")
        outfiles.append(f"{outdir}/common/4_elements/{sample_id}/{sample_id}.msp.bed.gz")
        outfiles.append(f"{outdir}/common/4_elements/{sample_id}/{sample_id}.cpg.bed.gz")
        outfiles.append(f"{outdir}/common/4_elements/{sample_id}/{sample_id}.fire_peaks.bed")
        outfiles.append(f"{outdir}/common/4_elements/{sample_id}/{sample_id}.fire_fdr_table.tsv")
        outfiles.append(f"{outdir}/common/4_elements/{sample_id}/{sample_id}.fire.bed.gz")
        outfiles.append(f"{outdir}/common/4_elements/{sample_id}/{sample_id}.qc.tsv")

        # Step 4: downstream (gated by Params.*.enabled)
        if enable_annotate:
            outfiles.append(f"{outdir}/fiberseq_analysis/annotate/{sample_id}/{sample_id}.peaks_annotated.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/annotate/{sample_id}/{sample_id}.peaks_summary.tsv")
        if want_motif:
            outfiles.append(f"{outdir}/fiberseq_analysis/motif/{sample_id}/{sample_id}.motif_sites.bed")    
        if enable_motif:
            outfiles.append(f"{outdir}/fiberseq_analysis/motif/{sample_id}/{sample_id}.motif_profile.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/motif/{sample_id}/{sample_id}.motif_footprint_scores.tsv")
            _motif_fmts = (datajson.get("Params", {}).get("motif", {})
                           .get("fiberseq_motif_profile", {}).get("formats") or ["pdf"])
            if isinstance(_motif_fmts, str):
                _motif_fmts = [_motif_fmts]
            outfiles += _figs(
                f"{outdir}/fiberseq_analysis/motif/{sample_id}/{sample_id}.motif_profile",
                _motif_fmts)
        if enable_footprint:
            outfiles.append(f"{outdir}/fiberseq_analysis/footprint/{sample_id}/{sample_id}.footprint_categories.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/footprint/{sample_id}/{sample_id}.footprint_profiles.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/footprint/{sample_id}/{sample_id}.footprint_scores.tsv")
            _footprint_fmts = (datajson.get("Params", {}).get("footprint", {})
                               .get("fiberseq_footprint_run", {}).get("formats") or ["pdf"])
            if isinstance(_footprint_fmts, str):
                _footprint_fmts = [_footprint_fmts]
            outfiles += _figs(
                f"{outdir}/fiberseq_analysis/footprint/{sample_id}/{sample_id}.footprint_profile",
                _footprint_fmts)
        if enable_nucpos:
            outfiles.append(f"{outdir}/fiberseq_analysis/nucpos/{sample_id}/{sample_id}.nuc_offset_per_anchor.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/nucpos/{sample_id}/{sample_id}.nuc_offset_summary.tsv")
            _nucpos_fmts = (datajson.get("Params", {}).get("nucpos", {})
                            .get("fiberseq_nucpos_run", {}).get("formats") or ["pdf"])
            if isinstance(_nucpos_fmts, str):
                _nucpos_fmts = [_nucpos_fmts]
            outfiles += _figs(
                f"{outdir}/fiberseq_analysis/nucpos/{sample_id}/{sample_id}.nucpos",
                _nucpos_fmts)
        if enable_track:
            outfiles.append(f"{outdir}/fiberseq_analysis/track/{sample_id}/{sample_id}.accessibility.bw")
            outfiles.append(f"{outdir}/fiberseq_analysis/track/{sample_id}/{sample_id}.m6a_density.bw")
            outfiles.append(f"{outdir}/fiberseq_analysis/track/{sample_id}/trackDb.txt")
        if want_quant:
            outfiles.append(f"{outdir}/fiberseq_analysis/quant/{sample_id}/{sample_id}.region_quant.tsv")
        if enable_coactuation:
            outfiles.append(f"{outdir}/fiberseq_analysis/coactuation/{sample_id}/{sample_id}.coactuation.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/coactuation/{sample_id}/{sample_id}.coactuation_summary.tsv")
        if enable_censat:
            outfiles.append(f"{outdir}/fiberseq_analysis/censat/{sample_id}/{sample_id}.censat_regions.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/censat/{sample_id}/{sample_id}.censat_class_summary.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/censat/{sample_id}/{sample_id}.censat_profile.tsv")
            _censat_fmts = (datajson.get("Params", {}).get("censat", {})
                            .get("fiberseq_censat_run", {}).get("formats") or ["pdf"])
            if isinstance(_censat_fmts, str):
                _censat_fmts = [_censat_fmts]
            outfiles += _figs(
                f"{outdir}/fiberseq_analysis/censat/{sample_id}/{sample_id}.censat_profile",
                _censat_fmts)

        # Telomere length analysis (modules/telomere); method-gated.
        # NOTE: read_density on the aligned fire.bam is alignment-biased (module
        # README §3.1/§3.4) — ideally run on the raw HiFi BAM for unbiased ATL.
        if enable_telomere:
            _TEL_METHODS = ("telogator2", "assembly_scan", "read_density", "tidk")
            _tel_methods = (datajson.get("Params", {}).get("telomere", {}).get("method")
                            or list(_TEL_METHODS))
            if isinstance(_tel_methods, str):
                _tel_methods = [_tel_methods]
            _bad_tel = [m for m in _tel_methods if m not in _TEL_METHODS]
            if _bad_tel:
                raise ValueError(
                    f"Fiberseq: Params.telomere.method has unknown method(s) {_bad_tel}; "
                    f"valid methods: {list(_TEL_METHODS)}.")
            _tel_fmts = (datajson.get("Params", {}).get("telomere", {})
                         .get("formats") or ["png"])
            if isinstance(_tel_fmts, str):
                _tel_fmts = [_tel_fmts]
            _tb = f"{outdir}/telomere/length"
            if "telogator2" in _tel_methods:
                outfiles.append(f"{_tb}/telogator2/{sample_id}/tlens_by_allele.tsv")
                outfiles.append(f"{_tb}/telogator2/{sample_id}/all_final_alleles.png")
                outfiles.append(f"{_tb}/telogator2/{sample_id}/violin_atl.png")
            if "assembly_scan" in _tel_methods:
                outfiles.append(f"{_tb}/assembly_scan/{sample_id}/{sample_id}_assembly_telomere_stats.txt")
                outfiles += _figs(
                    f"{_tb}/assembly_scan/{sample_id}/{sample_id}_assembly_telomere_length", _tel_fmts)
            if "read_density" in _tel_methods:
                outfiles.append(f"{_tb}/read_density/{sample_id}/{sample_id}_read_telomere_stats.txt")
                outfiles += _figs(
                    f"{_tb}/read_density/{sample_id}/{sample_id}_read_telomere_content", _tel_fmts)
            if "tidk" in _tel_methods:
                outfiles.append(f"{_tb}/tidk/{sample_id}/{sample_id}_tidk_telomeres.tsv")
                outfiles += _figs(f"{_tb}/tidk/{sample_id}/{sample_id}_tidk_telomeres", _tel_fmts)
            if "read_density" in _tel_methods:
                for _f in _tel_fmts:
                    outfiles.append(f"{_tb}/read_density/combined.{_f}")
            if len(_tel_methods) >= 4:
                for _f in _tel_fmts:
                    outfiles.append(f"{_tb}/combined_methods.{_f}")
            _dq_bed = ((datajson.get("Params", {}).get("telomere", {})
                        .get("depth_qc", {}) or {}).get("telomere_bed"))
            if _dq_bed:
                outfiles.append(f"{_tb}/depth_qc/{sample_id}/{sample_id}_telomere_region_depth_stats.txt")
            outfiles.append(f"{outdir}/telomere/composition/telomere_summary.tsv")
            outfiles.append(f"{outdir}/telomere/composition/telomere_composition.tsv")
            for _f in _tel_fmts:
                outfiles.append(f"{outdir}/telomere/composition/telomere_length_violin.{_f}")
            if "assembly_scan" in _tel_methods:
                outfiles.append(f"{outdir}/telomere/regions/{sample_id}/{sample_id}.telomere.bed")

        # Step 5: haplotype chain
        if enable_haplotype:
            outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.vcf.gz")
            outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.vcf.gz.csi")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.svsig.gz")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.vcf.gz")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.vcf.gz.csi")
            # hiphase module writes to {outdir}/common/5_phased_bam (see
            # subworkflow/Fiberseq.smk hiphase_config), NOT {var}/phased.
            outfiles.append(f"{outdir}/common/5_phased_bam/{sample_id}/{sample_id}.phased.bam")
            phased_dir = f"{outdir}/common/5_phased_bam"
            for hap in ("hap1", "hap2"):
                outfiles.append(f"{phased_dir}/{sample_id}/{sample_id}.{hap}.bam")
                outfiles.append(f"{phased_dir}/{sample_id}/{sample_id}.{hap}.bam.bai")
                outfiles.append(f"{phased_dir}/{sample_id}/{sample_id}.{hap}.m6a.bed.gz")
                outfiles.append(f"{phased_dir}/{sample_id}/{sample_id}.{hap}.nuc.bed.gz")
                outfiles.append(f"{phased_dir}/{sample_id}/{sample_id}.{hap}.msp.bed.gz")
                outfiles.append(f"{phased_dir}/{sample_id}/{sample_id}.{hap}.cpg.bed.gz")
                outfiles.append(f"{outdir}/fiberseq_analysis/quant/{sample_id}/{sample_id}.{hap}.region_quant.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/haplotype/{sample_id}/{sample_id}.haplotype_diff.tsv")
            outfiles.append(f"{outdir}/fiberseq_analysis/haplotype/{sample_id}/{sample_id}.haplotype_summary.tsv")

        # Step 6: eccDNA reconstruction (requires genome GTF) 
        if want_eccdna:
            outfiles.append(f"{outdir}/eccdna/{sample_id}/decoil/reconstruct.bed")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/decoil/reconstruct.ecDNA.bed")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/decoil/summary.txt")
            # Visualization: per-sample stats + interactive IGV HTML
            # Same key and fallback as modules/decoil/decoil.smk ("formats",
            # plural, or ["png"]) so registered targets match rule outputs.
            # Normalize + validate: a str default or dirty JSON must not leak
            # through as a non-list ("for fmt in 'png'" iterates characters).
            formats = datajson.get("Params", {}).get("decoil", {}).get("structure", {}).get("formats") or ["png"]
            if isinstance(formats, str):
                formats = [formats]
            if not isinstance(formats, list) or not all(isinstance(f, str) for f in formats):
                raise ValueError("Params.decoil.structure.formats must be a list of strings")
            for fmt in formats:
                outfiles.append(f"{outdir}/eccdna/{sample_id}/visualization/decoil_stats.{fmt}")
                outfiles.append(f"{outdir}/eccdna/{sample_id}/visualization/decoil_structure.{fmt}")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/visualization/decoil_structure_genes.tsv")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/visualization/igv_ecdna.html")
            outfiles.append(f"{outdir}/eccdna/{sample_id}/visualization/ucsc_track.txt")
            # Cross-sample comparison (aggregated): the decoil module writes
            # these under {outdir}/compare/ (see modules/decoil/decoil.smk
            # rule decoil_compare), NOT {outdir}/visualization/.
            _decoil_cmp_fmts = (datajson.get("Params", {}).get("decoil", {})
                                .get("compare", {}).get("formats") or ["png"])
            if isinstance(_decoil_cmp_fmts, str):
                _decoil_cmp_fmts = [_decoil_cmp_fmts]
            outfiles += _figs(f"{outdir}/eccdna/compare/decoil_compare",
                              _decoil_cmp_fmts)
            outfiles.append(f"{outdir}/eccdna/compare/decoil_compare_shared_unique.tsv")
        if enable_sv:
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.svsig.gz")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.vcf.gz")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.vcf.gz.csi")

    # Comparison figure modules (not per-sample): telocen_norm (fig2) and
    # epi_drift (fig9). Both need exactly 2 groups in Params.<step>.groups.
    for _step in ("telocen_norm", "epi_drift"):
        if _enabled(datajson, _step):
            _gr = (datajson.get("Params", {}).get(_step, {}) or {}).get("groups") or {}
            if len(_gr) != 2 or not all(_gr.values()):
                raise ValueError(
                    f"Fiberseq: Params.{_step}.groups must define exactly 2 "
                    f"non-empty groups (label -> [samples]); got: {_gr!r}. "
                    f"Set Params.{_step}.enabled=false to skip this analysis.")
    if enable_telocen_norm:
        _tn_fmts = (datajson.get("Params", {}).get("telocen_norm", {})
                    .get("fiberseq_telocen_norm_fig", {}).get("formats") or ["png"])
        if isinstance(_tn_fmts, str):
            _tn_fmts = [_tn_fmts]
        outfiles.append(f"{outdir}/fiberseq_analysis/telocen_norm/region_normalized.tsv")
        outfiles.append(f"{outdir}/fiberseq_analysis/telocen_norm/class_diff_normalized.tsv")
        outfiles += _figs(f"{outdir}/fiberseq_analysis/telocen_norm/telocen_normalized",
                          _tn_fmts)
    if enable_epi_drift:
        _ed_fmts = (datajson.get("Params", {}).get("epi_drift", {})
                    .get("fiberseq_epi_drift_fig", {}).get("formats") or ["png"])
        if isinstance(_ed_fmts, str):
            _ed_fmts = [_ed_fmts]
        outfiles.append(f"{outdir}/fiberseq_analysis/epi_drift/epigenetic_drift_windows.tsv")
        outfiles.append(f"{outdir}/fiberseq_analysis/epi_drift/epigenetic_drift.json")
        outfiles += _figs(f"{outdir}/fiberseq_analysis/epi_drift/epigenetic_drift",
                          _ed_fmts)

    if len(genome_versions) != 1:
        raise ValueError(f"meta don't support multiple organsim temporarily, please check your meta file, found: {genome_versions}")
    genome_version = next(iter(genome_versions))
    order_genome_version = datajson["genome"]["default"]
    if genome_version == "UNKNOWN":
        # meta organism column is empty/unrecognized; the literal "UNKNOWN"
        # must never reach Params.sv.exp_specific.species/assembly: vep_install
        # would run with -s UNKNOWN and die with exit 255 "No matching
        # species". Fall back to the config genome.default, which the
        # reference table resolves.
        logger.warning(
            f"meta organism column is UNKNOWN; using genome.default from config.json "
            f"({order_genome_version}) for Params.sv.exp_specific instead of the literal UNKNOWN."
        )
        genome_version = order_genome_version
    if genome_version != order_genome_version:
        logger.warning(
            f"user defined genome versions in meta file ({genome_version}) do not match the genome version defined in config.json ({order_genome_version}) it may come from default value or command-line arguments. "
            f"Using the genome version from meta file ({genome_version}) for downstream analysis."
        )
    ### inject outfiles
    sv_params = datajson.get("Params", {}).get("sv", {})
    _validate_te_lib(datajson, "Fiberseq")
    _fill_mei_rm_species(datajson, "Fiberseq")
    _fill_telomere_species(datajson, "Fiberseq")
    _fill_telomere_subtel_ref(datajson, "Fiberseq")
    _fill_telomere_n_chrom_arms(datajson, "Fiberseq")
    comparisons = dict(sv_params.get("comparisons") or {})
    if enable_sv and not comparisons:
        for gp in group_pairs:
            name = f"{gp.ctr_group_name}_vs_{gp.exp_group_name}"
            comparisons[name] = {
                "control_group_name": gp.ctr_group_name,
                "experimental_group_name": gp.exp_group_name,
                "control_samples": gp.ctr_sample_ids,
                "experimental_samples": gp.exp_sample_ids,
            }
    datajson["Params"].setdefault("sv", {})["comparisons"] = comparisons
    if enable_sv and not comparisons:
        logger.warning(
            "Fiberseq: Params.sv.enabled=true but no SV comparisons are defined; "
            "fill Params.sv.comparisons in the config JSON or the meta design "
            "column (ctrl_TAG / exp_TAG). Skipping SV comparison outputs."
        )
    if enable_sv and comparisons:
        sv_outdir = f"{outdir}/variation/germline_sv_analysis"
        sv_genes = (datajson.get("Params", {}).get("sv", {})
                    .get("sv_gene_model", {}).get("genes") or [])
        sv_oncoprint = sv_params.get("oncoprint") or {}
        known_samples = set(samples_info_dict)
        for name, c in comparisons.items():
            for key in ("control_samples", "experimental_samples"):
                unknown = [s for s in (c.get(key) or []) if s not in known_samples]
                if unknown:
                    logger.warning(
                        f"Fiberseq: Params.sv.comparisons['{name}'].{key} lists "
                        f"samples not in the meta: {unknown}"
                    )
            base = f"{sv_outdir}/{name}"
            outfiles.append(f"{base}/{name}_merged.vcf")
            outfiles.append(f"{base}/{name}_only.vcf")
            outfiles.append(f"{base}/{name}_annotated.tab")
            # TE enrichment (RepeatMasker subfamily test) + per-SV stats from
            # sv_exp_specific — registered so the rule's declared directory
            # outputs are requested by the DAG
            outfiles.append(f"{base}/repeatmasker")
            outfiles.append(f"{base}/statistics")
            outfiles.append(f"{base}/circos")
            # MEI suite (sv_mei rule): classified insertions + 6-figure
            # suite; gate via Params.sv.mei.enabled (default on)
            if (datajson.get("Params", {}).get("sv", {}).get("mei", {}) or {}).get("enabled", True):
                outfiles.append(f"{base}/mei/mei.done")
            if sv_genes:
                outfiles.append(f"{base}/gene_model")
        # SV-affected genes GO/KEGG (modules/function/sv_go; macaque-aware
        # clusterProfiler route). Registered only when the gate is on.
        svgokegg_params = datajson.get("Params", {}).get("function", {}).get("sv_go_kegg", {}) or {}
        if svgokegg_params.get("enabled"):
            svgokegg_species = resolve_genome(genome_version,"alias")
            if svgokegg_species not in ("human", "mouse", "macaque"):
                raise ValueError(
                    "Fiberseq: Params.function.sv_go_kegg.enabled=true requires "
                    "Params.function.sv_go_kegg.species to be one of "
                    f"human/mouse/macaque (got: {svgokegg_species!r})."
                )
            datajson["Params"]["function"]["sv_go_kegg"]["species"] = svgokegg_species
            for name in comparisons:
                outfiles.append(f"{sv_outdir}/{name}/sv_go_kegg/sv_go_kegg.done")
        outfiles.append(f"{sv_outdir}/sv_diff_analysis")
        if sv_oncoprint.get("cosmic_file"):
            outfiles.append(f"{sv_outdir}/sv_oncoprint/oncoprint_matrix.csv")

    # Cross-sample outputs
    if enable_quant or enable_diff:
        for metric in ("percent_accessible", "nuc_occupancy", "m6a_per_fiber", "fire_frac"):
            outfiles.append(f"{outdir}/fiberseq_analysis/quant/matrix/{metric}_matrix.tsv")

    # Modification significance (modules/fiberseq/mods + sigtest): m6A/CpG
    # per-fiber density over FIRE peaks, tested per contrast; significant
    # regions land in a regionplot-ready BED under fiberseq_analysis/mods.
    # Groups/contrasts come from the meta-derived group_pairs (same source
    # sv uses) — the meta is the single source of grouping, no design file.
    # single config name (no compat alias per repo convention: a rename that
    # unifies a key ships no fallback read — regenerate raw.json instead)
    enable_mods = _enabled(datajson, "fiberseq_mods")
    if enable_mods:
        mods_params = (datajson.setdefault("Params", {}).get("fiberseq_mods")
                       or datajson["Params"].setdefault("fiberseq_mods", {}))
        # 参数精确到 rule:mod_type/region_axis 是各 rule 独立键
        # (Params.fiberseq_mods.<rule>.<key>;长表文件名由
        # fiberseq_mods_longtable 产出,以其子键为锚;config 全量路径读取,
        # 无旧扁平键/跨 rule 回退)
        mod_type = ((datajson.get("Params", {}).get("fiberseq_mods", {})
                     .get("fiberseq_mods_longtable", {}) or {})
                    .get("mod_type") or "m6a")
        # region_axis is a discriminant object ({default, options}); the
        # axis id selects the engine, the long-table schema, and which
        # axis-specific params (min_fibers/loci_top_n vs metric/top_n)
        # apply. Flat strings are NOT accepted (closed-set parse error).
        _ra = ((datajson.get("Params", {}).get("fiberseq_mods", {})
                .get("fiberseq_mods_longtable", {}) or {})
               .get("region_axis") or {})
        if not isinstance(_ra, dict):
            raise ValueError(
                "Params.fiberseq_mods.<rule>.region_axis must be a discriminant "
                'object {default, options} (genome-shaped); flat strings '
                "are not accepted")
        region_axis = _ra.get("default") or "fire_peaks"
        axis_opts = (_ra.get("options", {}) or {}).get(region_axis) or {}
        if region_axis not in (_ra.get("options") or {}):
            raise ValueError(
                f"Params.fiberseq_mods.<rule>.region_axis.default must name a key "
                f"of its options (got {region_axis!r}, "
                f"options={sorted(_ra.get('options') or {})})")
        metric = axis_opts.get("metric", f"{mod_type}_per_fiber")
        if region_axis == "te":
            # tested TE levels belong to sigtest (they choose its feature
            # column); read there (config 全量路径,无旧扁平键回退)
            sigtest_params = datajson["Params"].setdefault("sigtest", {})
            te_levels = ([lv for lv in sigtest_params.get("te_levels", ["subfamily", "family"])
                          if lv in ("subfamily", "family")])
            if not te_levels:
                raise ValueError(
                    "Params.sigtest.te_levels must name at least one of "
                    "subfamily/family."
                )
            long_name = f"te_{mod_type}_long.tsv"
        else:
            te_levels = []
            long_name = f"{metric}_long.tsv"
        # grouping: meta design ONLY (group_pairs, same source sv uses)
        if not group_pairs:
            raise ValueError(
                "Params.fiberseq_mods.enabled=true needs the meta "
                "design/group columns filled (they drive group_pairs)."
            )
        mods_contrasts = _contrasts_from_group_pairs(group_pairs)
        # TE axis: one sigtest contrast per (group pair x te level) named
        # "<contrast>.<level>"; fire_peaks axis: bare contrast names.
        contrast_names = sorted(
            f"{c}.{lv}" if lv else c
            for c in mods_contrasts
            for lv in (te_levels or [None])
        )
        # group -> sample lists for the module (long-table group column +
        # top-loci group specs); parse-time expand needs the contrast list
        # contrasts 是接线键(subworkflow 转发为模块顶层 contrasts),非 rule 参数
        mods_params["contrasts"] = contrast_names
        # 写回精确到 rule:group_samples 消费于 fiberseq_mods_longtable 与
        # fiberseq_mods_top_loci;groups_by_contrast 仅 top_loci 消费(其规则
        # 只读自己的子键——多组设计否则每个 contrast 都回退"前两组")
        _group_samples = {
            g: samples for ctr in mods_contrasts.values() for g, samples in ctr.items()
        }
        _groups_by_contrast = {
            f"{c}.{lv}" if lv else c: sorted(groups)
            for c, groups in mods_contrasts.items()
            for lv in (te_levels or [None])
        }
        for _mr in ("fiberseq_mods_longtable", "fiberseq_mods_top_loci"):
            mods_params.setdefault(_mr, {})["group_samples"] = _group_samples
        mods_params.setdefault("fiberseq_mods_top_loci", {})["groups_by_contrast"] = _groups_by_contrast
        mods_base = f"{outdir}/fiberseq_analysis/mods"
        # sigtest module wiring (subworkflow use rule): its input table is
        # our long table; group pair per contrast comes from the same map,
        # feature column per contrast from the level suffix
        sigtest_params = datajson["Params"].setdefault("sigtest", {})
        sigtest_params["sigtest_table"] = f"{mods_base}/{long_name}"
        sigtest_params["contrasts"] = contrast_names
        sigtest_params["groups_by_contrast"] = {
            f"{c}.{lv}" if lv else c: sorted(groups)
            for c, groups in mods_contrasts.items()
            for lv in (te_levels or [None])
        }
        if te_levels:
            sigtest_params["feature_col_by_contrast"] = {
                f"{c}.{lv}": lv for c in mods_contrasts for lv in te_levels
            }
        # sigtest plot formats ([] = no figure; absorbs the old no_plot)
        _st_plot = (sigtest_params.get("plot", {}) or {})
        _plot_fmts = _st_plot.get("formats")
        if _plot_fmts is None:
            _plot_fmts = ["png"]
        if isinstance(_plot_fmts, str):
            _plot_fmts = [_plot_fmts]
        outfiles.append(f"{mods_base}/{long_name}")
        # sigtest subtree: significance tables/figures AND the mods
        # postprocess products (significant_objects/top_loci) all live
        # under sigtest/<contrast>/ — postprocess outputs are classified
        # by their DATA (sigtest results), not by the producing module
        for contrast in contrast_names:
            sig_base = f"{outdir}/fiberseq_analysis/sigtest/{contrast}"
            outfiles.append(f"{sig_base}/significance_test.tsv")
            for fmt in _plot_fmts:
                # per-format figures live under significance_bars/ +
                # optional significance_bars_overview.{fmt}; all tracked
                # by the done marker the rule declares
                pass
            if _plot_fmts:
                outfiles.append(f"{sig_base}/significance_bars.done")
            outfiles.append(f"{sig_base}/significant_objects.tsv")
            if region_axis != "te":
                # BED only on the fire_peaks axis (coordinate features);
                # TE features are family names -- no BED
                outfiles.append(f"{sig_base}/significant_objects.bed")
            outfiles.append(f"{sig_base}/nonsignificant_objects.tsv")
            if region_axis == "te":
                outfiles.append(f"{sig_base}/top_loci.tsv")
                outfiles.append(f"{sig_base}/top_loci.bed")

    if enable_diff:
        # grouping: meta design ONLY (group_pairs) — same source mods uses
        if not group_pairs:
            raise ValueError(
                "Params.diff.enabled=true needs the meta "
                "design/group columns filled (they drive group_pairs)."
            )
        contrast_groups = _contrasts_from_group_pairs(group_pairs)
        # contrast spec string for the engine (diff_accessibility --group-spec)
        # group_spec 精确到 rule:写回 fiberseq_diff_run 子键(消费 rule 经
        # rule_param 读取)
        datajson["Params"].setdefault("diff", {}).setdefault(
            "fiberseq_diff_run", {})["group_spec"] = ",".join(
            f"{c}:" + "|".join(f"{g}:{'+'.join(ss)}" for g, ss in gs.items())
            for c, gs in contrast_groups.items()
        )
        for contrast, groups in sorted(contrast_groups.items()):
            base = f"{outdir}/fiberseq_analysis/diff/{contrast}"
            outfiles.append(f"{base}/diff_{contrast}.tsv")
            _diff_fmts = (datajson.get("Params", {}).get("diff", {})
                          .get("fiberseq_diff_run", {}).get("formats") or ["pdf"])
            if isinstance(_diff_fmts, str):
                _diff_fmts = [_diff_fmts]
            outfiles += _figs(f"{base}/volcano_{contrast}", _diff_fmts)
            outfiles.append(f"{base}/class_summary_{contrast}.tsv")
            outfiles.append(f"{base}/TEcount/{contrast}.TEcount_Gene.name.tsv")

        # function module needs group_pairs keyed by contrast; derive when absent
        # (group order matches diff_accessibility: sorted group names -> control,
        # experimental)
        if not (datajson.get("Params", {}).get("function", {}) or {}).get("group_pairs"):
            gp = {}
            for contrast, groups in contrast_groups.items():
                (g_ctrl, s_ctrl), (g_exp, s_exp) = sorted(groups.items())
                gp[contrast] = {
                    "control_group_name": g_ctrl,
                    "experimental_group_name": g_exp,
                    "control_samples": s_ctrl,
                    "experimental_samples": s_exp,
                }
            datajson["Params"].setdefault("function", {})["group_pairs"] = gp

        # GO/KEGG + GSEA outputs from the shared function module
        if enable_function:
            for contrast in sorted(contrast_groups):
                base = f"{outdir}/function/{contrast}"
                _fn_go_fmts = (datajson.get("Params", {}).get("function", {})
                               .get("function_go_kegg", {}).get("formats") or ["png"])
                if isinstance(_fn_go_fmts, str):
                    _fn_go_fmts = [_fn_go_fmts]
                outfiles += _figs(f"{base}/go_back_to_back", _fn_go_fmts)
                outfiles += _figs(f"{base}/kegg_back_to_back", _fn_go_fmts)
                outfiles.append(f"{base}/go_up.csv")
                outfiles.append(f"{base}/go_down.csv")
                outfiles.append(f"{base}/kegg_up.csv")
                outfiles.append(f"{base}/kegg_down.csv")
                outfiles.append(f"{base}/up_genes.txt")
                outfiles.append(f"{base}/down_genes.txt")
            if ((datajson.get("Params", {}).get("function", {})
                 .get("function_gsea", {}) or {}).get("gmt")
                    and default_ref.get("geneIDAnno")):
                for contrast in sorted(contrast_groups):
                    base = f"{outdir}/function/{contrast}/GSEA"
                    _fn_gsea_fmts = (datajson.get("Params", {}).get("function", {})
                                     .get("function_gsea", {}).get("formats") or ["jpeg"])
                    if isinstance(_fn_gsea_fmts, str):
                        _fn_gsea_fmts = [_fn_gsea_fmts]
                    outfiles += _figs(f"{base}/TEcount_Gene_GSEA", _fn_gsea_fmts)
                    outfiles.append(f"{base}/TEcount_Gene_GSEA.csv")

    if enable_lsdc:
        prefix = ((datajson.get("Params", {}).get("ldsc", {})
                   .get("ldsc_annot", {}) or {}).get("annot_prefix") or "peaks")
        outfiles.append(f"{outdir}/ldsc/annot/{prefix}.annot_manifest.tsv")
        for trait in sorted((datajson.get("Params", {}).get("ldsc", {}).get("gwas") or {}).keys()):
            outfiles.append(f"{outdir}/ldsc/{trait}/{trait}.sumstats.gz")
            outfiles.append(f"{outdir}/ldsc/{trait}/{trait}.cts.results")

    datajson["Params"]["sv"]["exp_specific"]["species"] = genome_version
    datajson["Params"]["sv"]["exp_specific"]["assembly"] = genome_version
    datajson["Params"]["workflow"]["samples"] = samples
    datajson["Params"]["workflow"]["raw_files"] = raw_files
    outfiles += _register_regionplot(datajson, samples_info_dict, samples, outdir, "Fiberseq")
    outfiles += _register_fiberbrowser(datajson, samples_info_dict, samples, outdir, "Fiberseq")
    # Fiber pattern clustering + report package (modules/fibercluster):
    # per-setting done markers + report files; gated by Params.fibercluster.
    outfiles += _register_fibercluster(datajson, outdir, "Fiberseq")
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runLRtranscriptome(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ) -> str:
    """Prepare input JSON for LRtranscriptome (Long-Read Transcriptome) workflow.

    Long-read full-length transcriptome analysis pipeline:
    1. minimap2 alignment (PacBio HiFi / ONT)
    2. StringTie assembly (long-read mode, -L)
    3. StringTie merge (merge transcripts across samples)
    4. StringTie quantification (-e -B -L, re-estimate with merged GTF)

    Supports both PacBio (map-hifi) and ONT (map-ont) sequencing data.
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    wf["indir"] = indir
    wf["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    wf["logdir"] = logdir
    tmp_dir = datajson.get("tmp_dir") or os.path.join(outdir, "tmp")
    if not datajson.get("tmp_dir"):
        os.makedirs(tmp_dir, exist_ok=True)
    datajson["tmp_dir"] = tmp_dir

    outfiles = []
    samples = []
    skip_quant = datajson.get("Params", {}).get("skip_quant", False)

    for sample_id, sample_info in samples_info_dict.items():
        samples.append(sample_id)
        # Alignment output
        outfiles.append(f"{outdir}/alignment/{sample_id}/{sample_id}.sorted.bam")
        outfiles.append(f"{outdir}/alignment/{sample_id}/{sample_id}.sorted.bam.bai")
        # Assembly output
        outfiles.append(f"{outdir}/stringtie/raw/{sample_id}/{sample_id}.gtf")

    # Merged transcriptome (across all samples)
    outfiles.append(f"{outdir}/stringtie/stringtie_merged.gtf")

    # Quantification output
    if not skip_quant:
        for sample_id in samples:
            outfiles.append(f"{outdir}/quantification/{sample_id}/{sample_id}.gene_abund.tab")
            outfiles.append(f"{outdir}/quantification/{sample_id}/{sample_id}.quant.gtf")

    wf["samples"] = samples
    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the alignment BAMs; registered for every sample BAM.
    outfiles.append(f"{outdir}/QC/sexcheck/sexcheck.tsv")
    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.html")
    datajson["outfiles"] = outfiles
    return _write_raw_json(datajson, outdir)


def runDR(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        designPairs: List[DesignPair],
        group_pairs: List[CompareGroupPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ) -> str:
    """Prepare input JSON for the DR workflow (DNA x RNA joint evidence).

    End-to-end pipeline (subworkflow/DR.smk):
    0. DNA chain (data_type=DNA rows) from raw fastq (dna_indir) mirroring the
       full Mutation.smk analysis set: fastqc -> cutadapt -> bwa-mem2 -> samtools
       -> gatk_prepare markdup -> targets -> Mutect2+Filter -> msi ->
       gatk_germline -> somatic spectrum/signature -> fragment_size ->
       Manta (germline + somatic SV) -> CNVkit -> burden -> sexcheck ->
       qc_summary -> Mutation_report. Outputs under {outdir}/DNA/.
       Assay-aware via the meta assay column (WGS / Exome / Panel): it drives
       CNVkit --seq-method and the exome capture targets BED; nothing is
       WES-specific by default.
    1. RNA chain (data_type=RNA rows) from raw fastq (indir), replicating the
       full RNAseq.smk analysis set: trimmer -> TE quantification (TEtranscripts)
       -> featureCounts -> DESeq2 -> function (GO/KEGG/GSEA/GSVA) ->
       StringTie -> fusion (RmrRNA/bowtie2/STAR-chimera/arriba) -> RNA SNP
       (star -> XenofilteR -> markdup -> gatk_RNAseq) -> RNAseq_report ->
       sexcheck -> qc_summary. Outputs under {outdir}/RNA/.
    2. DRclass module: anchors -> bcftools pileup -> state-machine
       classification -> per-sample tables + summary ({outdir}/DRclass/).

    Meta conventions:
      - Chain routing is by the meta data_type column (required on every row,
        case-insensitive): DNA rows run chain 0, RNA rows run chain 1.
        MetaUtil stages each row's fastq with <ID>_{1,2}.fq.gz naming into
        {outdir}/DNA/common/1_raw_data/ or {outdir}/RNA/common/1_raw_data/.
      - design keeps the project-wide ctr_<tag>/exp_<tag> semantics
        (control/experiment, Mutation parity): DNA design pairs are matched
        normal/tumor somatic pairs (control = normal). Design pairs are never
        used for chain routing; mixed DNA-RNA pairs are ignored with a warning
        (DRclass evidence pairing is expressed by group or Params.workflow
        .dr_pairs instead).
      - RNA differential expression contrasts: hand-filled
        Params.DESeq2.group_pairs wins, else RNA-only design group pairs,
        else all pairs of the RNA rows' group column (alphabetical first =
        control).
      - DRclass evidence pairing: a DNA sample and a RNA sample with the SAME
        sample name (different data_type) are the two assays of one biological
        sample and form one DRclass pairing; hand-filled Params.workflow
        .dr_pairs ({dna_name: rna_name}) overrides. Entries need a somatic VCF,
        so only somatic tumors of the DNA design pairs get classified.
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files, workflow="DR")
    _resolve_genome_ref(datajson, samples_info_dict, "DR")

    # --- Chain routing by meta data_type (the ONLY routing axis) ---
    dna_samples: List[str] = []
    rna_samples: List[str] = []
    for sample_id, info in samples_info_dict.items():
        dt = (getattr(info, "data_type", None) or "").strip().upper()
        if dt == "DNA":
            dna_samples.append(sample_id)
        elif dt == "RNA":
            rna_samples.append(sample_id)
        else:
            raise ValueError(
                f"DR: sample {sample_id} has data_type={getattr(info, 'data_type', None)!r}; "
                "every DR meta row must carry data_type=DNA or data_type=RNA "
                "(case-insensitive) — it is the only chain routing axis")
    if not dna_samples or not rna_samples:
        raise ValueError(
            f"DR: meta must contain both data_type=DNA rows (got {len(dna_samples)}) "
            f"and data_type=RNA rows (got {len(rna_samples)})")

    # samples_dict keys are unique sample keys; modules/paths use the real
    # sample name (display). The same display name under two data_types is one
    # biological sample assayed twice — the DRclass pairing key.
    def _disp(key: str) -> str:
        return samples_info_dict[key].sample_id
    dna_info = {_disp(k): samples_info_dict[k] for k in dna_samples}
    rna_info_by_name = {_disp(k): samples_info_dict[k] for k in rna_samples}

    outfiles = []
    dna_paired_samples: List[str] = []
    dna_single_samples: List[str] = []
    genome_paired_samples = {}
    genome_single_samples = {}
    paired_design = []          # [{"control": normal_dna, "experimental": tumor_dna}]
    pair_strains = {}
    drclass_samples = {}
    strain_vcfs = datajson.get("Params", {}).get("workflow", {}).get("drclass_strain_vcfs", {})

    # --- DNA design pairs (standard ctr_/exp_ semantics, Mutation parity) ---
    # Both endpoints DNA -> matched normal/tumor somatic pair. Cross-type or
    # RNA-only pairs are not chain-routing devices; cross-type ones are ignored.
    dna_set, rna_set = set(dna_samples), set(rna_samples)
    dna_design_pairs = []
    for designPair in designPairs:
        ctr_id, exp_id = designPair.ctr_sample_id, designPair.exp_sample_id
        if ctr_id in dna_set and exp_id in dna_set:
            dna_design_pairs.append(designPair)
        elif ctr_id in rna_set and exp_id in rna_set:
            continue  # RNA contrast pair; DE uses group_pairs below
        else:
            logger.warning(
                f"DR: ignoring mixed DNA/RNA design pair {ctr_id} vs {exp_id}; "
                "design pairs never route chains (data_type does) and DRclass "
                "evidence pairing comes from group or Params.workflow.dr_pairs")
    for designPair in dna_design_pairs:
        paired_design.append({"control": _disp(designPair.ctr_sample_id),
                              "experimental": _disp(designPair.exp_sample_id)})
        strain = getattr(samples_info_dict.get(designPair.exp_sample_id), "strain", None)
        if strain:
            pair_strains[_disp(designPair.exp_sample_id)] = strain
    if pair_strains:
        logger.info(f"DR: strain-aware germline filtering enabled for pairs: {pair_strains}")

    # --- Raw staging roots (MetaUtil data_type routing) ---
    # Derived from where the fastq were actually linked: MetaUtil is the
    # single authority on the staging layout. Falls back to the shared root
    # when nothing was staged (the per-sample verification below reports).
    def _raw_staging_root(sample_ids: List[str]) -> str:
        for sid in sample_ids:
            info = samples_info_dict.get(sid)
            if info and info.fastq_1:
                return os.path.dirname(os.path.dirname(info.fastq_1))
        return os.path.join(outdir, "common", "1_raw_data")

    dna_raw_staging = _raw_staging_root(dna_samples)
    rna_raw_staging = _raw_staging_root(rna_samples)
    # Chain consumption: chain 1 reads indir (RNA root), chain 0 reads the
    # dna_indir key (DNA root) — both rewritten here after staging.
    datajson["Params"]["workflow"]["indir"] = rna_raw_staging
    datajson["Params"]["workflow"]["dna_indir"] = dna_raw_staging

    # --- DNA chain: verify staged fastq, classify layouts ---
    for key in dna_samples:
        dna_id = _disp(key)
        info = samples_info_dict[key]
        if info.layout == "PE":
            dna_paired_samples.append(dna_id)
        elif info.layout == "SE":
            dna_single_samples.append(dna_id)
        else:
            raise ValueError(f"DR: unknown layout for DNA sample {dna_id}: {info.layout}")
        staged = [info.fastq_1] + ([info.fastq_2] if info.layout == "PE" else [])
        if not all(p and os.path.exists(p) for p in staged):
            raise ValueError(f"DR: DNA fastq not staged for {dna_id} under {dna_raw_staging} "
                             f"(expected {staged}; meta must carry fastq_1/fastq_2 for DNA rows)")

    # --- RNA chain sample grouping ---
    for key in rna_samples:
        rna_id = _disp(key)
        rna_info = samples_info_dict[key]
        organism = rna_info.organism
        if rna_info.layout == "PE":
            genome_paired_samples.setdefault(organism, []).append(rna_id)
        elif rna_info.layout == "SE":
            genome_single_samples.setdefault(organism, []).append(rna_id)
        else:
            raise ValueError(f"DR: unknown layout for RNA sample {rna_id}: {rna_info.layout}")

    # --- RNA chain analysis products (RNAseq.smk parity) ---
    # Counts / DE / enrichment / assembly / fusion / report / QC registration
    # mirrored from runRNAseq, every path scoped under {outdir}/RNA. DE
    # contrasts (DESeq2.group_pairs): config wins, else the RNA-only design
    # group pairs, else all pairs of the RNA rows' group column (alphabetical
    # first = control). Function enrichment stays gated by
    # Params.function.<organism>.enabled like runRNAseq.
    enable_te = _enabled(datajson, "TEtranscripts")
    enable_fusion = _enabled(datajson, "arriba")
    enable_deseq2 = _enabled(datajson, "DESeq2")
    enable_TEcount = _enabled(datajson, "TEcount")
    enable_TElocal = _enabled(datajson, "TElocal")
    enable_report = _enabled(datajson, "RNAseq_report")
    enable_StringTie = _enabled(datajson, "StringTie")

    _deseq_params = datajson["Params"].setdefault("DESeq2", {})
    if not _deseq_params.get("group_pairs"):
        _derived = {}
        # 1) RNA-only design group pairs (runRNAseq semantics)
        for gp in group_pairs:
            if set(gp.ctr_sample_ids) <= rna_set and set(gp.exp_sample_ids) <= rna_set:
                _derived.setdefault(gp.organism, {})[f"{gp.ctr_group_name}_vs_{gp.exp_group_name}"] = {
                    "control_group_name": gp.ctr_group_name,
                    "experimental_group_name": gp.exp_group_name,
                    "control_samples": [_disp(k) for k in gp.ctr_sample_ids],
                    "experimental_samples": [_disp(k) for k in gp.exp_sample_ids],
                }
        _source = "RNA design group pairs"
        # 2) fall back to all pairs of the RNA rows' group column
        if not _derived:
            _rna_groups = {}
            for sid in rna_samples:
                _g = (samples_info_dict[sid].group or "").strip()
                if _g:
                    _rna_groups.setdefault(samples_info_dict[sid].organism, {}).setdefault(_g, []).append(_disp(sid))
            for _org, _gmap in _rna_groups.items():
                _names = sorted(_gmap)
                for _i, _ctr in enumerate(_names):
                    for _exp in _names[_i + 1:]:
                        _derived.setdefault(_org, {})[f"{_ctr}_vs_{_exp}"] = {
                            "control_group_name": _ctr,
                            "experimental_group_name": _exp,
                            "control_samples": list(_gmap[_ctr]),
                            "experimental_samples": list(_gmap[_exp]),
                        }
            _source = "group column"
        _deseq_params["group_pairs"] = _derived
        if _derived:
            logger.info(f"DR: RNA DE contrasts derived from {_source}: "
                        f"{ {o: list(v) for o, v in _derived.items()} } "
                        "(override via Params.DESeq2.group_pairs)")
    rna_group_pairs = _deseq_params.get("group_pairs") or {}
    if not rna_group_pairs:
        logger.info("DR: no RNA DE contrasts (need Params.DESeq2.group_pairs, RNA "
                    "design pairs, or >=2 groups in the RNA rows' group column); "
                    "DE/function pair outputs skipped")

    # StringTie sample grouping from the RNA rows' group column
    _sample_groups = {}
    for sid in rna_samples:
        _g = (samples_info_dict[sid].group or "").strip()
        if _g:
            _sample_groups.setdefault(_g, []).append(_disp(sid))
    datajson["Params"].setdefault("StringTie", {})["sample_groups"] = _sample_groups

    # XenofilteR host/contaminant map (runRNAseq parity; contaminant defaults
    # to the host genome when the meta leaves contaminated_organism empty)
    _sample_contamination = {}
    for sid in rna_samples:
        _info = samples_info_dict[sid]
        _sample_contamination[_disp(sid)] = {
            "host": _info.organism,
            "contaminant": _info.contaminated_organism or _info.organism,
        }
    datajson["Params"].setdefault("XenofilteR", {})["sample_contamination"] = _sample_contamination

    # per-contrast DE + function outputs
    for _org, _pairs in rna_group_pairs.items():
        for _pair_name in _pairs:
            if enable_deseq2 and enable_TEcount:
                outfiles.append(f"{outdir}/RNA/diff_expression/{_org}/{_pair_name}/TEcount/DESeq2.done")
            if enable_deseq2 and enable_TElocal:
                _te_dir = f"{outdir}/RNA/diff_expression/{_org}/{_pair_name}/TElocal"
                outfiles.append(f"{_te_dir}/DESeq2.done")
                _m_fmts = (datajson.get("Params", {}).get("DESeq2", {})
                           .get("DESeq2_TElocal_manhattan", {}).get("formats") or ["png"])
                if isinstance(_m_fmts, str):
                    _m_fmts = [_m_fmts]
                outfiles += _figs(f"{_te_dir}/manhattan/{_pair_name}.manhattan", _m_fmts)
                # locus boxplot / subfamily bar only run when their per-contrast
                # family patterns are configured (Params.DESeq2.<rule>.families)
                for _te_fig_rule, _te_fig_dir in (("DESeq2_TElocal_locus_boxplot", "locus_box"),
                                                  ("DESeq2_TElocal_subfamily_bar", "subfamily_bar")):
                    _fam = (datajson.get("Params", {}).get("DESeq2", {}).get(_te_fig_rule, {}) or {}).get("families", {}) or {}
                    if _fam.get(_org, {}).get(_pair_name):
                        _f_fmts = (datajson.get("Params", {}).get("DESeq2", {}).get(_te_fig_rule, {}).get("formats") or ["png"])
                        if isinstance(_f_fmts, str):
                            _f_fmts = [_f_fmts]
                        outfiles += _figs(f"{_te_dir}/{_te_fig_dir}/{_pair_name}.{_te_fig_dir}", _f_fmts)
            if (datajson.get("Params", {}).get("function", {}).get(_org, {}) or {}).get("enabled", False):
                _pair_dir = f"{outdir}/RNA/function/{_org}/{_pair_name}"
                _go_fmts = (datajson.get("Params", {}).get("function", {}).get(_org, {})
                            .get("function_go_kegg", {}).get("formats") or ["png"])
                if isinstance(_go_fmts, str):
                    _go_fmts = [_go_fmts]
                outfiles += _figs(f"{_pair_dir}/go_back_to_back", _go_fmts)
                outfiles += _figs(f"{_pair_dir}/kegg_back_to_back", _go_fmts)
                _gsea_fmts = (datajson.get("Params", {}).get("function", {}).get(_org, {})
                              .get("function_gsea", {}).get("formats") or ["jpeg"])
                if isinstance(_gsea_fmts, str):
                    _gsea_fmts = [_gsea_fmts]
                outfiles += _figs(f"{_pair_dir}/GSEA/TEcount_Gene_GSEA", _gsea_fmts)
                outfiles.append(f"{_pair_dir}/go_up.csv")
                outfiles.append(f"{_pair_dir}/go_down.csv")
                outfiles.append(f"{_pair_dir}/kegg_up.csv")
                outfiles.append(f"{_pair_dir}/kegg_down.csv")
                outfiles.append(f"{_pair_dir}/up_genes.txt")
                outfiles.append(f"{_pair_dir}/down_genes.txt")
                outfiles.append(f"{_pair_dir}/GSEA/TEcount_Gene_GSEA.csv")

    # per-sample RNA products (the variation VCF + BAM are registered with the
    # DRclass map below as chain evidence)
    for key in rna_samples:
        sid = _disp(key)
        _org = samples_info_dict[key].organism
        if enable_StringTie:
            outfiles.append(f"{outdir}/RNA/transcripts/{_org}/raw/{sid}/{sid}_TE_chimeric_transcripts.txt")
        if enable_fusion:
            outfiles.append(f"{outdir}/RNA/fusion/{_org}/{sid}/{sid}_passed_fusions.tsv")

    # per-organism RNA products
    for _org in dict.fromkeys(samples_info_dict[sid].organism for sid in rna_samples):
        # function species mirror runRNAseq (checked only when function is
        # enabled for the organism; unsupported species hard-error like there)
        _fn_org = datajson["Params"].setdefault("function", {}).setdefault(_org, {})
        if _fn_org.get("enabled", False):
            if _org in ("Homo sapiens", "human", "GRCh38"):
                _fn_species = "human"
            elif _org in ("Mus musculus", "mouse", "GRCm39"):
                _fn_species = "mouse"
            else:
                raise ValueError(f"function module don't support {_org}, only support human or mouse(Homo sapiens or Mus musculus)")
            for _fn_rule in ("function_go_kegg", "function_gsva"):
                _fn_org.setdefault(_fn_rule, {})["species"] = _fn_species
            outfiles.append(f"{outdir}/RNA/function/{_org}/gsva/gsva_scores.tsv")
            _gsva_fmts = (_fn_org.get("function_gsva", {}).get("formats") or ["png"])
            if isinstance(_gsva_fmts, str):
                _gsva_fmts = [_gsva_fmts]
            outfiles += _figs(f"{outdir}/RNA/function/{_org}/gsva/gsva_heatmap", _gsva_fmts)
        if enable_fusion:
            outfiles.append(f"{outdir}/RNA/fusion/{_org}/arriba_report/arriba_fusion_report.html")
        if enable_StringTie:
            outfiles.append(f"{outdir}/RNA/transcripts/{_org}/stringtie_merged.gtf")
            if rna_group_pairs.get(_org):
                _st_fmts = (datajson.get("Params", {}).get("StringTie", {}).get("formats") or ["png"])
                if isinstance(_st_fmts, str):
                    _st_fmts = [_st_fmts]
                for _st_fig in ("TE_chimeric_group_stacked",
                                "TE_chimeric_te_type_top",
                                "TE_chimeric_te_type_by_group"):
                    outfiles += _figs(f"{outdir}/RNA/transcripts/{_org}/TE_chimeric/{_st_fig}", _st_fmts)
                outfiles.append(f"{outdir}/RNA/transcripts/{_org}/TE_chimeric/TE_chimeric_sample_summary.tsv")
                outfiles.append(f"{outdir}/RNA/transcripts/{_org}/TE_chimeric/TE_chimeric_group_summary.tsv")
                outfiles.append(f"{outdir}/RNA/transcripts/{_org}/TE_chimeric/TE_chimeric_te_type_counts.tsv")
        if enable_te or enable_deseq2:
            outfiles.append(f"{outdir}/RNA/counts/TEtranscripts/{_org}/{_org}_TEcount.tsv")
            outfiles.append(f"{outdir}/RNA/counts/TEtranscripts/{_org}/{_org}_TEcount_name.tsv")
        if enable_TElocal:
            outfiles.append(f"{outdir}/RNA/counts/TEtranscripts/{_org}/{_org}_TElocal.tsv")
            outfiles.append(f"{outdir}/RNA/counts/TEtranscripts/{_org}/{_org}_TElocal_name.tsv")
        # gene expression quantification (always delivered)
        outfiles.append(f"{outdir}/RNA/counts/featureCounts/{_org}/{_org}_featureCounts.tsv")
        if enable_report:
            outfiles.append(f"{outdir}/RNA/results/{_org}/RNAseq_report.pptx")
            outfiles.append(f"{outdir}/RNA/results/{_org}/figures/manifest.tsv")

    # RNA report Params (Params.RNAseq_report 段; runRNAseq 同款 genome/date 写回)
    if enable_report:
        _rna_report = datajson["Params"].setdefault("RNAseq_report", {}).setdefault("generate_report", {})
        _rna_report["genome"] = "".join(dict.fromkeys(
            samples_info_dict[s].organism for s in rna_samples))
        _rna_report["date"] = time.strftime("%Y-%m-%d", time.localtime())

    # Sex check + QC summary (modules/sexcheck, modules/qc_summary): mirrored
    # from runRNAseq, registered unconditionally for the delivered output set.
    outfiles.append(f"{outdir}/RNA/QC/sexcheck/sexcheck.tsv")
    outfiles.append(f"{outdir}/RNA/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/RNA/QC/qc_summary/qc_summary.html")

    # --- DRclass per-sample explicit file map (chain outputs) ---
    # Evidence pairing rule: a DNA sample and a RNA sample sharing the same
    # sample name (different data_type) are the two assays of one biological
    # sample -> one DRclass pairing. Hand-filled Params.workflow.dr_pairs
    # ({dna_name: rna_name}) overrides. Entries need a somatic VCF, so only
    # somatic tumors of the DNA design pairs get classified.
    _paired_map = {
        rec["experimental"]: rec["control"]
        for rec in paired_design
        if rec.get("experimental") and rec.get("control")
    }
    dr_pairs_cfg = datajson.get("Params", {}).get("workflow", {}).get("dr_pairs") or {}
    if dr_pairs_cfg:
        dr_pairs = dict(dr_pairs_cfg)
        for _dna_name, _rna_name in dr_pairs.items():
            if _dna_name not in dna_info or _rna_name not in rna_info_by_name:
                raise ValueError(
                    f"DR: Params.workflow.dr_pairs entry {_dna_name} -> {_rna_name} must map a "
                    f"data_type=DNA sample name to a data_type=RNA sample name")
        logger.info(f"DR: DRclass evidence pairing taken from Params.workflow.dr_pairs: {dr_pairs}")
    else:
        dr_pairs = {_n: _n for _n in sorted(dna_info) if _n in rna_info_by_name}
        if dr_pairs:
            logger.info(f"DR: DRclass pairing by same sample name across data_type: {dr_pairs}")
        else:
            logger.info("DR: no DNA/RNA sample-name matches, no DRclass pair (give the two "
                        "assays of one sample the same name, or set Params.workflow.dr_pairs)")
    _tumors_without_rna = [t for t in sorted(_paired_map) if t not in dr_pairs]
    if _tumors_without_rna:
        logger.warning(f"DR: somatic tumors without same-name RNA evidence, no DRclass "
                       f"entry: {_tumors_without_rna}")

    for tumor_name, rna_name in dr_pairs.items():
        if tumor_name not in _paired_map:
            logger.info(f"DR: {tumor_name} is not a somatic tumor of a DNA design pair "
                        "(no somatic VCF), skipping DRclass entry")
            continue
        rna_info = rna_info_by_name[rna_name]
        organism = rna_info.organism
        control_id = _paired_map[tumor_name]
        somatic_vcf = (f"{outdir}/DNA/variation/somatic_snv_indel/"
                       f"{control_id}_vs_{tumor_name}/{control_id}_vs_{tumor_name}.filtered.vcf.gz")
        germline_vcf = (f"{outdir}/DNA/variation/germline_snv_indel/"
                        f"{tumor_name}/{tumor_name}.raw.vcf.gz")
        rna_vcf = (f"{outdir}/RNA/variation/germline_snv_indel_RNAseq/"
                   f"{organism}/{rna_name}/{rna_name}.filtered.vcf.gz")
        rna_bam = (f"{outdir}/RNA/common/10_markdup_deContaimation_bam/"
                   f"{organism}/{rna_name}/{rna_name}.sorted_markdup.bam")
        drclass_samples[tumor_name] = {
            "somatic_vcf": somatic_vcf,
            "germline_vcf": germline_vcf,
            "rna_vcf": rna_vcf,
            "rna_bam": rna_bam,
            "rna_bam_index": rna_bam.replace(".bam", ".bai"),
            "strain_vcf": strain_vcfs.get(tumor_name) or strain_vcfs.get("*") or [],
            "splice_bed": [],
            "editing_bed": [],
        }
        outfiles += [
            f"{outdir}/DRclass/{tumor_name}/{tumor_name}.dna_anchored.tsv",
            f"{outdir}/DRclass/{tumor_name}/{tumor_name}.rna_anchored.tsv",
        ]
        # chain-1 RNA products (so the DAG builds them)
        outfiles += [rna_vcf, rna_bam]

    # --- chain-0 DNA products: full Mutation.smk analysis set ---
    # Sequencing assay -> CNVkit --seq-method (Mutation parity; judged over the
    # DNA rows only — RNA rows share the meta). The meta `assay` column wins
    # over the config value; all DNA rows must share one assay. DNA is not
    # assumed to be WES: an absent assay column keeps the config value as-is.
    _ref_entry = datajson["genome"]["references"][datajson["genome"]["default"]]
    _ASSAY_TO_CNVKIT_METHOD = {
        "wgs": "wgs", "exome": "hybrid", "wes": "hybrid",
        "panel": "amplicon", "amplicon": "amplicon",
    }
    _dna_assays = sorted({
        (getattr(samples_info_dict.get(sid), "assay", "") or "").strip().lower()
        for sid in dna_samples
    } - {""})
    if len(_dna_assays) > 1:
        raise ValueError(
            f"DR: meta assay column holds mixed assay types {_dna_assays}; "
            "CNVkit needs a single sequencing method per run."
        )
    _dna_assay_method = (_ASSAY_TO_CNVKIT_METHOD.get(_dna_assays[0])
                         if _dna_assays else None)
    if _dna_assays and _dna_assay_method is None:
        raise ValueError(
            f"DR: unknown meta assay '{_dna_assays[0]}'; expected one of "
            f"{sorted(_ASSAY_TO_CNVKIT_METHOD)}"
        )
    # Capture region (Mutation parity): an explicit genome.references.<g>.targets
    # BED wins; exome/hybrid assays derive it from markdup BAM coverage
    # (DR.smk derive_targets), registered here for GATK -L, CNVkit --targets
    # and burden. wgs/amplicon assays and an absent assay column have no
    # capture region (DNA is not assumed to be WES).
    if _ref_entry.get("targets"):
        _targets = _ref_entry["targets"]
    elif _dna_assay_method == "hybrid":
        _targets = f"{outdir}/DNA/common/5_targets_bed/targets.bed"
        outfiles.append(_targets)
        _ref_entry["targets"] = _targets
    else:
        _targets = None
        _ref_entry["targets"] = None

    # Somatic (Mutect2) outputs per design pair + the sample -> file maps the
    # spectrum/burden modules consume.
    sample_somatic_vcf_dict = {}
    sample_group_dict = {}
    for rec in paired_design:
        somatic_file = (f"{outdir}/DNA/variation/somatic_snv_indel/"
                        f"{rec['control']}_vs_{rec['experimental']}/"
                        f"{rec['control']}_vs_{rec['experimental']}.filtered.vcf.gz")
        sample_somatic_vcf_dict[rec["experimental"]] = somatic_file
        sample_group_dict[rec["experimental"]] = (
            getattr(dna_info.get(rec["experimental"]), "group", "") or "")
        outfiles.append(somatic_file)
    # Germline (HaplotypeCaller) raw + filtered VCF for every DNA sample:
    # raw is the DRclass anchor, filtered completes the Mutation germline
    # analysis for delivery.
    for key in dna_samples:
        dna_id = _disp(key)
        outfiles.append(f"{outdir}/DNA/variation/germline_snv_indel/"
                        f"{dna_id}/{dna_id}.raw.vcf.gz")
        outfiles.append(f"{outdir}/DNA/variation/germline_snv_indel/"
                        f"{dna_id}/{dna_id}.filtered.vcf.gz")

    # Alignment QC (samtools flagstat) + hybrid-selection QC (CollectHsMetrics
    # over the capture region; Mutation parity — also the qc_summary anchors).
    for key in dna_samples:
        dna_id = _disp(key)
        outfiles.append(f"{outdir}/DNA/common/3_sorted_bam/{dna_id}/{dna_id}_flagstat.txt")
        if _targets:
            outfiles.append(f"{outdir}/DNA/common/4_markdup_bam/{dna_id}/{dna_id}.hs_metrics.txt")

    # Somatic mutation spectrum (modules/spectrum) — wired whenever design
    # pairs exist; signature analysis follows Params.somatic_spectrum
    # .signature.enabled (default off), mode selects which outputs register.
    if paired_design:
        _spectrum_fmts = (datajson.get("Params", {}).get("somatic_spectrum", {})
                          .get("plot", {}).get("formats") or ["png"])
        if isinstance(_spectrum_fmts, str):
            _spectrum_fmts = [_spectrum_fmts]
        outfiles += _figs(f"{outdir}/DNA/stability/spectrum/somatic_spectrum_stacked_bar",
                          _spectrum_fmts)
        outfiles += _figs(f"{outdir}/DNA/stability/spectrum/somatic_spectrum_heatmap",
                          _spectrum_fmts)
        outfiles.append(f"{outdir}/DNA/stability/spectrum/somatic_spectrum.csv")
        outfiles.append(f"{outdir}/DNA/stability/spectrum/somatic_spectrum_per_sample.csv")
        outfiles.append(f"{outdir}/DNA/stability/spectrum/somatic_spectrum_sample_stats.csv")
        _sig_cfg = datajson["Params"].get("somatic_spectrum", {}).get("signature", {}) or {}
        if _sig_cfg.get("enabled"):
            _sig_mode = str(_sig_cfg.get("mode", "both")).lower()
            if _sig_mode in ("refit", "both"):
                _sig_fmts = (datajson.get("Params", {}).get("somatic_spectrum", {})
                             .get("signature", {}).get("formats") or ["png"])
                if isinstance(_sig_fmts, str):
                    _sig_fmts = [_sig_fmts]
                for _f in ("signature_refit_exposures.csv",
                           "signature_refit_exposure_fractions.csv",
                           "signature_refit_family_exposures.csv",
                           "signature_refit_summary.csv"):
                    outfiles.append(f"{outdir}/DNA/stability/spectrum/{_f}")
                for _b in ("signature_refit_exposure_stacked_bar",
                           "signature_refit_exposure_heatmap",
                           "signature_refit_signature_cosine"):
                    outfiles += _figs(f"{outdir}/DNA/stability/spectrum/{_b}",
                                      _sig_fmts)
            if _sig_mode in ("nmf", "both"):
                _sig_fmts = (datajson.get("Params", {}).get("somatic_spectrum", {})
                             .get("signature", {}).get("formats") or ["png"])
                if isinstance(_sig_fmts, str):
                    _sig_fmts = [_sig_fmts]
                for _f in ("signature_nmf_signatures.csv",
                           "signature_nmf_exposures.csv",
                           "signature_nmf_rank_diagnostics.csv",
                           "signature_nmf_cosmic_similarity.csv",
                           "signature_nmf_cosmic_best_match.csv"):
                    outfiles.append(f"{outdir}/DNA/stability/spectrum/{_f}")
                for _b in ("signature_nmf_rank_diagnostics",
                           "signature_nmf_exposures_stacked_bar",
                           "signature_nmf_cosmic_match"):
                    outfiles += _figs(f"{outdir}/DNA/stability/spectrum/{_b}",
                                      _sig_fmts)

    # Fragment size (cfDNA) gate: Params.fragment_size.enabled.
    if _enabled(datajson, "fragment_size"):
        outfiles.append(f"{outdir}/DNA/results/fragment_size/fragment/FragmentSize.txt")
        _frag_fmts = (datajson.get("Params", {}).get("fragment_size", {})
                      .get("plotFragmentSize", {}).get("formats") or ["png"])
        if isinstance(_frag_fmts, str):
            _frag_fmts = [_frag_fmts]
        outfiles += _figs(f"{outdir}/DNA/results/fragment_size/fragment/FragmentSize",
                          _frag_fmts)

    # Structural variants (Manta) gate: Params.manta.enabled. Matched-normal
    # design pairs get somatic SV (mirroring the Mutect2 pair judgment);
    # samples outside the pairs get germline SV.
    _pair_members = ({rec["control"] for rec in paired_design}
                     | {rec["experimental"] for rec in paired_design})
    if _enabled(datajson, "manta"):
        for rec in paired_design:
            outfiles.append(f"{outdir}/DNA/variation/somatic_sv/"
                            f"{rec['control']}_vs_{rec['experimental']}/"
                            f"results/variants/somaticSV.vcf.gz")
        for key in dna_samples:
            dna_id = _disp(key)
            if dna_id in _pair_members:
                logger.info(f"Sample {dna_id} is involved in manta somatic analysis, skipping germline SV workflow for this sample.")
                continue
            outfiles.append(f"{outdir}/DNA/variation/germline_sv/{dna_id}/results/variants/candidateSV.vcf.gz")

    # MSI (msisensorpro) gate: Params.msi.enabled, tumor of each design pair.
    if _enabled(datajson, "msi") and paired_design:
        for rec in paired_design:
            outfiles.append(
                f"{outdir}/DNA/stability/msi/{rec['experimental']}/"
                f"{rec['experimental']}_msi_pro.tsv")
        outfiles.append(f"{outdir}/DNA/stability/msi/msi_input.txt")
        # formats 精确到 rule:Params.msi.msi_plot.formats(config 全量路径读取,
        # 无 module 级共享回退)。
        _msi_fmts = (datajson.get("Params", {}).get("msi", {})
                     .get("msi_plot", {}).get("formats") or ["jpg"])
        if isinstance(_msi_fmts, str):
            _msi_fmts = [_msi_fmts]
        outfiles += _figs(f"{outdir}/DNA/stability/msi/msi_barplot", _msi_fmts)

    # CNVkit gate: Params.cnvkit.enabled. The assay-resolved method computed
    # with the targets gating above is applied here (Mutation parity; judged
    # over the DNA rows only — RNA rows share the meta).
    _cnv_cfg = datajson["Params"].setdefault("cnvkit", {})
    if _dna_assays:
        configured = _cnv_cfg.get("method")
        if configured != _dna_assay_method:
            logger.info(
                f"DR: meta assay '{_dna_assays[0]}' -> Params.cnvkit.method "
                f"'{configured}' overridden to '{_dna_assay_method}'"
            )
        _cnv_cfg["method"] = _dna_assay_method
    # Calling mode resolution (Mutation parity): matched normals in the
    # design pairs -> paired; else a normal panel
    # (Params.workflow.control_samples) -> pooled; else flat. An explicit
    # mode wins over "auto"; impossible combos hard-error below.
    _cnv_mode = _cnv_cfg.get("mode") or "auto"
    if _cnv_mode == "auto":
        _cnv_mode = "paired" if _paired_map else (
            "pooled" if (datajson["Params"]["workflow"].get("control_samples") or []) else "flat"
        )
        _cnv_cfg["mode"] = _cnv_mode
        logger.info(f"DR: cnvkit mode auto -> {_cnv_mode} "
                    f"(design pairs: {len(_paired_map)}, "
                    f"control_samples: {len(datajson['Params']['workflow'].get('control_samples') or [])})")
    if _cnv_mode == "paired" and not _paired_map:
        raise ValueError(
            "DR: Params.cnvkit.mode 'paired' requires matched normal "
            "samples in the meta design; none found. Use mode 'pooled' (with "
            "control_samples) or 'flat'."
        )
    if _cnv_mode == "pooled" and not (datajson["Params"]["workflow"].get("control_samples") or []):
        raise ValueError(
            "DR: Params.cnvkit.mode 'pooled' requires "
            "Params.workflow.control_samples; none found. Use mode 'flat' "
            "or provide a normal panel."
        )
    if _cnv_mode == "paired":
        # Paired mode must not also build a pooled reference alongside the
        # per-pair references (conflicting reference semantics).
        datajson["Params"]["workflow"]["control_samples"] = []
    # 参数精确到 rule:resolve 结果(以及 meta assay 推导的 method)镜像写入
    # rule 子段——modules/cnvkit 只读 Params.cnvkit.<rule>.mode/.method
    # (rule_param,无 module 级回退),扁平键仅留作 node.py 自身门控/兼容。
    for _cnv_rule in ("cnvkit_reference", "cnvkit_batch", "cnvkit_paired"):
        _cnv_sub = _cnv_cfg.setdefault(_cnv_rule, {})
        _cnv_sub["mode"] = _cnv_mode
        if _cnv_cfg.get("method"):
            _cnv_sub["method"] = _cnv_cfg["method"]
    if _enabled(datajson, "cnvkit"):
        # Mode-aware registration: paired -> per-pair tumors only (no pooled
        # reference); pooled -> all samples + shared reference.cnn; flat ->
        # all samples, no reference.
        if _cnv_mode == "paired":
            for sample_id in _paired_map:
                outfiles.append(f"{outdir}/DNA/variation/cnv/{sample_id}/{sample_id}.cnr")
                outfiles.append(f"{outdir}/DNA/variation/cnv/{sample_id}/{sample_id}.cns")
        else:
            for sample_id in [_disp(k) for k in dna_samples]:
                outfiles.append(f"{outdir}/DNA/variation/cnv/{sample_id}/{sample_id}.cnr")
                outfiles.append(f"{outdir}/DNA/variation/cnv/{sample_id}/{sample_id}.cns")
            if _cnv_mode == "pooled":
                outfiles.append(f"{outdir}/DNA/variation/cnv/reference/reference.cnn")

    # Mutation burden (modules/burden): SNV/InDel TMB + SV + CNV burden for
    # each tumor sample (Mutation parity). The three sample -> file maps
    # derive from the somatic and CNV outputs above (hand-filled Params.burden
    # entries win); the CNV sample set follows the resolved CNVkit mode.
    if _enabled(datajson, "burden"):
        burden_cfg = datajson["Params"].setdefault("burden", {})
        snv_vcf_dict = dict(sample_somatic_vcf_dict)  # {tumor: somatic SNV VCF}
        sv_vcf_dict = {
            rec["experimental"]:
                f"{outdir}/DNA/variation/somatic_sv/"
                f"{rec['control']}_vs_{rec['experimental']}/"
                f"results/variants/somaticSV.vcf.gz"
            for rec in paired_design
        }
        if _cnv_mode == "paired":
            burden_cnv_samples = list(_paired_map)
        else:
            burden_cnv_samples = [_disp(k) for k in dna_samples]
        cnv_dict = {
            sample_id: f"{outdir}/DNA/variation/cnv/{sample_id}/{sample_id}.cns"
            for sample_id in burden_cnv_samples
        }
        for key, derived in (("sample_snv_vcf_dict", snv_vcf_dict),
                             ("sample_sv_vcf_dict", sv_vcf_dict),
                             ("sample_cnv_dict", cnv_dict)):
            if not burden_cfg.get(key):
                burden_cfg[key] = derived
        # VEP annotation settings (Params.burden.burden_snv_indel.annotate,
        # default on; config 全量路径读取,无 module 级共享回退):
        # species/assembly are derived from the genome; an explicit config
        # value wins and an unmapped genome hard-errors instead of guessing.
        _burden_snvi = burden_cfg.setdefault("burden_snv_indel", {})
        if _burden_snvi.get("annotate", True):
            if not _burden_snvi.get("vep_species") \
                    or not _burden_snvi.get("vep_assembly"):
                _genome_key = str(datajson["genome"]["default"])
                mapped = GENOME_VEP_SPECIES.get(_genome_key)
                if not mapped:
                    raise ValueError(
                        f"DR: genome '{_genome_key}' has no known VEP "
                        f"species/assembly for the burden annotation; set "
                        f"'Params.burden.burden_snv_indel.vep_species' and "
                        f"'Params.burden.burden_snv_indel.vep_assembly' in the config JSON, or "
                        f"disable annotation with Params.burden.burden_snv_indel.annotate=false."
                    )
                _burden_snvi.setdefault("vep_species", mapped[0])
                _burden_snvi.setdefault("vep_assembly", mapped[1])
                logger.info(
                    f"DR: burden VEP settings auto-filled for genome "
                    f"'{_genome_key}': species={mapped[0]} assembly={mapped[1]}"
                )
        for sample_id in burden_cfg["sample_snv_vcf_dict"]:
            outfiles.append(f"{outdir}/DNA/stability/burden/{sample_id}/{sample_id}.snv_indel_burden.tsv")
        for sample_id in burden_cfg["sample_sv_vcf_dict"]:
            outfiles.append(f"{outdir}/DNA/stability/burden/{sample_id}/{sample_id}.sv_burden.tsv")
        for sample_id in burden_cfg["sample_cnv_dict"]:
            outfiles.append(f"{outdir}/DNA/stability/burden/{sample_id}/{sample_id}.cnv_burden.tsv")
        if (burden_cfg["sample_snv_vcf_dict"] or burden_cfg["sample_sv_vcf_dict"]
                or burden_cfg["sample_cnv_dict"]):
            outfiles.append(f"{outdir}/DNA/stability/burden/burden_summary.tsv")

    datajson["Params"].setdefault("somatic_spectrum", {})["sample_somatic_vcf_dict"] = sample_somatic_vcf_dict
    datajson["Params"]["somatic_spectrum"]["sample_group_dict"] = sample_group_dict

    # Sex check (modules/sexcheck): per-sample sex call from sex-chromosome
    # coverage over the markdup BAMs (Mutation parity).
    outfiles.append(f"{outdir}/DNA/QC/sexcheck/sexcheck.tsv")

    # QC summary module (modules/qc_summary): cross-sample metric table +
    # HTML summary; always part of the delivered output set.
    outfiles.append(f"{outdir}/DNA/QC/qc_summary/qc_summary.tsv")
    outfiles.append(f"{outdir}/DNA/QC/qc_summary/qc_summary.html")

    # Final PPT report + results workbook (modules/Mutation_report), gated
    # by Params.report.enabled like the other report modules.
    if _enabled(datajson, "report"):
        # 写回精确到 rule:Params.report.generate_report.date(与
        # runMutation/runPeakCalling/runRIPseq/runRNAseq 的写回同一键形)
        datajson["Params"].setdefault("report", {}).setdefault(
            "generate_report", {}).setdefault(
            "date", time.strftime("%Y-%m-%d", time.localtime()))
        outfiles.append(f"{outdir}/DNA/Mutation_report.pptx")
        outfiles.append(f"{outdir}/DNA/Mutation_report.xlsx")
        # 报告自有图导出清单(报告模块规范:figures/<name>.<fmt> 条件产出,
        # 以 manifest.tsv 锚定 DAG 与磁盘一致)
        outfiles.append(f"{outdir}/DNA/figures/manifest.tsv")

    outfiles.append(f"{outdir}/DRclass/DRclass_summary.tsv")

    datajson["Params"]["workflow"]["dr_pairs"] = dr_pairs
    datajson["Params"]["workflow"]["drclass_samples"] = drclass_samples
    datajson["Params"]["workflow"]["genome_paired_samples"] = genome_paired_samples
    datajson["Params"]["workflow"]["genome_single_samples"] = genome_single_samples
    datajson["Params"]["workflow"]["dna_paired_samples"] = dna_paired_samples
    datajson["Params"]["workflow"]["dna_single_samples"] = dna_single_samples
    datajson["Params"]["workflow"]["paired_design"] = paired_design
    datajson["Params"]["workflow"]["pair_strains"] = pair_strains
    datajson["Params"]["workflow"]["control_samples"] = datajson["Params"]["workflow"].get("control_samples", []) or []
    datajson["outfiles"] = outfiles
    return _write_raw_json(datajson, outdir)

def runPopulation_genomics(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ) -> str:
    """Prepare input JSON for the Population_genomics workflow.

    Population-level variant analysis chain: gatk HaplotypeCaller (per-sample
    gVCF) -> joint genotyping -> bcftools filtering -> plink2 / vcftools /
    admixture / easySFS / PopLDdecay analyses gated by Params.analysis.

    Params.workflow.sample_bams and Params.workflow.populations are
    user-supplied maps; empty ones fall back to the meta (per-sample ``bam``
    column for sample_bams, ``group`` column for populations).
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    analysis = datajson.get("Params", {}).get("analysis", {}) or {}

    wf["samples"] = list(samples_info_dict) or sorted(wf.get("sample_bams") or {})
    # sample_bams: explicit config entries win; empty slots fall back to the
    # meta's per-sample BAM column.
    sample_bams = dict(wf.get("sample_bams") or {})
    for sid, info in samples_info_dict.items():
        bam = getattr(info, "bam", None)
        if bam and not sample_bams.get(sid):
            sample_bams[sid] = bam
    wf["sample_bams"] = sample_bams
    missing = [sid for sid in wf["samples"] if sid not in sample_bams]
    if missing:
        raise ValueError(
            f"Population_genomics: no BAM for sample(s) {missing}; fill "
            f"Params.workflow.sample_bams or record the meta 'bam' column"
        )
    # populations: explicit config map wins; otherwise group samples by the
    # meta's ``group`` column.
    populations = wf.get("populations") or {}
    if not populations:
        by_group: Dict[str, List[str]] = {}
        for sid, info in samples_info_dict.items():
            group = getattr(info, "group", None)
            if group:
                by_group.setdefault(group, []).append(sid)
        populations = by_group
    wf["populations"] = populations

    # outfiles mirror the fallback list in subworkflow/Population_genomics.smk
    pop_names = sorted(populations)
    pop_pairs = list(combinations(pop_names, 2))
    filtered_vcf = f"{outdir}/variants/filtered/population_snps.vcf.gz"
    plink_prefix = f"{outdir}/analysis/plink/population_snps"
    outfiles = [filtered_vcf, f"{filtered_vcf}.tbi"]
    if analysis.get("pca", True):
        outfiles += [f"{plink_prefix}.eigenvec", f"{plink_prefix}.eigenval"]
    if analysis.get("structure", True):
        outfiles += [f"{outdir}/analysis/structure/admixture.K{int(k)}.Q"
                     for k in analysis.get("admixture_k", [2, 3, 4, 5, 6])]
    if analysis.get("diversity", True):
        for population in pop_names:
            outfiles += [f"{outdir}/analysis/diversity/{population}.windowed.pi",
                         f"{outdir}/analysis/diversity/{population}.Tajima.D"]
    if analysis.get("fst", True):
        outfiles += [f"{outdir}/analysis/fst/{a}_vs_{b}.weir.fst"
                     for a, b in pop_pairs]
    if analysis.get("history", False):
        outfiles.append(f"{outdir}/analysis/history/easySFS.done")
    if analysis.get("ld_decay", False):
        outfiles += [f"{outdir}/analysis/ld_decay/{population}.stat.gz"
                     for population in pop_names]
    if analysis.get("gwas", False):
        outfiles.append(f"{outdir}/analysis/gwas/gwas.done")
    datajson["outfiles"] = outfiles
    return _write_raw_json(datajson, outdir)

def runspatial_transcriptomics(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ) -> str:
    """Prepare input JSON for the spatial_transcriptomics workflow.

    Scanpy-based Visium analysis: QC -> clustering -> spatial
    autocorrelation (advanced, gated by
    Params.spatial.spatial_scanpy_advanced.spatial_autocorrelation) over one
    h5ad input (Params.workflow.input_h5ad or Params.workflow.visium_h5).
    """
    root_dir = _init_raw_json(datajson, indir, outdir, raw_files)
    wf = datajson.setdefault("Params", {}).setdefault("workflow", {})
    if not wf.get("input_h5ad") and not wf.get("visium_h5"):
        raise ValueError(
            "spatial_transcriptomics: Params.workflow.input_h5ad or "
            "Params.workflow.visium_h5 is required"
        )
    spatial_cfg = datajson.get("Params", {}).get("spatial", {}) or {}
    advanced_enabled = (spatial_cfg.get("spatial_scanpy_advanced") or {}).get(
        "spatial_autocorrelation", True)

    # outfiles mirror the fallback list in
    # subworkflow/spatial_transcriptomics.smk
    outfiles = [
        f"{outdir}/spatial/qc/filtered.h5ad",
        f"{outdir}/spatial/qc/qc_metrics.tsv",
        f"{outdir}/spatial/cluster/clustered.h5ad",
        f"{outdir}/spatial/cluster/spatial_markers.tsv",
    ]
    if advanced_enabled:
        outfiles.append(
            f"{outdir}/spatial/advanced/spatial_autocorrelation.h5ad")
    datajson["outfiles"] = outfiles
    return _write_raw_json(datajson, outdir)

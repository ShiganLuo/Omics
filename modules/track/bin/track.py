from __future__ import annotations

import argparse
import html
import json
import logging
import os
import re
from typing import Dict, List, Sequence, Tuple, Union
from urllib.parse import urlparse


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def validate_igv_config(config_data: Dict) -> Tuple[bool, List[str]]:
    errors: List[str] = []

    if not isinstance(config_data, dict):
        return False, ["IGV config must be a JSON object/dictionary"]

    for key in ["id", "name", "tracks"]:
        if key not in config_data:
            errors.append(f"Missing required top-level key: {key}")

    has_top_level_fasta = "fastaURL" in config_data and "indexURL" in config_data
    has_reference = isinstance(config_data.get("reference"), dict) and (
        "fastaURL" in config_data["reference"] and "indexURL" in config_data["reference"]
    )
    if not (has_top_level_fasta or has_reference):
        errors.append("Either top-level 'fastaURL'/'indexURL' or a 'reference' object with those keys must be present")

    for key in ["id", "name"]:
        value = config_data.get(key)
        if value is not None and not isinstance(value, str):
            errors.append(f"Top-level key '{key}' must be a string")
        elif isinstance(value, str) and not value.strip():
            errors.append(f"Top-level key '{key}' cannot be empty")

    def _validate_url_string(value, context_name):
        if value is None:
            return
        if not isinstance(value, str) or not value.strip():
            errors.append(f"{context_name} must be a non-empty string")
            return
        parsed = urlparse(value)
        if parsed.scheme and parsed.scheme not in {"http", "https", "file"}:
            errors.append(f"{context_name} has unsupported URL scheme '{parsed.scheme}'")

    if has_top_level_fasta:
        _validate_url_string(config_data.get("fastaURL"), "Top-level 'fastaURL'")
        _validate_url_string(config_data.get("indexURL"), "Top-level 'indexURL'")
    if has_reference:
        reference = config_data.get("reference", {})
        _validate_url_string(reference.get("fastaURL"), "reference.fastaURL")
        _validate_url_string(reference.get("indexURL"), "reference.indexURL")

    tracks = config_data.get("tracks")
    if not isinstance(tracks, list):
        errors.append("Top-level key 'tracks' must be a list")
        tracks = []

    allowed_formats = {"bigwig", "wig", "bedgraph", "gtf", "bed", "bam", "vcf", "junction", "snp"}
    for index, track in enumerate(tracks):
        if not isinstance(track, dict):
            errors.append(f"Track entry at index {index} must be an object/dictionary")
            continue

        for key in ["name", "format", "url"]:
            if key not in track:
                errors.append(f"Track entry at index {index} is missing required key: {key}")

        track_name = track.get("name")
        track_format = track.get("format")
        track_url = track.get("url")

        if track_name is not None and (not isinstance(track_name, str) or not track_name.strip()):
            errors.append(f"Track entry at index {index} has an invalid 'name'")
        if track_format is not None:
            if not isinstance(track_format, str) or not track_format.strip():
                errors.append(f"Track entry at index {index} has an invalid 'format'")
            elif track_format.lower() not in allowed_formats:
                errors.append(f"Track entry at index {index} uses unsupported format '{track_format}'")
        if track_url is not None:
            if not isinstance(track_url, str) or not track_url.strip():
                errors.append(f"Track entry at index {index} has an invalid 'url'")
            else:
                parsed = urlparse(track_url)
                if parsed.scheme and parsed.scheme not in {"http", "https", "file"}:
                    errors.append(f"Track entry at index {index} has unsupported URL scheme '{parsed.scheme}'")

    return len(errors) == 0, errors


def _track_type(file_path: str) -> str:
    suffix = os.path.splitext(file_path)[1].lower()
    if suffix in {".bw", ".bigwig"}:
        return "bigWig"
    if suffix in {".bedgraph", ".bdg"}:
        return "bedGraph"
    return "bigWig"


def _track_name(file_path: str) -> str:
    return os.path.splitext(os.path.basename(file_path))[0]


def _parse_track_identity(track_name: str) -> Dict[str, str]:
    raw_name = track_name.strip()
    tokens = [token for token in re.split(r"[._\-]+", raw_name) if token]

    strand_map = {
        "plus": "plus",
        "p": "plus",
        "pos": "plus",
        "positive": "plus",
        "forward": "plus",
        "fwd": "plus",
        "minus": "minus",
        "m": "minus",
        "neg": "minus",
        "negative": "minus",
        "reverse": "minus",
        "rev": "minus",
    }
    drop_tokens = {"bigwig", "bw", "bedgraph", "track", "signal", "coverage", "sorted", "dedup", "unique", "rmdup", "wig"}

    strand = "unknown"
    cleaned_tokens: List[str] = []
    for token in tokens:
        lowered = token.lower()
        if lowered in strand_map:
            strand = strand_map[lowered]
            continue
        if lowered in drop_tokens:
            continue
        cleaned_tokens.append(token)

    if not cleaned_tokens:
        cleaned_tokens = tokens or [raw_name]

    group_tokens: List[str] = []
    for token in cleaned_tokens:
        lowered = token.lower()
        if re.fullmatch(r"rep\d+", lowered):
            continue
        if re.fullmatch(r"r\d+", lowered):
            continue
        if re.fullmatch(r"lane\d+", lowered):
            continue
        if re.fullmatch(r"l\d+", lowered):
            continue
        if re.fullmatch(r"read\d+", lowered):
            continue
        group_tokens.append(token)

    sample = "_".join(cleaned_tokens) if cleaned_tokens else raw_name
    group = "_".join(group_tokens) if group_tokens else sample
    return {"sample": sample, "group": group, "strand": strand}


def _ensure_parent_dir(file_path: str) -> None:
    parent_dir = os.path.dirname(file_path)
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)


def _load_json_file(file_path: str) -> Dict:
    with open(file_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _as_path_map(path_map: Union[Dict, List]) -> List[Tuple[str, str]]:
    if isinstance(path_map, dict):
        items = list(path_map.items())
    elif isinstance(path_map, list):
        items = []
        for entry in path_map:
            if isinstance(entry, (list, tuple)) and len(entry) == 2:
                items.append((str(entry[0]), str(entry[1])))
    else:
        items = []

    normalized_items = ((str(source), str(target)) for source, target in items)
    return sorted(normalized_items, key=lambda item: len(item[0]), reverse=True)


def _rewrite_resource_url(resource_url: str, output_dir: str, path_map: List[Tuple[str, str]]) -> str:
    if not isinstance(resource_url, str) or not resource_url.strip():
        return resource_url

    parsed = urlparse(resource_url)
    if parsed.scheme in {"http", "https"}:
        return resource_url

    if resource_url.startswith("/"):
        for source_prefix, public_prefix in path_map:
            if resource_url.startswith(source_prefix):
                suffix = resource_url[len(source_prefix):].lstrip("/")
                return f"{public_prefix.rstrip('/')}/{suffix}" if suffix else public_prefix.rstrip("/")

    try:
        return os.path.relpath(resource_url, start=output_dir)
    except Exception:
        return resource_url


def _normalize_browser_config(browser_config: Dict, output_dir: str) -> Dict:
    path_map = _as_path_map(browser_config.get("publicPathMap", []))

    for key in ["fastaURL", "indexURL", "cytobandURL"]:
        if key in browser_config:
            browser_config[key] = _rewrite_resource_url(browser_config[key], output_dir, path_map)

    reference = browser_config.get("reference")
    if not isinstance(reference, dict):
        reference = {}

    for key in ["id", "name", "fastaURL", "indexURL", "cytobandURL"]:
        if key in browser_config and key not in reference:
            reference[key] = browser_config[key]

    if reference:
        browser_config["reference"] = reference
        browser_config["genome"] = reference.get("id") or reference.get("name") or browser_config.get("id")

    for key in ("fastaURL", "indexURL", "cytobandURL"):
        if key in browser_config and isinstance(browser_config.get("reference"), dict):
            try:
                if browser_config.get(key) == browser_config["reference"].get(key):
                    del browser_config[key]
            except Exception:
                pass

    for track in browser_config.get("tracks", []):
        if isinstance(track, Dict) and "url" in track:
            track["url"] = _rewrite_resource_url(track["url"], output_dir, path_map)
        if isinstance(track, Dict) and "indexURL" in track:
            track["indexURL"] = _rewrite_resource_url(track["indexURL"], output_dir, path_map)

    return browser_config


def _track_color(track_name: str, index: int) -> str:
    lowered_name = track_name.lower()
    if "plus" in lowered_name:
        return "#1f77b4"
    if "minus" in lowered_name:
        return "#d62728"
    return "#1f77b4" if index % 2 == 0 else "#d62728"


def _build_igv_track_record(track_file: str, output_dir: str, index: int) -> Dict:
    track_type = _track_type(track_file)
    if track_type != "bigWig":
        raise ValueError(f"IGV.js recommended input type is bigWig: {track_file}")

    track_name = _track_name(track_file)
    return {
        "name": track_name,
        "url": track_file,
        "type": "wig",
        "format": "bigwig",
        "color": _track_color(track_name, index),
        "displayMode": "COLLAPSED",
        "autoscale": True,
        "height": 25,
    }


def _build_igv_tracks(track_files: Sequence[str], output_dir: str) -> List[Dict]:
    return [_build_igv_track_record(track_file, output_dir, index) for index, track_file in enumerate(track_files)]


def _build_grouped_track_manifest(tracks: Sequence[Dict]) -> List[Dict]:
    grouped: Dict[str, Dict] = {}
    for index, track in enumerate(tracks):
        track_name = str(track.get("name", f"track_{index + 1}"))
        identity = _parse_track_identity(track_name)
        group_name = identity["group"]
        sample_name = identity["sample"]
        strand = identity["strand"]

        group_entry = grouped.setdefault(group_name, {"group": group_name, "samples": {}})
        sample_entry = group_entry["samples"].setdefault(sample_name, {"sample": sample_name, "tracks": []})

        sample_entry["tracks"].append(
            {
                "key": f"generated_track_{index + 1}",
                "name": track_name,
                "strand": strand,
                "config": track,
            }
        )

    grouped_list: List[Dict] = []
    for group_name in sorted(grouped):
        sample_list: List[Dict] = []
        samples = grouped[group_name]["samples"]
        for sample_name in sorted(samples):
            sample_tracks = sorted(samples[sample_name]["tracks"], key=lambda item: item["name"].lower())
            sample_list.append({"sample": sample_name, "tracks": sample_tracks})
        grouped_list.append({"group": group_name, "samples": sample_list})
    return grouped_list


def _split_eager_and_deferred_tracks(tracks: Sequence[Dict]) -> Tuple[List[Dict], List[Dict]]:
    deferred_formats = {"gtf", "gff", "gff3"}
    eager: List[Dict] = []
    deferred: List[Dict] = []

    for track in tracks:
        if not isinstance(track, dict):
            continue
        track_format = str(track.get("format", "")).lower()
        if track_format in deferred_formats:
            deferred.append(track)
        else:
            eager.append(track)

    return eager, deferred


def _build_reference_track_group(tracks: Sequence[Dict]) -> List[Dict]:
    if not tracks:
        return []

    reference_items: List[Dict] = []
    for index, track in enumerate(tracks):
        track_name = str(track.get("name", f"reference_track_{index + 1}"))
        reference_items.append(
            {
                "key": f"reference_track_{index + 1}",
                "name": track_name,
                "strand": "annotation",
                "config": track,
            }
        )

    reference_items = sorted(reference_items, key=lambda item: item["name"].lower())
    return [
        {
            "group": "Reference_annotations",
            "samples": [
                {
                    "sample": "reference",
                    "tracks": reference_items,
                }
            ],
        }
    ]


def _render_igv_html(template_path: str, title: str, igv_js: str, browser_config: Dict, grouped_manifest: List[Dict]) -> str:
    browser_json = json.dumps(browser_config, ensure_ascii=False, separators=(",", ":"))
    grouped_json = json.dumps(grouped_manifest, ensure_ascii=False, separators=(",", ":"))

    if os.path.exists(template_path):
        with open(template_path, "r", encoding="utf-8") as handle:
            template = handle.read()
        return (
            template.replace("{{TITLE}}", html.escape(title))
            .replace("{{IGV_JS}}", igv_js)
            .replace("{{BROWSER_CONFIG}}", browser_json)
            .replace("{{GROUPED_MANIFEST}}", grouped_json)
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>{html.escape(title)}</title>
<script src="{igv_js}"></script>
</head>
<body>
<div id="igv-browser" style="height:100vh"></div>
<script>
const options = {browser_json};
const groupedManifest = {grouped_json};
console.log('IGV fallback template loaded', options, groupedManifest);
</script>
</body>
</html>
"""


def ucsc_track_format(track_files: Sequence[str], output: str) -> List[str]:
    output_dir = os.path.dirname(output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    track_lines = []
    for track_file in track_files:
        track_type = _track_type(track_file)
        track_name = _track_name(track_file)
        track_line = f'track type={track_type} name="{track_name}" description="{track_name}" bigDataUrl={track_file}'
        track_lines.append(track_line)

    track_lines = sorted(track_lines)
    with open(output, "w", encoding="utf-8") as handle:
        handle.write("\n".join(track_lines) + ("\n" if track_lines else ""))
    return track_lines


def igvjs_track_format(
    track_files: Sequence[str],
    output: str,
    config: str,
    title: str = "IGV Browser",
) -> Dict:
    output_dir = os.path.dirname(output)
    _ensure_parent_dir(output)
    track_files = sorted(track_files)

    browser_config = _load_json_file(config)
    path_map = _as_path_map(browser_config.get("publicPathMap", []))

    igv_js = browser_config.get("js", "https://cdn.jsdelivr.net/npm/igv/dist/igv.min.js")
    igv_js = _rewrite_resource_url(igv_js, output_dir, path_map)

    generated_tracks = _build_igv_tracks(track_files, output_dir)
    for track in generated_tracks:
        if "url" in track:
            track["url"] = _rewrite_resource_url(track["url"], output_dir or ".", path_map)

    browser_config["tracks"] = browser_config.get("tracks", [])
    browser_config = _normalize_browser_config(browser_config, output_dir or ".")

    normalized_base_tracks = browser_config.get("tracks", [])
    eager_base_tracks, deferred_base_tracks = _split_eager_and_deferred_tracks(normalized_base_tracks)
    # Load both eager and deferred (e.g. GTF) base tracks into the browser by default
    # so reference/annotation tracks are visible on initial load. Assign a high
    # order for reference tracks so they appear below sample tracks.
    REFERENCE_TRACK_ORDER_BASE = 1000000
    for idx, t in enumerate(deferred_base_tracks):
        if isinstance(t, dict) and t.get("order") is None:
            try:
                t["order"] = REFERENCE_TRACK_ORDER_BASE + idx
            except Exception:
                pass

    browser_config["tracks"] = eager_base_tracks + deferred_base_tracks

    grouped_manifest = _build_grouped_track_manifest(generated_tracks) + _build_reference_track_group(deferred_base_tracks)

    is_valid, errors = validate_igv_config(browser_config)
    if not is_valid:
        logging.warning("IGV config validation failed:")
        for error in errors:
            logging.warning("  - %s", error)
        raise ValueError("Invalid IGV config, see log for details")

    template_path = os.path.join(os.path.dirname(__file__), "igv_template.html")
    html_text = _render_igv_html(template_path, title, igv_js, browser_config, grouped_manifest)

    with open(output, "w", encoding="utf-8") as handle:
        handle.write(html_text)

    return browser_config


def build_track_files(inputs: Sequence[str], output: str, mode: str, config: str = "") -> None:
    if mode == "ucsc":
        ucsc_track_format(inputs, output)
    elif mode == "igv":
        if not config:
            raise ValueError("IGV mode requires a config JSON file")
        igvjs_track_format(inputs, output, config=config)
    else:
        raise ValueError(f"Unsupported mode: {mode}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate UCSC track files or IGV.js HTML from coverage tracks.")
    parser.add_argument("-i", "--input", required=True, nargs="+", help="Input track files.")
    parser.add_argument("-o", "--output", required=True, help="Output file path.")
    parser.add_argument("--mode", choices=["ucsc", "igv"], default="ucsc", help="Output format to generate.")
    parser.add_argument(
        "--config",
        default="",
        help="IGV config JSON file containing genome and other settings (required for igv mode).",
    )
    args = parser.parse_args()

    logging.info("Generating %s output from %d tracks.", args.mode, len(args.input))
    build_track_files(args.input, args.output, args.mode, args.config)
    logging.info("Conversion complete.")


if __name__ == "__main__":
    main()

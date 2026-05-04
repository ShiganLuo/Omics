from __future__ import annotations

import argparse
import html
import json
import logging
import os
from typing import Dict, List, Sequence, Tuple, Union
from urllib.parse import urlparse


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def validate_igv_config(config_data: Dict) -> Tuple[bool, List[str]]:
    """Validate an IGV.js browser config dictionary.

    Parameters
    ----------
    config_data : dict
        IGV.js browser configuration loaded from JSON.

    Returns
    -------
    tuple[bool, list[str]]
        A tuple containing the validation result and a list of error messages.
    """
    errors: List[str] = []

    if not isinstance(config_data, dict):
        return False, ["IGV config must be a JSON object/dictionary"]

    # Require id, name, and tracks. Accept either top-level fasta/index or a
    # `reference` object that provides those URLs.
    for key in ["id", "name", "tracks"]:
        if key not in config_data:
            errors.append(f"Missing required top-level key: {key}")

    # Accept reference object as alternative to top-level fasta/index
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

    # Validate fasta/index presence in either top-level or reference and ensure strings
    def _validate_url_string(v, context_name):
        if v is None:
            return
        if not isinstance(v, str) or not v.strip():
            errors.append(f"{context_name} must be a non-empty string")
        else:
            parsed = urlparse(v)
            if parsed.scheme and parsed.scheme not in {"http", "https", "file"}:
                errors.append(f"{context_name} has unsupported URL scheme '{parsed.scheme}'")

    if has_top_level_fasta:
        _validate_url_string(config_data.get("fastaURL"), "Top-level 'fastaURL'")
        _validate_url_string(config_data.get("indexURL"), "Top-level 'indexURL'")
    if has_reference:
        ref = config_data.get("reference", {})
        _validate_url_string(ref.get("fastaURL"), "reference.fastaURL")
        _validate_url_string(ref.get("indexURL"), "reference.indexURL")

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
                    errors.append(
                        f"Track entry at index {index} has unsupported URL scheme '{parsed.scheme}'"
                    )

    return len(errors) == 0, errors


def _track_type(file_path: str) -> str:
    """Infer the IGV track type from a file extension.

    Parameters
    ----------
    file_path : str
        Track file path.

    Returns
    -------
    str
        IGV track type string.
    """
    suffix = os.path.splitext(file_path)[1].lower()
    if suffix in {".bw", ".bigwig"}:
        return "bigWig"
    if suffix in {".bedgraph", ".bdg"}:
        return "bedGraph"
    return "bigWig"


def _track_name(file_path: str) -> str:
    """Derive a display name from a track file path.

    Parameters
    ----------
    file_path : str
        Track file path.

    Returns
    -------
    str
        Base name without the final extension.
    """
    return os.path.basename(file_path).split(".")[0]


def _ensure_parent_dir(file_path: str) -> None:
    """Create the parent directory for a file path when needed.

    Parameters
    ----------
    file_path : str
        Target file path.
    """
    parent_dir = os.path.dirname(file_path)
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)


def _load_json_file(file_path: str) -> Dict:
    """Load a JSON file into a dictionary.

    Parameters
    ----------
    file_path : str
        JSON file path.

    Returns
    -------
    dict
        Parsed JSON object.
    """
    with open(file_path, "r") as f:
        return json.load(f)


def _as_path_map(path_map: Union[Dict, List]) -> List[Tuple[str, str]]:
    """Normalize a public path mapping into a sorted list of prefix pairs.

    Parameters
    ----------
    path_map : dict or list
        Mapping from local prefixes to public URL prefixes.

    Returns
    -------
    list[tuple[str, str]]
        Path prefix pairs sorted from longest source prefix to shortest.
    """
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
    """Rewrite a resource URL to a public or relative path.

    Parameters
    ----------
    resource_url : str
        Original resource URL or file path.
    output_dir : str
        Directory containing the generated HTML output.
    path_map : list[tuple[str, str]]
        Public path mappings for absolute local prefixes.

    Returns
    -------
    str
        Rewritten URL suitable for use in the generated HTML.
    """
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
    """Rewrite IGV browser config URLs to match the deployment layout.

    Parameters
    ----------
    browser_config : dict
        IGV browser configuration.
    output_dir : str
        Directory containing the generated HTML output.

    Returns
    -------
    dict
        Updated browser configuration.
    """
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

    # Remove duplicate top-level fasta/index/cytoband keys once moved into `reference`.
    # IGV prefers a single `reference` object; keeping both forms is redundant and
    # may confuse downstream validators or consumers.
    for _key in ("fastaURL", "indexURL", "cytobandURL"):
        if _key in browser_config and isinstance(browser_config.get("reference"), Dict):
            # Only delete if the same value exists in reference (safe no-op otherwise)
            try:
                if browser_config.get(_key) == browser_config["reference"].get(_key):
                    del browser_config[_key]
            except Exception:
                # If anything unexpected happens, leave the top-level key intact.
                pass

    for track in browser_config.get("tracks", []):
        if isinstance(track, dict) and "url" in track:
            track["url"] = _rewrite_resource_url(track["url"], output_dir, path_map)

    return browser_config


def _track_color(track_name: str, index: int) -> str:
    """Select a display color for a track.

    Parameters
    ----------
    track_name : str
        Track name.
    index : int
        Track index in the generated list.

    Returns
    -------
    str
        Hex color string.
    """
    lowered_name = track_name.lower()
    if "plus" in lowered_name:
        return "#1f77b4"
    if "minus" in lowered_name:
        return "#d62728"
    return "#1f77b4" if index % 2 == 0 else "#d62728"


def _build_igv_track_record(track_file: str, output_dir: str, index: int) -> Dict:
    """Build a single IGV track record from an input file.

    Parameters
    ----------
    track_file : str
        Input track file path.
    output_dir : str
        Directory containing the generated HTML output.
    index : int
        Track index in the generated list.

    Returns
    -------
    dict
        IGV track configuration entry.

    Raises
    ------
    ValueError
        If the input file is not a supported bigWig track.
    """
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
        "autoscale": True,
        "height": 50,
    }


def _build_igv_tracks(track_files: Sequence[str], output_dir: str) -> List[Dict]:
    """Build IGV track records for all input files.

    Parameters
    ----------
    track_files : sequence of str
        Input track files.
    output_dir : str
        Directory containing the generated HTML output.

    Returns
    -------
    list[dict]
        IGV track records.
    """
    return [_build_igv_track_record(track_file, output_dir, index) for index, track_file in enumerate(track_files)]


def ucsc_track_format(track_files: Sequence[str], output: str) -> List[str]:
    """Write a UCSC track list for the provided files.

    Parameters
    ----------
    track_files : sequence of str
        Input track files.
    output : str
        Output text file path.

    Returns
    -------
    list[str]
        Sorted UCSC track lines.
    """
    output_dir = os.path.dirname(output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    track_lines = []
    for track_file in track_files:
        track_type = _track_type(track_file)
        track_name = _track_name(track_file)
        relative_url = os.path.relpath(track_file, start=output_dir or ".")
        track_line = f'track type={track_type} name="{track_name}" description="{track_name}" bigDataUrl={relative_url}'
        track_lines.append(track_line)

    track_lines = sorted(track_lines)
    with open(output, "w") as f:
        f.write("\n".join(track_lines) + ("\n" if track_lines else ""))
    return track_lines


def igvjs_track_format(
    track_files: Sequence[str],
    output: str,
    config: str,
    title: str = "IGV Browser",
) -> Dict:
    """Generate an IGV.js HTML page from track files and a browser config.

    Parameters
    ----------
    track_files : sequence of str
        Input track files.
    output : str
        Output HTML file path.
    config : str
        JSON file containing the base IGV browser configuration.
    title : str, optional
        Title shown in the HTML page.

    Returns
    -------
    dict
        Final IGV browser configuration used to render the page.

    Raises
    ------
    ValueError
        If the resulting IGV configuration is invalid.
    """
    output_dir = os.path.dirname(output)
    _ensure_parent_dir(output)
    track_files = sorted(track_files)
    browser_config = _load_json_file(config)
    igv_js = browser_config.get("js", "https://cdn.jsdelivr.net/npm/igv/dist/igv.min.js")
    browser_config["tracks"] = browser_config.get("tracks", []) + _build_igv_tracks(track_files, output_dir)
    browser_config = _normalize_browser_config(browser_config, output_dir or ".")

    is_valid, errors = validate_igv_config(browser_config)
    if not is_valid:
        logging.warning("IGV config validation failed:")
        for error in errors:
            logging.warning("  - %s", error)
        raise ValueError("Invalid IGV config, see log for details")

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>{html.escape(title)}</title>
<script src="{igv_js}"></script>
<style>
html, body {{
  margin: 0;
  padding: 0;
  height: 100%;
    background: #ffffff;
    color: #111111;
  font-family: Arial, sans-serif;
}}
header {{
  padding: 16px 20px 8px;
    border-bottom: 1px solid #ddd;
}}
h1 {{
  margin: 0;
  font-size: 18px;
  font-weight: 600;
}}
#igv-browser {{
  height: calc(100vh - 58px);
}}
</style>
</head>
<body>
<div id="igv-browser"></div>

<script>
const options = {json.dumps(browser_config, ensure_ascii=False, indent=2)};
igv.createBrowser(document.getElementById('igv-browser'), options)
  .then(browser => console.log("IGV loaded"));
</script>
</body>
</html>
"""

    with open(output, "w", encoding="utf-8") as f:
        f.write(html_text)

    return browser_config


def build_track_files(inputs: Sequence[str], output: str, mode: str, config: str = "") -> None:
    """Dispatch track generation to the requested output mode.

    Parameters
    ----------
    inputs : sequence of str
        Input track files.
    output : str
        Output file path.
    mode : str
        Output mode, either ``ucsc`` or ``igv``.
    config : str
        IGV browser JSON config file path.

    Raises
    ------
    ValueError
        If ``mode`` is unsupported.
    """
    if mode == "ucsc":
        ucsc_track_format(inputs, output)
    elif mode == "igv":
        if not config:
            raise ValueError("IGV mode requires a config JSON file")
        igvjs_track_format(inputs, output, config=config)
    else:
        raise ValueError(f"Unsupported mode: {mode}")


def main() -> None:
    """Run the command-line interface for track generation.

    Returns
    -------
    None
    """
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
"""Shared plotting helpers for the telomere module (multi-format output).

Self-contained helper (modules.md 脚本撰写规范 #0): lives next to the bin
scripts so the module directory runs when copied wholesale. Provides the
``PlotFormat`` type, the matplotlib format set, and a multi-format savefig
helper per 规范 #9 (default ``["png"]``, CLI ``-f/--format`` append).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")

PlotFormat = str

VALID_FORMATS = {
    "png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg", "pgf",
    "raw", "rgba",
}

_RASTER_FORMATS = {"png", "tif", "tiff", "jpg", "jpeg", "raw", "rgba"}


def normalize_formats(formats: Optional[Sequence[str]]) -> List[str]:
    """Return a validated format list; default ``["png"]`` when empty/None.

    Raises ValueError listing the valid set for any unknown format (fast fail).
    """
    if not formats:
        return ["png"]
    out: List[str] = []
    for f in formats:
        fl = str(f).lower()
        if fl not in VALID_FORMATS:
            raise ValueError(
                f"Unknown image format '{f}'. Valid formats: {sorted(VALID_FORMATS)}")
        if fl not in out:
            out.append(fl)
    return out


def save_fig(fig, out_base: str, formats: Optional[Sequence[str]] = None,
             dpi: int = 150) -> List[str]:
    """Save ``fig`` to ``out_base.<fmt>`` for each format; return written paths.

    dpi is applied only to raster formats; vector formats ignore it (规范 #9).
    """
    fmts = normalize_formats(formats)
    written: List[str] = []
    for fmt in fmts:
        path = f"{out_base}.{fmt}"
        save_kwargs: Dict[str, Any] = {"bbox_inches": "tight"}
        if fmt in _RASTER_FORMATS:
            save_kwargs["dpi"] = dpi
        fig.savefig(path, format=fmt, **save_kwargs)
        written.append(path)
    return written

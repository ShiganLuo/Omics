#!/usr/bin/env python
"""De novo mutational-signature extraction by NMF on SBS96 spectra.

Consumes the per-sample count matrix written by spectrum.py (96 contexts x
samples) and factorizes it as

    M (96 x N) ~= W (96 x K) * H (K x N)

with W, H >= 0 under the KL divergence (Poisson likelihood; Brunet et al.
2004 multiplicative updates). For each candidate rank the extraction is
repeated --restarts times with random initializations; the restart with the
lowest reconstruction error is kept. Rank quality is summarized by the
reconstruction error, the cophenetic correlation coefficient of the sample
consensus clustering (Brunet 2004: high and flat before the true rank, drops
after it), and the consensus stability. The chosen rank is reported, not
silently applied -- review the diagnostics before interpreting the signatures.

Extracted signatures can be annotated against a COSMIC matrix via cosine
similarity (--signatures); cross-species matches (human COSMIC on mouse data
or vice versa) are informative but must be reported as cross-species
reference, not homologous assignment.

Outputs (<prefix> = --output-prefix):
  <prefix>_signatures.csv             96 x K extracted signatures (fractions)
  <prefix>_exposures.csv              samples x K exposures (counts)
  <prefix>_rank_diagnostics.csv       per-rank error / cophenetic / stability
  <prefix>_rank_diagnostics.png       diagnostics curves
  <prefix>_signatures.png             96-channel profiles of extracted signatures
  <prefix>_exposures_stacked_bar.png  per-sample exposure stacked bar
  <prefix>_cosmic_similarity.csv      extracted x COSMIC cosine matrix (opt.)
  <prefix>_cosmic_best_match.csv      best COSMIC match per signature (opt.)

Usage:
    python signature_nmf.py \
        --counts somatic_spectrum_per_sample.csv \
        --output-prefix signature_nmf \
        --max-rank 8 --restarts 20 --seed 42 \
        --signatures COSMIC_v3.6_SBS_GRCh38.txt
"""

import argparse
import logging
import sys
from typing import Dict, List, Optional, TextIO, Tuple

import numpy as np

# setup_logger copied verbatim from src/common/util/LogUtil.py (module bin
# scripts are self-contained; see modules.md 脚本撰写规范 #0).
def setup_logger(
    logger_name: str,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    stream: Optional[TextIO] = None,
) -> logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


logger = setup_logger("SignatureNMF", stream=sys.stdout)


# ---------------------------------------------------------------------------
# NMF core (KL divergence, multiplicative updates; Brunet et al. 2004)
# ---------------------------------------------------------------------------

def kl_divergence(m: np.ndarray, wh: np.ndarray) -> float:
    """Generalized KL divergence D(m || wh) = sum(m*log(m/wh) - m + wh).

    Zero entries of m contribute only wh (the 0*log0 convention).
    """
    eps = 1e-10
    wh = np.maximum(wh, eps)
    mask = m > 0
    return float(
        np.sum(m[mask] * np.log(m[mask] / wh[mask]) - m[mask] + wh[mask])
        + np.sum(wh[~mask])
    )


def nmf_kl(m: np.ndarray, k: int, rng: np.random.Generator,
           max_iter: int = 500, tol: float = 1e-6) -> Tuple[np.ndarray, np.ndarray, float]:
    """Factorize m (96 x N) ~= W * H under KL divergence.

    Parameters
    ----------
    m : numpy.ndarray
        Non-negative count matrix (96 x N).
    k : int
        Rank (number of signatures).
    rng : numpy.random.Generator
        Random source for the initialization.
    max_iter : int
        Maximum multiplicative-update iterations.
    tol : float
        Relative KL improvement below which iteration stops.

    Returns
    -------
    (W, H, kl) : signatures (96 x K, columns sum to 1), exposures (K x N),
    final KL divergence.
    """
    eps = 1e-10
    n_features, n_samples = m.shape
    w = rng.random((n_features, k)) + 0.1
    h = rng.random((k, n_samples)) + 0.1
    w /= w.sum(axis=0, keepdims=True)
    prev = np.inf
    for _ in range(max_iter):
        wh = np.maximum(w @ h, eps)
        # H update: H *= (W^T (M / WH)) / (W^T 1)
        h *= (w.T @ (m / wh)) / np.maximum(w.sum(axis=0)[:, None], eps)
        wh = np.maximum(w @ h, eps)
        # W update: W *= ((M / WH) H^T) / (1 H^T)
        w *= ((m / wh) @ h.T) / np.maximum(h.sum(axis=1)[None, :], eps)
        w /= np.maximum(w.sum(axis=0, keepdims=True), eps)
        cur = kl_divergence(m, w @ h)
        if prev < np.inf and abs(prev - cur) / max(prev, eps) < tol:
            break
        prev = cur
    return w, h, kl_divergence(m, w @ h)


def sample_assignments(h: np.ndarray) -> np.ndarray:
    """Assign each sample to the component with its largest exposure."""
    return np.argmax(h, axis=0)


def cophenetic_coefficient(h: np.ndarray) -> float:
    """Cophenetic correlation of the sample consensus clustering.

    Samples are clustered by largest-exposure assignment; the cophenetic
    coefficient of the resulting connectivity matrix is computed with average
    linkage. High values indicate a clean, stable cluster structure.
    """
    from scipy.cluster.hierarchy import cophenet, linkage
    from scipy.spatial.distance import pdist, squareform

    labels = sample_assignments(h)
    n = labels.shape[0]
    if n < 3 or len(set(labels.tolist())) < 2:
        return float("nan")
    connectivity = (labels[:, None] == labels[None, :]).astype(float)
    distance = 1.0 - connectivity
    condensed = pdist(distance)
    if not condensed.any():
        return float("nan")
    z = linkage(condensed, method="average")
    # scipy versions differ in cophenet's return shape ((distances, corrcoef)
    # vs distances only); pick whichever element matches the condensed size.
    res = cophenet(z, condensed)
    candidates = res if isinstance(res, tuple) else (res,)
    coph = None
    for cand in candidates:
        arr = np.asarray(cand, dtype=float).ravel()
        if arr.size == condensed.size:
            coph = arr
            break
    if coph is None:
        return float("nan")
    return float(np.corrcoef(condensed, coph)[0, 1])


def consensus_stability(m: np.ndarray, k: int, rng: np.random.Generator,
                        restarts: int, max_iter: int) -> float:
    """Mean pairwise cosine agreement of signatures across independent restarts.

    Restarts are matched to the first restart's signatures by greedy cosine
    matching; the mean matched cosine is the stability score (1.0 = every
    restart recovered the same signatures).
    """
    if restarts < 2:
        return float("nan")
    runs = []
    for _ in range(min(restarts, 10)):
        w, _, _ = nmf_kl(m, k, rng, max_iter=max_iter)
        runs.append(w)
    ref = runs[0]
    scores = []
    for w in runs[1:]:
        sim = cosine_matrix(ref.T, w.T)  # k x k
        matched = sim.max(axis=1)
        scores.append(float(matched.mean()))
    return float(np.mean(scores)) if scores else float("nan")


def cosine_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pairwise cosine similarity between row vectors of a and b."""
    na = np.linalg.norm(a, axis=1, keepdims=True)
    nb = np.linalg.norm(b, axis=1, keepdims=True)
    na = np.maximum(na, 1e-12)
    nb = np.maximum(nb, 1e-12)
    return (a / na) @ (b / nb).T


# ---------------------------------------------------------------------------
# Analysis driver
# ---------------------------------------------------------------------------

def run_nmf(counts, args, sig_ref=None) -> None:
    """Extract signatures across ranks, pick the best, and write outputs."""
    import pandas as pd
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m = counts.to_numpy(dtype=float)
    contexts = counts.index
    samples = list(counts.columns)
    max_rank = min(args.max_rank, m.shape[0], m.shape[1])
    logger.info(f"Matrix {m.shape[0]} x {m.shape[1]}; scanning ranks 1..{max_rank}")

    rows = []
    best_per_rank: Dict[int, Tuple[np.ndarray, np.ndarray, float]] = {}
    for k in range(1, max_rank + 1):
        best = None
        for r in range(args.restarts):
            rng = np.random.default_rng(args.seed + 1000 * k + r)
            w, h, kl = nmf_kl(m, k, rng, max_iter=args.max_iter)
            if best is None or kl < best[2]:
                best = (w, h, kl)
        assert best is not None
        best_per_rank[k] = best
        w, h, kl = best
        rng_stab = np.random.default_rng(args.seed + 99991 + k)
        stab = consensus_stability(m, k, rng_stab, args.restarts, args.max_iter)
        coph = cophenetic_coefficient(h)
        rows.append([k, kl, coph, stab])
        logger.info(f"rank {k}: KL={kl:.1f} cophenetic={coph:.4f} stability={stab:.4f}")

    diag = pd.DataFrame(rows, columns=["rank", "kl_divergence", "cophenetic", "stability"])

    # Choose rank: the KL elbow (first K where the improvement from K to K+1
    # collapses relative to the improvement that got K here), restricted to
    # ranks meeting the cophenetic floor. Fall back to best stability, then
    # rank 1. The choice is always reported, never silent.
    chosen = args.rank
    if chosen is None:
        kl = diag["kl_divergence"].to_numpy()
        coph = diag["cophenetic"].to_numpy()
        chosen = None
        for idx in range(1, len(kl) - 1):          # K = idx+1, needs K+1 to exist
            gain_here = kl[idx - 1] - kl[idx]
            gain_next = kl[idx] - kl[idx + 1]
            if gain_here > 0 and gain_next < args.elbow_ratio * gain_here \
                    and not np.isnan(coph[idx]) and coph[idx] >= args.min_cophenetic:
                chosen = idx + 1
                break
        if chosen is None:
            coph_ok = diag[diag["cophenetic"] >= args.min_cophenetic]
            if len(coph_ok) and bool(coph_ok["stability"].notna().any()):
                chosen = int(coph_ok.loc[coph_ok["stability"].idxmax(), "rank"])
            elif bool(diag["stability"].notna().any()):
                chosen = int(diag.loc[diag["stability"].idxmax(), "rank"])
            else:
                chosen = 1
    logger.info(f"Chosen rank: {chosen} (override with --rank)")
    diag["selected"] = diag["rank"] == chosen
    diag.to_csv(f"{args.output_prefix}_rank_diagnostics.csv", index=False)

    w, h, kl = best_per_rank[chosen]
    order = np.argsort(-h.sum(axis=1))
    w = w[:, order]
    h = h[order, :]
    names = [f"Sig{i + 1}" for i in range(chosen)]

    sig_df = pd.DataFrame(w, index=contexts, columns=names)
    sig_df.to_csv(f"{args.output_prefix}_signatures.csv")
    exp_df = pd.DataFrame(h.T, index=samples, columns=names)
    exp_df.to_csv(f"{args.output_prefix}_exposures.csv")

    # Rank diagnostics plot.
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(diag["rank"], diag["kl_divergence"], "o-", color="tab:blue", label="KL divergence")
    ax1.set_xlabel("Rank K")
    ax1.set_ylabel("KL divergence", color="tab:blue")
    ax2 = ax1.twinx()
    ax2.plot(diag["rank"], diag["cophenetic"], "s-", color="tab:orange", label="cophenetic")
    ax2.plot(diag["rank"], diag["stability"], "^-", color="tab:green", label="stability")
    ax2.set_ylabel("Coefficient")
    ax2.set_ylim(0, 1.05)
    ax1.axvline(chosen, color="red", ls="--", lw=1, label=f"chosen K={chosen}")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [l.get_label() for l in lines], fontsize=8, loc="center right")
    plt.title("NMF rank diagnostics")
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_rank_diagnostics.png", dpi=300)
    plt.close()

    # Signature profiles: 96-channel bars colored by substitution class.
    class_colors = {
        "C>A": "#4DBBD5", "C>G": "#000000", "C>T": "#E64B35",
        "T>A": "#7E9AC8", "T>C": "#4DAF4A", "T>G": "#CC79A7",
    }
    fig, axes = plt.subplots(chosen, 1, figsize=(12, 2.6 * chosen), squeeze=False)
    x = np.arange(96)
    for idx, name in enumerate(names):
        ax = axes[idx][0]
        vals = sig_df[name].to_numpy()
        colors = [class_colors[label.split("[")[1].split("]")[0]] for label in contexts]
        ax.bar(x, vals, color=colors, width=0.85)
        ax.set_ylabel("Fraction")
        ax.set_title(f"{name} (KL rank {chosen})", fontsize=9)
        ax.set_xlim(-1, 96)
    axes[-1][0].set_xticks(x)
    axes[-1][0].set_xticklabels(contexts, rotation=90, fontsize=5)
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_signatures.png", dpi=300)
    plt.close()

    # Exposure stacked bar per sample.
    ax = exp_df.plot(kind="bar", stacked=True, figsize=(max(6, 1.2 * len(samples)), 6))
    ax.set_ylabel("Exposure (mutations)")
    ax.set_xlabel("Sample")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_exposures_stacked_bar.png", dpi=300)
    plt.close()

    # COSMIC annotation (optional).
    if sig_ref is not None:
        sim = cosine_matrix(sig_df.T.to_numpy(), sig_ref.T.to_numpy())
        sim_df = pd.DataFrame(sim, index=names, columns=list(sig_ref.columns))
        sim_df.to_csv(f"{args.output_prefix}_cosmic_similarity.csv")
        best_rows = []
        for name in names:
            order_idx = np.argsort(-sim_df.loc[name].to_numpy())
            i1, i2 = int(order_idx[0]), int(order_idx[1]) if sim.shape[1] > 1 else None
            best_rows.append([
                name,
                sim_df.columns[i1], float(sim_df.loc[name, sim_df.columns[i1]]),
                sim_df.columns[i2] if i2 is not None else "",
                float(sim_df.loc[name, sim_df.columns[i2]]) if i2 is not None else "",
            ])
            logger.info(
                f"{name}: best match {sim_df.columns[i1]} "
                f"(cosine {sim_df.loc[name, sim_df.columns[i1]]:.3f})"
            )
        pd.DataFrame(best_rows, columns=[
            "extracted_signature", "best_match", "cosine",
            "second_match", "second_cosine",
        ]).to_csv(f"{args.output_prefix}_cosmic_best_match.csv", index=False)
        logger.info(
            f"Wrote {args.output_prefix}_cosmic_similarity.csv, "
            f"{args.output_prefix}_cosmic_best_match.csv"
        )

    logger.info(
        f"Wrote {args.output_prefix}_signatures.csv, "
        f"{args.output_prefix}_exposures.csv, "
        f"{args.output_prefix}_rank_diagnostics.csv, "
        f"{args.output_prefix}_rank_diagnostics.png, "
        f"{args.output_prefix}_signatures.png, "
        f"{args.output_prefix}_exposures_stacked_bar.png"
    )


def load_signatures(path: str, contexts):
    """Load a COSMIC signature matrix and align rows to *contexts*."""
    import pandas as pd

    sig = pd.read_csv(path, sep="\t", index_col=0)
    if set(sig.index) != set(contexts):
        raise ValueError(
            f"{path}: context labels do not match the count matrix "
            f"({len(set(sig.index) ^ set(contexts))} labels differ)."
        )
    return sig.loc[contexts]


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="De novo mutational-signature extraction by NMF")
    parser.add_argument("-c", "--counts", required=True,
                        help="Per-sample count matrix CSV from spectrum.py "
                             "(96 contexts x samples)")
    parser.add_argument("-o", "--output-prefix", required=True,
                        help="Output prefix")
    parser.add_argument("-k", "--max-rank", type=int, default=8,
                        help="Highest rank to scan (default: 8)")
    parser.add_argument("-n", "--restarts", type=int, default=20,
                        help="Random restarts per rank (default: 20)")
    parser.add_argument("-d", "--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("-i", "--max-iter", type=int, default=500,
                        help="Max NMF iterations per restart (default: 500)")
    parser.add_argument("-r", "--rank", type=int, default=None,
                        help="Force this rank instead of auto-selection")
    parser.add_argument("-p", "--min-cophenetic", type=float, default=0.8,
                        help="Cophenetic floor for auto rank selection "
                             "(default: 0.8)")
    parser.add_argument("-E", "--elbow-ratio", type=float, default=0.25,
                        help="KL elbow detector: stop when the next rank's "
                             "improvement is below this ratio of the current "
                             "one (default: 0.25)")
    parser.add_argument("-s", "--signatures", default=None,
                        help="Optional COSMIC matrix TSV for cosine annotation")
    parser.add_argument("-l", "--log", default=None,
                        help="Optional log file (in addition to stdout)")
    return parser.parse_args()


def main():
    """Entry point: load matrix, scan ranks, extract and annotate signatures."""
    args = parse_args()
    if args.log:
        from pathlib import Path

        _log_path = Path(args.log)
        _log_path.parent.mkdir(parents=True, exist_ok=True)
        _lg = setup_logger("SignatureNMF", log_file=str(_log_path), stream=sys.stdout)
        logger.handlers = _lg.handlers
    import pandas as pd

    counts = pd.read_csv(args.counts, index_col=0)
    if counts.shape[0] != 96:
        raise ValueError(
            f"{args.counts}: expected 96 context rows, found {counts.shape[0]}; "
            f"pass the *_per_sample.csv written by spectrum.py."
        )
    if (counts.sum(axis=0) == 0).all():
        raise ValueError(f"{args.counts}: all samples have zero SNVs; nothing to decompose.")
    sig_ref = load_signatures(args.signatures, counts.index) if args.signatures else None
    run_nmf(counts, args, sig_ref)


if __name__ == "__main__":
    main()

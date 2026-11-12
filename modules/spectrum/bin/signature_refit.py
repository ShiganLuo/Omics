#!/usr/bin/env python
"""Reference-signature refitting: NNLS decomposition of SBS96 spectra onto COSMIC signatures.

Consumes the per-sample count matrix written by spectrum.py (96 contexts x
samples) and a COSMIC signature matrix (e.g. COSMIC_v3.6_SBS_GRCh38.txt or
COSMIC_v3.6_SBS_mm39.txt), and for every sample solves the non-negative
least-squares problem

    min_{h >= 0} || v - W h ||^2

where v is the sample's 96-channel spectrum and W the signature matrix.
Signatures whose exposure fraction falls below --min-exposure-frac are pruned
and the fit repeated (deconstructSigs-style iterative pruning) to keep the
attribution sparse. Per sample the output records the exposure vector, the
reconstruction cosine (goodness of fit) and a reliability flag driven by the
mutation count (samples below --min-mutations are statistically unreliable;
see the module design doc, section 5.1).

Outputs (<prefix> = --output-prefix):
  <prefix>_exposures.csv            samples x signatures, exposure counts
  <prefix>_exposure_fractions.csv   row-normalized exposure fractions
  <prefix>_family_exposures.csv     samples x signature families, fractions
  <prefix>_summary.csv              per-sample fit statistics and flags
  <prefix>_exposure_stacked_bar.png  family-fraction stacked bar per sample
  <prefix>_exposure_heatmap.png     sample x signature fraction heatmap
  <prefix>_signature_cosine.png     signature x signature cosine heatmap of
                                    the signatures used anywhere (degeneracy
                                    evidence: near-1 pairs make per-signature
                                    attribution non-unique)

Usage:
    python signature_refit.py \
        --counts somatic_spectrum_per_sample.csv \
        --signatures COSMIC_v3.6_SBS_GRCh38.txt \
        --output-prefix signature_refit \
        --min-mutations 200 --min-exposure-frac 0.02
"""

import argparse
import logging
import sys
from typing import Dict, List, Optional, TextIO, Tuple

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


logger = setup_logger("SignatureRefit", stream=sys.stdout)

# Built-in signature families. Collinear family members (SBS2/13, SBS1/5/40)
# are aggregated because single members are not separable at low mutation
# counts (design doc section 4.3). A --family-map TSV overrides this default.
DEFAULT_FAMILIES: Dict[str, List[str]] = {
    "Clock-like": ["SBS1", "SBS5", "SBS40a", "SBS40b", "SBS40c"],
    "APOBEC": ["SBS2", "SBS13"],
    "MMR": ["SBS6", "SBS15", "SBS20", "SBS26", "SBS44"],
    "POLE_POLD1": ["SBS10a", "SBS10b", "SBS10c", "SBS10d", "SBS14"],
    "HRD": ["SBS3"],
    "UV": ["SBS7a", "SBS7b", "SBS7c", "SBS7d"],
    "Tobacco": ["SBS4"],
    "Temozolomide": ["SBS11"],
    "Platinum": ["SBS31", "SBS35"],
    "Aristolochic_acid": ["SBS22a", "SBS22b", "SBS22c"],
    "Aflatoxin": ["SBS24"],
}


def load_counts(path: str):
    """Load the 96 x samples count matrix (index = context labels)."""
    import pandas as pd

    mat = pd.read_csv(path, index_col=0)
    if mat.shape[0] != 96:
        raise ValueError(
            f"{path}: expected 96 context rows, found {mat.shape[0]}; "
            f"pass the *_per_sample.csv written by spectrum.py."
        )
    return mat


def load_signatures(path: str, contexts):
    """Load a COSMIC signature matrix and align rows to *contexts*.

    Parameters
    ----------
    path : str
        COSMIC matrix TSV (first column ``Type`` = context label).
    contexts : pandas.Index
        Expected context labels (96 rows).

    Returns
    -------
    pandas.DataFrame
        96 x K signature matrix with rows in *contexts* order.
    """
    import pandas as pd

    sig = pd.read_csv(path, sep="\t", index_col=0)
    if set(sig.index) != set(contexts):
        raise ValueError(
            f"{path}: context labels do not match the count matrix "
            f"({len(set(sig.index) ^ set(contexts))} labels differ); the signature "
            f"matrix must use the same 96 SBS channels as the spectrum output."
        )
    return sig.loc[contexts]


def cosine(a, b) -> float:
    """Cosine similarity of two flat arrays."""
    import numpy as np

    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na == 0.0 or nb == 0.0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def refit_sample(v, w_mat, names: List[str], min_exposure_frac: float,
                 max_rounds: int = 500) -> Tuple:
    """NNLS-fit one spectrum with iterative low-exposure pruning.

    Parameters
    ----------
    v : numpy.ndarray
        96-channel observed spectrum (counts).
    w_mat : numpy.ndarray
        96 x K signature matrix (columns are signatures).
    names : list of str
        Signature names, aligned with w_mat columns.
    min_exposure_frac : float
        Drop signatures whose exposure fraction falls below this value and
        refit; stop when the active set is stable.

    Returns
    -------
    (exposures, reconstruction, active_names) where exposures is a dict of
    name -> exposure count over the kept signatures.
    """
    import numpy as np
    from scipy.optimize import nnls

    active = list(range(w_mat.shape[1]))
    h_full = np.zeros(w_mat.shape[1])
    for _ in range(max_rounds):
        h, _ = nnls(w_mat[:, active], v)
        total = h.sum()
        if total <= 0 or len(active) == 1:
            break
        fracs = h / total
        min_idx = int(np.argmin(fracs))
        if fracs[min_idx] >= min_exposure_frac:
            break
        # Remove only the weakest signature per round and refit, so the
        # remaining signatures can absorb its mass before the next cut.
        # Batch pruning can drop a genuinely active signature whose mass
        # the initial fit had spread across correlated columns.
        active.pop(min_idx)
    h, _ = nnls(w_mat[:, active], v)
    h_full[active] = h
    reconstruction = w_mat @ h_full
    exposures = {names[i]: float(hi) for i, hi in zip(active, h) if hi > 0}
    return exposures, reconstruction, [names[i] for i in active]


def read_family_map(path: Optional[str]) -> Dict[str, List[str]]:
    """Read a family map TSV (``family<TAB>SBS1,SBS2``) or the built-in map."""
    if not path:
        return dict(DEFAULT_FAMILIES)
    families: Dict[str, List[str]] = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            family, members = line.split("\t", 1)
            families[family] = [m.strip() for m in members.split(",") if m.strip()]
    return families


def aggregate_families(exposure_dict: Dict[str, float],
                       families: Dict[str, List[str]]) -> Dict[str, float]:
    """Sum exposures into families; ungrouped signatures form their own family."""
    member_to_family = {
        member: family for family, members in families.items() for member in members
    }
    agg: Dict[str, float] = {}
    for sig, value in exposure_dict.items():
        family = member_to_family.get(sig, sig)
        agg[family] = agg.get(family, 0.0) + value
    return agg


def run_refit(counts, sig, args) -> None:
    """Refit every sample and write the exposure/summary outputs and plots."""
    import numpy as np
    import pandas as pd
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(sig.columns)
    w_mat = sig.to_numpy(dtype=float)
    samples = list(counts.columns)
    families = read_family_map(args.family_map)

    exposures: Dict[str, Dict[str, float]] = {}
    fractions: Dict[str, Dict[str, float]] = {}
    family_frac: Dict[str, Dict[str, float]] = {}
    summary_rows: List[List] = []

    for sample in samples:
        v = counts[sample].to_numpy(dtype=float)
        n_mut = int(v.sum())
        exp, recon, active = refit_sample(v, w_mat, names, args.min_exposure_frac)
        total = sum(exp.values())
        frac = {k: (val / total if total > 0 else 0.0) for k, val in exp.items()}
        fam = aggregate_families(frac, families)
        exposures[sample] = exp
        fractions[sample] = frac
        family_frac[sample] = fam
        top = max(frac.items(), key=lambda kv: kv[1]) if frac else ("", 0.0)
        summary_rows.append([
            sample, n_mut, n_mut >= args.min_mutations,
            len(active), cosine(v, recon), top[0], top[1],
        ])
        logger.info(
            f"{sample}: n={n_mut} reliable={n_mut >= args.min_mutations} "
            f"signatures={active} recon_cosine={cosine(v, recon):.4f}"
        )

    pd.DataFrame(exposures).T.reindex(samples).to_csv(f"{args.output_prefix}_exposures.csv")
    pd.DataFrame(fractions).T.reindex(samples).to_csv(f"{args.output_prefix}_exposure_fractions.csv")
    fam_df = pd.DataFrame(family_frac).T.reindex(samples).fillna(0.0)
    fam_df.to_csv(f"{args.output_prefix}_family_exposures.csv")

    summary = pd.DataFrame(summary_rows, columns=[
        "sample_id", "n_mutations", "reliable", "n_signatures_used",
        "reconstruction_cosine", "top_signature", "top_signature_fraction",
    ])
    summary.to_csv(f"{args.output_prefix}_summary.csv", index=False)

    # Family-fraction stacked bar per sample.
    ax = fam_df.plot(kind="bar", stacked=True, figsize=(max(6, 1.2 * len(samples)), 6))
    ax.set_ylabel("Exposure fraction")
    ax.set_xlabel("Sample")
    ax.set_ylim(0, 1)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_exposure_stacked_bar.png", dpi=300)
    plt.close()

    # Sample x signature fraction heatmap (only signatures used somewhere).
    frac_df = pd.DataFrame(fractions).T.reindex(samples).fillna(0.0)
    used = [c for c in frac_df.columns if bool((frac_df[c] > 0).any())]
    if used:
        plt.figure(figsize=(max(6, 0.35 * len(used)), max(4, 0.4 * len(samples))))
        import seaborn as sns

        sns.heatmap(frac_df[used], cmap="Reds", annot=False)
        plt.ylabel("Sample")
        plt.xlabel("Signature")
        plt.tight_layout()
        plt.savefig(f"{args.output_prefix}_exposure_heatmap.png", dpi=300)
        plt.close()

    # Signature x signature cosine heatmap over the used signatures
    # (design doc §6 common outputs). The refit solution is only unique up
    # to collinear columns: two near-identical signatures can trade mass
    # with no change in reconstruction, so a high off-diagonal value marks
    # attribution that must be read at family level, not per signature.
    if used:
        cols = [names.index(u) for u in used]
        sub = w_mat[:, cols]
        norms = np.linalg.norm(sub, axis=0)
        cos_mat = (sub.T @ sub) / np.outer(norms, norms)
        import seaborn as sns

        side = max(5, 0.5 * len(used))
        plt.figure(figsize=(side, side * 0.85))
        sns.heatmap(
            pd.DataFrame(cos_mat, index=used, columns=used),
            cmap="Blues", vmin=0, vmax=1, annot=True, fmt=".2f",
            annot_kws={"fontsize": 7},
        )
        plt.title("Signature-signature cosine (used signatures)")
        plt.tight_layout()
        plt.savefig(f"{args.output_prefix}_signature_cosine.png", dpi=300)
        plt.close()

    logger.info(
        f"Wrote {args.output_prefix}_exposures.csv, "
        f"{args.output_prefix}_exposure_fractions.csv, "
        f"{args.output_prefix}_family_exposures.csv, "
        f"{args.output_prefix}_summary.csv, "
        f"{args.output_prefix}_exposure_stacked_bar.png, "
        f"{args.output_prefix}_exposure_heatmap.png, "
        f"{args.output_prefix}_signature_cosine.png"
    )


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="NNLS refitting of SBS96 spectra onto COSMIC signatures")
    parser.add_argument("-c", "--counts", required=True,
                        help="Per-sample count matrix CSV from spectrum.py "
                             "(96 contexts x samples)")
    parser.add_argument("-s", "--signatures", required=True,
                        help="COSMIC signature matrix TSV (e.g. "
                             "COSMIC_v3.6_SBS_GRCh38.txt / COSMIC_v3.6_SBS_mm39.txt)")
    parser.add_argument("-o", "--output-prefix", required=True,
                        help="Output prefix")
    parser.add_argument("-m", "--min-mutations", type=int, default=200,
                        help="Samples with fewer SNVs are flagged unreliable "
                             "(default: 200)")
    parser.add_argument("-e", "--min-exposure-frac", type=float, default=0.02,
                        help="Prune signatures below this exposure fraction "
                             "(default: 0.02)")
    parser.add_argument("-f", "--family-map", default=None,
                        help="Optional TSV 'family<TAB>SBS1,SBS2' overriding the "
                             "built-in family map")
    parser.add_argument("-l", "--log", default=None,
                        help="Optional log file (in addition to stdout)")
    return parser.parse_args()


def main():
    """Entry point: load matrices, refit every sample, write outputs."""
    args = parse_args()
    if args.log:
        from pathlib import Path

        _log_path = Path(args.log)
        _log_path.parent.mkdir(parents=True, exist_ok=True)
        _lg = setup_logger("SignatureRefit", log_file=str(_log_path), stream=sys.stdout)
        logger.handlers = _lg.handlers
    counts = load_counts(args.counts)
    sig = load_signatures(args.signatures, counts.index)
    run_refit(counts, sig, args)


if __name__ == "__main__":
    main()

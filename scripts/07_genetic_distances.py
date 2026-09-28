#!/usr/bin/env python3
"""
07_genetic_distances.py
-----------------------
Compute pairwise genetic distances from the concatenated supermatrix.

Method: p-distance (proportion of differing sites, gaps excluded).
Also provides Jukes-Cantor corrected distances.

Outputs:
  results/distance_matrix.csv   – full N×N distance matrix
  figures/distance_heatmap.png  – annotated heatmap

Usage:
    python scripts/07_genetic_distances.py
"""

from __future__ import annotations

import logging
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO

# ─── Constants ────────────────────────────────────────────────────────────────

SUPERMATRIX   = Path("data/concatenated/supermatrix.fasta")
RESULTS_DIR   = Path("results")
FIGURES_DIR   = Path("figures")
DIST_CSV      = RESULTS_DIR / "distance_matrix.csv"
HEATMAP_PNG   = FIGURES_DIR / "distance_heatmap.png"
LOG_FILE      = Path("logs/07_distances.log")

# ─── Logging ──────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_FILE, mode="w"),
    ],
)
log = logging.getLogger(__name__)


# ─── Distance functions ───────────────────────────────────────────────────────

def p_distance(seq1: str, seq2: str) -> float:
    """
    Compute the p-distance between two aligned sequences.
    Sites with a gap or ambiguity character in either sequence are excluded.
    """
    SKIP = set("-?NX")
    n_diff = n_comp = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a in SKIP or b in SKIP:
            continue
        n_comp += 1
        if a != b:
            n_diff += 1
    return n_diff / n_comp if n_comp > 0 else 0.0


def jukes_cantor(p: float) -> float:
    """
    Jukes-Cantor distance from p-distance.
    Returns nan if p >= 0.75 (formula undefined).
    """
    import math
    if p >= 0.75:
        return float("nan")
    return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)


# ─── Matrix computation ───────────────────────────────────────────────────────

def compute_distance_matrix(
    sequences: dict[str, str],
) -> tuple[np.ndarray, list[str]]:
    """
    Compute the symmetric N×N p-distance matrix.

    Args:
        sequences: {species_name: aligned_sequence_string}

    Returns:
        matrix    – numpy float64 array
        labels    – ordered list of species names
    """
    labels = sorted(sequences.keys())
    n = len(labels)
    matrix = np.zeros((n, n), dtype=np.float64)

    for i, j in combinations(range(n), 2):
        d = p_distance(sequences[labels[i]], sequences[labels[j]])
        matrix[i, j] = d
        matrix[j, i] = d

    return matrix, labels


# ─── Visualization ────────────────────────────────────────────────────────────

def plot_heatmap(matrix: np.ndarray, labels: list[str], out_path: Path) -> None:
    """
    Generate and save an annotated distance heatmap.
    """
    n = len(labels)
    # Clean species names for display: replace underscores with spaces, italicise
    display_labels = [sp.replace("_", " ") for sp in labels]

    fig_size = max(10, n * 0.6)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.85))

    # Mask the diagonal (self-distances = 0)
    mask = np.eye(n, dtype=bool)

    cmap = sns.color_palette("YlOrRd", as_cmap=True)
    sns.heatmap(
        matrix,
        mask=mask,
        xticklabels=display_labels,
        yticklabels=display_labels,
        cmap=cmap,
        annot=n <= 25,          # show values only for small matrices
        fmt=".3f",
        linewidths=0.4,
        linecolor="white",
        vmin=0,
        ax=ax,
        cbar_kws={"label": "p-distance", "shrink": 0.7},
    )

    ax.set_title(
        "Pairwise Genetic Distances — Mitochondrial and Nuclear Supermatrix",
        fontsize=13,
        fontweight="bold",
        pad=14,
    )
    plt.xticks(rotation=45, ha="right", fontsize=9, style="italic")
    plt.yticks(rotation=0, fontsize=9, style="italic")
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Heatmap saved → %s", out_path)


# ─── Distance summary ─────────────────────────────────────────────────────────

def distance_summary(matrix: np.ndarray, labels: list[str]) -> None:
    """Log closest and most-divergent species pairs."""
    n = len(labels)
    pairs: list[tuple[float, str, str]] = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((matrix[i, j], labels[i], labels[j]))

    pairs.sort()
    log.info("─── 5 closest pairs ───────────────────────────────")
    for d, sp1, sp2 in pairs[:5]:
        log.info("  %.4f  %s ↔ %s", d, sp1.replace("_", " "), sp2.replace("_", " "))

    log.info("─── 5 most divergent pairs ─────────────────────────")
    for d, sp1, sp2 in pairs[-5:][::-1]:
        log.info("  %.4f  %s ↔ %s", d, sp1.replace("_", " "), sp2.replace("_", " "))


# ─── Main ─────────────────────────────────────────────────────────────────────

def run() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    if not SUPERMATRIX.exists():
        log.error("Supermatrix not found: %s — run Step 5 first.", SUPERMATRIX)
        sys.exit(1)

    log.info("Loading supermatrix …")
    sequences: dict[str, str] = {}
    for rec in SeqIO.parse(str(SUPERMATRIX), "fasta"):
        sequences[rec.id] = str(rec.seq)

    if len(sequences) < 2:
        log.error("Need at least 2 sequences to compute distances.")
        sys.exit(1)

    log.info("Computing pairwise p-distances for %d species …", len(sequences))
    matrix, labels = compute_distance_matrix(sequences)

    # Also compute JC matrix and store alongside p-distance
    jc_matrix = np.vectorize(jukes_cantor)(matrix)

    # Save as CSV (p-distance)
    df = pd.DataFrame(matrix, index=labels, columns=labels)
    df.to_csv(str(DIST_CSV))
    log.info("Distance matrix → %s", DIST_CSV)

    # Save JC matrix
    jc_path = RESULTS_DIR / "jc_distance_matrix.csv"
    df_jc = pd.DataFrame(jc_matrix, index=labels, columns=labels)
    df_jc.to_csv(str(jc_path))
    log.info("JC distance matrix → %s", jc_path)

    distance_summary(matrix, labels)
    plot_heatmap(matrix, labels, HEATMAP_PNG)

    log.info("Step 7 complete.")


if __name__ == "__main__":
    run()

#!/usr/bin/env python3
"""
10_summary_figure.py
--------------------
Generate a single publication-quality summary figure combining:
  Panel A: Phylogenetic tree (from IQ-TREE ML result)
  Panel B: Pairwise distance heatmap
  Panel C: Gene coverage & bootstrap bar charts

Output: results/summary_figure.png  (300 dpi)
        results/summary_figure.pdf

Usage:
    python scripts/10_summary_figure.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo, SeqIO

# ─── Constants ────────────────────────────────────────────────────────────────
RESULTS_DIR  = Path("results")
FIGURES_DIR  = Path("figures")
TREE_PREFIX  = Path("results/tree/butterfly")
DIST_CSV     = RESULTS_DIR / "distance_matrix.csv"
ALIGN_QC_CSV = Path("data/aligned/alignment_qc_per_gene.csv")
SUPERMATRIX = Path("data/concatenated/supermatrix.fasta")

LOG_FILE = Path("logs/10_summary.log")

FAMILY_COLORS = {
    "Papilio":      "#E63946",
    "Pieris":       "#457B9D",
    "Gonepteryx":   "#457B9D",
    "Aporia":       "#457B9D",
    "Eurema":       "#457B9D",
    "Danaus":       "#2A9D8F",
    "Vanessa":      "#2A9D8F",
    "Junonia":      "#2A9D8F",
    "Melitaea":     "#2A9D8F",
    "Lycaena":      "#F4A261",
    "Plebejus":     "#F4A261",
    "Curetis":      "#F4A261",
    "Ampittia":     "#8338EC",
    "Ochlodes":     "#8338EC",
    "Parnara":      "#8338EC",
    "Heteropterus": "#8338EC",
    "Pyrgus":       "#8338EC",
    "Celaenorrhinus": "#8338EC",
    "Ctenoptilum":  "#8338EC",
    "Notocrypta":   "#8338EC",
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout),
                               logging.FileHandler(LOG_FILE, mode="w")])
log = logging.getLogger(__name__)


def get_tip_color(name: str) -> str:
    genus = name.split("_")[0].split("|")[0]
    return FAMILY_COLORS.get(genus, "#555555")


def draw_tree_panel(ax: plt.Axes) -> None:
    treefile = Path(str(TREE_PREFIX) + ".contree")
    if not treefile.exists():
        treefile = Path(str(TREE_PREFIX) + ".treefile")
    if not treefile.exists():
        ax.text(0.5, 0.5, "Tree file not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11)
        return

    tree = Phylo.read(str(treefile), "newick")
    tree.root_at_midpoint()
    tree.ladderize()

    for clade in tree.get_terminals():
        raw = clade.name or ""
        genus = raw.split("_")[0]
        species = raw.split("_")[1] if "_" in raw else ""
        clade.name = f"{genus} {species}".strip()

    Phylo.draw(
        tree, axes=ax, do_show=False,
        show_confidence=True,
        branch_labels=lambda c: (
            f"{int(c.confidence)}" if c.confidence is not None and c.confidence >= 70 else ""
        ),
    )

    for text_obj in ax.texts:
        txt = text_obj.get_text().strip()
        if txt and any(c.islower() for c in txt):
            text_obj.set_fontstyle("italic")
            text_obj.set_fontsize(8)
            genus = txt.split()[0]
            text_obj.set_color(FAMILY_COLORS.get(genus, "#111111"))

    # Family legend
    legend_elements = [
        mpatches.Patch(color="#E63946", label="Papilionidae"),
        mpatches.Patch(color="#457B9D", label="Pieridae"),
        mpatches.Patch(color="#2A9D8F", label="Nymphalidae"),
        mpatches.Patch(color="#F4A261", label="Lycaenidae"),
        mpatches.Patch(color="#8338EC", label="Hesperiidae"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=7,
              framealpha=0.85, title="Family", title_fontsize=7)
    ax.set_title("A   Maximum-Likelihood Phylogenetic Tree", fontweight="bold",
                 fontsize=10, loc="left", pad=6)
    ax.set_xlabel("Substitutions per site", fontsize=8)
    ax.set_ylabel("")


def draw_heatmap_panel(ax: plt.Axes) -> None:
    if not DIST_CSV.exists():
        ax.text(0.5, 0.5, "Distance matrix not found", ha="center", va="center",
                transform=ax.transAxes)
        return

    df = pd.read_csv(str(DIST_CSV), index_col=0)
    labels = [c.replace("_", " ") for c in df.columns]
    matrix = df.values
    n = len(labels)

    mask = np.eye(n, dtype=bool)
    vmax = max(0.25, min(float(np.nanmax(matrix)), 0.5))
    sns.heatmap(
        matrix, mask=mask,
        xticklabels=labels, yticklabels=labels,
        cmap="YlOrRd", ax=ax,
        annot=False,
        vmin=0, vmax=vmax,
        linewidths=0.3, linecolor="white",
        cbar_kws={"label": "p-distance", "shrink": 0.7},
    )
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=6.5, rotation=45,
                       ha="right", style="italic")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6.5, rotation=0,
                       style="italic")
    ax.set_title("B   Pairwise Genetic Distance Heatmap (p-distance)",
                 fontweight="bold", fontsize=10, loc="left", pad=6)


def draw_stats_panel(ax_top: plt.Axes, ax_bot: plt.Axes) -> None:
    # ── Gene alignment stats bar chart ──────────────────────────────────────
    if ALIGN_QC_CSV.exists():
        df = pd.read_csv(str(ALIGN_QC_CSV))
        genes = df["gene"].tolist()
        x = np.arange(len(genes))
        w = 0.28
        ax_top.bar(x - w, df["conserved_sites"],   width=w, label="Conserved",
                   color="#2196F3", alpha=0.85)
        ax_top.bar(x,      df["variable_sites"],    width=w, label="Variable",
                   color="#FF9800", alpha=0.85)
        ax_top.bar(x + w,  df["parsimony_informative_sites"], width=w,
                   label="Parsimony-inf.", color="#4CAF50", alpha=0.85)
        ax_top.set_xticks(x)
        ax_top.set_xticklabels(genes, fontsize=9)
        ax_top.set_ylabel("Sites (bp)", fontsize=8)
        ax_top.legend(fontsize=7, framealpha=0.8)
        ax_top.set_title("C   Alignment Site Composition per Gene",
                          fontweight="bold", fontsize=10, loc="left", pad=6)
        ax_top.tick_params(axis="y", labelsize=8)
        ax_top.grid(axis="y", alpha=0.3, linestyle="--")
    else:
        ax_top.text(0.5, 0.5, "Alignment QC data not found",
                    ha="center", va="center", transform=ax_top.transAxes)

    # ── Bootstrap distribution bar chart ────────────────────────────────────
    bins = ["<50", "50-69", "70-94", "≥95"]
    treefile = Path(str(TREE_PREFIX) + ".contree")
    if not treefile.exists():
        treefile = Path(str(TREE_PREFIX) + ".treefile")
    support = []
    if treefile.exists():
        tree = Phylo.read(str(treefile), "newick")
        support = [float(c.confidence) for c in tree.find_clades()
                   if c.confidence is not None]
    counts = [
        sum(v < 50 for v in support),
        sum(50 <= v < 70 for v in support),
        sum(70 <= v < 95 for v in support),
        sum(v >= 95 for v in support),
    ]
    colors = ["#F44336", "#FF9800", "#8BC34A", "#2196F3"]
    bars = ax_bot.bar(bins, counts, color=colors, alpha=0.9, edgecolor="white", linewidth=1.2)
    for bar, count in zip(bars, counts):
        ax_bot.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
                    str(count), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax_bot.set_ylabel("Number of nodes", fontsize=8)
    ax_bot.set_xlabel("Ultrafast Bootstrap (%)", fontsize=8)
    ax_bot.set_title(f"D   Bootstrap Support Distribution (n={len(support)} nodes)",
                      fontweight="bold", fontsize=10, loc="left", pad=6)
    ax_bot.tick_params(labelsize=8)
    ax_bot.set_ylim(0, max(counts, default=0) + 2)
    ax_bot.grid(axis="y", alpha=0.3, linestyle="--")

    # Annotation
    ax_bot.axhline(0, color="black", linewidth=0.5)
    mean_support = float(np.mean(support)) if support else float("nan")
    summary = f"Mean BS = {mean_support:.1f}%" if support else "Bootstrap support unavailable"
    ax_bot.text(0.97, 0.95, summary, transform=ax_bot.transAxes,
                ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#E3F2FD", alpha=0.8))


def run() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(20, 14))
    fig.patch.set_facecolor("#FAFAFA")

    gs = gridspec.GridSpec(
        2, 2,
        left=0.05, right=0.97,
        top=0.93, bottom=0.06,
        hspace=0.38, wspace=0.30,
        height_ratios=[1.6, 1],
    )

    ax_tree = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    ax_aln  = fig.add_subplot(gs[1, 0])
    ax_bs   = fig.add_subplot(gs[1, 1])

    draw_tree_panel(ax_tree)
    draw_heatmap_panel(ax_heat)
    draw_stats_panel(ax_aln, ax_bs)

    n_taxa = sum(1 for _ in SeqIO.parse(str(SUPERMATRIX), "fasta")) if SUPERMATRIX.exists() else 0
    fig.suptitle(
        f"Phylogenomics of {n_taxa} Butterfly Taxa — Mitochondrial and Nuclear Markers",
        fontsize=13, fontweight="bold", y=0.975,
    )

    for ext in ("png", "pdf"):
        out = RESULTS_DIR / f"summary_figure.{ext}"
        fig.savefig(str(out), dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
        log.info("Saved → %s", out)

    plt.close(fig)
    log.info("Summary figure generation complete.")


if __name__ == "__main__":
    run()

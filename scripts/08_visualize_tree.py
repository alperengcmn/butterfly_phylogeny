#!/usr/bin/env python3
"""
08_visualize_tree.py
--------------------
Generate publication-quality phylogenetic tree figures.

Two rendering backends are used:
  1. Biopython Phylo + Matplotlib  →  figures/tree_biopython.png / .pdf
  2. ETE3 (if available)           →  figures/tree_ete3.png / .pdf

Both figures include:
  - Bootstrap support labels on internal nodes
  - Italicised species name labels on tips
  - Scale bar

Usage:
    python scripts/08_visualize_tree.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo

# ─── Constants ────────────────────────────────────────────────────────────────

TREE_PREFIX  = Path("results/tree/butterfly")
FIGURES_DIR  = Path("figures")
LOG_FILE     = Path("logs/08_visualize.log")

# Expected IQ-TREE output (contree = consensus tree with bootstrap values)
TREEFILE     = Path(str(TREE_PREFIX) + ".contree")
TREEFILE_ML  = Path(str(TREE_PREFIX) + ".treefile")

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


# ─── Helpers ──────────────────────────────────────────────────────────────────

def find_treefile() -> Path:
    """Return the best available tree file."""
    for candidate in (TREEFILE, TREEFILE_ML):
        if candidate.exists():
            log.info("Using tree file: %s", candidate)
            return candidate
    raise FileNotFoundError(
        "No IQ-TREE output found. Run Step 6 first.\n"
        f"Looked for: {TREEFILE}, {TREEFILE_ML}"
    )


def clean_species_name(name: str) -> str:
    """Convert 'Genus_species|...' to 'Genus species' for display."""
    base = name.split("|")[0]
    return base.replace("_", " ")


# ─── Biopython tree drawing ───────────────────────────────────────────────────

def draw_biopython_tree(tree_path: Path) -> None:
    """
    Draw the tree using Biopython Phylo + Matplotlib.
    Saves PNG and PDF to figures/.
    """
    tree = Phylo.read(str(tree_path), "newick")
    tree.root_at_midpoint()
    tree.ladderize()

    # Rename tips to pretty species names
    for clade in tree.get_terminals():
        clade.name = clean_species_name(clade.name or "")

    # Count leaves to set figure height
    n_tips = len(tree.get_terminals())
    fig_h  = max(8, n_tips * 0.45)
    fig, ax = plt.subplots(figsize=(14, fig_h))

    # Phylo.draw sets x/y positions; we capture them
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=True,
        branch_labels=lambda c: (
            f"{int(c.confidence)}" if c.confidence is not None and c.confidence > 0 else ""
        ),
    )

    # Style tweaks
    ax.set_title(
        "Maximum-Likelihood Phylogenetic Tree of 20 Butterfly Species\n"
        "Six-locus Supermatrix (mtDNA + nuclear markers)",
        fontsize=12, fontweight="bold", pad=12,
    )
    ax.set_xlabel("Substitutions per site", fontsize=10)
    ax.set_ylabel("")

    # Italicise tip labels
    for text_obj in ax.texts:
        txt = text_obj.get_text().strip()
        if txt and len(txt.split()) >= 2:  # looks like a species name
            text_obj.set_fontstyle("italic")
            text_obj.set_fontsize(9)

    # Bootstrap legend
    legend_patch = mpatches.Patch(
        color="grey", label="Node labels = ultrafast bootstrap support (%)"
    )
    ax.legend(handles=[legend_patch], loc="lower right", fontsize=8)

    plt.tight_layout()

    for ext in ("png", "pdf"):
        out = FIGURES_DIR / f"tree_biopython.{ext}"
        plt.savefig(str(out), dpi=300, bbox_inches="tight")
        log.info("Saved → %s", out)

    plt.close(fig)


# ─── ETE3 tree drawing ────────────────────────────────────────────────────────

def draw_ete3_tree(tree_path: Path) -> None:
    """
    Draw the tree using ETE3.
    Falls back gracefully if ETE3 is not installed or headless rendering fails.
    """
    try:
        from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces
    except ImportError:
        log.warning("ETE3 not installed — skipping ETE3 figure. pip install ete3")
        return

    try:
        t = Tree(str(tree_path))
        t.set_outgroup(t.get_midpoint_outgroup())
        t.ladderize()

        # ── Tree style ──────────────────────────────────────────────────
        ts = TreeStyle()
        ts.mode              = "r"   # rectangular
        ts.show_leaf_name    = False # we add custom labels
        ts.show_branch_length= True
        ts.show_branch_support = True
        ts.scale             = 1200
        ts.branch_vertical_margin = 12
        ts.title.add_face(
            TextFace(
                "Butterfly Phylogenetic Tree — Mitochondrial Supermatrix",
                fsize=14, bold=True,
            ),
            column=0,
        )

        # ── Node styles ────────────────────────────────────────────────
        for node in t.traverse():
            ns = NodeStyle()
            ns["size"] = 0

            if node.is_leaf():
                # Italicised species name
                label = clean_species_name(node.name)
                face  = TextFace(label, fsize=10, fstyle="italic")
                face.margin_left = 4
                node.add_face(face, column=0, position="branch-right")
                ns["size"] = 0
            else:
                # Bootstrap label on internal nodes (only if ≥50)
                support = getattr(node, "support", None)
                if support is not None and float(support) >= 50:
                    bface = TextFace(f"{int(float(support))}", fsize=7, fgcolor="#555555")
                    node.add_face(bface, column=0, position="branch-top")
                ns["fgcolor"] = "#333333"
                ns["size"]    = 4

            node.set_style(ns)

        # ── Export ─────────────────────────────────────────────────────
        for ext in ("png", "pdf"):
            out = FIGURES_DIR / f"tree_ete3.{ext}"
            t.render(str(out), tree_style=ts, dpi=300, w=1600)
            log.info("ETE3 figure saved → %s", out)

    except Exception as exc:
        log.warning("ETE3 rendering failed: %s", exc)


# ─── Main ─────────────────────────────────────────────────────────────────────

def run() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    try:
        tree_path = find_treefile()
    except FileNotFoundError as exc:
        log.error(str(exc))
        sys.exit(1)

    log.info("Drawing Biopython tree …")
    draw_biopython_tree(tree_path)

    log.info("Drawing ETE3 tree …")
    draw_ete3_tree(tree_path)

    log.info("Step 8 complete.")


if __name__ == "__main__":
    run()

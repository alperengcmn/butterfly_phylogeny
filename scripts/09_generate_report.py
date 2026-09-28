#!/usr/bin/env python3
"""
09_generate_report.py
---------------------
Aggregate outputs from all pipeline steps into human-readable reports.

Outputs:
  results/project_summary.txt      – pipeline-wide statistics
  results/evolutionary_summary.txt – biological interpretation

Usage:
    python scripts/09_generate_report.py
"""

from __future__ import annotations

import logging
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO

# ─── Constants ────────────────────────────────────────────────────────────────

RESULTS_DIR   = Path("results")
FIGURES_DIR   = Path("figures")
LOG_FILE      = Path("logs/09_report.log")

PROJECT_SUMMARY     = RESULTS_DIR / "project_summary.txt"
EVO_SUMMARY         = RESULTS_DIR / "evolutionary_summary.txt"

# Input files
METADATA_CSV        = Path("data/raw/metadata.csv")
QC_CSV              = Path("data/cleaned/qc_stats.csv")
ALIGN_STATS_CSV     = Path("data/aligned/alignment_stats.csv")
ALIGN_QC_PER_GENE   = Path("data/aligned/alignment_qc_per_gene.csv")
SUPERMATRIX         = Path("data/concatenated/supermatrix.fasta")
PARTITIONS          = Path("data/concatenated/partitions.txt")
DIST_MATRIX_CSV     = RESULTS_DIR / "distance_matrix.csv"
IQTREE_LOG          = Path("results/tree/butterfly.iqtree")
TREEFILE            = Path("results/tree/butterfly.contree")
TREEFILE_ML         = Path("results/tree/butterfly.treefile")

SPECIES: list[str] = [
    "Papilio machaon", "Papilio xuthus", "Papilio glaucus",
    "Papilio polytes", "Papilio bianor", "Pieris rapae",
    "Pieris napi", "Gonepteryx rhamni", "Aporia crataegi",
    "Eurema hecabe", "Danaus plexippus", "Vanessa indica",
    "Vanessa cardui", "Junonia almana", "Melitaea cinxia",
    "Lycaena phlaeas", "Plebejus argus", "Curetis bulis",
    "Ampittia dioscorides", "Ochlodes venata", "Parnara guttata",
    "Heteropterus morpheus", "Pyrgus malvae",
    "Celaenorrhinus maculosus", "Ctenoptilum vasava", "Notocrypta curvifascia",
]

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


# ─── Data loaders ─────────────────────────────────────────────────────────────

def safe_read_csv(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(str(path))
    log.warning("File not found: %s", path)
    return pd.DataFrame()


def parse_best_model(log_path: Path) -> str:
    if not log_path.exists():
        return "N/A"
    text = log_path.read_text(errors="replace")
    m = re.search(r"Best-fit model(?: according to BIC)?:\s*(.+)", text)
    return m.group(1).strip() if m else "see IQ-TREE log"


def load_treefile() -> Path | None:
    for t in (TREEFILE, TREEFILE_ML):
        if t.exists():
            return t
    return None


def bootstrap_stats(tree_path: Path) -> dict:
    try:
        tree = Phylo.read(str(tree_path), "newick")
        vals = [
            float(c.confidence)
            for c in tree.find_clades()
            if c.confidence is not None
        ]
        if not vals:
            return {}
        return {
            "n_nodes":  len(vals),
            "min":      min(vals),
            "max":      max(vals),
            "mean":     sum(vals) / len(vals),
            "n_above_70": sum(1 for v in vals if v >= 70),
            "n_above_95": sum(1 for v in vals if v >= 95),
        }
    except Exception as exc:
        log.warning("Bootstrap stats failed: %s", exc)
        return {}


def species_pair_summary(dist_csv: Path) -> dict:
    if not dist_csv.exists():
        return {}
    df = pd.read_csv(str(dist_csv), index_col=0)
    labels = list(df.index)
    n = len(labels)
    pairs: list[tuple[float, str, str]] = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((df.iloc[i, j], labels[i], labels[j]))
    pairs.sort()
    return {
        "closest":     pairs[:5],
        "divergent":   pairs[-5:][::-1],
        "mean_dist":   float(np.mean([p[0] for p in pairs])),
        "max_dist":    float(pairs[-1][0]) if pairs else 0.0,
        "min_dist":    float(pairs[0][0])  if pairs else 0.0,
    }


# ─── Report builders ──────────────────────────────────────────────────────────

DIVIDER = "=" * 72


def build_project_summary() -> str:
    lines: list[str] = []

    def section(title: str) -> None:
        lines.append("")
        lines.append(DIVIDER)
        lines.append(f"  {title}")
        lines.append(DIVIDER)

    lines.append("BUTTERFLY PHYLOGENOMICS PROJECT — SUMMARY REPORT")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    # 1 — Species
    section("1. SPECIES")
    lines.append(f"  Total target species: {len(SPECIES)}")
    for sp in SPECIES:
        lines.append(f"    · {sp}")

    # 2 — Data collection
    section("2. DATA COLLECTION (Step 1)")
    meta = safe_read_csv(METADATA_CSV)
    if not meta.empty:
        lines.append(f"  Total sequences downloaded : {len(meta)}")
        for gene, grp in meta.groupby("gene"):
            lines.append(f"    {gene:6s}: {len(grp)} sequences")
    else:
        lines.append("  metadata.csv not found — re-run Step 1.")

    # 3 — QC & cleaning
    section("3. SEQUENCE QUALITY CONTROL (Step 2)")
    qc = safe_read_csv(QC_CSV)
    if not qc.empty:
        for _, row in qc.iterrows():
            lines.append(
                f"  {row['gene']:6s}: {row['input_seqs']} in → "
                f"{row['output_seqs']} out  "
                f"(-{row['removed_short']} short, "
                f"-{row['removed_ambig']} ambig, "
                f"-{row['removed_dup']} dup)"
            )
    else:
        lines.append("  qc_stats.csv not found.")

    # 4 — Alignments
    section("4. MULTIPLE SEQUENCE ALIGNMENT (Steps 3–4)")
    aln_qc = safe_read_csv(ALIGN_QC_PER_GENE)
    if not aln_qc.empty:
        header = f"  {'Gene':<8} {'Seqs':>5} {'Length':>8} {'Gap%':>6} {'Conserved':>10} {'Variable':>9} {'PI sites':>9}"
        lines.append(header)
        lines.append("  " + "-" * (len(header) - 2))
        for _, row in aln_qc.iterrows():
            lines.append(
                f"  {row['gene']:<8} {row['n_sequences']:>5} "
                f"{row['alignment_length']:>8} "
                f"{float(row['overall_gap_fraction']) * 100:>5.1f}% "
                f"{row['conserved_sites']:>10} "
                f"{row['variable_sites']:>9} "
                f"{row['parsimony_informative_sites']:>9}"
            )
    else:
        lines.append("  alignment_qc_per_gene.csv not found.")

    # 5 — Supermatrix
    section("5. CONCATENATED SUPERMATRIX (Step 5)")
    if SUPERMATRIX.exists():
        recs = list(SeqIO.parse(str(SUPERMATRIX), "fasta"))
        if recs:
            lines.append(f"  Species in supermatrix : {len(recs)}")
            lines.append(f"  Supermatrix length     : {len(recs[0].seq):,} bp")
    if PARTITIONS.exists():
        lines.append("  Partitions:")
        for line in PARTITIONS.read_text().splitlines():
            if "charset" in line.lower():
                lines.append("    " + line.strip())

    # 6 — IQ-TREE
    section("6. PHYLOGENETIC ANALYSIS — IQ-TREE (Step 6)")
    best_model = parse_best_model(IQTREE_LOG)
    lines.append(f"  Best-fit model         : {best_model}")
    lines.append(f"  Bootstrap replicates   : 1000 (ultrafast)")
    tree_path = load_treefile()
    if tree_path:
        bs = bootstrap_stats(tree_path)
        if bs:
            lines.append(f"  Internal nodes         : {bs['n_nodes']}")
            lines.append(f"  Bootstrap range        : {bs['min']:.0f} – {bs['max']:.0f}")
            lines.append(f"  Mean bootstrap support : {bs['mean']:.1f}")
            lines.append(f"  Nodes ≥ 70%            : {bs['n_above_70']}")
            lines.append(f"  Nodes ≥ 95%            : {bs['n_above_95']}")
    else:
        lines.append("  Tree file not found — run Step 6.")

    # 7 — Genetic distances
    section("7. GENETIC DISTANCES (Step 7)")
    dist = species_pair_summary(DIST_MATRIX_CSV)
    if dist:
        lines.append(f"  Mean pairwise p-distance : {dist['mean_dist']:.4f}")
        lines.append(f"  Min  pairwise p-distance : {dist['min_dist']:.4f}")
        lines.append(f"  Max  pairwise p-distance : {dist['max_dist']:.4f}")
    else:
        lines.append("  distance_matrix.csv not found.")

    # 8 — Output files
    section("8. OUTPUT FILES")
    for path in sorted(FIGURES_DIR.glob("*")) if FIGURES_DIR.exists() else []:
        lines.append(f"  figures/{path.name}")
    for path in sorted(RESULTS_DIR.glob("*")) if RESULTS_DIR.exists() else []:
        if path.is_file():
            lines.append(f"  results/{path.name}")
    inheritance_dir = RESULTS_DIR / "inheritance"
    if inheritance_dir.exists():
        for path in sorted(inheritance_dir.iterdir()):
            if path.is_file() and path.suffix in {".txt", ".treefile"}:
                lines.append(f"  results/inheritance/{path.name}")

    lines.append("")
    lines.append(DIVIDER)
    lines.append("END OF REPORT")
    lines.append(DIVIDER)

    return "\n".join(lines)


def build_evolutionary_summary() -> str:
    lines: list[str] = []

    def section(title: str) -> None:
        lines.append("")
        lines.append(DIVIDER)
        lines.append(f"  {title}")
        lines.append(DIVIDER)

    lines.append("BUTTERFLY PHYLOGENOMICS — EVOLUTIONARY SUMMARY")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    matrix_taxa = [r.id for r in SeqIO.parse(str(SUPERMATRIX), "fasta")] if SUPERMATRIX.exists() else []
    hesperiidae = {
        "Ampittia_dioscorides", "Ochlodes_venata", "Parnara_guttata",
        "Heteropterus_morpheus", "Pyrgus_malvae", "Celaenorrhinus_maculosus",
        "Ctenoptilum_vasava", "Notocrypta_curvifascia",
    }

    section("OVERVIEW")
    lines.append("  This analysis reconstructed the phylogenetic relationships among")
    lines.append(f"  {len(matrix_taxa)} taxa present in the supermatrix (from {len(SPECIES)} targets)")
    lines.append(f"  across 5 families; Hesperiidae contributes {len(set(matrix_taxa) & hesperiidae)} taxa")
    lines.append("  using four mitochondrial (COI, COII, CytB, ND5) and two nuclear")
    lines.append("  markers (EF1a and wingless), where available.")

    section("SPECIES PAIRS — CLOSEST (lowest p-distance)")
    dist = species_pair_summary(DIST_MATRIX_CSV)
    if dist.get("closest"):
        for d, sp1, sp2 in dist["closest"]:
            sp1_disp = sp1.replace("_", " ")
            sp2_disp = sp2.replace("_", " ")
            lines.append(f"  {d:.4f}  {sp1_disp} ↔ {sp2_disp}")
    else:
        lines.append("  Distance matrix not available.")

    section("SPECIES PAIRS — MOST DIVERGENT (highest p-distance)")
    if dist.get("divergent"):
        for d, sp1, sp2 in dist["divergent"]:
            sp1_disp = sp1.replace("_", " ")
            sp2_disp = sp2.replace("_", " ")
            lines.append(f"  {d:.4f}  {sp1_disp} ↔ {sp2_disp}")
    else:
        lines.append("  Distance matrix not available.")

    section("SUBSTITUTION MODEL")
    best_model = parse_best_model(IQTREE_LOG)
    lines.append(f"  Best-fit model selected by ModelFinder: {best_model}")
    lines.append("")
    lines.append("  GTR+F+I+G4 (if selected) indicates:")
    lines.append("    GTR — General Time-Reversible substitution matrix")
    lines.append("    +F  — empirical base frequencies")
    lines.append("    +I  — proportion of invariable sites")
    lines.append("    +G4 — gamma-distributed rate variation (4 categories)")

    section("BOOTSTRAP SUPPORT SUMMARY")
    tree_path = load_treefile()
    if tree_path:
        bs = bootstrap_stats(tree_path)
        if bs:
            lines.append(f"  Mean support       : {bs['mean']:.1f}%")
            lines.append(f"  Well-supported (≥70%): {bs['n_above_70']} / {bs['n_nodes']} nodes")
            lines.append(f"  Strongly supported (≥95%): {bs['n_above_95']} / {bs['n_nodes']} nodes")
            if bs["mean"] >= 85:
                lines.append("")
                lines.append("  Most internal nodes have high bootstrap support under this dataset.")
            else:
                lines.append("")
                lines.append("  ⚠ Several nodes have low support — consider adding")
                lines.append("    more loci or reviewing data quality.")
            lines.append("  Bootstrap support is node-level resampling support; it is not the")
            lines.append("  probability that the complete tree is correct.")
    else:
        lines.append("  Phylogenetic analysis results not yet available.")

    section("NOTES ON DATA COMPLETENESS")
    missing_csv = Path("data/aligned/missing_taxa_report.csv")
    if missing_csv.exists():
        df_missing = pd.read_csv(str(missing_csv))
        if df_missing.empty:
            lines.append(f"  All {len(SPECIES)} species represented across all loci.")
        else:
            lines.append(f"  {len(df_missing)} missing gene/species combinations:")
            for _, row in df_missing.iterrows():
                lines.append(f"    {row['missing_species'].replace('_', ' ')} — {row['gene']}")
            lines.append("")
            lines.append("  Missing data was filled with '?' in the supermatrix.")
    else:
        lines.append("  Missing taxa report not available.")

    if ALIGN_QC_PER_GENE.exists():
        alignment_qc = pd.read_csv(str(ALIGN_QC_PER_GENE))
        high_gap = alignment_qc[
            pd.to_numeric(alignment_qc["overall_gap_fraction"], errors="coerce") >= 0.20
        ]
        if not high_gap.empty:
            lines.append("")
            lines.append("  High-gap alignments requiring cautious interpretation:")
            for _, row in high_gap.iterrows():
                lines.append(
                    f"    {row['gene']}: {float(row['overall_gap_fraction']) * 100:.1f}% gaps"
                )

    comparison = Path("results/inheritance/comparison.txt")
    if comparison.exists():
        section("MITOCHONDRIAL VS NUCLEAR COMPARISON")
        lines.extend("  " + line for line in comparison.read_text(encoding="utf-8").splitlines())

    lines.append("")
    lines.append(DIVIDER)
    lines.append("END OF EVOLUTIONARY SUMMARY")
    lines.append(DIVIDER)

    return "\n".join(lines)


# ─── Main ─────────────────────────────────────────────────────────────────────

def run() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    log.info("Building project_summary.txt …")
    proj = build_project_summary()
    PROJECT_SUMMARY.write_text(proj, encoding="utf-8")
    log.info("Written → %s", PROJECT_SUMMARY)

    log.info("Building evolutionary_summary.txt …")
    evo = build_evolutionary_summary()
    EVO_SUMMARY.write_text(evo, encoding="utf-8")
    log.info("Written → %s", EVO_SUMMARY)

    log.info("Report generation complete.")


if __name__ == "__main__":
    run()

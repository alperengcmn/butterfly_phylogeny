#!/usr/bin/env python3
"""
04_alignment_qc.py
------------------
Detailed quality-control report on each aligned FASTA produced in Step 3.

Calculates per-gene and per-sequence metrics:
  - Alignment length
  - Gap fraction (overall and per-sequence)
  - Conserved / variable / parsimony-informative site counts
  - GC content per sequence
  - Missing-taxa report

Outputs:
  data/aligned/alignment_qc_per_gene.csv   – one row per gene
  data/aligned/alignment_qc_per_seq.csv    – one row per sequence × gene
  data/aligned/missing_taxa_report.csv     – species missing in each gene

Usage:
    python scripts/04_alignment_qc.py
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# ─── Constants ────────────────────────────────────────────────────────────────

GENES: list[str] = ["COI", "COII", "CytB", "ND5"]

ALIGNED_DIR = Path("data/aligned")
LOG_FILE    = Path("logs/04_qc.log")

QC_PER_GENE = Path("data/aligned/alignment_qc_per_gene.csv")
QC_PER_SEQ  = Path("data/aligned/alignment_qc_per_seq.csv")
QC_MISSING  = Path("data/aligned/missing_taxa_report.csv")

ALL_SPECIES: set[str] = {
    "Papilio_machaon", "Papilio_xuthus", "Papilio_glaucus",
    "Papilio_polytes", "Papilio_bianor", "Pieris_rapae",
    "Pieris_napi", "Gonepteryx_rhamni", "Delias_pasithoe",
    "Eurema_hecabe", "Danaus_plexippus", "Vanessa_indica",
    "Vanessa_cardui", "Junonia_almana", "Melitaea_cinxia",
    "Lycaena_phlaeas", "Arhopala_japonica", "Curetis_bulis",
    "Ampittia_dioscorides", "Lerema_accius",
}

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


# ─── Sequence-level metrics ───────────────────────────────────────────────────

def seq_gc_content(seq: str) -> float:
    """GC content ignoring gaps."""
    bases = [c for c in seq.upper() if c in "ACGT"]
    if not bases:
        return 0.0
    return (bases.count("G") + bases.count("C")) / len(bases)


def seq_gap_fraction(seq: str) -> float:
    return seq.count("-") / max(len(seq), 1)


def extract_species(record_id: str) -> str:
    """Pull first token (species name) from IDs like 'Papilio_machaon|ACC|COI'."""
    return record_id.split("|")[0]


# ─── Alignment-level metrics ──────────────────────────────────────────────────

def site_statistics(alignment: MultipleSeqAlignment) -> dict[str, int]:
    """
    Scan each column of the alignment to count:
      conserved, variable, parsimony_informative sites.
    """
    aln_len = alignment.get_alignment_length()
    conserved = variable = parsimony_inf = 0

    for i in range(aln_len):
        col   = alignment[:, i].upper()
        bases = [c for c in col if c not in ("-", "?", "N", "X")]
        unique = set(bases)

        if len(unique) <= 1:
            conserved += 1
        else:
            variable += 1
            freq = {b: bases.count(b) for b in unique}
            if sum(1 for f in freq.values() if f >= 2) >= 2:
                parsimony_inf += 1

    return {
        "conserved_sites":             conserved,
        "variable_sites":              variable,
        "parsimony_informative_sites": parsimony_inf,
    }


# ─── Per-gene QC ─────────────────────────────────────────────────────────────

def qc_gene(gene: str) -> tuple[dict, list[dict], list[str]]:
    """
    Run QC on a single gene alignment.

    Returns:
        gene_stats  – dict with aggregate metrics
        seq_rows    – list of per-sequence metric dicts
        missing     – list of species absent from this alignment
    """
    aln_path = ALIGNED_DIR / f"{gene}_aligned.fasta"
    empty_gene = {
        "gene": gene, "n_sequences": 0, "alignment_length": 0,
        "overall_gap_fraction": "N/A", "mean_seq_gap_fraction": "N/A",
        "mean_gc_content": "N/A", "conserved_sites": "N/A",
        "variable_sites": "N/A", "parsimony_informative_sites": "N/A",
    }

    if not aln_path.exists():
        log.warning("Aligned file not found: %s", aln_path)
        return empty_gene, [], list(ALL_SPECIES)

    try:
        alignment = AlignIO.read(str(aln_path), "fasta")
    except Exception as exc:
        log.error("Cannot parse %s: %s", aln_path, exc)
        return empty_gene, [], list(ALL_SPECIES)

    n_seq   = len(alignment)
    aln_len = alignment.get_alignment_length()

    # Per-sequence metrics
    seq_rows: list[dict] = []
    gap_fracs: list[float] = []
    gc_contents: list[float] = []
    present_species: set[str] = set()

    for rec in alignment:
        seq_str = str(rec.seq)
        gf  = seq_gap_fraction(seq_str)
        gc  = seq_gc_content(seq_str)
        sp  = extract_species(rec.id)
        gap_fracs.append(gf)
        gc_contents.append(gc)
        present_species.add(sp)

        seq_rows.append({
            "gene":          gene,
            "sequence_id":   rec.id,
            "species":       sp,
            "seq_length":    aln_len,
            "gap_fraction":  round(gf, 4),
            "gc_content":    round(gc, 4),
        })

    # Overall gap fraction
    all_chars   = n_seq * aln_len
    total_gaps  = sum(str(r.seq).count("-") for r in alignment)
    overall_gap = total_gaps / all_chars if all_chars > 0 else 0.0

    site_stats = site_statistics(alignment)

    gene_stats = {
        "gene":                        gene,
        "n_sequences":                 n_seq,
        "alignment_length":            aln_len,
        "overall_gap_fraction":        round(overall_gap, 4),
        "mean_seq_gap_fraction":       round(float(np.mean(gap_fracs)), 4),
        "mean_gc_content":             round(float(np.mean(gc_contents)), 4),
        **site_stats,
    }

    missing = sorted(ALL_SPECIES - present_species)

    log.info(
        "%s: %d seqs, %d bp, gap=%.2f%%, conserved=%d, variable=%d, PI=%d",
        gene, n_seq, aln_len,
        overall_gap * 100,
        site_stats["conserved_sites"],
        site_stats["variable_sites"],
        site_stats["parsimony_informative_sites"],
    )
    if missing:
        log.warning("  Missing in %s: %s", gene, ", ".join(missing))

    return gene_stats, seq_rows, missing


# ─── Main ─────────────────────────────────────────────────────────────────────

def run() -> None:
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    all_gene_stats: list[dict] = []
    all_seq_rows:   list[dict] = []
    missing_rows:   list[dict] = []

    for gene in GENES:
        gene_stats, seq_rows, missing = qc_gene(gene)
        all_gene_stats.append(gene_stats)
        all_seq_rows.extend(seq_rows)
        for sp in missing:
            missing_rows.append({"gene": gene, "missing_species": sp})

    # Write per-gene CSV
    if all_gene_stats:
        with open(QC_PER_GENE, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(all_gene_stats[0].keys()))
            writer.writeheader()
            writer.writerows(all_gene_stats)
        log.info("Per-gene QC → %s", QC_PER_GENE)

    # Write per-sequence CSV
    if all_seq_rows:
        with open(QC_PER_SEQ, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(all_seq_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_seq_rows)
        log.info("Per-sequence QC → %s", QC_PER_SEQ)

    # Write missing-taxa report
    if missing_rows:
        with open(QC_MISSING, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["gene", "missing_species"])
            writer.writeheader()
            writer.writerows(missing_rows)
        log.info("Missing taxa report → %s", QC_MISSING)
    else:
        log.info("All species represented in all gene alignments.")

    log.info("Step 4 complete.")


if __name__ == "__main__":
    run()

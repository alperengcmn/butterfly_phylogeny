#!/usr/bin/env python3
"""
05_concatenate.py
-----------------
Concatenate aligned gene FASTA files into a supermatrix for phylogenetic analysis.

Workflow:
  1. Read aligned FASTA for each mitochondrial and nuclear marker.
  2. Map sequences by species name.
  3. Fill gaps for species missing in a gene with an all-'?' placeholder.
  4. Write the concatenated supermatrix FASTA → data/concatenated/supermatrix.fasta
  5. Write partition file (RAxML/IQ-TREE format) → data/concatenated/partitions.txt

Usage:
    python scripts/05_concatenate.py
"""

from __future__ import annotations

import logging
import sys
from collections import defaultdict
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ─── Constants ────────────────────────────────────────────────────────────────

GENES: list[str] = ["COI", "COII", "CytB", "ND5", "EF1a", "wingless"]

ALIGNED_DIR = Path("data/aligned/trimmed")
CONCAT_DIR  = Path("data/concatenated")

SUPERMATRIX_FILE = CONCAT_DIR / "supermatrix.fasta"
PARTITION_FILE   = CONCAT_DIR / "partitions.txt"
LOG_FILE         = Path("logs/05_concat.log")

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

def extract_species(record_id: str) -> str:
    """Extract species component from ID 'Species_name|ACC|Gene'."""
    return record_id.split("|")[0].removeprefix("_R_")


def load_gene_alignment(gene: str) -> dict[str, str]:
    """
    Load an aligned FASTA and return {species: sequence_str}.
    If multiple sequences exist per species, keeps the first.
    """
    aln_path = ALIGNED_DIR / f"{gene}_trimmed.fasta"
    species_seqs: dict[str, str] = {}

    if not aln_path.exists():
        log.warning("Alignment not found: %s", aln_path)
        return species_seqs

    for rec in SeqIO.parse(str(aln_path), "fasta"):
        sp = extract_species(rec.id)
        if sp not in species_seqs:
            species_seqs[sp] = str(rec.seq)
        else:
            log.debug("Duplicate species %s in %s — keeping first.", sp, gene)

    log.info("Loaded %s: %d species, alignment length %d bp",
             gene, len(species_seqs),
             len(next(iter(species_seqs.values()), "")))
    return species_seqs


# ─── Supermatrix construction ─────────────────────────────────────────────────

def build_supermatrix() -> tuple[list[SeqRecord], list[tuple[str, int, int]]]:
    """
    Build the concatenated supermatrix.

    Returns:
        records    – list of SeqRecords (one per species)
        partitions – list of (gene_name, start_1based, end_1based)
    """
    # Load all gene alignments
    gene_data:   dict[str, dict[str, str]] = {}   # gene → {species: seq}
    gene_lengths: dict[str, int] = {}

    for gene in GENES:
        seqs = load_gene_alignment(gene)
        gene_data[gene] = seqs
        # Infer alignment length from first available sequence
        lengths = {len(s) for s in seqs.values()}
        if len(lengths) > 1:
            log.warning("%s has sequences of unequal length: %s", gene, lengths)
        gene_lengths[gene] = max(lengths, default=0)

    # Union of all species present across any gene
    all_species: set[str] = set()
    for seqs in gene_data.values():
        all_species.update(seqs.keys())
    all_species_sorted = sorted(all_species)

    log.info("Total unique species across all genes: %d", len(all_species_sorted))

    # Build concatenated sequences
    species_concat: dict[str, list[str]] = defaultdict(list)
    partitions: list[tuple[str, int, int]] = []
    cursor = 1

    for gene in GENES:
        glen = gene_lengths.get(gene, 0)
        if glen == 0:
            log.warning("Skipping gene %s (no data).", gene)
            continue

        seqs = gene_data[gene]
        for sp in all_species_sorted:
            if sp in seqs:
                fragment = seqs[sp]
            else:
                log.warning("Species %s missing in %s — filling with '?'", sp, gene)
                fragment = "?" * glen
            species_concat[sp].append(fragment)

        partitions.append((gene, cursor, cursor + glen - 1))
        cursor += glen

    # Build SeqRecord list
    records: list[SeqRecord] = []
    for sp in all_species_sorted:
        concat_seq = "".join(species_concat[sp])
        rec = SeqRecord(Seq(concat_seq), id=sp, description="")
        records.append(rec)

    return records, partitions


# ─── Main ─────────────────────────────────────────────────────────────────────

def run() -> None:
    CONCAT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    records, partitions = build_supermatrix()

    if not records:
        log.error("No sequences in supermatrix — aborting.")
        return

    # Write supermatrix FASTA
    SeqIO.write(records, str(SUPERMATRIX_FILE), "fasta")
    total_len = len(records[0].seq)
    log.info(
        "Supermatrix written → %s  (%d species × %d bp)",
        SUPERMATRIX_FILE, len(records), total_len,
    )

    # Write partitions.txt  (RAxML/IQ-TREE Nexus-lite format)
    with open(PARTITION_FILE, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for gene, start, end in partitions:
            fh.write(f"  charset {gene} = {start}-{end};\n")
        fh.write("end;\n")

    log.info("Partitions written → %s", PARTITION_FILE)

    # Print summary table to log
    log.info("─── Partition summary ───")
    for gene, start, end in partitions:
        log.info("  %-6s  %5d – %5d  (%d bp)", gene, start, end, end - start + 1)
    log.info("  TOTAL               %5d bp", total_len)

    log.info("Step 5 complete.")


if __name__ == "__main__":
    run()

#!/usr/bin/env python3
"""
03_align_sequences.py
---------------------
Multiple sequence alignment of cleaned FASTA files using MAFFT.

For each gene (COI, COII, CytB, ND5):
  - Runs MAFFT (--auto strategy, adjustdirection for strand issues).
  - Saves aligned FASTA to data/aligned/.
  - Computes per-gene alignment statistics.
  - Exports alignment_stats.csv.

Requires MAFFT installed and available on PATH.
  macOS: conda install -c bioconda mafft
         or: brew install mafft

Usage:
    python scripts/03_align_sequences.py
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

# ─── Constants ────────────────────────────────────────────────────────────────

GENES: list[str] = ["COI", "COII", "CytB", "ND5"]

CLEAN_DIR   = Path("data/cleaned")
ALIGNED_DIR = Path("data/aligned")
STATS_FILE  = Path("data/aligned/alignment_stats.csv")
LOG_FILE    = Path("logs/03_align.log")

# MAFFT parameters
MAFFT_ARGS: list[str] = [
    "--auto",           # auto-select strategy based on data size
    "--adjustdirection",# correct for reverse complement sequences
    "--thread", "-1",   # use all available CPU cores
    "--quiet",          # suppress progress output
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


# ─── Helpers ──────────────────────────────────────────────────────────────────

def check_mafft() -> str:
    """Return MAFFT executable path or raise RuntimeError."""
    path = shutil.which("mafft")
    if path is None:
        raise RuntimeError(
            "MAFFT not found on PATH.\n"
            "Install with:  conda install -c bioconda mafft\n"
            "           or: brew install mafft"
        )
    log.info("MAFFT found: %s", path)
    return path


def run_mafft(input_fasta: Path, output_fasta: Path, mafft_bin: str) -> None:
    """Run MAFFT and write the aligned output to *output_fasta*."""
    cmd = [mafft_bin] + MAFFT_ARGS + [str(input_fasta)]
    log.info("Running: %s", " ".join(cmd))

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        log.error("MAFFT stderr:\n%s", result.stderr)
        raise RuntimeError(f"MAFFT failed for {input_fasta.name} (exit {result.returncode})")

    output_fasta.write_text(result.stdout, encoding="utf-8")
    log.info("Alignment written → %s", output_fasta)


# ─── Alignment statistics ─────────────────────────────────────────────────────

def compute_alignment_stats(alignment: MultipleSeqAlignment, gene: str) -> dict:
    """
    Calculate basic alignment statistics.

    Returns a dict with:
      gene, n_sequences, alignment_length,
      gap_fraction, conserved_sites, variable_sites,
      parsimony_informative_sites
    """
    n_seq = len(alignment)
    aln_len = alignment.get_alignment_length()

    gap_count      = 0
    conserved      = 0
    variable       = 0
    parsimony_inf  = 0

    for col_i in range(aln_len):
        col = alignment[:, col_i].upper()  # string of characters at this column

        gap_count += col.count("-") + col.count("?")

        bases = [c for c in col if c not in ("-", "?", "N")]
        unique = set(bases)

        if len(unique) == 1:
            conserved += 1
        elif len(unique) > 1:
            variable += 1
            # Parsimony-informative: ≥2 different characters each present ≥2 times
            freq = {b: bases.count(b) for b in unique}
            if sum(1 for f in freq.values() if f >= 2) >= 2:
                parsimony_inf += 1

    total_cells  = n_seq * aln_len
    gap_fraction = gap_count / total_cells if total_cells > 0 else 0.0

    return {
        "gene":                       gene,
        "n_sequences":                n_seq,
        "alignment_length":           aln_len,
        "gap_fraction":               round(gap_fraction, 4),
        "conserved_sites":            conserved,
        "variable_sites":             variable,
        "parsimony_informative_sites": parsimony_inf,
    }


# ─── Main ─────────────────────────────────────────────────────────────────────

def run() -> None:
    ALIGNED_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    mafft_bin = check_mafft()
    all_stats: list[dict] = []

    for gene in GENES:
        in_path  = CLEAN_DIR   / f"{gene}_cleaned.fasta"
        out_path = ALIGNED_DIR / f"{gene}_aligned.fasta"

        if not in_path.exists():
            log.warning("Cleaned FASTA not found: %s — skipping.", in_path)
            continue

        # Count input sequences
        n_input = sum(1 for _ in SeqIO.parse(str(in_path), "fasta"))
        if n_input < 2:
            log.warning("%s: only %d sequence(s) — skipping alignment.", gene, n_input)
            continue

        log.info("Aligning %s (%d sequences) …", gene, n_input)

        try:
            run_mafft(in_path, out_path, mafft_bin)
        except RuntimeError as exc:
            log.error("Alignment failed for %s: %s", gene, exc)
            continue

        # Load alignment and compute stats
        try:
            alignment = AlignIO.read(str(out_path), "fasta")
            stats = compute_alignment_stats(alignment, gene)
            all_stats.append(stats)
            log.info(
                "  %s aligned: %d seqs × %d bp, gap=%.2f%%, conserved=%d, variable=%d",
                gene,
                stats["n_sequences"],
                stats["alignment_length"],
                stats["gap_fraction"] * 100,
                stats["conserved_sites"],
                stats["variable_sites"],
            )
        except Exception as exc:
            log.error("Could not read alignment for %s: %s", gene, exc)

    # Write stats CSV
    if all_stats:
        fieldnames = list(all_stats[0].keys())
        with open(STATS_FILE, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_stats)
        log.info("Alignment stats → %s", STATS_FILE)

    log.info("Step 3 complete.")


if __name__ == "__main__":
    run()

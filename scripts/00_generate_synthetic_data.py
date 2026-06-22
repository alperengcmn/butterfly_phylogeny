#!/usr/bin/env python3
"""
00_generate_synthetic_data.py
-----------------------------
Generates realistic synthetic mitochondrial gene sequences for 20 butterfly
species when NCBI is not accessible.

Simulation model:
  - Butterfly mtDNA is AT-rich (~75 % AT); reflected in base frequencies.
  - A root sequence is created for each gene.
  - Family-level mutations (substitution rate ~12 %) branch into 5 clades.
  - Species-level mutations (2-5 %) are added per terminal.
  - Realistic gene lengths (COI 1542 bp, COII 690 bp, CytB 1140 bp, ND5 1740 bp).
  - Produces metadata.csv compatible with Step 2 onwards.

Usage:
    python scripts/00_generate_synthetic_data.py
"""

from __future__ import annotations

import csv
import logging
import random
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# ─── Constants ────────────────────────────────────────────────────────────────

SEED = 42

GENE_LENGTHS: dict[str, int] = {
    "COI":  1542,
    "COII": 690,
    "CytB": 1140,
    "ND5":  1740,
}

# Butterfly mtDNA base frequencies (AT-rich)
BASE_FREQ: dict[str, float] = {"A": 0.38, "T": 0.37, "C": 0.14, "G": 0.11}

# Family → [species]
FAMILIES: dict[str, list[str]] = {
    "Papilionidae": [
        "Papilio_machaon", "Papilio_xuthus", "Papilio_glaucus",
        "Papilio_polytes", "Papilio_bianor",
    ],
    "Pieridae": [
        "Pieris_rapae", "Pieris_napi", "Gonepteryx_rhamni",
        "Delias_pasithoe", "Eurema_hecabe",
    ],
    "Nymphalidae": [
        "Danaus_plexippus", "Vanessa_indica", "Vanessa_cardui",
        "Junonia_almana", "Melitaea_cinxia",
    ],
    "Lycaenidae": [
        "Lycaena_phlaeas", "Arhopala_japonica", "Curetis_bulis",
    ],
    "Hesperiidae": [
        "Ampittia_dioscorides", "Lerema_accius",
    ],
}

# Family-level divergence from root (substitution fraction)
FAMILY_DIV   = 0.12
# Species-level divergence within family (substitution fraction)
SPECIES_DIV  = 0.035

RAW_DIR      = Path("data/raw")
LOG_FILE     = Path("logs/00_synthetic.log")
METADATA_FILE = RAW_DIR / "metadata.csv"

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


# ─── Sequence simulation ──────────────────────────────────────────────────────

def random_seq(length: int, rng: random.Random) -> str:
    """Generate a random DNA sequence with butterfly-like base composition."""
    bases = list(BASE_FREQ.keys())
    weights = list(BASE_FREQ.values())
    return "".join(rng.choices(bases, weights=weights, k=length))


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    """Apply random point substitutions at given rate."""
    bases = list("ACGT")
    result = list(seq)
    for i in range(len(result)):
        if rng.random() < rate:
            current = result[i]
            others  = [b for b in bases if b != current]
            result[i] = rng.choice(others)
    return "".join(result)


# ─── Main ─────────────────────────────────────────────────────────────────────

def run() -> None:
    rng = random.Random(SEED)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    metadata_rows: list[dict] = []
    gene_records:  dict[str, list[SeqRecord]] = {g: [] for g in GENE_LENGTHS}

    for gene, glen in GENE_LENGTHS.items():
        # Root ancestral sequence
        root_seq = random_seq(glen, rng)
        log.info("Gene %s root: %d bp", gene, glen)

        for family, species_list in FAMILIES.items():
            # Family clade ancestor
            family_seq = mutate(root_seq, FAMILY_DIV, rng)

            for species in species_list:
                # Species terminal
                sp_seq = mutate(family_seq, SPECIES_DIV, rng)
                acc    = f"SYN_{species[:3].upper()}{rng.randint(100000, 999999)}"
                rec_id = f"{species}|{acc}|{gene}"
                rec    = SeqRecord(Seq(sp_seq), id=rec_id, description="")
                gene_records[gene].append(rec)

                metadata_rows.append({
                    "species":         species.replace("_", " "),
                    "accession":       acc,
                    "gene":            gene,
                    "sequence_length": glen,
                    "ncbi_id":         acc,
                    "organism":        species.replace("_", " "),
                })

        # Write per-gene FASTA
        out = RAW_DIR / f"{gene}.fasta"
        SeqIO.write(gene_records[gene], str(out), "fasta")
        log.info("  Written %d seqs → %s", len(gene_records[gene]), out)

    # Write metadata CSV
    with open(METADATA_FILE, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(metadata_rows[0].keys()))
        writer.writeheader()
        writer.writerows(metadata_rows)
    log.info("Metadata → %s (%d rows)", METADATA_FILE, len(metadata_rows))
    log.info("Step 0 (synthetic data) complete.")


if __name__ == "__main__":
    run()

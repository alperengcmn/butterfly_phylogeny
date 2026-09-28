#!/usr/bin/env python3
"""
02_clean_sequences.py
---------------------
Quality-control and deduplication of raw FASTA files produced by Step 1.

Operations:
  1. Remove sequences shorter than a per-gene minimum length threshold.
  2. Remove exact-duplicate sequences (same sequence string, keep first).
  3. Standardise species names in FASTA headers.
  4. Write cleaned FASTA files to data/cleaned/.
  5. Export QC statistics to data/cleaned/qc_stats.csv.

Usage:
    python scripts/02_clean_sequences.py
"""

from __future__ import annotations

import csv
import hashlib
import logging
import re
import sys
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ─── Constants ────────────────────────────────────────────────────────────────

GENES: list[str] = ["COI", "COII", "CytB", "ND5", "EF1a", "wingless"]

# Minimum accepted sequence length per gene (bp)
MIN_LEN: dict[str, int] = {
    # Preserve the standard barcode-length threshold recommended for this dataset.
    "COI":  658,
    "COII": 600,
    "CytB": 900,
    "ND5":  1000,
    "EF1a": 500,
    "wingless": 200,
}

# Maximum tolerated ambiguous nucleotide fraction (N or ?)
MAX_AMBIGUOUS_FRAC: float = 0.05

RAW_DIR     = Path("data/raw")
CLEAN_DIR   = Path("data/cleaned")
QC_FILE     = Path("data/cleaned/qc_stats.csv")
LOG_FILE    = Path("logs/02_clean.log")
METADATA_FILE = Path("data/raw/metadata.csv")


@lru_cache(maxsize=1)
def _metadata_taxa_by_accession() -> dict[str, str]:
    """Map fetched accessions to the canonical query species in metadata.csv."""
    if not METADATA_FILE.exists():
        return {}
    with METADATA_FILE.open(newline="", encoding="utf-8") as handle:
        return {
            row["accession"]: row["species"].replace(" ", "_")
            for row in csv.DictReader(handle)
            if row.get("accession") and row.get("species")
        }

# Canonical species names → normalised form
SPECIES_NORM: dict[str, str] = {
    "papilio_machaon":       "Papilio_machaon",
    "papilio_xuthus":        "Papilio_xuthus",
    "papilio_glaucus":       "Papilio_glaucus",
    "papilio_polytes":       "Papilio_polytes",
    "papilio_bianor":        "Papilio_bianor",
    "pieris_rapae":          "Pieris_rapae",
    "pieris_napi":           "Pieris_napi",
    "gonepteryx_rhamni":     "Gonepteryx_rhamni",
    "aporia_crataegi":       "Aporia_crataegi",
    "eurema_hecabe":         "Eurema_hecabe",
    "danaus_plexippus":      "Danaus_plexippus",
    "vanessa_indica":        "Vanessa_indica",
    "vanessa_cardui":        "Vanessa_cardui",
    "junonia_almana":        "Junonia_almana",
    "melitaea_cinxia":       "Melitaea_cinxia",
    "lycaena_phlaeas":       "Lycaena_phlaeas",
    "plebejus_argus":        "Plebejus_argus",
    "curetis_bulis":         "Curetis_bulis",
    "ampittia_dioscorides":  "Ampittia_dioscorides",
    "ochlodes_venata":       "Ochlodes_venata",
    "parnara_guttata":       "Parnara_guttata",
    "heteropterus_morpheus": "Heteropterus_morpheus",
    "pyrgus_malvae":         "Pyrgus_malvae",
    "celaenorrhinus_maculosus": "Celaenorrhinus_maculosus",
    "ctenoptilum_vasava":    "Ctenoptilum_vasava",
    "notocrypta_curvifascia":"Notocrypta_curvifascia",
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


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _seq_hash(seq: str) -> str:
    """MD5 hash of a nucleotide sequence (upper-cased)."""
    return hashlib.md5(seq.upper().encode()).hexdigest()


def _ambiguous_fraction(seq: str) -> float:
    """Fraction of characters that are N or ?."""
    n = sum(1 for c in seq.upper() if c in ("N", "?", "-"))
    return n / max(len(seq), 1)


def _normalise_species(header_id: str) -> str:
    """
    Extract and normalise the species component of a FASTA ID.
    Expects IDs produced by Step 1: Species_name|AccessionNo|Gene
    Falls back gracefully if the format differs.
    """
    parts = header_id.split("|")
    if len(parts) > 1:
        canonical = _metadata_taxa_by_accession().get(parts[1])
        if canonical:
            return canonical
    raw_name = parts[0] if parts else header_id
    key = raw_name.lower().replace(" ", "_")
    return SPECIES_NORM.get(key, raw_name)


def iter_fasta(path: Path) -> Iterator[SeqRecord]:
    """Yield SeqRecords from a FASTA file; silently skip malformed entries."""
    try:
        yield from SeqIO.parse(str(path), "fasta")
    except Exception as exc:
        log.error("Could not parse %s: %s", path, exc)


# ─── Per-gene cleaning ────────────────────────────────────────────────────────

def clean_gene(gene: str) -> dict:
    """
    Clean a single gene FASTA.
    Returns a stats dict for the QC report.
    """
    in_path  = RAW_DIR  / f"{gene}.fasta"
    out_path = CLEAN_DIR / f"{gene}_cleaned.fasta"

    stats: dict = {
        "gene":            gene,
        "input_seqs":      0,
        "removed_short":   0,
        "removed_ambig":   0,
        "removed_dup":     0,
        "output_seqs":     0,
    }

    if not in_path.exists():
        log.warning("Input file not found: %s — skipping.", in_path)
        return stats

    min_len       = MIN_LEN.get(gene, 500)
    # Identical haplotypes in different species are valid observations; only
    # remove duplicate accessions within the same taxon.
    seen_hashes:  set[tuple[str, str]] = set()
    kept_records: list[SeqRecord] = []

    for rec in iter_fasta(in_path):
        stats["input_seqs"] += 1
        seq_str = str(rec.seq)

        # 1 — Length filter
        if len(seq_str) < min_len:
            log.debug("  REMOVED (short, %d bp): %s", len(seq_str), rec.id)
            stats["removed_short"] += 1
            continue

        # 2 — Ambiguity filter
        if _ambiguous_fraction(seq_str) > MAX_AMBIGUOUS_FRAC:
            log.debug("  REMOVED (ambiguous %.1f%%): %s",
                      _ambiguous_fraction(seq_str) * 100, rec.id)
            stats["removed_ambig"] += 1
            continue

        # 3 — Duplicate filter
        species_norm = _normalise_species(rec.id)
        h = _seq_hash(seq_str)
        duplicate_key = (species_norm, h)
        if duplicate_key in seen_hashes:
            log.debug("  REMOVED (duplicate): %s", rec.id)
            stats["removed_dup"] += 1
            continue
        seen_hashes.add(duplicate_key)

        # 4 — Normalise ID
        parts = rec.id.split("|")
        acc   = parts[1] if len(parts) > 1 else rec.id
        rec.id          = f"{species_norm}|{acc}|{gene}"
        rec.description = ""
        rec.name        = ""

        kept_records.append(rec)

    stats["output_seqs"] = len(kept_records)

    if kept_records:
        SeqIO.write(kept_records, str(out_path), "fasta")
        log.info(
            "%s: %d → %d sequences kept "
            "(-%d short, -%d ambig, -%d dup)",
            gene,
            stats["input_seqs"],
            stats["output_seqs"],
            stats["removed_short"],
            stats["removed_ambig"],
            stats["removed_dup"],
        )
    else:
        log.warning("%s: no sequences passed QC filters.", gene)

    return stats


# ─── Main ─────────────────────────────────────────────────────────────────────

def run() -> None:
    CLEAN_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    all_stats: list[dict] = []
    for gene in GENES:
        stats = clean_gene(gene)
        all_stats.append(stats)

    # Write QC report
    if all_stats:
        fieldnames = list(all_stats[0].keys())
        with open(QC_FILE, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_stats)
        log.info("QC stats written → %s", QC_FILE)

    # Summary
    total_in  = sum(s["input_seqs"]  for s in all_stats)
    total_out = sum(s["output_seqs"] for s in all_stats)
    log.info(
        "Step 2 complete. Total: %d input → %d output sequences (%.1f%% retained).",
        total_in,
        total_out,
        100 * total_out / max(total_in, 1),
    )


if __name__ == "__main__":
    run()

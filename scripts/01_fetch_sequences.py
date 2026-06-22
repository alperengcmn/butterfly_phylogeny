#!/usr/bin/env python3
"""
01_fetch_sequences.py
---------------------
Query NCBI Entrez for mitochondrial gene sequences of 20 butterfly species.
Downloads GenBank records, extracts COI / COII / CytB / ND5 coding sequences,
writes per-gene FASTA files, and produces a metadata CSV.

Usage:
    python scripts/01_fetch_sequences.py --email your@email.com
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# ─── Constants ────────────────────────────────────────────────────────────────

SPECIES: list[str] = [
    "Papilio machaon",
    "Papilio xuthus",
    "Papilio glaucus",
    "Papilio polytes",
    "Papilio bianor",
    "Pieris rapae",
    "Pieris napi",
    "Gonepteryx rhamni",
    "Delias pasithoe",
    "Eurema hecabe",
    "Danaus plexippus",
    "Vanessa indica",
    "Vanessa cardui",
    "Junonia almana",
    "Melitaea cinxia",
    "Lycaena phlaeas",
    "Arhopala japonica",
    "Curetis bulis",
    "Ampittia dioscorides",
    "Lerema accius",
]

GENES: list[str] = ["COI", "COII", "CytB", "ND5"]

# Alternative name patterns that NCBI uses in feature qualifiers
GENE_ALIASES: dict[str, list[str]] = {
    "COI":  ["COI", "COX1", "cox1", "COXI", "cytochrome oxidase subunit I",
              "cytochrome c oxidase subunit I"],
    "COII": ["COII", "COX2", "cox2", "COXII", "cytochrome oxidase subunit II",
              "cytochrome c oxidase subunit II"],
    "CytB": ["CytB", "CYTB", "cytb", "cytochrome b", "cob"],
    "ND5":  ["ND5", "nad5", "NADH5", "NADH dehydrogenase subunit 5"],
}

RAW_DIR        = Path("data/raw")
FASTA_DIR      = Path("data/raw")
METADATA_FILE  = Path("data/raw/metadata.csv")
LOG_FILE       = Path("logs/01_fetch.log")

NCBI_DELAY     = 0.4   # seconds between requests (NCBI policy: ≤3/s without API key)
MAX_RECORDS    = 5     # max GenBank hits per species×gene query

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


# ─── Entrez helpers ───────────────────────────────────────────────────────────

def search_ncbi(species: str, gene: str, retmax: int = MAX_RECORDS) -> list[str]:
    """Return a list of GenBank accession IDs for *species* + *gene*."""
    query = (
        f'"{species}"[Organism] AND ({gene}[Gene Name] OR {gene}[All Fields]) '
        f'AND mitochondrion[Filter] AND CDS[Feature Key]'
    )
    log.debug("Searching NCBI: %s", query)
    time.sleep(NCBI_DELAY)
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    ids: list[str] = record.get("IdList", [])
    log.info("  %s | %s → %d hit(s)", species, gene, len(ids))
    return ids


def fetch_genbank(ids: list[str]) -> list[SeqRecord]:
    """Download and parse GenBank records for the given list of NCBI IDs."""
    if not ids:
        return []
    time.sleep(NCBI_DELAY)
    handle = Entrez.efetch(
        db="nucleotide",
        id=",".join(ids),
        rettype="gb",
        retmode="text",
    )
    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    return records


# ─── Gene extraction ──────────────────────────────────────────────────────────

def _matches_gene(qualifiers: dict, gene_name: str) -> bool:
    """Check whether a CDS feature corresponds to *gene_name*."""
    aliases = [a.lower() for a in GENE_ALIASES.get(gene_name, [gene_name])]
    for key in ("gene", "product", "note"):
        for val in qualifiers.get(key, []):
            if any(alias in val.lower() for alias in aliases):
                return True
    return False


def extract_gene_from_record(
    record: SeqRecord,
    gene_name: str,
) -> Optional[SeqRecord]:
    """
    Scan all CDS features in *record* for a match to *gene_name*.
    Returns a trimmed SeqRecord or None.
    """
    for feature in record.features:
        if feature.type != "CDS":
            continue
        if _matches_gene(feature.qualifiers, gene_name):
            try:
                seq = feature.extract(record.seq)
            except Exception as exc:
                log.warning("Could not extract feature from %s: %s", record.id, exc)
                continue
            if len(seq) < 100:
                continue  # skip suspiciously short extractions
            species_clean = record.annotations.get("organism", record.id).replace(" ", "_")
            new_id = f"{species_clean}|{record.id}|{gene_name}"
            return SeqRecord(seq, id=new_id, description="")
    return None


# ─── Main workflow ────────────────────────────────────────────────────────────

def run(email: str) -> None:
    """Top-level orchestration: fetch → extract → write."""
    Entrez.email = email

    # Prepare output dirs
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Per-gene FASTA accumulator: gene → list[SeqRecord]
    gene_records: dict[str, list[SeqRecord]] = {g: [] for g in GENES}
    metadata_rows: list[dict] = []

    for species in SPECIES:
        log.info("Processing: %s", species)
        for gene in GENES:
            ids = search_ncbi(species, gene)
            if not ids:
                log.warning("  No hits for %s | %s", species, gene)
                continue

            gb_records = fetch_genbank(ids)
            found = False
            for gb in gb_records:
                extracted = extract_gene_from_record(gb, gene)
                if extracted is None:
                    continue

                gene_records[gene].append(extracted)
                metadata_rows.append(
                    {
                        "species":        species,
                        "accession":      gb.id,
                        "gene":           gene,
                        "sequence_length": len(extracted.seq),
                        "ncbi_id":        gb.id,
                        "organism":       gb.annotations.get("organism", ""),
                    }
                )
                log.info("    ✓ %s extracted (%d bp) from %s", gene, len(extracted.seq), gb.id)
                found = True
                break  # take the first good match per gene per species

            if not found:
                log.warning("  Gene %s not found in downloaded records for %s", gene, species)

    # Write per-gene FASTA files
    for gene, records in gene_records.items():
        out_path = FASTA_DIR / f"{gene}.fasta"
        if records:
            SeqIO.write(records, out_path, "fasta")
            log.info("Wrote %d sequences → %s", len(records), out_path)
        else:
            log.warning("No sequences collected for %s", gene)

    # Write metadata CSV
    if metadata_rows:
        df = pd.DataFrame(metadata_rows)
        df.to_csv(METADATA_FILE, index=False)
        log.info("Metadata written → %s (%d rows)", METADATA_FILE, len(df))
    else:
        log.error("No metadata collected — check NCBI queries and network access.")

    log.info("Step 1 complete.")


# ─── Entry point ──────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch butterfly mitochondrial gene sequences from NCBI.")
    parser.add_argument(
        "--email",
        required=True,
        help="E-mail address for NCBI Entrez (required by NCBI policy).",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=MAX_RECORDS,
        help=f"Max GenBank hits per query (default: {MAX_RECORDS}).",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=NCBI_DELAY,
        help=f"Delay between NCBI requests in seconds (default: {NCBI_DELAY}).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(email=args.email)

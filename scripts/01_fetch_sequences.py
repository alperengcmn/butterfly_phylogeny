#!/usr/bin/env python3
"""
01_fetch_sequences.py
---------------------
Query NCBI Entrez for mitochondrial and nuclear gene sequences of butterfly species.
Downloads GenBank records and extracts COI / COII / CytB / ND5 / EF1a / wingless,
writes per-gene FASTA files, and produces a metadata CSV.

Usage:
    python scripts/01_fetch_sequences.py --email your@email.com
"""

from __future__ import annotations

import argparse
import logging
import re
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
    "Aporia crataegi",
    "Eurema hecabe",
    "Danaus plexippus",
    "Vanessa indica",
    "Vanessa cardui",
    "Junonia almana",
    "Melitaea cinxia",
    "Lycaena phlaeas",
    "Plebejus argus",
    "Curetis bulis",
    "Ampittia dioscorides",
    "Ochlodes venata",
    "Parnara guttata",
    "Heteropterus morpheus",
    "Pyrgus malvae",
    "Celaenorrhinus maculosus",
    "Ctenoptilum vasava",
    "Notocrypta curvifascia",
]

MITOCHONDRIAL_GENES: list[str] = ["COI", "COII", "CytB", "ND5"]
NUCLEAR_GENES: list[str] = ["EF1a", "wingless"]
GENES: list[str] = MITOCHONDRIAL_GENES + NUCLEAR_GENES

# Alternative name patterns that NCBI uses in feature qualifiers
GENE_ALIASES: dict[str, list[str]] = {
    "COI":  ["COI", "COX1", "cox1", "COXI", "cytochrome oxidase subunit I",
              "cytochrome c oxidase subunit I"],
    "COII": ["COII", "COX2", "cox2", "COXII", "cytochrome oxidase subunit II",
              "cytochrome c oxidase subunit II"],
    "CytB": ["CytB", "CYTB", "cytb", "cytochrome b", "cob"],
    "ND5":  ["ND5", "nad5", "NADH5", "NADH dehydrogenase subunit 5"],
    "EF1a": ["EF1A", "EF-1alpha", "EF1-alpha", "elongation factor 1-alpha",
              "elongation factor 1 alpha"],
    "wingless": ["wingless", "Wnt1", "Wnt-1"],
}

RAW_DIR        = Path("data/raw")
FASTA_DIR      = Path("data/raw")
METADATA_FILE  = Path("data/raw/metadata.csv")
LOG_FILE       = Path("logs/01_fetch.log")

NCBI_DELAY     = 0.4   # seconds between requests (NCBI policy: ≤3/s without API key)
MAX_RECORDS    = 25    # search past short/poor records for a usable accession
MIN_FETCH_LENGTH = {
    "COI": 658, "COII": 600, "CytB": 900, "ND5": 1000,
    "EF1a": 500, "wingless": 300,
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


# ─── Entrez helpers ───────────────────────────────────────────────────────────

def search_ncbi(species: str, gene: str, retmax: int = MAX_RECORDS,
                delay: float = NCBI_DELAY) -> list[str]:
    """Return a list of GenBank accession IDs for *species* + *gene*."""
    aliases = GENE_ALIASES.get(gene, [gene])
    gene_query = " OR ".join(f'"{alias}"[All Fields]' for alias in aliases)
    query = f'"{species}"[Organism] AND ({gene_query}) AND CDS[Feature Key]'
    if gene in MITOCHONDRIAL_GENES:
        query += " AND mitochondrion[Filter]"
    log.debug("Searching NCBI: %s", query)
    time.sleep(delay)
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    ids: list[str] = record.get("IdList", [])
    log.info("  %s | %s → %d hit(s)", species, gene, len(ids))
    return ids


def fetch_genbank(ids: list[str], delay: float = NCBI_DELAY) -> list[SeqRecord]:
    """Download and parse GenBank records for the given list of NCBI IDs."""
    if not ids:
        return []
    time.sleep(delay)
    handle = Entrez.efetch(
        db="nucleotide",
        id=",".join(ids),
        rettype="gb",
        retmode="text",
    )
    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    return records


def enrich_mitochondrial_sequences(email: str, delay: float = NCBI_DELAY) -> None:
    """Reuse each taxon's selected mitochondrial records to recover every CDS.

    NCBI may return a short barcode for a COI query even when a complete
    mitogenome was found by the COII/CytB/ND5 queries. This pass extracts all
    four mitochondrial markers from those already downloaded accessions and
    keeps the longest qualifying sequence for each taxon and marker.
    """
    Entrez.email = email
    if not METADATA_FILE.exists():
        log.warning("Cannot enrich mitochondrial loci without metadata.csv")
        return

    metadata = pd.read_csv(METADATA_FILE)
    mito_meta = metadata[metadata["gene"].isin(MITOCHONDRIAL_GENES)]
    accessions_by_species = {
        species: sorted(set(group["accession"].astype(str)))
        for species, group in mito_meta.groupby("species")
    }
    current: dict[str, dict[str, SeqRecord]] = {gene: {} for gene in MITOCHONDRIAL_GENES}
    for gene in MITOCHONDRIAL_GENES:
        fasta = RAW_DIR / f"{gene}.fasta"
        if fasta.exists():
            current[gene] = {
                rec.id.split("|")[0]: rec for rec in SeqIO.parse(str(fasta), "fasta")
            }

    best: dict[tuple[str, str], tuple[int, SeqRecord, SeqRecord]] = {}
    for species, accessions in accessions_by_species.items():
        for gb in fetch_genbank(accessions, delay=delay):
            for gene in MITOCHONDRIAL_GENES:
                extracted = extract_gene_from_record(gb, gene)
                if extracted is None or len(extracted.seq) < MIN_FETCH_LENGTH[gene]:
                    continue
                key = (species, gene)
                old_length = best.get(key, (0, None, None))[0]
                current_record = current[gene].get(species.replace(" ", "_"))
                current_length = len(current_record.seq) if current_record else 0
                if len(extracted.seq) > max(old_length, current_length):
                    extracted.id = f"{species.replace(' ', '_')}|{gb.id}|{gene}"
                    extracted.description = ""
                    extracted.name = ""
                    best[key] = (len(extracted.seq), extracted, gb)

    replaced: dict[tuple[str, str], dict] = {}
    for (species, gene), (_, extracted, gb) in best.items():
        current[gene][species.replace(" ", "_")] = extracted
        replaced[(species, gene)] = {
            "species": species,
            "accession": gb.id,
            "gene": gene,
            "sequence_length": len(extracted.seq),
            "ncbi_id": gb.id,
            "organism": gb.annotations.get("organism", ""),
        }

    if not replaced:
        log.info("No longer mitochondrial CDSs found among downloaded records.")
        return
    for gene, records in current.items():
        SeqIO.write(list(records.values()), RAW_DIR / f"{gene}.fasta", "fasta")
    mito_meta = mito_meta[
        ~mito_meta.apply(lambda row: (row["species"], row["gene"]) in replaced, axis=1)
    ]
    metadata = pd.concat([metadata[~metadata["gene"].isin(MITOCHONDRIAL_GENES)],
                          mito_meta, pd.DataFrame(replaced.values())], ignore_index=True)
    metadata.to_csv(METADATA_FILE, index=False)
    log.info("Recovered longer mitochondrial CDSs for %d species/marker pairs.", len(replaced))


# ─── Gene extraction ──────────────────────────────────────────────────────────

def _matches_gene(qualifiers: dict, gene_name: str) -> bool:
    """Check whether a CDS feature corresponds to *gene_name*."""
    aliases = GENE_ALIASES.get(gene_name, [gene_name])
    for val in qualifiers.get("gene", []):
        token = re.sub(r"[^a-z0-9]", "", val.lower())
        if any(token == re.sub(r"[^a-z0-9]", "", alias.lower()) for alias in aliases):
            return True
    for key in ("product", "note"):
        for val in qualifiers.get(key, []):
            text = val.lower()
            for alias in aliases:
                alias = alias.lower()
                if len(alias) <= 5:
                    if re.search(rf"(?<![a-z0-9]){re.escape(alias)}(?![a-z0-9])", text):
                        return True
                elif alias in text:
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

def run(email: str, max_records: int = MAX_RECORDS,
        delay: float = NCBI_DELAY, only_missing: bool = False) -> None:
    """Top-level orchestration: fetch → extract → write."""
    Entrez.email = email

    # Prepare output dirs
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Per-gene FASTA accumulator: gene → list[SeqRecord]
    gene_records: dict[str, list[SeqRecord]] = {g: [] for g in GENES}
    existing_taxa: dict[str, set[str]] = {g: set() for g in GENES}
    metadata_rows: list[dict] = []
    if only_missing:
        for gene in GENES:
            raw_path = RAW_DIR / f"{gene}.fasta"
            if raw_path.exists():
                gene_records[gene] = list(SeqIO.parse(str(raw_path), "fasta"))
                existing_taxa[gene] = {rec.id.split("|")[0] for rec in gene_records[gene]}
        if METADATA_FILE.exists():
            metadata_rows = pd.read_csv(METADATA_FILE).to_dict(orient="records")

    for species in SPECIES:
        log.info("Processing: %s", species)
        for gene in GENES:
            query_taxon = species.replace(" ", "_")
            if only_missing and query_taxon in existing_taxa[gene]:
                continue
            ids = search_ncbi(species, gene, retmax=max_records, delay=delay)
            if not ids:
                log.warning("  No hits for %s | %s", species, gene)
                continue

            gb_records = fetch_genbank(ids, delay=delay)
            candidates: list[tuple[int, SeqRecord, SeqRecord]] = []
            for gb in gb_records:
                extracted = extract_gene_from_record(gb, gene)
                if extracted is None:
                    continue
                if len(extracted.seq) < MIN_FETCH_LENGTH[gene]:
                    log.info("    Skipping short %s record (%d bp) from %s",
                             gene, len(extracted.seq), gb.id)
                    continue

                # Keep the query taxon label, not a subspecies/synonym string
                # from GenBank, so all loci join on the same canonical taxon.
                extracted.id = f"{species.replace(' ', '_')}|{gb.id}|{gene}"
                extracted.description = ""
                extracted.name = ""

                candidates.append((len(extracted.seq), extracted, gb))

            if candidates:
                _, extracted, gb = max(candidates, key=lambda item: item[0])
                gene_records[gene].append(extracted)
                existing_taxa[gene].add(query_taxon)
                metadata_rows.append({
                    "species": species,
                    "accession": gb.id,
                    "gene": gene,
                    "sequence_length": len(extracted.seq),
                    "ncbi_id": gb.id,
                    "organism": gb.annotations.get("organism", ""),
                })
                log.info("    ✓ Selected longest %s record (%d bp) from %s",
                         gene, len(extracted.seq), gb.id)
            else:
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
        enrich_mitochondrial_sequences(email=email, delay=delay)
    else:
        log.error("No metadata collected — check NCBI queries and network access.")

    log.info("Step 1 complete.")


# ─── Entry point ──────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch mitochondrial and nuclear butterfly markers from NCBI.")
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
    parser.add_argument(
        "--only-missing", action="store_true",
        help="Keep current raw data and query only species/marker combinations without a sequence.",
    )
    parser.add_argument(
        "--enrich-mitochondrial", action="store_true",
        help="Reuse current mitochondrial accessions to recover longer CDSs for all four mtDNA markers.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.enrich_mitochondrial:
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        enrich_mitochondrial_sequences(email=args.email, delay=args.delay)
    else:
        run(email=args.email, max_records=args.max_records, delay=args.delay,
            only_missing=args.only_missing)

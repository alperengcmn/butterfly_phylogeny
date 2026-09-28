#!/usr/bin/env python3
"""Infer separate mitochondrial and nuclear trees and compare shared splits."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ALIGNED_DIR = Path("data/aligned/trimmed")
OUT_DIR = Path("results/inheritance")
GROUPS = {
    "mitochondrial": ["COI", "COII", "CytB", "ND5"],
    "nuclear": ["EF1a", "wingless"],
}
BOOTSTRAP = 1000
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def species_id(record_id: str) -> str:
    return record_id.split("|")[0].removeprefix("_R_")


def load_alignment(gene: str) -> dict[str, str]:
    path = ALIGNED_DIR / f"{gene}_trimmed.fasta"
    if not path.exists():
        return {}
    alignment = AlignIO.read(str(path), "fasta")
    result: dict[str, str] = {}
    for record in alignment:
        result.setdefault(species_id(record.id), str(record.seq))
    return result


def build_group(group: str, genes: list[str]) -> tuple[Path, Path, list[str]]:
    data = {gene: load_alignment(gene) for gene in genes}
    data = {gene: seqs for gene, seqs in data.items() if seqs}
    if not data:
        raise RuntimeError(f"No alignments available for {group} markers")

    taxa = sorted({taxon for seqs in data.values() for taxon in seqs})
    lengths = {gene: len(next(iter(seqs.values()))) for gene, seqs in data.items()}
    records = []
    for taxon in taxa:
        sequence = "".join(data[gene].get(taxon, "?" * lengths[gene]) for gene in data)
        records.append(SeqRecord(Seq(sequence), id=taxon, description=""))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fasta = OUT_DIR / f"{group}.fasta"
    partitions = OUT_DIR / f"{group}.nex"
    SeqIO.write(records, fasta, "fasta")
    cursor = 1
    with partitions.open("w", encoding="utf-8") as handle:
        handle.write("#nexus\nbegin sets;\n")
        for gene, length in lengths.items():
            handle.write(f"  charset {gene} = {cursor}-{cursor + length - 1};\n")
            cursor += length
        handle.write("end;\n")
    return fasta, partitions, taxa


def infer_tree(group: str, fasta: Path, partitions: Path) -> Path:
    prefix = OUT_DIR / group
    treefile = Path(str(prefix) + ".treefile")
    contree = Path(str(prefix) + ".contree")
    treefile.unlink(missing_ok=True)
    contree.unlink(missing_ok=True)
    iqtree = shutil.which("iqtree3") or shutil.which("iqtree2") or shutil.which("iqtree")
    if iqtree:
        cmd = [iqtree, "-s", str(fasta), "-p", str(partitions), "-m", "MFP",
               "-B", str(BOOTSTRAP), "-T", "AUTO", "--prefix", str(prefix), "--redo"]
        subprocess.run(cmd, check=True)
        return treefile

    log.warning("IQ-TREE unavailable; making an NJ tree without bootstrap support for %s", group)
    records = list(SeqIO.parse(str(fasta), "fasta"))
    alignment = MultipleSeqAlignment(records)
    distances = DistanceCalculator("identity").get_distance(alignment)
    tree = DistanceTreeConstructor().nj(distances)
    Phylo.write(tree, str(treefile), "newick")
    return treefile


def bipartitions(tree_path: Path, shared_taxa: set[str]) -> set[tuple[str, ...]]:
    tree = Phylo.read(str(tree_path), "newick")
    result: set[tuple[str, ...]] = set()
    for clade in tree.get_nonterminals():
        side = {tip.name for tip in clade.get_terminals()} & shared_taxa
        other = shared_taxa - side
        if len(side) < 2 or len(other) < 2:
            continue
        left, right = tuple(sorted(side)), tuple(sorted(other))
        result.add(min((left, right), key=lambda x: (len(x), x)))
    return result


def run() -> None:
    outputs: dict[str, tuple[Path, list[str]]] = {}
    for group, genes in GROUPS.items():
        fasta, partitions, taxa = build_group(group, genes)
        outputs[group] = (infer_tree(group, fasta, partitions), taxa)

    mt_tree, mt_taxa = outputs["mitochondrial"]
    nuc_tree, nuc_taxa = outputs["nuclear"]
    shared = set(mt_taxa) & set(nuc_taxa)
    mt_splits = bipartitions(mt_tree, shared)
    nuc_splits = bipartitions(nuc_tree, shared)
    common = mt_splits & nuc_splits
    union = mt_splits | nuc_splits
    similarity = 100 * len(common) / len(union) if union else float("nan")

    report = [
        "MITOCHONDRIAL VS NUCLEAR PHYLOGENY COMPARISON",
        f"Mitochondrial taxa: {len(mt_taxa)}",
        f"Nuclear taxa: {len(nuc_taxa)}",
        f"Shared taxa compared: {len(shared)}",
        f"Mitochondrial internal splits: {len(mt_splits)}",
        f"Nuclear internal splits: {len(nuc_splits)}",
        f"Shared splits: {len(common)}",
        f"Split similarity (shared splits / union): {similarity:.1f}%",
        "",
        "Different supported topologies may indicate mitonuclear discordance. This comparison",
        "is a screening result, not proof of introgression or heteroplasmy. Interpret node",
        "support, marker coverage, and taxon sampling together.",
    ]
    report_path = OUT_DIR / "comparison.txt"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    log.info("Inheritance comparison written to %s", report_path)


if __name__ == "__main__":
    run()

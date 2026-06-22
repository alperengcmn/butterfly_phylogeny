#!/usr/bin/env python3
"""
06_run_iqtree.py
----------------
Run IQ-TREE 2 on the concatenated supermatrix to infer a maximum-likelihood
phylogenetic tree with partition-aware model selection and ultrafast bootstraps.

IQ-TREE options used:
  -s  supermatrix.fasta
  -p  partitions.txt      (partition model)
  -m  MFP                 (ModelFinder Plus – select best model per partition)
  -B  1000                (1000 ultrafast bootstrap replicates)
  -T  AUTO                (auto-select number of CPU threads)
  --prefix results/tree/butterfly
  --redo                  (overwrite previous run if present)

Requires IQ-TREE 2 on PATH.
  macOS M1/M2: conda install -c bioconda iqtree

Usage:
    python scripts/06_run_iqtree.py
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path

# ─── Constants ────────────────────────────────────────────────────────────────

CONCAT_DIR   = Path("data/concatenated")
RESULTS_DIR  = Path("results/tree")
LOG_FILE     = Path("logs/06_iqtree.log")

SUPERMATRIX  = CONCAT_DIR / "supermatrix.fasta"
PARTITIONS   = CONCAT_DIR / "partitions.txt"
PREFIX       = RESULTS_DIR / "butterfly"

BOOTSTRAP    = 1000
THREADS      = "AUTO"

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

def check_iqtree() -> str | None:
    """Find IQ-TREE executable; prefer 'iqtree2' then 'iqtree'. Returns None if not found."""
    for name in ("iqtree2", "iqtree"):
        path = shutil.which(name)
        if path:
            log.info("IQ-TREE found: %s", path)
            return path
    log.warning("IQ-TREE not found on PATH — falling back to Biopython NJ tree.")
    return None


def parse_best_model(log_path: Path) -> str:
    """Extract the best-fit model from the IQ-TREE log file."""
    if not log_path.exists():
        return "unknown"
    text = log_path.read_text(encoding="utf-8", errors="replace")
    # Pattern for partition model: "Best-fit model: GTR+F+I+G4 chosen…"
    match = re.search(r"Best-fit model:\s+(\S+)", text)
    if match:
        return match.group(1)
    return "see log"


def parse_bootstrap_summary(treefile: Path) -> dict[str, float]:
    """
    Read the ML tree and compute simple bootstrap statistics.
    Returns dict with min, max, mean support values.
    """
    from Bio import Phylo
    import io

    if not treefile.exists():
        return {}

    try:
        tree = Phylo.read(str(treefile), "newick")
        supports = [
            float(c.confidence)
            for c in tree.find_clades()
            if c.confidence is not None
        ]
        if not supports:
            return {}
        return {
            "min_bootstrap":  min(supports),
            "max_bootstrap":  max(supports),
            "mean_bootstrap": sum(supports) / len(supports),
            "n_nodes":        len(supports),
        }
    except Exception as exc:
        log.warning("Could not parse tree for bootstrap summary: %s", exc)
        return {}


# ─── IQ-TREE run ─────────────────────────────────────────────────────────────

def run_iqtree(iqtree_bin: str) -> None:
    """Build and execute the IQ-TREE command."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = [
        iqtree_bin,
        "-s",  str(SUPERMATRIX),
        "-p",  str(PARTITIONS),
        "-m",  "MFP",
        "-B",  str(BOOTSTRAP),
        "-T",  THREADS,
        "--prefix", str(PREFIX),
        "--redo",
    ]

    log.info("IQ-TREE command:\n  %s", " ".join(cmd))

    result = subprocess.run(
        cmd,
        capture_output=False,   # let IQ-TREE print to terminal in real-time
        text=True,
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(f"IQ-TREE exited with code {result.returncode}")

    log.info("IQ-TREE finished successfully.")


# ─── Main ─────────────────────────────────────────────────────────────────────

def run_nj_fallback() -> None:
    """
    Biopython NJ tree fallback when IQ-TREE is not available.
    Builds a neighbour-joining tree from the supermatrix and writes a Newick file
    compatible with Steps 8 and 9.
    """
    from Bio import SeqIO, AlignIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio import Phylo

    log.info("Building NJ tree with Biopython (IQ-TREE fallback) …")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load alignment
    records = list(SeqIO.parse(str(SUPERMATRIX), "fasta"))
    alignment = MultipleSeqAlignment(records)

    calculator   = DistanceCalculator("identity")
    dm           = calculator.get_distance(alignment)
    constructor  = DistanceTreeConstructor()
    tree         = constructor.nj(dm)

    # Write as .treefile and .contree (same content, both expected by downstream)
    for suffix in (".treefile", ".contree"):
        out = Path(str(PREFIX) + suffix)
        Phylo.write(tree, str(out), "newick")
        log.info("NJ tree written → %s", out)

    # Write a stub .log so report generator doesn't crash
    stub_log = Path(str(PREFIX) + ".log")
    stub_log.write_text("Best-fit model: NJ (Biopython fallback — IQ-TREE not installed)\n")
    log.info("NJ fallback complete. Note: install IQ-TREE for ML tree with bootstrap support.")


def run() -> None:
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    if not SUPERMATRIX.exists():
        log.error("Supermatrix not found: %s — run Step 5 first.", SUPERMATRIX)
        sys.exit(1)
    if not PARTITIONS.exists():
        log.error("Partition file not found: %s — run Step 5 first.", PARTITIONS)
        sys.exit(1)

    iqtree_bin = check_iqtree()

    if iqtree_bin is None:
        run_nj_fallback()
        return

    run_iqtree(iqtree_bin)

    # Post-run summary
    iqtree_log  = Path(str(PREFIX) + ".log")
    treefile    = Path(str(PREFIX) + ".treefile")
    iqtree_best = Path(str(PREFIX) + ".best_model.nex")

    best_model = parse_best_model(iqtree_log)
    log.info("Best-fit model: %s", best_model)

    bs_stats = parse_bootstrap_summary(treefile)
    if bs_stats:
        log.info(
            "Bootstrap support — min: %.1f, max: %.1f, mean: %.1f (n=%d nodes)",
            bs_stats["min_bootstrap"],
            bs_stats["max_bootstrap"],
            bs_stats["mean_bootstrap"],
            bs_stats["n_nodes"],
        )

    log.info("Output files in: %s/", RESULTS_DIR)
    for f in sorted(RESULTS_DIR.iterdir()):
        log.info("  %s", f.name)

    log.info("Step 6 complete.")


if __name__ == "__main__":
    run()

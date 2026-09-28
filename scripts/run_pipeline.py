#!/usr/bin/env python3
"""
run_pipeline.py
---------------
Master orchestrator for the Butterfly Phylogenomics pipeline.
Runs all 11 steps in order with logging, timing, and graceful error handling.

Usage:
    python scripts/run_pipeline.py --email your@email.com

    # Run from a specific step (e.g. re-run from alignment):
    python scripts/run_pipeline.py --email your@email.com --start-step 3

    # Run only specific steps:
    python scripts/run_pipeline.py --email your@email.com --steps 7 8 9
"""

from __future__ import annotations

import argparse
import importlib
import logging
import sys
import time
from pathlib import Path
from typing import Callable

# ─── Constants ────────────────────────────────────────────────────────────────

LOG_FILE = Path("logs/pipeline.log")

PIPELINE_STEPS: list[tuple[int, str, str]] = [
    (1, "01_fetch_sequences",   "NCBI data collection"),
    (2, "02_clean_sequences",   "Sequence QC & deduplication"),
    (3, "03_align_sequences",   "MAFFT multiple sequence alignment"),
    (4, "04_alignment_qc",      "Alignment quality control"),
    (5, "05_concatenate",       "Supermatrix construction"),
    (6, "06_run_iqtree",        "IQ-TREE phylogenetic analysis"),
    (7, "07_genetic_distances", "Pairwise genetic distances"),
    (8, "08_visualize_tree",    "Tree visualisation"),
    (9, "11_compare_inheritance", "Mitochondrial vs nuclear trees"),
    (10, "09_generate_report",   "Report generation"),
    (11, "10_summary_figure",    "Summary figure generation"),
]

# ─── Logging ──────────────────────────────────────────────────────────────────

LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_FILE, mode="w"),
    ],
)
log = logging.getLogger(__name__)

# ─── Step runner ──────────────────────────────────────────────────────────────

def run_step(step_num: int, module_name: str, description: str, email: str) -> bool:
    """
    Dynamically import and execute a pipeline step module.
    Returns True on success, False on failure.
    """
    log.info("")
    log.info("━" * 60)
    log.info("  STEP %d: %s", step_num, description.upper())
    log.info("━" * 60)

    # Add scripts/ to sys.path so importlib can find the modules
    scripts_dir = str(Path(__file__).parent)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)

    t0 = time.perf_counter()
    try:
        mod = importlib.import_module(module_name)

        # Step 1 needs an email argument; others take none
        run_fn: Callable = getattr(mod, "run")
        import inspect
        sig = inspect.signature(run_fn)
        if "email" in sig.parameters:
            run_fn(email=email)
        else:
            run_fn()

        elapsed = time.perf_counter() - t0
        log.info("  ✓ Step %d completed in %.1f s", step_num, elapsed)
        return True

    except Exception as exc:
        elapsed = time.perf_counter() - t0
        log.error("  ✗ Step %d FAILED after %.1f s: %s", step_num, elapsed, exc, exc_info=True)
        return False


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Butterfly Phylogenomics — full pipeline runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full run
  python scripts/run_pipeline.py --email you@example.com

  # Start from step 3 (alignment)
  python scripts/run_pipeline.py --email you@example.com --start-step 3

  # Run only the final comparison, report, and figure steps
  python scripts/run_pipeline.py --email you@example.com --steps 9 10 11

  # Stop on first failure (default: continue)
  python scripts/run_pipeline.py --email you@example.com --stop-on-error
""",
    )
    parser.add_argument("--email",        required=True, help="E-mail for NCBI Entrez.")
    parser.add_argument("--start-step",   type=int, default=1,  help="Start from step N (1–11).")
    parser.add_argument("--steps",        type=int, nargs="+",  help="Run only these step numbers.")
    parser.add_argument("--stop-on-error",action="store_true",  help="Abort pipeline on first failure.")
    return parser.parse_args()


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    wall_t0 = time.perf_counter()
    log.info("═" * 60)
    log.info("  BUTTERFLY PHYLOGENOMICS PIPELINE")
    log.info("  %s", time.strftime("%Y-%m-%d %H:%M:%S"))
    log.info("═" * 60)

    # Determine which steps to run
    if args.steps:
        steps_to_run = [s for s in PIPELINE_STEPS if s[0] in args.steps]
    else:
        steps_to_run = [s for s in PIPELINE_STEPS if s[0] >= args.start_step]

    if not steps_to_run:
        log.error("No steps selected — check --steps / --start-step arguments.")
        sys.exit(1)

    log.info("Steps to run: %s", [s[0] for s in steps_to_run])

    results: dict[int, bool] = {}
    for step_num, module_name, description in steps_to_run:
        success = run_step(step_num, module_name, description, email=args.email)
        results[step_num] = success
        if not success and args.stop_on_error:
            log.error("Pipeline aborted at Step %d (--stop-on-error).", step_num)
            break

    # Final summary
    elapsed_total = time.perf_counter() - wall_t0
    log.info("")
    log.info("═" * 60)
    log.info("  PIPELINE SUMMARY  (total: %.1f s)", elapsed_total)
    log.info("═" * 60)
    for step_num, module_name, description in steps_to_run:
        status  = "✓" if results.get(step_num, False) else "✗"
        ran     = step_num in results
        skipped = "" if ran else "  [SKIPPED]"
        log.info("  %s  Step %d: %s%s", status, step_num, description, skipped)

    n_failed = sum(1 for ok in results.values() if not ok)
    if n_failed == 0:
        log.info("")
        log.info("All steps completed successfully.")
    else:
        log.warning("")
        log.warning("%d step(s) failed — check logs/ for details.", n_failed)
        sys.exit(1)


if __name__ == "__main__":
    main()

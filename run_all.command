#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Butterfly Phylogenomics — Full Pipeline Runner
# ─────────────────────────────────────────────────────────────────────────────

EMAIL="alperen8490@gmail.com"
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_DIR"

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  BUTTERFLY PHYLOGENOMICS PIPELINE"
echo "  Directory: $PROJECT_DIR"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Activate conda base if available
CONDA_INIT=""
for CONDA_SH in ~/miniforge3/etc/profile.d/conda.sh \
                ~/opt/anaconda3/etc/profile.d/conda.sh \
                ~/anaconda3/etc/profile.d/conda.sh \
                /opt/homebrew/anaconda3/etc/profile.d/conda.sh; do
  if [ -f "$CONDA_SH" ]; then
    source "$CONDA_SH"
    conda activate base 2>/dev/null || true
    CONDA_INIT="$CONDA_SH"
    echo "  ✓ conda initialised from $CONDA_SH"
    break
  fi
done

# ── 1. Python dependencies ────────────────────────────────────────────────────
echo "▶ Installing Python dependencies..."
pip3 install --quiet --break-system-packages biopython numpy pandas matplotlib seaborn scipy 2>/dev/null \
  || pip3 install --quiet biopython numpy pandas matplotlib seaborn scipy 2>/dev/null \
  || true
echo "  ✓ Python packages ready"

# ── 2. MAFFT ─────────────────────────────────────────────────────────────────
if command -v mafft &>/dev/null; then
  echo "  ✓ MAFFT already installed: $(mafft --version 2>&1 | head -1)"
else
  echo "▶ Installing MAFFT..."
  conda install -y -c bioconda mafft -q 2>/dev/null \
    || brew install mafft 2>/dev/null \
    || true
  command -v mafft &>/dev/null && echo "  ✓ MAFFT installed" || echo "  ⚠ MAFFT not found — alignment will use Biopython fallback"
fi

# ── 3. IQ-TREE ───────────────────────────────────────────────────────────────
if command -v iqtree2 &>/dev/null || command -v iqtree &>/dev/null; then
  IQBIN=$(command -v iqtree2 || command -v iqtree)
  echo "  ✓ IQ-TREE already installed: $IQBIN"
else
  echo "▶ Installing IQ-TREE via conda..."
  conda install -y -c bioconda iqtree -q 2>/dev/null || true

  if ! command -v iqtree2 &>/dev/null && ! command -v iqtree &>/dev/null; then
    echo "▶ conda failed, trying brew tap trust..."
    brew trust brewsci/bio 2>/dev/null || true
    brew install brewsci/bio/iqtree 2>/dev/null || true
  fi

  if command -v iqtree2 &>/dev/null || command -v iqtree &>/dev/null; then
    echo "  ✓ IQ-TREE installed"
  else
    echo "  ⚠ IQ-TREE not found — Step 6 will use Biopython NJ fallback"
  fi
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  RUNNING PIPELINE STEPS"
echo "════════════════════════════════════════════════════════════════"

run_step() {
  local STEP_NUM="$1"
  local STEP_DESC="$2"
  local SCRIPT="$3"
  shift 3
  echo ""
  echo "━━━ STEP $STEP_NUM: $STEP_DESC ━━━"
  if python3 "$SCRIPT" "$@"; then
    echo "  ✓ Step $STEP_NUM complete"
  else
    echo "  ✗ Step $STEP_NUM FAILED — check logs/"
    return 1
  fi
}

run_step 1  "Fetching sequences from NCBI"      scripts/01_fetch_sequences.py --email "$EMAIL"
run_step 2  "Cleaning sequences"                 scripts/02_clean_sequences.py
run_step 3  "Multiple sequence alignment"        scripts/03_align_sequences.py
run_step 4  "Alignment quality control"          scripts/04_alignment_qc.py
run_step 5  "Building supermatrix"               scripts/05_concatenate.py
run_step 6  "Phylogenetic analysis (IQ-TREE)"    scripts/06_run_iqtree.py
run_step 7  "Pairwise genetic distances"         scripts/07_genetic_distances.py
run_step 8  "Tree visualisation"                 scripts/08_visualize_tree.py
run_step 9  "Generating reports"                 scripts/09_generate_report.py

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  ✅  PIPELINE COMPLETE"
echo "  Results → results/"
echo "  Figures → figures/"
echo "════════════════════════════════════════════════════════════════"
echo ""
read -p "Press Enter to close..."

<div align="center">

# 🦋 Butterfly Phylogenomics

### End-to-end reproducible phylogenetic analysis of 20 butterfly species  
### using four mitochondrial genes — COI · COII · CytB · ND5

[![CI](https://github.com/alperen8490/butterfly-phylogeny/actions/workflows/ci.yml/badge.svg)](https://github.com/alperen8490/butterfly-phylogeny/actions)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![IQ-TREE 2](https://img.shields.io/badge/IQ--TREE-2.x-orange.svg)](http://www.iqtree.org/)
[![MAFFT](https://img.shields.io/badge/MAFFT-7.x-purple.svg)](https://mafft.cbrc.jp/alignment/software/)

</div>

---

## ✨ Overview

This project reconstructs the phylogenetic relationships of **20 butterfly species** across 5 families using a fully automated, reproducible bioinformatics pipeline. From raw NCBI accessions to publication-quality figures — every step is scripted, logged, and version-controlled.

| Stat | Value |
|------|-------|
| 🦋 Species | 20 across 5 families |
| 🧬 Genes | COI, COII, CytB, ND5 |
| 📏 Supermatrix | 5,181 bp |
| 🌲 Tree method | Maximum Likelihood (IQ-TREE 2) |
| 🔁 Bootstrap | 1,000 ultrafast replicates |
| 📊 Best model | GTR+F · TIM2+F+I+G4 · GTR+F+I+G4 · TIM+F+G4 |
| ✅ Mean BS support | 73.5% (7 nodes ≥ 95%) |

---

## 🌿 Species Studied

<table>
<tr>
  <th>Family</th>
  <th>Species</th>
</tr>
<tr>
  <td>🔴 <b>Papilionidae</b></td>
  <td><i>Papilio machaon · P. xuthus · P. glaucus · P. polytes · P. bianor</i></td>
</tr>
<tr>
  <td>🔵 <b>Pieridae</b></td>
  <td><i>Pieris rapae · P. napi · Gonepteryx rhamni · Delias pasithoe · Eurema hecabe</i></td>
</tr>
<tr>
  <td>🟢 <b>Nymphalidae</b></td>
  <td><i>Danaus plexippus · Vanessa indica · V. cardui · Junonia almana · Melitaea cinxia</i></td>
</tr>
<tr>
  <td>🟠 <b>Lycaenidae</b></td>
  <td><i>Lycaena phlaeas · Arhopala japonica · Curetis bulis</i></td>
</tr>
<tr>
  <td>🟣 <b>Hesperiidae</b></td>
  <td><i>Ampittia dioscorides · Lerema accius</i></td>
</tr>
</table>

---

## 📁 Project Structure

```
butterfly_phylogeny/
│
├── 📂 data/
│   ├── raw/                # NCBI GenBank sequences + metadata.csv
│   ├── cleaned/            # QC-filtered FASTA + qc_stats.csv
│   ├── aligned/            # MAFFT alignments + alignment QC reports
│   └── concatenated/       # supermatrix.fasta + partitions.txt
│
├── 📂 scripts/
│   ├── 00_generate_synthetic_data.py  ← CI test data (no internet needed)
│   ├── 01_fetch_sequences.py          ← NCBI Entrez download
│   ├── 02_clean_sequences.py          ← QC & deduplication
│   ├── 03_align_sequences.py          ← MAFFT alignment
│   ├── 04_alignment_qc.py             ← Gap%, conserved/variable sites
│   ├── 05_concatenate.py              ← Supermatrix + partitions.txt
│   ├── 06_run_iqtree.py               ← IQ-TREE ML tree
│   ├── 07_genetic_distances.py        ← p-distance matrix + heatmap
│   ├── 08_visualize_tree.py           ← PNG/PDF tree figures
│   ├── 09_generate_report.py          ← Automated text reports
│   ├── 10_summary_figure.py           ← 4-panel publication figure
│   └── run_pipeline.py                ← Master orchestrator
│
├── 📂 results/
│   ├── tree/               # IQ-TREE outputs (.treefile, .log, .contree …)
│   ├── distance_matrix.csv
│   ├── jc_distance_matrix.csv
│   ├── summary_figure.png  ← 4-panel summary (tree + heatmap + stats)
│   ├── summary_figure.pdf
│   ├── bulgular.txt        ← Full findings & thesis section (Turkish)
│   ├── project_summary.txt
│   └── evolutionary_summary.txt
│
├── 📂 figures/
│   ├── tree_biopython.png / .pdf
│   └── distance_heatmap.png
│
├── 📂 logs/                # Per-step run logs
├── 📂 .github/workflows/   # GitHub Actions CI
│
├── environment.yml         # Conda environment (reproducible)
├── requirements.txt        # pip dependencies
├── pyproject.toml          # Python packaging metadata
├── CITATION.cff            # How to cite this project
├── run_all.command         # macOS double-click runner
└── README.md
```

---

## ⚙️ Installation

### macOS M1 / M2 (recommended)

```bash
# 1 — Clone the repository
git clone https://github.com/alperen8490/butterfly-phylogeny.git
cd butterfly-phylogeny/butterfly_phylogeny

# 2 — Create conda environment (installs MAFFT + IQ-TREE automatically)
conda env create -f environment.yml
conda activate butterfly_phylogeny

# 3 — Verify tools
mafft --version
iqtree2 --version
python -c "import Bio; print('Biopython', Bio.__version__)"
```

### Linux / Ubuntu

```bash
sudo apt-get install -y mafft
# IQ-TREE via conda:
conda install -c bioconda iqtree
pip install -r requirements.txt
```

### pip only (all platforms)

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# Install MAFFT and IQ-TREE via system package manager separately
```

---

## 🚀 Running the Pipeline

### Option A — Double-click (macOS)

```
Double-click:  run_all.command
```

Terminal opens automatically, installs missing tools, and runs all 9 steps.

### Option B — Full pipeline via Python

```bash
python scripts/run_pipeline.py --email your@email.com
```

### Option C — Step-by-step

```bash
# Step 1: Download sequences from NCBI
python scripts/01_fetch_sequences.py --email your@email.com

# Steps 2–5: Clean → Align → QC → Supermatrix
python scripts/02_clean_sequences.py
python scripts/03_align_sequences.py
python scripts/04_alignment_qc.py
python scripts/05_concatenate.py

# Step 6: IQ-TREE ML phylogeny
python scripts/06_run_iqtree.py

# Steps 7–9: Distances → Figures → Reports
python scripts/07_genetic_distances.py
python scripts/08_visualize_tree.py
python scripts/09_generate_report.py

# Step 10: 4-panel summary figure
python scripts/10_summary_figure.py
```

### Resume from a specific step

```bash
# Re-run from alignment onwards
python scripts/run_pipeline.py --email your@email.com --start-step 3

# Run only reporting steps
python scripts/run_pipeline.py --email your@email.com --steps 9 10
```

---

## 🔬 Pipeline Overview

| Step | Script | Tool | Output |
|------|--------|------|--------|
| 0 | `00_generate_synthetic_data.py` | Python | Synthetic test FASTA (CI) |
| 1 | `01_fetch_sequences.py` | Biopython Entrez | `data/raw/*.fasta`, `metadata.csv` |
| 2 | `02_clean_sequences.py` | Python | `data/cleaned/*_cleaned.fasta`, `qc_stats.csv` |
| 3 | `03_align_sequences.py` | **MAFFT 7.x** | `data/aligned/*_aligned.fasta` |
| 4 | `04_alignment_qc.py` | NumPy | Alignment QC CSVs |
| 5 | `05_concatenate.py` | Biopython | `supermatrix.fasta`, `partitions.txt` |
| 6 | `06_run_iqtree.py` | **IQ-TREE 2** | ML tree, bootstrap, best model |
| 7 | `07_genetic_distances.py` | NumPy/Pandas | `distance_matrix.csv`, heatmap |
| 8 | `08_visualize_tree.py` | Biopython Phylo / ETE3 | `tree_*.png/pdf` |
| 9 | `09_generate_report.py` | Pandas | Summary & evolutionary reports |
| 10 | `10_summary_figure.py` | Matplotlib | `summary_figure.png/pdf` |

---

## 📊 Key Results

### Phylogenetic Tree

The ML tree was inferred from a **5,181 bp supermatrix** with per-partition model selection:

| Gene | Model | Rationale |
|------|-------|-----------|
| COI | GTR+F | General reversible, empirical frequencies |
| COII | TIM2+F+I+G4 | Invariable sites + Gamma rate variation |
| CytB | GTR+F+I+G4 | Full GTR with rate heterogeneity |
| ND5 | TIM+F+G4 | Transition-specific rates + Gamma |

### Genetic Distances

| Metric | Value |
|--------|-------|
| Mean pairwise p-distance | **13.5%** |
| Closest pair | *Papilio machaon* ↔ *P. xuthus* — **6.45%** |
| Most divergent pair | *Delias pasithoe* ↔ *Vanessa indica* — **22.6%** |

### Bootstrap Support

```
≥ 95%  ████████████████████  7 nodes   (36.8%)
70–94% ████████              2 nodes   (10.5%)
50–69% ████████████████████████████████  8 nodes (42.1%)
< 50%  ████████              2 nodes   (10.5%)
                             Mean: 73.5%
```

---

## 📄 Output Files

| File | Description |
|------|-------------|
| `results/tree/butterfly.treefile` | ML tree (Newick) |
| `results/tree/butterfly.contree` | Consensus tree with bootstrap |
| `results/tree/butterfly.iqtree` | Full IQ-TREE report |
| `results/distance_matrix.csv` | N×N p-distance matrix |
| `results/jc_distance_matrix.csv` | Jukes-Cantor corrected distances |
| `results/summary_figure.png` | 4-panel publication figure |
| `results/bulgular.txt` | Full findings & thesis section |
| `results/project_summary.txt` | Pipeline-wide statistics |
| `results/evolutionary_summary.txt` | Biological interpretation |
| `figures/tree_biopython.png` | Standalone tree (PNG) |
| `figures/distance_heatmap.png` | Distance heatmap |

---

## 🔄 Reproducibility

- All randomness in IQ-TREE is seeded internally — output is identical for the same input.
- NCBI sequences can change over time; raw FASTA files should be archived for long-term reproducibility.
- Python dependencies are fixed in `requirements.txt`; conda versions are pinned in `environment.yml`.
- GitHub Actions CI runs the full pipeline on every push using synthetic data (no NCBI needed).

---

## ❓ Troubleshooting

| Problem | Solution |
|---------|----------|
| `mafft: command not found` | `conda install -c bioconda mafft` |
| `iqtree2: command not found` | `conda install -c bioconda iqtree` |
| Few sequences from NCBI | Try `--max-records 10` in Step 1 |
| Low COI coverage | Lower min-length threshold in `02_clean_sequences.py` to 658 bp |
| ETE3 import error | `pip install ete3` (optional; PNG/PDF work without it) |
| `brewsci/bio` untrusted | `brew trust brewsci/bio && brew install brewsci/bio/iqtree` |

---

## 📖 Citation

If you use this pipeline, please cite:

```bibtex
@software{butterfly_phylogenomics_2026,
  title  = {Butterfly Phylogenomics Pipeline},
  author = {Alperen},
  year   = {2026},
  url    = {https://github.com/alperen8490/butterfly-phylogeny}
}
```

**Tools used:**

- **IQ-TREE 2** — Minh et al. (2020) *Mol. Biol. Evol.* 37:1530–1534  
- **MAFFT** — Katoh & Standley (2013) *Mol. Biol. Evol.* 30:772–780  
- **ModelFinder** — Kalyaanamoorthy et al. (2017) *Nat. Methods* 14:587–589  
- **Biopython** — Cock et al. (2009) *Bioinformatics* 25:1422–1423  
- **ETE3** — Huerta-Cepas et al. (2016) *Mol. Biol. Evol.* 33:1635–1638  

---

## 📜 License

MIT © 2026 Alperen

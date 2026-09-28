<div align="center">

# 🦋 Butterfly Phylogenomics

### Reproducible phylogenetic analysis of 26 target butterfly species
### using four mitochondrial and two nuclear markers

[![CI](https://github.com/alperengcmn/butterfly_phylogeny/actions/workflows/ci.yml/badge.svg)](https://github.com/alperengcmn/butterfly_phylogeny/actions)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![IQ-TREE](https://img.shields.io/badge/IQ--TREE-2%20%7C%203-orange.svg)](http://www.iqtree.org/)
[![MAFFT](https://img.shields.io/badge/MAFFT-7.x-purple.svg)](https://mafft.cbrc.jp/alignment/software/)

</div>

---

## ✨ Overview

This project targets **26 butterfly species** across five families, using four mitochondrial genes and two nuclear markers. Separate mitochondrial and nuclear trees provide an independent screen for mitonuclear discordance.

| Stat | Value |
|------|-------|
| 🦋 Taxa | 26 targets; 26 in the current matrix (8 Hesperiidae) |
| 🧬 Markers | COI, COII, CytB, ND5, EF1a, wingless |
| 📏 Supermatrix | 5,651 bp after complete-case trimming |
| 🌲 Tree method | Maximum Likelihood (IQ-TREE) |
| 🔁 Bootstrap | 1,000 ultrafast replicates |
| 📊 Best model | Selected per marker by IQ-TREE ModelFinder |
| ✅ Mean BS support | 92.6% across 23 internal nodes |

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
  <td><i>Pieris rapae · P. napi · Gonepteryx rhamni · Aporia crataegi · Eurema hecabe</i></td>
</tr>
<tr>
  <td>🟢 <b>Nymphalidae</b></td>
  <td><i>Danaus plexippus · Vanessa indica · V. cardui · Junonia almana · Melitaea cinxia</i></td>
</tr>
<tr>
  <td>🟠 <b>Lycaenidae</b></td>
  <td><i>Lycaena phlaeas · Plebejus argus · Curetis bulis</i></td>
</tr>
<tr>
  <td>🟣 <b>Hesperiidae</b></td>
  <td><i>Ampittia dioscorides · Ochlodes venata · Parnara guttata · Heteropterus morpheus · Pyrgus malvae · Celaenorrhinus maculosus · Ctenoptilum vasava · Notocrypta curvifascia</i></td>
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
│   ├── aligned/            # MAFFT alignments + QC reports
│   │   └── trimmed/        # 100%-occupancy per-locus FASTA files used downstream
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
│   ├── 11_compare_inheritance.py       ← mtDNA/nuclear comparison
│   └── run_pipeline.py                ← Master orchestrator
│
├── 📂 results/
│   ├── tree/               # IQ-TREE outputs (.treefile, .log, .contree …)
│   ├── distance_matrix.csv
│   ├── jc_distance_matrix.csv
│   ├── summary_figure.png  ← 4-panel summary (tree + heatmap + stats)
│   ├── inheritance/        # Separate mtDNA/nuclear trees and comparison
│   ├── summary_figure.pdf
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
├── bulgular.txt            # Current findings (Turkish)
└── README.md
```

---

## ⚙️ Installation

### macOS M1 / M2 (recommended)

```bash
# 1 — Clone the repository
git clone https://github.com/alperengcmn/butterfly_phylogeny.git
cd butterfly_phylogeny

# 2 — Create conda environment (installs MAFFT + IQ-TREE automatically)
conda env create -f environment.yml
conda activate butterfly_phylogeny

# 3 — Verify tools
mafft --version
iqtree --version  # or iqtree2, depending on the installed release
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

Terminal opens automatically, installs missing tools, and runs all 11 steps.

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

# Steps 7–11: Distances → Trees → Reports → Figures
python scripts/07_genetic_distances.py
python scripts/08_visualize_tree.py
python scripts/11_compare_inheritance.py
python scripts/09_generate_report.py
python scripts/10_summary_figure.py
```

### Resume from a specific step

```bash
# Re-run from alignment onwards
python scripts/run_pipeline.py --email your@email.com --start-step 3

# Run only reporting steps
python scripts/run_pipeline.py --email your@email.com --steps 9 10 11
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
| 6 | `06_run_iqtree.py` | **IQ-TREE** | ML tree, bootstrap, best model |
| 7 | `07_genetic_distances.py` | NumPy/Pandas | `distance_matrix.csv`, heatmap |
| 8 | `08_visualize_tree.py` | Biopython Phylo / ETE3 | `tree_*.png/pdf` |
| 9 | `11_compare_inheritance.py` | IQ-TREE / Biopython | Separate mitochondrial and nuclear trees + split comparison |
| 10 | `09_generate_report.py` | Pandas | Summary & evolutionary reports |
| 11 | `10_summary_figure.py` | Matplotlib | `summary_figure.png/pdf` |

---

## 📊 Key Results

### Phylogenetic Tree

The current run contains all 26 target taxa and 5,651 positions across six markers. Each locus was trimmed to columns present in all 26 taxa (100% site occupancy), so the concatenated matrix has no missing taxa or gap characters. This strict complete-case rule avoids gaps caused by partial overlap, while discarding 16.5–30.6% of some locus alignments. The exact matrix and model results are recorded in `results/project_summary.txt` and `results/tree/butterfly.iqtree`.

IQ-TREE ModelFinder Plus selects a model separately for each of COI, COII, CytB, ND5, EF1a, and wingless during the run.

### Genetic Distances

| Metric | Value |
|--------|-------|
| Mean pairwise p-distance | 0.1600 |
| Lowest pairwise p-distance | 0.0713 — *Pieris napi* / *P. rapae* |
| Highest pairwise p-distance | 0.2527 — *Papilio xuthus* / *Pieris rapae* |

### Bootstrap Support

| Statistic | Current result |
|-----------|----------------|
| Ultrafast replicates | 1,000 |
| Mean support | 92.6% across 23 internal nodes |
| Nodes with support ≥70% | 22 / 23 |
| Nodes with support ≥95% | 13 / 23 |

All 26 taxa have sequences for all six markers; the missing taxon/marker report is empty. Per-locus complete-case lengths are COI 1,512 bp, COII 614 bp, CytB 1,119 bp, ND5 1,714 bp, EF1a 354 bp, and wingless 338 bp. Separate mitochondrial and nuclear trees share 4 of 42 distinct internal splits (9.5% split similarity); this is a screening comparison and does not establish introgression or heteroplasmy. IQ-TREE composition/model assumptions and the limited nuclear sequence retained after trimming remain relevant when interpreting support. See [the current findings](bulgular.txt) for details.

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
| `results/project_summary.txt` | Pipeline-wide statistics |
| `results/evolutionary_summary.txt` | Biological interpretation |
| `results/inheritance/comparison.txt` | Mitochondrial/nuclear split comparison |
| `results/inheritance/mitochondrial.treefile` | Mitochondrial-marker ML tree |
| `results/inheritance/nuclear.treefile` | Nuclear-marker ML tree |
| `figures/tree_biopython.png` | Standalone tree (PNG) |
| `figures/distance_heatmap.png` | Distance heatmap |
| `bulgular.txt` | Updated Turkish findings and limitations |

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
| Low COI coverage | COI now retains sequences ≥658 bp; check `qc_stats.csv` for retained counts |
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
  url    = {https://github.com/alperengcmn/butterfly_phylogeny}
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

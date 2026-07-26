```markdown

# Evolutionary Deoptimization and Host-Mediated Genome Remodeling in Hepatitis A Virus (HAV)


## Overview
This repository contains the data, phylogenetic models, and interactive analysis pipelines for investigating evolutionary conflict and genome remodeling in **Hepatitis A Virus (HAV)**. 

Contrary to the traditional assumption of "evolutionary stasis" in acute hepatoviruses, this study demonstrates that perceived low substitution rates ($6.74 \times 10^{-4}$ substitutions/site/year) and effective population size ($N_e$) contractions are driven by pervasive purifying selection and host-mediated antiviral editing. By combining Bayesian phylodynamics with structural and dinucleotide motif analyses, we uncover how host deaminase activity (e.g., APOBEC3-mediated $C>T$ transitions) shapes viral codon deoptimization over time.

While the primary focus is on evolutionary mechanics, clinical datasets (such as Genotype III progression to Acute Liver Failure) are integrated as functional validation of genome remodeling.

---

## Repository Structure

```text
HAV_Evolution/
├── data/
│   ├── raw/                 # Unaligned sequences (.fasta), metadata, and raw BEAST logs (.log/.trees)*
│   └── processed/           # Curated alignments, thinned posterior trees, and skygrid exports
├── notebooks/               # Marimo interactive notebooks for statistical testing & figure generation
├── results/
│   ├── figures/             # Publication vector/raster graphics (.svg, .pdf, .png)
│   └── tables/              # Automated statistical reports and summary matrices (.txt, .csv)
├── .gitignore               # Excludes large MCMC posteriors (>100MB) and system cache
└── README.md                # Project documentation

```

**Note: Large posterior `.trees` files (>100 MB) are excluded via `.gitignore` to comply with GitHub file size limits and are archived separately.*

---

## Computational Workflow & Methodology

### 1. Bayesian Phylodynamics & Tree Annotation

* **Alignment & Filtering:** Sequences aligned via **MAFFT** and curated for temporal signal.
* **MCMC Sampling:** Evolutionary rates and demographic history modeled using **BEAST 1.x/2.x** under a Skygrid coalescent model.
* **Posterior Summarization:** MCMC runs converged with Effective Sample Sizes ($\text{ESS} > 200$). Posterior topologies are thinned (removing a 10% burn-in via `LogCombiner`) and summarized into Maximum Clade Credibility (MCC) trees using **TreeAnnotator** with common ancestor node heights.

### 2. Mutational Signatures & Dinucleotide Motifs

* **Deaminase Tracking:** Custom Python scripts scan pairwise sequence divergence against reference genomes (`NC_001489`) to quantify directional transitions, specifically isolating $C>T$ mutations indicative of APOBEC3 host defense mechanisms.
* **Motif Remodeling:** Dinucleotide Observed/Expected (O/E) ratios for CpG and GpA motifs are tracked across historical (<2018) and contemporary (>= 2018) epidemic cohorts using non-parametric Kruskal-Wallis and Dunn's post-hoc tests with Holm-Bonferroni correction.

### 3. Interactive Visualization & iTOL Mapping

* **Figure Generation:** Main manuscript visualizations (e.g., Figure 1 statistical summaries, regression models, violin distributions, and Skygrid overlays) are built natively using **Marimo**, **Seaborn**, and **Matplotlib**.
* **Phylogenetic Mapping:** Leaf-node mutational proportions are exported as structured `DATASET_SIMPLEBAR` formatting for direct overlay onto MCC trees in the **Interactive Tree Of Life (iTOL)**.

---

## Setup & Reproduction

### Prerequisites

Ensure you have Python 3.10+ installed along with the required scientific suites:

```bash
pip install pandas numpy scipy matplotlib seaborn scikit-posthocs biopython marimo

```

### Running the Interactive Notebooks

To interact with or regenerate the statistical plots and figures:

```bash
# Launch Marimo in edit/run mode
marimo edit notebooks/HAV_Evolutionary_Dynamics.py

```

### Generating iTOL Annotations

To recalculate $C>T$ transition proportions and output updated simple-bar mapping datasets for iTOL:

---

## Citation & Contact

If you adapt or utilize these pipelines for viral evolutionary analysis, please reference the upcoming publication (full citation pending review).

For technical inquiries or access to raw MCMC posterior tree distributions, please open an **Issue** in this repository.

```

---


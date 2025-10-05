# Aomine_Iijima_2025

# Repository Overview

This repository contains the analysis code used in our research paper.

All scripts and notebooks included here were used for processing, analyzing, and visualizing single-cell RNA-seq (scRNA-seq) and T-cell receptor sequencing (TCR-seq) data derived from human and mouse samples.

---

## File Descriptions

### scRNA-seq Analysis (Python and R Scripts)

- **`py_human_scRNAseq.ipynb`**  
  Code used for the analysis of **human scRNA-seq** datasets.

- **`py_mouse_snRNAseq.ipynb`**  
  Code used for the analysis of **mouse scRNA-seq and TCR-seq** datasets.

- **`r_human_scRNAseq_Conversion_to_h5ad.R`**  
  R script used to convert **human scRNA-seq data** into `.h5ad` format for downstream analysis with Python.

- **`yml_scRNA-seq_full_env.yml`**  
  Complete Python environment file used for all scRNA-seq and TCR-seq analyses (both human and mouse).  
  This file allows reproducible setup of the same computational environment via Conda.

---

### TCR-seq Analysis (R and Shell Scripts)

- **`r_mouse_TCRseq_Donut.R`**  
  R script for **TCR-seq analysis**, generating **donut plots**.

- **`r_mouse_TCRseq_Table.R`**  
  R script for **TCR-seq analysis**, generating **summary tables**.

- **`r_mouse_TCRseq_preanalysis.R`**  
  R script for **preprocessing** TCR-seq data prior to generating donut plots and tables.

- **`r_mouse_TCRseq_Ubuntu24.sh`**  
  Shell script for **raw TCR-seq data processing** on Ubuntu 24 environment.

---

## Usage

1. Clone this repository:
   ```bash
   git clone https://github.com/<username>/<repository>.git
   cd <repository>

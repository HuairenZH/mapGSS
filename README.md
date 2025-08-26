# mapGSS

## 1. Introduction

`mapGSS` is a computational tool designed to quantify the spatial specificity of gene expression in spatial transcriptomics datasets. It integrates gene spatial specificity scores (GSS) computed by gsMap with tissue or cell-type annotations, enabling systematic analysis and visualization of gene-specific spatial enrichment patterns across regions or cell types.

```
mapGSS
├── README.md
├── data
│   ├── Genes
│   │   ├── KCNQ3
│   │   │   ├── map_gene_KCNQ3.pdf
│   │   │   ├── map_gene_KCNQ3.png
│   │   │   ├── mouse_Adult_Mouse_brain_cell_bin_KCNQ3.csv
│   │   │   └── mouse_E16.5_E1S1_KCNQ3.csv
│   │   └──  SLC32A1
│   │       ├── map_gene_SLC32A1.pdf
│   │       ├── map_gene_SLC32A1.png
│   │       ├── mouse_Adult_Mouse_brain_cell_bin_SLC32A1.csv
│   │       └── mouse_E16.5_E1S1_SLC32A1.csv
│   └── ST
│       ├── E16.5_E1S1
│       │   ├── E16.5_E1S1_gene_marker_score.feather
│       │   └── mouse_E16.5_E1S1.txt
│       └── mouse_Adult_Mouse_brain_cell_bin
│           ├── mouse_Adult_Mouse_brain_cell_bin.txt
│           └── mouse_Adult_Mouse_brain_cell_bin_gene_marker_score.feather
├── pyproject.toml
└── src
    └── extract_gene_gss
        ├── __init__.py
        ├── cli.py
        ├── map_cli.py
        ├── map_gene.R
        └── map_gene.py

```



## 2. Key Features

- __Spatially-aware High-Resolution Trait Mapping__
- __Spatial Region Identification__
- __Putative Causal Genes Identification__

![Model Architecture](./data/Genes/SLC32A1/map_gene_SLC32A1.png)

## 3. Installation

**Before getting started, you will also need the R packages `ggplot2`, `dplyr`, `readr`, `ggpubr`, `ggridges`, `RColorBrewer`, `forcats`, `stringr`, `rlang`, `scales`, `patchwork`, and `extrafont`. Before running the code for the first time, you may need to execute the following commands in R to load the fonts:**

```R
library(extrafont)

# When running for the first time, import system fonts (this may take a few minutes)
font_import(prompt = FALSE)    # Do not show confirmation prompts; import all fonts by default
loadfonts()                    # Load the fonts to make them available for use
```



install from source:

```bash
conda create -n mapGSS python = 3.12.2
conda activate mapGSS
git clone https://github.com/HuairenZH/mapGSS
cd mapGSS
pip install -e .
```

Verify the installation by running the following command:

```bash
extract_gene_GSS --help
map_gene_GSS --help  
```

## 4. Usage

```bash
map_gene_GSS --gene SLC32A1
map_gene_GSS --gene SLC32A1
```



## 5. Citation

Please cite the paper and give us a STAR if you find mapGSS useful for your research.

# mapGSS

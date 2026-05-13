
# panNEC Single-Cell Atlas — Code & Analysis Repository

**Associated manuscript:** Aberrant brain-type neuronal programs in large-cell pancreatic neuroendocrine carcinoma  
**Zenodo DOI:** https://doi.org/10.5281/zenodo.20097468  
**GEO accession:** [GSE291816](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291816)

---

## Overview

This repository contains the analysis code and notebooks underlying the panNEC atlas. Fresh-frozen tumour tissue from 5 patients with high-grade large-cell pancreatic neuroendocrine carcinoma (panNEC) was profiled by single-nucleus RNA-sequencing (snRNA-seq), yielding 45,015 nuclei across different shared and unique cell states. Single-nucleus transcriptomics revealed two distinct neuroendocrine cell states: one highly proliferative with an aberrant brain-specific neuronal differentiation programme, and a stress-responsive state enriched for heat stress, hypoxia, and glycolysis.

Processed data (h5ad) & the list of 4,000 highly variable genes are available via Zenodo. Raw CellRanger count matrices are deposited at GEO (GSE291816).

---

## Contact

**Dr. Olivia Debnath** — olivia_debnath@dfci.harvard.edu  
**Dr. Hilmar Berger** — hilmar.berger@charite.de
**Dr. Naveed Ishaque** - naveed.ishaque@bih-charite.de

## Repository structure

```
panNEC_study2025/
├── Figure_1/                        # Main figure 1 notebooks
├── Figure_2/                        # Main figure 2 notebooks
├── Supplementary_figure_analyses/   # Notebooks per latest manuscript supplementary figure numbering
├── bioarxiv_2025_extended/          # Original bioRxiv extended data figures (archived)
│   ├── Extended_1/
│   ├── Extended_2/
│   ├── Extended_3/
│   └── Extended_4/
├── analysis/                        # Core atlas building pipeline
├── snRNAseq_preprocessing/          # Preprocessing pipeline (SoupX, alignment)
└── readme.md
```

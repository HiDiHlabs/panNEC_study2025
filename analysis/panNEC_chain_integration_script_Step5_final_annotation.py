"""
panNEC Single-Cell Atlas — Chain Integration Script (Step 5: Semifinal Annotation)
===============================================================================
Description:
    This script produces the semifinal cell annotations for the
    panNEC atlas. Starting from the subclustered h5ad (Step 4), it:
      1. Creates 'Cell states' — 12 fine-grained publication labels
      2. Creates 'Cell types' — 4 coarse-grained labels
      3. Creates 'PatientID' — clean de-identified patient labels
      4. Scores external NE/LCNEC/PDAC/pathway gene signatures
      5. Generates composition bar plots
      6. Saves the final annotated objects

Reference notebook:
    pNET_scRNAseq_hilmar_2022/notebooks_analysis/pNEC_updated_annotations_06082023.ipynb

Semifinal annotation layers (further refined in Figure notebooks):
    'Cell states'  — 12 fine-grained states (primary annotation)
    'Cell types'   — 4 coarse types: Neuroendocrine, Amphicrine acinar,
                     Stroma (normal), Immune
    'PatientID'    — P07, P018, P019, P021, P024 (acinar-like)

Author: Olivia Debnath (Current email: olivia_debnath@dfci.harvard.edu)

Associated manuscript:
    Aberrant brain-type neuronal programs in large-cell pancreatic
    neuroendocrine carcinoma
    Zenodo: https://zenodo.org/records/20097469
    GEO: GSE291816

Note: File paths are relative to the HPC project directory
    /dh-projects/ag-ishaque/analysis/debnatho/pNET_scRNAseq_hilmar_2022/
    Update paths accordingly before running.
"""

#── Libraries ──────────────────────────────────────────────────────────────────
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (4, 4)
plt.rcParams["axes.grid"] = False
sc.settings.set_figure_params(dpi=80, dpi_save=180, vector_friendly=True, transparent=True)

#── 1. Load subclustered object (Step 4 output) ───────────────────────────────
results = "./pNET_comparative_clustering/pNET_updated_annotations_180123.h5ad"
adata = sc.read_h5ad(results)

#Remove ambiguous L12 cluster (discussed in Step 3):
adata_sub = adata[adata.obs['res_1'] != '12']

#── 2. Cell states — fine-grained annotation (12 states) ─────────────────────
#First pass: intermediate labels
adata_sub.obs['cell_states'] = (
    adata_sub.obs["cell_type_semifinal_v4"].map(lambda x: {
        "Acinar_like_NET01":              "Amphicrine01 NE",
        "Acinar_like_NET02":              "Amphicrine01 NE",
        "BRCA1/FANCA_NETp":              "BRCA1/FANCA NEp",
        "BRCA1/FANCA_NETp (acinar like)":"BRCA1/FANCA NEp (amphicrine)",
        "Exo_glandular_mixed01":          "Amphicrine02 NE",
        "Exo_glandular_mixed02":          "Amphicrine03 NE",
        "HSPpos_NET":                     "HSP+ NE",
        "GIPR+ EEC like":                "GIPR-high EEC-like NE",
        "Smooth muscle like tumors (SMT)":"WNT-high EEC-like NE",
        "Stromal2":                       "Stromal (mixed)"
    }.get(x, x)).astype("category")
)

#Second pass: simplify per Philipp's feedback
adata_sub.obs['Cell states'] = (
    adata_sub.obs["cell_states"].map(lambda x: {
        "Amphicrine01 NE":           "Amphicrine acinar01",
        "Amphicrine02 NE":           "Amphicrine acinar02",
        "BRCA1/FANCA NEp":           "Neuroendocrine proliferating",
        "BRCA1/FANCA NEp (amphicrine)":"Amphicrine acinar proliferating",
        "Amphicrine03 NE":           "Amphicrine acinar03",
        "WNT-high EEC-like NE":      "Neuroendocrine EEC-like01",
        "HSP+ NE":                   "Neuroendocrine HSP+",
        "GIPR-high EEC-like NE":     "Neuroendocrine EEC-like02",
        "Stromal1":                  "Neuroendocrine stromal-like",
        "Stromal (mixed)":           "Stroma (normal)",
        "Lymphocytes (non-tumor)":   "Lymphocytes",
        "Macrophages (non-tumor)":   "Macrophages"
    }.get(x, x)).astype("category")
)

#Reorder for plotting:
adata_sub.obs['Cell states'] = adata_sub.obs['Cell states'].cat.reorder_categories([
    'Amphicrine acinar01', 'Amphicrine acinar02', 'Amphicrine acinar03',
    'Amphicrine acinar proliferating', 'Neuroendocrine EEC-like01',
    'Neuroendocrine EEC-like02', 'Neuroendocrine proliferating',
    'Neuroendocrine HSP+', 'Neuroendocrine stromal-like', 'Stroma (normal)',
    'Lymphocytes', 'Macrophages'
])

#Color palette:
adata_sub.uns['Cell states_colors'] = [
    '#b3ccff', '#cc99ff', '#e6ccff', '#3973ac', '#9999ff',
    '#00cccc', '#008080', '#004d4d', '#4dffb8', '#ff99cc', '#ff0000', '#ff704d'
]

sc.pl.umap(adata_sub, color=['Cell states'], save="_Cell_states_final_07082023.pdf")

#── 3. Cell types — coarse annotation (4 types) ───────────────────────────────
#Stromal1 = Stromal-like NE (aneuploid per InferCNV → classified as tumor)
adata_sub.obs['Cell types'] = (
    adata_sub.obs["cell_type_semifinal_v4"].map(lambda x: {
        "Acinar_like_NET01":              "Amphicrine acinar",
        "Acinar_like_NET02":             "Amphicrine acinar",
        "BRCA1/FANCA_NETp":             "Neuroendocrine",
        "BRCA1/FANCA_NETp (acinar like)":"Amphicrine acinar",
        "Exo_glandular_mixed01":         "Amphicrine acinar",
        "Exo_glandular_mixed02":         "Amphicrine acinar",
        "HSPpos_NET":                    "Neuroendocrine",
        "GIPR+ EEC like":               "Neuroendocrine",
        "Smooth muscle like tumors (SMT)":"Neuroendocrine",
        "Stromal1":                      "Neuroendocrine",
        "Stromal2":                      "Stroma (normal)",
        "Lymphocytes (non-tumor)":       "Immune",
        "Macrophages (non-tumor)":       "Immune"
    }.get(x, x)).astype("category")
)

adata_sub.obs['Cell types'] = adata_sub.obs['Cell types'].cat.reorder_categories([
    'Neuroendocrine', 'Amphicrine acinar', 'Stroma (normal)', 'Immune'
])

adata_sub.uns['Cell types_colors'] = ['#00cccc', '#e6ccff', '#ff99cc', '#991f00']

sc.pl.umap(adata_sub, color=['Cell types'], save="_celltypes_final_11072023.pdf")

#── 4. PatientID — clean de-identified labels ─────────────────────────────────
adata_sub.obs['PatientID'] = (
    adata_sub.obs["source_label"].map(lambda x: {
        "SP084_007": "P07",
        "SP084_018": "P018",
        "SP084_019": "P019",
        "SP084_021": "P021",
        "SP084_024": "P024 (acinar-like)"
    }.get(x, x)).astype("category")
)

sc.pl.umap(adata_sub, color=['PatientID', 'Cell states'], ncols=1)

#── 5. Patient composition bar plot ───────────────────────────────────────────
df = pd.crosstab(adata_sub.obs['PatientID'], adata_sub.obs['Cell states'],
                 normalize="index")

plt.rcParams["figure.figsize"] = (2, 6)
pl = df.plot(kind="bar", stacked=True, rot=90,
             color=adata_sub.uns['Cell states_colors'])
pl.grid(False)
pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('Cell_states_composition_v1_ExtdFig2C.pdf')

#── 6. Normalize for gene scoring ─────────────────────────────────────────────
adata_norm = adata_sub.raw.to_adata()
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)

#── 6. Save semifinal annotated object ────────────────────────────────────────
import os
os.makedirs("pnec_anndata_082023", exist_ok=True)

results_out = "./pnec_anndata_082023/pNEC_updated_annot_07082023.h5ad" #Input to 1_panNEC_analysis_Figure1.ipynb 
adata_sub.write(results_out)

print("Semifinal annotation layers:")
print(f"  Cell states : {list(adata_sub.obs['Cell states'].cat.categories)}")
print(f"  Cell types  : {list(adata_sub.obs['Cell types'].cat.categories)}")
print(f"  PatientID   : {list(adata_sub.obs['PatientID'].cat.categories)}")

#Note: downstream gene signature scoring (George LCNEC, Lazar, Alshalalfa,
#Metascape pathway modules) is performed in the Figure notebooks on adata_norm
#(normalized version of adata_sub). Those analyses do not modify cell annotations.

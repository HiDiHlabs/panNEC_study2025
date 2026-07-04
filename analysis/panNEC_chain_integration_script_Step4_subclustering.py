"""
panNEC Single-Cell Atlas — Chain Integration Script (Step 4: Subclustering & Final Annotation)
================================================================================================
Description:
    This script performs stromal subclustering, lymphocyte subclustering, and
    final cluster label refinement of the panNEC atlas. Starting from the
    filtered annotated h5ad (Step 3), it:
      1. Subclusters Stromal population into Stromal1 (shared) & Stromal2 (YAP/TAZ-high P7)
      2. Refines cluster labels based on December 2022 feedback (cell_type_semifinal_v4)
      3. Subclusters Lymphocytes into 3 subtypes (cell_type_semifinal_v5)
      4. Runs Wilcoxon marker detection on final labels
      5. Exports minimal obs metadata for Seurat/R handoff
      6. Saves the final annotated object

Reference notebooks:
    pNET_scRNAseq_hilmar_2022/notebooks_analysis/pNET_subclustering_11122022.ipynb
    pNET_scRNAseq_hilmar_2022/notebooks_analysis/pNET_manuscript_figure01_0123.ipynb
        (cluster refinement sections only; figure generation is in the Figure 1 notebook)

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
sc.settings.set_figure_params(dpi=80, dpi_save=150, facecolor='white')

#── 1. Load filtered annotated object (Step 3 output) ─────────────────────────
results = "./pNET_comparative_clustering/pNET_harmony_bbknn_filtered_021122.h5ad"
adata = sc.read_h5ad(results)
sc.pl.umap(adata, color=['source_label', 'cell_type_semifinal_v2'], ncols=1)

#── 2. EEC marker validation per patient ──────────────────────────────────────
#GIPR/TPH1 is robust in P19 (~60% GIPR+); P7 is ~10%
#For non-P19 patients, classification is more EEC-like or NET-like
sc.pl.umap(adata, color=['GIPR', 'TPH1'], ncols=1)

adata_eec = adata[adata.obs['cell_type_semifinal_v2'] == "GIPR_EEC_like"]
sc.pl.dotplot(adata_eec, ['GIPR', 'SYP', 'CHGA', 'CHGB', 'TPH1', 'TAC1'],
              groupby="sample_id_final")

for i in adata.obs['sample_id_final'].cat.categories:
    print(i)
    sc.pl.dotplot(adata[adata.obs['sample_id_final'] == i],
                  ['GIPR', 'SYP', 'CHGA', 'CHGB', 'TPH1', 'TAC1'],
                  groupby="cell_type_semifinal02", standard_scale="var", title=i)

#── 3. Stromal subclustering ──────────────────────────────────────────────────
#A portion of P7 stromal is YAP/TAZ-high with specific mesenchymal genes
#not found in other patients → subdivide into Stromal1 & Stromal2
sc.tl.leiden(adata, restrict_to=('cell_type_semifinal_v2', ['Stromal']),
             resolution=0.4, key_added='Stromal_sub')

sc.pl.umap(adata, color=['cell_type_semifinal_v2', 'Stromal_sub'], ncols=1)

adata.obs['cell_type_semifinal_v3'] = (
    adata.obs['Stromal_sub'].map(lambda x: {
        'Stromal,0': 'Stromal1',
        'Stromal,1': 'Stromal2'
    }.get(x, x)).astype("category")
)

adata.obs["cell_type_semifinal_v3"].value_counts()
sc.pl.umap(adata, color=['source_label', 'cell_type_semifinal_v3'], ncols=1)

#Stromal marker validation per patient:
for i in adata.obs['sample_id_final'].cat.categories:
    print(i)
    sc.pl.dotplot(adata[adata.obs['sample_id_final'] == i], [
        'COL1A2', 'CALD1', 'CASC15', 'FBXL7', 'ZFPM2', 'BICC1', 'SYNPO2',
        'ROBO2', 'DPYSL3', 'GPC6', 'COL6A3', 'PRKG1', 'DLG2', 'SOX5',
        'AEBP1', 'FN1', 'PDGFRA', 'YAP1', 'WWTR1'
    ], groupby="cell_type_semifinal_v3", standard_scale="var", title=i)

#── 4. Wilcoxon markers for v3 labels ────────────────────────────────────────
sc.tl.rank_genes_groups(adata, 'cell_type_semifinal_v3', method='wilcoxon')

sc.pl.rank_genes_groups_dotplot(adata, n_genes=6,
    values_to_plot="logfoldchanges", cmap='bwr',
    vmin=-4, vmax=4, min_logfoldchange=3,
    colorbar_title='LFC cut-off: 3')

#Export marker table:
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
     for group in groups
     for key in ['names', 'pvals', 'pvals_adj', 'logfoldchanges', 'scores']}
).head(50).to_csv("pNET_cluster_markers_Wilcoxon_12122022.csv")

#── 5. Final label refinement — cell_type_semifinal_v4 ────────────────────────
#Renamed based on Bertram/Katherina's December 2022 feedback:
#AcinarS_NET_like → Acinar_like_NET01 (relatively robust NET characteristics:
#    CERS4, CACNA1A etc.)
#AcinarS_like → Acinar_like_NET02
#Non-tumor populations explicitly labelled in brackets
adata.obs['cell_type_semifinal_v4'] = (
    adata.obs['cell_type_semifinal_v3'].map(lambda x: {
        'AcinarS_NET_like':              'Acinar_like_NET01',
        'AcinarS_like':                  'Acinar_like_NET02',
        'BRCA1/FANCA_prol01':            'BRCA1/FANCA_NETp',
        'BRCA1/FANCA_prol02_acinar_like':'BRCA1/FANCA_NETp (acinar like)',
        'HSP_pos':                       'HSPpos_NET',
        'Microglia_like/Macrophages':    'Macrophages (non-tumor)',
        'Lymphocytes':                   'Lymphocytes (non-tumor)',
        'SMC_like':                      'Smooth muscle like tumors (SMT)',
        'GIPR_EEC_like':                 'GIPR+ EEC like'
    }.get(x, x)).astype("category")
)

adata.obs["cell_type_semifinal_v4"].value_counts()
sc.pl.umap(adata, color=['cell_type_seafinal_v4'], ncols=1)

#Patient composition bar plot:
adata.obs['cell_type_semifinal_v4_ordered'] = (
    adata.obs['cell_type_semifinal_v4'].cat.reorder_categories([
        'SMT', 'BRCA1/FANCA_NETp (acinar like)', 'Exo_glandular_mixed01',
        'Exo_glandular_mixed02', 'Acinar_like_NET02', 'Acinar_like_NET01',
        'HSPpos_NET', 'BRCA1/FANCA_NETp', 'GIPR+ EEC like',
        'Stromal1', 'Stromal2', 'Lymphocytes (non-tumor)', 'Macrophages (non-tumor)'
    ], ordered=False)
)

df = pd.crosstab(adata.obs['cell_type_semifinal_v4_ordered'],
                 adata.obs['sample_id_final'], normalize="index")

plt.rcParams["figure.figsize"] = (6, 4)
pl = df.plot(kind="bar", stacked=True, rot=90)
pl.grid(False)
pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#── 6. Lymphocyte subclustering — cell_type_semifinal_v5 ──────────────────────
#Subdivide Lymphocytes (non-tumor) into 3 subtypes at low resolution
sc.tl.leiden(adata, restrict_to=('cell_type_semifinal_v4', ['Lymphocytes (non-tumor)']),
             resolution=0.3, key_added='Lympho_sub')

sc.pl.umap(adata, color=['Lympho_sub'], ncols=1)

adata.obs['cell_type_semifinal_v5'] = (
    adata.obs['Lympho_sub'].map(lambda x: {
        'Lymphocytes (non-tumor),0': 'Lymphocytes1 (non-tumor)',
        'Lymphocytes (non-tumor),1': 'Lymphocytes2 (non-tumor)',
        'Lymphocytes (non-tumor),2': 'Lymphocytes3 (non-tumor)'
    }.get(x, x)).astype("category")
)

sc.pl.umap(adata, color=['cell_type_semifinal_v5'], ncols=1)

#Wilcoxon markers for lymphocyte subtypes:
sc.tl.rank_genes_groups(adata, 'cell_type_semifinal_v5', method='wilcoxon',
                        key_added="Lymph_sub")

result = adata.uns['Lymph_sub']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
     for group in groups
     for key in ['names', 'pvals', 'pvals_adj', 'logfoldchanges', 'scores']}
).head(10).to_csv("Lymphocyte_subclust_180123.csv")

#── 7. Export annotation tables ───────────────────────────────────────────────
adata.obs['cell_type_semifinal_v4'].to_csv(
    './pNET_comparative_clustering/pNET_updated_annotation_v4_180123.csv')
adata.obs['cell_type_semifinal_v5'].to_csv(
    './pNET_comparative_clustering/pNET_updated_annotation_v5_180123.csv')

#── 8. Seurat/R export — minimal obs object ───────────────────────────────────
#Strip all embeddings, graphs, and var metadata; keep only essential obs columns
#for transfer to Seurat (Logistic Regression DEG analysis)
import os
os.makedirs("seurat_analysis", exist_ok=True)

adata_new = adata.copy()
del adata_new.obsm
del adata_new.varm
del adata_new.uns
del adata_new.obsp
del adata_new.obs
del adata_new.var

#Restore essential obs columns only:
adata_new.obs['cluster_label_single_MV']  = adata.obs['cluster_label_single_MV']
adata_new.obs['cell_type_semifinal_v3']   = adata.obs['cell_type_semifinal_v3']
adata_new.obs['cell_type_semifinal_v4']   = adata.obs['cell_type_semifinal_v4']
adata_new.obs['cell_type_semifinal_v5']   = adata.obs['cell_type_semifinal_v5']
adata_new.obs['sample_id_final']          = adata.obs['sample_id_final']
adata_new.obs['source_label']             = adata.obs['source_label']

#Expand to full gene space via raw slot before writing:
adata_new = adata_new.raw.to_adata()
del adata_new.var

results_seurat = "./seurat_analysis/pNET_seurat_12122022.h5ad"
adata_new.write(results_seurat)

#── 9. Save final annotated object ────────────────────────────────────────────
results_final = "./pNET_comparative_clustering/pNET_updated_annotations_180123.h5ad"
adata.write(results_final)

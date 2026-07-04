"""
panNEC Single-Cell Atlas — Chain Integration Script (Step 1)
=============================================================
Description:
    This script performs the core integration and clustering pipeline for the
    panNEC single-nucleus RNA-seq atlas. Starting from per-patient SoupX-corrected
    and normalized h5ad files, it:
      1. Concatenates all 5 patient datasets
      2. Removes MALAT1 and re-normalizes
      3. Computes QC metrics (MT, Ribo, HB genes) and cell cycle scores
      4. Identifies highly variable genes (HVGs) with and without batch correction
      5. Performs PCA and batch correction via Harmony
      6. Computes UMAP using BBKNN (batch-balanced k-nearest neighbours)
      7. Clusters cells using Leiden at multiple resolutions
      8. Runs Wilcoxon rank-sum marker gene detection
      9. Computes per-sample PCA correlation matrix

Author: Olivia Debnath (Current email: olivia_debnath@dfci.harvard.edu)

Associated manuscript:
    Aberrant brain-type neuronal programs in large-cell pancreatic
    neuroendocrine carcinoma
    Zenodo: https://doi.org/10.5281/zenodo.20097469
    GEO: GSE291816
"""

#── Libraries ──────────────────────────────────────────────────────────────────
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2

#── 1. Load per-patient normalized data ───────────────────────────────────────
#SoupX-corrected, integer count h5ad files per patient from Dr. Hilmar Berger 
adata_sp7  = sc.read_h5ad("/dh-projects/ag-ishaque/analysis/debnatho/pNET_scRNAseq_hilmar_2022/normalized_patient_data/NET_Soupx_integer_P7.h5ad")
adata_sp24 = sc.read_h5ad("/dh-projects/ag-ishaque/analysis/debnatho/pNET_scRNAseq_hilmar_2022/normalized_patient_data/NET_Soupx_integer_P24ab.h5ad")
adata_sp18 = sc.read_h5ad("/dh-projects/ag-ishaque/analysis/debnatho/pNET_scRNAseq_hilmar_2022/normalized_patient_data/NET_Soupx_integer_P18.h5ad")
adata_sp19 = sc.read_h5ad("/dh-projects/ag-ishaque/analysis/debnatho/pNET_scRNAseq_hilmar_2022/normalized_patient_data/NET_Soupx_integer_P19.h5ad")
adata_sp21 = sc.read_h5ad("/dh-projects/ag-ishaque/analysis/debnatho/pNET_scRNAseq_hilmar_2022/normalized_patient_data/NET_Soupx_integer_P21.h5ad")

#── 2. Concatenate all patients ────────────────────────────────────────────────
adata_net = adata_sp7.concatenate(adata_sp24, adata_sp21, adata_sp18, adata_sp19, join="outer")

#Merge both technical replicates of P24 under the same source label
adata_net.obs['source_label'] = adata_net.obs["sample_id_final"].map(
    lambda x: {"SP084_027": "SP084_024",
               "SP084_028": "SP084_024"}.get(x, x)
).astype("category")

adata_net.obs['source_label'].value_counts()

#── 3. Remove MALAT1 ───────────────────────────────────────────────────────────
#MALAT1 is a ubiquitous nuclear lncRNA that dominates snRNA-seq libraries
#and can obscure biologically meaningful variation; it is removed prior to
#normalization
malat1 = adata_net.var_names.str.startswith('MALAT1')
remove  = malat1
keep    = np.invert(remove)

adata_new = adata_net[:, keep]
adata_new.n_vars

adata_new.X.expm1().sum(axis=1)  # confirm not yet normalized after MALAT1 removal

#── 4. Normalize ──────────────────────────────────────────────────────────────
#Re-normalize after MALAT1 removal to ensure library-size correction is applied
#to the remaining gene set
sc.pp.normalize_total(adata_new, target_sum=1e4)
sc.pp.log1p(adata_new)

#Confirm normalization
adata_new.X.expm1().sum(axis=1)

#Preserve the normalized full data slot for downstream use (e.g. visualisation,
#differential expression on all genes)
adata_new.raw = adata_new

#── 5. QC metrics ─────────────────────────────────────────────────────────────
adata_new.var['MT_genes']   = adata_new.var_names.str.startswith('MT-')
adata_new.var['Ribo_genes'] = adata_new.var_names.str.startswith(("RPS", "RPL"))
adata_new.var['HB_genes']   = adata_new.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(adata_new, qc_vars=['MT_genes', 'Ribo_genes', 'HB_genes'],
                            percent_top=None, inplace=True)

#Mitochondrial fraction per cell
mito_genes = adata_new.var_names.str.startswith('MT-')
adata_new.obs['percent_mt2'] = (
    np.sum(adata_new[:, mito_genes].X, axis=1).A1 /
    np.sum(adata_new.X, axis=1).A1
)

#Total counts per cell
adata_new.obs['n_counts'] = adata_new.X.sum(axis=1).A1

#Ribosomal fraction per cell
ribo_genes = adata_new.var_names.str.startswith(("RPS", "RPL"))
adata_new.obs['percent_Ribo2'] = (
    np.sum(adata_new[:, ribo_genes].X, axis=1).A1 /
    np.sum(adata_new.X, axis=1).A1
)

#Haemoglobin fraction per cell
hb_genes = adata_new.var_names.str.contains(("^HB[^(P)]"))
adata_new.obs['percent_HB2'] = (
    np.sum(adata_new[:, hb_genes].X, axis=1).A1 /
    np.sum(adata_new.X, axis=1).A1
)

#── 6. Cell cycle scoring ─────────────────────────────────────────────────────
cell_cycle_genes = [x.strip() for x in open(
    '/dh-projects/preeclampsia_2022/analysis/herse4olivia/regev_lab_cell_cycle_genes.txt')]
print(len(cell_cycle_genes))

s_genes   = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata_new.var_names]
print(len(cell_cycle_genes))

sc.tl.score_genes_cell_cycle(adata_new, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pl.violin(adata_new, ['S_score', 'G2M_score'],
             jitter=0, groupby='sample_id_final', rotation=90)

#Add XIST expression to obs for sex-check
adata_new.obs["XIST-counts"] = adata_new.X[:, adata_new.var_names.str.match('XIST')].toarray()

sc.pl.violin(adata_new, ["XIST-counts"], jitter=0, groupby='source_label', rotation=90)

#── 7. Highly variable gene (HVG) selection ───────────────────────────────────
#Two strategies were evaluated:
#
#Case I — No batch correction during HVG selection:
#  Selects genes variable across the full dataset without accounting for
#  patient-specific effects. Risk: retains batch-specific genes (e.g. YAP1
#  for P7, acinar genes for P24).
#
#Case II — Light-weight batch correction (used in manuscript):
#  Genes are scored per sample (batch_key="source_label") and only genes
#  variable in a minimum of 2 samples are retained. This prevents the
#  selection of highly patient-specific genes while preserving biologically
#  meaningful shared variation. Note: major acinar genes are lost with this
#  approach, which was acceptable given our focus on NE cell states.

adata_hvg = adata_new.copy()

#Case I: No batch key
sc.pp.highly_variable_genes(adata_hvg, n_top_genes=4000, subset=True)
hvg_df = pd.DataFrame(adata_hvg.var_names)

#Case II: Light-weight batch correction (source_label has P24 TRs merged)
adata_hvg02 = adata_new.copy()
sc.pp.highly_variable_genes(adata_hvg02, n_top_genes=4000,
                             batch_key="source_label", subset=True)
hvg_df02 = pd.DataFrame(adata_hvg02.var_names)

#Compare HVG overlap between strategies
venn2([set(hvg_df[0]), set(hvg_df02[0])], ('HVG_no_BC', 'HVG_lightwt_BC'))

#── 8. PCA ────────────────────────────────────────────────────────────────────
sc.tl.pca(adata_hvg02)
sc.pl.pca_overview(adata_hvg02, color="source_label")

#Inspect top genes driving each PC
adata_hvg02.var['gene_ids'] = adata_hvg02.var.index
genes = adata_hvg02.var['gene_ids']

for pc in [1, 2, 3, 4, 5, 6, 7, 8]:
    g   = adata_hvg02.varm['PCs'][:, pc - 1]
    o   = np.argsort(g)
    sel = np.concatenate((o[:10], o[-10:])).tolist()
    emb = adata_hvg02.obsm['X_pca'][:, pc - 1]
    tempdata = adata_hvg02[np.argsort(emb),]
    sc.pl.heatmap(tempdata, var_names=genes[sel].index.tolist(),
                  groupby='source_label', swap_axes=True, use_raw=False)

#── 9. Initial UMAP (pre-integration) ─────────────────────────────────────────
sc.pp.neighbors(adata_hvg02)
sc.tl.umap(adata_hvg02, random_state=0, maxiter=200)
sc.pl.umap(adata_hvg02, color='source_label')  # samples cleanly separate pre-integration

#── 10. Harmony batch correction ──────────────────────────────────────────────
#Harmony corrects the PCA embedding by iteratively adjusting cluster centroids
#to be sample-invariant. It operates in PCA space (not gene space) and is
#therefore computationally efficient. Correction is performed per
#sample_id_final to account for individual patient effects.
#The Harmony-corrected embedding (X_pca_harmony) is used as input for BBKNN.
sce.pp.harmony_integrate(adata_hvg02, 'sample_id_final')

adata_corr = adata_hvg02.copy()
adata_corr.obsm['X_pca'] = adata_corr.obsm['X_pca_harmony']

#Inspect PCs after Harmony correction
sc.pl.pca_overview(adata_corr, color="source_label",
                   components=['1,2', '2,3', '1,3', '1,4'])

#── 11. BBKNN (used in manuscript) ───────────────────────────────────────────
#Batch-balanced k-nearest neighbours (BBKNN) builds a neighbourhood graph
#that is balanced across batches (patients). For each cell, it finds
#`neighbors_within_batch` nearest neighbours in each batch separately,
#then merges these into a single graph. This prevents dominant batches from
#monopolising the neighbourhood of any given cell.
#Parameters:
#  batch_key='sample_id_final'  — correct per individual patient
#  neighbors_within_batch=5     — 5 neighbours per patient per cell
#  use_rep='X_pca_harmony'      — use Harmony-corrected PCA as input
#  trim=80                      — trim edges beyond the 80th percentile
#                                 connectivity to reduce noise
sce.pp.bbknn(adata_hvg02, batch_key='sample_id_final',
             neighbors_within_batch=5, use_rep='X_pca_harmony', trim=80)

sc.tl.umap(adata_hvg02, maxiter=200, random_state=0)
sc.pl.umap(adata_hvg02, color='source_label')

#── 12. Leiden clustering ─────────────────────────────────────────────────────
sc.tl.leiden(adata_hvg02, resolution=1,   key_added="res_1")
sc.tl.leiden(adata_hvg02, resolution=0.5, key_added="res_0.5")
sc.pl.umap(adata_hvg02, color='res_0.5')

#── 13. Marker gene detection (Wilcoxon rank-sum) ────────────────────────────
#Initial marker gene detection per Leiden cluster using Wilcoxon rank-sum test.
#Note: final cell type annotation in the manuscript used Logistic Regression
#in Seurat for robustness.
sc.tl.rank_genes_groups(adata_hvg02, 'res_0.5', method='wilcoxon')
sc.tl.dendrogram(adata_hvg02, groupby="res_0.5")

sc.pl.rank_genes_groups_dotplot(
    adata_hvg02,
    n_genes=10,
    values_to_plot="logfoldchanges",
    cmap='bwr',
    vmin=-4,
    vmax=4,
    min_logfoldchange=3,
    colorbar_title='log fold change'
)

#Compare proliferative cluster markers across resolutions
sc.pl.dotplot(adata_hvg02, ['FANCA', 'FANCI', 'EZH2', 'BRIP1'], groupby="res_1")
#Leiden cluster 5 is proliferative

#── 14. Save integrated object ────────────────────────────────────────────────
adata_hvg02.write("./Harmony_bbknn_pNET/pNET_harmonization_4000hvg_BC_181022.h5ad")

#── 15. Per-sample PCA correlation matrix ────────────────────────────────────
#Computes the mean PCA coordinate per sample (in uncorrected PCA space) and
#calculates Pearson correlation between samples. This serves as a quality
#check for batch effects and sample similarity prior to integration.
annot_key = 'source_label'

for group in adata_hvg02.obs[annot_key].cat.categories:
    print(group)
    adata_hvg02_sub = adata_hvg02.copy()
    ix = np.in1d(adata_hvg02.obs[annot_key], group)
    adata_hvg02_sub = adata_hvg02.copy()[ix]
    adata_hvg02.uns['mean_pca_' + group] = np.squeeze(
        np.asarray(adata_hvg02_sub.obsm['X_pca'].mean(axis=0).T)
    )
    del adata_hvg02_sub

#Compile mean PCA vectors across samples
d = {
    'SP084_007': adata_hvg02.uns['mean_pca_SP084_007'],
    'SP084_018': adata_hvg02.uns['mean_pca_SP084_018'],
    'SP084_019': adata_hvg02.uns['mean_pca_SP084_019'],
    'SP084_021': adata_hvg02.uns['mean_pca_SP084_021'],
    'SP084_024': adata_hvg02.uns['mean_pca_SP084_024']
}

df = pd.DataFrame(data=d)

#Pearson correlation between samples based on mean PCA coordinates
df_cor_pearson       = np.corrcoef(df.T)
names                = df.columns
df_cor_named_pearson = pd.DataFrame(df_cor_pearson, index=names, columns=names)

#Hierarchically-clustered heatmap of inter-sample PCA correlations
#HVGs where genes are variable in at least 2 or more batches
ax = sns.clustermap(df_cor_named_pearson, cmap='RdBu_r',
                    method='average', figsize=(5, 5))

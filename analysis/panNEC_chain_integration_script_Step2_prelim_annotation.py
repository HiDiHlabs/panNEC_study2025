"""
panNEC Single-Cell Atlas — Chain Integration Script (Step 2: Annotation)
========================================================================
Description:
    This script performs cell type annotation of the panNEC atlas following
    Harmony/BBKNN integration (Step 1). Starting from the integrated h5ad,
    it:
      1. Explores clustering at multiple Leiden resolutions (0.5, 0.8, 1.0)
      2. Identifies cluster markers using Wilcoxon rank-sum test
      3. Assigns cell type labels (cell_type_semifinal01)
      4. Scores PDAC Regev lineage signatures against cell states
      5. Saves the annotated object

Reference notebook:
    pNET_scRNAseq_hilmar_2022/notebooks_analysis/pNET_integration_caseII_181022.ipynb

Author: Olivia Debnath (Current email: olivia_debnath@dfci.harvard.edu)

Associated manuscript:
    Aberrant brain-type neuronal programs in large-cell pancreatic
    neuroendocrine carcinoma
    Zenodo: https://doi.org/10.5281/zenodo.20097469
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

sc.settings.set_figure_params(dpi=90, dpi_save=150, facecolor='white')

#── 1. Load integrated object (Step 1 output) ─────────────────────────────────
adata = sc.read_h5ad("./Harmony_bbknn_pNET/pNET_harmonization_4000hvg_BC_181022.h5ad")

#── 2. Inspect PCA & UMAP ─────────────────────────────────────────────────────
#PCA overview (uncorrected):
sc.pl.pca_overview(adata, color="source_label")

#Inspect Harmony-corrected PCs:
adata_corr = adata.copy()
adata_corr.obsm['X_pca'] = adata_corr.obsm['X_pca_harmony']
sc.pl.pca(adata_corr, color="source_label", components=['1,2', '2,3', '1,3', '1,4'], ncols=1)

#UMAP at res_0.5 (computed in Step 1):
sc.pl.umap(adata, color='res_0.5')

#Proliferation markers to identify cycling cluster:
sc.pl.umap(adata, color=['TOP2A', 'MKI67', 'EZH2'])

#── 3. Per-patient UMAP inspection ────────────────────────────────────────────
#Check per source_label to assess batch mixing:
for i in adata.obs['source_label'].cat.categories:
    print(i)
    fig = sc.pl.umap(adata[adata.obs['source_label'] == i],
                     color='res_0.5', return_fig=True, title=i)

#── 4. Leiden clustering at multiple resolutions ──────────────────────────────
#Resolution 0.8:
sc.tl.leiden(adata, resolution=0.8, key_added="res_0.8")
sc.pl.umap(adata, color='res_0.8')
sc.pl.dotplot(adata, ['TOP2A', 'MKI67', 'DIAPH3', 'STMN1', 'EZH2'], groupby="res_0.8")

#Resolution 1.0 (used for final annotation):
sc.tl.leiden(adata, resolution=1, key_added="res_1")
sc.pl.umap(adata, color='res_1')
sc.pl.dotplot(adata, ['TOP2A', 'MKI67', 'DIAPH3', 'STMN1', 'EZH2'], groupby="res_1")

#── 5. Marker gene exploration ────────────────────────────────────────────────
#Stromal/mesenchymal markers:
sc.pl.dotplot(adata, ['COL1A2', 'YAP1', 'WWTR1', 'FAP', 'S100A4'], groupby="res_1")

#EEC & proliferating EEC markers (L4=EEC_prol, L0=EEC):
sc.pl.dotplot(adata, ['ADARB2', 'GIPR', 'COL26A1', 'CENPP', 'FANCA', 'MKI67', 'EZH2'], groupby="res_1")

#DNA repair / acinar-like markers:
sc.pl.dotplot(adata, ['BRCA1', 'FANCA', 'FANCI'], groupby="res_1")
sc.pl.dotplot(adata, ['BRCA1', 'FANCA', 'FANCI', 'EZH2'], groupby="source_label")

#── 6. Wilcoxon rank-sum marker detection (res_1) ────────────────────────────
sc.tl.rank_genes_groups(adata, 'res_1', method='wilcoxon')
sc.tl.dendrogram(adata, groupby="res_1")

sc.pl.rank_genes_groups_dotplot(adata, n_genes=6,
    values_to_plot="logfoldchanges", cmap='bwr',
    vmin=-4, vmax=4, min_logfoldchange=3,
    colorbar_title='LFC cut-off: 3')

sc.pl.rank_genes_groups_dotplot(adata, n_genes=6,
    values_to_plot="logfoldchanges", cmap='bwr',
    vmin=-4, vmax=4, min_logfoldchange=3.5,
    colorbar_title='LFC cut-off: 3.5')

#── 7. Cell type annotation ───────────────────────────────────────────────────
#L2 & L5 are very similar; L2 has no ductal markers like PROM1/ABI3BP
adata.obs['cell_type_semifinal01'] = (
    adata.obs['res_1'].map(lambda x: {
        '7':  'Lymphocytes',
        '10': 'Microglia_like',
        '1':  'Mesenchymal',
        '0':  'EEC',
        '4':  'EEC_prol',
        '11': 'BRCA1_FANCA_acinar_like',
        '3':  'GP2_acinar_like',
        '8':  'GP2_acinar_like',
        '5':  'Exo_glandular_mixed',
        '6':  'HSP_pos',
        '9':  'SMC_like',
        '2':  'Exo_glandular_acinar_like',
        '12': 'GP2_acinar_like'
    }.get(x, x)).astype("category")
)

adata.obs['cell_type_semifinal01'].value_counts()

#SMC-like vs acinar marker validation:
sc.pl.dotplot(adata, ['ENPP3', 'FREM2', 'NKD1', 'GP2', 'RBPJL',
                      'COL1A2', 'ABI3BP', 'PROM1'], groupby="res_1")

#── 8. Re-stratify markers by cell_type_semifinal01 ──────────────────────────
sc.tl.rank_genes_groups(adata, 'cell_type_semifinal01', method='wilcoxon')
sc.tl.dendrogram(adata, groupby="cell_type_semifinal01")

sc.pl.rank_genes_groups_dotplot(adata, n_genes=10,
    values_to_plot="logfoldchanges", cmap='bwr',
    vmin=-4, vmax=4, min_logfoldchange=3,
    colorbar_title='LFC cut-off: 3')

sc.pl.umap(adata, color='cell_type_semifinal01')

#Full marker panel dotplot:
sc.pl.dotplot(adata, [
    'APCDD1', 'NKD1', 'ENPP3', 'FREM2', 'SEZ6L', 'NOTUM', 'RASGRF1', 'ALK', 'DYNC1I1',
    'MKI67', 'DIAPH3', 'TOP2A', 'CENPP', 'CENPK', 'NUSAP1',
    'BRCA1', 'BRIP1', 'DTL', 'ATAD5', 'FANCI', 'FANCA',
    'RIMBP2', 'ADARB2', 'GSE1', 'LDLRAD4', 'GIPR', 'COL1A2', 'CALD1', 'CASC15', 'BICC1', 'CDH11',
    'HSPH1', 'HSPE1', 'HSPB1', 'PGK1', 'NDRG1', 'VEGFA', 'DNAJB1', 'RBPJL', 'GP2', 'NR5A2', 'MECOM', 'CEL',
    'PROM1', 'ABI3BP', 'FKBP9', 'IGFN1',
    'PTPRC', 'ARHGAP15', 'CELF2', 'ETS1', 'IKZF1', 'DOCK2', 'DOCK8', 'PLXDC2', 'ZEB2', 'DPYD', 'SAT1', 'CTSB'
], groupby="cell_type_semifinal01", dendrogram=True)

#── 9. Reorder categories & patient composition ───────────────────────────────
adata.obs['cell_type_semifinal_v2'] = adata.obs['cell_type_semifinal01'].cat.reorder_categories([
    'BRCA1_FANCA_acinar_like', 'GP2_acinar_like', 'Exo_glandular_acinar_like', 'Exo_glandular_mixed',
    'SMC_like', 'HSP_pos', 'EEC', 'EEC_prol', 'Lymphocytes', 'Mesenchymal', 'Microglia_like'
])

df = pd.crosstab(adata.obs['cell_type_semifinal_v2'], adata.obs['sample_id_final'], normalize="index")

plt.rcParams["figure.figsize"] = (8, 4)
pl = df.plot(kind="bar", stacked=True, rot=90)
pl.grid(False)
pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#── 10. Export marker table ───────────────────────────────────────────────────
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
     for group in groups
     for key in ['names', 'pvals', 'pvals_adj', 'logfoldchanges', 'scores']}
).head(30).to_csv("./CUP7_NEC_manuscript2022/Markers_Wilcoxon_241022.csv")

#── 11. Save annotated object ─────────────────────────────────────────────────
results = "./pNET_comparative_clustering/pNET_harmony_bbknn_241022.h5ad"
adata.write(results)

#── 12. PDAC Regev signature scoring ─────────────────────────────────────────
#Score malignant cell state & lineage programme signatures from Broad PDAC atlas
#(Raghavan et al.) against panNEC cell types
pdac_regev = pd.read_csv("./pNET_comparative_clustering/PDAC_Regev_signatures.csv")

#Score all signatures (top 30 genes per signature):
for i in pdac_regev.columns:
    pdac_sig = pdac_regev[i].dropna().values.tolist()[:30]
    sc.tl.score_genes(adata, pdac_sig, score_name=[i])

#Visualize scored signatures by cell type:
sig_cols = ['Cycling (S)', 'Cycling (G2/M)', 'MYC signaling', 'Adhesive', 'Ribosomal',
            'Interferon signaling', 'TNF-NFkB signaling', 'Acinar-like', 'Classical-like',
            'Basaloid', 'Squamoid', 'Mesenchymal', 'Neuroendocrine-like', 'Neural-like progenitor',
            'iPSC', 'NPC', 'fetal replicating', 'fetal quiescent', 'Oligodendrocyte precursors',
            'Neurons', 'Astrocytes', 'Oligodendrocytes', 'Microglia',
            'Highly-expressed in RB1-loss:', 'Low-expressed in RB1-loss:']

sc.pl.dotplot(adata, sig_cols, groupby="cell_type_semifinal_v2", dendrogram=True, vmax=0.4)

#Tumor cells only (exclude non-malignant populations):
non_tumor = ['3_Lymphocytes', '3_Macrophages', '3_Non-Tm S21', '4_Leukocytes NOS',
             '5_Endothelial', '5_Macrophages', '6_Leukocytes NOS', '2_Lymphocytes', '2_Non-Tm S21']

sc.pl.dotplot(adata[~adata.obs['cluster_label_single_MV'].isin(non_tumor)],
              sig_cols[:15], groupby="cell_type_semifinal_v2")

"""
panNEC Single-Cell Atlas — Chain Integration Script (Step 3: Interpatient Variability)
=======================================================================================
Description:
    This script performs refined annotation and interpatient variability analysis
    of the panNEC atlas following Step 2. Starting from the annotated h5ad, it:
      1. Refines cell type labels (cell_type_semifinal02)
      2. Identifies and filters L12 barcode contamination (P18/P19/P21)
      3. Recomputes UMAP post-filtering
      4. Scores HSP/hypoxia/AGE-RAGE/O-linked glycosylation gene signatures
      5. Runs Wilcoxon marker detection on refined labels
      6. Computes adjusted mutual information & silhouette batch scores
      7. Saves the filtered annotated object

Reference notebook:
    pNET_scRNAseq_hilmar_2022/notebooks_analysis/pNET_interpatient_variability_update_02112022.ipynb

Note: Integration variant checks (PC=10/15/20/30/40 BBKNN sweeps) are documented
    separately in panNEC_PCA_integration_variants_check_24102022.ipynb.
    The final atlas uses PC=50 (default).

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
import scanpy.external as sce
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.metrics.cluster import silhouette_samples, silhouette_score

sc.settings.set_figure_params(dpi=90, dpi_save=150, facecolor='white')

#── 1. Load annotated object (Step 2 output) ──────────────────────────────────
results = "./pNET_comparative_clustering/pNET_harmony_bbknn_241022.h5ad"
adata_annot = sc.read_h5ad(results)

#── 2. Inspect existing annotations ───────────────────────────────────────────
sc.pl.umap(adata_annot, color='cell_type_semifinal01')
sc.pl.umap(adata_annot, color='source_label')
sc.pl.umap(adata_annot, color='res_1')

#── 3. EEC marker validation ──────────────────────────────────────────────────
#Lacks major hormones (INS, GCG, SST); ~30% of EEC state expresses GIPR &
#voltage-gated gene sets (known EEC enrichment)
sc.pl.dotplot(adata_annot, [
    'INS', 'SST', 'GHRL', 'GCG', 'GIPR', 'CLDN4', 'CHGA', 'FFAR4', 'PCSK1', 'ARX',
    'CHGB', 'GFRA3', 'TRPA1', 'CPE', 'GAL', 'ENPP2', 'VEGFA', 'TPH1', 'ALCAM',
    'CACNA1A', 'CACNB2', 'CACNA1C', 'PDE4D', 'PLCB1', 'VAV2', 'SIPA1L3',
    'NRXN1', 'KCNJ3', 'KCNJ6', 'MKI67', 'DIAPH3'
], groupby="res_1")

#── 4. BRCA1/FANCA proliferative cluster (L4/L11) ────────────────────────────
#L11 separates from P24 owing to robust acinar marker expression
sc.pl.dotplot(adata_annot, [
    'BRCA1', 'BRCA2', 'FANCA', 'FANCI', 'DIAPH3', 'MKI67', 'CENPE',
    'GP2', 'MECOM', 'NR5A2', 'RBPJL'
], groupby="res_1")

#── 5. Acinar & stromal marker exploration ────────────────────────────────────
#Mesenchymal (L1): mainly P7 contributing; separate P7 analysis done elsewhere
sc.pl.dotplot(adata_annot, ['COL1A1', 'COL1A2', 'CALD1', 'CASC15', 'CDH11'], groupby="res_1")

#Leiden composition by source:
pd.crosstab(adata_annot.obs['source_label'], adata_annot.obs['res_1'])

#Acinar/stromal mixed populations — PROM1=ductal type; ABI3BP=stromal/fibroblasts
sc.pl.dotplot(adata_annot, [
    'TPH1', 'TAC1', 'GP2', 'NR5A2', 'MECOM', 'RBPJL', 'CEL',
    'PRSS1', 'PRSS2', 'AMY2A', 'ELANE', 'REG3A', 'REG3G', 'REG1A',
    'PROM1', 'ABI3BP', 'FKBP9', 'IGFN1', 'BICC1', 'ANXA4', 'SPP1', 'CFTR'
], groupby="res_1", dendrogram=True)

#── 6. Refined cell type annotation (cell_type_semifinal02) ───────────────────
#L2 & L5 are very similar; L2 has low ductal markers (PROM1/ABI3BP)
#but medium FKBP9/IGFN1 (ECM remodelling)
adata_annot.obs['cell_type_semifinal02'] = (
    adata_annot.obs['res_1'].map(lambda x: {
        '7':  'Lymphocytes',
        '10': 'Microglia_like/Macrophages',
        '1':  'Stromal',
        '0':  'GIPR_EEC_like',
        '4':  'BRCA1/FANCA_prol01',
        '11': 'BRCA1/FANCA_prol02_acinar_like',
        '3':  'AcinarS_NET_like',
        '8':  'AcinarS_like',
        '5':  'Exo_glandular_mixed01',
        '6':  'HSP_pos',
        '9':  'SMC_like',
        '2':  'Exo_glandular_mixed02',
        '12': 'AcinarS_like'
    }.get(x, x)).astype("category")
)

adata_annot.obs['cell_type_semifinal02'].value_counts()
sc.pl.umap(adata_annot, color='cell_type_semifinal02')

#── 7. Full marker panel dotplot ──────────────────────────────────────────────
sc.pl.dotplot(adata_annot, [
    'APCDD1', 'NKD1', 'ENPP3', 'FREM2', 'SEZ6L', 'NOTUM', 'RASGRF1', 'ALK', 'DYNC1I1',
    'RBPJL', 'GP2', 'NR5A2', 'MECOM', 'CEL', 'PRSS1', 'PRSS2',
    'MKI67', 'DIAPH3', 'TOP2A', 'CENPP', 'CENPK', 'NUSAP1', 'BRCA1', 'BRIP1',
    'DTL', 'ATAD5', 'FANCI', 'FANCA',
    'TPH1', 'PROM1', 'ABI3BP', 'FKBP9', 'IGFN1', 'NRXN3', 'DKK4',
    'CACNA1A', 'CERS4', 'TNS1', 'NEDD4L', 'PLXNA2', 'ACACB', 'RASA3',
    'RIMBP2', 'ADARB2', 'GSE1', 'LDLRAD4', 'CACNB2', 'KCNJ3', 'GIPR',
    'COL1A1', 'COL1A2', 'CALD1', 'CASC15', 'BICC1', 'CDH11',
    'HSPH1', 'HSPE1', 'HSPB1', 'PGK1', 'NDRG1', 'VEGFA', 'DNAJB1',
    'PTPRC', 'ARHGAP15', 'CELF2', 'ETS1', 'IKZF1', 'DOCK2', 'DOCK8',
    'PLXDC2', 'ZEB2', 'DPYD', 'SAT1', 'CTSB'
], groupby="cell_type_semifinal02", dendrogram=True, standard_scale="var")

#── 8. L12 barcode contamination filter ───────────────────────────────────────
#L12 contains P18/P19/P21 contamination — filter these barcodes
sc.tl.rank_genes_groups(adata_annot, groupby="res_1", groups=['12'],
                        key_added="L12 vs rest", method='wilcoxon')
sc.pl.rank_genes_groups_dotplot(adata_annot, key="L12 vs rest", n_genes=20)

filter_barcodes = adata_annot.obs.loc[
    (adata_annot.obs.res_1 == '12') &
    (adata_annot.obs.source_label.isin(["SP084_021", "SP084_019", "SP084_018"]))
].index.tolist()

adata_annot_filter = adata_annot[~adata_annot.obs_names.isin(filter_barcodes)]
print(f"Cells after filtering: {adata_annot_filter.n_obs}")

#Recompute UMAP post-filtering
sc.tl.umap(adata_annot_filter, maxiter=200, random_state=0)
sc.pl.umap(adata_annot_filter, color='cell_type_semifinal02')

#── 9. Reorder categories & patient composition ───────────────────────────────
adata_annot_filter.obs['cell_type_semifinal_v2'] = (
    adata_annot_filter.obs['cell_type_semifinal02'].cat.reorder_categories([
        'SMC_like', 'BRCA1/FANCA_prol02_acinar_like', 'Exo_glandular_mixed01',
        'Exo_glandular_mixed02', 'AcinarS_like', 'AcinarS_NET_like',
        'HSP_pos', 'BRCA1/FANCA_prol01', 'GIPR_EEC_like', 'Stromal',
        'Lymphocytes', 'Microglia_like/Macrophages'
    ])
)

df = pd.crosstab(adata_annot_filter.obs['cell_type_semifinal_v2'],
                 adata_annot_filter.obs['sample_id_final'], normalize="index")

plt.rcParams["figure.figsize"] = (6, 4)
pl = df.plot(kind="bar", stacked=True, rot=90)
pl.grid(False)
pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#── 10. Wilcoxon marker detection ─────────────────────────────────────────────
sc.tl.rank_genes_groups(adata_annot_filter, 'cell_type_semifinal02',
                        method='wilcoxon', pts=True)
sc.tl.dendrogram(adata_annot_filter, groupby="cell_type_semifinal02")

sc.pl.rank_genes_groups_dotplot(adata_annot_filter, n_genes=10,
    values_to_plot="logfoldchanges", cmap='bwr',
    vmin=-4, vmax=4, min_logfoldchange=3,
    colorbar_title='LFC cut-off: 3',
    save="_pNET_states_logFC_Wilcoxon_021122.pdf")

#Export marker table:
result = adata_annot_filter.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
     for group in groups
     for key in ['names', 'pvals', 'pvals_adj', 'logfoldchanges', 'scores']}
).head(100).to_csv("./pNET_comparative_clustering/pNET_Wilcoxon_marker_table_02112022.csv")

#── 11. Gene signature scoring ────────────────────────────────────────────────
#HSP stress response:
HSP_resp = ['HSPE1', 'HSPH1', 'HSPA9', 'HSPA1A', 'HSPA1B', 'HSPD1', 'DNAJB1']

#MSigDB Hallmark Hypoxia:
msigdb_hypoxia = [
    "ACKR3","ADM","ADORA2B","AK4","AKAP12","ALDOA","ALDOB","ALDOC","AMPD3","ANGPTL4",
    "ANKZF1","ANXA2","ATF3","ATP7A","B3GALT6","B4GALNT2","BCAN","BCL2","BGN","BHLHE40",
    "BNIP3L","BRS3","BTG1","CA12","CASP6","CAV1","CAVIN1","CAVIN3","CCN1","CCN2","CCN5",
    "CCNG2","CDKN1A","CDKN1B","CDKN1C","CHST2","CHST3","CITED2","COL5A1","CP","CSRP2",
    "CXCR4","DCN","DDIT3","DDIT4","DPYSL4","DTNA","DUSP1","EDN2","EFNA1","EFNA3","EGFR",
    "ENO1","ENO2","ENO3","ERO1A","ERRFI1","ETS1","EXT1","F3","FAM162A","FBP1","FOS",
    "FOSL2","FOXO3","GAA","GALK1","GAPDH","GAPDHS","GBE1","GCK","GCNT2","GLRX","GPC1",
    "GPC3","GPC4","GPI","GRHPR","GYS1","HAS1","HDLBP","HEXA","HK1","HK2","HMOX1",
    "HOXB9","HS3ST1","HSPA5","IDS","IER3","IGFBP1","IGFBP3","IL6","ILVBL","INHA","IRS2",
    "ISG20","JMJD6","JUN","KDELR3","KDM3A","KIF5A","KLF6","KLF7","KLHL24","LALBA",
    "LARGE1","LDHA","LDHC","LOX","LXN","MAFF","MAP3K1","MIF","MT1E","MT2A","MXI1",
    "MYH9","NAGK","NCAN","NDRG1","NDST1","NDST2","NEDD4L","NFIL3","NOCT","NR3C1",
    "P4HA1","P4HA2","PAM","PCK1","PDGFB","PDK1","PDK3","PFKFB3","PFKL","PFKP","PGAM2",
    "PGF","PGK1","PGM1","PGM2","PHKG1","PIM1","PKLR","PKP1","PLAC8","PLAUR","PLIN2",
    "PNRC1","PPARGC1A","PPFIA4","PPP1R15A","PPP1R3C","PRDX5","PRKCA","PYGM","RBPJ",
    "RORA","RRAGD","S100A4","SAP30","SCARB1","SDC2","SDC3","SDC4","SELENBP1","SERPINE1",
    "SIAH2","SLC25A1","SLC2A1","SLC2A3","SLC2A5","SLC37A4","SLC6A6","SRPX","STBD1",
    "STC1","STC2","SULT2B1","TES","TGFB3","TGFBI","TGM2","TIPARP","TKTL1","TMEM45A",
    "TNFAIP3","TPBG","TPD52","TPI1","TPST2","UGP2","VEGFA","VHL","VLDLR","WSB1",
    "XPNPEP1","ZFP36","ZNF292"
]

age_rage = ['BCL2', 'SMAD3', 'TGFBR1', 'VEGFA']
Olinked_glyco = ['LARGE1', 'B4GALT6', 'MUC13', 'ADAMTS1', 'THSD4', 'B3GNT5', 'GALNT13', 'GALNT18']

sc.tl.score_genes(adata_annot_filter, HSP_resp, score_name="HSP_response")
sc.tl.score_genes(adata_annot_filter, msigdb_hypoxia, score_name="MSigDB_Hallmark_Hypoxia")
sc.tl.score_genes(adata_annot_filter, age_rage, score_name="AGE_RAGE_signaling")
sc.tl.score_genes(adata_annot_filter, Olinked_glyco, score_name="Olinked_glycosylation")

#Tumor cells only (exclude non-malignant):
tumor_only = ~adata_annot_filter.obs['cell_type_semifinal02'].isin(
    ['Lymphocytes', 'Microglia_like/Macrophages'])

sc.pl.dotplot(adata_annot_filter[tumor_only],
    ['HSP_response', 'MSigDB_Hallmark_Hypoxia', 'AGE_RAGE_signaling', 'Olinked_glycosylation'],
    groupby="cell_type_semifinal02", dendrogram=False, standard_scale="var",
    save="_Hypoxia_pathway_plot_021122.pdf")

#Hypoxia-HSP gene-level dotplot:
sc.pl.dotplot(adata_annot_filter[tumor_only], [
    'IARS2', 'HSPE1', 'HSPH1', 'HSPA9', 'HSPA1A', 'HSPA1B', 'HSPD1', 'DNAJB1',
    'VEGFA', 'HIF1A', 'NDRG1', 'TPD52', 'MIF', 'PGAM1', 'LDHA', 'P4HA1', 'SLC2A1'
], groupby="cell_type_semifinal02", dendrogram=True, standard_scale="var",
    save="_HSPpos_hypoxia_tumor_021122.pdf")

#── 12. Batch integration quality metrics ─────────────────────────────────────
#Adjusted mutual information — measures dependency between cell type & patient labels
#(lower = better mixing across patients within cell types)
ami_score = adjusted_mutual_info_score(
    adata_annot_filter.obs['cell_type_semifinal_v2'],
    adata_annot_filter.obs['source_label']
)
print(f"Adjusted Mutual Information (cell type vs patient): {ami_score:.4f}")

#Silhouette batch score (adapted from Luecken et al. 2021, Nature Methods)
#Measures how well batch effects are removed within each cell type
#Score range: 0 (worst) to 1 (best mixing across batches)
def silhouette_batch(adata, batch_key, group_key, embed,
                     metric='euclidean', scale=True, verbose=True):
    """
    Absolute silhouette score of batch labels subsetted per cell type group.
    Adapted from Luecken et al. 2021, Nature Methods.
    Higher score = better batch mixing within cell types.
    """
    if embed not in adata.obsm.keys():
        raise KeyError(f'{embed} not in obsm')

    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])

    for group in adata.obs[group_key].unique():
        adata_group = adata[adata.obs[group_key] == group]
        n_batches = adata_group.obs[batch_key].nunique()
        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue
        sil_per_group = silhouette_samples(
            adata_group.obsm[embed], adata_group.obs[batch_key], metric=metric)
        sil_per_group = [abs(i) for i in sil_per_group]
        if scale:
            sil_per_group = [1 - i for i in sil_per_group]
        sil_all = sil_all.append(pd.DataFrame({
            'group': [group] * len(sil_per_group),
            'silhouette_score': sil_per_group
        }))

    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby('group').mean()
    asw = sil_means['silhouette_score'].mean()
    if verbose:
        print(f'Mean silhouette per cell type:\n{sil_means}')
    return asw

emb_key_ = "X_umap"

#PC=50 (final atlas):
print("ASW (PC=50, all groups):")
silhouette_batch(adata_annot_filter, 'source_label', 'cell_type_semifinal_v2', emb_key_)

#Excluding AcinarS_like (P24-dominant, expected low mixing):
print("ASW (PC=50, excl. AcinarS_like & SMC_like):")
silhouette_batch(
    adata_annot_filter[~adata_annot_filter.obs['cell_type_semifinal_v2'].isin(
        ["AcinarS_like", "SMC_like"])],
    'source_label', 'cell_type_semifinal_v2', emb_key_)

#── 13. PCA variance explained ────────────────────────────────────────────────
print("Observed variance explained (top 20 PCs):")
print(adata_annot_filter.uns['pca']['variance_ratio'][0:20])

plt.bar(range(len(adata_annot_filter.uns['pca']['variance_ratio'])),
        adata_annot_filter.uns['pca']['variance_ratio'])
plt.xlabel('#Principal Components', fontsize=8)
plt.ylabel('Explained Variance', fontsize=8)

#── 14. Save filtered annotated object ────────────────────────────────────────
results_new = "./pNET_comparative_clustering/pNET_harmony_bbknn_filtered_021122.h5ad"
adata_annot_filter.write(results_new)

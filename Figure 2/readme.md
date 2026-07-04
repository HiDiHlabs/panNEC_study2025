This notebook generates panels 2b, 2c & 2d.

2a — Transcription factor regulon activity heatmap generated using pySCENIC (script: analysis/run_PySCENIC.sh). GRN inference was run with GRNBoost2, followed by regulon enrichment using cisTarget (hg19, 500bp/5kb/10kb upstream motif databases) and AUCell scoring. Final analysis done by Dr. Hilmar Berger. 

2b — Dotplot of pancreas-upregulated vs brain-upregulated gene modules across cell states. Gene modules were derived from bulk RNA-seq differential expression between GTEx pancreas and brain samples (DESeq2 LRT; script: analysis/DGE_Brain_Pancreas.R). VST-normalised average expression per GTEx tissue subtype was used to annotate directionality.

2c — Bubble plot showing enrichment of PTF1A, PAX6 and NKX2.2 target gene sets (pancreas vs brain targets) across cell states. Enrichment was assessed using a hypergeometric test against cell state marker genes; scores capped at 95th percentile, FDR-corrected.

2d — Dotplot of proliferation-associated brain-type markers across NE cell states.
Reads pNEC_updated_annot_07082023.h5ad (same input as Figure 1, with identical Cell-states renames applied). Immune cells are excluded for all panels. Cell state marker gene lists are loaded from Supplementary Data 2 (internally named as Supplementary Table 23_panNEC_substate_markers.xlsx).

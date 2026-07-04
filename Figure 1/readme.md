This notebook generates the main Figure 1 panels. It reads the semi-final annotated object (pNEC_updated_annot_07082023.h5ad, output of the Step 5 annotation pipeline in analysis/) and applies three final cell states renames before figure generation:

Neuroendocrine EEC-like01 → Amphicrine progenitor-like
Neuroendocrine EEC-like02 → Neuroendocrine
Neuroendocrine HSP+ → Neuroendocrine HSP+ (hypoxic)

Panels generated: UMAP by PatientID (1b), UMAP by Cell types (1c), diagnostic marker dotplot per Cell type (1d), UMAP by Cell states (1e), patient composition bar plot (1f), and Metascape pathway module matrix plot (1g). Pathway enrichment results are tabulated in Supplementary Table 3. Font was adjusted to Arial in Adobe Illustrator post-export.

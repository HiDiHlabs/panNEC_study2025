# Integration & Annotation Pipeline

This folder documents the computational pipeline used to build the panNEC single-cell atlas — from raw CellRanger output to the annotated object underlying all manuscript figures.

**Manuscript:** Aberrant brain-type neuronal programs in large-cell pancreatic neuroendocrine carcinoma  
**Zenodo:** https://zenodo.org/records/20097469  
**GEO:** [GSE291816](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291816)  
**Contact:** olivia_debnath@dfci.harvard.edu

---

## Pipeline overview

```
CellRanger output (GEO: GSE291816)
        │
        ▼
   Step 0 ── SoupX ambient RNA removal (R)
              Seurat → AnnData conversion (per patient)
        │
        ▼  NET_Soupx_integer_P*.h5ad  (one per patient)
        │
   Step 1 ── Concatenation, QC, normalization, HVG selection,
              Harmony + BBKNN integration, Leiden clustering
        │
        ▼  pNET_harmonization_4000hvg_BC_181022.h5ad
        │
   Step 2 ── Cluster exploration, initial cell type annotation
              (cell_type_semifinal01)
        │
        ▼  pNET_harmony_bbknn_241022.h5ad
        │
   Step 3 ── Refined annotation (cell_type_semifinal02),
              L12 contamination filter, interpatient variability QC
        │
        ▼  pNET_harmony_bbknn_filtered_021122.h5ad
        │
   Step 4 ── Stromal & lymphocyte subclustering,
              iterative label refinement (v3 → v5),
              Seurat export for Logistic Regression DEGs
        │
        ▼  pNET_updated_annotations_180123.h5ad
        │
   Step 5 ── Final label harmonisation →
              'Cell states' (12), 'Cell types' (4), 'PatientID' (5)
        │
        ▼  pNEC_updated_annot_07082023.h5ad  ← input to Figure notebooks
        │
   Figure notebooks
              (3 final Cell state renames → manuscript labels)
        │
        ▼  Zenodo deposit: panNEC_atlas_dataset_raw_20260503.h5ad
```

---

## Step 0 — Preprocessing (SoupX + R → AnnData)

Raw FASTQ files were aligned using CellRanger. Per-patient count matrices were processed in R:

- Ambient RNA contamination was estimated and removed using **SoupX**
- The cleaned integer count matrices were exported from Seurat and converted to `.h5ad` format (one file per patient), with technical replicates of P24 kept separate at this stage

Raw CellRanger matrices are available at GEO (GSE291816); personally identifiable sequence data will be available via EGA.

---

## Steps 1–5 — Integration & annotation scripts

| Script | Input | Output | Key steps |
|--------|-------|--------|-----------|
| `panNEC_chain_integration_script_Step1.py` | Per-patient SoupX h5ad | `pNET_harmonization_4000hvg_BC_181022.h5ad` | Concatenation, MALAT1 removal, re-normalization, QC, cell cycle scoring, HVG selection (4,000 genes, per-batch), PCA, Harmony, BBKNN, Leiden clustering |
| `panNEC_chain_integration_script_Step2_annotation.py` | Step 1 h5ad | `pNET_harmony_bbknn_241022.h5ad` | Cluster exploration at multiple resolutions, initial cell type annotation (`cell_type_semifinal01`) |
| `panNEC_chain_integration_script_Step3_interpatient_variability.py` | Step 2 h5ad | `pNET_harmony_bbknn_filtered_021122.h5ad` | Refined annotation (`cell_type_semifinal02`), partial L12 filter (P18/P19/P21 barcodes), UMAP recompute, batch integration quality metrics |
| `panNEC_chain_integration_script_Step4_subclustering.py` | Step 3 h5ad | `pNET_updated_annotations_180123.h5ad` | Stromal subclustering (Stromal1/2), iterative label refinement (v3→v5), Lymphocyte subclustering, Seurat export |
| `panNEC_chain_integration_script_Step5_final_annotation.py` | Step 4 h5ad | `pNEC_updated_annot_07082023.h5ad` | L12 cluster fully excluded; final label harmonisation → `Cell states` (12), `Cell types` (4), `PatientID` (5); 3 `Cell states` labels receive a final rename in the Figure 1 notebook |

---

## Final annotation layers

Three `Cell states` labels receive a final rename in the Figure 1 notebook before figures are generated (noted below). All other labels are stable from Step 5 onward.

### `Cell states` — 12 states

| Cell state | Broad class |
|------------|-------------|
| Neuroendocrine ¹ | NE |
| Neuroendocrine proliferating | NE |
| Neuroendocrine HSP+ (hypoxic) ² | NE |
| Neuroendocrine stromal-like | NE |
| Amphicrine acinar01 | Amphicrine |
| Amphicrine acinar02 | Amphicrine |
| Amphicrine acinar03 | Amphicrine |
| Amphicrine acinar proliferating | Amphicrine |
| Amphicrine progenitor-like ³ | Amphicrine |
| Stroma (normal) | Stroma |
| Lymphocytes | Immune |
| Macrophages | Immune |

¹ Renamed from `Neuroendocrine EEC-like02` in Figure 1 notebook  
² Renamed from `Neuroendocrine HSP+` in Figure 1 notebook  
³ Renamed from `Neuroendocrine EEC-like01` in Figure 1 notebook  

### `Cell types` — 4 coarse types

Neuroendocrine / Amphicrine acinar / Stroma (normal) / Immune

### `PatientID`

P1 / P2 / P3 / P4 / P5 (acinar)

---

## Integration strategy

Rather than a single batch correction method, we combined two complementary approaches:

**Harmony** adjusts PCA coordinates iteratively until cluster centroids are patient-invariant, without touching the count matrix. This gives a corrected embedding (`X_pca_harmony`) while leaving the expression data intact for downstream analysis.

**BBKNN** then uses those corrected PCs to build a neighbourhood graph that is explicitly balanced across patients — each cell gets 5 nearest neighbours from each patient separately, preventing the largest patient (P24, ~19k nuclei) from monopolising the graph. The `trim=80` parameter removes weak long-range edges to reduce noise.

---

## Dependencies

```
scanpy >= 1.8
anndata
harmonypy
bbknn
pandas
numpy
matplotlib
seaborn
scikit-learn
```

For Step 0: `Seurat`, `SoupX`, `SeuratDisk`

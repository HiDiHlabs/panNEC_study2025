library(Seurat)
library(SeuratDisk)

setwd("Final_scripts/scRNA_preprocessing")

# source: https://www.cellphonedb.org/faq-and-troubleshooting
# take raw data and normalise it

all_samples = c("P7","P18","P19","P21","P24ab")

load(file.path("../Data/Processed", "NET_SoupX_downsampled_6best_and_all_Tosti_pancreas_BBKNN.RData" ))
scData = scData.combined
scData@meta.data = scData@meta.data[, c("orig_cell_id","seurat_clusters")]
SaveLoom(scData, filename="../Data/Processed/NET_all_NET_with_Tosti_complete_v2.loom", overwrite=T)

quit("no")


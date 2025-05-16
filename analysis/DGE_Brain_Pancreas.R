library(data.table)
library(DESeq2)
# GTEX download page: https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
# gene reads downloaded from here: 
# https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
d <- fread("/mnt/data/data_1/DataSets/GTEX/Expression/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
# Meta-data obtained from https://gtexportal.org/home/downloads/adult-gtex/metadata
# Download link: https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
sample_description <- fread("/mnt/data/data_1/DataSets/GTEX/Meta/annotations_v8_GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
sample_description <- subset(sample_description, SAMPID %in% colnames(d))

sel_samples = subset(sample_description, SMTS %in% c("Pancreas","Brain"))

expr_mat <- as.matrix(d[, sel_samples$SAMPID, with=F])
rownames(expr_mat) <- d$Name
colData <- sel_samples[, c("SAMPID","SMTS","SMTSD")]
rownames(colData) <- colData$SAMPID

anno <- d[, 1:2]
setkey(anno, "Name")
#row_sums <- rowSums(expr_mat)

dds <- DESeq2::DESeqDataSetFromMatrix(expr_mat, colData = colData, design = ~ SMTS)
dds_LRT <- DESeq(dds, fitType = "glmGamPoi", full = ~ SMTS, reduced = ~ 1, test = "LRT")
res <- results(dds_LRT, name = "SMTSPancreas") |> as.data.frame()
#resLFC <- lfcShrink(dds_LRT, coef = "SMTSPancreas", type="apeglm")
res$gene <- anno[rownames(res)]$Description
res$Ensembl_gene <- rownames(res)
res2 <- subset(res, abs(log2FoldChange) < 500 & !is.na(log2FoldChange))
setorder(res2, -log2FoldChange)

fwrite(res2, file = "/mnt/data/data_1/scSeq_NET/GTEX/Pancreas_vs_brain_DGE.tsv", sep="\t", quote = F)

norm <- vst(dds_LRT) |> assay()
rownames(norm) <- anno[rownames(norm)]$Description

expr_by_subtype <- limma::avearrays(norm, sel_samples$SMTSD)
tmp <- as.data.table(cbind(gene=rownames(expr_by_subtype), expr_by_subtype))
tmp$gene <- rownames(expr_by_subtype)
fwrite(tmp, file = "/mnt/data/data_1/scSeq_NET/GTEX/Pancreas_Brain_Subtype_Avg_VST.tsv", sep="\t", quote = F)


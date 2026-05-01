# ================================
# Load libraries
# ================================
library(DESeq2)
library(ggplot2)
library(pheatmap)

# ================================
# Load count matrix
# ================================
# genes as rows, samples as columns
counts <- read.csv("gene_count_matrix.csv", row.names = 1)
head(counts)
# ================================
# Sample information
# ================================
sample_info <- read.csv("metadata.csv",row.names=1, sep="\t")
head(sample_info)
#sample_info <- data.frame(
#  row.names = c("SRR21496993", "SRR21496994",
#                "SRR21496997", "SRR21496998"),
#  condition = c("diseased", "control",
#                 "control", "diseased")
#)

# ================================
# Create DESeq2 dataset
# ================================
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition
)

# ================================
# Manual DESeq2 steps
# ================================
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)

# ================================
# Extract results
# ================================
res <- results(dds)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

write.csv(res_df, "DESeq2_results_95.csv")
print(head(res, 20))

# ================================
# rlog transformation
# ================================
rld <- rlog(dds, blind = TRUE)

# ================================
# PCA Plot
# ================================
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("Plot_PCA_95.pdf", width = 6, height = 6)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA of Samples Based on Gene Expression")
dev.off()

# ================================
# Volcano Plot
# ================================
res_df <- res_df[!is.na(res_df$padj), ]

res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

pdf("Volcano_Plot_95.pdf", width = 7, height = 6)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
dev.off()

# ================================
# Heatmap of Top 50 DE Genes
# ================================
topGenes <- head(order(res$padj), 50)
mat <- assay(rld)[topGenes, ]

# Scale gene expression
mat <- t(scale(t(mat)))

# Sample annotation
annotation_col <- data.frame(
  condition = sample_info$condition
)
rownames(annotation_col) <- rownames(sample_info)

pdf("Heatmap_Top50_DEGs_95.pdf", width = 7, height = 9)
pheatmap(mat,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 6,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Top 50 Differentially Expressed Genes")
dev.off()


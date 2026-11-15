# Assignment 2: miRNA Differential Expression, SARS-CoV-2 vs Mock
library(DESeq2)

setwd("/N/scratch/vgutta/asg2")
dir.create("results", showWarnings = FALSE)

counts <- read.delim("counts/mirna_counts_clean.tsv", row.names = 1, check.names = FALSE)
meta   <- read.csv("counts/metadata.csv", row.names = 1)

cat("Samples match:", all(rownames(meta) == colnames(counts)), "\n")
cat("miRNAs annotated:", nrow(counts), "\n")

meta$condition <- factor(meta$condition, levels = c("Mock", "SARS"))
meta$time      <- factor(meta$time, levels = c("24h", "72h"))

# ---- Pre-filter: keep miRNAs with >=10 counts in at least 3 samples ----
# (3 = smallest group size; standard DESeq2 practice)
keep <- rowSums(counts >= 10) >= 3
cat("miRNAs after filtering:", sum(keep), "\n\n")
counts <- counts[keep, ]

# ---- Comparison 1: Mock vs SARS ----
dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ time + condition)
dds <- DESeq(dds)

res1 <- results(dds, contrast = c("condition", "SARS", "Mock"))
res1 <- res1[order(res1$padj), ]
deg1 <- subset(res1, padj < 0.05)

cat("=== Mock vs SARS ===\n")
cat("Significant DEGs:", nrow(deg1), "\n")
cat("  Up in SARS:  ", sum(deg1$log2FoldChange > 0), "\n")
cat("  Down in SARS:", sum(deg1$log2FoldChange < 0), "\n")
print(head(as.data.frame(deg1)[, c("baseMean","log2FoldChange","padj")], 10))
cat("\n")

write.csv(as.data.frame(res1), "results/DEG_Mock_vs_SARS_full.csv")
write.csv(as.data.frame(deg1), "results/DEG_Mock_vs_SARS_significant.csv")

# ---- Comparison 2: SARS 24h vs 72h ----
dds_s <- dds[, dds$condition == "SARS"]
dds_s$condition <- droplevels(dds_s$condition)
design(dds_s) <- ~ time
dds_s <- DESeq(dds_s)

res2 <- results(dds_s, contrast = c("time", "72h", "24h"))
res2 <- res2[order(res2$padj), ]
deg2 <- subset(res2, padj < 0.05)

cat("=== SARS 24h vs 72h ===\n")
cat("Significant DEGs:", nrow(deg2), "\n")
cat("  Up at 72h:  ", sum(deg2$log2FoldChange > 0), "\n")
cat("  Down at 72h:", sum(deg2$log2FoldChange < 0), "\n")
print(head(as.data.frame(deg2)[, c("baseMean","log2FoldChange","padj")], 10))
cat("\n")

write.csv(as.data.frame(res2), "results/DEG_SARS_24h_vs_72h_full.csv")
write.csv(as.data.frame(deg2), "results/DEG_SARS_24h_vs_72h_significant.csv")

# ---- Normalized expression matrix (Submission 1) ----
norm <- counts(dds, normalized = TRUE)
write.csv(norm, "results/miRNA_normalized_expression.csv")

# ---- PCA (using rlog, robust for small/low-count datasets) ----
rld <- rlog(dds, blind = FALSE)
pdf("results/PCA_plot.pdf", width = 7, height = 5)
print(DESeq2::plotPCA(rld, intgroup = c("condition", "time")))
dev.off()

cat("Done. Results in results/\n")

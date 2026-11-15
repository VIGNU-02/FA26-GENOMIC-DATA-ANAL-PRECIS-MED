library(ggplot2)
setwd("/N/scratch/vgutta/asg2")

make_volcano <- function(csv, title, outfile, up_lab, down_lab) {
  res <- read.csv(csv, row.names = 1)
  res <- res[!is.na(res$padj), ]
  res$sig <- "Not significant"
  res$sig[res$padj < 0.05 & res$log2FoldChange > 0] <- up_lab
  res$sig[res$padj < 0.05 & res$log2FoldChange < 0] <- down_lab
  res$label <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1,
                      rownames(res), NA)

  cols <- setNames(c("grey70", "#D62728", "#1F77B4"),
                   c("Not significant", up_lab, down_lab))

  p <- ggplot(res, aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_colour_manual(values = cols, name = NULL) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
    ggrepel::geom_text_repel(aes(label = label), size = 2.8,
                             max.overlaps = 20, show.legend = FALSE) +
    labs(title = title,
         x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"adjusted p-value")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  ggsave(outfile, p, width = 7, height = 6, dpi = 300)
  cat("Wrote", outfile, " | significant:",
      sum(res$padj < 0.05), "\n")
}

make_volcano("results/DEG_Mock_vs_SARS_full.csv",
             "miRNA differential expression: SARS-CoV-2 vs Mock",
             "results/volcano_Mock_vs_SARS.png",
             "Up in SARS-CoV-2", "Down in SARS-CoV-2")

make_volcano("results/DEG_SARS_24h_vs_72h_full.csv",
             "miRNA differential expression: SARS-CoV-2 72h vs 24h",
             "results/volcano_SARS_24h_vs_72h.png",
             "Up at 72h", "Down at 72h")

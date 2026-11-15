library(multiMiR); library(clusterProfiler); library(org.Hs.eg.db)
setwd("/N/scratch/vgutta/asg2")

# Background: targets of all miRNAs expressed in this dataset
expressed <- rownames(read.delim("counts/mirna_counts_clean.tsv", row.names=1, check.names=FALSE))
cnt <- read.delim("counts/mirna_counts_clean.tsv", row.names=1, check.names=FALSE)
expressed <- rownames(cnt)[rowSums(cnt >= 10) >= 3]
cat("Expressed miRNAs for background:", length(expressed), "\n")

bg <- get_multimir(org="hsa", mirna=expressed, table="mirtarbase", summary=FALSE)@data
bg <- bg[grepl("Functional MTI", bg$support_type) & !grepl("Weak", bg$support_type), ]
universe <- unique(bg$target_symbol); universe <- universe[universe!="" & !is.na(universe)]
cat("Background universe size:", length(universe), "\n")

strong <- function(mirs){
  if(!length(mirs)) return(character(0))
  d <- get_multimir(org="hsa", mirna=mirs, table="mirtarbase", summary=FALSE)@data
  d <- d[grepl("Functional MTI", d$support_type) & !grepl("Weak", d$support_type), ]
  t <- unique(d$target_symbol); t[t!="" & !is.na(t)]
}

run <- function(f, label){
  degs <- read.csv(f, row.names=1)
  for(dir in c("UP","DOWN")){
    mirs <- rownames(degs)[if(dir=="UP") degs$log2FoldChange>0 else degs$log2FoldChange<0]
    tg <- strong(mirs)
    cat("\n===", label, dir, "=== miRNAs:", length(mirs), " targets:", length(tg), "\n")
    if(length(tg) < 10){ cat("too few targets\n"); next }
    ego <- enrichGO(gene=tg, universe=universe, OrgDb=org.Hs.eg.db,
                    keyType="SYMBOL", ont="BP", pAdjustMethod="BH",
                    pvalueCutoff=0.05, qvalueCutoff=0.05)
    df <- as.data.frame(ego)
    cat("GO BP terms:", nrow(df), "\n")
    if(nrow(df)>0) print(head(df[,c("Description","GeneRatio","p.adjust")], 10))
    write.csv(df, paste0("results/GOfinal_", label, "_", dir, ".csv"), row.names=FALSE)
  }
}
run("results/DEG_Mock_vs_SARS_significant.csv","Mock_vs_SARS")
run("results/DEG_SARS_24h_vs_72h_significant.csv","SARS_24h_vs_72h")

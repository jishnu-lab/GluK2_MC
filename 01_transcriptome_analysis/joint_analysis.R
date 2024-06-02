### Load Libraries ----
library(biomaRt)
library(sleuth)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(matrixStats)
library(pheatmap)
library(VennDiagram)
library(stringr)
library(ReactomePA)
library(clusterProfiler)
# library(WebGestaltR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(tidyverse)
# library(futile.logger)
# flog.threshold(INFO)
# flog.layout(layout.format("[%.Time%] %.Level%: %.Message%"))

### Loading up DEG's and making list of common values ----
DEG_Up_MRG <- read.xlsx('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Shiqun/02-03-2024-DE-Analysis/sig_DE_MRG_vs_WT.xlsx',sheet = "sigDEG_Up")
DEG_Dn_SYM <- read.xlsx('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/sig_DE_GLU_vs_UNT.xlsx',sheet = "sigDEG_Dn")
DEG_Dn_GLU <- read.xlsx('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/sig_DE_SYM_vs_UNT.xlsx',sheet = "sigDEG_Dn")

DET_Up_MRG <- read.xlsx('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Shiqun/02-03-2024-DE-Analysis/sig_DE_MRG_vs_WT.xlsx',sheet = "sigDET_Up")
DET_Dn_SYM <- read.xlsx('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/sig_DE_GLU_vs_UNT.xlsx',sheet = "sigDET_Dn")
DET_Dn_GLU <- read.xlsx('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/sig_DE_SYM_vs_UNT.xlsx',sheet = "sigDET_Dn")

common_values<- Reduce(intersect, list(DEG_Up_MRG$gene_name, DEG_Dn_SYM$gene_name, DEG_Dn_GLU$gene_name))
write.csv(common_values, file = "/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/common_genes.csv",row.names = FALSE)

target_ids_SZ <- Reduce(union, list(DET_Up_MRG[DET_Up_MRG$gene_name %in% common_values,]$gene_name))
target_ids_YZ <-  Reduce(union, list(DET_Dn_SYM[DET_Dn_SYM$gene_name %in% common_values,]$target_id,DET_Dn_GLU[DET_Dn_GLU$gene_name %in% common_values,]$gene_name))

# common_values_SYM <- Reduce(intersect, list(DEG_Up_MRG$gene_name, DEG_Dn_SYM$gene_name))
# common_values_GLU <- Reduce(intersect, list(DEG_Up_MRG$gene_name, DEG_Dn_GLU$gene_name))
# common_values_weak <- Reduce(union, list(common_values_SYM, common_values_GLU))
# write.csv(common_values_weak, file = "/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/common_genes_weak.csv",row.names = FALSE)

### Functions used ----
z_score_TPM <- function(path, target_ids) {
  TPM <- read.csv(path)
  TPM_subset <- TPM[TPM$X %in% target_ids,]
  TPM_mean <- TPM_subset %>% group_by(gene_name) %>% summarise(across(c(names(TPM)[-c(1, 2)]),\(x) sum(x, na.rm = TRUE)))
  print(dim(TPM_mean))
  z_scores <- apply(TPM_mean[-c(1)], 1, function(x) (x - mean(x)) / sd(x))
  TPM_z_scores <- as.data.frame(t(z_scores))
  TPM_z_scores$gene_name <- TPM_mean$gene_name
  return (TPM_z_scores)
}

heatmap_plot <- function(df,custom_order,name,width = 4, height = 6) {
  TPM_Joint <- df
  TPM_Joint_matrix <- as.matrix(TPM_Joint[,custom_order])
  colnames(TPM_Joint_matrix) = custom_order
  rownames(TPM_Joint_matrix) = TPM_Joint$gene_name
  color_breaks <- seq(-2.5, 2.5, length.out = 101)
  pheatmap(TPM_Joint_matrix, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE, angle_col = "45",
          fontsize = 6, cellwidth = 15, cellheight = 10, main = name,color = colorRampPalette(c("blue", "white", "red"))(100),breaks = color_breaks,
          filename = paste0("/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/",name,"_heatmap_plot_targets_only.pdf"), width = width, height = height, units = "in", res = 300)
}

### Calculating matrices and plotting heatmap ----
TPM_SZ_zs_target <- z_score_TPM('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Shiqun/02-03-2024-DE-Analysis/tpm_names.csv',target_ids_SZ)
TPM_YZ_zs_target <- z_score_TPM('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/tpm_names.csv',target_ids_YZ)

common_values_with_transcripts<- Reduce(intersect, list(TPM_SZ_zs_target$gene_name, TPM_YZ_zs_target$gene_name))

TPM_SZ_zs <- TPM_SZ_zs_target[TPM_SZ_zs_target$gene_name %in% common_values_with_transcripts,]
TPM_YZ_zs <- TPM_YZ_zs_target[TPM_YZ_zs_target$gene_name %in% common_values_with_transcripts,]

custom_order <- c(names(TPM_SZ_zs)[c(5:8)],names(TPM_SZ_zs)[c(1:4)])
heatmap_plot(TPM_SZ_zs,custom_order,"TPM_SZ_zs")
custom_order <- c(names(TPM_YZ_zs)[c(9:12)],names(TPM_YZ_zs)[c(1:8)])
heatmap_plot(TPM_YZ_zs,custom_order,"TPM_YZ_zs")

# ### Calculating matrices and plotting heatmap ----
# TPM_SZ_zs_weak <- z_score_TPM('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Shiqun/02-03-2024-DE-Analysis/tpm_names.csv',common_values_weak)
# TPM_YZ_zs_weak <- z_score_TPM('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/tpm_names.csv',common_values_weak)
# custom_order <- c(names(TPM_SZ_zs_weak)[c(5:8)],names(TPM_SZ_zs_weak)[c(1:4)])
# heatmap_plot(TPM_SZ_zs_weak,custom_order,"TPM_SZ_zs_weak",width = 4, height = 25)
# custom_order <- c(names(TPM_YZ_zs_weak)[c(9:12)],names(TPM_YZ_zs_weak)[c(1:8)])
# heatmap_plot(TPM_YZ_zs_weak,custom_order,"TPM_YZ_zs_weak",width = 4, height = 25)

### Calculating matrices and plotting heatmap ----
lf_z45_run5 <- read.table("/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/slide_analysis/slide_runs/1k_features_deep/run_5/0.01_1_out/feature_list_Z45.txt", header = TRUE, sep = "\t")
common_values_lf45 <- Reduce(intersect, list(lf_z45_run5$names, lf_z45_run5$names))
TPM_SZ_zs_lfz45 <- z_score_TPM('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Shiqun/02-03-2024-DE-Analysis/tpm_names.csv',common_values_lf45)
TPM_YZ_zs_lfz45 <- z_score_TPM('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/tpm_names.csv',common_values_lf45)
custom_order <- c(names(TPM_SZ_zs_lfz45)[c(5:8)],names(TPM_SZ_zs_lfz45)[c(1:4)])
heatmap_plot(TPM_SZ_zs_lfz45,custom_order,"TPM_SZ_zs_lfz45",width = 4, height = 25)
custom_order <- c(names(TPM_YZ_zs_lfz45)[c(9:12)],names(TPM_YZ_zs_lfz45)[c(1:8)])
heatmap_plot(TPM_YZ_zs_lfz45,custom_order,"TPM_YZ_zs_lfz45",width = 4, height = 25)

### Calculating matrices and plotting heatmap ----
lf_z32_run5 <- read.table("/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/slide_analysis/slide_runs/1k_features_deep/run_5/0.01_1_out/feature_list_Z32.txt", header = TRUE, sep = "\t")
common_values_lf32 <- Reduce(intersect, list(lf_z32_run5$names, lf_z32_run5$names))
TPM_SZ_zs_lfz32 <- z_score_TPM('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Shiqun/02-03-2024-DE-Analysis/tpm_names.csv',common_values_lf32)
TPM_YZ_zs_lfz32 <- z_score_TPM('/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/tpm_names.csv',common_values_lf32)
custom_order <- c(names(TPM_SZ_zs_lfz32)[c(5:8)],names(TPM_SZ_zs_lfz32)[c(1:4)])
heatmap_plot(TPM_SZ_zs_lfz32,custom_order,"TPM_SZ_zs_lfz32",width = 4, height = 25)
custom_order <- c(names(TPM_YZ_zs_lfz32)[c(9:12)],names(TPM_YZ_zs_lfz32)[c(1:8)])
heatmap_plot(TPM_YZ_zs_lfz32,custom_order,"TPM_YZ_zs_lfz32",width = 4, height = 25)

# ### ORA Analysis ----
# ORA_analysis <- function(gene_list,g1){
#   genes <- gene_list
#   entrez_genes <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

#   # Reactome pathway over-representation analysis
#   reactome <- enrichPathway(gene = entrez_genes$ENTREZID, organism = "mouse", pvalueCutoff=pvalueCutoff, readable=TRUE)
#   ORA[[paste0(g1, "_Reactome")]] <- reactome

#   # # KEGG pathway over-representation analysis
#   # kegg <- enrichKEGG(gene = entrez_genes$ENTREZID, organism = "mmu",  pvalueCutoff=pvalueCutoff, use_internal_data = FALSE)
#   # kegg2 <- setReadable(kegg, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
#   # kegg2@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", kegg2@result$Description, fixed = T)    
#   # ORA[[paste0(g1, "_KEGG")]] <- kegg2

#   for(a1 in names(ORA)){
#     reacORkegg <- ORA[[a1]]
#     if(nrow(reacORkegg) !=0){
#       # Create directories for each result
#       ORA_path <- file.path(cp_path, a1)
#       dir.create(ORA_path, recursive=TRUE)
#       ORA_df <- as.data.frame(reacORkegg)
#       # Convert 'GeneRatio' and 'BgRatio' from character to numeric
#       ORA_df$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", ORA_df$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", ORA_df$GeneRatio, perl=T) )
#       ORA_df$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", ORA_df$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", ORA_df$BgRatio, perl=T) )
#       # Now perform the division to create 'RatioOfRatios'
#       ORA_df <- transform(ORA_df, EnrichmentRatio = GeneRatio / BgRatio)
#       write.xlsx(ORA_df, file.path(ORA_path, paste0(a1, ".xlsx")))
      
#       # Plot enrichment result - dotplot
#       fit <- dotplot(reacORkegg, showCategory = 40)
#       ggsave(file.path(ORA_path, paste0(a1, "_dotplot.pdf")), plot = fit, dpi = 300, width = 15, height = 12, units = "in", device="pdf")
#       plots_list[[a1]] <- fit
#     }
    
#     if (nrow(reacORkegg) > 40){
#       # Weighted Set Cover of geneSets
#       weightedPath <- file.path(cp_path, "weightedSetCover", a1)
#       dir.create(weightedPath, recursive=TRUE)
      
#       setCoverNum = abs(0.40*(nrow(ORA_df)))
#       nThreads = 4
#       idsInSet <- sapply(ORA_df$geneID, strsplit, split="/")
#       names(idsInSet) <- ORA_df$ID    
#       minusLogP <- -log(ORA_df$pvalue)
#       minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
#       wscRes <- weightedSetCover(idsInSet=idsInSet, costs=(1 / minusLogP), topN=setCoverNum, nThreads=nThreads)
      
#       wscRes_full <- ORA_df[c(match(wscRes$topSets, ORA_df$ID)),]
#       wscRes_full <- wscRes_full[order(wscRes_full$p.adjust), ]
#       write.xlsx(wscRes_full, file.path(weightedPath, paste0(a1, ".xlsx")))
      
#       # Plot enrichment result of weighted set cover - dotplot
#       reacORkegg2 <- reacORkegg
#       reacORkegg2@result <- wscRes_full
      
#       fit <- dotplot(reacORkegg2, showCategory = 20)
#       ggsave(file.path(weightedPath, paste0(a1, "_dotplot.pdf")), plot = fit, dpi = 300, width = 17, height = 12, units = "in", device="pdf")
#       plots_list[[paste0(a1, "weightSet")]] <- fit
#     }
#   }
# }
# cp_path <- file.path("/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/joint_analysis", "pathway_analysis_ORA")
# dir.create(cp_path, recursive=TRUE)
# pvalueCutoff <-0.05
# ORA <- list()
# plots_list <- list()
# ORA_analysis(common_values,"joint_analysis")
# ORA_analysis(common_values_weak,"joint_analysis_weak")


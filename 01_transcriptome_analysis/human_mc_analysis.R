# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap) 

setwd("/ix/djishnu/Swapnil/kaplanAnalysis/GluK2_MC/")

# Load the data
data <- read.csv("03_relevent_data_objects/01_transcriptome/human_mast_cells/CELLxGENE_gene_expression_082524.csv",comment.char = "#")

# Combine Tissue and Cell Type into a new column tissue_celltype
heatmap_data <- data  %>% 
        select(Tissue, Cell.Type , Gene.Symbol, Expression) %>%
        filter(Cell.Type != "aggregated") %>%
        mutate(Expression = replace_na(Expression, 0), tissue_celltype = paste(Tissue, Cell.Type, sep = "_")) %>%
        select(tissue_celltype, Gene.Symbol, Expression) %>%
        pivot_wider(names_from = Gene.Symbol, values_from = Expression)



TPM_Joint_matrix <- 10^as.matrix(heatmap_data[-1])
rownames(TPM_Joint_matrix) = heatmap_data$tissue_celltype
pheatmap(TPM_Joint_matrix, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE, angle_col = "90",
         fontsize = 6, cellwidth = 15, cellheight = 10, color = colorRampPalette(c("white", "red"))(100),
         units = "in", res = 300)
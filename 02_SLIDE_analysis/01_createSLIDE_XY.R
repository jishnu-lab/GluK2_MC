library(readr)
library(dplyr)

setwd(".")
### Loading up DEG's and making list of common values ----
TPM_YZ <- read_csv("/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/02-03-2024-DE-Analysis/tpm_names.csv")
TPM_SZ <- read_csv("/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Shiqun/02-03-2024-DE-Analysis/tpm_names.csv")
TPM_joint <- inner_join(TPM_YZ, TPM_SZ, by = "...1")
TPM_joint <- TPM_joint %>% rename(transcript_id = "...1", gene_name = "gene_name.x") %>% select(-gene_name.y)
TPM_joint <- TPM_joint[complete.cases(TPM_joint), ] # Dropping rows with NA values
TPM_joint <- TPM_joint %>% filter(!(grepl("^MT-", TPM_joint$gene_name, ignore.case = TRUE) | 
                                    grepl("^RPS", TPM_joint$gene_name, ignore.case = TRUE) | 
                                    grepl("^RPL", TPM_joint$gene_name, ignore.case = TRUE) ))
                                    # |grepl("^GM\\d+", TPM_joint$gene_name, ignore.case = TRUE) )) # Removing mitochondrial and ribosomal genes
TPM_joint_aggr <- aggregate(. ~ gene_name, data = TPM_joint[,-c(1)], FUN = mean, na.rm = TRUE)
rownames(TPM_joint_aggr) <- TPM_joint_aggr[, 1]
TPM_joint_aggr <- TPM_joint_aggr[, -1]

### Filter on variances across the samples ----
# Make sure while filtering if there are some genes of interest they still remain.
# This is not wrong since anyway filtering on variable genes is a heuristic.
# Add back the genes of interest if they are not in the top variable genes.
gene_variances <- apply(TPM_joint_aggr, 1, var)
top_genes <- sort(gene_variances, decreasing = TRUE)[1:1000]
tpm_matrix_top_genes <- TPM_joint_aggr[names(top_genes), ]

gene_joint_aggr <- as.data.frame(t(as.matrix(tpm_matrix_top_genes)))
colnames(gene_joint_aggr) <- rownames(tpm_matrix_top_genes)
gene_joint_aggr <- cbind(X = colnames(tpm_matrix_top_genes),gene_joint_aggr)
Y <- cbind(X = colnames(tpm_matrix_top_genes),Y=c(1,1,1,1,1,1,1,1,0,0,0,0,-1,-1,-1,-1,0,0,0,0))
write.csv(gene_joint_aggr, "../03_relevent_data_objects/02_SLIDE/X_1k.csv", row.names = FALSE)
write.csv(Y, "../03_relevent_data_objects/02_SLIDE/Y.csv", row.names = FALSE)


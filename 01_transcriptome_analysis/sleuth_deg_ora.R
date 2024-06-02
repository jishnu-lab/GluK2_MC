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

### Load the configuration from the YAML file ----
setwd(".") # Set working directory to the location of the script
args <- commandArgs(trailingOnly = TRUE)
configFilePath <- if(length(args) > 0) args[1] else stop("No YAML file path provided.")
# configFilePath <- file.path("/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/Youran/config.yaml")
config <- yaml::read_yaml(configFilePath)

# Basic file path setup
proj_path <- config$proj_path
dataDir <- config$dataDir
metaData <- config$metadataPath
current_date <- Sys.Date()
formatted_date <- format(current_date, "%d-%m-%Y")
experiment <- paste0(formatted_date, "-DE-Analysis")
experiment_path <- file.path(proj_path, experiment)
outPath <- experiment_path

# Variables for sleuth
skipLines <- config$skipLines
sampleCol <- config$sampleCol
factors <- config$factorCol
factorOfInt <- config$factorOfInt
refLevel <- config$refLevel
filter_target_id <- config$filter_target_id
sleuth_test <- config$sleuth_test
geneList <- config$geneList
drop_dupGenes <- config$drop_dupGenes

# Setting up cutoffs
qval <- config$qval # used for DEG analysis
bval <- config$bval # used for DEG analysis
pvalueCutoff <- config$pathway_pvalueCutoff # used when performing PATHWAY ANALYSIS

# Setting pathways to run
runReactome <- config$runReactome
runKEGG <- config$runKEGG

### DEG analysis : Preliminary variable setting based on inputs ----
# Function to analyze Up and Down Transcript/Genes in case of wald test
load(file = file.path(outPath, "sleuth_object.RData"))
load(file = file.path(outPath, "sleuth_results_noDEG.RData"))
normCounts <- result$normCounts
normCounts_names <- result$normCounts_names
sampleDF <- result$sampleDF

processRegulation <- function(df, baseName, results) {
  # Separate up- and down-regulated entries
  df_Up <- df[df$b > 0, ]
  df_Dn <- df[df$b < 0, ]
  # Remove entries with blank gene names
  df_Up <- df_Up[!(df_Up$gene_name == ""), ]
  df_Dn <- df_Dn[!(df_Dn$gene_name == ""), ]
  # Construct new names for up- and down-regulated groups
  name_Up <- paste(baseName, "_Up", sep="")
  name_Dn <- paste(baseName, "_Dn", sep="")
  # Save the processed DataFrames into the results list
  results[[name_Up]] <- df_Up
  results[[name_Dn]] <- df_Dn
  
  # Return a list containing both up and down regulated DataFrames
  return(list(up = df_Up, down = df_Dn, result=results))
}

# Creating empty data for saving the result for each factor
if(sleuth_test == "wt")
{ # Wald Test
  summaryDF <- data.frame(Comparisons=character(), DE_Genes=numeric(), Up_Reg_Genes=numeric(), Down_Reg_Genes=numeric(),
                          DE_Transcripts=numeric(), Up_Reg_Transcripts=numeric(), Down_Reg_Transcripts=numeric())
  all_sigDETs <- data.frame(matrix(nrow = 0, ncol = 13))
  all_sigDEGs <- data.frame(matrix(nrow = 0, ncol = 7))
} 
else 
{ # Likelihood Ratio Test (LRT)
  summaryDF <- data.frame(Comparison=character(), DE_Genes=numeric(), DE_Transcripts=numeric())
  all_sigDETs <- data.frame(matrix(nrow = 0, ncol = 14))
  all_sigDEGs <- data.frame(matrix(nrow = 0, ncol = 6))
}

# Create list of factors/levels to test
testsToRun <- colnames(so$fits$full$design_matrix)[2:ncol(so$fits$full$design_matrix)]
result[["tests"]] <- testsToRun
if(!is.null(factorOfInt))
{
  testsToRun <- colnames(so$fits$reduced$design_matrix)[2:ncol(so$fits$reduced$design_matrix)]
  result[["tests"]] <- testsToRun
}

### Running DEG analysis : For all possible combinations ----
for(i2 in 1:length(testsToRun))
{
  ## Common steps regardless of sleuth test type
  # Variables defining the writing path of files
  t1 <- testsToRun[i2]
  c1 <- gsub(factors, "", t1)
  if(!is.null(factorOfInt))
  {
    c1 <- gsub(factorOfInt, "", t1)
  }
  fr1 <- gsub(c1, "", t1)
  r1 <- refLevel[factors == fr1]
  compStr <- paste(c1, "_vs_", r1,sep="")
  
  ## Separate steps depending on sleuth test type
  # Running DEG test
  if(sleuth_test == "wt")
  {
    # Run Wald test for each factor-level listed
    so <- sleuth_wt(so, which_beta = t1, which_model="full")
    # Delete negative mean_obs values (for RIP-Seq analyses)
    mean_obs_pos <- so$tests$wt$full[[t1]][so$tests$wt$full[[t1]]$mean_obs > 0, ]
    so$tests$wt$full[[t1]] <- mean_obs_pos
    # Extract results from sleuth object
    deTrans <- sleuth_results(so, test = t1, test_type = "wt", which_model = "full", pval_aggregate = FALSE)
    deGenes <- sleuth_results(so, test = t1, test_type = "wt", which_model = "full", pval_aggregate = TRUE)
  } 
  else 
  { ## LRT test
    # Run LRT(likelihood ratio test (LRT)) test if specifies for each factor-level listed
    so <- sleuth_lrt(so, null_model = "reduced", alt_model = "full") # Run a likelihood ratio test (LRT) between the two models
    # Delete negative mean_obs values (for RIP-Seq analyses)
    mean_obs_pos <- so$tests$lrt$`reduced:full`[so$tests$lrt$`reduced:full`$mean_obs > 0, ]
    so$tests$lrt$`reduced:full` <- mean_obs_pos
    # Extract results from sleuth object
    deTrans <- sleuth_results(so, test = "reduced:full", test_type = "lrt", pval_aggregate = FALSE)
    deGenes <- sleuth_results(so, test = "reduced:full", test_type = "lrt", pval_aggregate = TRUE)
  } 
  
  ## Common steps regardless of sleuth test type
  # remove missing values (NA) in de_transcript data frame 
  deTrans2 <- deTrans[complete.cases(deTrans), ]
  deGenes2 <- deGenes[complete.cases(deGenes), ]
  # if specified, remove duplicate genes (genes w/ multiple transcripts mapped to them); keep transcript with lowest qval (most significant)
  if(drop_dupGenes)
  { 
    deTrans2 <- deTrans2[order(deTrans2[,'gene_name'], deTrans2[,'qval']), ]
    deTrans2 <- deTrans2[!duplicated(deTrans2$gene_name), ]
    deTrans2 <- deTrans2[order(deTrans2[,'qval']), ]    #re-sort by qval
  }
  
  deTrans2Name <- paste(compStr, "_allDET", sep="")
  result[[deTrans2Name]] <- deTrans2
  deGenes2Name <- paste(compStr, "_allDEG", sep="")
  result[[deGenes2Name]] <- deGenes2
  
  ## Separate steps depending on sleuth test type
  # Transcripts/genes with *significant* diff. expression
  if(sleuth_test == "wt")
  {
    sigDET <- deTrans2[deTrans2$qval < qval & abs(deTrans2$b) > bval, ]
    sigDETp <- deTrans2[deTrans2$pval < qval & abs(deTrans2$b) > bval, ]
    # Aggregating b-value per gene
    agg_bval <- aggregate(b ~ ensembl_gene, data = deTrans2[ ,c("ensembl_gene", "b")], sum)
    deGenes3 <- merge(deGenes2, agg_bval, by.x = "target_id", by.y = "ensembl_gene", all = FALSE, sort = FALSE)
    sigDEG2 <- deGenes3[deGenes3$qval < qval & abs(deGenes3$b) > bval, ]
    # Saving deGenes3
    deGenes3Name <- paste(compStr, "_allDEG_bvalue", sep="")
    result[[deGenes3Name]] <- deGenes3
  } 
  else 
  { ## LRT test
    sigDET <- deTrans2[deTrans2$qval < qval, ]
    sigDETp <- deTrans2[deTrans2$pval < qval, ]
    sigDEG2 <- deGenes2[deGenes2$qval < qval, ]
  }
  
  ## Common steps regardless of sleuth test type
  sigDetName <- paste("sigDET_", compStr, sep="")
  sigDegName <- paste("sigDEG_", compStr, sep="")
  result[[sigDetName]] <- sigDET
  result[[sigDegName]] <- sigDEG2
  numDET <- dim(sigDET)[1]
  numDEG <- dim(sigDEG2)[1]
  # Add (union of) all sigDET to all_sigDETs dataframe
  # used when performing union of PCA plot for "n" way (example n = 4 or 6)  comparisons
  all_sigDETs <- rbind(all_sigDETs, sigDET) 
  all_sigDEGs <- rbind(all_sigDEGs, sigDEG2)
  
  # Defining names for DEG analysis downstream
  signifFile <- file.path(outPath, paste("sig_DE_", compStr, ".xlsx", sep=""))

  ## Separate steps depending on sleuth test
  if(sleuth_test =="wt")
  {
    sigDET_UD <- processRegulation(sigDET, sigDetName,result); result <- sigDET_UD$result
    sigDEG_UD <- processRegulation(sigDEG2, sigDegName,result); result <- sigDEG_UD$result
    write.xlsx(x=list(sigDEG=sigDEG2, sigDEG_Up=sigDEG_UD$up, sigDEG_Dn=sigDEG_UD$down, 
                      sigDET=sigDET[,1:7], sigDET_Up=sigDET_UD$up[,1:7], sigDET_Dn=sigDET_UD$down[,1:7]), file=paste0(signifFile))
    
    # Add a row to summary dataframe
    numUpG <- dim(sigDEG_UD$up)[1]
    numDnG <- dim(sigDEG_UD$down)[1]
    numUpT <- dim(sigDET_UD$up)[1]
    numDnT <- dim(sigDET_UD$down)[1]
    newRow <- data.frame(Comparisons=compStr, DE_Genes=numDEG, Up_Reg_Genes=numUpG, Down_Reg_Genes=numDnG,
                         DE_Transcripts=numDET, Up_Reg_Transcripts=numUpT, Down_Reg_Transcripts=numDnT)
  } 
  else 
  {  ## LRT test
    write.xlsx(x=list(sigDET=sigDET, sigDEG=sigDEG2), file=paste0(signifFile))  
    # Add a row to summary dataframe
    newRow <- data.frame(Comparison=compStr, DE_Genes=numDEG, DE_Transcr=numDET)
  }

  ## Common steps regardless of sleuth test type
  # For each test pair adding the summary row to the final DF
  summaryDF <- rbind(summaryDF, newRow)
} # completes one iteration of DEG analysis for each factor-level pair

### Running DEG analysis : Saving results for all possible combinations ----
result[["summaryDF"]] <- summaryDF
result[["all_sigDETs"]] <- all_sigDEGs
result[["all_sigDETs"]] <- all_sigDETs
save(result, file = file.path(outPath, "sleuth_results.RData"))

### Plotting DEG analysis : Generating and Saving Plots (PCA/Volcano) ----
for(i2 in 1:length(testsToRun))
{
  ## Common steps regardless of sleuth test type
  # Variables defining the writing path of files
  t1 <- testsToRun[i2]
  c1 <- gsub(factors, "", t1)
  if(!is.null(factorOfInt)){
    c1 <- gsub(factorOfInt, "", t1)
  }
  fr1 <- gsub(c1, "", t1)
  r1 <- refLevel[factors == fr1]
  compStr <- paste(c1, "_vs_", r1,sep="")
  
  ### PCA plotting and Saving ----
  sigDET <- result[[paste("sigDET_", compStr, sep="")]]
  # Create subCounts for sigDET
  subCountsResults <- c("subCountsGenes")
  if(nrow(sigDET) > 10)
  {
    subCountsTrans <- normCounts[sigDET$target_id, ]
    subCountsResults <- c("subCountsTrans", subCountsResults)
  }
  # Create subCounts for sigDEG
  subCountsGenes <- result[["normCounts_names"]][result[["normCounts_names"]]$gene_name %in% sigDEG2$gene_name, ]
  subCountsGenes <- as.matrix(subCountsGenes[,-1])
  subCountsGenes <- subCountsGenes[rowSums(subCountsGenes) > 0,]
  # Create plot
  for(s1 in 1:length(subCountsResults))
  {
    subCounts <- get(subCountsResults[s1])
    pcaName <- sub("subCounts", "", subCountsResults[s1])

    if(nrow(sigDET) > 10)
    {
      inPCA <- as.data.frame(t(subCounts))
      inPCA <- inPCA[sampleDF$sample, ]   #put inPCA in same order as sampleDF
      keepSamples <- sampleDF[sampleDF[,fr1] %in% c(r1, c1), "sample"]
      inPCA <- inPCA[keepSamples,]    #keep only the samples corresp. to the reference & comparison levels
      inPCA <- inPCA[ , !apply(inPCA, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))] #Exclude columns that have constant/zero values
      inPCA <- as.data.frame(cbind(inPCA, sampleDF[sampleDF[,fr1] %in% c(c1, r1), fr1]))
      colnames(inPCA)[colnames(inPCA) == "sampleDF[sampleDF[, fr1] %in% c(c1, r1), fr1]"] <- fr1
      inPCA[ ,fr1] <- as.factor(inPCA[ ,fr1])
      samplePCA <- prcomp(inPCA[ , 1:ncol(inPCA)-1], center = TRUE, scale. = TRUE)
      pca_plot <- autoplot(samplePCA, data = inPCA, colour = fr1, size = 3) + theme(legend.position = "right", axis.text=element_text(size=12)) + geom_text(aes(label = rownames(inPCA)), nudge_y = 0.1)
      # Add PCA plot to result
      pca_name <- paste(compStr, "_pca_", pcaName, sep="")
      write.xlsx(pca_plot$data[, c("PC1", "PC2", "Groups")], file.path(outPath, paste0(pca_name, ".xlsx")))
      ggsave(filename=file.path(outPath, paste0(pca_name, ".pdf")),width = 10, height = 8, dpi = 300,plot = pca_plot)
      result[[pca_name]] <- print(pca_plot)
    }
  }
  ### Volcano Plotting and Saving ----
  # defining variables we need volcano plot for
  deTrans2Name <- paste(compStr, "_allDET", sep="")
  deGenes3Name <- paste(compStr, "_allDEG_bvalue", sep="")
  deResult <- c(deTrans2Name, deGenes3Name)
  for(d1 in deResult)
  {
    if (is.null(geneList)) 
    {
      selectLab <- NULL
    } 
    else 
    if (is.null(geneList[[t1]])) 
    {
      selectLab <- NULL
    } 
    else 
    {
      selectLab <- geneList[[t1]]
    }
    # Create plot
    volc_plot <- EnhancedVolcano::EnhancedVolcano(result[[d1]],lab=result[[d1]]$gene_name, x='b', y='qval',
                                                  pCutoff = qval, pointSize = 2.0,  
                                                  xlab = bquote(~Log[2] ~ "fold change"), axisLabSize = 11.0,
                                                  labCol= 'black', labFace = 'bold', labSize = 3.0,  
                                                  boxedLabels = FALSE, drawConnectors = FALSE, 
                                                  widthConnectors = 1.0, colConnectors = 'black',
                                                  legendLabSize = 8.0, legendIconSize = 3.0,
                                                  legendPosition = 'right', title = paste(c1, "vs", r1),
                                                  selectLab = c("NonExistentLabel"),
                                                  col = c("grey30", "grey30", "grey30", "red2"))
                                                
    volc_plot + theme(axis.text.x = element_text(color="black", size=15, family = "Helvetica"), 
                      axis.text.y = element_text(color="black", size=15, family = "Helvetica"),
                      panel.grid.minor = element_blank(),panel.grid.major = element_blank())
    # Add volcano plot to result
    volcName <- paste(d1, "_vPlot", sep="")
    write.xlsx(volc_plot$data[, c("target_id", "gene_name", "xvals", "yvals")], file.path(outPath, paste0(volcName, ".xlsx")))
    ggsave(filename=file.path(outPath, paste0(volcName, ".pdf")),width = 10, height = 8, dpi = 300,plot = volc_plot)
    result[[volcName]] <- volc_plot
  }
} # completes one iteration of DEG analysis for each factor-level pair

### Plotting DEG analysis : Generating and Saving Plots (PCA All Samples (Allgenes and DEG)) ----
sleuth_results <-result
sampleInfo <- sleuth_results$sampleDF
unionDEGs <-rbind(sleuth_results$sigDEG_GLU_vs_UNT,sleuth_results$sigDEG_SYM_vs_UNT)
normCounts_names_All <- sleuth_results$normCounts_names
normCounts_names_DEG <- normCounts_names[normCounts_names$gene_name %in% unionDEGs$gene_name, ]
subCountsResults <- c("normCounts_names_All", "normCounts_names_DEG")
for(s1 in 1:length(subCountsResults))
{
  subCounts <- get(subCountsResults[s1])
  pcaName <- gsub("normCounts_names", "pca", subCountsResults[s1])
  subCounts <- as.matrix(subCounts[,-1])
  subCounts <- subCounts[rowSums(subCounts) > 0,]

  inPCA <- as.data.frame(t(subCounts))
  inPCA <- inPCA[sampleInfo$sample , ]
  inPCA["Groups"] <- sampleInfo["Groups"]

  samplePCA <- prcomp(inPCA[ , 1:(ncol(inPCA)-1)], center = TRUE, scale. = TRUE)
  pca_plot <- autoplot(samplePCA, data = inPCA, colour = "Groups", size = 3) + theme(legend.position = "right", axis.text=element_text(size=12)) 

  write.xlsx(pca_plot$data[, c("PC1", "PC2", "Groups")], file.path(outPath, paste0(pcaName, ".xlsx")))
  ggsave(filename=file.path(outPath, paste0(pcaName, ".pdf")),width = 10, height = 8, dpi = 300,plot = pca_plot)
}

### Fethcing gene list for ORA ----
if (runReactome | runKEGG)
{
  cp_path <- file.path(outPath, "pathway_analysis_ORA")
  dir.create(cp_path, recursive=TRUE)
  ORA <- list()
  plots_list <- list()
  sigDE_UpOrDn <- grep("Up|Dn", names(sleuth_results), value=TRUE)
  sigDE_UpOrDn <- grep("^sigDET",sigDE_UpOrDn, value=TRUE) #change to sigDET or sigDEG if using transcript or genes to run ORA
  for(g1 in sigDE_UpOrDn)
  {
    # Get gene names and convert to entrez id
    genes <- sleuth_results[[g1]]$gene_name
    entrez_genes <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
    if(runReactome)
    {
      # Reactome pathway over-representation analysis
      reactome <- enrichPathway(gene = entrez_genes$ENTREZID, organism = "mouse", pvalueCutoff=pvalueCutoff, readable=TRUE)
      ORA[[paste0(g1, "_Reactome")]] <- reactome
    }
    else
    if(runKEGG)
    {
      # KEGG pathway over-representation analysis
      kegg <- enrichKEGG(gene = entrez_genes$ENTREZID, organism = "mmu",  pvalueCutoff=pvalueCutoff, use_internal_data = FALSE)
      kegg2 <- setReadable(kegg, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
      kegg2@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", kegg2@result$Description, fixed = T)    
      ORA[[paste0(g1, "_KEGG")]] <- kegg2
    }
    else
    {
      stop("Only Reactome or KEGG are supported.")
    }
  }
}

### Performing Over Representation Analysis (ORA) ----
for(a1 in names(ORA))
{
  reacORkegg <- ORA[[a1]]
  if(nrow(reacORkegg) !=0)
  {
    # Create directories for each result
    ORA_path <- file.path(cp_path, a1)
    dir.create(ORA_path, recursive=TRUE)
    ORA_df <- as.data.frame(reacORkegg)
    # Convert 'GeneRatio' and 'BgRatio' from character to numeric
    ORA_df$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", ORA_df$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", ORA_df$GeneRatio, perl=T) )
    ORA_df$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", ORA_df$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", ORA_df$BgRatio, perl=T) )
    # Now perform the division to create 'RatioOfRatios'
    ORA_df <- transform(ORA_df, EnrichmentRatio = GeneRatio / BgRatio)
    write.xlsx(ORA_df, file.path(ORA_path, paste0(a1, ".xlsx")))
    
    # Plot enrichment result - dotplot
    fit <- dotplot(reacORkegg, showCategory = 40)
    ggsave(file.path(ORA_path, paste0(a1, "_dotplot.pdf")), plot = fit, dpi = 300, width = 15, height = 12, units = "in", device="pdf")
    plots_list[[a1]] <- fit
  }

### Collapsing into Weighted Set Cover (ORA) ----  
  if (nrow(reacORkegg) > 40)
  {
    # Weighted Set Cover of geneSets
    weightedPath <- file.path(cp_path, "weightedSetCover", a1)
    dir.create(weightedPath, recursive=TRUE)
    
    setCoverNum = abs(0.40*(nrow(ORA_df)))
    nThreads = 4
    idsInSet <- sapply(ORA_df$geneID, strsplit, split="/")
    names(idsInSet) <- ORA_df$ID    
    minusLogP <- -log(ORA_df$pvalue)
    minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
    wscRes <- weightedSetCover(idsInSet=idsInSet, costs=(1 / minusLogP), topN=setCoverNum, nThreads=nThreads)
    
    wscRes_full <- ORA_df[c(match(wscRes$topSets, ORA_df$ID)),]
    wscRes_full <- wscRes_full[order(wscRes_full$p.adjust), ]
    write.xlsx(wscRes_full, file.path(weightedPath, paste0(a1, ".xlsx")))
    
    # Plot enrichment result of weighted set cover - dotplot
    reacORkegg2 <- reacORkegg
    reacORkegg2@result <- wscRes_full
    
    fit <- dotplot(reacORkegg2, showCategory = 20)
    ggsave(file.path(weightedPath, paste0(a1, "_dotplot.pdf")), plot = fit, dpi = 300, width = 17, height = 12, units = "in", device="pdf")
    plots_list[[paste0(a1, "weightSet")]] <- fit
  }
}

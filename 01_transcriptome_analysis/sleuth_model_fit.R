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
library(yaml)

### Load the configuration from the YAML file ----
args <- commandArgs(trailingOnly = TRUE)
configFilePath <- if(length(args) > 0) args[1] else stop("No YAML file path provided.")
config <- yaml::read_yaml(configFilePath)

# Basic file path setup
proj_path <- config$proj_path
dataDir <- config$dataDir
metadataPath <- config$metadataPath
current_date <- Sys.Date()
formatted_date <- format(current_date, "%d-%m-%Y")
experiment <- paste0(formatted_date, "-DE-Analysis")
experiment_path <- file.path(proj_path, experiment)
outPath <- experiment_path
dir.create(experiment_path, recursive=TRUE)

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

biomart <- config$biomart
dataset <- config$dataset

### Sleuth : Loading annotation and Setting up variable list ----
# Loading genome and genomic annotation
mart <- useMart(biomart, dataset=dataset) 
mappingDF <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"), mart = mart) 
colnames(mappingDF) <- c("ensembl_gene", "target_id", "gene_name")
aggCol <- "ensembl_gene"

### Running Sleuth : Preliminary variable setting based on inputs ----
# QC of inputs
if(length(factors) > 1 && is.null(factorOfInt))
{
  stop("factorOfInt cannot be null for multifactor analysis")
}
if(length(factors) == 1 && is.null(factorOfInt))
{
  factorOfInt == factors
}

# Return results
result <- list()

# Add kallisto file paths to sample info dataFrame 
sampleDF <- read.csv(metaData, skip = skipLines)
colnames(sampleDF)[names(sampleDF) == sampleCol] <- "sample"

# Converting factor columns to factor class
sampleDF[factors] <- lapply(sampleDF[factors], factor)  ## as.factor() could also be used
sampleDF <- sampleDF[ ,c("sample", factors)]            ## keep only "sample", factor columns

# Re-level factors 
for (i in 1:length(factors))
{  #loop thru factors by index
  fac <- factors[[i]]
  lev <- refLevel[[i]]
  sampleDF[,fac] <- relevel(sampleDF[,fac], ref = lev)
}

### Running Sleuth : inferring kallisto abundance files paths ----
# Getting abundance.h5 paths and adding path column to sampleDF
abundancePaths <- c()
for(a1 in dataDir)
{
  aPaths <- list.dirs(path = a1, recursive = FALSE)
  kallistoDirNames <- list.dirs(path = a1, full.names = FALSE, recursive = FALSE) ## gets directory names
  names(aPaths) <- kallistoDirNames
  if(any(kallistoDirNames %in% sampleDF$sample))
  {
    names(aPaths) <- kallistoDirNames
  }
  else
  {
    isPresent <- rowSums(sapply(sampleDF$sample, grepl, kallistoDirNames))
    names(isPresent) <- kallistoDirNames
    if(sum(isPresent) > 0)
    {
      keepKallistoDirNames <- names(isPresent)[isPresent == 1]
      aPaths <- aPaths[names(aPaths) %in% keepKallistoDirNames]
      rowOrder <- max.col(sapply(sampleDF$sample, grepl, keepKallistoDirNames))
      names(aPaths) <- sampleDF$sample[rowOrder]
    }
  }
  abundancePaths <- c(abundancePaths, aPaths)
}

sampleDF$path <- abundancePaths[match(sampleDF$sample, names(abundancePaths))]
result[["sampleDF"]] <- sampleDF

### Running Sleuth : Creating 'prepped' sleuth object ----
# Create full model formula out of list of factors
fmla <- as.formula(paste0("~", paste0(factors, collapse=" + ")))
result[["design"]] <- fmla

# Create reduced formula out of list of factors
if(!is.null(factorOfInt))
{
  for(f1 in factorOfInt)
  {
    redFmla <- as.formula(paste0("~", f1))
    result[["design"]] <- c(result[["design"]], redFmla)
  }
}

if(!is.null(filter_target_id))
{
  sPrep <- sleuth_prep(sampleDF, read_bootstraps = FALSE, normalize = FALSE)
  filter_sleuth <- sPrep$filter_df$target_id
  filter_target_id <- filter_target_id[which(filter_target_id %in% filter_sleuth)]
}

# Create 'prepped' sleuth object
so <- sleuth_prep(sampleDF, full_model = fmla, target_mapping = mappingDF, 
                  read_bootstrap_tpm = TRUE, filter_target_id = filter_target_id,
                  extra_bootstrap_summary = TRUE, num_cores = 1,
                  aggregation_column = aggCol, 
                  transformation_function = function(x) log2(x + 1))

### Running Sleuth : Saving TPM and Counts ----
namesDF <- so$target_mapping[,c('gene_name','target_id')]
matrix_save <-function(df,namesDF,savename)
{
  DF <- as.data.frame(df)
  DF_names <- merge(DF, namesDF, by.x=0, by.y='target_id', all.x = TRUE)
  rownames(DF_names) <- DF_names$Row.names
  DF_names <- subset(DF_names, select=-Row.names)
  DF_names <- DF_names[, c(which(colnames(DF_names)=="gene_name"), which(colnames(DF_names)!="gene_name"))]
  write.csv(DF_names,savename)
  return(DF_names)
}

Tpm <- sleuth_to_matrix(obj = so, which_df = "obs_raw", which_units = "tpm") #Get TPM counts
result[["TPM"]] = Tpm; result[["TPM_names"]] = matrix_save(Tpm, namesDF, file.path(outPath, "tpm_names.csv"));
normTpm <- sleuth_to_matrix(obj = so, which_df = "obs_norm", which_units = "tpm") # Get normalized TPM counts
normCounts <- sleuth_to_matrix(obj = so, which_df = "obs_norm", which_units = "est_counts")
write.csv(normCounts,file.path(outPath, "normalized_counts.csv"))
result[["normCounts"]] = normCounts; result[["normCounts_names"]] = matrix_save(normCounts, namesDF, file.path(outPath, "tpm_names.csv"));

### Running Sleuth : Fitting sleuth model for DEGs and Saving sleuth object ----
# Fit full model
so <- sleuth_fit(so) # fit full model (includes all terms from 'fmla')
if(!is.null(factorOfInt))
{ # Fit reduced model (only one factor)
  so <- sleuth_fit(so, redFmla, "reduced")
} 
else
{ # Fit intercept only model
  so <- sleuth_fit(obj = so, formula = ~1, fit_name = "reduced")
}

### Running Sleuth : Saving sleuth and result object ----
result[["sleuth_model_fit"]] <- so
result[["sleuth_models"]] <- models(so)
save(so,file=file.path(outPath,"sleuth_object.RData"))
save(result, file = file.path(outPath, "sleuth_results_noDEG.RData"))
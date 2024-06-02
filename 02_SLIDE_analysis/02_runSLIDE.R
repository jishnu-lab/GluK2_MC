#### Install SLIDE if not installed already ----
# library(devtools)
# devtools::install_github("jishnu-lab/SLIDE")

#### Running SLIDE ----
library(SLIDE)
library(yaml)
# library(rstudioapi)
require(doParallel)
require(dplyr)

# current_file_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(".")

set.seed(1234)
args <- commandArgs(trailingOnly = T)
yaml_path  <- args[1]
# yaml_path = "./slide.yaml"
input_params <- yaml::yaml.load_file(yaml_path)

SLIDE::checkDataParams(input_params)
SLIDE::optimizeSLIDE(input_params, sink_file = TRUE)
SLIDE::plotCorrelationNetworks(input_params)
SLIDE::SLIDEcv(yaml_path, nrep = input_params$CViter, k = input_params$sampleCV_K)

# If you want to save the boxplot data as csv file
write.csv(readRDS("SLIDECV_boxplot_data.rds"), "SLIDECV_boxplot.csv")

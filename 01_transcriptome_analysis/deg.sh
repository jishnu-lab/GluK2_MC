#!/bin/bash
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH --job-name=sleuth_main
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=swk25@pitt.edu

module purge
module load gcc/12.2.0
module load r/4.3.0
# Requires path to yaml as input as argument 1
Rscript sleuth_deg_ora.R $1
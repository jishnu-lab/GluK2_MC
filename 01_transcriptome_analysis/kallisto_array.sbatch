#!/bin/bash
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH --job-name=kallisto_array
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=swk25@pitt.edu
#SBATCH --array=0-19 # job array index. Please, change the max number. number of samples - 1

## Loading kallisto module
module purge
module load kallisto/0.46.1

DATAPATH=/ix/djishnu/Swapnil/kaplanAnalysis/kaplanData/data/YZ_Files/Youran_1_05_2023 ## Change the datapath
OUTPATH=/ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/mapped ## Change result path
mkdir -p $OUTPATH

names=($(cat /ix/djishnu/Swapnil/kaplanAnalysis/final_analysis/YZ_list_fastq.txt)) ## Read the text file with sample names
echo ${names[${SLURM_ARRAY_TASK_ID}]}

kallisto quant -i /ix/cigcore/utils/mouse_indexed/mm10_kallisto.idx \
               -b 100 -o $OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]} \
	             $DATAPATH/${names[${SLURM_ARRAY_TASK_ID}]}_R1_001.fastq.gz \
	             $DATAPATH/${names[${SLURM_ARRAY_TASK_ID}]}_R2_001.fastq.gz

crc-job-stats
#!/bin/bash

## Credit: The associated R script is adapted from: https://benjjneb.github.io/dada2/bigdata.html
## This script is to run the main script "02.Denoise.R" from a HPC cluster environemnt. 
## See the R script "02.Denoise.R" for details on script functions

printf "\nStarting at: %s\n" "$(date)"

## Activate the R environment
conda activate R

## Run R script in a clean R instance, output a logfile
R_SCRIPT_PATH="../bioinformatics/02.Denoise.R"
LOG_FILE="../bioinformatics/slurm/02.Denoise.sh.${SLURM_JOBID}.Rout"
OUTPUT_LOG="../bioinformatics/slurm/02.Denoise.sh.${SLURM_JOBID}.out"

Rscript --vanilla --verbose "$R_SCRIPT_PATH" > "$LOG_FILE" 2>&1 

## Append Rout log to the slurm out log
cat "$LOG_FILE" >> "$OUTPUT_LOG"

## Remove Rout log
rm "$LOG_FILE"

## Deactivate the R environment
conda deactivate

printf "\nFinished at: %s\n" "$(date)"

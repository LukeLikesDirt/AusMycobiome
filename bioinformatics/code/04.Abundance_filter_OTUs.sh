#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --cpus-per-task=40
#SBATCH --output=slurm/%x.%j.out

# Script:   Remove low abundance OTUs and samples 
# Author:   Luke Florence.
# Date:     28th Janurary 2024.

# Notes:
# --------------------------------------------------------------------------- #
#   -   See the associated R scripts for details on script purose.
#   -   File paths are handeled in the R scripts - this aspect needs updating.
#   -   The 'ASVs' outputs corresponds to taxonomically informed dynamic 
#       clustering developed for this project. The 'OTUs' output is used to 
#       compare results between conventional 97% OTU clustering and the dynamic
#       clustering.
#   -   The R code with 'thresholds' suffixes is for compareing results using
#       a range of filtering threcholds. 

echo "Starting at: $(date)"

source ~/.bashrc
conda activate dynamic_clustering

Rscript --vanilla --verbose abundance_filter_ASVs_thresholds.R
Rscript --vanilla --verbose abundance_filter_OTUs_thresholds.R
Rscript --vanilla --verbose abundance_filter_ASVs.R
Rscript --vanilla --verbose abundance_filter_OTUs.R

conda deactivate

echo "Finished at: $(date)"
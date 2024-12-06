#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --cpus-per-task=40
#SBATCH --output=slurm/%x.%j.out

# Script:   Classify OTUs.
# Author:   Luke Florence.
# Date:     18th November 2024.

# Notes:
# --------------------------------------------------------------------------- #
# - Similarity cutoffs were predicted using DNABarcoder.
# - Coverage cutoffs are applied at 90% for genus-level and 95% for 
#   species-level identifications.
# - Consensus cutoffs are arbitrarily set at 66% to ensure some degree of
#   precision when accepting best hit at a given taxanomic rank.
# - The 'ASVs' outputs corresponds to taxonomically informed dynamic 
#   clustering developed for this project. The 'OTUs' output is used to 
#   compare results between conventional 97% OTU clustering and the dynamic
#   clustering.

echo "Starting at: $(date)"

source ~/.bashrc
conda activate dynamic_clustering

Rscript --vanilla --verbose classify_ASVs.R
Rscript --vanilla --verbose classify_OTUs.R

conda deactivate

echo "Finished at: $(date)"
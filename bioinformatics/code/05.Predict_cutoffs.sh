#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=slurm/%x.%j.out

# Script: Predict taxonomically guided cutoffs using the dnabarcoder framework
# Credit: Adapted from https://github.com/vuthuyduong/dnabarcoder
# Author: Luke Florence
# Date: 21st October 2024
#
# Software:
# --------------------------------------------------------------------------- #
#  - BLAST v2.14.1 - https://blast.ncbi.nlm.nih.gov/Blast.cgi
#  - dnabarcoder v1.0.6 - https://github.com/vuthuyduong/dnabarcoder
#
# Notes:
# --------------------------------------------------------------------------- #
#  - Run in parallel to speed up predictions
#  - Should take about a week with the above slurm specifications
#  - 05a.Predict_cutoffs_dikarya_fungi.sh is the pinch point, taking 4-5 days

# Activate the conda environment
source ~/.bashrc
conda activate dynamic_clustering

echo "Starting at:" $(date)

# Predict cutoff for dikaryotic fungi
bash group_cutoff_predictions/05a.Predict_cutoffs_dikarya_fungi.sh
echo "05a.Predict_cutoffs_dikarya_fungi.sh completed at:" $(date)

# Predict cutoffs for (mostly) terrestrial basal fungi
bash group_cutoff_predictions/05b.Predict_cutoffs_terrestrial_basal_fungi.sh
echo "05b.Predict_cutoffs_terrestrial_basal_fungi.sh completed at:" $(date)

# Predict cutoffs for (mostly) single-celled zoosporic fungi
bash group_cutoff_predictions/05c.Predict_cutoffs_zoosporic_basal_fungi.sh
echo "05c.Predict_cutoffs_zoosporic_basal_fungi.sh completed at:" $(date)

# Predict global cutoffs for all fungi and the above groups individually
bash group_cutoff_predictions/05d.Predict_cutoffs_global_fungi.sh
echo "05d.Predict_cutoffs_global_fungi.sh completed at:" $(date)

# Merge files ready for taxonomic assignment and OTU clustering
bash group_cutoff_predictions/05e.Predict_cutoffs_merge_fungi.sh
echo "05e.Predict_cutoffs_merge_fungi.sh completed at:" $(date)

# Deactivate the conda environment
conda deactivate

echo "Finished at:" $(date)
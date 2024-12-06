#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --cpus-per-task=60
#SBATCH --output=slurm/%x.%j.out

# Script:   Denoise Illumina single-end reads.
# Purpose:  Prepare Illumina forwards reads for chimera detection in VSEARCH.
# Credit:   The associated R script is adapted from: https://benjjneb.github.io/dada2/bigdata.html
# Author:   Luke Florence.
# Date:     28th Janurary 2024.

echo "Starting at: $(date)"

source ~/.bashrc
conda activate dynamic_clustering

Rscript --vanilla --verbose denoise_dada2.R

conda deactivate

echo "Finished at: $(date)"
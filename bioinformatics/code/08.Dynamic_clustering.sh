#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=60
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=slurm/%x.%j.out

# Script: Cluster OTUs (i.e. species level cluster cores) and pseudo taxa across ranks using taxanomically informed cutoffs
# Credit: This approch has been adapted from a similar approch using other tools for cutoff predictions and clustering: see https://zenodo.org/records/11125610, https://github.com/brendanf/deadwood_restoration, https://github.com/brendanf/deadwood_priority_effects, and https://github.com/brendanf/bistorta_vivipara_translocation
# Author: Luke Florence
# Date:   24th November 2024
#
# Note on couverage thresholds for clustering:
#  - Requires >90% coverage on at least one sequence for phylum to family
#  - Requires >90% coverage on both sequences for genus
#  - Requires >95% coverage on both sequences for species
#  - blastn (reference clustering): Uses custom code to implement these rules
#  - blastclust (de novo clustering): -L X Length coverage threshold (default = 0.9), switches to 0.95 for species
#  - blastclust (de novo clustering): -b F Do not require coverage on both neighbors, switches to T for genus and species
#
# Notes on need improvments:
#  - To speed up clusering, once an ASV partitions into a single query move to an "OTU bin" and megere once clustering is complete
#  - Couverage threshold are arbitrarily assigned and future efforts could estimate length variation

# Activate the conda environment
echo "Starting at:" $(date)
source ~/.bashrc
conda activate dynamic_clustering

# Define constants, input files and output directory
TAXA_FILE="../data/annotated_asvs/taxonomy.txt"
TAXA_CUTOFFS="./dnabarcoder/ITS1_cutoffs.txt"
ASV_SEQUENCES="../data/annotated_asvs/sequences.fasta"
ASV_TABLE="../data/annotated_asvs/asv_table.txt"
OUTPUT="../../output/"
THREADS=60
MINLEN=50

# Make output directoiry
mkdir -p $OUTPUT

# Dynamic clustering in R
Rscript --vanilla dynamic_clustering.R "$TAXA_FILE" $TAXA_CUTOFFS "$ASV_SEQUENCES" "$ASV_TABLE" "$OUTPUT" "$THREADS" "$MINLEN"

conda deactivate
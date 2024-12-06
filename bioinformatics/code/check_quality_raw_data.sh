#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Script:   Quality check raw reads
# Purpose:  Inform truncation length
# Author:   Luke Florence
# Date:     16th September 2024

# Software:
# --------------------------------------------------------------------------- #
# FastQC v0.12.1: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# MultiQC v1.15: https://multiqc.info/
#
# Note:
# --------------------------------------------------------------------------- #
# The raw reads are expected in subdirectories named 'run1', 'run2', 'run3', etc.

# Raw data directory
RAW_DATA="../data/raw_data"

# Constants
readonly NUM_THREADS=40           # The number of threads to use for parallel processing
readonly FIRST_RUN=1              # The first run to process
readonly LAST_RUN=42              # The last run to process

# Function to generate read quality report for each sample using FastQC and MultiQC
generate_quality_report() {
    local run_number="$1"

    echo "Starting to generate the quality report for run $run_number"
    
    fastqc "${RAW_DATA}/run${run_number}"/*.fastq.gz -o "${RAW_DATA}/run${run_number}"
    multiqc "${RAW_DATA}/run${run_number}" -o "${RAW_DATA}/run${run_number}"
    
    # Remove intermediate files and 'multiqc_data' directory
    rm "${RAW_DATA}/run${run_number}"/*fastqc.zip "${RAW_DATA}/run${run_number}"/*fastqc.html
    rm -r "${RAW_DATA}/run${run_number}"/multiqc_data

    echo "Completed generating a quality report for run $run_number"
}

###############################################################################
### Main script ###############################################################
###############################################################################

echo "Starting at:" $(date)

### Exicute functions #########################################################

# Activate the conda environment
source ~/.bashrc
conda activate sequence_prep

# Process the fastq files
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    generate_quality_report "$run_number"
done

# Deactivate the conda environment
conda deactivate
echo "Finished at:" $(date)
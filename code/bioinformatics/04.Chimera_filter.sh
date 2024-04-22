#!/bin/bash

## Author: Luke Florence
## Date:   4th February 2024
## Credit: This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline.

# Software:
# --------------------------------------------------------------------------- #
# VSEARCH v2.22.1: https://github.com/torognes/vsearch
# Perl v5.32.1: https://www.perl.org/get.html

# Constants and subdirectories
readonly THREADS=8
readonly IDENTITY=0.97
readonly MAP_SCRIPT="../bioinformatics/map.pl" # see the following link for an example pipeline that does not require a mapping file: https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline
readonly REFERENCE_SEQS="../../data/bioinformatics/06.Reference_dataset/general_release/ITS1/ITS1.fasta"
CHIMERA_FILTERED_DIR=../../data/bioinformatics/07.Chimera_filtered
TRACK_REPSEQS_READS_DIR="../../data/bioinformatics/"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp" | tee -a "$LOG_FILE"
}
LOG_FILE="slurm/%x.%j.out"

## Function for dereplication across samples, de novo and referenced-based chimera filtering
chimera_filter() {
    local files_to_remove=("all.denovo.nonchimeras.fasta" "all.derep.fasta" "all.derep.uc" "all.nonchimeras.derep.fasta" "all.nonchimeras.fasta" "all.preclustered.fasta" "all.preclustered.uc" "all.ref.nonchimeras.fasta")

    for file in "${files_to_remove[@]}"; do
        rm -f "$CHIMERA_FILTERED_DIR/$file"
    done

    log 'Dereplicating across samples at:'
    vsearch \
        --derep_fulllength "$CHIMERA_FILTERED_DIR/all.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.derep.uc" \
        --output "$CHIMERA_FILTERED_DIR/all.derep.fasta"

    log 'Preclustering reads at:'
    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.derep.fasta" \
        --threads $THREADS \
        --id $IDENTITY \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.preclustered.uc" \
        --centroids "$CHIMERA_FILTERED_DIR/all.preclustered.fasta"

    log 'Starting de novo chimera detection at:'
    vsearch \
        --uchime3_denovo "$CHIMERA_FILTERED_DIR/all.preclustered.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta"

    log 'Starting reference-based chimera detection at:'
    vsearch \
        --uchime_ref "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta" \
        --threads "$THREADS" \
        --db "$REFERENCE_SEQS" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta"

    log 'Extracting all non-chimeric sequences at:'
    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.derep.fasta" "$CHIMERA_FILTERED_DIR/all.preclustered.uc" "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta"
    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.fasta" "$CHIMERA_FILTERED_DIR/all.derep.uc" "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta"
}

###############################################################################
### Main script ###############################################################
###############################################################################

log 'Starting at:'

# Activate the conda environment
conda activate shell

# Dereplicate across samples and remove chimeras in the background
chimera_filter

# Deactivate the conda environment
conda deactivate

log 'Finishing at:'
#!/bin/bash

## Script: Cluster OTUs
## Author: Luke Florence
## Date: 5th February 2024
## Software: VSEARCH v2.22.1: https://github.com/torognes/vsearch

## Constants and file paths
readonly THREADS=8                                                   ## Set the number of threads
readonly IDENTITY=0.97                                               ## Set the identity threshold for clustering
CHIMERA_FILTERED_DIR=../../data/bioinformatics/07.Chimera_filtered   ## Path to chimera filtered fasta file
CLUSTERED_OTUs_DIR="../../data/bioinformatics/08.OTUs"               ## Path to clustered fasta file and OTU table

# Define extensions for different methods used for denoising

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## Function for clustering OTUs and formatting the OTU table
cluster_OTUs() {

    # Create the directory if it doesn't exist
    mkdir -p "$CLUSTERED_OTUs_DIR"

    log "Clustering at 97% at:"

    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta" \
        --threads "$THREADS" \
        --id "$IDENTITY" \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --relabel_sha \
        --uc "$CHIMERA_FILTERED_DIR/OTUs.uc" \
        --centroids "$CLUSTERED_OTUs_DIR/OTUs.fasta" \
        --biomout "$CLUSTERED_OTUs_DIR/OTUs.biom" \
        --otutabout "$CLUSTERED_OTUs_DIR/OTUs.txt"

    # Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$CLUSTERED_OTUs_DIR/OTUs.txt"

    printf '\nNumber of unique sequences and OTUs\n'
    printf '    Unique non-chimeric sequence: %s\n' "$(grep -c "^>" "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta")"
    printf '    Clustered OTUs: %s\n' "$(grep -c "^>" "$CLUSTERED_OTUs_DIR/OTUs.fasta")"
}

# Function to build sequence matching table for post-clustering curation with LULU
match_list() {

    log "Building sequence matching table for at:"

    vsearch \
        --usearch_global "$CLUSTERED_OTUs_DIR/OTUs.fasta" \
        --db "$CLUSTERED_OTUs_DIR/OTUs.fasta" \
        --self --id .84 \
        --iddef 1 --userout "$CLUSTERED_OTUs_DIR/match_list.txt" \
        -userfields query+target+id \
        --maxaccepts 0 --query_cov 0.9 --maxhits 10
}

###############################################################################
## Main script ################################################################
###############################################################################

log 'Starting at:'

# Activate the conda environment
conda activate shell

# Run cluster_OTUs function
cluster_OTUs

# Run match_list function
match_list

# Deactivate the conda environment
conda deactivate

log 'Finished at:'
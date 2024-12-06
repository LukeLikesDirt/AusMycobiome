#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --output=slurm/%x.%j.out

# Author: Luke Florence
# Date:   4th February 2024
# Credit: This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline.

# Notes:
# --------------------------------------------------------------------------- #
# The 'OTUs_97' output is used to compare results between conventional 97% OTU 
# clustering and the contemporary taxonomically informed dynamic clustering 
# developed for this project. The 'ASVs' output corresponds to the output used
# for dynamic clustering. The 'match_list' output is designed for use with
# LULU (https://github.com/tobiasgf/lulu), an R package for distribution-based,
# post-clustering curation of amplicon data.

# Software:
# --------------------------------------------------------------------------- #
# VSEARCH v2.22.1: https://github.com/torognes/vsearch
# Perl v5.32.1: https://www.perl.org/get.html

# Constants and subdirectories
readonly THREADS=40
readonly IDENTITY=0.97
readonly MAP_SCRIPT="map.pl" # see the following link for an example pipeline that does not require a mapping file: https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline
readonly REFERENCE_SEQS="../data/reference_data/ITS.fasta" # Full ITS because amplicon tageted full ITS using ITS1F/ITS4
readonly CHIMERA_FILTERED_DIR="../data/chimera_filtered"
readonly CLUSTERED_DIR="../data/OTUs"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

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

## Function for generating an ASV table
generate_ASVs() {

    # Create the directory if it doesn't exist
    mkdir -p "$CLUSTERED_DIR"

    log "Generating ASVs at:"

    vsearch \
        --cluster_unoise "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta" \
        --threads $THREADS \
        --sizein --sizeout \
        --relabel_sha \
        --uc "$CLUSTERED_DIR/ASVs.uc" \
        --centroids "$CLUSTERED_DIR/ASVs.fasta" \
        --biomout "$CLUSTERED_DIR/ASVs.biom" \
        --otutabout "$CLUSTERED_DIR/ASVs.txt"

    # Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$CLUSTERED_DIR/ASVs.txt"

    printf '\nNumber of unique sequences and ASVs\n'
    printf '    Unique non-chimeric sequence: %s\n' "$(grep -c "^>" "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta")"
    printf '    Clustered ASVs: %s\n' "$(grep -c "^>" "$CLUSTERED_DIR/ASVs.fasta")"
}

generate_OTUs() {

    # Create the directory if it doesn't exist
    mkdir -p "$OTU_DIR"

    log "Generating OTUs at:"

    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta" \
        --threads $THREADS \
        --sizein --sizeout \
        --id $IDENTITY \
        --strand plus \
        --relabel_sha \
        --uc "$CLUSTERED_DIR/OTUs_97.uc" \
        --centroids "$CLUSTERED_DIR/OTUs_97.fasta" \
        --biomout "$CLUSTERED_DIR/OTUs_97.biom" \
        --otutabout "$CLUSTERED_DIR/OTUs_97.txt"

    # Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$CLUSTERED_DIR/ASVs.txt"

    printf '\nNumber of unique sequences and ASVs\n'
    printf '    Unique non-chimeric sequence: %s\n' "$(grep -c "^>" "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta")"
    printf '    Clustered OTUs at 97%: %s\n' "$(grep -c "^>" "$CLUSTERED_DIR/OTUs_97.fasta")"
}

# Function to build sequence matching table for post-clustering curation with LULU
match_list() {

    log "Building sequence matching table at:"

    vsearch \
        --usearch_global "$CLUSTERED_DIR/ASVs.fasta" \
        --db "$CLUSTERED_DIR/ASVs.fasta" \
        --self --id .84 \
        --iddef 1 --userout "$CLUSTERED_DIR/match_list.txt" \
        -userfields query+target+id \
        --maxaccepts 0 --query_cov 0.9 --maxhits 10
}

###############################################################################
### Main script ###############################################################
###############################################################################

log 'Starting at:'

# Activate the conda environment
source ~/.bashrc
conda activate dynamic_clustering

# Dereplicate across samples and remove chimeras in the background
chimera_filter

# Run cluster_OTUs for both methods concurrently
generate_ASVs
generate_OTUs

# Run match_list for both methods concurrently
match_list

# Deactivate the conda environment
conda deactivate

log 'Finishing at:'
#!/bin/bash

## Script: Taxonomic assignment with BLASTn
## Author: Luke Florence
## Date: 5th November 2023
## Software: BLAST v2.14.1 - https://blast.ncbi.nlm.nih.gov/Blast.cgi

# Constants and subdirectories
readonly THREADS=8
readonly REFERENCE_SEQUENCES="../../data/bioinformatics/06.Reference_dataset/ITS1/ITS1"
CLUSTERED_DIR="../../data/bioinformatics/08.OTUs"
TAXA_DIR="../../data/bioinformatics/09.Taxonomy"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## BLAST function

blast_best_ten_hits() {
    local OTU_FASTA="$CLUSTERED_DIR/OTUs.fasta"
    local TAXA_DIR="$TAXA_DIR"

    # Create the directory if it doesn't exist
    mkdir -p "$TAXA_DIR"

    # blast ten hits
    blastn \
         -task blastn \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -out "$TAXA_DIR/BLAST_best_10.txt" \
        -num_threads "$THREADS"

    
    sed -i '1s/^/OTU_ID;abundance\treference;kingdom;phylum;class;order;family;genus;species\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/' "$TAXA_DIR/BLAST_best_10.txt"
    sed 's/[[:space:]]\{1,\}/;/g' "$TAXA_DIR/BLAST_best_10.txt" > "$TAXA_DIR/BLAST_best_10.csv"
    rm "$TAXA_DIR/BLAST_best_10.txt"
        
}

###############################################################################
## Main script ################################################################
###############################################################################

log 'Starting at:'

conda activate shell

blast_best_ten_hits

conda deactivate

log 'Finished at:'
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=60
#SBATCH --partition=day
#SBATCH --time=1-00:00:00
#SBATCH --output="./slurm/%x.%j.out"

# Script: Taxonomic assignment with BLASTn
# Author: Luke Florence
# Date: 5th November 2024
# Software: BLAST v2.14.1 - https://blast.ncbi.nlm.nih.gov/Blast.cgi
#
# Notes:
# --------------------------------------------------------------------------- #
# - BLAST is performed on UNITE with species information and the full UNTE 
#   with prefernce given to the species only dataset when classifying OTUs.
# - The 'ASVs' outputs corresponds to taxonomically informed dynamic 
#   clustering developed for this project. The 'OTUs' output is used to 
#   compare results between conventional 97% OTU clustering and the dynamic
#   clustering.

# Constants and subdirectories
THREADS=60
REFERENCE_SEQUENCES_SPECIES="../data/reference_data/ITS1_species"
REFERENCE_SEQUENCES_GENUS="../data/reference_data/ITS1_genus"
REFERENCE_SEQUENCES_ALL="../data/reference_data/ITS1"
OTU_FASTA="../data/OTUs/OTUs_abundance_filtered.fasta"
ASV_FASTA="../data/OTUs/ASVs_abundance_filtered.fasta"
TAXA_DIR="../data/taxonomy"
mkdir -p $TAXA_DIR

blast_OTUs_species() {

    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES_SPECIES" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$TAXA_DIR/blast_output.txt" \
        -num_threads "$THREADS"

    # Define the header
    local header="OTU_ID\tread_abundance\treference_ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tspecies_hypothesis\tpident\tlength\tqlen\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

    # Process the file
    echo -e "$header" > "$TAXA_DIR/BLAST_OTUs_species.txt"
    sed -e 's/;size=/\t/g' -e 's/;/\t/g' -e 's/|/\t/g' \
        -e 's/k__//g' -e 's/p__//g' -e 's/c__//g' \
        -e 's/o__//g' -e 's/f__//g' -e 's/g__//g' -e 's/s__//g' \
        "$TAXA_DIR/blast_output.txt" >> "$TAXA_DIR/BLAST_OTUs_species.txt"
    
    # Remove temporary file
    rm -f "$TAXA_DIR/blast_output.txt"
        
}

blast_OTUs_genus() {

    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES_GENUS" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$TAXA_DIR/blast_output.txt" \
        -num_threads "$THREADS"

    # Define the header
    local header="OTU_ID\tread_abundance\treference_ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tspecies_hypothesis\tpident\tlength\tqlen\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

    # Process the file
    echo -e "$header" > "$TAXA_DIR/BLAST_OTUs_genus.txt"
    sed -e 's/;size=/\t/g' -e 's/;/\t/g' -e 's/|/\t/g' \
        -e 's/k__//g' -e 's/p__//g' -e 's/c__//g' \
        -e 's/o__//g' -e 's/f__//g' -e 's/g__//g' -e 's/s__//g' \
        "$TAXA_DIR/blast_output.txt" >> "$TAXA_DIR/BLAST_OTUs_genus.txt"
    
    # Remove temporary file
    rm -f "$TAXA_DIR/blast_output.txt"

}

blast_OTUs_all() {
    
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES_ALL" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$TAXA_DIR/blast_output.txt" \
        -num_threads "$THREADS"

    # Define the header
    local header="OTU_ID\tread_abundance\treference_ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tspecies_hypothesis\tpident\tlength\tqlen\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

    # Process the file
    echo -e "$header" > "$TAXA_DIR/BLAST_OTUs_all.txt"
    sed -e 's/;size=/\t/g' -e 's/;/\t/g' -e 's/|/\t/g' \
        -e 's/k__//g' -e 's/p__//g' -e 's/c__//g' \
        -e 's/o__//g' -e 's/f__//g' -e 's/g__//g' -e 's/s__//g' \
        "$TAXA_DIR/blast_output.txt" >> "$TAXA_DIR/BLAST_OTUs_all.txt"
    
    # Remove temporary file
    rm -f "$TAXA_DIR/blast_output.txt"
        
}

blast_ASVs_species() {
    
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$ASV_FASTA" \
        -db "$REFERENCE_SEQUENCES_SPECIES" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$TAXA_DIR/blast_output.txt" \
        -num_threads "$THREADS"

    # Define the header
    local header="OTU_ID\tread_abundance\treference_ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tspecies_hypothesis\tpident\tlength\tqlen\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

    # Process the file
    echo -e "$header" > "$TAXA_DIR/BLAST_ASVs_species.txt"
    sed -e 's/;size=/\t/g' -e 's/;/\t/g' -e 's/|/\t/g' \
        -e 's/k__//g' -e 's/p__//g' -e 's/c__//g' \
        -e 's/o__//g' -e 's/f__//g' -e 's/g__//g' -e 's/s__//g' \
        "$TAXA_DIR/blast_output.txt" >> "$TAXA_DIR/BLAST_ASVs_species.txt"
    
    # Remove temporary file
    rm -f "$TAXA_DIR/blast_output.txt"
        
}

blast_ASVs_genus() {

    # BLAST 5 best hits
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$ASV_FASTA" \
        -db "$REFERENCE_SEQUENCES_GENUS" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$TAXA_DIR/blast_output.txt" \
        -num_threads "$THREADS"
        
    # Define the header
    local header="OTU_ID\tread_abundance\treference_ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tspecies_hypothesis\tpident\tlength\tqlen\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

    # Process the file
    echo -e "$header" > "$TAXA_DIR/BLAST_ASVs_genus.txt"
    sed -e 's/;size=/\t/g' -e 's/;/\t/g' -e 's/|/\t/g' \
        -e 's/k__//g' -e 's/p__//g' -e 's/c__//g' \
        -e 's/o__//g' -e 's/f__//g' -e 's/g__//g' -e 's/s__//g' \
        "$TAXA_DIR/blast_output.txt" >> "$TAXA_DIR/BLAST_ASVs_genus.txt"
    
    # Remove temporary file
    rm -f "$TAXA_DIR/blast_output.txt"

}

blast_ASVs_all() {
    
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$ASV_FASTA" \
        -db "$REFERENCE_SEQUENCES_ALL" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$TAXA_DIR/blast_output.txt" \
        -num_threads "$THREADS"

    # Define the header
    local header="OTU_ID\tread_abundance\treference_ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tspecies_hypothesis\tpident\tlength\tqlen\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

    # Process the file
    echo -e "$header" > "$TAXA_DIR/BLAST_ASVs_all.txt"
    sed -e 's/;size=/\t/g' -e 's/;/\t/g' -e 's/|/\t/g' \
        -e 's/k__//g' -e 's/p__//g' -e 's/c__//g' \
        -e 's/o__//g' -e 's/f__//g' -e 's/g__//g' -e 's/s__//g' \
        "$TAXA_DIR/blast_output.txt" >> "$TAXA_DIR/BLAST_ASVs_all.txt"
    
    # Remove temporary file
    rm -f "$TAXA_DIR/blast_output.txt"
        
}

###############################################################################
## Main script ################################################################
###############################################################################

echo 'Starting at:' $(date)

# Activate the conda environment
source ~/.bashrc
conda activate dynamic_clustering

blast_OTUs_species
blast_OTUs_genus
blast_OTUs_all
blast_ASVs_species
blast_ASVs_genus
blast_ASVs_all

conda deactivate

echo 'Finished at:' $(date)

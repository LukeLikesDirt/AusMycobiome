#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=slurm/%x.%j.out

# Script:   Quality truncate, extract the ITS and quality filter Illumina
#           single-end reads targeting the ITS region.
# Purpose:  Prepare Illumina single-end reads for denoising.
# Author:   Luke Florence
# Date:     28th Janurary 2024

# Software:
# --------
# ITSxpress v2.0.0: https://github.com/USDA-ARS-GBRU/itsxpress
# VSEARCH v2.22.1: https://github.com/torognes/vsearch
# Trimmomatic v0.36: http://www.usadellab.org/cms/?page=trimmomatic
# FastQC v0.12.1: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# MultiQC v1.15: https://multiqc.info/

# Script Overview:
# --------------------------------------------------------------------------- #
# This script performs the following tasks:
#   (1) Quality truncate reads with Trimmomatic
#   (2) Extract the ITS region with ITSxpress
#   (3) Quality filter reads with VSEARCH
#   (4) Quality check read quality with FastQC and MultiQC
#   (5) Track reads across the pipeline using a library-wise approach

# Pre-requisites:
# --------------------------------------------------------------------------- #
# The raw reads should be demultiplexed and separated into unique 
# subdirectories based on sequencing run. The sequencing run subdirectories
# should be named 'run1', 'run2', 'run3', etc.

# Notes:
# --------------------------------------------------------------------------- #
# The chosen PCR primers of the Australian Microbiome, ITS1F and ITS4, targeted
# the full ITS region, leading to two amplicons (i.e. ITS1 for forward 
# sequences and ITS2 for reverse sequences) from which reads were generally
# too short to be merged. We therefore focus on forward reads that targt the 
# ITS1.
#
# Ideally, partially trimmed ITS1 reads (i.e., those containing the SSU region
# but with the 5.8S region removed) would be retained to capture fungi with long
# (>230 bp) ITS1 sequences. However, ASV calling on partially trimmed ITS reads
# is not recommended. To address underlying quality issues in the Australian 
# Microbiome data, we have chosen to use DADA2 and therefore do not retiain 
# partially trimmed ITS1. For guidance on retaining partially trimmed ITS1 reads,
# see: https://github.com/USDA-ARS-GBRU/itsxpress/issues/50
#
# ITS extraction for each libray can take in excess of one day depending on your 
# system – parallelise the ITSxpress step as needed.

# Constants
readonly NUM_THREADS=64           # The number of threads to use for parallel processing
readonly FIRST_RUN=1              # The first run to process
readonly LAST_RUN=42              # The last run to process
readonly FILE_EXT=".fastq.gz"     # The file extension of the raw reads
readonly ITS_REGION="ITS1"        # ITSxpress: The ITS region to extract. Options: ITS1, ITS2 or All
readonly TAXA="All"               # ITSxpress: The target taxonomic group: Alveolata, Bryophyta, Bacillariophyta, Amoebozoa, Euglenozoa, Fungi, Chlorophyta, Rhodophyta, Phaeophyceae, Marchantiophyta, Metazoa, Oomycota, Haptophyceae, Raphidophyceae, Rhizaria, Synurophyceae, Tracheophyta, Eustigmatophyceae, All. Default Fungi.
readonly CLUSTER=1.0              # ITSxpress: The identity threshold for clustering reads. Range 0.99-1.0. Set to 1 for exact de-replication. Default 1.0.
readonly WINDOW=4                 # Trimmomatic: The sliding window size for averaging quality scores
readonly QUAL=13                  # Trimmomatic: The quality threshold for sliding window trimming
readonly MAXEE=1                  # VSEARCH: Maximum expected error rate
readonly MAXN=0                   # VSEARCH: Maximum number of mismatches
readonly QMAX=41                  # VSEARCH: Maximum quality score
readonly MINLEN=50                # Various functions: Minimum length of reads

# Subdirectories
RAW_DATA="../data/bioinformatics/raw_data"
QUALITY_TRIMMED_DIR="../data/bioinformatics/quality_trimmed"
ITS_EXTRACTED_DIR="../data/bioinformatics/ITS_extracted"
QUALITY_FILTERED_DIR="../data/bioinformatics/quality_filtered"

# Define summary files for tracking reads
TRACK_READS_TRIMMED="../data/bioinformatics/summary_trimmed.txt"
TRACK_READS_ITS="../data/bioinformatics/summary_its.txt"
TRACK_READS_QUAL_FILTERED="../data/bioinformatics/summary_quality_filtered.txt"

# Log function
log() {
    local timestamp
    timestamp=$(date)
    echo -e "\n$1 $timestamp\n"
}

# Function to check if a directory exists, and if not, create it
create_directory() {
    local directory="$1"
    if [ ! -d "$directory" ]; then
        mkdir -p "$directory"
    fi
}

# Function to perform quality trimming with Trimmomatic
quality_trim() {
    local run_number="$1"
    log "Starting quality trimming with Trimmomatic for run $run_number"

    local input_dir="$RAW_DATA/run$run_number"
    local output_dir="$QUALITY_TRIMMED_DIR/run$run_number"

    create_directory "$output_dir"

    for file in "$input_dir"/*"$FILE_EXT"; do
        local filename=$(basename "$file")
        log "Trimming file $filename"
        
        trimmomatic SE -threads "$NUM_THREADS" "${file}" "$output_dir/${filename}" \
            SLIDINGWINDOW:"$WINDOW":"$QUAL" MINLEN:"$MINLEN"
        
        log "Completed trimming file $filename"
    done

    log "Completed quality trimming for run $run_number"
}

# Function to extract ITS region with ITSxpress
extract_its() {
    local run_number="$1"
    log "Starting ITS extraction with ITSxpress for run $run_number"

    local input_dir="$QUALITY_TRIMMED_DIR/run$run_number"
    local output_dir="$ITS_EXTRACTED_DIR/run$run_number"

    create_directory "$output_dir"

    for file in "$input_dir"/*"$FILE_EXT"; do
        local filename=$(basename "$file")
        log "Extracting ITS region from file $filename"

        itsxpress \
            --single_end \
            --fastq "${file}" \
            --cluster_id "$CLUSTER" \
            --region "$ITS_REGION" \
            --taxa "$TAXA" \
            --log "$output_dir/logfile.txt" \
            --outfile "$output_dir/${filename}" \
            --threads "$NUM_THREADS"

        log "Completed ITS extraction from file $filename"
    done

    log "Completed ITS extraction with ITSxpress for run $run_number"
}

# Function to perform quality filtering with VSEARCH
quality_filter() {
    local run_number="$1"
    log "Starting quality filtering with VSEARCH for run $run_number"

    local input_dir="$ITS_EXTRACTED_DIR/run$run_number"
    local output_dir="$QUALITY_FILTERED_DIR/run$run_number"

    create_directory "$output_dir"

    for file in "$input_dir"/*.fastq.gz; do
        local filename=$(basename "$file")

        vsearch \
            --fastq_filter "${file}" \
            --fastq_maxee $MAXEE \
            --fastq_minlen $MINLEN \
            --fastq_maxns $MAXN \
            --fastq_qmax $QMAX \
            --fastqout "$output_dir/${filename%$FILE_EXT}.fastq"

    done

    log "Completed quality filtering for run $run_number"

    log "Compressing the quality-filtered reads"

    for file in "$output_dir"/*"fastq"; do

        gzip -c "$file" > "${file}.gz"
        rm "$file"

    done

    log "Finised compressing the quality-filtered reads"

}

# Function to generate read quality report for each sample using FastQC and MultiQC
generate_quality_report() {
    local run_number="$1"

    log "Starting to generate the quality report for run $run_number"
    
    fastqc "$ITS_EXTRACTED_DIR/run$run_number"/*FILE_EXT -o "$ITS_EXTRACTED_DIR/run$run_number"
    multiqc "$ITS_EXTRACTED_DIR/run$run_number" -o "$ITS_EXTRACTED_DIR/run$run_number"
    
    # Remove intermediate files and 'multiqc_data' directory
    rm "$ITS_EXTRACTED_DIR/run$run_number"/*fastqc.zip "$ITS_EXTRACTED_DIR/run$run_number"/*fastqc.html
    rm -r "$ITS_EXTRACTED_DIR/run$run_number"/multiqc_data

    log "Completed generating a quality report for run $run_number"
}

###############################################################################
### Main script ###############################################################
###############################################################################

log "Starting at:"

### Exicute functions #########################################################

# Activate Conda environment
source ~/.bashrc
conda activate dynamic_clustering

# Process the fastq files
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    quality_trim "$run_number"
    extract_its "$run_number"
    quality_filter "$run_number"
    generate_quality_report "$run_number"
done

# Deactivate Conda environment
conda deactivate

### Track the reads across the pipeline #######################################

# Open the file for writing
output_file="read_track_summary.txt"

# Write the headers to the file
echo -e "Run Number\tRaw Sequences\tQuality Trimmed\tITS Extracted\tQuality Filtered" > "$output_file"

# Loop through the runs and count sequences
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    # Count raw sequences
    raw_seq_count=$(($(zcat "$RAW_DATA/run$run_number"/*"$FILE_EXT" | wc -l) / 4))

    # Count quality trimmed sequences
    trimmed_seq_count=$(($(zcat "$QUALITY_TRIMMED_DIR/run$run_number"/*"$FILE_EXT" | wc -l) / 4))

    # Count ITS extracted sequences
    its_extracted_seq_count=$(($(zcat "$ITS_EXTRACTED_DIR/run$run_number"/*"$FILE_EXT" | wc -l) / 4))

    # Count quality filtered sequences
    quality_filtered_seq_count=$(($(zcat "$QUALITY_FILTERED_DIR/run$run_number"/*"$FILE_EXT" | wc -l) / 4))

    # Append the row for the current run to the file
    echo -e "$run_number\t$raw_seq_count\t$trimmed_seq_count\t$its_extracted_seq_count\t$quality_filtered_seq_count" >> "$output_file"
done

log "Finished at:"
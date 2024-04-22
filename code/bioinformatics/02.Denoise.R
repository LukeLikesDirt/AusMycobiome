
# Script:   Denoise Illumina single-end reads.
# Purpose:  Prepare Illumina forwards reads for chimera detection in VSEARCH.
# Credit:   https://benjjneb.github.io/dada2/bigdata.html
# Author:   Luke Florence.
# Date:     28th Janurary 2024.
#
# Script Purpose
# --------------------------------------------------------------------------- #
# This R script is designed for the denoising of Illumina forward reads using 
# DADA2. The script is designed for "big data" projects, and therefore, performs
# quality filtering and denoising on sequencing run individually to improve 
# error rate estimations. This script converts the DADA2 'rds' sequence table 
# output to a 'fasta' file for chimera detection and removal in VSEARCH.
#
# Script Contents:
# --------------------------------------------------------------------------- #
#   (1) Denoise the reads
#   (2) Merges the denoised sequence tables
#   (3) Convert the sequence table to a FASTA file formatted for VSEARCH

# Required packages
library(dada2)
library(seqinr)
library(here)
library(tidyverse)

# Define the number of sequencing runs to be processed
num_runs <- 42

# Define the path to the ITS data directory
path <- here("data/bioinformatics/")

# Create subdirectories
dir.create(file.path(path, '05.Denoised'))
dir.create(file.path(path, '07.Chimera_filtered'))

# Create an empty data frame to store the read tracking summary
summary_track <- data.frame()

### 1. Denoise the reads #######################################################

# Denoise each run individually
for (run in 1:num_runs) {

  ## Input and output directories
  qualFilt_dir <- file.path(path, "04.Quality_filtered",
                            paste0("run", run))
  denoised_dir <- file.path(path, "05.Denoised")

  ## Denoising
  filts <- list.files(qualFilt_dir, pattern = 'fastq.gz', full.names = TRUE)
  # Extract prefix and paste the run number
  sample.names <- sapply(strsplit(basename(filts), '_'), function(x) {
    prefix <- paste(x[1], sep = "_")  # Grab the filename prefix
    result <- paste(prefix, run, sep = "_")  # Append run ID to the prefix
    return(result)
  })
  names(filts) <- sample.names

  ## Learn error rates
  set.seed(1986)
  err <- learnErrors(filts, nbases = 1e8, multithread = TRUE,
                     randomize = TRUE)

  ## Infer sequence variants
  dds <- vector('list', length(sample.names))
  names(dds) <- sample.names
  for (sam in sample.names) {
    cat('Processing:', sam, '\n')
    derep <- derepFastq(filts[[sam]])
    dds[[sam]] <- dada(derep, err = err, multithread = TRUE,
                       DETECT_SINGLETONS = TRUE)
  }

  ## Build sequence table
  seqtab <- makeSequenceTable(dds)
  saveRDS(seqtab, file.path(denoised_dir, paste0(run, '_seqtab.rds')))

  ## Track reads through the DADA2 pipeline for the current run
  getN <- function(x) sum(getUniques(x))

  ## Check if there are rows in the "stats.csv" file
  track <- data.frame(Run = run,
                      sample = sample.names,
                      output = sapply(dds, getN)) %>%
    filter(substr(sample, 1, 1) == "s") %>%
    select(-sample) %>%
    as.data.frame()

  ## Append the track information for the current run to the summary_track file
  summary_track <- rbind(summary_track, track)

}

summary_track %>%
  as.data.frame() %>%
  write.table(file.path(path, "summary_denoised.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

### 2. Merge the seuqnce tables ################################################

# Merge all sequence tables and convert rds to fasta format for chimera
# detection in VSEARCH

# NOTE: Some samples have been re-sequenced across multiple sequencing runs.
# Therefore, I sum read counts and merge re-sequenced samples when executing the
# mergeSequenceTables()function.

# Initialise an empty list to store the sequence tables
seqtab_list <- list()

# Loop through the sequencing runs to read and merge the sequence tables
for (run in 1:num_runs) {
  seqtab_file <- file.path(path, paste0('05.Denoised/',
                                        run, '_seqtab.rds'))
  seqtab <- readRDS(seqtab_file)
  seqtab_list[[run]] <- seqtab
}

# Merge all sequence tables with repeats = "sum" to combine abundance across
# samples that have been sequenced in multiple runs
all.seqtab <- mergeSequenceTables(tables = seqtab_list)

### 3. Convert the sequence table to a FASTA file formatted for VSEARCH ########

## Convert rds to fasta for chimera detection in VSEARCH

## Format data frame
fasta.tab <- as.data.frame(all.seqtab) %>%
  rownames_to_column('sample.names') %>%
  as_tibble() %>%
  pivot_longer(-sample.names, names_to = 'seq', values_to = 'size') %>%
  filter(size > 0) %>%
  ungroup() %>%
  mutate(seq.name = paste0(sample.names, '.fasta', '.', row_number(),
                           ';size=', size)) %>%
  select(seq.name, seq)

## Save the rds file for the merged sequence table
saveRDS(all.seqtab, file.path(path, '05.Denoised/all_seqtab.rds'))

## Write and save the fasta file for the merged sequence table
write.fasta(as.list(fasta.tab$seq), fasta.tab$seq.name,
            file.path(path, '07.Chimera_filtered/all.fasta'),
                      open = 'w', nbchar = 60, as.string = FALSE)
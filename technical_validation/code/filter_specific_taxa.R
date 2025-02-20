
# Required packages
require(Biostrings)
require(data.table)
require(tidyverse)

# Path to input files
taxonomy_file <- "../output/taxonomy_pseudo.txt"
otu_table_file <- "../output/otu_table.txt"
sample_metadata_file <- "../output/sample_metadata.txt"
fasta_file <- "../output/sequences.fasta"

# Path to the output directory
output_dir <- "../output/"

# Taxon of interest and rank
taxon <- "Mycena"
rank <- "genus"

# Filter the taxonomy file to the taxon of interest
taxonomy_filtered <- fread(taxonomy_file) %>%
  filter(.data[[rank]] == taxon)

# Filter the OTU table to the taxon of interest
otu_table_filtered <- fread(otu_table_file) %>%
  filter(OTU_ID %in% taxonomy_filtered$OTU_ID) %>%
  # Remove redundant samples
  pivot_longer(
    -OTU_ID,
    names_to = "sample_id",
    values_to = "abundance"
    ) %>%
  filter(abundance > 0) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = abundance,
    values_fill = 0
    )

# Read in the fasta file
fasta <- readDNAStringSet(fasta_file)

# Filter the fasta file to the taxon of interest
fasta_filtered <- fasta[names(fasta) %in% otu_table_filtered$OTU_ID]

# Filter the sample metadata to the taxon of interest
sample_metadata_filtered <- fread(sample_metadata_file) %>%
  filter(sample_id %in% names(otu_table_filtered))

# Write the filtered files to the output directory
fwrite(
  taxonomy_filtered,
  file = paste0(output_dir, "taxonomy_", taxon, ".txt"),
  sep = "\t"
  )
fwrite(
  otu_table_filtered,
  file = paste0(output_dir, "otu_table_", taxon, ".txt"),
  sep = "\t"
  )
writeXStringSet(
  fasta_filtered,
  file = paste0(output_dir, "sequences_", taxon, ".fasta")
  )
fwrite(
  sample_metadata_filtered,
  file = paste0(output_dir, "sample_metadata_", taxon, ".txt"),
  sep = "\t"
  )


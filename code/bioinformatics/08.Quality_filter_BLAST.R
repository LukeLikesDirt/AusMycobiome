
# Script:   Remove poor quality BLAST hits.
# Author:   Luke Florence
# Date:     31st January 2024
#
# Contents:
#   (1) Filter taxa to OTUs 
#   (2) Remove poor quality BLAST hits
#   (3) Subset taxa tables based on phylum for manual inspection of top hits
#   (4) Save the quality filtered taxa table
#
# Required packages and functions
require(data.table)
require(tidyverse)
source("code/bioinformatics/functions.R")

### (1) Filter taxa to OTUs ###################################################

# Read in OTUs
OTUs <- fread(
  "data/bioinformatics/08.OTUs/OTUs_abundance_filtered.csv"
  ) %>%
  select(OTU_ID)

# Read in the taxa tables and filter to OTUs
taxa <- left_join(
  OTUs,
  fread("data/bioinformatics/09.Taxonomy/BLAST_best_10.csv", sep = ";"),
  by = "OTU_ID"
  ) %>%
  # Rename taxa levels - remove prefixes
  mutate_all(~gsub("^.*__", "", .)) %>%
  # Abundance, similarity, coverage (to percent) and e-value as numeric
  mutate(
    similarity = as.numeric(pident),
    coverage = (as.numeric(length) / as.numeric(slen)) * 100,
    coverage = ifelse(coverage > 100, 100, coverage),
    evalue = as.numeric(evalue),
    species_hypothesis = sapply(strsplit(species, "\\|"), `[`, 2),
    species = sapply(strsplit(species, "\\|"), `[`, 1)
  ) %>%
  # Select values of interest
  select(
    OTU_ID, kingdom, phylum, class, order, family, genus, species,
    species_hypothesis, similarity, coverage, evalue
  ) %>%
  glimpse(.)

### (2) Remove poor quality hits ##############################################

# Summaries and manually inspect phylum level coverage and similarity across
# abundant taxa:

# Summaries similarity statistics for phylum
similarity_phylum <- taxa %>%
  filter(!is.na(phylum)) %>%
  group_by(phylum) %>%
  summarise(
    mean_similarity = mean(similarity),
    median_similarity = median(similarity),
    sd_similarity = sd(similarity),
    q25_similarity = quantile(similarity, 0.25),
    q75_similarity = quantile(similarity, 0.75),
    lower_limit_outlier = max(q25_similarity - 1.5 * IQR(coverage), 0),
    upper_limit_outlier = min(q75_similarity + 1.5 * IQR(coverage), 100),
    n = n()
  ) %>%
  print(n = Inf)

# Summaries coverage statistics for phylum
coverage_phylum <- taxa %>%
  filter(!is.na(phylum)) %>%
  group_by(kingdom, phylum) %>%
  summarise(
    mean_coverage = mean(coverage),
    median_coverage = median(coverage),
    sd_coverage = sd(coverage),
    q25_coverage = quantile(coverage, 0.25),
    q75_coverage = quantile(coverage, 0.75),
    lower_limit_outlier = max(q25_coverage - 1.5 * IQR(coverage), 0),
    upper_limit_outlier = min(q75_coverage + 1.5 * IQR(coverage), 100),
    n = n()
  ) %>%
  print(n = Inf)

# Which fungal phylum have q25 coverage < 70%
coverage_phylum %>%
  filter(
    kingdom == "Fungi",
    q25_coverage < 70
  )

# I will filter taxa table to phylum based on a similarity and coverage
# thresholds of 70%. I include a coverage exception where q25_coverage is < 70%,
# in which case I will use the q25_coverage values as the threshold. I also
# remove similarity outliers. Hits within this "quality filtered" taxa table
# will later be used to annotate OTUs to a single consensus hit.

# Create a quality filter at level phylum
phylum_filter <- taxa %>%
  filter(!is.na(phylum)) %>%
  group_by(phylum) %>%
  mutate(
    q25_coverage = quantile(coverage, 0.25),
    q25_similarity = quantile(similarity, 0.25),
    low_outlier_similarity = max(q25_similarity - 1.5 * IQR(similarity), 0)
  ) %>%
  ungroup() %>%
  select(phylum, q25_coverage, low_outlier_similarity)

# Quality filter the taxa table and save the abundance and quality filtered taxa
# table
taxa_phylum <- taxa %>%
  filter(!is.na(phylum)) %>%
  filter(
    # Coverage threshold of 70% or q25 values for phylum with poor coverage due
    # to the combination of long ITS regions and short-read sequencing
    (coverage >= 70 | 
       (phylum %in% phylum_filter$phylum &
          coverage > phylum_filter$q25_coverage)) &
      # Similarity threshold of of 70% or similarity outliers which indicate 
      # poor exceptionally quality hits
      (similarity > 70 | 
         (phylum %in% phylum_filter$phylum &
            similarity > phylum_filter$low_outlier_similarity))
  )

# See what we have retain: 89.3% of the original hits have been retained
nrow(taxa_phylum) / nrow(taxa)

### (3) Subset taxa tables based on phylum ####################################

# Subset taxa table based on phylum for manual inspection.
# The phylum sub-setting is based on OTU best hit.
# Create a list of unique fungal phyla names
fungal_phyla <- taxa_phylum %>%
  select(kingdom, phylum) %>%
  filter(
    kingdom == "Fungi"
  ) %>%
  select(phylum) %>%
  unique() %>%
  print(n = Inf)

# Initialise an empty list to store the data frames
phyla_dataframes <- list()

# Loop through each phylum, subset by phyla and save phyla-specific taxa tables
# The subsets are based on the best hit, but the corresponding subset will
# contain all the hits associated to the given OTU

# Loop through each phylum and save the corresponding data frames
for (phyla in fungal_phyla$phylum) {
  
  # Grab the best hit for each OTU
  first_otu_ids <- taxa_phylum %>%
    group_by(OTU_ID) %>%
    slice_head(n = 1) %>%
    filter(phylum == phyla) %>%
    pull(OTU_ID)
  
  # Subset taxa tables based on best hits, but include all hits for each OTU
  phyla_df <- taxa_phylum %>%
    filter(OTU_ID %in% first_otu_ids)
  
  # Save the data frame to a CSV
  csv_file <- paste0(
    "data/bioinformatics/09.Taxonomy/phylum_",
    gsub(" ", "_", phyla), ".csv"
  )
  
  fwrite(phyla_df, csv_file)
  
  # Append the data frame to the list
  phyla_dataframes[[phyla]] <- phyla_df
}

### (4) Save the quality and abundance filtered taxa table #####################

# We re-annotate "GS01_phy_Incertae_sedis" to "GS01" for interoperability.
# We also re-annotate the beetle genus "Metazoa; Arthropoda; Insecta; 
# Coleoptera; Alexiidae; Sphaerosoma" to the ectomycorrrhiza genus "Fungi; 
# Ascomycota; Pezizomycetes; Pezizales; Pyronemataceae; Sphaerosoma". This is
# a bug in the UNITE+INDC v.9.0 databases that we found during our our analysis.
# This bug was confirmed  by Kessy Abarenkov (a UNITE developer) and will be 
# fixed in UNITE+INDC v10.0.

# Reannotate and save the quality and abundance filtered taxa table
taxa_phylum %>%
  mutate(
    phylum = case_when(
      phylum == "GS01_phy_Incertae_sedis" ~ "GS01",
      TRUE ~ phylum
    ),
    kingdom = case_when(
      genus == "Sphaerosoma" ~ "Fungi",
      TRUE ~ kingdom
    ),
    phylum = case_when(
      genus == "Sphaerosoma" ~ "Ascomycota",
      TRUE ~ phylum
    ),
    class = case_when(
      genus == "Sphaerosoma" ~ "Pezizomycetes",
      TRUE ~ class
    ),
    order = case_when(
      genus == "Sphaerosoma" ~ "Pezizales",
      TRUE ~ order
    ),
    family = case_when(
      genus == "Sphaerosoma" ~ "Pyronemataceae",
      TRUE ~ family
    )
  ) %>%
  fwrite(
    ., "data/bioinformatics/09.Taxonomy/taxa_abundance_quality_filtered.csv"
  )

# Clear all objects from the environment
rm(list = ls())

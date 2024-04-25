
# Script:   Abundance filter for OTUs and samples.
# Author:   Luke Florence
# Date:     31st January 2024
#
# Contents:
#   (1) Sample-wise abundance filter @ 0.1%
#   (2) Library-wise abundance filter @ 0.5%
#   (3) Filter positive control OTUs library-wise @ 3%
#   (3) Low abundance sample filter @ 5,000 reads
#
# Note that the filtering thresholds used here have been established based on 
# in-depth analyses on this dataset and are likely to not be appropriate for
# other datasets. Specifically, the thresholds are very stringent because of 
# the poor quality of the data.

# Required packages and functions
require(data.table)
require(ggpubr)
require(terra)
require(tidyverse)
source("code/bioinformatics/functions.R")

# Metadata
metadata <- fread("data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

### 1. Sample-wise abundance filter ##########################################

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "data/bioinformatics/08.OTUs/OTUs.txt", header = TRUE, sep = "\t"
  ) %>%
  filter_low_abundance_otus(., threshold = 0.1)

### 2. Library-wise abundance filter ########################################

# Filter OTUs library-wise using an OTU abundance threshold 0.5%.

# Controls are not in the metadata so I will grab control sample names from the
# OTU table
flow_id_data <- 
  inner_join_colnames(
    
    OTUs %>%
      select(starts_with(c("Mock", "Neg", "Pos"))) %>%
      rownames_to_column(var = "OTU_ID") %>%
      pivot_longer(-OTU_ID, names_to = "sample_id") %>%
      select(sample_id) %>%
      mutate(
        flow_id = str_extract(sample_id, "(?<=_)\\d+"),
      ) %>%
      unique(),
    
    metadata %>%
      select(sample_id, flow_id) %>%
      as_tibble()
    
  ) %>%
  glimpse()

# Unique flow ID to vector
flow_ids <- unique(flow_id_data$flow_id)

# List to store filtered OTU tables for each flow_id
filtered_otu_tables <- list()

# Loop over each flow_id and...
for (lib in flow_ids) {
  # Step 1: Subset sample_IDs by flow_id
  subset_metadata <- flow_id_data %>%
    filter(flow_id == lib) %>%
    select(sample_id)
  
  # Step 2: Subset OTU table based on unique sample_IDs within each flow_ID
  subset_otu_table <- OTUs %>%
    select(any_of(subset_metadata$sample_id))
  
  # Step 3: Remove rows where row sums are equal to 0
  filtered_otu_table <- subset_otu_table %>%
    filter(rowSums(.) > 0)
  
  # Step 4: Apply the library-wise abundance filter function to each library
  filtered_otu_table <- library_wise_filter(filtered_otu_table, 0.5)
  
  # Step 5: Generate a filtered OTU table for each library
  assign(paste("otu_", lib, sep = ""), filtered_otu_table)
  
  # Step 6: Store filtered OTU tables in a list for later use in the join 
  # function
  filtered_otu_tables[[lib]] <- filtered_otu_table
  
}

# Initialise the combined OTU table
OTUs_filtered <- data.frame(OTU_ID = character())

# Loop over each flow_id and join tables
for (lib in flow_ids) {
  OTUs_filtered <- full_join(OTUs_filtered,
                             filtered_otu_tables[[lib]],
                             by = "OTU_ID") %>%
    mutate_all(~replace_na(., 0))
}

# Check that we still have the same number of samples: Looking good
ncol(OTUs_filtered)
ncol(OTUs)

### 3. Remove positive control OTUs ###########################################

# Above found that OTU prevalence was associated with the relative abundance of
# OTUs in positive control samples after the above filtering. However, we also
# found that OTUs in positive control samples were true positives in some test
# samples, potentially due to over clustering. Rather than removing all positive
# control OTUs I will remove those that are present at 3% within a sample on a 
# library-by-library basis.

# Define control OTUs
control_OTUs_IDs <- OTUs_filtered %>%
  select(OTU_ID, starts_with("Mock"), starts_with("Pos")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  pull(OTU_ID) %>%
  unique()

# List to store filtered OTU tables for each flow_id
filtered_otu_tables_control <- list()

# Loop over each flow_id and...
for (lib in flow_ids) {
  # Step 1: Subset sample_IDs by flow_id
  subset_flow_id <- flow_id_data %>%
    filter(flow_id == lib) %>%
    select(sample_id)
  
  # Step 2: Subset OTU table based on unique sample_IDs within each flow_ID
  subset_otu_table <- OTUs_filtered %>%
    column_to_rownames(var = "OTU_ID") %>%
    select(any_of(subset_flow_id$sample_id))
  
  # Step 3: Remove rows where row sums are equal to 0
  filtered_otu_table <- subset_otu_table %>%
    filter(rowSums(.) > 0)
  
  # Step 4: Apply the library-wise abundance filter function to each library for
  # the control OTUs only
  filtered_otu_table <- library_wise_filter_controls(
    filtered_otu_table, 3, control_OTUs_IDs
  )
  
  # Step 5: Generate a filtered OTU table for each library
  assign(paste("otu_", lib, sep = ""), filtered_otu_table)
  
  # Step 6: Store filtered OTU tables in a list for later use in the join 
  # function
  filtered_otu_tables_control[[lib]] <- filtered_otu_table
  
}

# Initialise the combined OTU table
OTUs_filtered_control <- data.frame(OTU_ID = character())

# Loop over each flow_id and join tables
for (lib in flow_ids) {
  OTUs_filtered_control <- full_join(OTUs_filtered_control,
                                     filtered_otu_tables_control[[lib]],
                                     by = "OTU_ID") %>%
    mutate_all(~replace_na(., 0))
}

### 4. Low abundance sample filter ############################################

# Here I will remove samples with a sequencing depth of less than 5,000 reads.
# I also dereplicate re-sequenced samples by retaining the sample with the
# lowest richness-to-abundance ratio. This approach is atypical and conservative
# but we consider this to be the most appropriate approach for this dataset.

# Remove low abundance samples and summaries sequencing depth per sample
seq_depth <- OTUs_filtered_control %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance >= 5000) %>%
  arrange(abundance) %>%
  print(n = 200)

# Grab the re-sequenced samples
resequenced_samples <- seq_depth %>%
  mutate(
    sample_flow_id = sample_id,
    flow_id = str_extract(sample_id, "(?<=_)\\d+"),
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  filter(duplicated(sample_id) | duplicated(sample_id, fromLast = TRUE)) %>%
  as_tibble() %>%
  arrange(sample_id) %>%
  print(n = Inf)

# Retain the sample with the lowest richness-to-abundance ratio
resequenced_uniques <- resequenced_samples %>%
  mutate(
    ratio = richness/abundance
  ) %>%
  arrange(ratio) %>%
  group_by(sample_id) %>% 
  slice(1) %>%
  ungroup() %>%
  # Rebuild the sample_id name so to contain the flow_id suffix
  mutate(sample_flow_id = paste(sample_id, flow_id, sep = "_")) %>%
  select(-c(flow_id, ratio)) %>%
  print(n = Inf)

# Create a filter for the samples that I'll remove
resequenced_filter <- resequenced_samples %>%
  filter(!sample_flow_id %in% resequenced_uniques$sample_flow_id) %>%
  select(sample_id = sample_flow_id)

# Filter the sequencing depth summary to selected samples
seq_depth_filtered <- seq_depth %>%
  filter(!sample_id %in% resequenced_filter$sample_id)

# Check to see if we retained the 2101 unique samples
nrow(seq_depth_filtered)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_samples <- OTUs_filtered_control %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_samples,
  "data/bioinformatics/08.OTUs/OTUs_abundance_filtered.csv"
)

# Clear all objects from the environment
rm(list = ls())

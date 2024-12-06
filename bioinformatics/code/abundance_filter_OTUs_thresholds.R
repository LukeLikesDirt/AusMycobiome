
# Required packages and functions
require(data.table)
require(Biostrings)
require(ggpubr)
require(tidyverse)

### Threshold 01 ##############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()


# OTUs
OTUs <- fread(
  "../data/OTUs/OTUs.txt"
  )

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- OTUs %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
  ) %>%
  filter_low_abundance_otus(., threshold = 0.01)

# What have we lost?
# At 0.01% we removed 55,816 OTUs or 23.6% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 1,346,518 or 0.66% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_01.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_01.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 59,076 or 25% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 25,743,865 or 12.6% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_01.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 02 ##############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.02% (2 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 2 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.02)

# What have we lost?
# At 0.1% we removed 91,280 OTUs or 38.6% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 3,327,232 or 1.6% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_02.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_02.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 93,895 or 39.7% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 27,427,785 or 13.4% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_02.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 05 ##############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 5 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.05)

# What have we lost?
# At 0.1% we removed 139,129 OTUs or 58.9% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 9,119,769 or 4.4% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_05.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_05.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 140,765 or 59.6% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 32,255,730 or 15.8% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_05.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 10 ##############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.1)

# What have we lost?
# At 0.1% we removed 169,850 OTUs or 71.8% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 16,853,596 or 8.2% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_10.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_10.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 170,904 or 72.3% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 39,031,186 or 19.1% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_10.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 20 ##############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.2% (20 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 20 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.2)

# What have we lost?
# At 0.1% we removed 193,284 OTUs or 81.79% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 27,893,650 or 13.66% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_20.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_20.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 194,065 or 82.1% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 48,931,356 or 23.4% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_20.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 50 ##############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.5% (50 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 50 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.5)

# What have we lost?
# At 0.5% we removed 213,567 OTUs or 90.37% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 47,660,184 or 23.3% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_50.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_50.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 214,123 or 90.6% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 67,241,510 or 32.9% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_50.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 01_1 #############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.01)

# What have we lost?
# At 0.1% we removed 113,672 OTUs or 74.9% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 11,682,397 or 6.6% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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
  filtered_otu_table <- library_wise_filter(filtered_otu_table, 1)
  
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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_01_1.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_01_1.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 114,693 or 75.6% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 32,541,646 or 18.4% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_01_1.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 02_1 #############################################################


source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.02)

# What have we lost?
# At 0.1% we removed 113,672 OTUs or 74.9% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 11,682,397 or 6.6% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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
  filtered_otu_table <- library_wise_filter(filtered_otu_table, 1)
  
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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_02_1.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_02_1.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 114,693 or 75.6% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 32,541,646 or 18.4% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_02_1.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 05_1 #############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.05)

# What have we lost?
# At 0.1% we removed 113,672 OTUs or 74.9% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 11,682,397 or 6.6% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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
  filtered_otu_table <- library_wise_filter(filtered_otu_table, 1)
  
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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_05_1.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_05_1.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 114,693 or 75.6% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 32,541,646 or 18.4% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_05_1.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 10_1 #############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.1)

# What have we lost?
# At 0.1% we removed 113,672 OTUs or 74.9% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 11,682,397 or 6.6% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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
  filtered_otu_table <- library_wise_filter(filtered_otu_table, 1)
  
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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_10_1.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_10_1.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 114,693 or 75.6% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 32,541,646 or 18.4% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_10_1.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 20_1 #############################################################

source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.2)

# What have we lost?
# At 0.1% we removed 113,672 OTUs or 74.9% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 11,682,397 or 6.6% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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
  filtered_otu_table <- library_wise_filter(filtered_otu_table, 1)
  
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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_20_1.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_20_1.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 114,693 or 75.6% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 32,541,646 or 18.4% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_20_1.csv"
)

# Clear all objects from the environment
rm(list = ls())

### Threshold 50_1 #############################################################
source("abundance_filter_functions.R")

### (1) Sample-wise abundance filter ##########################################

# Filter OTUs sample-wise using threshold of 0.1% (10 in 10,000).

# Metadata
metadata <- fread("../../technical_validation/data/metadata.csv") %>%
  mutate(flow_id = as.character(flow_id)) %>%
  glimpse()

# Check control sample richness: If samples are clean richness will = 10
mock_controls <- fread(
  "../data/OTUs/OTUs_97.txt") %>%
  select(OTU_ID, starts_with("Mock")) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    richness = n()
  ) %>%
  print(n = Inf)

# # Calculate the total number of OTUs before filtering: 236,303
total_OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  nrow(.) %>%
  print(.)

# Calculate the total abundance: 204,134,208
total_abundance <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  sum(.) %>%
  print(.)

# Calculate the total number of samples: 2,536
total_samples <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  column_to_rownames(var = "OTU_ID") %>%
  ncol(.) %>%
  print(.)

# Read in the OTU table and filter OTUs at 10 in 10,000
OTUs <- fread(
  "../data/OTUs/OTUs_97.txt"
) %>%
  filter_low_abundance_otus(., threshold = 0.5)

# What have we lost?
# At 0.1% we removed 113,672 OTUs or 74.9% compared to the unfiltered table
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# This corresponds to 11,682,397 or 6.6% of the total abundance
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Check mock sample richness after 0.1% sample-wise abundance filter
mock_controls_filtered_samples <- inner_join(
  mock_controls,
  (OTUs %>%
     rownames_to_column(var = "OTU_ID") %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_samples = n()
     )),
  by = "sample") %>%
  print(n = Inf)

### (2) Library-wise abundance filter ########################################

# Here I will remove OTUs across each library using an OTU abundance threshold
# 0.5%.

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
  filtered_otu_table <- library_wise_filter(filtered_otu_table, 1)
  
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

# Check that we still have the same number of samples:
# The extra sample is the row names occupying the first column
ncol(OTUs_filtered)
ncol(OTUs)

# Check mock community richness after library-wise filtering:
mock_controls_filtered_lib <- inner_join(
  mock_controls_filtered_samples,
  (OTUs_filtered %>%
     select(OTU_ID, starts_with("Mock")) %>%
     pivot_longer(-OTU_ID, names_to = "sample") %>%
     filter(value > 0) %>%
     group_by(sample) %>%
     summarise(
       richness_filtered_lib = n()
     )),
  by = "sample") %>%
  print(n = Inf)
# If the filtering has worked well, I expect 10 OTUs per mock sample
median(mock_controls_filtered_lib$richness_filtered_lib)
mean(mock_controls_filtered_lib$richness_filtered_lib)

# Save mock community data
mock_controls <- mock_controls_filtered_lib %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5)
  ) %>%
  arrange(flow_id) %>%
  group_by(flow_id) %>%
  mutate(replicate = row_number()) %>%
  select(
    flow_id, flow_id_verbatim, replicate,
    richness_no_filter = richness,
    richness_sample_wise_filter = richness_filtered_samples,
    richness_lib_wise_filter = richness_filtered_lib) %>%
  print(n = Inf) %>%
  write_csv(., "../../output/supplementary/mock_samples_fungi_50_1.csv")

# Summarise richness and sequencing depth of controls
all_controls <- OTUs_filtered %>%
  select(OTU_ID, starts_with(c("Neg", "Pos", "Mock"))) %>%
  pivot_longer(-OTU_ID, names_to = "sample") %>%
  filter(value > 0) %>%
  group_by(sample) %>%
  summarise(
    seq_depth = sum(value),
    richness = n()
  ) %>%
  ungroup() %>%
  mutate(
    flow_id = str_extract(sample, "(?<=_)\\d+"),
    sample_id = str_extract(sample, "(.*?)(?=_)"),
    flow_id_verbatim = str_sub(sample_id, start = -5),
    sample_id = str_sub(sample_id, end = -6),
  ) %>%
  select(flow_id, flow_id_verbatim, sample_id, seq_depth, richness) %>%
  arrange(seq_depth) %>%
  arrange(flow_id) %>%
  print(n = 100) %>%
  fwrite(., "../../output/supplementary/control_samples_fungi_50_1.csv")

### (3) Low abundance sample filter ############################################

# Remove low abundance samples and summarise sequencing depth per sample
seq_depth <- OTUs_filtered %>%
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
seq_depth  %>%
  print(n = 200)
# Calculate the number of unique samples after low abundance filtering
unique_samples <- seq_depth %>%
  select(sample_id) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
  ) %>%
  unique()

# There are 2109 unique samples
n_samples <- unique(unique_samples$sample_id) %>%
  length() %>%
  print()

# Take a look at what we lost, mainly for shits and giggles:
# Most in the 1,000 to 5,000 range are from Antarctica
seq_depth_low <- OTUs_filtered %>%
  select(OTU_ID, starts_with("s")) %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n()
  ) %>%
  # Define abundance threshold
  filter(abundance < 5000) %>%
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

# Calculate the number of re-sequenced samples: 93
n_dups <- unique(resequenced_samples$sample_id) %>%
  length() %>%
  print()

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

# Check to see if we retained the 2109 unique samples
nrow(seq_depth_filtered)

seq_depth_filtered %>% print(n = 30)

# Filter the the OTU table to unique test sample and remove rows where OTU
# abundance is equal to zero
OTUs_filtered_1 <- OTUs_filtered %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(
    sample_id %in% seq_depth_filtered$sample_id,
    value > 0
  ) %>%
  pivot_wider(names_from = "sample_id", values_from = "value", values_fill = 0)

# What have we lost?
# Remember this loss includes controls and the removal of resequenced samples
# In total we lost 114,693 or 75.6% of OTUs from the out of the box OUT table 
total_OTUs - nrow(OTUs_filtered_1)
100 - (nrow(OTUs_filtered_1) / total_OTUs * 100)
# Corresponding to 32,541,646 or 18.4% of the total abundance
total_abundance - sum(
  OTUs_filtered_1 %>%
    column_to_rownames(var = "OTU_ID")
)
100 - (sum(
  OTUs_filtered_1 %>% 
    column_to_rownames(var = "OTU_ID")) / total_abundance * 100
)

# Compared to what we lost at the start:
total_OTUs - nrow(OTUs)
100 - (nrow(OTUs) / total_OTUs * 100)
# Corresponding to:
total_abundance - sum(OTUs)
100 - (sum(OTUs) / total_abundance * 100)

# Visualise the relationship between sequencing depth and richness
seq_depth %>%
  ggplot(aes(abundance, richness)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "top", label.x.npc = "left") +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = ",")
  ) +
  labs(y = "Observed OTU richness", x = "Logartium of sequencing depth") +
  theme_bw()

# Save the filtered OTU table
fwrite( 
  OTUs_filtered_1,
  "../data/OTUs/OTUs_abundance_filtered_50_1.csv"
)

# Clear all objects from the environment
rm(list = ls())

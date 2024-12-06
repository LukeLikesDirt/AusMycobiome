
# Script: Calculate ECM occurrences in Antarctica

# Required packages
require(data.table)
require(tidyverse)

# ECM filter based on FungalTraits
ecm_filter <- fread(
  "data//fungal_traits/fungal_traits.csv"
) %>%
  filter(
    primary_lifestyle == "ectomycorrhizal"
  ) %>%
  select(genus) %>%
  unique(.)

### 1. My ECM fungi ############################################################

# Filter my dataset to ECM fungi
my_ecm <- fread(
  "../output/otu_table.txt"
  ) %>%
  pivot_longer(
    cols = -OTU_ID,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(
    abundance > 0
  ) %>%
  select(OTU_ID, sample_id) %>%
  left_join(
    fread("../output/taxonomy.txt"),
    by = "OTU_ID"
  ) %>%
  select(OTU_ID, sample_id, genus, species) %>%
  filter(
    genus %in% ecm_filter$genus
    ) %>%
  left_join(
    fread("../output/sample_metadata.txt"),
    by = "sample_id"
  ) %>%
  select(OTU_ID, sample_id, genus, latitude, longitude) %>%
  glimpse(.)

# How many ECM occurrences are there: 11,833
my_ecm %>%
  nrow(.)

# How many ECM OTUs are there: 1960
my_ecm %>%
  select(OTU_ID) %>%
  unique(.) %>%
  nrow(.)

# How many samples have ECM fungi: 1,586
my_ecm %>%
  select(sample_id) %>%
  unique(.) %>%
  nrow(.)

# How many genera are there: 93
my_ecm %>%
  select(genus) %>%
  unique(.) %>%
  nrow(.)

# How many ECM OTUs are there in Antarctica: 0
my_ecm %>%
  filter(
    latitude < -65
  ) %>%
  select(OTU_ID) %>%
  unique(.) %>%
  print(n = Inf) %>%
  nrow(.)

### 2. ALA ECM fungi ###########################################################

ala_ecm <- fread(
  "data/australian_microbiome_ALA/observations.csv"
  ) %>%
  select(
    genus, species, taxonConceptID, verbatimDepth, occurrenceID, recordID,
    longitude = decimalLongitude, latitude = decimalLatitude
  ) %>%
  filter(
    grepl("FungiITS1", occurrenceID) & 
      verbatimDepth == "0 cm" &
      genus %in% ecm_filter$genus) %>%
  glimpse(.)

# How many sites are there on ALA: 1,723
fread(
  "data/australian_microbiome_ALA/observations.csv"
) %>%
  select(
    decimalLongitude, decimalLatitude
  ) %>%
  unique(.) %>%
  nrow(.)

# How many plots have ECM: 1,503
ala_ecm %>%
  select(longitude, latitude) %>%
  unique(.) %>%
  nrow(.)

# Proportion of plots with ECM fungi in ALA:
1503/1723 * 100

# How many genera are there: 89
ala_ecm %>%
  select(genus) %>%
  unique(.) %>%
  nrow(.)

# How many ECM occurrences are there in Antarctica: 818
ala_ecm %>%
  filter(
    latitude < -65
  ) %>%
  nrow(.)

# How many ECM genera in Antarctica: 38
ala_ecm %>%
  filter(
    latitude < -65
  ) %>%
  select(genus) %>%
  unique(.) %>%
  nrow(.)

# How many occurrences of Cortinarius are there: 102
ala_ecm %>%
  filter(
    latitude < -65 &
    genus == "Cortinarius"
  ) %>%
  nrow(.)

# List of genera in Antarctica
ala_ecm %>%
  filter(
    latitude < -65
  ) %>%
  group_by(genus) %>%
  summarise(
    n = n()
  ) %>%
  arrange(desc(n)) %>%
  print(n = Inf)
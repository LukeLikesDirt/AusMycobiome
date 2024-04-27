
# Script: Calculate ECM occurrences in Antarctica

# Required packages
require(data.table)
require(tidyverse)

# ECM filter based on FungalTraits
ecm_filter <- fread(
  "data/bioinformatics/fungal_traits/fungal_traits.csv"
) %>%
  filter(
    primary_lifestyle == "ectomycorrhizal"
  ) %>%
  select(genus) %>%
  unique(.)

### 1 My ECM fungi ############################################################

# Filter my dataset to ECM fungi
my_ecm <- fread(
  "output/OTUs.csv"
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
    fread("output/taxonomy.csv"),
    by = "OTU_ID"
  ) %>%
  select(OTU_ID, sample_id, genus, species) %>%
  filter(
    genus %in% ecm_filter$genus
  ) %>%
  left_join(
    fread("output/sample_metadata.csv"),
    by = "sample_id"
  ) %>%
  select(OTU_ID, sample_id, genus, latitude, longitude) %>%
  glimpse(.)

# How many ECM occurrences are there: 11,197
my_ecm %>%
  nrow(.)

# How many ECM OTUs are there: 2079
my_ecm %>%
  select(OTU_ID) %>%
  unique(.) %>%
  nrow(.)

# How many samples have ECM fungi: 1411
my_ecm %>%
  select(sample_id) %>%
  unique(.) %>%
  nrow(.)

# How many genera are there: 75
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
  nrow(.)

### 2 ALA ECM fungi ###########################################################

ala_ecm <- fread(
  "data/technical_validation/australian_microbiome_ALA/observations.csv"
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
  "data/technical_validation/australian_microbiome_ALA/observations.csv"
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

# Proportion of plots with ECM fungi in ALA: 0.0003
1503/1723 * 100

# How many genera are there: 89
ala_ecm %>%
  select(genus) %>%
  unique(.) %>%
  nrow(.)

# How many ECM occurances are there in Antarctica: 818
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

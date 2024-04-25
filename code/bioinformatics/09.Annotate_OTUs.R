
# Script:   Annotate BLAST hits to OTUs using taxa-specific similarity
#           thresholds and a 66.67 consensus cutoff across ther top three hits
# Author:   Luke Florence
# Date:     1st March 2024
#
# Contents:
#   (1) Consensus estimates of BLAST hits
#   (2) Annotate OTUs based similarity thresholds
#
# Note: For both processes I subset the taxa tables from species to phylum. I 
# apply the thresholds from low to high levels, and only the taxa that do not 
# meet the criteria at a given level are assessed on the subsequent level. I 
# then I rejoin the taxa tables subsets at the end of each process. There's
# probably a way to do this process more effectively, but this is the best I can
# do for now.
#
# Required functions
require(data.table)
require(ggpubr)
require(patchwork)
require(tidyverse)

### (1) Consensus estimates of BLAST hits ######################################

# I apply 66.7% consensus cutoff across the top three hits in case the top hit
# is achieved by chance.

# Required functions
source("code/bioinformatics/functions.R")

# Quality and abundance filtered taxa table
taxa <- fread(
  "data/bioinformatics/09.Taxonomy/taxa_abundance_quality_filtered.csv"
  ) %>%
  # Generate a unique ID for each hit for the final inner join
  mutate(
    hit_ID = paste(OTU_ID, row_number(), sep = "_")
  ) %>%
  glimpse(.)

# How many OTUs do I have: 34,132
n_OTUs <- unique(taxa$OTU_ID) %>% 
  length(.) %>%
  print(.)

#### (1a) Species consensus ######################################################

species_consensus_1 <- taxa %>%
  select(OTU_ID, hit_ID, species) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[species == first(species)]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(species_consensus_1$OTU_ID))

# How many OTUs have been annotated to species: 29,570 or 86.6%
nrow(species_consensus_1)
nrow(species_consensus_1)/n_OTUs * 100

# Check second best hits
species_consensus_2 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus_1$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, species) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[species == species[2]]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(2) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(species_consensus_2$OTU_ID))

# Based on the second best hit: 3,332 or 9.8 %
nrow(species_consensus_2)
nrow(species_consensus_2)/n_OTUs * 100

# Join the species taxa tables
# In addition, I will add fungal traits based on genera
species_consensus <- inner_join_colnames(
  species_consensus_1, species_consensus_2
) %>%
  left_join(
    taxa, by = c("OTU_ID", "hit_ID", "species")
  ) %>%
  select(-hit_ID) %>%
  glimpse(.)

# Check for duplicated OTU_IDs
any(duplicated(species_consensus$OTU_ID))

# Number of OTUs annotated to species: 32,902 or 96.4 %
nrow(species_consensus)
nrow(species_consensus)/n_OTUs * 100

#### (1b) Genus consensus ######################################################

genus_consensus_1 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, genus) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[genus == first(genus)]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(genus_consensus_1$OTU_ID))

# How many OTUs have been annotated to genus: 359 or 1.1%
nrow(genus_consensus_1)
nrow(genus_consensus_1)/n_OTUs * 100

# Check second best hits
genus_consensus_2 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus_1$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, genus) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[genus == genus[2]]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(2) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(genus_consensus_2$OTU_ID))

# Based on the second best hit: 87 or 0.3 %
nrow(genus_consensus_2)
nrow(genus_consensus_2)/n_OTUs * 100

# Join the genus taxa tables
# In addition, I will add fungal traits based on genera
genus_consensus <- inner_join_colnames(
  genus_consensus_1, genus_consensus_2
) %>%
  left_join(
    taxa, by = c("OTU_ID", "hit_ID", "genus")
  ) %>%
  select(-hit_ID) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", genus),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", genus),
      paste0(genus, "_sp_unidentified")
    )
  ) %>%
  glimpse(.)

# Check for duplicated OTU_IDs
any(duplicated(genus_consensus$OTU_ID))

# Number of OTUs annotated to genus: 446 or 1.3 %
nrow(genus_consensus)
nrow(genus_consensus)/n_OTUs * 100

#### (1c) Family consensus ######################################################

family_consensus_1 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, family) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[family == first(family)]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(family_consensus_1$OTU_ID))

# How many OTUs have been annotated to family: 191 or 0.56%
nrow(family_consensus_1)
nrow(family_consensus_1)/n_OTUs * 100

# Check second best hits
family_consensus_2 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus_1$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, family) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[family == family[2]]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(2) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(family_consensus_2$OTU_ID))

# Based on the second best hit: 72 or 0.2 %
nrow(family_consensus_2)
nrow(family_consensus_2)/n_OTUs * 100

# Join the family taxa tables
family_consensus <- inner_join_colnames(
  family_consensus_1, family_consensus_2
) %>%
  left_join(
    taxa, by = c("OTU_ID", "hit_ID", "family")
  ) %>%
  select(-hit_ID) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", family),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", family),
      paste0(family, "_sp_unidentified")
    ),
    genus = ifelse(
      grepl("_Incertae_sedis$", family),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", family),
      paste0(family, "_gen_unidentified")
    )
  ) %>%
  glimpse(.)

# Check for duplicated OTU_IDs
any(duplicated(family_consensus$OTU_ID))

# Number of OTUs annotated to family: 263 or 0.8 %
nrow(family_consensus)
nrow(family_consensus)/n_OTUs * 100

#### (1d) Order consensus ######################################################

order_consensus_1 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
  ) %>%
  select(OTU_ID, hit_ID, order) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[order == first(order)]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(order_consensus_1$OTU_ID))

# How many OTUs have been annotated to order: 143 or 0.4%
nrow(order_consensus_1)
nrow(order_consensus_1)/n_OTUs * 100

# Check second best hits
order_consensus_2 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
    !OTU_ID %in% order_consensus_1$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, order) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[order == order[2]]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(2) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(order_consensus_2$OTU_ID))

# Based on the second best hit: 52 or 0.15 %
nrow(order_consensus_2)
nrow(order_consensus_2)/n_OTUs * 100

# Join the order taxa tables
order_consensus <- inner_join_colnames(
  order_consensus_1, order_consensus_2
) %>%
  left_join(
    taxa, by = c("OTU_ID", "hit_ID", "order")
  ) %>%
  select(-hit_ID) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", order),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", order),
      paste0(order, "_sp_unidentified")
    ),
    genus = ifelse(
      grepl("_Incertae_sedis$", order),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", order),
      paste0(order, "_gen_unidentified")
    ),
    family = ifelse(
      grepl("_Incertae_sedis$", order),
      sub("_\\w{3}_Incertae_sedis$", "_fam_unidentified", order),
      paste0(order, "_fam_unidentified")
    )
  ) %>%
  glimpse(.)

# Check for duplicated OTU_IDs
any(duplicated(order_consensus$OTU_ID))

# Number of OTUs annotated to order: 195 or 0.6 %
nrow(order_consensus)
nrow(order_consensus)/n_OTUs * 100

#### (1e) Class consensus ######################################################

class_consensus_1 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
    !OTU_ID %in% order_consensus$OTU_ID,
  ) %>%
  select(OTU_ID, hit_ID, class) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[class == first(class)]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(class_consensus_1$OTU_ID))

# How many OTUs have been annotated to class: 48 or 0.1%
nrow(class_consensus_1)
nrow(class_consensus_1)/n_OTUs * 100

# Check second best hits
class_consensus_2 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
    !OTU_ID %in% order_consensus$OTU_ID,
    !OTU_ID %in% class_consensus_1$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, class) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[class == class[2]]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(2) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(class_consensus_2$OTU_ID))

# Based on the second best hit: 20 or 0.1 %
nrow(class_consensus_2)
nrow(class_consensus_2)/n_OTUs * 100

# Join the class taxa tables
class_consensus <- inner_join_colnames(
  class_consensus_1, class_consensus_2
) %>%
  left_join(
    taxa, by = c("OTU_ID", "hit_ID", "class")
  ) %>%
  select(-hit_ID) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", class),
      paste0(class, "_sp_unidentified")
    ),
    genus = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", class),
      paste0(class, "_gen_unidentified")
    ),
    family = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_fam_unidentified", class),
      paste0(class, "_fam_unidentified")
    ),
    order = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_ord_unidentified", class),
      paste0(class, "_ord_unidentified")
    )
  ) %>%
  glimpse(.)

# Check for duplicated OTU_IDs
any(duplicated(class_consensus$OTU_ID))

# Number of OTUs annotated to class: 68 or 0.2 %
nrow(class_consensus)
nrow(class_consensus)/n_OTUs * 100

#### (1f) Phylum consensus #####################################################

phylum_consensus_1 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
    !OTU_ID %in% order_consensus$OTU_ID,
    !OTU_ID %in% class_consensus$OTU_ID,
  ) %>%
  select(OTU_ID, hit_ID, phylum) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[phylum == first(phylum)]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(phylum_consensus_1$OTU_ID))

# How many OTUs have been annotated to phylum: 50 or 0.1%
nrow(phylum_consensus_1)
nrow(phylum_consensus_1)/n_OTUs * 100

# Check second best hits
phylum_consensus_2 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
    !OTU_ID %in% order_consensus$OTU_ID,
    !OTU_ID %in% class_consensus$OTU_ID,
    !OTU_ID %in% phylum_consensus_1$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, phylum) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[phylum == phylum[2]]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(2) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(phylum_consensus_2$OTU_ID))

# Based on the second best hit: 19 or 0.1 %
nrow(phylum_consensus_2)
nrow(phylum_consensus_2)/n_OTUs * 100

# Join the phylum taxa tables
phylum_consensus <- inner_join_colnames(phylum_consensus_1, phylum_consensus_2) %>%
  left_join(taxa, by = c("OTU_ID", "hit_ID", "phylum")) %>%
  select(-hit_ID) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", phylum),
      paste0(phylum, "_sp_unidentified")
    ),
    genus = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", phylum),
      paste0(phylum, "_gen_unidentified")
    ),
    family = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_fam_unidentified", phylum),
      paste0(phylum, "_fam_unidentified")
    ),
    order = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_ord_unidentified", phylum),
      paste0(phylum, "_ord_unidentified")
    ),
    class = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_cls_unidentified", phylum),
      paste0(phylum, "_cls_unidentified")
    )
  ) %>%
  glimpse(.)

# Check for duplicated OTU_IDs
any(duplicated(phylum_consensus$OTU_ID))

# Number of OTUs annotated to phylum: 69 or 0.2 %
nrow(phylum_consensus)
nrow(phylum_consensus)/n_OTUs * 100

#### (1g) Kingdom consensus ######################################################

kingdom_consensus_1 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
    !OTU_ID %in% order_consensus$OTU_ID,
    !OTU_ID %in% class_consensus$OTU_ID,
    !OTU_ID %in% phylum_consensus$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, kingdom) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[kingdom == first(kingdom)]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(kingdom_consensus_1$OTU_ID))

# How many OTUs have been annotated to kingdom: 80 or 0.2%
nrow(kingdom_consensus_1)
nrow(kingdom_consensus_1)/n_OTUs * 100

# Check second best hits
kingdom_consensus_2 <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
    !OTU_ID %in% order_consensus$OTU_ID,
    !OTU_ID %in% class_consensus$OTU_ID,
    !OTU_ID %in% phylum_consensus$OTU_ID,
    !OTU_ID %in% kingdom_consensus_1$OTU_ID
  ) %>%
  select(OTU_ID, hit_ID, kingdom) %>%
  mutate(hit_count = 1) %>%
  group_by(OTU_ID) %>%
  # Based on best three hits
  slice(1:3) %>%
  mutate(
    total_hits = sum(hit_count),
    same_taxa_hits = sum(hit_count[kingdom == kingdom[2]]),
    taxa_consensus = same_taxa_hits / total_hits * 100
  ) %>%
  slice(2) %>%
  ungroup() %>%
  filter(taxa_consensus > 66) %>%
  select(-c(hit_count, total_hits, same_taxa_hits)) %>%
  glimpse()

# See if I have retained any duplicated OTU_IDs
any(duplicated(kingdom_consensus_2$OTU_ID))

# Based on the second best hit: 42 or 0.1 %
nrow(kingdom_consensus_2)
nrow(kingdom_consensus_2)/n_OTUs * 100

# Join the kingdom taxa tables
kingdom_consensus <- inner_join_colnames(
  kingdom_consensus_1, kingdom_consensus_2
) %>%
  left_join(
    taxa, by = c("OTU_ID", "hit_ID", "kingdom")
  ) %>%
  select(-hit_ID) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", kingdom),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", kingdom),
      paste0(kingdom, "_sp_unidentified")
    ),
    genus = ifelse(
      grepl("_Incertae_sedis$", kingdom),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", kingdom),
      paste0(kingdom, "_gen_unidentified")
    ),
    family = ifelse(
      grepl("_Incertae_sedis$", kingdom),
      sub("_\\w{3}_Incertae_sedis$", "_fam_unidentified", kingdom),
      paste0(kingdom, "_fam_unidentified")
    ),
    order = ifelse(
      grepl("_Incertae_sedis$", kingdom),
      sub("_\\w{3}_Incertae_sedis$", "_ord_unidentified", kingdom),
      paste0(kingdom, "_ord_unidentified")
    ),
    class = ifelse(
      grepl("_Incertae_sedis$", kingdom),
      sub("_\\w{3}_Incertae_sedis$", "_cls_unidentified", kingdom),
      paste0(kingdom, "_cls_unidentified")
    ),
    phylum = ifelse(
      grepl("_Incertae_sedis$", kingdom),
      sub("_\\w{3}_Incertae_sedis$", "_phy_unidentified", kingdom),
      paste0(kingdom, "_phy_unidentified")
    )
  )

# Check for duplicated OTU_IDs
any(duplicated(kingdom_consensus$OTU_ID))

# Number of OTUs annotated to kingdom: 122 or 0.4%
nrow(kingdom_consensus)
nrow(kingdom_consensus)/n_OTUs * 100

#### (1h) Unidentified taxa ####################################################

# Check n_OTUs that will remain unidentified at level kingdom
unidentified_OTUs <- taxa %>%
  filter(
    !OTU_ID %in% species_consensus$OTU_ID,
    !OTU_ID %in% genus_consensus$OTU_ID,
    !OTU_ID %in% family_consensus$OTU_ID,
    !OTU_ID %in% order_consensus$OTU_ID,
    !OTU_ID %in% class_consensus$OTU_ID,
    !OTU_ID %in% phylum_consensus$OTU_ID,
    !OTU_ID %in% kingdom_consensus$OTU_ID
  ) %>%
  group_by(OTU_ID) %>%
  slice(1) %>%
  ungroup()

# Check for duplicated OTU_IDs
any(duplicated(unidentified_OTUs$OTU_ID))

# Number of OTUs that remain unidentified: 67 or 0.2 %
nrow(unidentified_OTUs)
nrow(unidentified_OTUs)/n_OTUs * 100

# Check to see if the numbers match up:
# The number of OTUs minus the number of unidentified OTUs: 34,065
n_OTUs - nrow(unidentified_OTUs)
# Sum of identified and unassigned OTUs: 34,065
sum(
  nrow(species_consensus),
  nrow(genus_consensus), nrow(family_consensus), nrow(order_consensus),
  nrow(class_consensus), nrow(phylum_consensus), nrow(kingdom_consensus)
)
# Looks good!

#### (1i) Join the taxa tables and save #######################################

# Fungi taxa table
fungi_consensus <- inner_join_colnames(
  species_consensus, genus_consensus
) %>%
  inner_join_colnames(family_consensus) %>%
  inner_join_colnames(order_consensus) %>%
  inner_join_colnames(class_consensus) %>%
  inner_join_colnames(phylum_consensus) %>%
  inner_join_colnames(kingdom_consensus) %>%
  filter(kingdom == "Fungi") %>%
  # I'll mutate "GS01_phy_Incertae_sedis" to GS01 for interoperability
  mutate(phylum = case_when(
    phylum == "GS01_phy_Incertae_sedis" ~ "GS01",
    TRUE ~ as.character(phylum)
  )) %>%
  select(
    OTU_ID, phylum, class, order, family, genus, species, species_hypothesis,
    similarity, coverage, evalue
  ) %>%
  mutate(across(everything(), ~ifelse(
    . == "", if(is.integer(.)) NA_integer_ else NA_character_, .)
  )) %>%
  glimpse(.)

# How many fungal OTUs do we have: 31,458 or 92.3% of all OTUs were fungi
nrow(fungi_consensus)
nrow(fungi_consensus) / 
  (sum(
    nrow(species_consensus),
    nrow(genus_consensus), nrow(family_consensus), nrow(order_consensus),
    nrow(class_consensus), nrow(phylum_consensus), nrow(kingdom_consensus)
  )) * 100

# Save the taxa table with one hit per OTU_ID
fwrite(
  fungi_consensus,
  "data/bioinformatics/09.Taxonomy/fungi_consensus.csv"
)

# Clear the environment
rm(list = ls())

### (2) Annotate OTUs based similarity thresholds #############################

# Here I annotate taxa according to similarity thresholds defined by Tedersoo et
# al. (2022). Best practices in metabarcoding of fungi. Molecular Ecology, 
# 31(10), 2769-2795.

# Similarity thresholds for genera are defined at level phylum or order. The 
# Phylum level thresholds are set at 70% similarity. For intermediate levels
# (class, order, family), the thresholds are set at intermediate points between
# the phylum level thresholds genus level thresholds.

# Required functions
source("code/bioinformatics/functions.R")

# Consensus hit taxa table
taxa <- fread("data/bioinformatics/09.Taxonomy/fungi_consensus.csv")

# How many fungal OTUs do I have: 31,511
n_OTUs <- unique(taxa$OTU_ID) %>% 
  length(.) %>%
  print(.)

# Read in similarity thresholds

# Order threshold exceptions
order_level_threshold <- read.csv(
  'data/bioinformatics/similarity_thresholds/taxon_thresholds.csv'
) %>%
  filter(
    taxon_rank == "order"
  ) %>%
  select(
    order, class_threshold, order_threshold, family_threshold, genus_threshold
  ) %>%
  glimpse(.)

# Phylum thresholds
phylum_level_threshold <- read.csv(
  'data/bioinformatics/similarity_thresholds/taxon_thresholds.csv'
) %>%
  filter(
    taxon_rank == "phylum"
  ) %>%
  select(
    phylum, class_threshold, order_threshold, family_threshold, genus_threshold
  ) %>%
  glimpse(.)

#### (2a) Annotate species ####################################################

# Apply species threshold
species_threshold <- taxa %>%
  filter(
    coverage >= 95, similarity >= 97
  ) %>%
  as_tibble() %>%
  glimpse(.)

# See if I have retained any duplicated OTU_IDs
any(duplicated(species_threshold$OTU_ID))

# How many OTUs have been annotated to genus: 12,926 or 41.1%
nrow(species_threshold)
nrow(species_threshold)/n_OTUs * 100

#### (2b) Annotate genus ######################################################

# Apply order level thresholds
genus_threshold_1 <- taxa %>%
  filter(
    coverage >= 90,
    !OTU_ID %in% species_threshold$OTU_ID,
    order %in% order_level_threshold$order
  ) %>%
  left_join(order_level_threshold, by = "order") %>%
  filter(similarity >= genus_threshold) %>%
  select(
    -c(class_threshold, order_threshold, family_threshold, genus_threshold)
  ) %>%
  as_tibble(.) %>%
  glimpse(.)

# Filter by phylum level thresholds
genus_threshold_2 <- taxa %>%
  filter(
    coverage >= 90,
    !OTU_ID %in% species_threshold$OTU_ID,
    !order %in% order_level_threshold$order,
    phylum %in% phylum_level_threshold$phylum
  ) %>%
  left_join(phylum_level_threshold, by = "phylum") %>%
  filter(similarity >= genus_threshold) %>%
  select(
    -c(class_threshold, order_threshold, family_threshold, genus_threshold)
  ) %>%
  as_tibble(.) %>%
  glimpse(.)

# Join filtered taxa tables
genus_threshold <- inner_join_colnames(
  genus_threshold_1,
  genus_threshold_2
  ) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", genus),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", genus),
      ifelse(
        grepl("_unidentified$", genus),
        sub("_\\w{3}_unidentified$", "_sp_unidentified", genus),
        ifelse(
          grepl("_sp_unidentified$", species),
          species,
          paste0(genus, "_sp_unidentified")
        )
      )
    )
  ) %>%
  glimpse(.)

# Check the mutate function: Looks good!
genus_threshold %>%
  select(genus, species) %>%
  print(n = 500)

# See if I have retained any duplicated OTU_IDs
any(duplicated(genus_threshold$OTU_ID))

# How many OTUs have been annotated to genus: 10,372 or 32.9%
nrow(genus_threshold)
nrow(genus_threshold)/n_OTUs * 100

#### (2c) Annotate family ######################################################

# Apply order level thresholds
family_threshold_1 <- taxa %>%
  filter(
    coverage >= 85,
    !OTU_ID %in% species_threshold$OTU_ID,
    !OTU_ID %in% genus_threshold$OTU_ID,
    order %in% order_level_threshold$order
  ) %>%
  left_join(order_level_threshold, by = "order") %>%
  filter(similarity >= family_threshold) %>%
  select(
    -c(class_threshold, order_threshold, family_threshold, genus_threshold)
  ) %>%
  as_tibble(.) %>%
  glimpse(.)

# Filter by phylum level thresholds
family_threshold_2 <- taxa %>%
  filter(
    coverage >= 85,
    !OTU_ID %in% species_threshold$OTU_ID,
    !OTU_ID %in% genus_threshold$OTU_ID,
    !order %in% order_level_threshold$order,
    phylum %in% phylum_level_threshold$phylum
  ) %>%
  left_join(phylum_level_threshold, by = "phylum") %>%
  filter(similarity >= family_threshold) %>%
  select(
    -c(class_threshold, order_threshold, family_threshold, genus_threshold)
  ) %>%
  as_tibble(.) %>%
  glimpse(.)

# Join filtered taxa tables
family_threshold <- inner_join_colnames(
  family_threshold_1,
  family_threshold_2
  ) %>%
  mutate(
    genus = ifelse(
      grepl("_Incertae_sedis$", family),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", family),
      ifelse(
        grepl("_unidentified$", family),
        sub("_\\w{3}_unidentified$", "_gen_unidentified", family),
        ifelse(
          grepl("_gen_unidentified$", genus),
          genus,
          paste0(family, "_gen_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", family),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", family),
      ifelse(
        grepl("_unidentified$", family),
        sub("_\\w{3}_unidentified$", "_sp_unidentified", family),
        ifelse(
          grepl("_sp_unidentified$", species),
          species,
          paste0(family, "_sp_unidentified")
        )
      )
    )
  ) %>%
  glimpse(.)

# Check the mutate function: Looks good!
family_threshold %>%
  select(family, genus, species) %>%
  print(n = 500)

# See if I have retained any duplicated OTU_IDs
any(duplicated(family_threshold$OTU_ID))

# How many OTUs have been annotated to family: 3,468 or 11.0%
nrow(family_threshold)
nrow(family_threshold)/n_OTUs * 100

#### (2d) Annotate order ######################################################

# Apply order level thresholds
order_threshold_1 <- taxa %>%
  filter(
    coverage >= 80,
    !OTU_ID %in% species_threshold$OTU_ID,
    !OTU_ID %in% genus_threshold$OTU_ID,
    !OTU_ID %in% family_threshold$OTU_ID,
    order %in% order_level_threshold$order
  ) %>%
  left_join(order_level_threshold, by = "order") %>%
  filter(similarity >= order_threshold) %>%
  select(
    -c(class_threshold, order_threshold, family_threshold, genus_threshold)
  ) %>%
  as_tibble(.) %>%
  glimpse(.)

# Filter by phylum level thresholds
order_threshold_2 <- taxa %>%
  filter(
    coverage >= 80,
    !OTU_ID %in% species_threshold$OTU_ID,
    !OTU_ID %in% genus_threshold$OTU_ID,
    !OTU_ID %in% family_threshold$OTU_ID,
    !order %in% order_level_threshold$order,
    phylum %in% phylum_level_threshold$phylum
  ) %>%
  left_join(phylum_level_threshold, by = "phylum") %>%
  filter(similarity >= order_threshold) %>%
  select(
    -c(class_threshold, order_threshold, family_threshold, genus_threshold)
  ) %>%
  as_tibble(.) %>%
  glimpse(.)

# Join filtered taxa tables
order_threshold <- inner_join_colnames(
  order_threshold_1,
  order_threshold_2
  ) %>%
  mutate(
    family = ifelse(
      grepl("_Incertae_sedis$", order),
      sub("_\\w{3}_Incertae_sedis$", "_fam_unidentified", order),
      ifelse(
        grepl("_unidentified$", order),
        sub("_\\w{3}_unidentified$", "_fam_unidentified", order),
        ifelse(
          grepl("_fam_unidentified$", family),
          family,
          paste0(order, "_fam_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    genus = ifelse(
      grepl("_Incertae_sedis$", order),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", order),
      ifelse(
        grepl("_unidentified$", order),
        sub("_\\w{3}_unidentified$", "_gen_unidentified", order),
        ifelse(
          grepl("_gen_unidentified$", genus),
          genus,
          paste0(order, "_gen_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", order),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", order),
      ifelse(
        grepl("_unidentified$", order),
        sub("_\\w{3}_unidentified$", "_sp_unidentified", order),
        ifelse(
          grepl("_sp_unidentified$", species),
          species,
          paste0(order, "_sp_unidentified")
        )
      )
    )
  ) %>%
  glimpse(.)

# Check the mutate function: Looks good!
order_threshold %>%
  select(order, family, genus, species) %>%
  print(n = 500)

# See if I have retained any duplicated OTU_IDs
any(duplicated(order_threshold$OTU_ID))

# How many OTUs have been annotated to order: 1,868 or 5.9%
nrow(order_threshold)
nrow(order_threshold)/n_OTUs * 100

#### (2e) Annotate class ######################################################

# Apply order level thresholds
class_threshold <- taxa %>%
  filter(
    coverage >= 75,
    !OTU_ID %in% species_threshold$OTU_ID,
    !OTU_ID %in% genus_threshold$OTU_ID,
    !OTU_ID %in% family_threshold$OTU_ID,
    !OTU_ID %in% order_threshold$OTU_ID,
    phylum %in% phylum_level_threshold$phylum
  ) %>%
  left_join(phylum_level_threshold, by = "phylum") %>%
  filter(similarity >= class_threshold) %>%
  select(
    -c(class_threshold, order_threshold, family_threshold, genus_threshold)
  ) %>%
  as_tibble(.) %>%
  mutate(
    order = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_ord_unidentified", class),
      ifelse(
        grepl("_unidentified$", class),
        sub("_\\w{3}_unidentified$", "_ord_unidentified", class),
        ifelse(
          grepl("_ord_unidentified$", order),
          order,
          paste0(class, "_ord_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    family = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_fam_unidentified", class),
      ifelse(
        grepl("_unidentified$", class),
        sub("_\\w{3}_unidentified$", "_fam_unidentified", class),
        ifelse(
          grepl("_fam_unidentified$", family),
          family,
          paste0(class, "_fam_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    genus = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", class),
      ifelse(
        grepl("_unidentified$", class),
        sub("_\\w{3}_unidentified$", "_gen_unidentified", class),
        ifelse(
          grepl("_gen_unidentified$", genus),
          genus,
          paste0(class, "_gen_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", class),
      ifelse(
        grepl("_unidentified$", class),
        sub("_\\w{3}_unidentified$", "_sp_unidentified", class),
        ifelse(
          grepl("_sp_unidentified$", species),
          species,
          paste0(class, "_sp_unidentified")
        )
      )
    )
  ) %>%
  glimpse(.)

# Check the mutate function: Looks good!
class_threshold %>%
  select(class, order, family, genus, species) %>%
  print(n = 500)

# See if I have retained any duplicated OTU_IDs
any(duplicated(class_threshold$OTU_ID))

# How many OTUs have been annotated to class: 1,107 or 3.5%
nrow(class_threshold)
nrow(class_threshold)/n_OTUs * 100

#### (2f) Annotate phyla ######################################################

# Apply order level thresholds
phylum_threshold <- taxa %>%
  filter(
    !OTU_ID %in% species_threshold$OTU_ID,
    !OTU_ID %in% genus_threshold$OTU_ID,
    !OTU_ID %in% family_threshold$OTU_ID,
    !OTU_ID %in% order_threshold$OTU_ID,
    !OTU_ID %in% class_threshold$OTU_ID
  ) %>%
  as_tibble(.) %>%
  mutate(
    class = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_cls_unidentified", phylum),
      ifelse(
        grepl("_unidentified$", phylum),
        sub("_\\w{3}_unidentified$", "_cls_unidentified", phylum),
        ifelse(
          grepl("_cls_unidentified$", class),
          class,
          paste0(phylum, "_cls_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    order = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_ord_unidentified", phylum),
      ifelse(
        grepl("_unidentified$", phylum),
        sub("_\\w{3}_unidentified$", "_ord_unidentified", phylum),
        ifelse(
          grepl("_ord_unidentified$", order),
          order,
          paste0(phylum, "_ord_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    family = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_fam_unidentified", phylum),
      ifelse(
        grepl("_unidentified$", phylum),
        sub("_\\w{3}_unidentified$", "_fam_unidentified", phylum),
        ifelse(
          grepl("_fam_unidentified$", family),
          family,
          paste0(phylum, "_fam_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    genus = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", phylum),
      ifelse(
        grepl("_unidentified$", phylum),
        sub("_\\w{3}_unidentified$", "_gen_unidentified", phylum),
        ifelse(
          grepl("_gen_unidentified$", genus),
          genus,
          paste0(phylum, "_gen_unidentified")
        )
      )
    )
  ) %>%
  mutate(
    species = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_sp_unidentified", phylum),
      ifelse(
        grepl("_unidentified$", phylum),
        sub("_\\w{3}_unidentified$", "_sp_unidentified", phylum),
        ifelse(
          grepl("_sp_unidentified$", species),
          species,
          paste0(phylum, "_sp_unidentified")
        )
      )
    )
  ) %>%
  glimpse(.)

# Check the mutate function: Looks good!
phylum_threshold %>%
  select(phylum, class, order, family, genus, species) %>%
  print(n = 500)

# See if I have retained any duplicated OTU_IDs
any(duplicated(phylum_threshold$OTU_ID))

# How many OTUs have been annotated to phylum: 1,744 or 5.5%
nrow(phylum_threshold)
nrow(phylum_threshold)/n_OTUs * 100

#### (2g) Annotate functional traits ###################################################

# Join the data frames
fungi_threshold <- inner_join_colnames(
  species_threshold, genus_threshold
) %>%
  inner_join_colnames(family_threshold) %>%
  inner_join_colnames(order_threshold) %>%
  inner_join_colnames(class_threshold) %>%
  inner_join_colnames(phylum_threshold) %>%
  # Add trait information
  left_join(., fread(
    "data/bioinformatics/fungal_traits/fungal_traits.csv"
  ), by = "genus") %>%
  # Glomeromycota traits are defined at the phylum-level
  mutate(
    primary_lifestyle = case_when(
      phylum == "Glomeromycota" ~ "arbuscular_mycorrhizal",
      TRUE ~ primary_lifestyle
    ),
    secondary_lifestyle = case_when(
      phylum == "Glomeromycota" ~ "root-associated",
      TRUE ~ primary_lifestyle
    ),
    plant_endophytic_capability = case_when(
      phylum == "Glomeromycota" ~ "root_endophyte",
      TRUE ~ primary_lifestyle
    ),
    growth_form = case_when(
      phylum == "Glomeromycota" ~ "filamentous_mycelium",
      TRUE ~ primary_lifestyle
    ),
    fruitbody_type = case_when(
      phylum == "Glomeromycota" ~ "none",
      TRUE ~ primary_lifestyle
    ),
    hymenium_type = case_when(
      phylum == "Glomeromycota" ~ "none",
      TRUE ~ primary_lifestyle
    )
  ) %>%
  mutate(across(everything(), ~ifelse(
    . == "", if(is.integer(.)) NA_integer_ else NA_character_, .)
  )) %>%
  glimpse(.)

# Check the mutate function:
fungi_threshold %>% 
  filter(phylum == "Glomeromycota") %>%
  select(genus, primary_lifestyle, secondary_lifestyle, plant_endophytic_capability, growth_form, fruitbody_type, hymenium_type) %>%
  print(n = 500)

# See if I have any duplicated OTU_IDs
any(duplicated(fungi_threshold$OTU_ID))

# Check to see if the numbers match up:
# The number of OTUs minus the number of unidentified OTUs: 31,458
n_OTUs
# Sum of identified and unassigned OTUs: 31,458
fungi_threshold %>%
  nrow(.)
# Looks good!

### (3) Save the final taxa and OTU tables ####################################

# Filter the OTU table to fungi and remove run IDs from sample names
OTUs_abundance_filtered <- fread(
  "data/bioinformatics/08.OTUs/OTUs_abundance_filtered.csv"
  ) %>%
  # Filter OTUs to taxa
  filter(
    OTU_ID %in% fungi_threshold$OTU_ID
  ) %>%
  # Remove library run suffix from sample names
  pivot_longer(
    cols = -OTU_ID,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  mutate(
    sample = sub("_.*", "", sample)
  ) %>%
  pivot_wider(
    names_from = "sample",
    values_from = "abundance"
  ) %>%
  # Removed because < 5000 reads: s13482, s13514, s7051, s13556, s19412
  # Removed because these are identified later identified as richness outlines
  # and are misassignments from mainland Australia: s13619, s13561
  select(-c(s13482, s13514, s7051, s13556, s19412, s13619, s13561, s13607)) %>%
  # Remove any redundant OTUs
  filter(
    rowSums(select(., -OTU_ID)) > 0
  )

# How many samples and OTUs in the final OTU table: 31,455 fungal OTUs in 2,093
# samples
nrow(OTUs_abundance_filtered)
ncol(OTUs_abundance_filtered) - 1

# Calculate OTU abundances
OTUs_abundances <- OTUs_abundance_filtered %>%
  pivot_longer(
    -OTU_ID,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  select(OTU_ID, abundance) %>%
  group_by(OTU_ID) %>%
  summarise(
    abundance = sum(abundance)
  ) %>%
  ungroup() %>%
  print()

# Save the filtered OTU table
fwrite(
  OTUs_abundance_filtered,
  "output/OTUs.csv"
  )

# Save the taxa table with one hit per OTU_ID
fungi_threshold %>%
  inner_join(OTUs_abundances, by = "OTU_ID") %>%
  select(
    OTU_ID, phylum, class, order, family, genus, species, species_hypothesis,
    abundance, everything(.)
  ) %>%
  arrange(desc(abundance)) %>%
  fwrite(
  "output/taxonomy.csv"
  )

# Clear the environment
rm(list = ls())

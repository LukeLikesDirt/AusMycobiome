
cutoff_file_path <- "./dnabarcoder/ITS1_cutoffs.txt"
require(Biostrings)
require(data.table)
require(tidyverse)
source("read_blast.R")

# NOTE: I have mannually moved "Global Fungi" cutoffs to the cutoff file to make
# the cutoff file complete for ASVs identifie to fungi. I need to automate this
# in the DNAbarcoder pipeline in later iterations

# Read in the dnabarcoder cutoff file
cutoff_file <- fread(cutoff_file_path) %>%
  select(rank = Rank, taxa = Dataset, cutoff = "cut-off")

# (1) Cutoff files for each rank ###############################################

### (1a) Species cutoffs ####
species_cutoffs <- fread("../data/taxonomy/BLAST_ASVs_species.txt") %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  mutate(rank = "species") %>%
  select(OTU_ID, rank, genus, family, order, class, phylum, kingdom) %>%
  pivot_longer(
    -c(OTU_ID, rank), 
    names_to = "level",
    values_to = "taxa"
  ) %>%
  left_join(
    .,
    cutoff_file,
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  group_by(OTU_ID) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(OTU_ID, cutoff)

# Check for NA values
species_cutoffs %>%
  filter(
    is.na(cutoff)
  )

### (1b) Genus cutoffs ####
genus_cutoffs <- fread("../data/taxonomy/BLAST_ASVs_genus.txt") %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  mutate(rank = "genus") %>%
  select(OTU_ID, rank, family, order, class, phylum, kingdom) %>%
  pivot_longer(
    -c(OTU_ID, rank), 
    names_to = "level",
    values_to = "taxa"
  ) %>%
  left_join(
    .,
    cutoff_file,
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  group_by(OTU_ID) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(OTU_ID, cutoff)

# Check for NA values
genus_cutoffs %>%
  filter(
    is.na(cutoff)
  )

### (1c) Family cutoffs ####
family_cutoffs <- fread("../data/taxonomy/BLAST_ASVs_genus.txt") %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  mutate(rank = "family") %>%
  select(OTU_ID, rank, order, class, phylum, kingdom) %>%
  pivot_longer(
    -c(OTU_ID, rank), 
    names_to = "level",
    values_to = "taxa"
  ) %>%
  #select(-level) %>%
  unique(.) %>%
  left_join(
    .,
    cutoff_file,
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  group_by(OTU_ID) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(OTU_ID, cutoff)

# Check for NA values
family_cutoffs %>%
  filter(
    is.na(cutoff)
  )

### (1d) Order cutoffs ####
order_cutoffs <- fread("../data/taxonomy/BLAST_ASVs_all.txt") %>%
  # Format the BLAST data
  read_blast(., minlen = 50) %>%
  mutate(rank = "order") %>%
  select(OTU_ID, rank, class, phylum, kingdom) %>%
  pivot_longer(
    -c(OTU_ID, rank), 
    names_to = "level",
    values_to = "taxa"
  ) %>%
  #select(-level) %>%
  unique(.) %>%
  left_join(
    .,
    cutoff_file,
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  group_by(OTU_ID) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(OTU_ID, cutoff)

# Check for NA values
order_cutoffs %>%
  filter(
    is.na(cutoff)
  )

### (1e) Class cutoffs ####
class_cutoffs <- fread("../data/taxonomy/BLAST_ASVs_all.txt") %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  mutate(rank = "class") %>%
  select(OTU_ID, rank, phylum, kingdom) %>%
  pivot_longer(
    -c(OTU_ID, rank), 
    names_to = "level",
    values_to = "taxa"
  ) %>%
  #select(-level) %>%
  unique(.) %>%
  left_join(
    .,
    cutoff_file,
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  group_by(OTU_ID) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(OTU_ID, cutoff)

# Check for NA values
class_cutoffs %>%
  filter(
    is.na(cutoff)
  )

phylum_cutoff <- cutoff_file %>%
  filter(rank == "phylum") %>%
  pull(cutoff)
print(phylum_cutoff)

# (2) Assign taxonomy ##########################################################

### (2a) Annotate species ####
species_hits <- fread("../data/taxonomy/BLAST_ASVs_species.txt") %>%
  # Format the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter using coverage cutoffs
  filter(
    subject_coverage > 0.9 & query_coverage > 0.9
  ) %>%
  # Filter based on species cutoffs
  left_join(
    .,
    species_cutoffs,
    by = "OTU_ID"
  ) %>%
  filter(
    score >= cutoff
  ) %>%
  # Filter based on consensus cutoffs
  group_by(OTU_ID, species) %>%
  mutate(
    hits = n()
  ) %>%
  ungroup() %>%
  group_by(OTU_ID) %>%
  mutate(
    total_hits = n(),
    consensus = hits / total_hits
  ) %>%
  filter(
    consensus > 0.66
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(
    kingdom %in% c("Fungi", "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi") 
  ) %>%
  mutate(
    rank = "species"
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus))

### (2b) Annotate genus ####
genus_hits <- fread("../data/taxonomy/BLAST_ASVs_genus.txt") %>%
  # Remove the OTUs identified to species
  filter(
    !OTU_ID %in% species_hits$OTU_ID
  ) %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter using coverage cutoffs
  # filter(
  #   subject_coverage > 0.90 & query_coverage > 0.90
  # ) %>%
  # Filter based on genus cutoffs
  left_join(
    .,
    genus_cutoffs,
    by = "OTU_ID"
  ) %>%
  filter(
    score >= cutoff
  ) %>%
  # Filter based on consensus cutoffs
  group_by(OTU_ID, genus) %>%
  mutate(
    hits = n()
  ) %>%
  ungroup() %>%
  group_by(OTU_ID) %>%
  mutate(
    total_hits = n(),
    consensus = hits / total_hits
  ) %>%
  filter(
    consensus > 0.66
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(
    kingdom %in% c("Fungi", "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi") 
  ) %>%
  mutate(
    rank = "genus",
    species = "unidentified"
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus))

### (2c) Annotate family ####
family_hits <- fread("../data/taxonomy/BLAST_ASVs_all.txt") %>%
  # Remove the OTUs identified to species
  filter(
    !OTU_ID %in% species_hits$OTU_ID,
    !OTU_ID %in% genus_hits$OTU_ID
  ) %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter based on family cutoffs
  left_join(
    .,
    family_cutoffs,
    by = "OTU_ID"
  ) %>%
  filter(
    score >= cutoff
  ) %>%
  # Filter based on consensus cutoffs
  group_by(OTU_ID, family) %>%
  mutate(
    hits = n()
  ) %>%
  ungroup() %>%
  group_by(OTU_ID) %>%
  mutate(
    total_hits = n(),
    consensus = hits / total_hits
  ) %>%
  filter(
    consensus > 0.66
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(
    kingdom %in% c("Fungi", "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi") 
  ) %>%
  mutate(
    rank = "family",
    species = "unidentified",
    genus = "unidentified"
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2d) Annotate order ####
order_hits <- fread("../data/taxonomy/BLAST_ASVs_all.txt") %>%
  # Remove the OTUs identified to species
  filter(
    !OTU_ID %in% species_hits$OTU_ID,
    !OTU_ID %in% genus_hits$OTU_ID,
    !OTU_ID %in% family_hits$OTU_ID
  ) %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter based on order cutoffs
  left_join(
    .,
    order_cutoffs,
    by = "OTU_ID"
  ) %>%
  filter(
    score >= cutoff
  ) %>%
  # Filter based on consensus cutoffs
  group_by(OTU_ID, order) %>%
  mutate(
    hits = n()
  ) %>%
  ungroup() %>%
  group_by(OTU_ID) %>%
  mutate(
    total_hits = n(),
    consensus = hits / total_hits
  ) %>%
  filter(
    consensus > 0.66
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(
    kingdom %in% c("Fungi", "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi") 
  ) %>%
  mutate(
    rank = "order",
    species = "unidentified",
    genus = "unidentified",
    family = "unidentified"
    ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2e) Annotate class ####
class_hits <- fread("../data/taxonomy/BLAST_ASVs_all.txt") %>%
  # Remove the OTUs identified to species, genus, family, and order
  filter(
    !OTU_ID %in% species_hits$OTU_ID,
    !OTU_ID %in% genus_hits$OTU_ID,
    !OTU_ID %in% family_hits$OTU_ID,
    !OTU_ID %in% order_hits$OTU_ID
  ) %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter based on order cutoffs
  left_join(
    .,
    class_cutoffs,
    by = "OTU_ID"
  ) %>%
  filter(
    score >= cutoff
  ) %>%
  # Filter based on consensus cutoffs
  group_by(OTU_ID, class) %>%
  mutate(
    hits = n()
  ) %>%
  ungroup() %>%
  group_by(OTU_ID) %>%
  mutate(
    total_hits = n(),
    consensus = hits / total_hits
  ) %>%
  filter(
    consensus > 0.66
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(
    kingdom %in% c("Fungi", "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi") 
  ) %>%
  mutate(
    rank = "class",
    species = "unidentified",
    genus = "unidentified",
    family = "unidentified",
    order = "unidentified"
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2f) Annotate phylum ####
phylum_hits <- fread("../data/taxonomy/BLAST_ASVs_all.txt") %>%
  # Remove the OTUs identified to species, genus, family, order, and class
  filter(
    !OTU_ID %in% species_hits$OTU_ID,
    !OTU_ID %in% genus_hits$OTU_ID,
    !OTU_ID %in% family_hits$OTU_ID,
    !OTU_ID %in% order_hits$OTU_ID,
    !OTU_ID %in% class_hits$OTU_ID
  ) %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter based on phylum cutoffs
  mutate(
    cutoff = phylum_cutoff
  ) %>%
  filter(
    score >= cutoff
  ) %>%
  # Filter based on consensus cutoffs
  group_by(OTU_ID, phylum) %>%
  mutate(
    hits = n()
  ) %>%
  ungroup() %>%
  group_by(OTU_ID) %>%
  mutate(
    total_hits = n(),
    consensus = hits / total_hits
  ) %>%
  filter(
    consensus > 0.66
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(
    kingdom %in% c("Fungi", "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi") 
  ) %>%
  mutate(
    rank = "phylum",
    species = "unidentified",
    genus = "unidentified",
    family = "unidentified",
    order = "unidentified",
    class = "unidentified"
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2g) Annotate kingdom ####
kingdom_hits <- fread("../data/taxonomy/BLAST_ASVs_all.txt") %>%
  # Remove the OTUs identified to species, genus, family, order, class, and phylum
  filter(
    !OTU_ID %in% species_hits$OTU_ID,
    !OTU_ID %in% genus_hits$OTU_ID,
    !OTU_ID %in% family_hits$OTU_ID,
    !OTU_ID %in% order_hits$OTU_ID,
    !OTU_ID %in% class_hits$OTU_ID,
    !OTU_ID %in% phylum_hits$OTU_ID
  ) %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter based on phylum cutoff
  mutate(
    cutoff = phylum_cutoff
  ) %>%
  filter(
    score >= cutoff
  ) %>%
  # Filter based on consensus cutoffs
    group_by(OTU_ID, kingdom) %>%
  mutate(
    hits = n()
  ) %>%
  ungroup() %>%
  group_by(OTU_ID) %>%
  mutate(
    total_hits = n(),
    consensus = hits / total_hits
  ) %>%
  filter(
    consensus > 0.66
  ) %>%
  slice(1) %>%
  ungroup() %>%
  filter(
    kingdom %in% c("Fungi", "dikarya_fungi", "terrestrial_basal_fungi", "zoosporic_basal_fungi") 
  ) %>%
  mutate(
    rank = "kingdom",
    species = "unidentified",
    genus = "unidentified",
    family = "unidentified",
    order = "unidentified",
    class = "unidentified",
    phylum = "unidentified"
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2h) Bind the data ####

n_OTUs <- rbind(
  species_hits,
  genus_hits,
  family_hits,
  order_hits,
  class_hits,
  phylum_hits,
  kingdom_hits
) %>% nrow()

taxonomy <- rbind(
  species_hits,
  genus_hits,
  family_hits,
  order_hits,
  class_hits,
  phylum_hits,
  kingdom_hits
  ) %>%
  mutate(
    # Fungi to unidentified across lower ranks
    phylum = case_when(
      kingdom == "Fungi" ~ "unidentified",
      TRUE ~ phylum
    ),
    class = case_when(
      kingdom == "Fungi" ~ "unidentified",
      TRUE ~ class
    ),
    order = case_when(
      kingdom == "Fungi" ~ "unidentified",
      TRUE ~ order
    ),
    family = case_when(
      kingdom == "Fungi" ~ "unidentified",
      TRUE ~ family
    ),
    genus = case_when(
      kingdom == "Fungi" ~ "unidentified",
      TRUE ~ genus
    ),
    species = case_when(
      kingdom == "Fungi" ~ "unidentified",
      TRUE ~ species
    ),
    # Update the rank of classification
    rank = case_when(
      kingdom != "unidentified" & phylum == "unidentified" & class == "unidentified" & order == "unidentified" & family == "unidentified" & genus == "unidentified" & species == "unidentified" ~ "kingdom",
      phylum != "unidentified" & class == "unidentified" & order == "unidentified" & family == "unidentified" & genus == "unidentified" & species == "unidentified" ~ "phylum",
      class != "unidentified" & order == "unidentified" & family == "unidentified" & genus == "unidentified" & species == "unidentified" ~ "class",
      order != "unidentified" & family == "unidentified" & genus == "unidentified" & species == "unidentified" ~ "order",
      family != "unidentified" & genus == "unidentified" & species == "unidentified" ~ "family",
      genus != "unidentified" & species == "unidentified" ~ "genus",
      species != "unidentified" ~ "species",
      TRUE ~ "NA"
    )
  )
filter(taxonomy, rank == "NA") %>% nrow()

# (3) Save data ####################################################################

### (3a) Taxonomy ####
taxonomy %>% 
  select(-species_hypothesis) %>%
fwrite(., "../data/annotated_asvs/taxonomy.txt", sep = "\t")

### (3b) OTUs ####
fread("../data/OTUs/ASVs_abundance_filtered.txt") %>%
  mutate(
    OTU_ID = str_extract(OTU_ID, "^[^;]+")
  ) %>%
  filter(
    OTU_ID %in% taxonomy$OTU_ID
  ) %>%
  pivot_longer(
    -OTU_ID,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  filter(
    abundance > 0
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = abundance,
    values_fill = list(abundance = 0)
  ) %>%
  fwrite(., "../data/annotated_asvs/asv_table.txt", sep = "\t")
    
### (3c) fasta ####

fasta <- readDNAStringSet("../data/OTUs/ASVs_abundance_filtered.fasta")

# Filter fasta to OTU IDs in taxonomy
fasta_tibble <- tibble(
  seq = as.character(fasta), 
  names = names(fasta)
  ) %>%
  mutate(
    names = str_extract(names, "^[^;]+")
  ) %>%
  filter(
    names %in% taxonomy$OTU_ID
  )

# Convert filtered tibble back to a DNAStringSet 
filtered_fasta <- DNAStringSet(fasta_tibble$seq) 
names(filtered_fasta) <- fasta_tibble$names

# Save the fasta
writeXStringSet(filtered_fasta, "../data/annotated_asvs/sequences.fasta")


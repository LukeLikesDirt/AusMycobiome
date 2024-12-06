
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

# Species cutoffs
species_cutoffs <- fread("../data/taxonomy/BLAST_OTUs_species.txt") %>%
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

genus_cutoffs <- fread("../data/taxonomy/BLAST_OTUs_genus.txt") %>%
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

family_cutoffs <- fread("../data/taxonomy/BLAST_OTUs_genus.txt") %>%
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

order_cutoffs <- fread("../data/taxonomy/BLAST_OTUs_all.txt") %>%
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

class_cutoffs <- fread("../data/taxonomy/BLAST_OTUs_all.txt") %>%
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

phylum_cutoff <- cutoff_file %>%
  filter(rank == "phylum") %>%
  pull(cutoff)

# (2) Assign taxonomy ##########################################################

### (2a) Annotate species ####
species_hits <- fread("../data/taxonomy/BLAST_OTUs_species.txt") %>%
  # Format the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter using coverage cutoffs
  filter(
    subject_coverage > 0.95 & query_coverage > 0.95
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
    kingdom = "Fungi",
    rank = "species"
    )%>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse()

### (2b) Annotate genus ####
genus_hits <- fread("../data/taxonomy/BLAST_OTUs_genus.txt") %>%
  # Remove the OTUs identified to species
  filter(
    !OTU_ID %in% species_hits$OTU_ID
  ) %>%
  # Formt the BLAST data
  read_blast(., minlen = 50) %>%
  # Filter using coverage cutoffs
  filter(
    subject_coverage > 0.90 & query_coverage > 0.90
  ) %>%
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
    kingdom = "Fungi",
    rank = "genus",
    species = ifelse(
      grepl("_Incertae_sedis$", genus),
      sub("_\\w{3}_Incertae_sedis$", "_sp", genus),
      paste0(genus, "_sp")
    )
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse()

### (2c) Annotate family ####
family_hits <- fread("../data/taxonomy/BLAST_OTUs_all.txt") %>%
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
    kingdom = "Fungi",
    rank = "family",
    species = ifelse(
      grepl("_Incertae_sedis$", family),
      sub("_\\w{3}_Incertae_sedis$", "_sp", family),
      paste0(family, "_sp")
    ),
    genus = ifelse(
      grepl("_Incertae_sedis$", family),
      sub("_\\w{3}_Incertae_sedis$", "_gen_unidentified", family),
      paste0(family, "_gen_unidentified")
    )
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2d) Annotate order ####
order_hits <- fread("../data/taxonomy/BLAST_OTUs_all.txt") %>%
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
    kingdom = "Fungi",
    rank = "order",
    species = ifelse(
      grepl("_Incertae_sedis$", order),
      sub("_\\w{3}_Incertae_sedis$", "_sp", order),
      paste0(order, "_sp")
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
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2e) Annotate class ####
class_hits <- fread("../data/taxonomy/BLAST_OTUs_all.txt") %>%
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
    kingdom = "Fungi",
    rank = "class",
    species = ifelse(
      grepl("_Incertae_sedis$", class),
      sub("_\\w{3}_Incertae_sedis$", "_sp", class),
      paste0(class, "_sp")
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
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2f) Annotate phylum ####
phylum_hits <- fread("../data/taxonomy/BLAST_OTUs_all.txt") %>%
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
    kingdom = "Fungi",
    rank = "phylum",
    species = ifelse(
      grepl("_Incertae_sedis$", phylum),
      sub("_\\w{3}_Incertae_sedis$", "_sp", phylum),
      paste0(phylum, "_sp")
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
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2g) Annotate kingdom ####
kingdom_hits <- fread("../data/taxonomy/BLAST_OTUs_all.txt") %>%
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
    kingdom = "Fungi",
    rank = "kingdom",
    species = ifelse(
      grepl("_Incertae_sedis$", kingdom),
      sub("_\\w{3}_Incertae_sedis$", "_sp", kingdom),
      paste0(kingdom, "_sp")
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
  ) %>%
  select(-c(subject_coverage, query_coverage, hits, total_hits, consensus)) %>%
  glimpse(.)

### (2h) Annotate functional traits ####
taxonomy <- rbind(
  species_hits,
  genus_hits,
  family_hits,
  order_hits,
  class_hits,
  phylum_hits,
  kingdom_hits
) %>%
  # Add trait information
  left_join(., fread(
    "../../technical_validation/data/fungal_traits/fungal_traits.csv"
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
  ))

# (3) Save data ####################################################################

### (3a) Taxonomy ####
fwrite(taxonomy, "../../output/taxonomy.txt", sep = "\t")

### (3b) OTUs ####
fread("../data/OTUs/OTUs_abundance_filtered.txt") %>%
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
  fwrite(., "../../output/otu_table.txt", sep = "\t")
    
### (3c) fasta ####

fasta <- readDNAStringSet("../data/OTUs/OTUs_abundance_filtered.fasta")

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
writeXStringSet(filtered_fasta, "../../output/sequences.fasta")

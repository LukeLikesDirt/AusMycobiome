suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

# Mutate GS01_phy_Incertae_sedis to GS01 to identify GS01 as a valid phylum considering it has at least been placed within fungi in UNITE
fread("data/unite2024ITS1.unique.classification") %>%
  mutate(
    phylum = str_replace(phylum, "GS01_phy_Incertae_sedis", "GS01"),
    class = str_replace(class, "GS01_phy_Incertae_sedis", "GS01"),
    order = str_replace(order, "GS01_phy_Incertae_sedis", "GS01"),
    family = str_replace(family, "GS01_phy_Incertae_sedis", "GS01"),
    genus = str_replace(genus, "GS01_phy_Incertae_sedis", "GS01"),
  ) %>%
  fwrite("data/unite2024ITS1.unique.classification", sep = "\t")

# Read in the sequences and classifications
classifications <- fread("data/unite2024ITS1.unique.classification")
sequences <- readDNAStringSet("data/unite2024ITS1.unique.fasta")

# Generate a list of IDs for each subkingdom or subkingdom group
dikarya_fungi_ids <- classifications %>%
  filter(phylum %in% c("Ascomycota", "Basidiomycota", "Entorrhizomycota")) %>%
  pull(id) %>%
  unique(.)

terrestrial_basal_fungi_ids <- classifications %>%
  filter(phylum %in% c("Glomeromycota", "Mortierellomycota", "Mucoromycota", "Calcarisporiellomycota", "Basidiobolomycota", "Entomophthoromycota", "Kickxellomycota", "Zoopagomycota")) %>%
  pull(id) %>%
  unique(.)

zoosporic_basal_fungi_ids <- classifications %>%
  filter(phylum %in% c("Aphelidiomycota", "Blastocladiomycota", "Chytridiomycota", "Monoblepharomycota", "Neocallimastigomycota", "Olpidiomycota", "Rozellomycota", "Sanchytriomycota", "GS01")) %>%
  pull(id) %>%
  unique(.)

Fungi_incertae_sedis_ids <- classifications %>%
  filter(
    phylum == "Fungi_phy_Incertae_sedis"
  ) %>%
  pull(id) %>%
  unique(.)

# Check that all IDs are accounted for
unclassified_rows <- classifications %>%
  filter(
    !id %in% c(dikarya_fungi_ids, terrestrial_basal_fungi_ids, zoosporic_basal_fungi_ids),
    phylum != "Fungi_phy_Incertae_sedis"
  )

# Issue a warning if there are unclassified rows
if (nrow(unclassified_rows) > 0) {
  warning(paste("There are", nrow(unclassified_rows), "sequences that have not been defined by a group!"))
} else {
  message("All classified sequences are accounted for!")
}

# Generate classification files for predictions and taxanomic assignments
classifications %>%
  mutate(
    kingdom = case_when(
      phylum %in% c("Ascomycota", "Basidiomycota", "Entorrhizomycota") ~ "dikarya_fungi",
      phylum %in% c("Glomeromycota", "Mortierellomycota", "Mucoromycota", "Calcarisporiellomycota", "Basidiobolomycota", "Entomophthoromycota", "Kickxellomycota", "Zoopagomycota") ~ "terrestrial_basal_fungi",
      phylum %in% c("Aphelidiomycota", "Blastocladiomycota", "Chytridiomycota", "Monoblepharomycota", "Neocallimastigomycota", "Olpidiomycota", "Rozellomycota", "Sanchytriomycota", "GS01") ~ "zoosporic_basal_fungi",
      TRUE ~ kingdom
    )
  ) %>%
  fwrite("data/unite2024ITS1.groups.classification", sep = "\t")
classifications <- fread("data/unite2024ITS1.groups.classification")

# Write out the sequences and classifications for groups
classifications %>%
  filter(id %in% dikarya_fungi_ids) %>%
  fwrite("data/dikarya_fungi.classification", sep = "\t")
writeXStringSet(sequences[dikarya_fungi_ids], "data/dikarya_fungi.fasta")

classifications %>%
  filter(id %in% terrestrial_basal_fungi_ids) %>%
  fwrite("data/terrestrial_basal_fungi.classification", sep = "\t")
writeXStringSet(sequences[terrestrial_basal_fungi_ids], "data/terrestrial_basal_fungi.fasta")

classifications %>%
  filter(id %in% zoosporic_basal_fungi_ids) %>%
  fwrite("data/zoosporic_basal_fungi.classification", sep = "\t")
writeXStringSet(sequences[zoosporic_basal_fungi_ids], "data/zoosporic_basal_fungi.fasta")
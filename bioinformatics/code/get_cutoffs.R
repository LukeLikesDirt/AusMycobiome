
# Function to get species cutoff
get_species_cutoff <- function(taxa_file, unique_taxa_cutoffs, cutoff_file, superranks_species) {
  species_cutoffs_all <- taxa_file %>%
    select(rank, all_of(superranks_species)) %>%
    filter(rank == "species") %>%
    pivot_longer(-rank, values_to = "taxa") %>%
    left_join(unique_taxa_cutoffs %>% filter(rank == "species"), by = c("rank", "taxa")) %>%
    filter(!is.na(cutoff)) %>%
    select(-rank, rank = name) %>%
    unique()
  
  species_cutoffs <- list()
  for (superrank in superranks_species) {
    species_cutoffs[[superrank]] <- species_cutoffs_all %>% filter(rank == superrank) %>% select(-rank)
  }
  
  taxa_file %>%
    mutate(
      species_cutoff = case_when(
        genus %in% species_cutoffs[["genus"]][["taxa"]] ~ species_cutoffs[["genus"]][["cutoff"]][match(genus, species_cutoffs[["genus"]][["taxa"]])],
        family %in% species_cutoffs[["family"]][["taxa"]] ~ species_cutoffs[["family"]][["cutoff"]][match(family, species_cutoffs[["family"]][["taxa"]])],
        order %in% species_cutoffs[["order"]][["taxa"]] ~ species_cutoffs[["order"]][["cutoff"]][match(order, species_cutoffs[["order"]][["taxa"]])],
        class %in% species_cutoffs[["class"]][["taxa"]] ~ species_cutoffs[["class"]][["cutoff"]][match(class, species_cutoffs[["class"]][["taxa"]])],
        phylum %in% species_cutoffs[["phylum"]][["taxa"]] ~ species_cutoffs[["phylum"]][["cutoff"]][match(phylum, species_cutoffs[["phylum"]][["taxa"]])],
        kingdom %in% species_cutoffs[["kingdom"]][["taxa"]] ~ species_cutoffs[["kingdom"]][["cutoff"]][match(kingdom, species_cutoffs[["kingdom"]][["taxa"]])],
        TRUE ~ ifelse(nrow(cutoff_file %>% filter(rank == "species" & taxa == "All")) > 0,
                      cutoff_file %>% filter(rank == "species" & taxa == "All") %>% pull(cutoff),
                      NA_real_)
      )
    )
}

# Function to get genus cutoff
get_genus_cutoff <- function(taxa_file, unique_taxa_cutoffs, cutoff_file, superranks_genus) {
  genus_cutoffs_all <- taxa_file %>%
    select(rank, all_of(superranks_genus), kingdom) %>%
    filter(rank == "genus") %>%
    pivot_longer(-rank, values_to = "taxa") %>%
    left_join(unique_taxa_cutoffs %>% filter(rank == "genus"), by = c("rank", "taxa")) %>%
    filter(!is.na(cutoff)) %>%
    select(-rank, rank = name) %>%
    unique()
  
  genus_cutoffs <- list()
  for (superrank in superranks_genus) {
    genus_cutoffs[[superrank]] <- genus_cutoffs_all %>% filter(rank == superrank) %>% select(-rank)
  }
  
  taxa_file %>%
    mutate(
      genus_cutoff = case_when(
        family %in% genus_cutoffs[["family"]][["taxa"]] ~ genus_cutoffs[["family"]][["cutoff"]][match(family, genus_cutoffs[["family"]][["taxa"]])],
        order %in% genus_cutoffs[["order"]][["taxa"]] ~ genus_cutoffs[["order"]][["cutoff"]][match(order, genus_cutoffs[["order"]][["taxa"]])],
        class %in% genus_cutoffs[["class"]][["taxa"]] ~ genus_cutoffs[["class"]][["cutoff"]][match(class, genus_cutoffs[["class"]][["taxa"]])],
        phylum %in% genus_cutoffs[["phylum"]][["taxa"]] ~ genus_cutoffs[["phylum"]][["cutoff"]][match(phylum, genus_cutoffs[["phylum"]][["taxa"]])],
        kingdom %in% genus_cutoffs[["kingdom"]][["taxa"]] ~ genus_cutoffs[["kingdom"]][["cutoff"]][match(kingdom, genus_cutoffs[["kingdom"]][["taxa"]])],
        TRUE ~ ifelse(nrow(cutoff_file %>% filter(rank == "genus" & taxa == "All")) > 0,
                      cutoff_file %>% filter(rank == "genus" & taxa == "All") %>% pull(cutoff),
                      NA_real_)
      )
    )
}

# Function to get family cutoff
get_family_cutoff <- function(taxa_file, unique_taxa_cutoffs, cutoff_file, superranks_family) {
  family_cutoffs_all <- taxa_file %>%
    select(rank, all_of(superranks_family), kingdom) %>%
    filter(rank == "family") %>%
    pivot_longer(-rank, values_to = "taxa") %>%
    left_join(unique_taxa_cutoffs %>% filter(rank == "family"), by = c("rank", "taxa")) %>%
    filter(!is.na(cutoff)) %>%
    select(-rank, rank = name) %>%
    unique()
  
  family_cutoffs <- list()
  for (superrank in superranks_family) {
    family_cutoffs[[superrank]] <- family_cutoffs_all %>% filter(rank == superrank) %>% select(-rank)
  }
  
  taxa_file %>%
    mutate(
      family_cutoff = case_when(
        order %in% family_cutoffs[["order"]][["taxa"]] ~ family_cutoffs[["order"]][["cutoff"]][match(order, family_cutoffs[["order"]][["taxa"]])],
        class %in% family_cutoffs[["class"]][["taxa"]] ~ family_cutoffs[["class"]][["cutoff"]][match(class, family_cutoffs[["class"]][["taxa"]])],
        phylum %in% family_cutoffs[["phylum"]][["taxa"]] ~ family_cutoffs[["phylum"]][["cutoff"]][match(phylum, family_cutoffs[["phylum"]][["taxa"]])],
        kingdom %in% family_cutoffs[["kingdom"]][["taxa"]] ~ family_cutoffs[["kingdom"]][["cutoff"]][match(kingdom, family_cutoffs[["kingdom"]][["taxa"]])],
        TRUE ~ ifelse(nrow(cutoff_file %>% filter(rank == "family" & taxa == "All")) > 0,
                      cutoff_file %>% filter(rank == "family" & taxa == "All") %>% pull(cutoff),
                      NA_real_)
      )
    )
}

# Function to get order cutoff
get_order_cutoff <- function(taxa_file, unique_taxa_cutoffs, cutoff_file, superranks_order) {
  order_cutoffs_all <- taxa_file %>%
    select(rank, all_of(superranks_order), kingdom) %>%
    filter(rank == "order") %>%
    pivot_longer(-rank, values_to = "taxa") %>%
    left_join(unique_taxa_cutoffs %>% filter(rank == "order"), by = c("rank", "taxa")) %>%
    filter(!is.na(cutoff)) %>%
    select(-rank, rank = name) %>%
    unique()
  
  order_cutoffs <- list()
  for (superrank in superranks_order) {
    order_cutoffs[[superrank]] <- order_cutoffs_all %>% filter(rank == superrank) %>% select(-rank)
  }
  
  taxa_file %>%
    mutate(
      order_cutoff = case_when(
        class %in% order_cutoffs[["class"]][["taxa"]] ~ order_cutoffs[["class"]][["cutoff"]][match(class, order_cutoffs[["class"]][["taxa"]])],
        phylum %in% order_cutoffs[["phylum"]][["taxa"]] ~ order_cutoffs[["phylum"]][["cutoff"]][match(phylum, order_cutoffs[["phylum"]][["taxa"]])],
        kingdom %in% order_cutoffs[["kingdom"]][["taxa"]] ~ order_cutoffs[["kingdom"]][["cutoff"]][match(kingdom, order_cutoffs[["kingdom"]][["taxa"]])],
        TRUE ~ ifelse(nrow(cutoff_file %>% filter(rank == "order" & taxa == "All")) > 0,
                      cutoff_file %>% filter(rank == "order" & taxa == "All") %>% pull(cutoff),
                      NA_real_)
      )
    )
}

# Function to get class cutoff
get_class_cutoff <- function(taxa_file, unique_taxa_cutoffs, cutoff_file, superranks_class) {
  class_cutoffs_all <- taxa_file %>%
    select(rank, all_of(superranks_class), kingdom) %>%
    filter(rank == "class") %>%
    pivot_longer(-rank, values_to = "taxa") %>%
    left_join(unique_taxa_cutoffs %>% filter(rank == "class"), by = c("rank", "taxa")) %>%
    filter(!is.na(cutoff)) %>%
    select(-rank, rank = name) %>%
    unique()
  
  class_cutoffs <- list()
  for (superrank in superranks_class) {
    class_cutoffs[[superrank]] <- class_cutoffs_all %>% filter(rank == superrank) %>% select(-rank)
  }
  
  taxa_file %>%
    mutate(
      class_cutoff = case_when(
        phylum %in% class_cutoffs[["phylum"]][["taxa"]] ~ class_cutoffs[["phylum"]][["cutoff"]][match(phylum, class_cutoffs[["phylum"]][["taxa"]])],
        kingdom %in% class_cutoffs[["kingdom"]][["taxa"]] ~ class_cutoffs[["kingdom"]][["cutoff"]][match(kingdom, class_cutoffs[["kingdom"]][["taxa"]])],
        TRUE ~ ifelse(nrow(cutoff_file %>% filter(rank == "class" & taxa == "All")) > 0,
                      cutoff_file %>% filter(rank == "class" & taxa == "All") %>% pull(cutoff),
                      NA_real_)
      )
    )
}

# Function to get phylum cutoff
get_phylum_cutoff <- function(taxa_file, unique_taxa_cutoffs, cutoff_file, superranks_phylum) {
  phylum_cutoffs_all <- taxa_file %>%
    select(rank, all_of(superranks_phylum), kingdom) %>%
    filter(rank == "phylum") %>%
    pivot_longer(-rank, values_to = "taxa") %>%
    left_join(unique_taxa_cutoffs %>% filter(rank == "phylum"), by = c("rank", "taxa")) %>%
    filter(!is.na(cutoff)) %>%
    select(-rank, rank = name) %>%
    unique()
  
  phylum_cutoffs <- list()
  for (superrank in superranks_phylum) {
    phylum_cutoffs[[superrank]] <- phylum_cutoffs_all %>% filter(rank == superrank) %>% select(-rank)
  }
  
  taxa_file %>%
    mutate(
      phylum_cutoff = case_when(
        kingdom %in% phylum_cutoffs[["kingdom"]][["taxa"]] ~ phylum_cutoffs[["kingdom"]][["cutoff"]][match(kingdom, phylum_cutoffs[["kingdom"]][["taxa"]])],
        TRUE ~ ifelse(nrow(cutoff_file %>% filter(rank == "phylum" & taxa == "All")) > 0,
                      cutoff_file %>% filter(rank == "phylum" & taxa == "All") %>% pull(cutoff),
                      NA_real_)
      )
    )
}

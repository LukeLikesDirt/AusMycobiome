
# Function for phylum rank
taxa_match_phylum <- function(this_superrank) {
  inner_join(
    fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
      select(-all_of(superranks_phylum)),
    fread("tmp_clusters/kingdom_clusters.txt") %>%
      select(OTU_ID, kingdom),
    by = "OTU_ID"
  )
}

# Function for class rank
taxa_match_class <- function(this_superrank) {
  inner_join(
    fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
      select(-all_of(superranks_class)),
    fread("tmp_clusters/kingdom_clusters.txt") %>%
      select(OTU_ID, kingdom),
    by = "OTU_ID"
  ) %>%
    inner_join(
      .,
      fread("tmp_clusters/phylum_clusters.txt") %>%
        select(OTU_ID, phylum),
      by = "OTU_ID"
    )
}

# Function for order rank
taxa_match_order <- function(this_superrank) {
  inner_join(
    fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
      select(-all_of(superranks_order)),
    fread("tmp_clusters/kingdom_clusters.txt") %>%
      select(OTU_ID, kingdom),
    by = "OTU_ID"
  ) %>%
    inner_join(
      .,
      fread("tmp_clusters/phylum_clusters.txt") %>%
        select(OTU_ID, phylum),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/class_clusters.txt") %>%
        select(OTU_ID, class),
      by = "OTU_ID"
    )
}

# Function for family rank
taxa_match_family <- function(this_superrank) {
  inner_join(
    fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
      select(-all_of(superranks_family)),
    fread("tmp_clusters/kingdom_clusters.txt") %>%
      select(OTU_ID, kingdom),
    by = "OTU_ID"
  ) %>%
    inner_join(
      .,
      fread("tmp_clusters/phylum_clusters.txt") %>%
        select(OTU_ID, phylum),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/class_clusters.txt") %>%
        select(OTU_ID, class),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/order_clusters.txt") %>%
        select(OTU_ID, order),
      by = "OTU_ID"
    )
}

# Function for genus rank
taxa_match_genus <- function(this_superrank) {
  inner_join(
    fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
      select(-all_of(superranks_genus)),
    fread("tmp_clusters/kingdom_clusters.txt") %>%
      select(OTU_ID, kingdom),
    by = "OTU_ID"
  ) %>%
    inner_join(
      .,
      fread("tmp_clusters/phylum_clusters.txt") %>%
        select(OTU_ID, phylum),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/class_clusters.txt") %>%
        select(OTU_ID, class),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/order_clusters.txt") %>%
        select(OTU_ID, order),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/family_clusters.txt") %>%
        select(OTU_ID, family),
      by = "OTU_ID"
    )
}

# Function for species rank
taxa_match_species <- function(this_superrank) {
  inner_join(
    fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
      select(-all_of(superranks_species)),
    fread("tmp_clusters/kingdom_clusters.txt") %>%
      select(OTU_ID, kingdom),
    by = "OTU_ID"
  ) %>%
    inner_join(
      .,
      fread("tmp_clusters/phylum_clusters.txt") %>%
        select(OTU_ID, phylum),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/class_clusters.txt") %>%
        select(OTU_ID, class),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/order_clusters.txt") %>%
        select(OTU_ID, order),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/family_clusters.txt") %>%
        select(OTU_ID, family),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/genus_clusters.txt") %>%
        select(OTU_ID, genus),
      by = "OTU_ID"
    )
}


merge_clusters <- function() {
  inner_join(
    fread("tmp_clusters/species_clusters.txt") %>%
      select(-all_of(superranks_species)),
    fread("tmp_clusters/kingdom_clusters.txt") %>%
      select(OTU_ID, kingdom),
    by = "OTU_ID"
  ) %>%
    inner_join(
      .,
      fread("tmp_clusters/phylum_clusters.txt") %>%
        select(OTU_ID, phylum),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/class_clusters.txt") %>%
        select(OTU_ID, class),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/order_clusters.txt") %>%
        select(OTU_ID, order),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/family_clusters.txt") %>%
        select(OTU_ID, family),
      by = "OTU_ID"
    ) %>%
    inner_join(
      .,
      fread("tmp_clusters/genus_clusters.txt") %>%
        select(OTU_ID, genus),
      by = "OTU_ID"
    ) %>%
    # Rename pseudo taxa
    group_by(phylum) %>%
    mutate(
      phylum = case_when(
        str_detect(phylum, "_pseudo_") ~ paste0(str_extract(phylum, "^[^_]+"), "_pseudo_phylum_", sprintf("%04d", cur_group_id())),
        TRUE ~ phylum
      )) %>%
    ungroup() %>%
    group_by(class) %>%
    mutate(
      class = case_when(
        str_detect(class, "_pseudo_") ~ paste0(str_extract(class, "^[^_]+"), "_pseudo_class_", sprintf("%04d", cur_group_id())),
        TRUE ~ class
      )) %>%
    ungroup() %>%
    group_by(order) %>%
    mutate(
      order = case_when(
        str_detect(order, "_pseudo_") ~ paste0(str_extract(order, "^[^_]+"), "_pseudo_order_", sprintf("%04d", cur_group_id())),
        TRUE ~ order
      )) %>%
    ungroup() %>%
    group_by(family) %>%
    mutate(
      family = case_when(
        str_detect(family, "_pseudo_") ~ paste0(str_extract(family, "^[^_]+"), "_pseudo_family_", sprintf("%04d", cur_group_id())),
        TRUE ~ family
      )) %>%
    ungroup() %>%
    group_by(genus) %>%
    mutate(
      genus = case_when(
        str_detect(genus, "_pseudo_") ~ paste0(str_extract(genus, "^[^_]+"), "_pseudo_genus_", sprintf("%04d", cur_group_id())),
        TRUE ~ genus
      )) %>%
    ungroup() %>%
    group_by(species) %>%
    mutate(
      species = case_when(
        str_detect(species, "_pseudo_") ~ paste0(str_extract(species, "^[^_]+"), "_pseudo_species_", sprintf("%04d", cur_group_id())),
        TRUE ~ species
      )) %>%
    ungroup() %>%
    # Generate unique OTU_IDs and sum the abundances
    mutate(
      ASV_ID = OTU_ID,
      ASV_abundance = str_extract(OTU_ID, "(?<=size=)\\d+") %>% as.numeric()
    ) %>%
    group_by(species) %>%
    mutate(
      OTU_abundance = sum(ASV_abundance)
    ) %>%
    # Select the representative ASV sequence for the OTU cluster based on 
    # score and abundance (i.e. when scores are equal take the
    # ID with the highest score followed by the most abundant sequence)
    arrange(desc(score), desc(ASV_abundance)) %>%
    mutate(
      OTU_ID = first(OTU_ID)
    ) %>%
    ungroup() %>%
    arrange(desc(OTU_abundance)) %>%
    select(
      OTU_ID, ASV_ID, reference_ID, OTU_abundance, ASV_abundance,
      kingdom, phylum, class, order, family, genus, species, rank, cutoff, score
    ) %>%
    print(n = 100)
}

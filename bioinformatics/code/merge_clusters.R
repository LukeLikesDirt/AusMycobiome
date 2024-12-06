merge_clusters <- function() {
  
  asv_abundance <- fread(asv_table_path) %>%
    pivot_longer(
      -OTU_ID,
      names_to = "sample_ID",
      values_to = "abundance"
    ) %>%
    group_by(OTU_ID) %>%
    summarise(
      ASV_abundance = sum(abundance)
    )
  
  fread("tmp_clusters/species_clusters.txt") %>%
    unique() %>%
    select(-all_of(superranks_species)) %>%
    inner_join(
      fread("tmp_clusters/kingdom_clusters.txt") %>%
        unique() %>%
        select(OTU_ID, kingdom),
      by = "OTU_ID"
    ) %>%
    inner_join(
      fread("tmp_clusters/phylum_clusters.txt") %>%
        unique() %>%
        select(OTU_ID, phylum),
      by = "OTU_ID"
    ) %>%
    inner_join(
      fread("tmp_clusters/class_clusters.txt") %>%
        unique() %>%
        select(OTU_ID, class),
      by = "OTU_ID"
    ) %>%
    inner_join(
      fread("tmp_clusters/order_clusters.txt") %>%
        unique() %>%
        select(OTU_ID, order),
      by = "OTU_ID"
    ) %>%
    inner_join(
      fread("tmp_clusters/family_clusters.txt") %>%
        unique() %>%
        select(OTU_ID, family),
      by = "OTU_ID"
    ) %>%
    inner_join(
      fread("tmp_clusters/genus_clusters.txt") %>%
        unique() %>%
        select(OTU_ID, genus),
      by = "OTU_ID"
    ) %>%
    # Rename pseudo taxa for each rank
    group_by(phylum) %>%
    mutate(
      phylum = case_when(
        str_detect(phylum, "_pseudo_") ~ paste0(
          str_extract(phylum, "^[^_]+"), 
          "_pseudo_phylum_", 
          sprintf("%04d", cur_group_id())
        ),
        str_detect(phylum, "unidentified") ~ paste0(
          str_extract(kingdom, "^[^_]+"), 
          "_pseudo_phylum_", 
          sprintf("%04d", cur_group_id())
        ),
        (str_detect(phylum, "Incerate") & str_detect(genus, "_pseudo_")) |
          (str_detect(phylum, "Incerate") & str_detect(phylum, "unidentified")) ~ paste0(
            str_extract(kingdom, "^[^_]+"), 
            "_pseudo_phylum_", 
            sprintf("%04d", cur_group_id())
          ),
        TRUE ~ phylum
      )
    ) %>%
    ungroup() %>%
    group_by(class) %>%
    mutate(
      class = case_when(
        str_detect(class, "_pseudo_") ~ paste0(
          str_extract(class, "^[^_]+"), 
          "_pseudo_class_", 
          sprintf("%04d", cur_group_id())
        ),
        str_detect(class, "unidentified") ~ paste0(
          str_extract(phylum, "^[^_]+"), 
          "_pseudo_class_", 
          sprintf("%04d", cur_group_id())
        ),
        (str_detect(class, "Incerate") & str_detect(genus, "_pseudo_")) |
          (str_detect(class, "Incerate") & str_detect(class, "unidentified")) ~ paste0(
            str_extract(phylum, "^[^_]+"), 
            "_pseudo_class_", 
            sprintf("%04d", cur_group_id())
          ),
        TRUE ~ class
      )
    ) %>%
    ungroup() %>%
    group_by(order) %>%
    mutate(
      order = case_when(
        str_detect(order, "_pseudo_") ~ paste0(
          str_extract(order, "^[^_]+"), 
          "_pseudo_order_", 
          sprintf("%04d", cur_group_id())
        ),
        str_detect(order, "unidentified") ~ paste0(
          str_extract(class, "^[^_]+"), 
          "_pseudo_order_", 
          sprintf("%04d", cur_group_id())
        ),
        (str_detect(order, "Incerate") & str_detect(genus, "_pseudo_")) |
          (str_detect(order, "Incerate") & str_detect(order, "unidentified")) ~ paste0(
            str_extract(class, "^[^_]+"), 
            "_pseudo_order_", 
            sprintf("%04d", cur_group_id())
          ),
        TRUE ~ order
      )
    ) %>%
    ungroup() %>%
    group_by(family) %>%
    mutate(
      family = case_when(
        str_detect(family, "_pseudo_") ~ paste0(
          str_extract(family, "^[^_]+"), 
          "_pseudo_family_", 
          sprintf("%04d", cur_group_id())
        ),
        str_detect(family, "unidentified") ~ paste0(
          str_extract(order, "^[^_]+"), 
          "_pseudo_family_", 
          sprintf("%04d", cur_group_id())
        ),
        (str_detect(family, "Incerate") & str_detect(genus, "_pseudo_")) |
          (str_detect(family, "Incerate") & str_detect(family, "unidentified")) ~ paste0(
            str_extract(order, "^[^_]+"), 
            "_pseudo_family_", 
            sprintf("%04d", cur_group_id())
          ),
        TRUE ~ family
      )
    ) %>%
    ungroup() %>%
    group_by(genus) %>%
    mutate(
      genus = case_when(
        str_detect(genus, "_pseudo_") ~ paste0(
          str_extract(genus, "^[^_]+"), 
          "_pseudo_genus_", 
          sprintf("%04d", cur_group_id())
        ),
        str_detect(genus, "unidentified") ~ paste0(
          str_extract(family, "^[^_]+"), 
          "_pseudo_genus_", 
          sprintf("%04d", cur_group_id())
        ),
        (str_detect(genus, "Incerate") & str_detect(genus, "_pseudo_")) |
          (str_detect(genus, "Incerate") & str_detect(genus, "unidentified")) ~ paste0(
            str_extract(family, "^[^_]+"), 
            "_pseudo_genus_", 
            sprintf("%04d", cur_group_id())
          ),
        TRUE ~ genus
      )
    ) %>%
    ungroup() %>%
    group_by(species) %>%
    mutate(
      species = case_when(
        str_detect(species, "_pseudo_") ~ paste0(
          str_extract(species, "^[^_]+"), 
          "_pseudo_species_", 
          sprintf("%04d", cur_group_id())
        ),
        str_detect(species, "unidentified") ~ paste0(
          str_extract(genus, "^[^_]+"), 
          "_pseudo_species_", 
          sprintf("%04d", cur_group_id())
        ),
        (str_detect(species, "Incerate") & str_detect(species, "_pseudo_")) |
          (str_detect(species, "Incerate") & str_detect(species, "unidentified")) ~ paste0(
            str_extract(genus, "^[^_]+"), 
            "_pseudo_species_", 
            sprintf("%04d", cur_group_id())
          ),
        TRUE ~ species
      )
    ) %>%
    ungroup() %>%
    # Add the ASV abundance
    left_join(
      asv_abundance,
      by = "OTU_ID"
    ) %>%
    # Generate unique OTU_IDs and sum the abundances
    mutate(
      ASV_ID = OTU_ID
    ) %>%
    group_by(species) %>%
    mutate(
      OTU_abundance = sum(ASV_abundance),
      # Prioritise cluster core sequences over reference-based clusters.
      # Reference-based clusters have reference_IDs starting lowercase letters 
      # or numbers, while UNITE reference_IDs start with uppercase letters.
      # Assign a score of 0 to reference_IDs starting with a number to 
      # prevent them from being selected as representative sequences.
      score = case_when(
        str_detect(reference_ID, "^[0-9a-z]") ~ 0,
        TRUE ~ score
      )
    ) %>%
    # Select representative ASV sequence based on the best score, then ASV abundance
    arrange(desc(score), desc(ASV_abundance)) %>%
    mutate(
      OTU_ID = first(OTU_ID),
      reference_ID = first(reference_ID),
      score = first(score),
      cutoff = first(cutoff)
    ) %>%
    ungroup() %>%
    arrange(desc(OTU_abundance)) %>%
    # !!! UPDATE THE CLUSTERING LOGIC TO MAKE THIS STEP UNEEDED !!!
    # Update the reference_ID and score for the pseudo taxa generated form the
    # reference clusters to their orginal reference_ID and score for the most
    # abundant ASV
    select(-reference_ID, -score) %>%
    left_join(
      fread(taxa_file_path) %>%
        select(OTU_ID, reference_ID, score),
      by = "OTU_ID"
    ) %>%
    select(
      OTU_ID, ASV_ID, reference_ID, 
      kingdom, phylum, class, order, family, genus, species,
      rank, score, cutoff, abundance = OTU_abundance
    ) %>%
    print(n = 100)
  
}

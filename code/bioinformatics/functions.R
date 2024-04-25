
#### Bind rows by matching column names ########################################
inner_join_colnames <- function(df1, df2) {
  common_names <- intersect(names(df1), names(df2))
  
  # Subset and bind rows based on common column names
  if (length(common_names) > 0) {
    rows_df1 <- df1[, common_names, drop = FALSE]
    rows_df2 <- df2[, common_names, drop = FALSE]
    result <- bind_rows(rows_df1, rows_df2)
  } else {
    warning("Column names do not match! Check column names.")
    result <- NULL
  }
  
  return(result)
}

#### Remove low abundance OTUs using a specified threshold #####################

# Removes OTUs sample-wise with an abundance <= threshold

filter_low_abundance_otus <- function(otu_tab, threshold) {
  # Filter out zero values and calculate relative abundance
  otu_tab %>%
    pivot_longer(-OTU_ID, names_to = 'sample') %>%
    filter(value > 0) %>%
    group_by(sample) %>%
    mutate(rel_abund = value / sum(value) * 100) %>%
    ungroup() %>%
    
    # Set values below the threshold to zero
    mutate(value = replace(value, rel_abund <= threshold, 0)) %>%
    
    # Remove rows where all values are zero 
    select(-rel_abund) %>%
    filter(value > 0) %>%
    
    # Reshape the table
    pivot_wider(names_from = "sample", values_from = "value", values_fill = 0) %>%
    column_to_rownames(var = "OTU_ID")
  
}

#### Filter rare occurrences of abundant OTUs ###################################

# Proportional version
library_wise_filter <- function(library_specific_OTU_table, threshold) {

  # Remove rare occurrences of abundant OTUs using a 0.5% threshold.
  filtered_otu_table <- library_specific_OTU_table %>%
    rownames_to_column(var = "OTU_ID") %>%
    pivot_longer(-OTU_ID, names_to = "sample",
                 values_to = "OTU_sample_abundance") %>%
    group_by(OTU_ID) %>%
    mutate(
      OTU_library_abundance = sum(OTU_sample_abundance),
      rel_abundance = (OTU_sample_abundance / OTU_library_abundance) * 100
    ) %>%
    mutate(OTU_sample_abundance = replace(
      OTU_sample_abundance, rel_abundance <= threshold, 0)) %>%
    select(-c(OTU_library_abundance, rel_abundance)) %>%
    pivot_wider(names_from = "sample", values_from = "OTU_sample_abundance")

  return(filtered_otu_table)

}

#### Filter positive control OTUs ##############################################

# Proportional version
library_wise_filter_controls <- function(
    library_specific_OTU_table, threshold, control_OTUs
    ) {
  
  # Remove positive control OTUs using a 3% threshold.
  filtered_otu_table <- library_specific_OTU_table %>%
    rownames_to_column(var = "OTU_ID") %>%
    pivot_longer(-OTU_ID, names_to = "sample",
                 values_to = "OTU_sample_abundance") %>%
    group_by(OTU_ID) %>%
    mutate(
      OTU_library_abundance = sum(OTU_sample_abundance),
      rel_abundance = (OTU_sample_abundance / OTU_library_abundance) * 100
    ) %>%
    mutate(
      OTU_sample_abundance = case_when(
        OTU_ID %in% control_OTUs & rel_abundance <= threshold ~ 0,
        TRUE ~ OTU_sample_abundance
        )
      )%>%
    select(-c(OTU_library_abundance, rel_abundance)) %>%
    pivot_wider(names_from = "sample", values_from = "OTU_sample_abundance")
  
  return(filtered_otu_table)
  
}

#### Calculate margin of error #################################################

# Function to calculate margin of error
calculate_margin_of_error <- function(
    confidence_level = 0.95, standard_deviation, sample_size
) {
  # Find the Z-score
  z_score <- qnorm((confidence_level + 1) / 2)
  
  # Calculate sample size
  moe <- round(
    z_score * ((standard_deviation) / sqrt(sample_size)), digits = 2
    )
  
  return(moe)
}

#### Annotate blast hits based on thresholds ###################################

annotate_blast_hits <- function(data) {
  
  data %>%
    mutate(
      genus = case_when(
          coverage < 90 | similarity < genus_threshold & 
            !grepl("_Incertae_sedis$", family) ~ 
          paste0(family, "_gen_unassigned"),
          coverage < 90 | similarity < genus_threshold & 
            grepl("_Incertae_sedis$", family) ~ 
          sub("_\\w{3}_Incertae_sedis$", "_gen_unassigned", family),
        TRUE ~ genus
      ),
      family = case_when(
        coverage < 85 | similarity < family_threshold & 
          !grepl("_Incertae_sedis$", order) ~ 
          paste0(order, "_fam_unassigned"),
        coverage < 85 | similarity < family_threshold & 
          grepl("_Incertae_sedis$", order) ~ 
          sub("_\\w{3}_Incertae_sedis$", "_fam_unassigned", order),
        TRUE ~ family
      ),
      order = case_when(
        coverage < 80 | similarity < order_threshold & 
          !grepl("_Incertae_sedis$", class) ~ 
          paste0(class, "_ord_unassigned"),
        coverage < 80 | similarity < order_threshold & 
          grepl("_Incertae_sedis$", class) ~ 
          sub("_\\w{3}_Incertae_sedis$", "_ord_unassigned", class),
        TRUE ~ order
      ),
      class = case_when(
        coverage < 75 | similarity < class_threshold & 
          !grepl("_Incertae_sedis$", phylum) ~ 
          paste0(phylum, "_cls_unassigned"),
        coverage < 75 | similarity < class_threshold & 
          grepl("_Incertae_sedis$", phylum) ~ 
          sub("_\\w{3}_Incertae_sedis$", "_cls_unassigned", phylum),
        TRUE ~ class
      ),
      phylum = case_when(
        similarity < phylum_threshold & 
          !grepl("_Incertae_sedis$", kingdom) ~ 
          paste0(kingdom, "_phy_unassigned"),
        similarity < phylum_threshold & 
          grepl("_Incertae_sedis$", kingdom) ~ 
          sub("_\\w{3}_Incertae_sedis$", "_phy_unassigned", kingdom),
        TRUE ~ phylum
      )
    ) %>%
    select(
      OTU_ID, kingdom, phylum, class, order, family, genus, species, similarity,
      coverage, evalue, abundance, hit_ID
      )
}

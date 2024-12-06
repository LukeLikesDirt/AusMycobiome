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
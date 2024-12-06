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

#### MY GGPLOT THEME ###########################################################

MyTheme <- function() {
  
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(0.9)),
    legend.position = "none"
  )
  
}

#### Maps ######################################################################

# Global map
world_map <- rnaturalearth::ne_countries(scale = 'medium', type = 'map_units',
                                         returnclass = 'sf')

aus_map <- world_map %>%
  filter(name == 'Australia')

ant_map <- world_map %>%
  filter(name == 'Antarctica')

# Basic Australia layer in GGPLOT
aus_plot <-
  ggplot() +
  geom_sf(data = aus_map) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    limits = c(113, 154),
    breaks = c(120, 130, 140, 150)
  ) +
  scale_y_continuous(
    limits = c(-43.5, -10),
    breaks = c(-40, -30, -20, -10)
  ) +
  MyTheme() +
  theme(
    axis.text = element_text(size = 8)
  )

# Basic Antarctica layer in GGPLOT
ant_plot <-
  ggplot() +
  geom_sf(data = ant_map) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    limits = c(-160, 160),
    breaks = c(-150, 0, 150)
  ) +
  scale_y_continuous(
    limits = c(-83, -57),
    breaks = c(-80, -70, -60)
  ) +
  MyTheme() +
  theme(
    axis.text = element_text(size = 8)
  )

rm(world_map)



### Summary statistics functions ###############################################

summarise_abundance_richness <- function(data) {
  data %>%
    group_by(threshold) %>%
    summarise(
      mean_total_abundance = mean(total_abundance),
      mean_relative_abundance = mean(relative_abundance),
      mean_total_richness = mean(total_richness),
      mean_relative_richness = mean(relative_richness),
      std_err_total_abundance = sd(total_abundance) / sqrt(n()),
      std_err_relative_abundance = sd(relative_abundance) / sqrt(n()),
      std_err_total_richness = sd(total_richness) / sqrt(n()),
      std_err_relative_richness = sd(relative_richness) / sqrt(n()),
      ci_low_total_abundance = mean_total_abundance - qt(0.975, df = n() - 1) * std_err_total_abundance,
      ci_high_total_abundance = mean_total_abundance + qt(0.975, df = n() - 1) * std_err_total_abundance,
      ci_low_relative_abundance = mean_relative_abundance - qt(0.975, df = n() - 1) * std_err_relative_abundance,
      ci_high_relative_abundance = mean_relative_abundance + qt(0.975, df = n() - 1) * std_err_relative_abundance,
      ci_low_total_richness = mean_total_richness - qt(0.975, df = n() - 1) * std_err_total_richness,
      ci_high_total_richness = mean_total_richness + qt(0.975, df = n() - 1) * std_err_total_richness,
      ci_low_relative_richness = mean_relative_richness - qt(0.975, df = n() - 1) * std_err_relative_richness,
      ci_high_relative_richness = mean_relative_richness + qt(0.975, df = n() - 1) * std_err_relative_richness
    )
}

### Calculate abundance and richness for each group ############################

caculate_abundance_richness <- function(group_name) {
  combined_df <- bind_rows(
    
    # Threshold 0.01
    get(paste0("OTUs_", group_name, "_01")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_01")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_01",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.02
    get(paste0("OTUs_", group_name, "_02")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_02")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_02",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.05
    get(paste0("OTUs_", group_name, "_05")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_05")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_05",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.10
    get(paste0("OTUs_", group_name, "_10")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_10")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_10",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.20
    get(paste0("OTUs_", group_name, "_20")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_20")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_20",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.50
    get(paste0("OTUs_", group_name, "_50")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_50")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_50",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.01_1
    get(paste0("OTUs_", group_name, "_01_1")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_01_1")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_01_1",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.02_1
    get(paste0("OTUs_", group_name, "_02_1")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_02_1")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_02_1",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.05_1
    get(paste0("OTUs_", group_name, "_05_1")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_05_1")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_05_1",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.10_1
    get(paste0("OTUs_", group_name, "_10_1")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_10_1")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_10_1",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.20_1
    get(paste0("OTUs_", group_name, "_20_1")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_20_1")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_20_1",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0)),
    
    # Threshold 0.50_1
    get(paste0("OTUs_", group_name, "_50_1")) %>%
      group_by(sample_id) %>%
      summarise(
        total_richness = n(),
        total_abundance = sum(abundance)
      ) %>%
      full_join(get(paste0("samples_50_1")), by = "sample_id") %>%
      mutate(
        threshold = "threshold_50_1",
        relative_richness = total_richness / sample_richness * 100,
        relative_abundance = total_abundance / sample_abundance * 100
      ) %>%
      select(
        threshold,
        total_abundance, relative_abundance,
        total_richness, relative_richness
      ) %>%
      replace_na(list(total_abundance = 0, relative_abundance = 0, total_richness = 0, relative_richness = 0))
    
    # Add other thresholds similarly...
  )
  
  return(combined_df)
}

#### Calculate mold stats ######################################################

# Create a function to calculate mold stats for a given mold dataframe
calculate_mould_stats <- function(otu_tab, mould_df, mould_name) {
  otu_tab %>%
    filter(OTU_ID %in% mould_df$OTU_ID) %>%
    pivot_longer(-OTU_ID, names_to = "sample_id") %>%
    filter(value > 0) %>%
    group_by(sample_id) %>%
    summarise(
      !!sym(paste(mould_name, "_abundance", sep = "")) := sum(value),
      !!sym(paste(mould_name, "_richness", sep = "")) := n()
    )
}

#### Maps ######################################################################

# Basic Australia layer in GGPLOT
aus_sf <- rnaturalearth::ne_countries(
  scale = 'medium',
  type = 'map_units',
  returnclass = 'sf') %>%
  filter(
    name == 'Australia'
  )

aus_plot <-
  ggplot() +
  geom_sf(data = aus_sf) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    limits = c(113, 154),
    breaks = c(120, 130, 140, 150)
  ) +
  scale_y_continuous(
    limits = c(-43.5, -10),
    breaks = c(-40, -30, -20, -10)
  ) +
  MyTheme() +
  theme(
    axis.text = element_text(size = 8)
  )

# Basic Antarctica layer in GGPLOT
ant_sf <- rnaturalearth::ne_countries(
  scale = 'medium',
  type = 'map_units',
  returnclass = 'sf') %>%
  filter(
    name == 'Antarctica'
  )

ant_plot <-
  ggplot() +
  geom_sf(data = ant_sf) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    #limits = c(-160, 160),
    #breaks = c(-150, 0, 150)
    limits = c(35, 165),
    breaks = c(50, 100, 150)
  ) +
  scale_y_continuous(
    limits = c(-83, -57),
    breaks = c(-80, -70, -60)
  ) +
  MyTheme() +
  theme(
    axis.text = element_text(size = 8)
  )

# Xmas isand study extent
xmas_sf <- giscoR::gisco_get_countries(
  country = "Christmas I.", resolution = 1)

# Basic Christmas Island layer in GGPLOT
xmas_plot <-
  ggplot() +
  geom_sf(data = xmas_sf) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    breaks = c(105.55, 105.65)
  ) +
  scale_y_continuous(
    breaks = c(-10.55, -10.50, -10.45)
  ) +
  MyTheme() +
  theme(
    axis.text = element_text(size = 8)
  )

# Full study extent
aus_micro_sf <- rnaturalearth::ne_countries(
  scale = 'medium',
  type = 'map_units',
  returnclass = 'sf') %>%
  filter(
    name %in% c('Antarctica', 'Christmas I.', 'Australia')
  )

aus_micro_map <-
  ggplot() +
  geom_sf(data = aus_micro_sf) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    limits = c(78, 154),
    breaks = c(80, 90, 100, 110, 120, 130, 140, 150)
  ) +
  scale_y_continuous(
    limits = c(-70, -10),
    breaks = c(-70, -60, -50, -40, -30, -20, -10)
  ) +
  MyTheme() +
  theme(
    axis.text = element_text(size = 8)
  )

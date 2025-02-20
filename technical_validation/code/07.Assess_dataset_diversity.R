require(data.table)
require(tidyverse)

# KRONA CHARTS #################################################################

# Data
data <- fread("../output/sample_metadata.txt")
OTUs <- fread("../output/otu_table.txt")
taxa <- fread("../output/taxonomy.txt")

# Samples by continent
australian_samples <- data %>%
  filter(continent == "Australia")
antarctic_samples <- data %>%
  filter(continent == "Antarctica")
asian_samples <- data %>%
  filter(continent == "Asia")

# OTUs by continent
australian_OTUs <- OTUs %>%
  pivot_longer(
    cols = -OTU_ID,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(sample_id %in% australian_samples$sample_id) %>%
  filter(abundance > 0) %>%
  group_by(OTU_ID) %>%
  summarise(
    abundance = sum(abundance)
  )
antarctic_OTUs <- OTUs %>%
  pivot_longer(
    cols = -OTU_ID,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(sample_id %in% antarctic_samples$sample_id) %>%
  filter(abundance > 0) %>%
  group_by(OTU_ID) %>%
  summarise(
    abundance = sum(abundance)
  )

# Mean sample richness
mean_richness <- OTUs %>%
  pivot_longer(
    cols = -OTU_ID,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  summarise(
    richness = n_distinct(OTU_ID),
    abundance = sum(abundance),
    .groups = "drop"
  ) %>%
  summarise(
    mean_richness = mean(richness),
    mean_abundance = mean(abundance),
    median_richness = median(richness),
    median_abundance = median(abundance)
  )

# Proportion of rank and genus annotations
taxa %>%
  mutate(n_OTUs = n_distinct(OTU_ID)) %>%
  group_by(rank) %>%
  summarise(
    n_OTUs = first(n_OTUs),
    n_annotations = n_distinct(OTU_ID),
    proportion = n_annotations / n_OTUs
  )

# Mean richness and abundnace per sample
OTUs %>%
  pivot_longer(
    cols = -OTU_ID,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  summarise(
    richness = n_distinct(OTU_ID),
    abundance = sum(abundance),
    .groups = "drop"
  ) %>%
  summarise(
    mean_richness = mean(richness),
    sd_richness = sd(richness),
    mean_abundance = mean(abundance),
    sd_abundance = sd(abundance)
  )

# All Krona chart
taxa %>%
  select(phylum, class, order, family, genus, OTU_ID, abundance) %>%
  mutate(
    # Replace all strings containing "pseudo" with "unidentified"
    across(c(phylum, class, order, family, genus), ~ if_else(str_detect(.x, "pseudo"), "unidentified", .x)),
    
    # Standardize "Incertae_sedis" as "Incertae sedis"
    across(c(phylum, class, order, family, genus), ~ if_else(str_detect(.x, "Incertae_sedis"), "Incertae sedis", .x))
  ) %>%
  rowwise() %>%
  mutate(
    # Correct "Incertae sedis" for each rank based on all lower ranks
    class = if_else(
      class == "Incertae sedis" & all(c(order, family, genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      class
    ),
    order = if_else(
      order == "Incertae sedis" & all(c(family, genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      order
    ),
    family = if_else(
      family == "Incertae sedis" & all(c(genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      family
    )
  ) %>%
  ungroup() %>%
  group_by(phylum, class, order, family, genus) %>%
  summarise(
    richness = n_distinct(OTU_ID),
    abundance = sum(abundance),
  ) %>%
  fwrite("../output/krona_contemporary.txt", sep = "\t")

# Australian Krona chart
taxa <- select(phylum, class, order, family, genus, OTU_ID, abundance) %>%
  mutate(
    # Replace all strings containing "pseudo" with "unidentified"
    across(c(phylum, class, order, family, genus), ~ if_else(str_detect(.x, "pseudo"), "unidentified", .x)),
    
    # Standardize "Incertae_sedis" as "Incertae sedis"
    across(c(phylum, class, order, family, genus), ~ if_else(str_detect(.x, "Incertae_sedis"), "Incertae sedis", .x))
  ) %>%
  rowwise() %>%
  mutate(
    # Correct "Incertae sedis" for each rank based on all lower ranks
    class = if_else(
      class == "Incertae sedis" & all(c(order, family, genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      class
    ),
    order = if_else(
      order == "Incertae sedis" & all(c(family, genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      order
    ),
    family = if_else(
      family == "Incertae sedis" & all(c(genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      family
    )
  ) %>%
  ungroup() %>%
  group_by(phylum, class, order, family, genus) %>%
  summarise(
    richness = n_distinct(OTU_ID),
    abundance = sum(abundance),
  ) %>%
  fwrite("../output/krona_australian.txt", sep = "\t")

# Antarctic Krona chart
taxa %>% select(phylum, class, order, family, genus, OTU_ID, abundance) %>%
  mutate(
    # Replace all strings containing "pseudo" with "unidentified"
    across(c(phylum, class, order, family, genus), ~ if_else(str_detect(.x, "pseudo"), "unidentified", .x)),
    
    # Standardize "Incertae_sedis" as "Incertae sedis"
    across(c(phylum, class, order, family, genus), ~ if_else(str_detect(.x, "Incertae_sedis"), "Incertae sedis", .x))
  ) %>%
  rowwise() %>%
  mutate(
    # Correct "Incertae sedis" for each rank based on all lower ranks
    class = if_else(
      class == "Incertae sedis" & all(c(order, family, genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      class
    ),
    order = if_else(
      order == "Incertae sedis" & all(c(family, genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      order
    ),
    family = if_else(
      family == "Incertae sedis" & all(c(genus) %in% c("unidentified", "Incertae sedis")),
      "unidentified",
      family
    )
  ) %>%
  ungroup() %>%
  group_by(phylum, class, order, family, genus) %>%
  summarise(
    richness = n_distinct(OTU_ID),
    abundance = sum(abundance),
  ) %>%
  fwrite("../output/krona_antarctic.txt", sep = "\t")

# Summarise most species rich genera
genus_diversity <- OTUs %>%
  pivot_longer(
    -OTU_ID,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  filter(
    abundance > 0
  ) %>%
  inner_join(
  .,
  taxa %>% select(-abundance),
  by = "OTU_ID"
  ) %>%
  mutate(
    total_abundance = sum(abundance),
    total_richness = n_distinct(OTU_ID),
    n_samples = n_distinct(sample)
  ) %>%
  group_by(genus) %>%
  summarise(
    total_abundance = first(total_abundance),
    n_samples = first(n_samples),
    total_richness = first(total_richness),
    absolute_abundance = sum(abundance),
    relative_abundance = (absolute_abundance / total_abundance) * 100,
    richness = n_distinct(OTU_ID),
    relative_richness = (richness / total_richness) * 100,
    prevalence_count = n_distinct(sample),
    prevalence = (prevalence_count / n_samples) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(richness)) %>%
  select(-total_abundance, -n_samples, -total_richness) %>%
  filter(
    !grepl("unidentified", genus),
    !grepl("Incertae_sedis", genus)
  ) %>%
  print(n = 20)

# Plot genus diversity
genus_diversity %>%
  slice(1:10) %>%
  inner_join(
    .,
    taxa %>% select(genus, primary_lifestyle) %>% distinct(),
    by = "genus"
  ) %>%
  select(
    genus,
    `Guild` = primary_lifestyle,
    Richness = richness, 
    `Relative Abundance (%)` =  relative_abundance,
    `Sample Prevalence (%)` = prevalence,
    ) %>%
  pivot_longer(
    cols = -c(genus, `Guild`),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("Richness", "Relative Abundance (%)", "Sample Prevalence (%)")),
    genus = factor(genus, levels = genus_diversity$genus),
    `Guild` = case_when(
      `Guild` == "plant_pathogen" ~ "Plant Pathogen",
      `Guild` == "arbuscular_mycorrhizal" ~ "Arbuscular Mycorrhizal",
      `Guild` == "ectomycorrhizal" ~ "Ectomycorrhizal",
      `Guild` == "soil_saprotroph" ~ "Soil Saprotroph",
      TRUE ~ "Unclassified"
    )
  ) %>%
  ggplot(aes(x = genus, y = value, fill = `Guild`)) +
  geom_col() +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)),
    strip.text = element_text(face = "bold"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    aspect.ratio = 1
  )

# Save plot
ggsave(
  filename = "../output/plots/genus_diversity.png",
  width = 15,
  height = 9,
  units = "cm",
  dpi = 300
)
ggsave(
  filename = "../output/plots/genus_diversity.tiff",
  width = 15,
  height = 9,
  units = "cm",
)
ggsave(
  filename = "../output/plots/genus_diversity.svg",
  width = 15,
  height = 9,
  units = "cm",
)

# Summarise most species rich genera
guild_diversity <- OTUs %>%
  pivot_longer(
    -OTU_ID,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  filter(
    abundance > 0
  ) %>%
  inner_join(
    .,
    taxa %>% select(-abundance),
    by = "OTU_ID"
  ) %>%
  mutate(
    total_abundance = sum(abundance),
    total_richness = n_distinct(OTU_ID),
    n_samples = n_distinct(sample)
  ) %>%
  group_by(primary_lifestyle) %>%
  summarise(
    total_abundance = first(total_abundance),
    n_samples = first(n_samples),
    total_richness = first(total_richness),
    absolute_abundance = sum(abundance),
    relative_abundance = (absolute_abundance / total_abundance) * 100,
    richness = n_distinct(OTU_ID),
    relative_richness = (richness / total_richness) * 100,
    prevalence_count = n_distinct(sample),
    prevalence = (prevalence_count / n_samples) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(richness)) %>%
  select(-total_abundance, -n_samples, -total_richness) %>%
  print(n = 20)

# Phylum diversity
phylum_diversity <- OTUs %>%
  pivot_longer(
    -OTU_ID,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  filter(
    abundance > 0
  ) %>%
  inner_join(
    .,
    taxa %>% select(-abundance),
    by = "OTU_ID"
  ) %>%
  mutate(
    total_abundance = sum(abundance),
    total_richness = n_distinct(OTU_ID),
    n_samples = n_distinct(sample)
  ) %>%
  group_by(phylum) %>%
  summarise(
    total_abundance = first(total_abundance),
    n_samples = first(n_samples),
    total_richness = first(total_richness),
    absolute_abundance = sum(abundance),
    relative_abundance = (absolute_abundance / total_abundance) * 100,
    richness = n_distinct(OTU_ID),
    relative_richness = (richness / total_richness) * 100,
    prevalence_count = n_distinct(sample),
    prevalence = (prevalence_count / n_samples) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(richness)) %>%
  select(-total_abundance, -n_samples, -total_richness) %>%
  filter(
    !grepl("unidentified", phylum),
    !grepl("Incertae_sedis", phylum)
  ) %>%
  print(n = 20)
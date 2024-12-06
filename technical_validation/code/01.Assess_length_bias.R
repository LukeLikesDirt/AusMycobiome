
# Required packages
require(ggpubr)
require(data.table)
require(Biostrings)
require(tidyverse)

# Names for described species in Australia need to be updated
taxa_match <- fread("data/fungi_names_mycobank.txt")
# Fungal trait information
fungal_traits <- fread("data/fungal_traits/fungal_traits.csv")

# (1) Organise data ############################################################

### (1a) Described species in Australia ####
taxa_australia_described <- fread("data/FungI-taxon-2024-07-25-4053.csv") %>%
  filter(
    kingdom == "Fungi",
    taxonomicStatus == "accepted",
    #nameType == "scientific",
    taxonRank == "Species"
  ) %>%
  mutate(
    # Extract the phylum name by trimming the domain and kingdom 
    phylum = str_extract(higherClassification, "(?:[^|]*\\|){2}([^|]*)") %>% str_trim(), 
    # Further trim everything before the second '|' 
    phylum = str_replace(phylum, "^[^|]*\\|[^|]*\\|(.*)", "\\1"),
    class = class,
    species = canonicalName,
    genus = word(species, 1)
  ) %>%
  select(phylum, class, family, genus, species) %>%
  unique() %>%
  left_join(taxa_match, by = c("species" = "old_name")) %>%
  mutate(
    species = case_when(
      !is.na(new_name) ~ new_name,  # When species matches old_name, take new_name
      TRUE ~ species  # Otherwise, keep the original species name
    ),
    genus = word(species, 1)
  ) %>%
  select(-new_name)

# How many described species in Australia?
n_distinct(taxa_australia_described$species)

# Phyla and classes of described species in Australia are not up to date
# (see below) so I use genus names for filtering purposes
unique(taxa_australia_described$phylum)
unique(taxa_australia_described$class)

### (1b) Australian Microbiome data ####
sequences_ausmic <- readDNAStringSet("../output/sequences.fasta")
taxa_ausmic <- fread("../output/taxonomy.txt") %>%
  select(OTU_ID, phylum, class, genus, species) %>%
  left_join(taxa_match, by = c("species" = "old_name")) %>%
  mutate(
    species = str_replace(species, "_", " "),
    species = case_when(
      !is.na(new_name) ~ new_name,  # When species matches old_name, take new_name
      TRUE ~ species  # Otherwise, keep the original species name
    ),
    genus = word(species, 1)
  ) %>%
  select(-new_name)

### (1c) UNITE data with unique sequences ####
sequences_unite <- readDNAStringSet("data/UNITE_2024/unite2024ITS1.unique.fasta")
taxa_unite <- fread("data/UNITE_2024/unite2024ITS1.classification") %>%
  filter(
    id %in% names(sequences_unite)
  )

# How many species in Australia have sequences in unite?
taxa_australia_described %>%
  inner_join(
    .,
    taxa_unite %>% select(species) %>% unique(),
    by = "species"
  ) %>%
  select(species) %>%
  unique()

### (1d) All fungi ####
length_data_all <- bind_rows(
  # Australian Microbiome data
  tibble(
    dataset = "Aust. Microbiome",
    OTU_ID = names(sequences_ausmic),
    length = width(sequences_ausmic)
  ) %>%
    left_join(
      .,
      taxa_ausmic %>% select(OTU_ID, phylum, class, genus) %>% unique(.),
      by = "OTU_ID"
    ) %>%
    select(dataset, length, phylum, class, genus),
  # UNITE data
  tibble(
    dataset = "UNITE+INSD",
    id = names(sequences_unite),
    length = width(sequences_unite)
  ) %>%
    left_join(
      .,
      taxa_unite %>% select(id, phylum, class, genus),
      by = "id"
    ) %>%
    select(dataset, length, phylum, class, genus)
) %>%
  # Remove outliers and retain classes in Asutralia
  filter(
    length <= 400,
    length >= 50
  ) %>%
  mutate(
    group = "All fungi"
  )

### (1d) ECM fungi ####
length_data_ecm <- length_data_all %>%
  left_join(
    .,
    fungal_traits %>% select(genus, primary_lifestyle) %>% unique(),
    by = "genus") %>%
  filter(
    primary_lifestyle == "ectomycorrhizal"
  ) %>%
  select(-primary_lifestyle) %>%
  mutate(
    group = "ECM fungi"
  )

### (1e) non-ECM Agaromycetes ####
length_data_agaricomycetes <- length_data_all %>%
  left_join(
    .,
    fungal_traits %>% select(genus, primary_lifestyle) %>% unique(),
    by = "genus") %>%
  filter(
    class == "Agaricomycetes",
    primary_lifestyle != "ectomycorrhizal"
  ) %>%
  select(-primary_lifestyle) %>%
  mutate(
    group = "non-ECM Agricomycetes"
  )

### (1f) non-ECM Mucoromyceta ####
length_data_mucoromycota <- length_data_all %>%
  filter(phylum %in% c("Mucoromycota")) %>%
  mutate(
    group = "Mucoromycota"
  )

### (1g) Join the data ####
length_data_combined <- bind_rows(
  length_data_all,
  length_data_ecm,
  length_data_agaricomycetes,
  length_data_mucoromycota
)

### (1f) Plot the data ####

length_distributions_plots <- ggplot(length_data_combined, aes(x = length, fill = dataset)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 230, colour = "black", linetype = "dotted") +
  scale_fill_manual(values = c("Aust. Microbiome" = "#0057B7", "UNITE+INSD" = "#FFD700")) +
  labs(x = "Sequence length (base pairs)", y = "Density") +
  scale_x_continuous(limits = c(50, 350), breaks = c(100, 200, 300)) +
  scale_y_continuous(
    limits = c(0, 0.02),
    breaks = seq(0, 0.02, by = 0.01)
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    #legend.position = c(1, 1), 
    #legend.justification = c(1.025, 1.05),
    legend.background = element_rect(fill = "NA"),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = rel(0.8)),
    legend.margin = margin(t=0, r=0, b=0, l=0),
    aspect.ratio = 1,
    strip.text = element_text(face = "bold"),
    strip.background = element_blank()
  ) + 
  facet_wrap(~group)

ggsave("../output/plots/length_distributions.png", length_distributions_plots, width = 15.75, height = 12, units = "cm")
ggsave("../output/plots/length_distributions.tiff", length_distributions_plots, width = 15.75, height = 12, units = "cm")

# (3) Length bias plots ########################################################

### (3a) Classes #### 

# Merge length data for classes with described species in Australia
class_filter <- taxa_unite %>%
  # Filter to classes with described species in Australia using genus as a proxy
  # because not all class names are available in the Australia dataset
  filter(
    genus %in% taxa_australia_described$genus,
    !grepl("Incertae_sedis", class)
  ) %>%
  pull(class) %>%
  unique(.)

# Mean lengths for Australian classes
length_data_australian_classes <- tibble(
  id = names(sequences_unite),
  length = width(sequences_unite)
) %>%
  left_join(
    .,
    taxa_unite %>% select(id, class),
    by = "id"
  ) %>%
  select(length, class) %>%
  # Remove outliers and retain classes in Asutralia
  filter(
    length <= 400,
    length >= 50,
    class %in% class_filter
  ) %>%
  mutate(
    count = 1,
    total_count = sum(count)
  ) %>%
  group_by(class) %>%
  summarise(
    mean_length = mean(length) %>% na.omit(),
    sd = sd(length),
    se = sd / sqrt(n()),
    lower_95 = mean_length - 1.96 * se,
    upper_95 = mean_length + 1.96 * se
  ) %>%
  # Add detection status
  mutate(
    detection_status = ifelse(class %in% taxa_ausmic$class, "Detected", "Not detected")
  ) %>%
  print(n = Inf)

# Plot the data
length_bias_classes <- length_data_australian_classes %>%
  select(class, mean_length, sd, detection_status) %>%
  ggplot(aes(x = mean_length, y = class, colour = detection_status)) +
  geom_point(size = 2) +
  geom_errorbar(aes(
    xmin = mean_length - sd,
    xmax = mean_length + sd), 
    width = 0.2
  ) +
  geom_vline(xintercept = 230, colour = "black", linetype = "dotted") +
  scale_color_manual( 
    values = c("Detected" = "black", "Not detected" = "#de2d26")
  ) +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(
    limits = c(50, 350),
    breaks = c(100, 200, 300)) +
  ylab(NULL) +
  xlab(expression("Mean sequence length (±"*italic(SD)*")")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none"
    # legend.position = c(0.975, 0.975),
    # legend.justification = c("right", "top"),
    # legend.key = element_blank(),
    # legend.background = element_blank(),
    # legend.title = element_blank()
    # legend.text = element_text(size = rel(0.8)),
    # legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )  +
  annotate("text", x = 50, y = 49, label = "a", fontface = "bold", size = 5, hjust = 0, vjust = 1)


# Print the plot
print(length_bias_classes)

### (3b) Length bias ECM fungi ####

# Merge length data for genus with described species in Australia
ecm_filter <- taxa_australia_described %>%
  inner_join(
    .,
    fungal_traits %>% 
      select(genus, primary_lifestyle) %>% 
      filter(primary_lifestyle == "ectomycorrhizal") %>%
      unique(),
    by = "genus"
  ) %>%
  group_by(genus) %>%
  summarise(
    class = dplyr::first(class),
    richness = n_distinct(species)
  ) %>%
  ungroup()

# Mean lengths for Australian ECM
length_data_australian_ecm <- tibble(
  id = names(sequences_unite),
  length = width(sequences_unite)
) %>%
  left_join(
    .,
    taxa_unite %>% select(id, class, genus),
    by = "id"
  ) %>%
  select(length, class, genus) %>%
  inner_join(
    .,
    ecm_filter,
    by = c("genus", "class")
  ) %>%
  # Remove outliers and retain classes in Asutralia
  filter(
    length <= 400,
    length >= 50
  ) %>%
  group_by(genus) %>%
  summarise(
    class = dplyr::first(class),
    richness = dplyr::first(richness),
    mean_length = mean(length) %>% na.omit(),
    sd = sd(length),
    se = sd / sqrt(n()),
    lower_95 = mean_length - 1.96 * se,
    upper_95 = mean_length + 1.96 * se
  ) %>%
  arrange(desc(richness)) %>%
  # Add detection status
  mutate(
    detection_status = ifelse(genus %in% taxa_ausmic$genus, "Detected", "Not detected")
  ) %>%
  print(n = Inf)

# Plot top 20 Agaricomycete ECM from Australia
length_bias_agaricomycetes_ecm <- length_data_australian_ecm %>%
  filter(class == "Agaricomycetes") %>%
  dplyr::slice(1:20) %>%
  ggplot(aes(x = mean_length, y = genus, colour = detection_status)) +
  geom_point(size = 2) +
  geom_errorbar(aes(
    xmin = mean_length - sd,
    xmax = mean_length + sd), 
    width = 0.2
  ) +
  geom_vline(xintercept = 230, colour = "black", linetype = "dotted") +
  scale_color_manual( 
    values = c("Detected" = "black", "Not detected" = "#de2d26")
  ) +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(
    limits = c(50, 350),
    breaks = c(100, 200, 300)) +
  ylab(NULL) +
  xlab(expression("Mean sequence length (±"*italic(SD)*")")) +
  theme(
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title = element_text(size = 12),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  annotate("text", x = 50, y = 20, label = "b", fontface = "bold", size = 5, hjust = 0, vjust = 1)

# Plot top 10 non-ascomycota ECM from Australia
length_bias_ascomycota_ecm <- length_data_australian_ecm %>%
  filter(class != "Agaricomycetes") %>%
  dplyr::slice(1:10) %>%
  ggplot(aes(x = mean_length, y = genus, colour = detection_status)) +
  geom_point(size = 2) +
  geom_errorbar(aes(
    xmin = mean_length - sd,
    xmax = mean_length + sd), 
    width = 0.2
  ) +
  geom_vline(xintercept = 230, colour = "black", linetype = "dotted") +
  scale_color_manual( 
    values = c("Detected" = "black", "Not detected" = "#de2d26")
  ) +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(
    limits = c(50, 350),
    breaks = c(100, 200, 300)) +
  ylab(NULL) +
  xlab(expression("Mean sequence length (±"*italic(SD)*")")) +
  theme(
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title = element_text(size = 12),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  annotate("text", x = 50, y = 10, label = "c", fontface = "bold", size = 5, hjust = 0, vjust = 1)

### (3c) Length bias AM fung i####

# Merge length data for AM genera with described species in Australia
am_filter <- taxa_australia_described %>%
  inner_join(
    .,
    fungal_traits %>% 
      select(genus, primary_lifestyle) %>% 
      filter(primary_lifestyle == "arbuscular_mycorrhizal") %>%
      unique(),
    by = "genus"
  ) %>%
  group_by(genus) %>%
  summarise(
    class = dplyr::first(class),
    richness = n_distinct(species)
  ) %>%
  ungroup()

# Mean lengths for Australian AM
length_data_australian_am <- tibble(
  id = names(sequences_unite),
  length = width(sequences_unite)
) %>%
  left_join(
    .,
    taxa_unite %>% select(id, class, genus),
    by = "id"
  ) %>%
  select(length, class, genus) %>%
  inner_join(
    .,
    am_filter,
    by = c("genus", "class")
  ) %>%
  # Remove outliers and retain classes in Asutralia
  filter(
    length <= 400,
    length >= 50
  ) %>%
  group_by(genus) %>%
  summarise(
    class = dplyr::first(class),
    richness = dplyr::first(richness),
    mean_length = mean(length) %>% na.omit(),
    sd = sd(length),
    se = sd / sqrt(n()),
    lower_95 = mean_length - 1.96 * se,
    upper_95 = mean_length + 1.96 * se
  ) %>%
  arrange(desc(richness)) %>%
  # Add detection status
  mutate(
    detection_status = ifelse(genus %in% taxa_ausmic$genus, "Detected", "Not detected")
  ) %>%
  print(n = Inf)

# Plot top 5 AM genera from Australia
length_bias_am <- length_data_australian_am %>%
  dplyr::slice(1:5) %>%
  ggplot(aes(x = mean_length, y = genus, colour = detection_status)) +
  geom_point(size = 2) +
  geom_errorbar(aes(
    xmin = mean_length - sd,
    xmax = mean_length + sd), 
    width = 0.2
  ) +
  geom_vline(xintercept = 230, colour = "black", linetype = "dotted") +
  scale_color_manual( 
    values = c("Detected" = "black", "Not detected" = "#de2d26")
  ) +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(
    limits = c(50, 350),
    breaks = c(100, 200, 300)) +
  ylab(NULL) +
  xlab(expression("Mean sequence length (±"*italic(SD)*")")) +
  theme(
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title = element_text(size = 12),
    legend.position = "none"
    # legend.position = c(0.975, 0.025),
    # legend.justification = c("right", "bottom"),
    # legend.key = element_blank(),
    # legend.background = element_blank(),
    # legend.title = element_blank(),
    # legend.text = element_text(size = rel(1.2)),
    # legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  annotate("text", x = 50, y = 5, label = "d", fontface = "bold", size = 5, hjust = 0, vjust = 1)

### (3d) Join the plots ####

# Get the legend
legend <- get_legend(
  length_bias_agaricomycetes_ecm +
    theme(
      legend.position = "bottom",
      legend.justification = "right",
      legend.text = element_text(size = rel(1)),
      legend.title = element_blank()
    ))

# Arrange mycorrhizal plots and add legend below
length_bias_mycorrhizal <- cowplot::plot_grid(
  length_bias_agaricomycetes_ecm, length_bias_ascomycota_ecm, length_bias_am, legend,
  ncol = 1,
  rel_heights = c(4, 2.25, 1.45, 0.2)  # Adjust legend height as needed
)


# Combine all the plots
length_bias_unite <- gridExtra::arrangeGrob(
  grobs = list(length_bias_classes, length_bias_mycorrhizal),
  ncol = 2,
  widths = c(1.2, 1)
)
grid::grid.draw(length_bias_unite)
# Save the class plot
ggsave("../output/plots/length_bias_unite.png", length_bias_unite, width = 20, height = 24, units = "cm")
ggsave("../output/plots/length_bias_unite.tiff", length_bias_unite, width = 20, height = 24, units = "cm")


# Map Amanita distribution in Australia

# Required packages
require(terra)
require(data.table)
require(tidyverse)
source("code/data_curation/functions.R")

### 1. Aminita data ############################################################

#### 1.1 My OTUs ####

my_amanita <- fread(
  "output/OTUs.csv"
) %>%
  pivot_longer(
    cols = -OTU_ID,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(
    abundance > 0
  ) %>%
  select(OTU_ID, sample_id) %>%
  left_join(
    fread("output/taxonomy.csv"),
    by = "OTU_ID"
  ) %>%
  select(sample_id, genus, species) %>%
  filter(
    genus == "Amanita"
  ) %>%
  left_join(
    fread("output/sample_metadata.csv"),
    by = "sample_id"
  )

#### 1.2 Hao Amanita observations ####

hao_amanita <- fread(
  "data/technical_validation/Hao_observations/output/AAFver1.csv"
) %>%
  select(species = SPECIES, longitude = LONGITUDE, latitude = LATITUDE) %>%
  # Mutate SPECIES column to contain only the genus name
  mutate(genus = str_extract(species, "^[[:alpha:]]+")) %>%
  filter(
    genus == "Amanita"
  )

#### 1.3 UNITE amanita reference sequences ####
unite_amanita <- fread("data/technical_validation/amanita/amanita_unite.csv")

### 2. Aminita coordinates #####################################################

#### 2.1 My Aminita species ####

# Identify the Amanita species in my data
my_amanita %>%
  select(species) %>%
  arrange(species) %>%
  distinct()

# Amanita djarilmari
my_coords_djarilmari <- my_amanita %>%
  filter(
    species == "Amanita_djarilmari"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita fuscosquamosa
my_coords_fuscosquamosa <- my_amanita %>%
  filter(
    species == "Amanita_fuscosquamosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita hemibapha
my_coords_hemibapha <- my_amanita %>%
  filter(
    species == "Amanita_hemibapha"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)


# Amanita luteolovelata
my_coords_luteolovelata <- my_amanita %>%
  filter(
    species == "Amanita_luteolovelata"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita marmorata
my_coords_marmorata <- my_amanita %>%
  filter(
    species == "Amanita_marmorata"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita muscaria
my_coords_muscaria <- my_amanita %>%
  filter(
    species == "Amanita_muscaria"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita oleosa
my_coords_oleosa <- my_amanita %>%
  filter(
    species == "Amanita_oleosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita peltigera
my_coords_peltigera <- my_amanita %>%
  filter(
    species == "Amanita_peltigera"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita pseudoinculta
my_coords_pseudoinculta <- my_amanita %>%
  filter(
    species == "Amanita_pseudoinculta"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita sabulosa
my_coords_sabulosa <- my_amanita %>%
  filter(
    species == "Amanita_sabulosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  na.omit() %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita wadjukiorum
my_coords_wadjukiorum <- my_amanita %>%
  filter(
    species == "Amanita_wadjukiorum"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita xanthocephala
my_coords_xanthocephala <- my_amanita %>%
  filter(
    species == "Amanita_xanthocephala"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

#### 2.2 Hao Aminita species ####

# Amanita djarilmari
hao_coords_djarilmari <- hao_amanita %>%
  filter(
    species == "Amanita djarilmari"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita fuscosquamosa
hao_coords_fuscosquamosa <- hao_amanita %>%
  filter(
    species == "Amanita fuscosquamosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita hemibapha
hao_coords_hemibapha <- hao_amanita %>%
  filter(
    species == "Amanita hemibapha"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita luteolovelata
hao_coords_luteolovelata <- hao_amanita %>%
  filter(
    species == "Amanita luteolovelata"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita marmorata
hao_coords_marmorata <- hao_amanita %>%
  filter(
    species == "Amanita marmorata"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita muscaria
hao_coords_muscaria <- hao_amanita %>%
  filter(
    species == "Amanita muscaria"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita oleosa
hao_coords_oleosa <- hao_amanita %>%
  filter(
    species == "Amanita oleosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita peltigera
hao_coords_peltigera <- hao_amanita %>%
  filter(
    species == "Amanita peltigera"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita pseudoinculta
hao_coords_pseudoinculta <- hao_amanita %>%
  filter(
    species == "Amanita pseudoinculta"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita sabulosa
hao_coords_sabulosa <- hao_amanita %>%
  filter(
    species == "Amanita sabulosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita wadjukiorum
hao_coords_wadjukiorum <- hao_amanita %>%
  filter(
    species == "Amanita wadjukiorum"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita xanthocephala
hao_coords_xanthocephala <- hao_amanita %>%
  filter(
    species == "Amanita xanthocephala"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

#### 2.3 UNITE Amanita species ####

# Amanita djarilmari
unite_coords_djarilmari <- unite_amanita %>%
  filter(
    species == "Amanita_djarilmari"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita fuscosquamosa
unite_coords_fuscosquamosa <- unite_amanita %>%
  filter(
    species == "Amanita_fuscosquamosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita hemibapha
unite_coords_hemibapha <- unite_amanita %>%
  filter(
    species == "Amanita_hemibapha"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita luteolovelata
unite_coords_luteolovelata <- unite_amanita %>%
  filter(
    species == "Amanita_luteolovelata"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita marmorata
unite_coords_marmorata <- unite_amanita %>%
  filter(
    species == "Amanita_marmorata"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita muscaria
unite_coords_muscaria <- unite_amanita %>%
  filter(
    species == "Amanita_muscaria"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita oleosa
unite_coords_oleosa <- unite_amanita %>%
  filter(
    species == "Amanita_oleosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita peltigera
unite_coords_peltigera <- unite_amanita %>%
  filter(
    species == "Amanita_peltigera"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita pseudoinculta
unite_coords_pseudoinculta <- unite_amanita %>%
  filter(
    species == "Amanita_pseudoinculta"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita sabulosa
unite_coords_sabulosa <- unite_amanita %>%
  filter(
    species == "Amanita_sabulosa"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita wadjukiorum
unite_coords_wadjukiorum <- unite_amanita %>%
  filter(
    species == "Amanita_wadjukiorum"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

# Amanita xanthocephala
unite_coords_xanthocephala <- unite_amanita %>%
  filter(
    species == "Amanita_xanthocephala"
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

### 3. Map Aminita species #####################################################

#### 3.1 Amanita djarilmari #####
map_djarilmari <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_djarilmari,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_djarilmari,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_djarilmari,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "Amanita djarilmari", fontface = "italic")

#### 3.2 Amanita fuscosquamosa #####
map_fuscosquamosa <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_fuscosquamosa,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_fuscosquamosa,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_fuscosquamosa,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.075, vjust = 1.5,
           label = "Amanita fuscosquamosa", fontface = "italic")

#### 3.3 Amanita hemibapha #####
map_hemibapha <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_hemibapha,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_hemibapha,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_hemibapha,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "Amanita hemibapha", fontface = "italic")

#### 3.4 Amanita luteolovelata #####
map_luteolovelata <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_luteolovelata,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_luteolovelata,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_luteolovelata,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "Amanita luteolovelata", fontface = "italic")

#### 3.5 Amanita marmorata #####
map_marmorata <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_marmorata,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_marmorata,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_marmorata,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "Amanita marmorata", fontface = "italic")

#### 3.6 Amanita muscaria #####
map_muscaria <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_muscaria,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_muscaria,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_muscaria,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "Amanita muscaria", fontface = "italic")

#### 3.7 Amanita oleosa #####
map_oleosa <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_oleosa,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_oleosa,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_oleosa,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "Amanita oleosa", fontface = "italic")

#### 3.8 Amanita peltigera #####
map_peltigera <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_peltigera,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_peltigera,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_peltigera,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "Amanita peltigera", fontface = "italic")

#### 3.9 Amanita pseudoinculta #####
map_pseudoinculta <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_pseudoinculta,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_pseudoinculta,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_pseudoinculta,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5,
           label = "Amanita pseudoinculta", fontface = "italic")

#### 3.10 Amanita sabulosa ####
map_sabulosa <- ggplot() +
  geom_sf(data = aus_map) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    limits = c(113, 154),
    breaks = c(120, 130, 140, 150)
  ) +
  scale_y_continuous(
    limits = c(-43.5, -8),
    breaks = c(-40, -30, -20, -10)
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_sabulosa,
    aes(fill = "OTUs"),
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_sabulosa,
    aes(fill = "Source spec."),
    colour = "#fbaa00",
    size = 2
  ) +
  # Make fake point to include in the legend
  tidyterra::geom_spatvector(
    data = data.frame(x = 100, y = 100) %>%
      as.matrix() %>%
      vect(., crs = '+proj=longlat') %>%
      project(., aus_map),
    aes(fill = "Obs. & spec."),
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  ylim(-43.5, -8) +
  scale_fill_manual(
    limits = c("OTUs", "Source spec.", "Obs. & spec."),
    values = c("OTUs" = "#1b9e77", "Source spec." = "#fbaa00", "Obs. & spec." = "#7570b3"),
    labels = c("OTUs", "Source spec.", "Obs. & spec.")
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(0.9)),
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm"),
    # # Modify legend
    legend.position = c(0, 0),
    legend.justification = c(-0.025,-0.05),
    legend.title = element_blank(),
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.25, "cm"),
    legend.background = element_rect(fill="NA",
                                     color="black",
                                     linewidth = 0.3,
                                     linetype = 'dotted'),
    legend.margin = margin(t=1, r=3, b=2, l=3),
    legend.key = element_blank()
  ) + 
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "Amanita sabulosa", fontface = "italic")
map_sabulosa

#### 3.11 Amanita wadjukiorum #####
map_wadjukiorum <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_wadjukiorum,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_wadjukiorum,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_wadjukiorum,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.09, vjust = 1.5,
           label = "Amanita wadjukiorum", fontface = "italic")

#### 3.12 Amanita xanthocephala #####
map_xanthocephala <- aus_plot +
  tidyterra::geom_spatvector(
    data = hao_coords_xanthocephala,
    colour = alpha("#7570b3", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = my_coords_xanthocephala,
    colour = alpha("#1b9e77", 0.66),
    size = 2
  ) +
  tidyterra::geom_spatvector(
    data = unite_coords_xanthocephala,
    colour = "#fbaa00",
    size = 2
  ) +
  ylim(-43.5, -8) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.075, vjust = 1.5,
           label = "Amanita xanthocephala", fontface = "italic")

### 4. Join maps ###############################################################

amanita_maps <- patchwork::wrap_plots(
  map_djarilmari, map_fuscosquamosa, map_hemibapha, map_luteolovelata,
  map_marmorata, map_muscaria, map_oleosa, map_peltigera,
  map_pseudoinculta, map_sabulosa, map_wadjukiorum, map_xanthocephala,
  ncol = 3
) +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm"))

ggsave("output/map_amanita.jpg", amanita_maps,
       width = 16, height = 20.5, units = "cm", dpi = 600)
ggsave("output/map_amanita.tiff", amanita_maps,
       width = 16, height = 20.5, units = "cm")

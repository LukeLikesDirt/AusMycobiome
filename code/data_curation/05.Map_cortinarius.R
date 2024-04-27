
# Map Cortinarius distribution in Australia

# Required packages
require(terra)
require(data.table)
require(tidyverse)
source("code/data_curation/functions.R")

### 1. Cortinaius coordinates ##################################################

#### 1.1 My coords ####

aus_my_coords_cortinarius <- fread(
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
    genus == "Cortinarius"
  ) %>%
  left_join(
    fread("output/sample_metadata.csv"),
    by = "sample_id"
  ) %>%
  select(genus, species,  sample_id, latitude, longitude)  %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

#### 1.2 ALA coords ####

aus_ala_coords_cortinarius <- fread(
  "data/technical_validation/australian_microbiome_ALA/observations.csv"
  ) %>%
  select(
    genus, depth = verbatimDepth,
    longitude = decimalLongitude, latitude = decimalLatitude
    ) %>%
  filter(
    genus == "Cortinarius",
    # Limit the dataset to the topsoil layer
    depth == "0 cm",
    # Exclude Antarctic
    latitude > -65
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

ant_ala_coords_cortinarius <- fread(
  "data/technical_validation/australian_microbiome_ALA/observations.csv"
  ) %>%
  select(
    genus, depth = verbatimDepth,
    longitude = decimalLongitude, latitude = decimalLatitude
  ) %>%
  filter(
    genus == "Cortinarius",
    # Limit the dataset to the topsoil layer
    depth == "0 cm",
    # Exclude Antarctic
    latitude < -65
  ) %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

#### 2.3 Hao et al. (2021) coords ####

aus_hao_coords_cortinarius <- fread(
  "data/technical_validation/Hao_observations/output/AAFver1.csv"
  ) %>%
  select(SPECIES, longitude = LONGITUDE, latitude = LATITUDE) %>%
  # Mutate SPECIES column to contain only the genus name
  mutate(genus = str_extract(SPECIES, "^[[:alpha:]]+")) %>%
  filter(
    genus == "Cortinarius"
  ) %>%
  select(longitude, latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

#### 2.4 Global Fungi coords ####

aus_gf_coords_cortinarius <- fread(
  "data/technical_validation/australian_microbiome_GF/cortinarius.csv"
  ) %>%
  filter(
    # Exclude Antarctic
    latitude > -65
  ) %>%
  select(longitude, latitude) %>%
  unique(.) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

ant_gf_coords_cortinarius <- fread(
  "data/technical_validation/australian_microbiome_GF/cortinarius.csv"
  ) %>%
  filter(
    # Exclude Australia
    latitude < -65
  ) %>%
  select(longitude, latitude) %>%
  unique(.) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., aus_map)

### 2. Cortinaius maps ##################################################

#### 2.1 My map ####
plot_my_cortinarius <- aus_plot +
  tidyterra::geom_spatvector(
    data = aus_my_coords_cortinarius,
    colour = alpha("#1b9e77", 0.5),
    size = 1
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5,
           label = "c", fontface = "bold") 

#### 2.1 ALA map ####
plot_ala_cortinarius <- aus_plot +
  tidyterra::geom_spatvector(
    data = aus_ala_coords_cortinarius,
    colour = alpha("#f03b20", 0.5),
    size = 1
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5,
           label = "b", fontface = "bold") 

#### 2.1 Hao map ####
plot_hao_cortinarius <- aus_plot +
  tidyterra::geom_spatvector(
    data = aus_hao_coords_cortinarius,
    colour = alpha("#7570b3", 0.5),
    size = 1
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5,
           label = "a", fontface = "bold") 

#### 2.1 GF map ####
plot_gf_cortinarius <- aus_plot +
  tidyterra::geom_spatvector(
    data = aus_gf_coords_cortinarius,
    colour = alpha("#2C85B2", 0.5),
    size = 1
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm")
  ) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5,
           label = "d", fontface = "bold") 

### 3. Join maps ##################################################

# Wrap the maps
cortinarius_maps <- patchwork::wrap_plots(
  plot_hao_cortinarius,
  plot_ala_cortinarius,
  plot_my_cortinarius) +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "mm"))

# Save the map
ggsave("output/map_cortinarius.jpg", cortinarius_maps,
       width = 15.75, height = 5.1, units = "cm", dpi = 600)
ggsave("output/map_cortinarius.tiff", cortinarius_maps,
       width = 15.75, height = 5.1, units = "cm")



# Script:   Collect georeferneced predictors, and estimate mould contamination
#           and alpha diversity metrics
# Author:   Luke Florence
# Date:     April 21st 2024
#
# Contents:
#   (1) Organise metadata 
#   (2) Estimate alpha diversity
#   (3) Estimate mould abundance
#   (4) Collect georeferenced predictors
#   (5) Save metadata

# Required packages and functions
require(data.table)
require(ggpubr)
require(terra)
require(tidyverse)
source("code/data_curation/functions.R")

### 1. Organise metadata ######################################################

# Read in the OTU table
OTUs <- fread("output/OTUs.csv")

# Read in the taxa file
taxa <- fread("output/taxonomy.csv")

# Read in the metadata and dereplicate samples
data <- fread("data/metadata.csv") %>% 
  filter(
    sample_id %in% colnames(fread(
      "data/bioinformatics/08.OTUs/OTUs_abundance_filtered.csv")[,-1]
      )
    ) %>%
  mutate(
    sample_id = str_extract(sample_id, "(.*?)(?=_)")
    ) %>%
  filter(
    sample_id %in% colnames(OTUs[,-1])
    ) %>%
  select(
    sample_id, everything(), -c(reads, flow_id)
    ) %>%
  rename(
    flow_id = flow_id_verbatim
  ) %>%
  distinct(
    sample_id, .keep_all = TRUE
    ) %>%
  mutate(
    continent = case_when(
      latitude > -65 & latitude < -11 ~ "Australia",
      latitude < -65  ~ "Antarctica",
      latitude > -11 ~ "Asia")
  ) %>%
  glimpse(.)

# Read in the original metadata with to add URLs to source files
data_AM <- inner_join_colnames(
  read.csv(
    "data/bioplatforms/bpa_4ce7fa3d_20230713T2218/package_metadata/package_metadata_bpa_4ce7fa3d_20230713T2218_amdb-genomics-amplicon.csv"
    ) %>% select(sample_id, flow_id, url),
  read.csv(
    "data/bioplatforms/bpa_4ce7fa3d_20230713T2218/package_metadata/package_metadata_bpa_4ce7fa3d_20230713T2218_base-genomics-amplicon.csv"
    ) %>% select(sample_id, flow_id, url)
  ) %>%
  inner_join_colnames(., read.csv(
    "data/bioplatforms/bpa_4f999bb9_20230717T0504/package_metadata/package_metadata_bpa_4f999bb9_20230717T0504_amdb-genomics-amplicon.csv"
    ) %>% select(sample_id, flow_id, url)
  ) %>%
  inner_join_colnames(., read.csv(
    "data/bioplatforms/bpa_10ad26d2_20230717T0504/package_metadata/package_metadata_bpa_10ad26d2_20230717T0504_amdb-genomics-amplicon.csv"
    ) %>% select(sample_id, flow_id, url)
  ) %>%
  filter(
    flow_id %in% data$flow_id
  ) %>%
  mutate(
    sample_id_verbatim = sample_id,
    sample_id = str_extract(sample_id, "\\d+$"),
    sample_id = paste0("s", sample_id)
  ) %>%
  select(
    sample_id,
    sample_id_verbatim,
    flow_id,
    dataset_url = url
  ) %>%
  unique(.) %>%
  as_tibble()

### 2. Estimate alpha diversity ################################################

# Calculate the richness, Shannon diversity and sequence depth
seq_depth_alpha <- OTUs %>%
  pivot_longer(-OTU_ID, names_to = "sample_id") %>%
  filter(value > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(value),
    richness = n(),
    shannon_diversity = -sum(value/abundance * log(value/abundance))
  )

### 3. Estimate mould abundance ################################################

# Our detailed analysis on the impact mould (Mucorales, Mortierellales,
# Umbelopsidales, Aspergillaceae, Trichocomaceae, Bifiguratus and Trichoderma)
# abundance on alpha diversity in this dataset showed that abundance of
# Mortierellales, Umbelopsidales and Trichoderma had a particularly strong
# negative impact on alpha diversity. We will therefore calculate the relative
# abundance of these taxa as proxy to mould contamination unsing 30% threshold 

# Calculate mould abundances:
mould_rel_abundance <- taxa %>%
  filter(
    order == "Mortierellales" |
    order == "Umbelopsidales" |
    genus == "Trichoderma"
  )

# Calculate mould stats
mould_abundance <- full_join(
  seq_depth_alpha,
  calculate_mould_stats(OTUs, mould_rel_abundance, "mould_contamination"),
  by = "sample_id"
  ) %>%
  mutate(
    mould_contamination = 
      mould_contamination_abundance / abundance * 100) %>%
  select(
    sample_id,
    mould_contamination
  ) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  print()

### 4. Georeferenced predictors ################################################

#### 4.1. Climate data #####

# Read in the rasters
bio01_rast <- rast("data/explanatory_variables/bioclim/bio1.tif") 
bio02_rast <- rast("data/explanatory_variables/bioclim/bio2.tif")
bio03_rast <- rast("data/explanatory_variables/bioclim/bio3.tif")
bio04_rast <- rast("data/explanatory_variables/bioclim/bio4.tif")
bio05_rast <- rast("data/explanatory_variables/bioclim/bio5.tif")
bio06_rast <- rast("data/explanatory_variables/bioclim/bio6.tif")
bio07_rast <- rast("data/explanatory_variables/bioclim/bio7.tif")
bio08_rast <- rast("data/explanatory_variables/bioclim/bio8.tif")
bio09_rast <- rast("data/explanatory_variables/bioclim/bio9.tif")
bio10_rast <- rast("data/explanatory_variables/bioclim/bio10.tif")
bio11_rast <- rast("data/explanatory_variables/bioclim/bio11.tif")
bio12_rast <- rast("data/explanatory_variables/bioclim/bio12.tif")
bio13_rast <- rast("data/explanatory_variables/bioclim/bio13.tif")
bio14_rast <- rast("data/explanatory_variables/bioclim/bio14.tif")
bio15_rast <- rast("data/explanatory_variables/bioclim/bio15.tif")
bio16_rast <- rast("data/explanatory_variables/bioclim/bio16.tif")
bio17_rast <- rast("data/explanatory_variables/bioclim/bio17.tif")
bio18_rast <- rast("data/explanatory_variables/bioclim/bio18.tif")
bio19_rast <- rast("data/explanatory_variables/bioclim/bio19.tif")
bio20_rast <- rast("data/explanatory_variables/bioclim/bio20.tif")
bio21_rast <- rast("data/explanatory_variables/bioclim/bio21.tif")
bio22_rast <- rast("data/explanatory_variables/bioclim/bio22.tif")
bio23_rast <- rast("data/explanatory_variables/bioclim/bio23.tif")
bio24_rast <- rast("data/explanatory_variables/bioclim/bio24.tif")
bio25_rast <- rast("data/explanatory_variables/bioclim/bio25.tif")
bio26_rast <- rast("data/explanatory_variables/bioclim/bio26.tif")
bio27_rast <- rast("data/explanatory_variables/bioclim/bio27.tif")
bio28_rast <- rast("data/explanatory_variables/bioclim/bio28.tif")
bio29_rast <- rast("data/explanatory_variables/bioclim/bio29.tif")
bio30_rast <- rast("data/explanatory_variables/bioclim/bio30.tif")
bio31_rast <- rast("data/explanatory_variables/bioclim/bio31.tif")
bio32_rast <- rast("data/explanatory_variables/bioclim/bio32.tif")
bio33_rast <- rast("data/explanatory_variables/bioclim/bio33.tif")
bio34_rast <- rast("data/explanatory_variables/bioclim/bio34.tif")
bio35_rast <- rast("data/explanatory_variables/bioclim/bio35.tif")

# Project the coordinates
coords_bio <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., bio01_rast)

# Extract the values
bio_values <- full_join(
  terra::extract(bio01_rast, coords_bio),
  terra::extract(bio02_rast, coords_bio),
  by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio03_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio04_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio05_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio06_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio07_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio08_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio09_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio10_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio11_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio12_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio13_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio14_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio15_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio16_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio17_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio18_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio19_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio20_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio21_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio22_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio23_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio24_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio25_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio26_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio27_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio28_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio29_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio30_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio31_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio32_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio33_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio34_rast, coords_bio),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bio35_rast, coords_bio),
    by = "ID"
  ) %>%
  glimpse(.)
  

#### 4.2. Soil data #####

# Read in the rasters
SBIO1_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO1_15cm.tif")
SBIO1_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO1_5cm.tif")
SBIO10_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO10_15cm.tif")
SBIO10_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO10_5cm.tif")
SBIO11_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO11_15cm.tif")
SBIO11_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO11_5cm.tif")
SBIO2_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO2_15cm.tif")
SBIO2_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO2_5cm.tif")
SBIO3_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO3_15cm.tif")
SBIO3_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO3_5cm.tif")
SBIO4_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO4_15cm.tif")
SBIO4_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO4_5cm.tif")
SBIO5_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO5_15cm.tif")
SBIO5_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO5_5cm.tif")
SBIO6_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO6_15cm.tif")
SBIO6_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO6_5cm.tif")
SBIO7_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO7_15cm.tif")
SBIO7_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO7_5cm.tif")
SBIO8_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO8_15cm.tif")
SBIO8_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO8_5cm.tif")
SBIO9_15cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO9_15cm.tif")
SBIO9_5cm_rast <- rast("data/explanatory_variables/soil_climate/SBIO9_5cm.tif")
AWC_15cm_rast <- rast("data/explanatory_variables/soil_grid/AWC_15cm.tif") %>%
  rename(AWC_15cm = AWC_005_015_EV_N_P_AU_TRN_N_20210614)
AWC_5cm_rast <- rast(
  "data/explanatory_variables/soil_grid/AWC_5cm.tif"
  ) %>%
  rename(AWC_5cm = AWC_000_005_EV_N_P_AU_TRN_N_20210614)
bulk_density_15cm <- rast(
  "data/explanatory_variables/soil_grid/bulk_density_15cm.tif"
  ) %>%
  rename(bulk_density_15cm = BDW_005_015_EV_N_P_AU_TRN_N_20230607)
bulk_density_5cm <- rast(
  "data/explanatory_variables/soil_grid/bulk_density_5cm.tif"
  ) %>%
  rename(bulk_density_5cm = BDW_000_005_EV_N_P_AU_TRN_N_20230607)
CEC_15cm <- rast(
  "data/explanatory_variables/soil_grid/CEC_15cm.tif"
  ) %>%
  rename(CEC_15cm = CEC_005_015_EV_N_P_AU_TRN_N_20220826)
CEC_5cm <- rast(
  "data/explanatory_variables/soil_grid/CEC_5cm.tif"
  ) %>%
  rename(CEC_5cm = CEC_000_005_EV_N_P_AU_TRN_N_20220826)
clay_15cm <- rast(
  "data/explanatory_variables/soil_grid/clay_15cm.tif"
  ) %>%
  rename(clay_15cm = CLY_005_015_EV_N_P_AU_TRN_N_20210902)
clay_5cm <- rast(
  "data/explanatory_variables/soil_grid/clay_5cm.tif"
  ) %>%
  rename(clay_5cm = CLY_000_005_EV_N_P_AU_TRN_N_20210902)
ECE_15cm <- rast(
  "data/explanatory_variables/soil_grid/ECE_15cm.tif"
  )
ECE_5cm <- rast(
  "data/explanatory_variables/soil_grid/ECE_5cm.tif"
  )
MAOC_density_15cm <- rast(
  "data/explanatory_variables/soil_grid/MAOC_density_15cm.tif"
  ) %>%
  rename(MAOC_density_15cm = SOF_005_015_EV_N_P_AU_TRN_N_20221006_Fraction_Density_MAOC)
MAOC_density_5cm <- rast(
  "data/explanatory_variables/soil_grid/MAOC_density_5cm.tif"
  ) %>%
  rename(MAOC_density_5cm = SOF_000_005_EV_N_P_AU_TRN_N_20221006_Fractions_Density_MAOC)
MAOC_proportion_15cm <- rast(
  "data/explanatory_variables/soil_grid/MAOC_proportion_15cm.tif"
  ) %>%
  rename(MAOC_proportion_15cm = SOF_015_030_05_N_P_AU_TRN_N_20221006_Proportion_MAOC)
MAOC_proportion_5cm <- rast(
  "data/explanatory_variables/soil_grid/MAOC_proportion_5cm.tif"
  ) %>%
  rename(MAOC_proportion_5cm = SOF_000_005_EV_N_P_AU_TRN_N_20221006_Proportion_MAOC)
N_total_15cm <- rast(
  "data/explanatory_variables/soil_grid/N_total_15cm.tif"
  )
N_total_5cm <- rast(
  "data/explanatory_variables/soil_grid/N_total_5cm.tif"
  )
P_available_15cm <- rast(
  "data/explanatory_variables/soil_grid/P_available_15cm.tif"
  ) %>%
  rename(P_available_15cm = AVP_005_015_EV_N_P_AU_TRN_N_20220826)
P_available_5cm <- rast(
  "data/explanatory_variables/soil_grid/P_available_5cm.tif"
  ) %>%
  rename(P_available_5cm = AVP_000_005_EV_N_P_AU_TRN_N_20220826)
P_total_15cm <- rast(
  "data/explanatory_variables/soil_grid/P_total_15cm.tif"
  )
P_total_5cm <- rast(
  "data/explanatory_variables/soil_grid/P_total_5cm.tif"
  )
pH_CaCl2_15cm <- rast(
  "data/explanatory_variables/soil_grid/pH_CaCl2_15cm.tif"
  )
pH_CaCl2_5cm <- rast(
  "data/explanatory_variables/soil_grid/pH_CaCl2_5cm.tif"
  )
pH_water_15cm <- rast(
  "data/explanatory_variables/soil_grid/pH_H2O_15cm.tif"
  ) %>%
  rename(pH_water_15cm = PHW_005_015_EV_N_P_AU_TRN_N_20220520)
pH_water_5cm <- rast(
  "data/explanatory_variables/soil_grid/pH_H2O_5cm.tif"
  ) %>%
  rename(pH_water_5cm = PHW_000_005_EV_N_P_AU_TRN_N_20222005)
POC_density_15cm <- rast(
  "data/explanatory_variables/soil_grid/POC_density_15cm.tif"
  ) %>%
  rename(POC_density_15cm = SOF_005_015_EV_N_P_AU_TRN_N_20221006_Fraction_Density_POC)
POC_density_5cm <- rast(
  "data/explanatory_variables/soil_grid/POC_density_5cm.tif"
  ) %>%
  rename(POC_density_5cm = SOF_000_005_EV_N_P_AU_TRN_N_20221006_Fractions_Density_POC)
POC_proportion_15cm <- rast(
  "data/explanatory_variables/soil_grid/POC_proportion_15cm.tif"
  ) %>%
  rename(POC_proportion_15cm = SOF_005_015_EV_N_P_AU_TRN_N_20221006_Proportion_POC)
POC_proportion_5cm <- rast(
  "data/explanatory_variables/soil_grid/POC_proportion_5cm.tif"
  ) %>%
  rename(POC_proportion_5cm = SOF_000_005_EV_N_P_AU_TRN_N_20221006_Proportion_POC)
PyOC_density_15cm <- rast(
  "data/explanatory_variables/soil_grid/PyOC_density_15cm.tif"
  ) %>%
  rename(PyOC_density_15cm = SOF_005_015_EV_N_P_AU_TRN_N_20221006_Fraction_Density_PyOC)
PyOC_density_5cm <- rast(
  "data/explanatory_variables/soil_grid/PyOC_density_5cm.tif"
  ) %>%
  rename(PyOC_density_5cm = SOF_000_005_EV_N_P_AU_TRN_N_20221006_Fractions_Density_PyOC)
PyOC_proportion_15cm <- rast(
  "data/explanatory_variables/soil_grid/PyOC_proportion_15cm.tif"
  ) %>%
  rename(PyOC_proportion_15cm = SOF_005_015_EV_N_P_AU_TRN_N_20221006_Proportion_PyOC)
PyOC_proportion_5cm <- rast(
  "data/explanatory_variables/soil_grid/PyOC_proportion_5cm.tif"
  ) %>%
  rename(PyOC_proportion_5cm = SOF_000_005_EV_N_P_AU_TRN_N_20221006_Proportion_PyOC)
sand_15cm <- rast(
  "data/explanatory_variables/soil_grid/sand_15cm.tif"
  ) %>%
  rename(sand_15cm = SND_005_015_EV_N_P_AU_TRN_N_20210902)
sand_5cm <- rast(
  "data/explanatory_variables/soil_grid/sand_5cm.tif"
  ) %>%
  rename(sand_5cm = SND_000_005_EV_N_P_AU_TRN_N_20210902)
silt_15cm <- rast(
  "data/explanatory_variables/soil_grid/silt_15cm.tif"
  ) %>%
  rename(silt_15cm = SLT_005_015_95_N_P_AU_TRN_N_20210902)
silt_5cm <- rast(
  "data/explanatory_variables/soil_grid/silt_5cm.tif"
  ) %>%
  rename(silt_5cm = SLT_000_005_EV_N_P_AU_TRN_N_20210902)
SOC_15cm <- rast(
  "data/explanatory_variables/soil_grid/SOC_15cm.tif"
  ) %>%
  rename(SOC_15cm = SOC_005_015_EV_N_P_AU_NAT_N_20120404)
SOC_5cm <- rast(
  "data/explanatory_variables/soil_grid/SOC_5cm.tif"
  ) %>%
  rename(SOC_5cm = SOC_000_005_EV_N_P_AU_NAT_N_20120404)
soil_depth_rast <- rast(
  "data/explanatory_variables/soil_grid/soil_depth.tif"
  ) %>%
  rename(soil_depth = DES_000_200_EV_N_P_AU_NAT_C_20190901)

# Project the coordinates
coords_soil_climate <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., SBIO1_5cm_rast)
coords_soil_grid <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., AWC_5cm_rast)

# Extract the values
soil_climate_values <- full_join(
  terra::extract(SBIO1_5cm_rast, coords_soil_climate),
  terra::extract(SBIO1_15cm_rast, coords_soil_climate),
  by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO2_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO2_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO3_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO3_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO4_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO4_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO5_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO5_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO6_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO6_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO7_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO7_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO8_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO8_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO9_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO9_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO10_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO10_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO11_5cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SBIO11_15cm_rast, coords_soil_climate),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(AWC_5cm_rast, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(AWC_15cm_rast, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bulk_density_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(bulk_density_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(CEC_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(CEC_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(clay_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(clay_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(ECE_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(ECE_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(MAOC_density_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(MAOC_density_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(MAOC_proportion_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(MAOC_proportion_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(N_total_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(N_total_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(P_available_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(P_available_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(P_total_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(P_total_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(pH_CaCl2_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(pH_CaCl2_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(pH_water_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(pH_water_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(POC_density_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(POC_density_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(POC_proportion_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(POC_proportion_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(PyOC_density_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(PyOC_density_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(PyOC_proportion_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(PyOC_proportion_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(sand_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(sand_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(silt_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(silt_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SOC_5cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(SOC_15cm, coords_soil_grid),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(soil_depth_rast, coords_soil_grid),
    by = "ID"
  ) %>%
  glimpse(.)

#### 4.3. Vegetation data #####

# Read in the rasters
woody_veg_cover_rast <- 
  rast("data/explanatory_variables/vegetation_structure/vegetation_cover.tif") %>%
  tidyterra::rename(woody_veg_cover = fcov_totalc1000)
woody_veg_height_rast <- 
  rast("data/explanatory_variables/vegetation_structure/vegetation_peak_cover_height.tif") %>%
  tidyterra::rename(woody_vege_height = hModec1000satcor)
woody_veg_max_height_rast <- 
  rast("data/explanatory_variables/vegetation_structure/vegetation_max_height.tif") %>%
  tidyterra::rename(woody_veg_max_height = hp95c1000satcor)
forest_structure_rast <-
  rast("data/explanatory_variables/vegetation_structure/structural_classification.tif") %>%
  rename(forest_structure = structural_classification)
plant_comp_uniqueness_rast <- 
  rast("data/explanatory_variables/plant_diversity/plant_div_importance.tif") %>%
  rename(plant_comp_uniqueness = richness_prediction_Aus9s_gda94)
plant_div_importance_rast <- 
  rast("data/explanatory_variables/plant_diversity/plant_comp_uniqueness.tif") %>%
  rename(plant_div_importance = MeanSimilarity_res10km)
plant_richness_rast <- 
  rast("data/explanatory_variables/plant_diversity/plant_richness.tif") %>%
  rename(plant_richness = richness_prediction_Aus3s_gda94)

# Project the coordinates
coords_vegetation_structure <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., forest_structure_rast)

coords_plant_diversity <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., plant_richness_rast)

# Extract the values
vegetation_values <- full_join(
  terra::extract(forest_structure_rast, coords_vegetation_structure),
  terra::extract(woody_veg_cover_rast, coords_vegetation_structure),
  by = "ID"
  ) %>%
  full_join(
    ., terra::extract(woody_veg_height_rast, coords_vegetation_structure),
    by = "ID") %>%
  full_join(
    ., terra::extract(woody_veg_max_height_rast, coords_vegetation_structure),
    by = "ID") %>%
  full_join(
    ., terra::extract(plant_richness_rast, coords_plant_diversity),
    by = "ID") %>%
  full_join(
    ., terra::extract(plant_comp_uniqueness_rast, coords_plant_diversity),
    by = "ID") %>%
  full_join(
    ., terra::extract(plant_div_importance_rast, coords_plant_diversity),
    by = "ID") %>%
  glimpse(.)
    
#### 4.4. Ecosystem data #####

# Read in the rasters and vectors
ecosystem_site_condition_rast <- rast(
    "data/explanatory_variables/habitat_condition/ecosystem_site_condition.tif"
    ) 
habitat_condition_rast <- rast(
  "data/explanatory_variables/habitat_condition/habitat_condition.tif"
  )
local_pressures_rast <- rast(
  "data/explanatory_variables/habitat_condition/local_pressures.tif"
  )
elevation_rast <- rast("data/explanatory_variables/elevation.tif")
bioregions_vect <- vect(
  "data/explanatory_variables/IBRA7_regions/ibra7_regions.shp"
  ) %>%
  select(bioregion_name = REG_NAME_7, bioregion_code = REG_CODE_7)

# Project the coordinates
coords_habitat_condition <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., habitat_condition_rast)
coords_elevation <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., elevation_rast)
coords_bioregions <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat') %>%
  project(., bioregions_vect)

# Not sure hoe to geab information from vectors in terra so here's a workaround
# for now

# Convert coordinates to a spatial points object
coords_sf <- sf::st_as_sf(coords_bioregions, coords = c("x", "y"))

# Create a spatial points data frame with the bioregion names
bioregions_sf <- sf::st_as_sf(bioregions_vect, coords = c("x", "y"))

# Validate and simplify geometries
coords_sf <- sf::st_make_valid(coords_sf)
bioregions_sf <- sf::st_make_valid(bioregions_sf)

# Extract the values
ecosystem_values <- full_join(
  terra::extract(ecosystem_site_condition_rast, coords_habitat_condition),
  terra::extract(habitat_condition_rast, coords_habitat_condition),
  by = "ID"
  ) %>%
  full_join(
    ., terra::extract(local_pressures_rast, coords_habitat_condition),
    by = "ID"
  ) %>%
  full_join(
    ., terra::extract(elevation_rast, coords_elevation),
    by = "ID"
  ) %>%
  full_join(
    ., sf::st_join(coords_sf, bioregions_sf) %>%
      as_tibble() %>%
      select(bioregion_name, bioregion_code) %>%
      mutate(ID = as.numeric(row.names(.))),
    by = "ID"
  ) %>%
  glimpse(.)

#### 4.5. Combine predictors #####

predictor_variables <- full_join(
  data %>%
    select(sample_id) %>%
    mutate(ID = as.numeric(row.names(.))),
  bio_values,
  by = "ID"
  ) %>%
  full_join(
    ., soil_climate_values,
    by = "ID"
  ) %>%
  full_join(
  ., vegetation_values,
  by = "ID"
  ) %>%
  full_join(
  ., ecosystem_values,
  by = "ID"
  ) %>%
  glimpse(.)

# Sample metadata
sample_metadata <- data %>%
  inner_join(., data_AM, by = c("sample_id", "flow_id")) %>%
  inner_join(., seq_depth_alpha, by = "sample_id") %>%
  inner_join(., mould_abundance, by = "sample_id") %>%
  inner_join(., predictor_variables, by = "sample_id") %>%
  select(
    sample_id,
    sample_id_verbatim,
    dataset_url,
    date,
    latitude,
    longitude,
    continent,
    bioregion_name,
    bioregion_code,
    abundance,
    richness,
    shannon_diversity,
    mould_contamination,
    ammonium_verbatim,
    clay_verbatim,
    conductivity_verbatim,
    elevation_verbatim = elev_verbatim,
    exc_aluminium_verbatim,
    exc_calcium_verbatim,
    exc_magnesium_verbatim,
    exc_potassium_verbatim,
    exc_sodium_verbatim,
    nitrate_verbatim,
    organic_carbon_verbatim,
    ph_verbatim,
    phosphorus_colwell_verbatim,
    potassium_colwell_verbatim,
    sand_verbatim,
    silt_verbatim,
    sulphur_verbatim,
    texture_verbatim,
    vegetation_type_verbatim,
    bio1,
    bio2,
    bio3,
    bio4,
    bio5,
    bio6,
    bio7,
    bio8,
    bio9,
    bio10,
    bio11,
    bio12,
    bio13,
    bio14,
    bio15,
    bio16,
    bio17,
    bio18,
    bio19,
    bio20,
    bio21,
    bio22,
    bio23,
    bio24,
    bio25,
    bio26,
    bio27,
    bio28,
    bio29,
    bio30,
    bio31,
    bio32,
    bio33,
    bio34,
    bio35,
    SBIO1_5cm,
    SBIO1_15cm,
    SBIO2_5cm,
    SBIO2_15cm,
    SBIO3_5cm,
    SBIO3_15cm,
    SBIO4_5cm,
    SBIO4_15cm,
    SBIO5_5cm,
    SBIO5_15cm,
    SBIO6_5cm,
    SBIO6_15cm,
    SBIO7_5cm,
    SBIO7_15cm,
    SBIO8_5cm,
    SBIO8_15cm,
    SBIO9_5cm,
    SBIO9_15cm,
    SBIO10_5cm,
    SBIO10_15cm,
    SBIO11_5cm,
    SBIO11_15cm,
    AWC_5cm,
    AWC_15cm,
    bulk_density_5cm,
    bulk_density_15cm,
    CEC_5cm,
    CEC_15cm,
    clay_5cm,
    clay_15cm,
    ECE_5cm,
    ECE_15cm,
    MAOC_density_5cm,
    MAOC_density_15cm,
    MAOC_proportion_5cm,
    MAOC_proportion_15cm,
    N_total_5cm,
    N_total_15cm,
    P_available_5cm,
    P_available_15cm,
    P_total_5cm,
    P_total_15cm,
    pH_CaCl2_5cm,
    pH_CaCl2_15cm,
    pH_water_5cm,
    pH_water_15cm,
    POC_density_5cm,
    POC_density_15cm,
    POC_proportion_5cm,
    POC_proportion_15cm,
    PyOC_density_5cm,
    PyOC_density_15cm,
    PyOC_proportion_5cm,
    PyOC_proportion_15cm,
    sand_5cm,
    sand_15cm,
    silt_5cm,
    silt_15cm,
    SOC_5cm,
    SOC_15cm,
    soil_depth,
    forest_structure,
    woody_veg_cover,
    woody_vege_height,
    woody_veg_max_height,
    plant_richness,
    plant_comp_uniqueness,
    plant_div_importance,
    ecosystem_site_condition,
    habitat_condition,
    local_pressures,
    elevation
    ) %>%
  glimpse(.)

# Save the metadata file
fwrite(sample_metadata, "output/sample_metadata.csv")

# Clear the environment
rm(list = ls())


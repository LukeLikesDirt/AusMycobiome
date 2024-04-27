
# Script:   Normalise OTUs to account for differences in sequencing depth across
#           samples
# Author:   Luke Florence
# Date:     February 2nd 2024

# Required packages and functions
require(data.table)
require(SRS)
require(ggpubr)
require(patchwork)
require(scales)
require(vegan)
require(tidyverse)
source("code/data_curation/functions.R")

# Read in the OTU table
OTUs <- fread(
  "output/OTUs.csv"
  ) %>%
  column_to_rownames(var = "OTU_ID") %>%
  filter(rowSums(.) > 0)
# Grab the OTU IDs
OTU_IDs <- rownames(OTUs)

# Read in the metadata
data <- fread("output/sample_metadata.csv")

# Normalise OTUs using the scaling with ranked subsampling method
OTUs_SRS <- SRS(OTUs, Cmin = 5000, set_seed = TRUE, seed = 1986) %>%
  as.data.frame(.) %>%
  mutate(OTU_ID = OTU_IDs) %>%
  column_to_rownames(var = "OTU_ID") %>%
  filter(rowSums(.) > 0) %>%
  rownames_to_column(var = "OTU_ID")

# Save the normalised OTUs
fwrite(OTUs_SRS, "output/OTUs_SRS.csv")

# Traditional rarefaction of OTUs
OTUs_rare <- rarefy(
  OTUs %>%
    t(.) %>%
    as.data.frame(.),
  sample = 5000) %>%
  as_tibble(
    rownames = "sample_id"
  ) %>%
  rename(richness_rare = value)

# Non-normalised OTUs with OTU IDs
OTUs <- OTUs %>%
  rownames_to_column(var = "OTU_ID")

# Check that OTUs have values greater than 0
min(rowSums(OTUs %>% select(-OTU_ID)))
min(rowSums(OTUs_SRS %>% select(-OTU_ID)))

# Richness and abundance of OTUs
seq_depth_richness <- data %>%
  select(sample_id, continent, abundance, richness_raw = richness) %>%
  filter(continent != "Asia") %>%
  inner_join(
    OTUs_SRS %>%
      pivot_longer(-OTU_ID, names_to = "sample_id") %>%
      filter(value > 0) %>%
      group_by(sample_id) %>%
      summarise(
        richness_SRS = n()
      ),
    by = "sample_id"
  ) %>%
  inner_join(
    OTUs_rare,
    by = "sample_id"
  ) 

# Gram max and min values for plots
max_richness <- max(seq_depth_richness$richness_raw)
min_richness <- min(seq_depth_richness$richness_raw)
max_abundance <- max(seq_depth_richness$abundance)
min_abundance <- min(seq_depth_richness$abundance)

### Effect of Sequencing depth on OTU richness #################################

seq_depth_richness_plot <- seq_depth_richness %>%
  ggplot(
    aes(abundance, richness_raw, colour = continent)
  ) +
  geom_point(
    size = 1, alpha = 0.33
  ) +
  scale_colour_manual(values = c("#7570b3", "#d95f02")) +
  stat_smooth(method = lm, data = seq_depth_richness %>% 
                filter(continent == "Australia"),
              colour = "#d95f02") +
  stat_smooth(method = lm, data = seq_depth_richness %>% 
                filter(continent == "Antarctica"),
              colour = "#7570b3") +
  stat_cor(data = seq_depth_richness %>% 
             filter(continent == "Australia"),
           aes(label = after_stat(rr.label)), method = "pearson",
           label.y = (max_richness * 1.075), label.x = log10(150000),
           colour = "#d95f02", r.accuracy = 0.0001) +
  stat_cor(data = seq_depth_richness %>% 
             filter(continent == "Antarctica"),
           aes(label = after_stat(rr.label)), method = "pearson",
           label.y = (max_richness * 0.975), label.x = log10(150000),
           colour = "#7570b3", r.accuracy = 0.001) +
  scale_x_log10(
    breaks = c(5000, 20000, 100000, 500000),
    labels = scales::comma_format(scale = 1, big.mark = " ")
  ) +
  scale_y_continuous(limits = c(0, (max_richness * 1.1)),
                     labels = comma_format(big.mark = ' ')) +
  labs(y = "Observed OTU richness", x = "Logarithm of read abundance") +
  annotate("text", x = 5000, y = (max_richness * 1.075), label = "a",
           size = 5, fontface = "bold") +
  annotate("text", x = 120000, y = (max_richness * 1.075),
           label = "Australia", size = 3.75, colour = "#d95f02", 
           vjust = 0.7, hjust = 1) +
  annotate("text", x = 120000, y = (max_richness * 0.975),
           label = "Antarctica", size = 3.75, colour = "#7570b3",
           vjust = 0.7, hjust = 1) +
  MyTheme()
  
### Effect of rarefaction and normalisation on richness #######################

max_antarctica_richness <- seq_depth_richness %>%
  filter(continent == "Antarctica") %>%
  summarise(max(richness_raw)) %>%
  pull()

# Antarctica rarefied OTU richness
seq_depth_richness %>%
  filter(continent == "Antarctica") %>%
  ggplot(aes(richness_raw, richness_rare)) +
  geom_point(size = 3, alpha = 0.33, colour = "#7570b3") +
  stat_smooth(method = lm, colour = "black", linewidth = 0.5) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "bottom", label.x.npc = 0.75,
           colour = "#7570b3") +
  scale_y_continuous(
    limits = c(0, max_antarctica_richness)
  ) +
  scale_x_continuous(
    limits = c(0, max_antarctica_richness)
  ) +
  labs(y = 'Observed OTU richness', x = "Rarefied OTU richness") +
  MyTheme()

# Antarctica SRS OTU richness
seq_depth_richness %>%
  filter(continent == "Antarctica") %>%
  ggplot(aes(richness_raw, richness_SRS)) +
  geom_point(size = 3, alpha = 0.33, colour = "#7570b3") +
  stat_smooth(method = lm, colour = "black", linewidth = 0.5) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "bottom", label.x.npc = 0.75,
           colour = "#7570b3") +
  scale_y_continuous(
    limits = c(0, max_antarctica_richness)
  ) +
  scale_x_continuous(
    limits = c(0, max_antarctica_richness)
  ) +
  labs(y = 'Observed OTU richness', x = "SRS OTU richness") +
  MyTheme()

# Australia rarefied OTU richness
seq_depth_richness %>%
  filter(continent == "Australia") %>%
  ggplot(aes(richness_raw, richness_rare)) +
  geom_point(size = 3, alpha = 0.33, colour = "#d95f02") +
  stat_smooth(method = lm, colour = "black", linewidth = 0.5) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "bottom", label.x.npc = 0.75,
           colour = "#d95f02") +
  labs(y = 'Observed OTU richness', x = "Rarefied OTU richness") +
  theme(axis.text.y = element_blank()) +
  MyTheme()

# Australia SRS OTU richness
seq_depth_richness %>%
  filter(continent == "Australia") %>%
  ggplot(aes(richness_raw, richness_rare)) +
  geom_point(size = 3, alpha = 0.33, colour = "#d95f02") +
  stat_smooth(method = lm, colour = "black", linewidth = 0.5) +
  stat_cor(aes(label = after_stat(rr.label)), method = "pearson",
           label.y.npc = "bottom", label.x.npc = 0.75,
           colour = "#d95f02") +
  labs(y = 'Observed OTU richness', x = "SRS OTU richness") +
  theme(axis.text.y = element_blank()) +
  MyTheme()

### Rarefaction curve ########################################################

# Read in the OTU table
OTUs <- fread(
  "output/OTUs.csv"
  ) %>%
  pivot_longer(-OTU_ID, names_to = 'sample') %>%
  pivot_wider(names_from = 'OTU_ID', values_from = 'value') %>%
  column_to_rownames(var = 'sample')

rare_curve <- rarecurve(OTUs, step = 100)

rare_curve_plot <- map_dfr(rare_curve, bind_rows) %>%
  bind_cols(sample = rownames(OTUs),.) %>%
  pivot_longer(-sample) %>%
  drop_na() %>%
  mutate(depth = as.numeric(str_replace(name, 'N', ''))) %>%
  # Adjusting depth: mutate depths greater than 20,000 to 19,999 for readability
  mutate(depth = ifelse(depth > 20000, 19999, depth)) %>% 
  select(-name) %>%
  ggplot(aes(x = depth, y = value, sample = sample)) +
  geom_vline(xintercept = min_abundance, colour = 'grey') +
  geom_line(
    linewidth = 0.1,
  ) +
  scale_x_continuous(labels = comma_format(big.mark = ' '),
                     limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000)) +
  scale_y_continuous(limits = c(0, (max_richness * 1.1)),
                     labels = comma_format(big.mark = ' ')) +
  xlab('Read abundance') +
  ylab(NULL) +
  annotate("text", x = 10, y = (max_richness * 1.075), label = "b",
           size = 5, fontface = "bold") +
  MyTheme()

# Join and save the plots
patchwork::wrap_plots(seq_depth_richness_plot, rare_curve_plot, ncol = 2)
ggsave("output/sequencing_depth_richness.jpg", width = 6, height = 3.15,
       dpi = 300)
ggsave("output/sequencing_depth_richness.tiff", width = 6, height = 3.15)

# Clear the environment
rm(list = ls())

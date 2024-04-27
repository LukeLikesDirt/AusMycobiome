
# Script:   Quality check: Impact of mould contamination on richness.
# Author:   Luke Florence
# Date:     21st April 2024

# Required functions
require(data.table)
require(ggpubr)
require(tidyverse)
source("code/data_curation/functions.R")

# Sample metadata
data <- fread("output/sample_metadata.csv")

# Our detailed analysis on the impact mould (Mucorales, Mortierellales,
# Umbelopsidales, Aspergillaceae, Trichocomaceae, Bifiguratus and Trichoderma)
# abundance on alpha diversity in this dataset showed that abundance of
# Mortierellales, Umbelopsidales and Trichoderma had a particularly strong
# negative impact on alpha diversity. We therefore calculated the relative
# abundance of these taxa as proxy to "mould contamination" unsing 30% threshold.
# The 30% threshold is based on mean mould relative abundance + 3 standard
# deviations, as per Tedersoo et al. (2021) The Global Soil Mycobiome consortium
# dataset for boosting fungal diversity research. Fungal Diversity, 111, 573-588. 

### 1. Calculate mould statistics #############################################

# Calculate mould relative abundance
mould_data <- data %>%
  filter(continent == "Australia") %>%
  select(
    sample_id, richness, mould_contamination
  ) %>%
  filter(mould_contamination > 0)

# Calculate median mould relative abundance:
Q2_mould <- mould_data %>%
  summarise(median(mould_contamination))  %>%
  as.numeric()

# Calculate mean and standard deviation
mean_mould <- mould_data %>%
  summarise(mean = mean(mould_contamination)) %>%
  pull(mean)
sd_mould <- mould_data %>%
  summarise(sd = sd(mould_contamination)) %>%
  pull(sd)

# Grab putative contaminated samples
mould_outliers <- mould_data %>%
  filter(mould_contamination > (mean_mould + (3 * sd_mould))) %>%
  select(sample_id, richness, mould_contamination)

# Find the minimum value of upper limit outliers:
min_mould <- min(mould_outliers$mould_contamination)

### 2. Visualise the impact of mould contamination #############################

# Impact of mould relative abundance on OTU richness
plot_with_mould <- mould_data %>%
  mutate(mould_contamination = mould_contamination + 1) %>%
  ggplot(aes(mould_contamination, richness)) +
  geom_point(size = 1, color = "black", alpha = 0.5) +
  stat_smooth(
    data = mould_data %>%
      mutate(mould_contamination = mould_contamination + 1),
    method = loess, se = FALSE, colour = "#de2d26"
  ) +
  stat_smooth(
    method = lm, data = mould_data %>%
      mutate(mould_contamination = mould_contamination + 1)
    ) +
  stat_smooth(
    method = lm, data = mould_data %>%
      filter(mould_contamination > Q2_mould) %>%
      mutate(mould_contamination = mould_contamination + 1),
    colour = "#31a354"
    ) +
  stat_cor(
    data = mould_data %>%
      mutate(mould_contamination = mould_contamination + 1),
           aes(label = after_stat(rr.label)), method = "pearson",
           label.y = 200, label.x = log10(33),
           colour = "#3366FF", r.digits = 1
    ) +
  stat_cor(
    data = mould_data %>%
             filter(mould_contamination > Q2_mould) %>%
             mutate(mould_contamination = mould_contamination + 1),
           aes(label = after_stat(rr.label)), method = "pearson",
           label.y = 187.5, label.x = log10(33),
           colour = "#31a354", r.digits = 3
    ) +
  scale_x_log10(
    breaks = c(1,  5, 20, 80),
    limits = c(1, (max(mould_data$mould_contamination) + 1))
  ) +
  labs(
    y = "Observed OTU richness",
    x = "Logarithm of mould relative abundance (%)") +
  MyTheme() +
  theme(
    axis.title = element_text(size = rel(1.2)),
    axis.text = element_text(size = rel(1.2))
  )

# Contaminated sample: 65 in total
mouldy_samples <- mould_data %>%
  filter(
    mould_contamination >= 30
  ) %>%
  arrange(desc(mould_contamination)) %>%
  print(n = 100)

# Effect of mould contamination on the dataset after filtering
sample_not_mouldy <- mould_data %>%
  filter(!sample_id %in% mouldy_samples$sample_id)

# Calculate mould median value:
Q2_no_mould <- sample_not_mouldy %>%
  select(mould_contamination) %>%
  summarise(median(mould_contamination))  %>%
  as.numeric()

# Without moulds plot
plot_without_mould <- sample_not_mouldy %>%
  mutate(mould_contamination = mould_contamination + 1) %>%
  ggplot(aes(mould_contamination, richness)) +
  geom_point(size = 1, color = "black", alpha = 0.5) +
  stat_smooth(
    data = sample_not_mouldy %>%
      mutate(mould_contamination = mould_contamination + 1),
    method = loess, se = FALSE, colour = "#de2d26"
  ) +
  stat_smooth(
   data = sample_not_mouldy %>%
     mutate(mould_contamination = mould_contamination + 1),
   method = lm
    ) +
  stat_smooth(
    data = sample_not_mouldy %>%
      filter(mould_contamination > Q2_no_mould) %>%
      mutate(mould_contamination = mould_contamination + 1),
    method = lm, colour = "#31a354"
    ) +
  stat_cor(
    data = sample_not_mouldy %>%
      mutate(mould_contamination = mould_contamination + 1),
    aes(label = after_stat(rr.label)), method = "pearson",
    label.y = 200, label.x = log10(15),
    colour = "#3366FF", r.digits = 2
    ) +
  stat_cor(
    data = sample_not_mouldy %>%
      filter(mould_contamination > Q2_mould) %>%
      mutate(mould_contamination = mould_contamination + 1),
    aes(label = after_stat(rr.label)), method = "pearson",
    label.y = 187.5, label.x = log10(15),
    colour = "#31a354", r.digits = 1
    ) +
  scale_x_log10(
    breaks = c(1,  5, 20, 40),
    labels = scales::comma_format(scale = 1, big.mark = ","),
    limits = c(1, (max(sample_not_mouldy$mould_contamination) + 1))
  ) +
  labs(
    y = "Observed OTU richness",
    x = "Logarithm of mould relative abundance (%)") +
  MyTheme() +
  theme(
    axis.title = element_text(size = rel(1.2)),
    axis.text = element_text(size = rel(1.2)),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )

# Annotate plots
plot_with_mould_A <- plot_with_mould +
  annotate("text", x = 1, y = 200, label = "a", size = 5,
           fontface = "bold")
plot_without_mould_B <- plot_without_mould +
  annotate("text", x = 1, y = 200, label = "b", size = 5,
           fontface = "bold")

# Save plot
mould_summary <- patchwork::wrap_plots(plot_with_mould_A, plot_without_mould_B)
ggsave("output/mould_summary.jpg", width = 8.5, height = 4.25, dpi = 300)
ggsave("output/mould_summary.tif", width = 8.5, height = 4.25)

# Clear the environment
rm(list = ls())

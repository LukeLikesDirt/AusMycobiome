
# Difference in ITS extraction method on diversity of distinct groups

# (1) ITS method comparison ########################################################

require(tidyverse)

# Combined sample data with updated and reverted values
combined_rank_data <- data.frame(
  method = c("ITSxpress", "ITSxpress", "ITSxpress", "ITSxpress", "ITSxpress", "ITSxpress", "ITSx", "ITSx", "ITSx", "ITSx", "ITSx", "ITSx"),
  genus = c("Glomus", "Rhizophagus", "Ruhlandiella", "Sphaerosoma", "Inocybe", "Cortinarius", "Sphaerosoma", "Inocybe", "Ruhlandiella", "Cortinarius", "Glomus", "Rhizophagus"),
  absolute_abundance = c(1849794, 398741, 2611577, 3803524, 1235021, 314165, 1960561, 1423392, 1042463, 946239, 489702, 132288),
  relative_abundance = c(1.4, 0.302, 1.98, 2.88, 0.934, 0.238, 2.24, 1.63, 1.19, 1.08, 0.56, 0.151),
  absolute_richness = c(536, 231, 164, 118, 109, 89, 189, 230, 266, 220, 483, 189),
  relative_richness = c(1.68, 0.724, 0.514, 0.37, 0.341, 0.279, 0.526, 0.641, 0.741, 0.613, 1.35, 0.526),
  absolute_prevalence = c(804, 567, 699, 580, 431, 125, 353, 491, 410, 302, 514, 282),
  relative_prevalence = c(38.2, 27, 33.2, 27.6, 20.5, 5.94, 17.1, 23.8, 19.8, 14.6, 24.9, 13.6)
)

# Melting the data for faceting
melted_data <- combined_rank_data %>%
  pivot_longer(
    cols = starts_with("relative_"),
    names_to = "metric",
    values_to = "relative_value"
  ) %>%
  pivot_longer(
    cols = starts_with("absolute_"),
    names_to = "absolute_metric",
    values_to = "absolute_value"
  ) %>%
  filter(substring(metric, 9) == substring(absolute_metric, 9)) %>%
  mutate(
    metric = factor(metric, levels = c("relative_richness", "relative_abundance", "relative_prevalence")),
    genus = factor(genus, levels = c("Cortinarius", "Inocybe", "Ruhlandiella", "Sphaerosoma", "Glomus", "Rhizophagus"))
  )

# Create the faceted plot
facet_plot <- ggplot(melted_data, aes(x = genus, y = relative_value, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(
    label = scales::comma(absolute_value),
    y = 0
  ),
  position = position_dodge(0.9),
  hjust = -0.1,  # Adjust horizontal justification for better alignment
  size = 2.75
  ) +
  labs(
    x = NULL,  # Remove x-axis title
    y = "Relative Measure (%)",
    fill = "Method"
  ) +
  facet_wrap(~ metric, scales = "free_x", labeller = as_labeller(c(
    "relative_abundance" = "Abundance",
    "relative_richness" = "Richness",
    "relative_prevalence" = "Sample prevalence"
  ))) +
  theme_bw() +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(as.factor(melted_data$genus)))) +
  scale_fill_manual(values = c("ITSxpress" = "#1b9e77", "ITSx" = "#d95f02")) +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.title.x = element_text(size = rel(0.9)),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    aspect.ratio = 1
  )

# Print the plot
print(facet_plot)

# Save the plot
ggsave("../output/plots/ITS_method_comparison.png", facet_plot, width = 15.75, height = 7.5, units = "cm")
ggsave("../output/plots/ITS_method_comparison.tiff", facet_plot, width = 15.75, height = 7.5, units = "cm")

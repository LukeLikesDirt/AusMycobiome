require(emmeans)
require(performance)
require(parameters)
require(data.table)
require(tidyverse)

# Read in the taxonomy data
taxa <- fread("../output/taxonomy.txt") %>%
  filter(
    genus %in% c("Ruhlandiella", "Sphaerosoma"),
  ) %>%
  select(OTU_ID, genus) %>%
  distinct() %>%
  glimpse()

# Read in the nitrogen data
data <- fread("../output/otu_table_srs.txt") %>%
  pivot_longer(
    -OTU_ID,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  inner_join(
    taxa,
    by = "OTU_ID"
  ) %>%
  inner_join(
    fread("../output/sample_metadata.txt") %>%
      select(
        sample_id, N_total_5cm, P_total_5cm, ph_verbatim, pH_CaCl2_5cm,
        ecoregion, bioregion_code, longitude, latitude
      ),
    by = "sample_id"
  ) %>%
  # Randomly select 1 sample per georeferenced coordinates
  group_by(longitude, latitude) %>%
  sample_n(1) %>%
  ungroup() %>%
  print(.)

# How many OTUs do we have?
data %>%
  group_by(genus) %>%
  summarise(n_OTUs = n_distinct(OTU_ID))

# Data distribution
ggplot(data, aes(x = genus, y = log(N_total_5cm+1))) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_minimal()
ggplot(data, aes(x = genus, y = log(P_total_5cm+1))) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_minimal()
ggplot(data, aes(x = genus, y = ph_verbatim)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_minimal()
ggplot(data, aes(x = genus, y = pH_CaCl2_5cm)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_minimal()

# (1) N_total_5cm models ########################################################

### (1.1) Compute niche position and width ####
data_N <- data %>%
  select(OTU_ID, genus, N_total_5cm, bioregion_code, ecoregion, longitude, latitude) %>%
  na.omit() %>%
  # Take species taht occur in at least 5 samples
  group_by(OTU_ID) %>%
  filter(n() >= 5) %>%
  # Compute niche position and breadth
  summarise(
    genus = first(genus),
    niche_breadth = sd(N_total_5cm, na.rm = TRUE),
    niche_position = mean(N_total_5cm, na.rm = TRUE)
  ) %>%
  ungroup()

### (1.2) Fit the models ####
lm_pos_N <- lm(log(niche_position) ~ genus, data = data_N)
lm_breadth_N <- lm(log(niche_breadth) ~ genus, data = data_N)

# Check model performance
check_model(lm_pos_N)
check_model(lm_breadth_N)

# Model parameters
parameters(lm_pos_N)
parameters(lm_breadth_N)

# Extract the marginal means
emmeans_pos_N <- emmeans(lm_pos_N, specs = ~genus) %>%
  as_tibble() %>%
  mutate(across(where(is.numeric), exp))
emmargins_breadth_N <- emmeans(lm_breadth_N, specs = ~genus) %>%
  as_tibble() %>%
  mutate(across(where(is.numeric), exp))

### (1.3) Create the plots ####
plot_pos_N <- data_N %>%
  mutate(
    `Total nitrogen (%)` = case_when(
      genus == "Ruhlandiella" ~ marginal_means_N$emmean[1],
      genus == "Sphaerosoma" ~ marginal_means_N$emmean[2]
    ),
    upper.CL = case_when(
      genus == "Ruhlandiella" ~ marginal_means_N$upper.CL[1],
      genus == "Sphaerosoma" ~ marginal_means_N$upper.CL[2]
    ),
    lower.CL = case_when(
      genus == "Ruhlandiella" ~ marginal_means_N$lower.CL[1],
      genus == "Sphaerosoma" ~ marginal_means_N$lower.CL[2]
    )
  ) %>%
  ggplot() +
  # Raw data points
  geom_beeswarm(
    aes(
      x = genus,
      y = niche_position,
      colour = genus
    ),
    alpha = 0.25
  ) +
  # Marginal means error bars
  geom_errorbar(
    aes(
      x = genus,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.1,
    colour = "black"
    #position = position_nudge(x = -0.2)  # Align with marginal mean points
  ) +
  # Estimated means points
  geom_point(
    aes(x = genus, y = `Total nitrogen (%)`, fill = genus),
    #position = position_nudge(x = -0.2),
    size = 3,
    shape = 21,
    colour = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
    aspect.ratio = 1
  ) +
  labs(
    y = "Niche position",
    title = "Total nitrogen (%)"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "a",
    hjust = -0.5,
    vjust = 1.5,
    fontface = "bold",
    size = 6
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(beta)[1] == -0.77 ~ "[" ~ -1.05 ~ "," ~ -0.49 ~ "]"),
    hjust = 1.05, 
    vjust = 1.5, 
    size = 2.5
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(t) ~ (47) == -5.50),
    hjust = 1.05, 
    vjust = 2.5, 
    size = 2.5
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(p) <= ".001"),
    hjust = 1.05, 
    vjust = 4.75, 
    size = 2.5
  )

plot_breadth_N <- data_N %>%
  mutate(
    `Total nitrogen (%)` = case_when(
      genus == "Ruhlandiella" ~ emmargins_breadth_N$emmean[1],
      genus == "Sphaerosoma" ~ emmargins_breadth_N$emmean[2]
    ),
    upper.CL = case_when(
      genus == "Ruhlandiella" ~ emmargins_breadth_N$upper.CL[1],
      genus == "Sphaerosoma" ~ emmargins_breadth_N$upper.CL[2]
    ),
    lower.CL = case_when(
      genus == "Ruhlandiella" ~ emmargins_breadth_N$lower.CL[1],
      genus == "Sphaerosoma" ~ emmargins_breadth_N$lower.CL[2]
    )
  ) %>%
  ggplot() +
  # Raw data points
  geom_beeswarm(
    aes(
      x = genus,
      y = niche_breadth,
      colour = genus
    ),
    alpha = 0.25
  ) +
  # Marginal means error bars
  geom_errorbar(
    aes(
      x = genus,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.1,
    colour = "black"
    #position = position_nudge(x = -0.2)  # Align with marginal mean points
  ) +
  # Estimated means points
  geom_point(
    aes(x = genus, y = `Total nitrogen (%)`, fill = genus),
    #position = position_nudge(x = -0.2),
    size = 3,
    shape = 21,
    colour = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "italic"),
    aspect.ratio = 1
  ) +
  labs(
    y = "Nitrogen breadth"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "c",
    hjust = -0.5,
    vjust = 1.5,
    fontface = "bold",
    size = 6
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(beta)[1] == -0.88 ~ "[" ~ -1.33 ~ "," ~ -0.42 ~ "]"),
    hjust = 1.05, 
    vjust = 1.5, 
    size = 2.5
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(t) ~ (47) == -3.87),
    hjust = 1.05, 
    vjust = 2.5, 
    size = 2.5
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(p) <= ".001"),
    hjust = 1.05, 
    vjust = 4.75, 
    size = 2.5
  )
  
# (2) P_total_5cm models ########################################################

### (2.1) Compute niche position and width ####
data_P <- data %>%
  select(OTU_ID, genus, P_total_5cm, bioregion_code, ecoregion, longitude, latitude) %>%
  na.omit() %>%
  # Take species taht occur in at least 5 samples
  group_by(OTU_ID) %>%
  filter(n() >= 5) %>%
  # Compute niche position and breadth
  summarise(
    genus = first(genus),
    niche_breadth = sd(P_total_5cm, na.rm = TRUE),
    niche_position = mean(P_total_5cm, na.rm = TRUE)
  ) %>%
  ungroup()

### (2.2) Fit the models ####
lm_pos_P <- lm(log(niche_position) ~ genus, data = data_P)
lm_breadth_P <- lm(log(niche_breadth) ~ genus, data = data_P)

# Check model performance
check_model(lm_pos_P)
check_model(lm_breadth_P)

# Model parameters
parameters(lm_pos_P)
parameters(lm_breadth_P)  

# Extract the marginal means
emmeans_pos_P <- emmeans(lm_pos_P, specs = ~genus) %>%
  as_tibble() %>%
  mutate(across(where(is.numeric), exp))
emmargins_breadth_P <- emmeans(lm_breadth_P, specs = ~genus) %>%  
  as_tibble() %>%
  mutate(across(where(is.numeric), exp))

### (2.3) Create the plots ####
plot_pos_P <- data_P %>%
  mutate(
    `Total phosphorus (%)` = case_when(
      genus == "Ruhlandiella" ~ emmeans_pos_P$emmean[1],
      genus == "Sphaerosoma" ~ emmeans_pos_P$emmean[2]
    ),
    upper.CL = case_when(
      genus == "Ruhlandiella" ~ emmeans_pos_P$upper.CL[1],
      genus == "Sphaerosoma" ~ emmeans_pos_P$upper.CL[2]
    ),
    lower.CL = case_when(
      genus == "Ruhlandiella" ~ emmeans_pos_P$lower.CL[1],
      genus == "Sphaerosoma" ~ emmeans_pos_P$lower.CL[2]
    )
  ) %>%
  ggplot() +
  # Raw data points
  geom_beeswarm(
    aes(
      x = genus,
      y = niche_position,
      colour = genus
    ),
    alpha = 0.25
  ) +
  # Marginal means error bars
  geom_errorbar(
    aes(
      x = genus,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.1,
    colour = "black"
    #position = position_nudge(x = -0.2)  # Align with marginal mean points
  ) +
  # Estimated means points
  geom_point(
    aes(x = genus, y = `Total phosphorus (%)`, fill = genus),
    #position = position_nudge(x = -0.2),
    size = 3,
    shape = 21,
    colour = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
    aspect.ratio = 1
  ) +
  labs(
    y = "Niche position",
    title = "Total phosphorus (%)"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "b",
    hjust = -0.5,
    vjust = 1.5,
    fontface = "bold",
    size = 6
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(beta)[1] == -0.32 ~ "[" ~ -0.45 ~ "," ~ -0.19 ~ "]"),
    hjust = 1.05, 
    vjust = 1.5, 
    size = 2.5
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(t) ~ (47) == -4.84),
    hjust = 1.05, 
    vjust = 2.5, 
    size = 2.5
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(p) <= ".001"),
    hjust = 1.05, 
    vjust = 4.75, 
    size = 2.5
  )

plot_breadth_P <- data_P %>%
  mutate(
    `Total phosphorus (%)` = case_when(
      genus == "Ruhlandiella" ~ emmargins_breadth_P$emmean[1],
      genus == "Sphaerosoma" ~ emmargins_breadth_P$emmean[2]
    ),
    upper.CL = case_when(
      genus == "Ruhlandiella" ~ emmargins_breadth_P$upper.CL[1],
      genus == "Sphaerosoma" ~ emmargins_breadth_P$upper.CL[2]
    ),
    lower.CL = case_when(
      genus == "Ruhlandiella" ~ emmargins_breadth_P$lower.CL[1],
      genus == "Sphaerosoma" ~ emmargins_breadth_P$lower.CL[2]
    )
  ) %>%
  ggplot() +
  # Raw data points
  geom_beeswarm(
    aes(
      x = genus,
      y = niche_breadth,
      colour = genus
    ),
    alpha = 0.25
  ) +
  # Marginal means error bars
  geom_errorbar(
    aes(
      x = genus,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.1,
    colour = "black"
    #position = position_nudge(x = -0.2)  # Align with marginal mean points
  ) +
  # Estimated means points
  geom_point(
    aes(x = genus, y = `Total phosphorus (%)`, fill = genus),
    #position = position_nudge(x = -0.2),
    size = 3,
    shape = 21,
    colour = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_text(face = "italic"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    aspect.ratio = 1
  ) +
  labs(
    y = "Niche breadth"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "d",
    hjust = -0.5,
    vjust = 1.5,
    fontface = "bold",
    size = 6
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(beta)[1] == -0.45 ~ "[" ~ -0.73 ~ "," ~ -0.16 ~ "]"),
    hjust = 1.05, 
    vjust = 1.5, 
    size = 2.5
    ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(t) ~ (47) == -3.18),
    hjust = 1.05, 
    vjust = 2.5, 
    size = 2.5
    ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = expression(italic(p) == ".003"),
    hjust = 1.05, 
    vjust = 4.75, 
    size = 2.5
      )
plot_breadth_P
  

# (3) Wrap and save plots ######################################################

# Wrap the plots
patchwork::wrap_plots(
  plot_pos_N + plot_pos_P + plot_breadth_N + plot_breadth_P
)

# Save the plots
ggsave(
  filename = "../output/niche_analysis.png",
  width = 15.75,
  height = 15.75,
  units = "cm"
)
ggsave(
  filename = "../output/niche_analysis.tiff",
  width = 15.75,
  height = 15.75,
  units = "cm"
)

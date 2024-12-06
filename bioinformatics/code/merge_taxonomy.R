
taxa_dir <- "../data/taxonomy/"
# Load packages
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(tidyverse))

# Step 1: Species annotations
species <- fread(paste0(taxa_dir, "species_classification.txt"))
species_ID <- species$ID

# Step 2: Genus annotations
genus <- fread(paste0(taxa_dir, "genus_classification.txt")) %>%
  filter(!ID %in% species_ID)
genus_ID <- genus$ID

# Step 3: Family and order annotations
family <- fread(paste0(taxa_dir, "family_classification.txt")) %>%
  filter(!ID %in% c(species_ID, genus_ID))
family_ID <- family$ID

# Step 4: Order annotations
order <- fread(paste0(taxa_dir, "order_classification.txt")) %>%
  filter(!ID %in% c(species_ID, genus_ID, family_ID))
order_ID <- order$ID

# Sterp 5: Class annotations
class <- fread(paste0(taxa_dir, "class_classification.txt")) %>%
  filter(!ID %in% c(species_ID, genus_ID, family_ID, order_ID))
class_ID <- class$ID

# Step 6: Phylum and kingdom annotations
all <- fread(paste0(taxa_dir, "phylum_classification.txt")) %>%
  filter(!ID %in% c(
      species_ID, genus_ID, family_ID, order_ID, class_ID
    ),
    # Remove the "unclassified" phyla with cutofff lower than min cutoff
    score >= fread(paste0(taxa_dir, "phylum_classification.txt")) %>% 
      select(cutoff) %>%
      filter(cutoff != "N/A") %>%
      pull(.) %>%
      unique(.) %>%
      as.numeric(.)
  )

# Build the taxonomy table
taxonomy <- rbind(
  species, genus, family, order, class, all
) %>%
  mutate(
    OTU_ID = ID,
    confidence = as.numeric(confidence),
    cutoff = as.numeric(cutoff),
    score = as.numeric(score)
  ) %>%
  select(OTU_ID, everything(), -ID) %>%
  arrange(
    desc(confidence),
    desc(score)
  )

# Check if any duplicate OTU_IDs
duplicates <- taxonomy %>%
  group_by(OTU_ID) %>%
  filter(n() > 1) %>%
  ungroup()

# Print duplicates
print(duplicates, n = Inf)

# Exit with an error if duplicates are found
if (nrow(duplicates) > 0) {
  stop("Error: Duplicate OTU_IDs found.")
}

# Save the taxonomy table
taxonomy %>%
  # Replace spaces with underscores in all columns
  mutate_all(~ gsub(" ", "_", .)) %>%
  fwrite(., paste0(taxa_dir, "/taxa_ASVs.txt"), sep = "\t")


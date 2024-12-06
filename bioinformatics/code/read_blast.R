
# Define the read_blast function
read_blast <- function(data, minlen = 50) {
  data %>%
    mutate(
      # Replace GS01 references
      phylum = str_replace(phylum, "GS01_phy_Incertae_sedis", "GS01"),
      class = str_replace(class, "GS01_phy_Incertae_sedis", "GS01"),
      order = str_replace(order, "GS01_phy_Incertae_sedis", "GS01"),
      family = str_replace(family, "GS01_phy_Incertae_sedis", "GS01"),
      genus = str_replace(genus, "GS01_phy_Incertae_sedis", "GS01"),
      # Adjust for Sphaerosoma genus: Previously Sphaerosoma fungi were classified as a to a bettel genus
      kingdom = case_when(
        genus == "Sphaerosoma" ~ "Fungi",
        TRUE ~ kingdom
      ),
      phylum = case_when(
        genus == "Sphaerosoma" ~ "Ascomycota",
        TRUE ~ phylum
      ),
      class = case_when(
        genus == "Sphaerosoma" ~ "Pezizomycetes",
        TRUE ~ class
      ),
      order = case_when(
        genus == "Sphaerosoma" ~ "Pezizales",
        TRUE ~ order
      ),
      family = case_when(
        genus == "Sphaerosoma" ~ "Pyronemataceae",
        TRUE ~ family
      ),
      # Mutate kingdom for cutoff categories
      kingdom = case_when(
        phylum %in% c("Ascomycota", "Basidiomycota", "Entorrhizomycota") ~ "dikarya_fungi",
        phylum %in% c("Glomeromycota", "Mortierellomycota", "Mucoromycota", "Calcarisporiellomycota", "Basidiobolomycota", "Entomophthoromycota", "Kickxellomycota", "Zoopagomycota") ~ "terrestrial_basal_fungi",
        phylum %in% c("Aphelidiomycota", "Blastocladiomycota", "Chytridiomycota", "Monoblepharomycota", "Neocallimastigomycota", "Olpidiomycota", "Rozellomycota", "Sanchytriomycota", "GS01") ~ "zoosporic_basal_fungi",
        TRUE ~ kingdom
      )
    ) %>%
    # Compute additional metrics
    mutate(
      subject_coverage = (length / slen), # Normalise subject_coverage to a fraction
      query_coverage = (length / qlen), # Normalise query_coverage to a fraction
      score = pident / 100, # Normalise pident to a fraction
      # Adjust score for short alignments
      score = case_when(
        length < minlen ~ (score * length) / minlen,
        TRUE ~ score
      )
    ) %>%
    select(OTU_ID, reference_ID, kingdom, phylum, class, order, family, genus, species, species_hypothesis, score, subject_coverage, query_coverage)
}

###############################################################################
# Supplementary Host Metadata Table (16S + ITS combined)
# One row per host species; columns include sample counts per dataset,
# number of paired individuals (present in both 16S and ITS), and
# host ecological/biological metadata.
# Alexander Rurik
###############################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
})

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate

# ========================== BASIC PARAMETERS ==================================
set.seed(52325)

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_tables    <- here("output", "tables")

# ========================== LOAD DATA =========================================
ps_16S <- readRDS(file.path(dir_processed, "16S_absolute_all_final.rds"))
ps_ITS <- readRDS(file.path(dir_processed, "ITS_absolute_all_final.rds"))

# ========================== EXTRACT METADATA ==================================
metadata_16S <- data.frame(sample_data(ps_16S))
metadata_ITS <- data.frame(sample_data(ps_ITS))

# ========================== STANDARDIZE SITE NAMES ===========================
# Shorten verbose captive site names for cleaner table display
clean_sites <- function(metadata) {
  metadata %>%
    mutate(site = as.character(site),
           site = case_when(
             site == "Amphibian Foundation captive collection" ~ "Amphibian Foundation",
             site == "Captive population at SELU"             ~ "SELU captive population",
             TRUE                                             ~ site
           ),
           # For captive samples use institution name; for wild samples use EPA Level III ecoregion
           Ecoregion_or_Site = ifelse(env_broad_scale == "Captive", site, ecoregion_III))
}

metadata_16S <- clean_sites(metadata_16S)
metadata_ITS <- clean_sites(metadata_ITS)

# ========================== BUILD PER-DATASET SUMMARIES =======================
# Summarize by host species for each dataset independently before joining,
# so that sample counts, sample types, and locations are accurate per dataset

summarize_dataset <- function(metadata, count_col) {
  metadata %>%
    group_by(host_taxon) %>%
    summarise(
      !!count_col                               := n(),
      `Wild/Captive`                             = paste(sort(unique(env_broad_scale)), collapse = ", "),
      `Sample Type`                              = paste(sort(unique(env_medium)),       collapse = ", "),
      `EPA Level III Ecoregion / Captive Site`   = paste(sort(unique(Ecoregion_or_Site)), collapse = ", "),
      `Primary Habitat`                          = paste(sort(unique(animal_ecomode)),   collapse = ", "),
      `Host Order`                               = paste(sort(unique(Clade_Order)),      collapse = ", "),
      `Host Family`                              = paste(sort(unique(Family)),            collapse = ", "),
      `Diet`                                     = paste(sort(unique(Diet)),              collapse = ", "),
      .groups = "drop"
    ) %>%
    rename(`Host Species` = host_taxon)
}

host_summary_16S <- summarize_dataset(metadata_16S, "# Samples (16S)")
host_summary_ITS <- summarize_dataset(metadata_ITS, "# Samples (ITS)")

# ========================== CALCULATE PAIRED INDIVIDUALS ======================
# Paired individuals = unique animal_number values present in BOTH datasets.
# Reported per host species.

animals_16S    <- unique(metadata_16S$animal_number)
animals_ITS    <- unique(metadata_ITS$animal_number)
paired_animals <- intersect(animals_16S, animals_ITS)

cat("Animals in 16S:", length(animals_16S), "\n")
cat("Animals in ITS:", length(animals_ITS), "\n")
cat("Animals in both datasets (paired):", length(paired_animals), "\n\n")

# Count paired individuals per host species
paired_counts <- metadata_16S %>%
  filter(animal_number %in% paired_animals) %>%
  group_by(host_taxon) %>%
  summarise(`# Paired Individuals` = n_distinct(animal_number), .groups = "drop") %>%
  rename(`Host Species` = host_taxon)

# ========================== JOIN INTO COMBINED TABLE ==========================
# Full join so species present in only one dataset are still represented.
# Shared metadata columns (Wild/Captive, Sample Type, etc.) are resolved by
# combining values from both datasets where they differ.

host_summary_combined <- full_join(
  host_summary_16S,
  host_summary_ITS,
  by = "Host Species",
  suffix = c(".16S", ".ITS")
) %>%
  # For shared metadata columns, merge values from both datasets,
  # deduplicating and sorting to handle cases where 16S and ITS differ
  # (e.g., a species with both wild and captive samples in one dataset only)
  mutate(
    `Wild/Captive` = map2_chr(`Wild/Captive.16S`, `Wild/Captive.ITS`,
                              ~ paste(sort(unique(c(str_split(.x, ", ")[[1]], str_split(.y, ", ")[[1]]))), collapse = ", ")),
    `Sample Type` = map2_chr(`Sample Type.16S`, `Sample Type.ITS`,
                             ~ paste(sort(unique(c(str_split(.x, ", ")[[1]], str_split(.y, ", ")[[1]]))), collapse = ", ")),
    `EPA Level III Ecoregion / Captive Site` = map2_chr(
      `EPA Level III Ecoregion / Captive Site.16S`,
      `EPA Level III Ecoregion / Captive Site.ITS`,
      ~ paste(sort(unique(c(str_split(.x, ", ")[[1]], str_split(.y, ", ")[[1]]))), collapse = ", ")),
    `Primary Habitat` = coalesce(`Primary Habitat.16S`, `Primary Habitat.ITS`),
    `Host Order`   = coalesce(`Host Order.16S`,   `Host Order.ITS`),
    `Host Family`  = coalesce(`Host Family.16S`,  `Host Family.ITS`),
    `Diet`         = coalesce(`Diet.16S`,          `Diet.ITS`)
  ) %>%
  # Replace NAs in sample count columns with 0 (species absent from one dataset)
  mutate(
    `# Samples (16S)` = replace_na(`# Samples (16S)`, 0),
    `# Samples (ITS)` = replace_na(`# Samples (ITS)`, 0)
  ) %>%
  # Add paired individual counts (0 if species has no paired samples)
  left_join(paired_counts, by = "Host Species") %>%
  mutate(`# Paired Individuals` = replace_na(`# Paired Individuals`, 0)) %>%
  # Select and order final columns — sample counts immediately after Host Species
  select(
    `Host Species`,
    `# Samples (16S)`,
    `# Samples (ITS)`,
    `# Paired Individuals`,
    `Host Order`,
    `Host Family`,
    `Diet`,
    `Primary Habitat`,
    `Wild/Captive`,
    `Sample Type`,
    `EPA Level III Ecoregion / Captive Site`
  ) %>%
  # Sort by Order → Family → Species for readability
  arrange(`Host Order`, `Host Family`, `Host Species`)

# ========================== SUMMARY STATS =====================================
cat("========== DATASET SUMMARY ==========\n")
cat("16S: ", nrow(metadata_16S), "samples |",
    n_distinct(metadata_16S$host_taxon), "species |",
    n_distinct(metadata_16S$Family), "families\n")
cat("ITS: ", nrow(metadata_ITS), "samples |",
    n_distinct(metadata_ITS$host_taxon), "species |",
    n_distinct(metadata_ITS$Family), "families\n")
cat("Combined:", n_distinct(c(metadata_16S$host_taxon, metadata_ITS$host_taxon)), "unique species |",
    n_distinct(c(metadata_16S$Family, metadata_ITS$Family)), "unique families\n")
cat("Paired individuals (both datasets):", length(paired_animals), "\n\n")

cat("Wild/Captive breakdown:\n")
cat("  16S — Wild:", sum(metadata_16S$env_broad_scale == "Wild"),
    sprintf("(%.1f%%)", 100 * mean(metadata_16S$env_broad_scale == "Wild")),
    "| Captive:", sum(metadata_16S$env_broad_scale == "Captive"),
    sprintf("(%.1f%%)\n", 100 * mean(metadata_16S$env_broad_scale == "Captive")))
cat("  ITS — Wild:", sum(metadata_ITS$env_broad_scale == "Wild"),
    sprintf("(%.1f%%)", 100 * mean(metadata_ITS$env_broad_scale == "Wild")),
    "| Captive:", sum(metadata_ITS$env_broad_scale == "Captive"),
    sprintf("(%.1f%%)\n\n", 100 * mean(metadata_ITS$env_broad_scale == "Captive")))

cat("Ecoregions represented:\n")
cat("  16S:", n_distinct(metadata_16S$Ecoregion_or_Site), "\n")
cat("  ITS:", n_distinct(metadata_ITS$Ecoregion_or_Site), "\n")
cat("  Combined:", n_distinct(c(metadata_16S$Ecoregion_or_Site, metadata_ITS$Ecoregion_or_Site)), "\n\n")

cat("Table dimensions:", nrow(host_summary_combined), "species rows x",
    ncol(host_summary_combined), "columns\n")
cat("Species with paired samples:",
    sum(host_summary_combined$`# Paired Individuals` > 0), "\n")
cat("Species with 16S only:",
    sum(host_summary_combined$`# Samples (16S)` > 0 & host_summary_combined$`# Samples (ITS)` == 0), "\n")
cat("Species with ITS only:",
    sum(host_summary_combined$`# Samples (ITS)` > 0 & host_summary_combined$`# Samples (16S)` == 0), "\n\n")

# ========================== EXPORT ============================================
outfile <- file.path(dir_tables, "combined_16S_ITS_supp_host_summary.csv")
write.csv(host_summary_combined, file = outfile, row.names = FALSE)
cat("Table saved:", outfile, "\n")

###############################################################################
###############################################################################
###############################################################################
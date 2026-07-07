###############################################################################
# Phyloseq Metadata Cleaning — 16S & ITS
#
# Applies standardization and corrections to sample metadata in the
# post-decontam phyloseq objects before downstream analyses. Changes made:
#
#   1. Initial env_medium text standardization (trim + recode raw variants)
#   2. Correct mislabeled env_medium: samples with spiked_volume == 1 must
#      be labeled "Cloacal swab" (key correction — see notes below)
#   3. Standardize env_broad_scale to "Wild" / "Captive"
#   4. Re-run env_medium standardization (catch-all pass after Correction 2)
#   5. Correct animal_ecomode entries for specific taxa
#   6. Standardize ecoregion_III naming ("Blue Ridge Mountains" → "Blue Ridge")
#   7. Update host species/genus names (Hyla → Dryophytes; Geochelone → Aldabrachelys)
#   8. Fix Zoo Knoxville / Knoxville Zoo site name discrepancy
#   9. Fix missing GPS coordinates for two ITS samples
#
# Alexander Rurik
###############################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
})

# Explicitly prioritize dplyr functions
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
# Load post-decontam phyloseq objects produced by the decontam scripts (Script 02)
ps_16S     <- readRDS(file.path(dir_processed, "16S_physeq_absolute_decontam_output.rds"))
ps_16S_all <- readRDS(file.path(dir_processed, "16S_physeq_absolute_PRE_depth_filter.rds"))

ps_ITS     <- readRDS(file.path(dir_processed, "ITS_physeq_absolute_decontam_output.rds"))
ps_ITS_all <- readRDS(file.path(dir_processed, "ITS_physeq_absolute_PRE_depth_filter.rds"))

# ========================== EXTRACT METADATA ==================================
metadata_16S <- data.frame(sample_data(ps_16S))
metadata_16S <- tibble::rownames_to_column(metadata_16S, var = "sample_name")

metadata_ITS <- data.frame(sample_data(ps_ITS))
metadata_ITS <- tibble::rownames_to_column(metadata_ITS, var = "sample_name")

cat("Metadata extracted:\n")
cat("  16S samples:", nrow(metadata_16S), "\n")
cat("  ITS samples:", nrow(metadata_ITS), "\n\n")

# Extract metadata for the ALL samples phyloseq object
metadata_16S_all <- data.frame(sample_data(ps_16S_all))
metadata_16S_all <- tibble::rownames_to_column(metadata_16S_all, var = "sample_name")

metadata_ITS_all <- data.frame(sample_data(ps_ITS_all))
metadata_ITS_all <- tibble::rownames_to_column(metadata_ITS_all, var = "sample_name")

cat("Metadata extracted:\n")
cat("  16S samples:", nrow(metadata_16S), "| 16S all-species:", nrow(metadata_16S_all), "\n")
cat("  ITS samples:", nrow(metadata_ITS), "| ITS all-species:", nrow(metadata_ITS_all), "\n\n")

# Export raw metadata CSVs for manual inspection if needed (uncomment to use)
# write.csv(metadata_16S, file = file.path(dir_tables, "16S_metadata_phyloseq_TO_CLEAN.csv"), row.names = FALSE)
# write.csv(metadata_ITS, file = file.path(dir_tables, "ITS_metadata_phyloseq_TO_CLEAN.csv"), row.names = FALSE)

##------------------------------------------------------------------------------
## Correction 1: Initial env_medium text standardization ####
##------------------------------------------------------------------------------
# Trim whitespace and recode raw variant spellings to clean standard labels
# before any downstream logic runs against env_medium values. This ensures
# the spiked_volume check in Correction 2 compares against consistent strings.

metadata_16S <- metadata_16S %>%
  mutate(env_medium = env_medium %>%
           str_trim() %>%
           recode(
             "cloacal swab - fecal" = "Cloacal swab",
             "cloacal swab"         = "Cloacal swab",
             "stomach"              = "Lower GI",
             "lower GI"             = "Lower GI",
             "fecal"                = "Fecal"
           ))

metadata_ITS <- metadata_ITS %>%
  mutate(env_medium = env_medium %>%
           str_trim() %>%
           recode(
             "cloacal swab - fecal" = "Cloacal swab",
             "cloacal swab"         = "Cloacal swab",
             "stomach"              = "Lower GI",
             "lower GI"             = "Lower GI",
             "fecal"                = "Fecal"
           ))

cat("--- Correction 1: Initial env_medium standardization ---\n")
cat("  16S:\n"); print(table(metadata_16S$env_medium, useNA = "ifany"))
cat("  ITS:\n"); print(table(metadata_ITS$env_medium, useNA = "ifany"))
cat("\n")

##------------------------------------------------------------------------------
## Correction 2: Fix mislabeled env_medium for cloacal swab samples ####
##------------------------------------------------------------------------------
# CRITICAL CORRECTION: Samples with spiked_volume == 1 are cloacal swabs.
# Some of these were previously mislabeled as "Fecal" or other sample types.
# Runs after Correction 1 so all comparisons are against standardized strings
# (i.e., "Cloacal swab" with correct capitalization is already the standard
# after Correction 1 — so any remaining mismatch here is a true labeling error,
# not a capitalization artifact).
#
# Background: This error was identified during metadata review (April 2026)
# and traced to mislabeled entries in the May 2023 turtle/snake field trip
# and July 2024 Louisiana alligator field trip. The master database has been
# corrected; this step ensures the phyloseq metadata reflects those corrections.

# Snapshot env_medium BEFORE correction so we can identify what actually changed.
# This is taken after Correction 1, so capitalization is already standardized —
# any sample not labeled "Cloacal swab" here is a genuine mislabeling error.
pre_correction_16S <- metadata_16S %>%
  select(sample_name, env_medium_before = env_medium)

pre_correction_ITS <- metadata_ITS %>%
  select(sample_name, env_medium_before = env_medium)

# --- Apply correction: all spiked_volume == 1 samples → "Cloacal swab" ---
metadata_16S <- metadata_16S %>%
  mutate(env_medium = if_else(spiked_volume == 1, "Cloacal swab", env_medium))

metadata_ITS <- metadata_ITS %>%
  mutate(env_medium = if_else(spiked_volume == 1, "Cloacal swab", env_medium))

# --- Verification ---
verify_correction <- function(metadata, pre_snap, dataset_label) {
  n_spike1       <- sum(metadata$spiked_volume == 1, na.rm = TRUE)
  n_now_cloacal  <- sum(metadata$spiked_volume == 1 &
                          metadata$env_medium == "Cloacal swab", na.rm = TRUE)
  n_changed      <- metadata %>%
    left_join(pre_snap, by = "sample_name") %>%
    filter(spiked_volume == 1, env_medium_before != "Cloacal swab") %>%
    nrow()
  
  cat("--- Correction 2:", dataset_label, "---\n")
  cat("  spiked_volume == 1 samples:         ", n_spike1, "\n")
  cat("  Already correctly labeled before:   ", n_spike1 - n_changed, "\n")
  cat("  Genuinely corrected (mislabeled):   ", n_changed, "\n")
  cat("  Now labeled 'Cloacal swab':         ", n_now_cloacal, "\n")
  if (n_spike1 == n_now_cloacal) {
    cat("  PASS\n\n")
  } else {
    cat("  FAIL: mismatch — investigate before proceeding\n\n")
  }
}

verify_correction(metadata_16S, pre_correction_16S, "16S")
verify_correction(metadata_ITS, pre_correction_ITS, "ITS")

cat("  16S env_medium distribution after correction:\n")
print(table(metadata_16S$env_medium, useNA = "ifany"))
cat("\n  ITS env_medium distribution after correction:\n")
print(table(metadata_ITS$env_medium, useNA = "ifany"))
cat("\n")

##------------------------------------------------------------------------------
## Export: Genuinely mislabeled sample list (Correction 2) ####
##------------------------------------------------------------------------------
# Captures only samples that were CHANGED by Correction 2 — i.e., had
# spiked_volume == 1 but were NOT already labeled "Cloacal swab" after
# Correction 1. Samples that were already correctly labeled are excluded.
# The pre-correction label (env_medium_before) reflects the state after
# Correction 1 capitalization standardization, so any label other than
# "Cloacal swab" here is a confirmed mislabeling error, not a formatting issue.

build_mislabeled_df <- function(metadata, pre_snap, dataset_label) {
  metadata %>%
    left_join(pre_snap, by = "sample_name") %>%
    filter(
      spiked_volume == 1,
      env_medium_before != "Cloacal swab"   # was genuinely mislabeled
    ) %>%
    mutate(
      label_before = env_medium_before,
      label_after  = env_medium,             # always "Cloacal swab" post-correction
      dataset      = dataset_label
    ) %>%
    select(sample_name, host_taxon, host_genus, Clade_Order,
           env_broad_scale, site, spiked_volume,
           label_before, label_after, dataset)
}

mislabeled_16S <- build_mislabeled_df(metadata_16S, pre_correction_16S, "bacteria_16S")
mislabeled_ITS <- build_mislabeled_df(metadata_ITS, pre_correction_ITS, "fungi_ITS")

mislabeled_all <- bind_rows(mislabeled_16S, mislabeled_ITS) %>%
  arrange(dataset, Clade_Order, host_taxon, sample_name)

# Species-level summary of genuinely corrected samples
mislabeled_species_summary <- mislabeled_all %>%
  group_by(dataset, Clade_Order, host_taxon) %>%
  summarise(
    n_samples_corrected = n(),
    label_before        = paste(sort(unique(label_before)), collapse = "; "),
    sites               = paste(sort(unique(site)),         collapse = "; "),
    env_broad_scale     = paste(sort(unique(env_broad_scale)), collapse = "; "),
    .groups             = "drop"
  ) %>%
  arrange(dataset, Clade_Order, host_taxon)

cat("--- Mislabeled sample export (genuinely corrected only) ---\n")
cat("  Samples genuinely corrected — 16S:", nrow(mislabeled_16S),
    "| ITS:", nrow(mislabeled_ITS), "\n")
cat("  Unique species affected    — 16S:", n_distinct(mislabeled_16S$host_taxon),
    "| ITS:", n_distinct(mislabeled_ITS$host_taxon), "\n\n")
cat("  Species summary:\n")
print(mislabeled_species_summary, n = Inf)

#write.csv(mislabeled_all, file.path(dir_tables, "mislabeled_cloacal_samples_corrected.csv"), row.names = FALSE)
#write.csv(mislabeled_species_summary, file.path(dir_tables, "mislabeled_cloacal_species_summary.csv"), row.names = FALSE)
cat("\nSaved mislabeled sample and species summary CSVs.\n\n")

##------------------------------------------------------------------------------
## Correction 3: Standardize env_broad_scale ("Wild" / "Captive") ####
##------------------------------------------------------------------------------
metadata_16S <- metadata_16S %>%
  mutate(env_broad_scale = if_else(tolower(env_broad_scale) == "wild", "Wild", "Captive"))

metadata_ITS <- metadata_ITS %>%
  mutate(env_broad_scale = if_else(tolower(env_broad_scale) == "wild", "Wild", "Captive"))

cat("--- Correction 2: env_broad_scale standardization ---\n")
cat("  16S:", table(metadata_16S$env_broad_scale), "\n")
cat("  ITS:", table(metadata_ITS$env_broad_scale), "\n\n")

##------------------------------------------------------------------------------
## Correction 4: Re-run env_medium standardization (catch-all pass) ####
##------------------------------------------------------------------------------
# Re-applies the same recode after the spiked_volume correction in case any
# values were modified or introduced by Correction 2. Also serves as a
# verification that no non-standard values remain.

recode_env_medium <- function(x) {
  x %>%
    str_trim() %>%
    recode(
      "cloacal swab - fecal" = "Cloacal swab",
      "cloacal swab"         = "Cloacal swab",
      "stomach"              = "Lower GI",
      "lower GI"             = "Lower GI",
      "fecal"                = "Fecal"
    )
}

metadata_16S <- metadata_16S %>% mutate(env_medium = recode_env_medium(env_medium))
metadata_ITS <- metadata_ITS %>% mutate(env_medium = recode_env_medium(env_medium))

cat("--- Correction 4: env_medium catch-all standardization ---\n")
cat("  16S final distribution:\n"); print(table(metadata_16S$env_medium, useNA = "ifany"))
cat("  ITS final distribution:\n"); print(table(metadata_ITS$env_medium, useNA = "ifany"))
cat("\n")

##------------------------------------------------------------------------------
## Correction 5: Fix animal_ecomode entries ####
##------------------------------------------------------------------------------
# Corrections based on natural history literature and expert judgment.
# Rationale for each change documented inline.

fix_ecomode <- function(metadata) {
  metadata %>%
    mutate(animal_ecomode = as.character(animal_ecomode)) %>%
    mutate(animal_ecomode = case_when(
      # Agkistrodon piscivorus: semi-aquatic, not fully aquatic
      host_taxon == "Agkistrodon piscivorus"                                       ~ "Aquatic-Terrestrial",
      # Osteolaemus tetraspis + Alligator mississippiensis: strongly aquatic but uses terrestrial burrows
      host_taxon %in% c("Osteolaemus tetraspis", "Alligator mississippiensis")     ~ "Aquatic-Terrestrial",
      # Anolis carolinensis + A. distichus: trunk-ground generalists depending on population
      host_taxon %in% c("Anolis carolinensis", "Anolis distichus")                 ~ "Arboreal-Terrestrial",
      # Bombina orientalis: semi-aquatic rather than terrestrial
      host_taxon == "Bombina orientalis"                                            ~ "Aquatic-Terrestrial",
      # Ambystoma spp.: most adult life fossorial/underground, short forest-floor activity periods
      host_genus == "Ambystoma"                                                     ~ "Fossorial-Terrestrial",
      TRUE                                                                          ~ animal_ecomode
    ))
}

metadata_16S <- fix_ecomode(metadata_16S)
metadata_ITS <- fix_ecomode(metadata_ITS)

cat("--- Correction 5: animal_ecomode ---\n")
cat("  16S:\n"); print(table(metadata_16S$animal_ecomode))
cat("  ITS:\n"); print(table(metadata_ITS$animal_ecomode))
cat("\n")

##------------------------------------------------------------------------------
## Correction 6: Standardize ecoregion_III naming ####
##------------------------------------------------------------------------------
# EPA Level III Ecoregion #66 is "Blue Ridge" on the continental and NC maps
# but "Blue Ridge Mountains" on the TN map. Standardizing to "Blue Ridge"
# per the continental map (https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/
# ORD/Ecoregions/us/Eco_Level_III_US.pdf) for consistency.

metadata_16S <- metadata_16S %>%
  mutate(ecoregion_III = as.character(ecoregion_III),
         ecoregion_III = if_else(ecoregion_III == "Blue Ridge Mountains",
                                 "Blue Ridge", ecoregion_III))

metadata_ITS <- metadata_ITS %>%
  mutate(ecoregion_III = as.character(ecoregion_III),
         ecoregion_III = if_else(ecoregion_III == "Blue Ridge Mountains",
                                 "Blue Ridge", ecoregion_III))

cat("--- Correction 6: ecoregion_III ---\n")
cat("  16S:\n"); print(table(metadata_16S$ecoregion_III))
cat("  ITS:\n"); print(table(metadata_ITS$ecoregion_III))
cat("\n")

##------------------------------------------------------------------------------
## Correction 7: Update host species and genus names ####
##------------------------------------------------------------------------------
# Hyla → Dryophytes: North American Hyla reclassified to Dryophytes by Duellman
# et al. (2016); followed by Amphibian Species of the World and TimeTree.
# Geochelone gigantea → Aldabrachelys gigantea: currently accepted name per
# Turtle Taxonomy Working Group; consistent with phylogenetic resources used here.

recode_taxa <- function(metadata) {
  metadata %>%
    mutate(
      host_taxon = recode(host_taxon,
                          "Hyla cinerea"       = "Dryophytes cinereus",
                          "Hyla avivoca"       = "Dryophytes avivoca",
                          "Hyla chrysoscelis"  = "Dryophytes chrysoscelis",
                          "Geochelone gigantea" = "Aldabrachelys gigantea"),
      host_genus = recode(host_genus,
                          "Hyla"       = "Dryophytes",
                          "Geochelone" = "Aldabrachelys")
    )
}

metadata_16S <- recode_taxa(metadata_16S)
metadata_ITS <- recode_taxa(metadata_ITS)

cat("--- Correction 7: host taxon/genus name updates ---\n")
cat("  Dryophytes spp. in 16S:",
    sum(metadata_16S$host_genus == "Dryophytes", na.rm = TRUE), "samples\n")
cat("  Dryophytes spp. in ITS:",
    sum(metadata_ITS$host_genus == "Dryophytes", na.rm = TRUE), "samples\n")
cat("  Aldabrachelys gigantea in 16S:",
    sum(metadata_16S$host_taxon == "Aldabrachelys gigantea", na.rm = TRUE), "samples\n")
cat("  Aldabrachelys gigantea in ITS:",
    sum(metadata_ITS$host_taxon == "Aldabrachelys gigantea", na.rm = TRUE), "samples\n\n")

##------------------------------------------------------------------------------
## Correction 8: Fix Zoo Knoxville / Knoxville Zoo site name ####
##------------------------------------------------------------------------------
metadata_16S <- metadata_16S %>%
  mutate(site = if_else(site == "Knoxville Zoo", "Zoo Knoxville", as.character(site)))

metadata_ITS <- metadata_ITS %>%
  mutate(site = if_else(site == "Knoxville Zoo", "Zoo Knoxville", as.character(site)))

cat("--- Correction 8: Zoo Knoxville site name ---\n")
cat("  'Knoxville Zoo' entries remaining in 16S:",
    sum(metadata_16S$site == "Knoxville Zoo", na.rm = TRUE), "\n")
cat("  'Knoxville Zoo' entries remaining in ITS:",
    sum(metadata_ITS$site == "Knoxville Zoo", na.rm = TRUE), "\n\n")

##------------------------------------------------------------------------------
## Correction 9: Fix missing GPS coordinates for two ITS samples ####
##------------------------------------------------------------------------------
# UHM481-20499 and UHM483-20699 were missing coordinates in the ITS metadata.
#
# NOTE ON COORDINATE PRECISION (public repository version):
# Coordinates below are truncated to whole-degree precision (~100 km) rather
# than the exact decimal values used internally. Several host species in this
# dataset are subject to poaching/collection pressure or are threatened/endangered,
# so precise site coordinates are withheld from the public code/data release.
# Contact the corresponding author for exact coordinates for legitimate
# research requests.

metadata_ITS <- metadata_ITS %>%
  mutate(
    gps_n = if_else(sample_name %in% c("UHM481-20499", "UHM483-20699"), 35.0,  gps_n),
    gps_w = if_else(sample_name %in% c("UHM481-20499", "UHM483-20699"), -86.0, gps_w)
  )

cat("--- Correction 9: GPS coordinates for UHM481-20499 and UHM483-20699 ---\n")
cat("  Coordinates after fix (precision truncated for public release):\n")
print(metadata_ITS %>%
        filter(sample_name %in% c("UHM481-20499", "UHM483-20699")) %>%
        select(sample_name, gps_n, gps_w))
cat("\n")

##------------------------------------------------------------------------------
## Paired sample count ####
##------------------------------------------------------------------------------
animals_16S    <- unique(metadata_16S$animal_number)
animals_ITS    <- unique(metadata_ITS$animal_number)
paired_animals <- intersect(animals_16S, animals_ITS)

cat("--- Paired 16S / ITS samples (shared animal_number) ---\n")
cat("  Unique animals in 16S:", length(animals_16S), "\n")
cat("  Unique animals in ITS:", length(animals_ITS), "\n")
cat("  Paired animals (both datasets):", length(paired_animals), "\n\n")

##------------------------------------------------------------------------------
## Final dataset summary ####
##------------------------------------------------------------------------------
cat("========== FINAL DATASET SUMMARY ==========\n")
cat("16S:\n")
cat("  Samples:", nrow(metadata_16S), "\n")
cat("  Host species:", n_distinct(metadata_16S$host_taxon), "\n")
cat("  env_medium distribution:\n")
print(table(metadata_16S$env_medium))
cat("  env_broad_scale distribution:\n")
print(table(metadata_16S$env_broad_scale))

cat("\nITS:\n")
cat("  Samples:", nrow(metadata_ITS), "\n")
cat("  Host species:", n_distinct(metadata_ITS$host_taxon), "\n")
cat("  env_medium distribution:\n")
print(table(metadata_ITS$env_medium))
cat("  env_broad_scale distribution:\n")
print(table(metadata_ITS$env_broad_scale))
cat("============================================\n\n")

##------------------------------------------------------------------------------
## Export final phyloseq objects ####
##------------------------------------------------------------------------------
# Apply cleaned metadata back to phyloseq objects and save

rownames(metadata_16S) <- metadata_16S$sample_name
sample_data(ps_16S)    <- sample_data(metadata_16S)
saveRDS(ps_16S, file = file.path(dir_processed, "16S_absolute_all_final.rds"))
cat("16S final phyloseq saved:", file.path(dir_processed, "16S_absolute_all_final.rds"), "\n")

rownames(metadata_ITS) <- metadata_ITS$sample_name
sample_data(ps_ITS)    <- sample_data(metadata_ITS)
saveRDS(ps_ITS, file = file.path(dir_processed, "ITS_absolute_all_final.rds"))
cat("ITS final phyloseq saved:", file.path(dir_processed, "ITS_absolute_all_final.rds"), "\n")

##------------------------------------------------------------------------------
## Apply all corrections to pre-depth-filter objects & export ####
##------------------------------------------------------------------------------
# These objects share the same metadata structure as ps_16S / ps_ITS and
# need identical cleaning. All correction functions defined above are reused
# directly — no logic changes, just applied to the _all metadata frames.

metadata_16S_all <- metadata_16S_all %>%
  mutate(env_medium = recode_env_medium(env_medium))

metadata_ITS_all <- metadata_ITS_all %>%
  mutate(env_medium = recode_env_medium(env_medium))

# Correction 2: cloacal swab fix
metadata_16S_all <- metadata_16S_all %>%
  mutate(env_medium = if_else(spiked_volume == 1, "Cloacal swab", env_medium))
metadata_ITS_all <- metadata_ITS_all %>%
  mutate(env_medium = if_else(spiked_volume == 1, "Cloacal swab", env_medium))

# Correction 3: env_broad_scale
metadata_16S_all <- metadata_16S_all %>%
  mutate(env_broad_scale = if_else(tolower(env_broad_scale) == "wild", "Wild", "Captive"))
metadata_ITS_all <- metadata_ITS_all %>%
  mutate(env_broad_scale = if_else(tolower(env_broad_scale) == "wild", "Wild", "Captive"))

# Correction 4: env_medium catch-all
metadata_16S_all <- metadata_16S_all %>% mutate(env_medium = recode_env_medium(env_medium))
metadata_ITS_all <- metadata_ITS_all %>% mutate(env_medium = recode_env_medium(env_medium))

# Correction 5: animal_ecomode
metadata_16S_all <- fix_ecomode(metadata_16S_all)
metadata_ITS_all <- fix_ecomode(metadata_ITS_all)

# Corrections 6–8: ecoregion, taxon names, zoo name
metadata_16S_all <- metadata_16S_all %>%
  mutate(ecoregion_III = as.character(ecoregion_III),
         ecoregion_III = if_else(ecoregion_III == "Blue Ridge Mountains",
                                 "Blue Ridge", ecoregion_III)) %>%
  recode_taxa() %>%
  mutate(site = if_else(site == "Knoxville Zoo", "Zoo Knoxville", as.character(site)))

metadata_ITS_all <- metadata_ITS_all %>%
  mutate(ecoregion_III = as.character(ecoregion_III),
         ecoregion_III = if_else(ecoregion_III == "Blue Ridge Mountains",
                                 "Blue Ridge", ecoregion_III)) %>%
  recode_taxa() %>%
  mutate(site = if_else(site == "Knoxville Zoo", "Zoo Knoxville", as.character(site)))

# Correction 9: GPS fix (ITS only)
metadata_ITS_all <- metadata_ITS_all %>%
  mutate(
    gps_n = if_else(sample_name %in% c("UHM481-20499", "UHM483-20699"), 35.0,  gps_n),
    gps_w = if_else(sample_name %in% c("UHM481-20499", "UHM483-20699"), -86.0, gps_w)
  )

# Apply cleaned metadata back to pre-depth-filter phyloseq objects and save
rownames(metadata_16S_all) <- metadata_16S_all$sample_name
sample_data(ps_16S_all)    <- sample_data(metadata_16S_all)
saveRDS(ps_16S_all, file = file.path(dir_processed, "16S_absolute_all_final_ALL_SPECIES.rds"))
cat("16S all-species phyloseq saved:",
    file.path(dir_processed, "16S_absolute_all_final_ALL_SPECIES.rds"), "\n")

rownames(metadata_ITS_all) <- metadata_ITS_all$sample_name
sample_data(ps_ITS_all)    <- sample_data(metadata_ITS_all)
saveRDS(ps_ITS_all, file = file.path(dir_processed, "ITS_absolute_all_final_ALL_SPECIES.rds"))
cat("ITS all-species phyloseq saved:",
    file.path(dir_processed, "ITS_absolute_all_final_ALL_SPECIES.rds"), "\n")


cat("\nAll outputs saved. Proceed to downstream analyses.\n")

###############################################################################
###############################################################################
###############################################################################
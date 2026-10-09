################################################################################
# MRM Analysis: Phylosymbiosis Testing in Amphibian & Reptile Gut Microbiomes
#
# Approach:  Bootstrap 100 iterations, 1 sample per species per iteration
# Analyses:  Amphibians & Reptiles | 16S & ITS | Bray-Curtis & Jaccard
#
# NOTE ON CLADE DEFINITIONS:
#   Reptilia: Squamata only. Testudines and Crocodilia were sampled exclusively
#   via cloacal swab and are therefore excluded from fecal-only analyses.
#   Amphibia: Caudata + Anura.
#
# ANALYSIS ORDER:
#   Section 1  — Setup: load data, harmonize species names
#   Section 2  — Diet distance matrix: TetrapodTraits integration (Reviewer 1 response)
#   Section 3  — Primary analysis setup (wild samples only)
#   Section 4  — Helper functions
#   Section 5  — Predictor label lookup
#   Section 6  — Primary MRM (wild samples only)
#   Section 7  — Wild + Managed MRM
#   Section 8  — Null phylogeny tests (Bray-Curtis and Jaccard)
#   Section 9  — Figures and supplementary table
#
# REVISION NOTE (Reviewer 1): Diet is now represented using the 10-column
# EltonTraits-style binary diet scheme from TetrapodTraits 3.0.0 (Moura et al.
# 2024; taxonomic backbone per Pyron et al. 2025), combined via
# asymmetric-binary Gower dissimilarity, replacing the original coarse
# 4-category Diet classification (Carnivore/Insectivore/Herbivore/Omnivore)
# See Section 2 for full documentation of this change.
#
# OUTPUT STRUCTURE (all under output/tables/mrm/, via here()):
#   mrm/                     <- main results (figures, summary table)
#   mrm/extended_outputs/    <- raw bootstrap CSVs/RDS per analysis
#   mrm/null_tests/          <- null phylogeny test outputs
#   mrm/phylo_range_diagnostic/  <- Section 11 outputs
#   mrm/range_subsampling/       <- Section 12 outputs
#
# Alexander Rurik
################################################################################

# ========== Load packages ===================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(vegan)
  library(ape)
  library(ecodist)
  library(geosphere)
  library(cluster)
  library(stringr)
  library(phytools)
  library(patchwork)
})

set.seed(52325)
n_bootstrap <- 100

# ========== Project directories =============================================
library(here)
dir_raw       <- here("data", "raw")
dir_processed <- here("data", "processed")
dir_trees     <- here("data", "raw", "trees")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")
dir_mrm       <- file.path(dir_tables, "mrm")

# ── Color palette ─────────────────────────────────────────────────────────────
# 16S = darker shades; ITS = lighter/brighter shades
dataset_colors <- c(
  "16S Amphibians" = "#0052B2",
  "16S Reptiles"   = "#009E73",
  "ITS Amphibians" = "#11E4E9",
  "ITS Reptiles"   = "#69FDAA"
)

################################################################################
# SECTION 1: SETUP — Load data and harmonize species names
################################################################################
cat("=== LOADING DATA ===\n")

ps_16S <- readRDS(file.path(dir_processed, "16S_abs_final_fecal.rds"))
ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_fecal.rds"))

cat("16S samples loaded:", nsamples(ps_16S), "\n")
cat("ITS samples loaded:", nsamples(ps_ITS), "\n")

# ── Output directories ─────────────────────────────────────────────────────────
dir_extended     <- file.path(dir_mrm, "extended_outputs")
dir_null         <- file.path(dir_mrm, "null_tests")

for (d in c(dir_mrm, dir_extended, dir_null)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("Output directories created:\n")
cat("  Extended outputs →", dir_extended,     "\n")
cat("  Null tests       →", dir_null,         "\n")

# ── Recode management status: "Captive" → "Managed" ──────────────────────────
recode_management_status <- function(ps, dataset_name) {
  meta      <- data.frame(sample_data(ps))
  n_recoded <- sum(meta$env_broad_scale == "Captive", na.rm = TRUE)
  meta$env_broad_scale <- recode(meta$env_broad_scale, "Captive" = "Managed")
  sample_data(ps) <- meta
  cat(dataset_name, ": recoded", n_recoded, "samples from 'Captive' to 'Managed'\n")
  cat("  Status breakdown —",
      "Wild:",    sum(meta$env_broad_scale == "Wild",    na.rm = TRUE),
      "| Managed:", sum(meta$env_broad_scale == "Managed", na.rm = TRUE), "\n")
  return(ps)
}

cat("\n=== RECODING MANAGEMENT STATUS ===\n")
ps_16S <- recode_management_status(ps_16S, "16S")
ps_ITS <- recode_management_status(ps_ITS, "ITS")

# ── Species name harmonization ────────────────────────────────────────────────
cat("\n=== HARMONIZING SPECIES NAMES ===\n")

replacements_16S <- c(
  "Sceloporus occidentalis bocourtii"  = "Sceloporus occidentalis",
  "Pituophis catenifer pumilus"        = "Pituophis catenifer",
  "Thamnophis elegans terrestris"      = "Thamnophis elegans",
  "Coluber constrictor mormon"         = "Coluber constrictor",
  "Thamnophis atratus atratus"         = "Thamnophis atratus",
  "Agama picticauda"                   = "Agama atra",
  "Terrapene carolina"                 = "Terrapene ornata",
  "Desmognathus adatsihi"              = "Desmognathus ocoee",
  "Nerodia rhombifer"                  = "Nerodia erythrogaster",
  "Thamnophis proximus"                = "Thamnophis sirtalis",
  "Anolis distichus"                   = "Anolis cristatellus",
  "Atelopus balios"                    = "Atelopus longirostris",
  "Dendrobates tinctorius azureus"     = "Dendrobates tinctorius",
  "Ambystoma annulatum"                = "Ambystoma opacum",
  "Heloderma exasperatum"              = "Heloderma horridum",
  "Hyla cinerea"                       = "Dryophytes cinereus",
  "Hyla avivoca"                       = "Dryophytes avivoca",
  "Hyla chrysoscelis"                  = "Dryophytes chrysoscelis",
  "Geochelone gigantea"                = "Aldabrachelys gigantea",
  "Litoria caerulea"                   = "Ranoidea caerulea"
)

replacements_ITS <- c(
  replacements_16S,
  "Crotalus willardi silus"   = "Crotalus willardi",
  "Crotalus lepidus klauberi" = "Crotalus lepidus"
)

apply_replacements <- function(ps, replacements, dataset_name) {
  meta <- data.frame(sample_data(ps))
  meta$host_taxon <- as.character(meta$host_taxon)
  cat("\n--- Applying to", dataset_name, "---\n")
  changes <- 0
  for (old in names(replacements)) {
    idx <- meta$host_taxon == old
    if (sum(idx, na.rm = TRUE) > 0) {
      cat(" ", old, "->", replacements[old], "(", sum(idx), "samples)\n")
      meta$host_taxon[idx] <- replacements[old]
      changes <- changes + sum(idx)
    }
  }
  sample_data(ps) <- meta
  cat("Total samples updated:", changes, "\n")
  return(ps)
}

ps_16S <- apply_replacements(ps_16S, replacements_16S, "16S")
ps_ITS <- apply_replacements(ps_ITS, replacements_ITS, "ITS")

cat("\nSpecies names harmonized\n")
cat("16S unique species:", length(unique(sample_data(ps_16S)$host_taxon)), "\n")
cat("ITS unique species:", length(unique(sample_data(ps_ITS)$host_taxon)), "\n")


################################################################################
# SECTION 2: DIET DISTANCE MATRIX — TetrapodTraits integration
#            (Response to Reviewer 1)
#
# REVIEWER CONCERN: Our original coarse 4-category diet classification
# (Carnivore/Insectivore/Herbivore/Omnivore) was flagged as lacking rigor --
# e.g., Diadophis and Thamnophis were both coded "Carnivore" despite very
# different diets. Reviewer 1 pointed us to TetrapodTraits (Moura et al. 2024); 
# as a finer-resolution alternative.
#
# WHAT THIS SECTION DOES:
#   1. Loads TetrapodTraits 3.0.0 and its 10 binary EltonTraits-style diet
#      component columns (DietInv, DietVend, DietVect, DietVfish, DietVunk,
#      DietScav, DietFruit, DietNect, DietSeed, DietPlant), following the
#      detailed diet classification scheme used by Youngblut et al. (2019
#      Nature Communications; 2021 Nature Microbiology), who ran MRM on Gower
#      distances calculated from this type of detailed, multi-component diet
#      data rather than a single coarse categorical label.
#   2. Cross-references our harmonized host_taxon names against TetrapodTraits
#      (translating a small number of species to TetrapodTraits' older
#      taxonomic backbone ONLY at the point of the diet lookup -- this never
#      touches host_taxon itself, so phylogenetic tree matching elsewhere in
#      this script is unaffected).
#   3. Fills remaining coverage gaps with (a) Tennessee Wildlife Resources
#      Agency (TWRA) species accounts for 10 Plethodontidae salamanders plus
#      Lithobates palustris, none of which have diet data in TetrapodTraits
#      itself, and (b) targeted manual curation for 5 additional species
#      present only in the wild+managed datasets.
#   4. Joins the resulting per-species diet table onto the full sample_data
#      of ps_16S and ps_ITS (i.e., before the wild-only filter below), so both
#      the primary wild-only MRM and the wild+managed MRM draw on the same
#      enriched metadata.
#
# This diet table is later used in place of the coarse "Diet" column when
# building the diet dissimilarity matrix in calc_all_distances() /
# calc_all_distances_managed() (see helper functions and asymmetric-binary
# Gower notes below). The original coarse "Diet" column is left untouched in
# the metadata (retained for reference/descriptive tables) but is no longer
# used in the MRM predictor matrices.
################################################################################
cat("\n\n################################################################\n")
cat("### SECTION 2: TETRAPODTRAITS DIET INTEGRATION (Reviewer 1 response)\n")
cat("################################################################\n")

# ── TetrapodTraits-specific name mapping (diet lookup ONLY) ──────────────────
# TetrapodTraits 3.0.0 uses an older taxonomic backbone for some species than
# our harmonized host_taxon names. These species have real diet data recorded
# under their OLDER name (e.g. "Hyla avivoca") even though our metadata/MRM
# pipeline correctly uses the updated name (e.g. "Dryophytes avivoca"). Rather
# than reverting our canonical host_taxon names (which would break
# phylogenetic tree matching and downstream MRM code), we translate ONLY at
# the moment of joining against TetrapodTraits, using this lookup: harmonized
# host_taxon -> TetrapodTraits' Scientific.Name. Species not listed here are
# looked up under their unmodified host_taxon name.
tt_diet_name_map <- c(
  "Dryophytes avivoca"        = "Hyla avivoca",
  "Dryophytes chrysoscelis"   = "Hyla chrysoscelis",
  "Dryophytes cinereus"       = "Hyla cinerea",
  "Lithobates catesbeianus"   = "Rana catesbeiana",
  "Lithobates clamitans"      = "Rana clamitans",
  "Lithobates sphenocephalus" = "Rana sphenocephala",
  "Lithobates sylvaticus"     = "Rana sylvatica",
  "Ranoidea caerulea"         = "Litoria caerulea",
  "Agalychnis dacnicolor"     = "Pachymedusa dacnicolor",
  "Testudo horsfieldii"       = "Agrionemys horsfieldii"
)

# ── Load TetrapodTraits 3.0.0 ─────────────────────────────────────────────────
cat("\n=== Loading TetrapodTraits 3.0.0 (Moura et al. 2024) ===\n")

# Third-party dataset (Moura et al. 2024), not redistributed in this repo --
# download TetrapodTraits_v3.0.0.csv from the original source (see README)
# and place it in data/raw/ before running this section.
tt_path <- file.path(dir_raw, "TetrapodTraits_v3.0.0.csv")
tt <- read_csv(tt_path, show_col_types = FALSE)

diet_cols <- c("Scientific.Name", "DietInv", "DietVend", "DietVect", "DietVfish",
               "DietVunk", "DietScav", "DietFruit", "DietNect", "DietSeed",
               "DietPlant", "DietBreadth", "SourceDiet")

tt_diet <- tt %>% select(any_of(diet_cols))

# The 10 binary EltonTraits-style diet component columns. These become the
# input matrix for the asymmetric-binary Gower diet distance (Section 4:
# calc_all_distances() / calc_all_distances_managed()), replacing the single
# coarse "Diet" categorical column used in the original v5 script.
diet_binary_cols <- c("DietInv", "DietVend", "DietVect", "DietVfish",
                      "DietVunk", "DietScav", "DietFruit", "DietNect",
                      "DietSeed", "DietPlant")

# ── Build the full herpetofauna species list to match (union across BOTH
# 16S and ITS phyloseq objects, using ALL samples -- not just wild-only --
# since the wild+managed analyses require diet data for additional species
# not present in the wild-only subset). Restricting to wild-only species
# happens later via the existing clade/order filters; this join simply needs
# to cover every species that could appear in either analysis. ───────────────
all_spp_16S <- unique(as.character(data.frame(sample_data(ps_16S))$host_taxon))
all_spp_ITS <- unique(as.character(data.frame(sample_data(ps_ITS))$host_taxon))
herp_spp_all <- sort(unique(c(all_spp_16S, all_spp_ITS)))

cat("Total unique host species across 16S + ITS datasets:", length(herp_spp_all), "\n")

# ── Match against TetrapodTraits ──────────────────────────────────────────────
diet_match <- tibble(host_taxon = herp_spp_all) %>%
  mutate(
    tt_lookup_name = ifelse(
      host_taxon %in% names(tt_diet_name_map),
      tt_diet_name_map[host_taxon],
      host_taxon
    )
  ) %>%
  left_join(tt_diet, by = c("tt_lookup_name" = "Scientific.Name")) %>%
  mutate(DietDataSource = ifelse(!is.na(DietBreadth), "TetrapodTraits", NA_character_))

# ── Supplemental diet data #1: TWRA (species not covered by TetrapodTraits) ──
# Source: Tennessee Wildlife Resources Agency (TWRA)
#   Salamanders:   https://www.tn.gov/twra/wildlife/amphibians/salamanders
#   Frogs & toads: https://www.tn.gov/twra/wildlife/amphibians/frogs-and-toads.html
# These 10 species (9 Plethodontidae salamanders + Lithobates palustris) have
# no diet data in TetrapodTraits 3.0.0. TWRA species accounts describe
# exclusively invertebrate diets (earthworms, arthropods, mites, snails,
# aquatic and terrestrial insect life stages, etc.) for each, with no mention
# of vertebrate prey, fruit, or other plant material. Coded here as "Inv"
# (DietInv = 1), consistent with how congeners already present in
# TetrapodTraits are coded. Applied ONLY to species still missing after the
# TetrapodTraits join, so any existing valid TetrapodTraits data is never
# overwritten.
twra_supplement <- c(
  "Desmognathus aeneus"      = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Desmognathus conanti"     = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Desmognathus marmoratus"  = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Desmognathus ocoee"       = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Desmognathus santeetlah"  = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Desmognathus wrighti"     = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Eurycea cirrigera"        = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Eurycea wilderae"         = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Plethodon metcalfi"       = "https://www.tn.gov/twra/wildlife/amphibians/salamanders",
  "Lithobates palustris"     = "https://www.tn.gov/twra/wildlife/amphibians/frogs-and-toads.html"
)

diet_match <- diet_match %>%
  mutate(
    supplement_applies = is.na(DietBreadth) & host_taxon %in% names(twra_supplement),
    DietInv     = ifelse(supplement_applies, 1, DietInv),
    DietVend    = ifelse(supplement_applies, 0, DietVend),
    DietVect    = ifelse(supplement_applies, 0, DietVect),
    DietVfish   = ifelse(supplement_applies, 0, DietVfish),
    DietVunk    = ifelse(supplement_applies, 0, DietVunk),
    DietScav    = ifelse(supplement_applies, 0, DietScav),
    DietFruit   = ifelse(supplement_applies, 0, DietFruit),
    DietNect    = ifelse(supplement_applies, 0, DietNect),
    DietSeed    = ifelse(supplement_applies, 0, DietSeed),
    DietPlant   = ifelse(supplement_applies, 0, DietPlant),
    DietBreadth = ifelse(supplement_applies, 1, DietBreadth),
    SourceDiet  = ifelse(supplement_applies, paste0("TWRA, ", twra_supplement[host_taxon]), SourceDiet),
    DietDataSource = ifelse(supplement_applies, "TWRA_supplemental", DietDataSource)
  ) %>%
  select(-supplement_applies)

# ── Supplemental diet data #2: targeted manual curation ──────────────────────
# These 5 species (present only in the wild+managed datasets) have no
# TetrapodTraits or TWRA coverage. Coded from named sources below, following
# the same "Inv/Vend/Vfish/.../Plant" binary scheme as TetrapodTraits/TWRA,
# and flagged with DietDataSource = "Manual_curation" for full transparency
# and reproducibility.
#
#   Agalychnis lemur   -- invertebrates only.
#     Source: https://nationalzoo.si.edu (Lemur Leaf Frog fact sheet)
#   Rhaebo guttatus    -- invertebrates and vertebrate endotherms.
#     Source: https://nationalzoo.si.edu (Smooth-sided Toad fact sheet)
#   Anolis ricordii    -- no species-specific diet information available.
#     Coded identically to Anolis equestris (Knight Anole), a congeneric
#     trunk-crown anole with a broadly similar ecomorphology/foraging guild,
#     as a documented proxy.
#   Cuora flavomarginata -- invertebrates, plants, fruit; opportunistic
#     scavenger (omnivorous box turtle diet, per husbandry/species accounts).
#   Cuora picturata    -- diet ecology closely mirrors C. flavomarginata
#     (congeneric Asian box turtles with equivalent omnivorous diets); coded
#     identically.
manual_diet_supplement <- tribble(
  ~host_taxon,               ~DietInv, ~DietVend, ~DietVect, ~DietVfish, ~DietVunk, ~DietScav, ~DietFruit, ~DietNect, ~DietSeed, ~DietPlant, ~source_note,
  "Agalychnis lemur",               1,         0,         0,          0,         0,         0,          0,         0,         0,          0, "https://nationalzoo.si.edu (Lemur Leaf Frog)",
  "Rhaebo guttatus",                1,         1,         0,          0,         0,         0,          0,         0,         0,          0, "https://nationalzoo.si.edu (Smooth-sided Toad)",
  "Cuora flavomarginata",           1,         0,         0,          0,         0,         1,          1,         0,         0,          1, "Species accounts: omnivorous, invertebrates/plants/fruit/scavenging",
  "Cuora picturata",                1,         0,         0,          0,         0,         1,          1,         0,         0,          1, "Congener of C. flavomarginata; equivalent omnivorous diet"
)

# Anolis ricordii: pull A. equestris' binary diet values dynamically from the
# TetrapodTraits table already loaded above (tt_diet), rather than hardcoding.
equestris_diet <- tt_diet %>% filter(Scientific.Name == "Anolis equestris")

if (nrow(equestris_diet) == 1) {
  anolis_ricordii_row <- tibble(
    host_taxon = "Anolis ricordii",
    !!!equestris_diet[, diet_binary_cols],
    source_note = "Proxy: coded identically to congener Anolis equestris (Knight Anole; no species-specific data available)"
  )
  manual_diet_supplement <- bind_rows(manual_diet_supplement, anolis_ricordii_row)
} else {
  warning("Anolis equestris not found in TetrapodTraits -- Anolis ricordii diet proxy NOT applied. ",
          "Check TetrapodTraits taxonomy/spelling before proceeding.")
}

diet_match <- diet_match %>%
  left_join(manual_diet_supplement, by = "host_taxon", suffix = c("", "_manual")) %>%
  mutate(
    supplement_applies = is.na(DietBreadth) & host_taxon %in% manual_diet_supplement$host_taxon,
    DietInv     = ifelse(supplement_applies, DietInv_manual,   DietInv),
    DietVend    = ifelse(supplement_applies, DietVend_manual,  DietVend),
    DietVect    = ifelse(supplement_applies, DietVect_manual,  DietVect),
    DietVfish   = ifelse(supplement_applies, DietVfish_manual, DietVfish),
    DietVunk    = ifelse(supplement_applies, DietVunk_manual,  DietVunk),
    DietScav    = ifelse(supplement_applies, DietScav_manual,  DietScav),
    DietFruit   = ifelse(supplement_applies, DietFruit_manual, DietFruit),
    DietNect    = ifelse(supplement_applies, DietNect_manual,  DietNect),
    DietSeed    = ifelse(supplement_applies, DietSeed_manual,  DietSeed),
    DietPlant   = ifelse(supplement_applies, DietPlant_manual, DietPlant),
    DietBreadth = ifelse(supplement_applies,
                          DietInv_manual + DietVend_manual + DietVect_manual + DietVfish_manual +
                            DietVunk_manual + DietScav_manual + DietFruit_manual + DietNect_manual +
                            DietSeed_manual + DietPlant_manual,
                          DietBreadth),
    SourceDiet     = ifelse(supplement_applies, source_note, SourceDiet),
    DietDataSource = ifelse(supplement_applies, "Manual_curation", DietDataSource)
  ) %>%
  select(-supplement_applies, -source_note,
         -ends_with("_manual"))

# ── Coverage summary ──────────────────────────────────────────────────────────
cat("\n=== DIET COVERAGE SUMMARY (all species, 16S + ITS combined) ===\n")
n_total   <- nrow(diet_match)
n_missing <- sum(is.na(diet_match$DietBreadth))
cat("Total unique host species:", n_total, "\n")
cat("Matched (TetrapodTraits):",  sum(diet_match$DietDataSource == "TetrapodTraits",   na.rm = TRUE), "\n")
cat("Matched (TWRA supplemental):", sum(diet_match$DietDataSource == "TWRA_supplemental", na.rm = TRUE), "\n")
cat("Matched (Manual curation):",  sum(diet_match$DietDataSource == "Manual_curation",  na.rm = TRUE), "\n")
cat("Still missing diet data:", n_missing, "\n")
if (n_missing > 0) {
  cat("  Missing species:", paste(diet_match$host_taxon[is.na(diet_match$DietBreadth)], collapse = ", "), "\n")
  cat("  NOTE: these species will drop out of any MRM analysis that includes them,\n")
  cat("  via the same NA-handling already in place for the original 'Diet' column.\n")
}

# ── Build readable diet_code from binary columns (for reference/QC only;
# NOT used directly in the distance calculation, which operates on the raw
# binary columns -- see Section 4) ────────────────────────────────────────────
diet_match <- diet_match %>%
  rowwise() %>%
  mutate(
    diet_code = if (is.na(DietBreadth)) {
      NA_character_
    } else {
      paste(
        gsub("^Diet", "", diet_binary_cols[c_across(all_of(diet_binary_cols)) == 1]),
        collapse = "+"
      )
    }
  ) %>%
  ungroup()

# ── Save the full per-species diet match table (mirrors "herp_diet_match.csv"
# from the coverage-check script; the other two exploratory outputs from that
# script -- metadata-vs-TetrapodTraits comparison, diet-code group sizes --
# are not needed here) ────────────────────────────────────────────────────────
f_diet_match <- file.path(dir_mrm, "herp_diet_match_MRM.csv")
write_csv(diet_match, f_diet_match)
cat("\nSaved diet match table:", f_diet_match, "\n")

# ── VALIDATION: species-level diet clustering (cf. Youngblut et al. 2019) ────
# Independent sanity check of the diet distance matrix, separate from the MRM
# pipeline itself. Mirrors the diet-clustering approach in Youngblut et al.
# (2019, Nat Commun), who visually inspected their EltonTraits-derived binary
# diet matrix via hierarchical clustering prior to running MRM. One row per
# SPECIES (not per sample), using the same Jaccard distance used in the MRM
# diet term (see calc_all_distances()/calc_all_distances_managed()).
cat("\n=== Validation: diet dendrogram (cf. Youngblut et al. 2019) ===\n")

diet_species_val <- diet_match %>%
  filter(!is.na(DietBreadth)) %>%
  distinct(host_taxon, .keep_all = TRUE)

diet_mat_val <- diet_species_val %>% select(all_of(diet_binary_cols)) %>% as.data.frame()
rownames(diet_mat_val) <- diet_species_val$host_taxon

diet_dist_val <- vegan::vegdist(diet_mat_val, method = "jaccard", binary = TRUE)
hc_val <- hclust(diet_dist_val)

f_diet_dendro <- file.path(dir_figures, "diet_dendrogram_validation.pdf")
cairo_pdf(f_diet_dendro, width = 10, height = 14)
plot(hc_val, cex = 0.5,
     main = paste0("Diet dendrogram validation (TetrapodTraits, Jaccard; n = ",
                   nrow(diet_mat_val), " species)"))
dev.off()
cat("Saved:", f_diet_dendro, "\n")

# ── Join diet data onto phyloseq sample_data (16S and ITS) ───────────────────
# Joins the full diet_match table (host_taxon, tt_lookup_name, 10 binary diet
# columns, DietBreadth, SourceDiet, DietDataSource, diet_code) onto each
# phyloseq object's sample_data, keyed on host_taxon. The original coarse
# "Diet" column is left in place, untouched, for reference; it is simply not
# used in calc_all_distances()/calc_all_distances_managed() going forward.
join_diet_to_ps <- function(ps, diet_table, dataset_name) {
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  meta$sample_name_join <- rownames(meta)  # preserve sample order/identity through the join
  meta_joined <- meta %>%
    left_join(diet_table, by = "host_taxon") %>%
    column_to_rownames("sample_name_join")
  sample_data(ps) <- sample_data(meta_joined)
  n_missing <- sum(is.na(meta_joined$DietBreadth))
  cat(dataset_name, ": diet data joined --", nrow(meta_joined) - n_missing,
      "/", nrow(meta_joined), "samples with complete diet data\n")
  return(ps)
}

cat("\n=== Joining TetrapodTraits diet data onto phyloseq metadata ===\n")
ps_16S <- join_diet_to_ps(ps_16S, diet_match, "16S")
ps_ITS <- join_diet_to_ps(ps_ITS, diet_match, "ITS")


################################################################################
# SECTION 3: PRIMARY ANALYSIS SETUP — Wild samples only
################################################################################
cat("\n=== Filtering to wild samples only ===\n")

filter_to_wild <- function(ps, dataset_name) {
  meta     <- data.frame(sample_data(ps))
  cat("\n---", dataset_name, "---\n")
  cat("Total:", nrow(meta),
      "| Wild:", sum(meta$env_broad_scale == "Wild", na.rm = TRUE),
      "| Managed (removed):", sum(meta$env_broad_scale == "Managed", na.rm = TRUE), "\n")
  wild_ids <- meta$sample_name[meta$env_broad_scale == "Wild"]
  ps_wild  <- prune_samples(wild_ids, ps)
  cat("After filter:", nsamples(ps_wild), "samples\n")
  return(ps_wild)
}

ps_16S_wild <- filter_to_wild(ps_16S, "16S")
ps_ITS_wild <- filter_to_wild(ps_ITS, "ITS")

# ── Fix GPS coordinate errors ─────────────────────────────────────────────────
cat("\n=== Fixing GPS coordinates ===\n")

fix_gps <- function(ps, dataset_name) {
  meta    <- data.frame(sample_data(ps))
  pos_lon <- meta$gps_w > 0 & !is.na(meta$gps_w)
  if (sum(pos_lon) > 0) {
    cat(dataset_name, ": converting", sum(pos_lon),
        "positive longitudes to negative (Western Hemisphere)\n")
    meta$gps_w[pos_lon] <- -1 * meta$gps_w[pos_lon]
  }
  sample_data(ps) <- meta
  cat("Longitude range:", min(meta$gps_w, na.rm = TRUE), "to",
      max(meta$gps_w, na.rm = TRUE), "\n")
  missing <- sum(is.na(meta$gps_n) | is.na(meta$gps_w))
  if (missing > 0) cat("WARNING:", missing, "samples still missing GPS!\n") else
    cat("All wild samples have GPS coordinates\n")
  return(ps)
}

ps_16S_wild <- fix_gps(ps_16S_wild, "16S")
ps_ITS_wild <- fix_gps(ps_ITS_wild, "ITS")

# ── Load host trees ───────────────────────────────────────────────────────────
cat("\n=== Loading host trees ===\n")

host_tree_16S_full <- read.tree(file.path(dir_trees, "host_tree_16S.nwk"))
host_tree_ITS_full <- read.tree(file.path(dir_trees, "host_tree_ITS.nwk"))

for (tr in c("host_tree_16S_full", "host_tree_ITS_full")) {
  t <- get(tr)
  t$tip.label <- gsub("_", " ", t$tip.label)
  assign(tr, t)
}

cat("16S tree tips:", length(host_tree_16S_full$tip.label), "\n")
cat("ITS tree tips:", length(host_tree_ITS_full$tip.label), "\n")

# ── Define clades and analysis combinations ───────────────────────────────────
# Reptilia: Squamata only. Testudines and Crocodilia were sampled exclusively
# via cloacal swab and are excluded from fecal-only primary analyses.
amphibian_orders <- c("Caudata", "Anura")
reptile_orders   <- c("Squamata")

clades <- list(
  list(name = "Amphibian", orders = amphibian_orders),
  list(name = "Reptile",   orders = reptile_orders)
)

dissim_methods <- list(
  list(name = "BrayCurtis", method = "bray"),
  list(name = "Jaccard",    method = "jaccard")
)

# ── Validate species coverage ─────────────────────────────────────────────────
cat("\n=== Validating species coverage ===\n")

meta_16S_wild <- data.frame(sample_data(ps_16S_wild))
meta_ITS_wild <- data.frame(sample_data(ps_ITS_wild))

check_coverage <- function(meta, tree, label) {
  spp  <- unique(meta$host_taxon)
  miss <- setdiff(spp, tree$tip.label)
  cat(label, "— species in data:", length(spp),
      "| in tree:", length(tree$tip.label),
      "| missing from tree:", length(miss), "\n")
  if (length(miss) > 0) cat("  Missing:", paste(miss, collapse = ", "), "\n")
  return(miss)
}

miss_16S <- check_coverage(meta_16S_wild, host_tree_16S_full, "16S wild")
miss_ITS <- check_coverage(meta_ITS_wild, host_tree_ITS_full, "ITS wild")

if (length(miss_16S) > 0 | length(miss_ITS) > 0) {
  stop("Species in data not present in tree. Update harmonization before proceeding.")
}
cat("All wild species present in trees. Proceeding...\n\n")


################################################################################
# SECTION 4: HELPER FUNCTIONS
################################################################################

# ── Rescale distance matrix to [0, 1] ─────────────────────────────────────────
rescale_dist <- function(d) {
  m <- as.matrix(d)
  as.dist((m - min(m)) / (max(m) - min(m)))
}

# ── Bootstrap summary ─────────────────────────────────────────────────────────
# NOTE (v4): prop_sig_0.05 / prop_sig_0.01 and mean_pval are now computed from
# pval_BH (BH-FDR adjusted within each bootstrap iteration across the 4
# simultaneously-fit predictors; see bootstrap_mrm()/bootstrap_mrm_managed()),
# following Youngblut et al. 2019/2021. Raw (unadjusted) values are retained
# as mean_pval_raw for reference. Column names used downstream (prop_sig_0.05,
# mean_pval) are unchanged, so no other function needs to be modified.
summarize_bootstrap <- function(boot_df, has_management = FALSE) {
  base_summary <- boot_df %>%
    group_by(predictor) %>%
    summarise(
      n_iterations    = n(),
      mean_coef       = mean(dists.comm,    na.rm = TRUE),
      sd_coef         = sd(dists.comm,      na.rm = TRUE),
      median_coef     = median(dists.comm,  na.rm = TRUE),
      q025_coef       = quantile(dists.comm, 0.025, na.rm = TRUE),
      q975_coef       = quantile(dists.comm, 0.975, na.rm = TRUE),
      mean_pval_raw   = mean(pval,          na.rm = TRUE),
      mean_pval       = mean(pval_BH,       na.rm = TRUE),
      prop_sig_0.05   = mean(pval_BH < 0.05, na.rm = TRUE),
      prop_sig_0.01   = mean(pval_BH < 0.01, na.rm = TRUE),
      mean_partial_r2 = mean(partial_r2,    na.rm = TRUE),
      mean_r2_full    = mean(r2_full,       na.rm = TRUE),
      mean_n_species  = mean(n_species,     na.rm = TRUE),
      sd_n_species    = sd(n_species,       na.rm = TRUE),
      .groups = "drop"
    )
  if (has_management && all(c("n_managed", "n_wild") %in% names(boot_df))) {
    mgmt_cols <- boot_df %>%
      group_by(predictor) %>%
      summarise(
        mean_n_managed = mean(n_managed, na.rm = TRUE),
        mean_n_wild    = mean(n_wild,    na.rm = TRUE),
        .groups = "drop"
      )
    base_summary <- left_join(base_summary, mgmt_cols, by = "predictor")
  }
  return(base_summary)
}

# ── Distance matrices: PRIMARY (wild-only) ────────────────────────────────────
# Geographic distances computed via geosphere::distm() (vectorized); distances
# are in km. All predictor matrices are rescaled to [0, 1] before MRM to place
# coefficients on a common scale.
calc_all_distances <- function(ps, sample_ids, tree, dissim_method = "bray") {
  
  ps_sub   <- prune_samples(sample_ids, ps)
  meta_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
  
  # Community dissimilarity
  # binary = TRUE only applies when dissim_method == "jaccard" (ignored by vegdist
  # for "bray"); ensures Jaccard is true presence/absence, not the abundance-
  # weighted quantitative version (2*Bray/(1+Bray); see vegan::vegdist docs).
  comm_dist <- phyloseq::distance(ps_sub, method = dissim_method,
                                  binary = (dissim_method == "jaccard"))
  
  # Host phylogenetic distance (cophenetic distances from TimeTree topology)
  spp_list   <- meta_sub$host_taxon
  tree_sub   <- keep.tip(tree, unique(spp_list))
  phylo_mat  <- cophenetic.phylo(tree_sub)[spp_list, spp_list]
  dimnames(phylo_mat) <- list(sample_ids, sample_ids)
  phylo_dist <- rescale_dist(as.dist(phylo_mat))
  
  # Geographic distance (great-circle km; vectorized via distm)
  geo_df <- meta_sub %>%
    mutate(lon = as.numeric(as.character(gps_w)),
           lat = as.numeric(as.character(gps_n))) %>%
    select(lon, lat)
  rownames(geo_df) <- sample_ids
  
  if (any(is.na(geo_df))) stop("Missing GPS in subsample!")
  
  coord_mat <- as.matrix(geo_df[, c("lon", "lat")])
  g_mat     <- geosphere::distm(coord_mat, fun = geosphere::distGeo) / 1000
  dimnames(g_mat) <- list(sample_ids, sample_ids)
  geo_dist <- rescale_dist(as.dist(g_mat))
  
  # Diet dissimilarity (Jaccard; 10 binary EltonTraits-style diet components
  # from TetrapodTraits -- see Section 2). Jaccard distance on pure
  # presence/absence data is mathematically equivalent to Gower distance with
  # variables declared asymmetric binary (both exclude double-zero/shared-
  # absence "agreement" from the calculation), so this preserves the intended
  # treatment while avoiding a daisy()+type=asymm incompatibility with
  # ecodist::MRM() confirmed during testing (daisy() produced a well-formed,
  # correctly-sized dist object every time, yet MRM() still failed on every
  # bootstrap iteration -- root cause isolated to that interaction, not to
  # the diet data itself). This also directly mirrors Youngblut et al. (2019
  # Nat Comms) detailed-diet MRM code, which uses
  # vegan::vegdist(..., method = 'jaccard') on binary diet data.
  diet_df   <- meta_sub %>% select(all_of(diet_binary_cols))
  rownames(diet_df) <- sample_ids
  diet_dist <- vegan::vegdist(diet_df, method = "jaccard", binary = TRUE)
  
  # Habitat (ecomode) dissimilarity (Gower; categorical)
  ecomode_df   <- meta_sub %>% select(animal_ecomode) %>%
    mutate(animal_ecomode = as.factor(as.character(animal_ecomode)))
  rownames(ecomode_df) <- sample_ids
  ecomode_dist <- daisy(ecomode_df, metric = "gower")
  
  list(comm    = comm_dist,
       phylo   = phylo_dist,
       geo     = geo_dist,
       diet    = diet_dist,
       ecomode = ecomode_dist)
}

# ── Distance matrices: MANAGED + WILD ─────────────────────────────────────────
# Geography is replaced by a binary management-status distance (0 = same status,
# 1 = different status) because managed samples lack meaningful GPS coordinates.
calc_all_distances_managed <- function(ps, sample_ids, tree, dissim_method = "bray") {
  
  ps_sub   <- prune_samples(sample_ids, ps)
  meta_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
  
  # Community dissimilarity
  # binary = TRUE only applies when dissim_method == "jaccard" (ignored by vegdist
  # for "bray"); ensures Jaccard is true presence/absence, not the abundance-
  # weighted quantitative version (2*Bray/(1+Bray); see vegan::vegdist docs).
  comm_dist <- phyloseq::distance(ps_sub, method = dissim_method,
                                  binary = (dissim_method == "jaccard"))
  
  # Host phylogenetic distance
  spp_list  <- meta_sub$host_taxon
  tree_sub  <- keep.tip(tree, unique(spp_list))
  phylo_mat <- cophenetic.phylo(tree_sub)[spp_list, spp_list]
  dimnames(phylo_mat) <- list(sample_ids, sample_ids)
  phylo_dist <- rescale_dist(as.dist(phylo_mat))
  
  # Management status distance (binary: 0 = same, 1 = different)
  status    <- meta_sub$env_broad_scale
  mgmt_mat  <- outer(status, status, FUN = function(x, y) as.numeric(x != y))
  dimnames(mgmt_mat) <- list(sample_ids, sample_ids)
  mgmt_dist <- as.dist(mgmt_mat)
  
  # Diet dissimilarity (Gower; 10 binary EltonTraits-style diet components
  # from TetrapodTraits, treated as asymmetric binary -- see Section 2 and the
  # matching comment in calc_all_distances() above for full rationale.)
  diet_df   <- meta_sub %>% select(all_of(diet_binary_cols))
  rownames(diet_df) <- sample_ids
  diet_dist <- vegan::vegdist(diet_df, method = "jaccard", binary = TRUE)
  
  # Habitat (ecomode) dissimilarity
  ecomode_df <- meta_sub %>% select(animal_ecomode) %>%
    mutate(animal_ecomode = as.factor(as.character(animal_ecomode)))
  rownames(ecomode_df) <- sample_ids
  ecomode_dist <- daisy(ecomode_df, metric = "gower")
  
  list(comm       = comm_dist,
       phylo      = phylo_dist,
       management = mgmt_dist,
       diet       = diet_dist,
       ecomode    = ecomode_dist)
}

# ── Bootstrap function: PRIMARY (wild-only) ───────────────────────────────────
# Model: community ~ phylogeny + geography + diet + habitat
# Partial R² = R²_full - R²_model-without-predictor (computed with nperm = 0
# to avoid redundant permutations). Values can be negative due to collinearity.
bootstrap_mrm <- function(meta, ps, tree, dissim_method = "bray") {
  
  sampled    <- meta %>% group_by(host_taxon) %>% slice_sample(n = 1) %>% ungroup()
  sample_ids <- sampled$sample_name
  
  gps_check <- sampled[, c("sample_name", "gps_n", "gps_w")]
  if (any(is.na(gps_check$gps_n) | is.na(gps_check$gps_w)))
    stop("NA GPS values in bootstrap sample!")
  if (any(abs(gps_check$gps_n) > 90))
    stop("Invalid latitude in bootstrap sample")
  
  dists <- calc_all_distances(ps, sample_ids, tree, dissim_method)
  
  mrm_result <- MRM(dists$comm ~ dists$phylo + dists$geo + dists$diet + dists$ecomode,
                    nperm = 999)
  
  coef_df           <- data.frame(mrm_result$coef)
  coef_df$predictor <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # BH-FDR correction across the 4 simultaneously-fit predictors within this
  # bootstrap iteration (excludes intercept), following Youngblut et al. 2019/2021
  coef_df$pval_BH <- NA_real_
  non_int <- coef_df$predictor != "Int"
  coef_df$pval_BH[non_int] <- p.adjust(coef_df$pval[non_int], method = "BH")
  
  r2_full <- mrm_result$r.squared[1]
  coef_df$partial_r2 <- c(
    NA,
    r2_full - MRM(dists$comm ~ dists$geo   + dists$diet + dists$ecomode, nperm = 0)$r.squared[1],
    r2_full - MRM(dists$comm ~ dists$phylo + dists$diet + dists$ecomode, nperm = 0)$r.squared[1],
    r2_full - MRM(dists$comm ~ dists$phylo + dists$geo  + dists$ecomode, nperm = 0)$r.squared[1],
    r2_full - MRM(dists$comm ~ dists$phylo + dists$geo  + dists$diet,    nperm = 0)$r.squared[1]
  )
  coef_df$r2_full   <- r2_full
  coef_df$n_species <- length(unique(sampled$host_taxon))
  
  return(coef_df)
}

# ── Bootstrap function: MANAGED + WILD ────────────────────────────────────────
# Model: community ~ phylogeny + management status + diet + habitat
bootstrap_mrm_managed <- function(meta, ps, tree, dissim_method = "bray") {
  
  sampled    <- meta %>% group_by(host_taxon) %>% slice_sample(n = 1) %>% ungroup()
  sample_ids <- sampled$sample_name
  
  dists <- calc_all_distances_managed(ps, sample_ids, tree, dissim_method)
  
  mrm_result <- MRM(dists$comm ~ dists$phylo + dists$management + dists$diet + dists$ecomode,
                    nperm = 999)
  
  coef_df           <- data.frame(mrm_result$coef)
  coef_df$predictor <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # BH-FDR correction across the 4 simultaneously-fit predictors within this
  # bootstrap iteration (excludes intercept), following Youngblut et al. 2019/2021
  coef_df$pval_BH <- NA_real_
  non_int <- coef_df$predictor != "Int"
  coef_df$pval_BH[non_int] <- p.adjust(coef_df$pval[non_int], method = "BH")
  
  r2_full <- mrm_result$r.squared[1]
  coef_df$partial_r2 <- c(
    NA,
    r2_full - MRM(dists$comm ~ dists$management + dists$diet + dists$ecomode, nperm = 0)$r.squared[1],
    r2_full - MRM(dists$comm ~ dists$phylo + dists$diet + dists$ecomode,       nperm = 0)$r.squared[1],
    r2_full - MRM(dists$comm ~ dists$phylo + dists$management + dists$ecomode, nperm = 0)$r.squared[1],
    r2_full - MRM(dists$comm ~ dists$phylo + dists$management + dists$diet,    nperm = 0)$r.squared[1]
  )
  coef_df$r2_full   <- r2_full
  coef_df$n_species <- length(unique(sampled$host_taxon))
  coef_df$n_managed <- sum(sampled$env_broad_scale == "Managed")
  coef_df$n_wild    <- sum(sampled$env_broad_scale == "Wild")
  
  return(coef_df)
}

# ── Generic bootstrap runner ──────────────────────────────────────────────────
# Seeds each iteration as set.seed(52325 + i) for full per-iteration
# reproducibility, regardless of execution order across Parts.
run_bootstrap <- function(boot_fn, n = n_bootstrap, label = "", ...) {
  results <- vector("list", n)
  for (i in seq_len(n)) {
    set.seed(52325 + i)
    if (i %% 10 == 0) cat("  Iteration", i, "/", n, "\n")
    tryCatch({
      results[[i]]           <- boot_fn(...)
      results[[i]]$iteration <- i
    }, error = function(e) {
      cat("  Error iteration", i, "(", label, "):", e$message, "\n")
    })
  }
  bind_rows(results)
}


################################################################################
# SECTION 5: PREDICTOR LABEL LOOKUP
################################################################################
rename_predictors <- function(df, analysis_type = "primary") {
  if (analysis_type == "primary") {
    df %>% mutate(predictor = recode(predictor,
                                     "dists$phylo"   = "Phylogeny",
                                     "dists$geo"     = "Geography",
                                     "dists$diet"    = "Diet",
                                     "dists$ecomode" = "Habitat"
    ))
  } else if (analysis_type == "managed_wild") {
    df %>% mutate(predictor = recode(predictor,
                                     "dists$phylo"       = "Phylogeny",
                                     "dists$management"  = "Management status",
                                     "dists$diet"        = "Diet",
                                     "dists$ecomode"     = "Habitat"
    ))
  }
}


################################################################################
# SECTION 6: PRIMARY MRM ANALYSES — Wild samples only
# Model: community dissimilarity ~ phylogeny + geography + diet + habitat
################################################################################
cat("\n=== PART 2: PRIMARY MRM ANALYSES (wild samples only) ===\n")

datasets_primary <- list(
  list(name = "16S", ps = ps_16S_wild, meta = meta_16S_wild, tree = host_tree_16S_full),
  list(name = "ITS", ps = ps_ITS_wild, meta = meta_ITS_wild, tree = host_tree_ITS_full)
)

results_primary <- list()

for (dataset in datasets_primary) {
  for (clade in clades) {
    for (dissim in dissim_methods) {
      
      key <- paste(dataset$name, clade$name, dissim$name, sep = "_")
      cat("\n--- PRIMARY:", key, "---\n")
      
      meta_clade <- dataset$meta %>% filter(Clade_Order %in% clade$orders)
      cat("Pool:", nrow(meta_clade), "samples |",
          length(unique(meta_clade$host_taxon)), "species\n")
      
      if (nrow(meta_clade) < 10) { cat("Skipping — too few samples\n"); next }
      
      ps_clade <- prune_samples(meta_clade$sample_name, dataset$ps)
      cat("Running", n_bootstrap, "bootstrap iterations...\n")
      
      boot_df <- run_bootstrap(
        boot_fn       = bootstrap_mrm,
        label         = key,
        meta          = meta_clade,
        ps            = ps_clade,
        tree          = dataset$tree,
        dissim_method = dissim$method
      )
      
      write.csv(boot_df,
                file.path(dir_extended, paste0(key, "_bootstrap.csv")),
                row.names = FALSE)
      saveRDS(boot_df,
              file.path(dir_extended, paste0(key, "_bootstrap.rds")))
      
      summary_df <- summarize_bootstrap(boot_df)
      write.csv(summary_df,
                file.path(dir_extended, paste0(key, "_summary.csv")),
                row.names = FALSE)
      
      results_primary[[key]] <- list(bootstrap = boot_df, summary = summary_df)
      cat("Saved:", key, "\n")
    }
  }
}

cat("\nPRIMARY ANALYSES COMPLETE —", length(results_primary), "analyses\n")
saveRDS(results_primary,file.path(dir_extended, "primary_all_results.rds"))


################################################################################
# SECTION 7: WILD + MANAGED MRM
# Model: community ~ phylogeny + management status + diet + habitat
# Geography replaced by management status (binary: Wild vs Managed) because
# managed samples lack meaningful GPS coordinates.
################################################################################
cat("\n=== PART 3: WILD + MANAGED MRM ===\n")

meta_16S_all <- data.frame(sample_data(ps_16S))
meta_ITS_all <- data.frame(sample_data(ps_ITS))

datasets_managed <- list(
  list(name = "16S", ps = ps_16S, meta = meta_16S_all, tree = host_tree_16S_full),
  list(name = "ITS", ps = ps_ITS, meta = meta_ITS_all, tree = host_tree_ITS_full)
)

results_managed <- list()

for (dataset in datasets_managed) {
  for (clade in clades) {
    for (dissim in dissim_methods) {
      
      key <- paste(dataset$name, clade$name, dissim$name, sep = "_")
      cat("\n--- MANAGED+WILD:", key, "---\n")
      
      meta_clade <- dataset$meta %>% filter(Clade_Order %in% clade$orders)
      cat("Pool:", nrow(meta_clade), "samples |",
          length(unique(meta_clade$host_taxon)), "species\n")
      cat("  Wild:", sum(meta_clade$env_broad_scale == "Wild"),
          "| Managed:", sum(meta_clade$env_broad_scale == "Managed"), "\n")
      
      if (nrow(meta_clade) < 10) { cat("Skipping — too few samples\n"); next }
      
      ps_clade <- prune_samples(meta_clade$sample_name, dataset$ps)
      
      boot_df <- run_bootstrap(
        boot_fn       = bootstrap_mrm_managed,
        label         = paste("managed+wild", key),
        meta          = meta_clade,
        ps            = ps_clade,
        tree          = dataset$tree,
        dissim_method = dissim$method
      )
      
      write.csv(boot_df,
                file.path(dir_extended, paste0("managed_wild_", key, "_bootstrap.csv")),
                row.names = FALSE)
      saveRDS(boot_df,
              file.path(dir_extended, paste0("managed_wild_", key, "_bootstrap.rds")))
      
      summary_df <- summarize_bootstrap(boot_df, has_management = TRUE)
      write.csv(summary_df,
                file.path(dir_extended, paste0("managed_wild_", key, "_summary.csv")),
                row.names = FALSE)
      
      results_managed[[key]] <- list(bootstrap = boot_df, summary = summary_df)
      cat("Saved:", key, "\n")
    }
  }
}

cat("\nMANAGED+WILD ANALYSES COMPLETE —", length(results_managed), "analyses\n")
saveRDS(results_managed,file.path(dir_extended, "managed_wild_all_results.rds"))


################################################################################
# SECTION 8: NULL PHYLOGENY TESTS
# Tests whether the observed phylogenetic partial R² exceeds chance expectation.
# 1,000 random Yule phylogenies (pbtree) are generated with tip labels randomly
# assigned from the observed species pool. For each null tree, 20 bootstrap
# iterations are run (vs. 100 for the observed) to reduce computational burden;
# with 1,000 null trees this asymmetry does not meaningfully inflate variance of
# the null distribution. Empirical p = proportion of null partial R² >= observed.
################################################################################
cat("\n=== PART 4: NULL PHYLOGENY TESTS ===\n")

run_mrm_null_test <- function(ps, meta, tree, dissim_method = "bray",
                              n_null = 1000, n_boot_per_null = 20) {
  
  n_spp <- length(unique(meta$host_taxon))
  cat("  Species:", n_spp, "| Null trees:", n_null,
      "| Bootstrap per null:", n_boot_per_null, "\n")
  
  # Observed: run n_boot_per_null iterations with the real tree
  obs_r2_phylo <- obs_r2_full <- numeric(n_boot_per_null)
  for (i in seq_len(n_boot_per_null)) {
    set.seed(52325 + i)
    r                <- bootstrap_mrm(meta, ps, tree, dissim_method)
    phylo_row        <- grep("phylo", r$predictor, ignore.case = TRUE)
    obs_r2_phylo[i]  <- r$partial_r2[phylo_row]
    obs_r2_full[i]   <- r$r2_full[1]
  }
  mean_obs_phylo <- mean(obs_r2_phylo, na.rm = TRUE)
  mean_obs_full  <- mean(obs_r2_full,  na.rm = TRUE)
  cat("  Observed phylo partial R² =", round(mean_obs_phylo, 4), "\n")
  
  # Null: generate random topologies and compute phylo partial R² under each
  null_r2_phylo <- null_r2_full <- numeric(n_null)
  for (i in seq_len(n_null)) {
    set.seed(52325 + i)
    if (i %% 100 == 0) cat("    Null", i, "/", n_null, "\n")
    null_tree           <- pbtree(n = n_spp)
    null_tree$tip.label <- sample(unique(meta$host_taxon), n_spp, replace = FALSE)
    tryCatch({
      r                  <- bootstrap_mrm(meta, ps, null_tree, dissim_method)
      phylo_row          <- grep("phylo", r$predictor, ignore.case = TRUE)
      null_r2_phylo[i]  <- r$partial_r2[phylo_row]
      null_r2_full[i]   <- r$r2_full[1]
    }, error = function(e) {
      null_r2_phylo[i] <<- NA
      null_r2_full[i]  <<- NA
    })
  }
  null_r2_phylo <- null_r2_phylo[!is.na(null_r2_phylo)]
  null_r2_full  <- null_r2_full[!is.na(null_r2_full)]
  
  # Empirical p-value: proportion of null values >= observed mean
  p_phylo <- mean(null_r2_phylo >= mean_obs_phylo)
  p_full  <- mean(null_r2_full  >= mean_obs_full)
  cat("  p (phylo):", round(p_phylo, 4),
      "| p (full):", round(p_full,  4), "\n")
  
  list(obs_r2_phylo      = obs_r2_phylo,
       obs_r2_full       = obs_r2_full,
       mean_obs_r2_phylo = mean_obs_phylo,
       mean_obs_r2_full  = mean_obs_full,
       null_r2_phylo     = null_r2_phylo,
       null_r2_full      = null_r2_full,
       p_phylo           = p_phylo,
       p_full            = p_full,
       n_null            = length(null_r2_phylo),
       n_species         = n_spp,
       dissim_method     = dissim_method)
}

null_combos <- list(
  list(label = "16S_Amphibian", ps = ps_16S_wild,
       meta  = meta_16S_wild %>% filter(Clade_Order %in% amphibian_orders),
       tree  = host_tree_16S_full),
  list(label = "16S_Reptile",   ps = ps_16S_wild,
       meta  = meta_16S_wild %>% filter(Clade_Order %in% reptile_orders),
       tree  = host_tree_16S_full),
  list(label = "ITS_Amphibian", ps = ps_ITS_wild,
       meta  = meta_ITS_wild %>% filter(Clade_Order %in% amphibian_orders),
       tree  = host_tree_ITS_full),
  list(label = "ITS_Reptile",   ps = ps_ITS_wild,
       meta  = meta_ITS_wild %>% filter(Clade_Order %in% reptile_orders),
       tree  = host_tree_ITS_full)
)

null_results_list <- list()
for (combo in null_combos) {
  cat("\nNULL TEST:", combo$label, "\n")
  nr <- run_mrm_null_test(
    ps              = prune_samples(combo$meta$sample_name, combo$ps),
    meta            = combo$meta,
    tree            = combo$tree,
    dissim_method   = "bray",
    n_null          = 1000,
    n_boot_per_null = 20
  )
  null_results_list[[combo$label]] <- nr
  saveRDS(nr, file.path(dir_null,
                        paste0("null_", combo$label, "_BrayCurtis.rds")))
}

# ── Null test figure: 2×2 multipanel ─────────────────────────────────────────
plot_null_panel <- function(null_res, panel_label, analysis_name) {
  h <- hist(null_res$null_r2_phylo, breaks = 30, col = "lightblue",
            main = "", xlab = expression("Phylogenetic partial " * R^2),
            xlim = c(min(c(null_res$null_r2_phylo, 0)),
                     max(c(null_res$null_r2_phylo,
                           null_res$mean_obs_r2_phylo)) * 1.15),
            cex.lab = 1.2, cex.axis = 1.1)
  mtext(panel_label,   side = 3, line = 2,   adj = 0, cex = 1.5, font = 1)
  mtext(analysis_name, side = 3, line = 0.5, cex = 1.1)
  abline(v = null_res$mean_obs_r2_phylo, col = "red", lwd = 3)
  x_range <- par("usr")[2] - par("usr")[1]
  rel_pos <- (null_res$mean_obs_r2_phylo - par("usr")[1]) / x_range
  txt_pos <- if (rel_pos > 0.6) 2L else 4L
  text(x      = null_res$mean_obs_r2_phylo,
       y      = max(h$counts) * 0.9,
       labels = paste0("Observed: ", round(null_res$mean_obs_r2_phylo, 3),
                       "\np = ",     round(null_res$p_phylo, 4)),
       col    = "red", pos = txt_pos, cex = 1.1, font = 2)
}

cairo_pdf(file.path(dir_null, "null_phylogeny_multipanel.pdf"), width = 14, height = 12)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
plot_null_panel(null_results_list[["16S_Amphibian"]], "A)", "16S Amphibians")
plot_null_panel(null_results_list[["16S_Reptile"]],   "B)", "16S Reptiles")
plot_null_panel(null_results_list[["ITS_Amphibian"]], "C)", "ITS Amphibians")
plot_null_panel(null_results_list[["ITS_Reptile"]],   "D)", "ITS Reptiles")
dev.off()
cat("Saved:", file.path(dir_null, "null_phylogeny_multipanel.pdf"), "\n")


################################################################################
# SECTION 8b: NULL PHYLOGENY TESTS — JACCARD
# Same procedure as the Bray-Curtis null test above (1,000 random Yule
# phylogenies via pbtree, tip labels randomly reassigned from the observed
# species pool, 20 bootstrap iterations per null tree vs. 100 for the
# observed), run on true presence/absence Jaccard distances instead of
# Bray-Curtis. Reuses run_mrm_null_test(),ull_combos, and plot_null_panel() 
# defined above; only dissim_method and output filenames differ.
################################################################################
cat("\n=== PART 4b: NULL PHYLOGENY TESTS (JACCARD) ===\n")

null_results_list_jaccard <- list()
for (combo in null_combos) {
  cat("\nNULL TEST (Jaccard):", combo$label, "\n")
  nr <- run_mrm_null_test(
    ps              = prune_samples(combo$meta$sample_name, combo$ps),
    meta            = combo$meta,
    tree            = combo$tree,
    dissim_method   = "jaccard",
    n_null          = 1000,
    n_boot_per_null = 20
  )
  null_results_list_jaccard[[combo$label]] <- nr
  saveRDS(nr, file.path(dir_null,
                        paste0("null_", combo$label, "_Jaccard.rds")))
}

# ── Null test figure: 2×2 multipanel (Jaccard) ───────────────────────────────
cairo_pdf(file.path(dir_null, "null_phylogeny_multipanel_Jaccard.pdf"), width = 14, height = 12)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
plot_null_panel(null_results_list_jaccard[["16S_Amphibian"]], "A)", "16S Amphibians")
plot_null_panel(null_results_list_jaccard[["16S_Reptile"]],   "B)", "16S Reptiles")
plot_null_panel(null_results_list_jaccard[["ITS_Amphibian"]], "C)", "ITS Amphibians")
plot_null_panel(null_results_list_jaccard[["ITS_Reptile"]],   "D)", "ITS Reptiles")
dev.off()
cat("Saved:", file.path(dir_null, "null_phylogeny_multipanel_Jaccard.pdf"), "\n")


################################################################################
# SECTION 9: FIGURES AND SUPPLEMENTARY TABLE
################################################################################
cat("\n=== PART 6: FIGURES AND SUPPLEMENTARY TABLE ===\n")

# ── Helper: prepare summary for plotting ─────────────────────────────────────
prep_plot_data <- function(results_list, analysis_label, analysis_type = "primary") {
  bind_rows(
    lapply(names(results_list), function(key) {
      parts  <- strsplit(key, "_")[[1]]
      marker <- parts[1]
      clade  <- parts[2]
      dissim <- paste(parts[3:length(parts)], collapse = "_")
      
      results_list[[key]]$summary %>%
        filter(predictor != "Int") %>%
        rename_predictors(analysis_type = analysis_type) %>%
        mutate(
          dataset  = paste(marker, paste0(clade, "s")),
          dissim   = dissim,
          analysis = analysis_label,
          # Proportion of iterations significant — primary inference metric
          prop_sig_label = paste0(round(prop_sig_0.05 * 100), "%"),
          # Flag predictors meeting the ≥95% bootstrap-significance threshold,
          # so labels can be rendered bold in the figure exactly as they were
          # previously bolded manually post-hoc.
          sig_face = ifelse(prop_sig_0.05 >= 0.95, "bold", "plain")
        )
    })
  )
}

plot_primary      <- prep_plot_data(results_primary, "Wild only",      "primary")
plot_managed_wild <- prep_plot_data(results_managed, "Wild + Managed", "managed_wild")

split_dissim <- function(df) {
  list(
    bc  = df %>% filter(dissim == "BrayCurtis"),
    jac = df %>% filter(dissim == "Jaccard")
  )
}

pd_primary <- split_dissim(plot_primary)
pd_mw      <- split_dissim(plot_managed_wild)

# ── Shared y-axis limits ──────────────────────────────────────────────────────
# Ceiling computed from data + error bar + label clearance (12% headroom).
# Floor at 0 or below if any mean - SD is negative.
compute_ylims <- function(...) {
  dfs  <- list(...)
  ymin <- min(sapply(dfs, function(d) min(d$mean_coef - d$sd_coef, na.rm = TRUE)))
  ymax <- max(sapply(dfs, function(d) max(d$mean_coef + d$sd_coef, na.rm = TRUE)))
  ymin <- min(ymin, 0)
  c(ymin, ymax * 1.12)
}

ylim_bc  <- compute_ylims(pd_primary$bc,  pd_mw$bc)
ylim_jac <- compute_ylims(pd_primary$jac, pd_mw$jac)

# ── Sample size annotation ────────────────────────────────────────────────────
build_n_label <- function(df) {
  df %>%
    group_by(dataset, dissim, analysis) %>%
    summarise(n_label = paste0("n=", round(unique(mean_n_species)), " spp"),
              .groups = "drop")
}

n_labels_primary <- build_n_label(plot_primary)
n_labels_mw      <- build_n_label(plot_managed_wild)

# ── Figure constructor ────────────────────────────────────────────────────────
# Bars show mean MRM coefficient ± SD across 100 bootstrap iterations.
# Labels above bars show the proportion of iterations with BH-adjusted p < 0.05
# (e.g. "94%"), providing a measure of effect consistency across species
# subsamples. Predictor order: Diet, Habitat, Geography/Management status, Phylogeny.
# Labels are bolded automatically when the predictor meets the ≥95% bootstrap-
# significance threshold (previously done manually post-hoc; see sig_face in
# prep_plot_data()).
make_mrm_plot <- function(df, ylimits, title_str, n_df = NULL, show_legend = TRUE) {
  
  pred_levels <- c("Diet", "Habitat", "Geography", "Management status", "Phylogeny")
  df <- df %>%
    mutate(
      predictor = factor(predictor,
                         levels = intersect(pred_levels, unique(predictor))),
      # Force 'dataset' to an explicit factor with fixed levels (matching
      # dataset_colors) so geom_col() and geom_text() compute IDENTICAL
      # dodge positions. Relying on implicit grouping via fill= alone can
      # cause the two layers to dodge slightly differently when any
      # predictor is missing a dataset, silently shifting labels onto the
      # wrong bar.
      dataset = factor(dataset, levels = names(dataset_colors))
    )
  
  # Prop-sig label y-position: just above the top of the error bar
  sig_offset <- diff(ylimits) * 0.025
  
  p <- ggplot(df, aes(x = predictor, y = mean_coef, fill = dataset, group = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_coef - sd_coef, ymax = mean_coef + sd_coef),
                  position = position_dodge(width = 0.8), width = 0.25,
                  linewidth = 0.5) +
    # Proportion of iterations significant, displayed instead of asterisks
    geom_text(aes(label = prop_sig_label,
                  y     = mean_coef + sd_coef + sig_offset,
                  fontface = sig_face,
                  group = dataset),
              position = position_dodge(width = 0.8),
              vjust    = 0,
              size     = 3.2,
              color    = "grey20") +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = dataset_colors) +
    coord_cartesian(ylim = ylimits, clip = "off") +
    labs(x = NULL, y = "MRM coefficient (mean ± SD)",
         fill = "Dataset", title = title_str) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 40, hjust = 1, size = 10),
      legend.position  = if (show_legend) "right" else "none",
      plot.title       = element_text(size = 11, face = "bold"),
      plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
    )
  
  if (!is.null(n_df)) {
    n_sub <- n_df %>%
      filter(analysis == unique(df$analysis),
             dissim   == unique(df$dissim)) %>%
      arrange(dataset)
    n_str <- paste(paste0(n_sub$dataset, ": ", n_sub$n_label), collapse = "  |  ")
    p <- p + labs(subtitle = n_str) +
      theme(plot.subtitle = element_text(size = 8, color = "grey30"))
  }
  
  return(p)
}

# ── Build individual figures ──────────────────────────────────────────────────
fig_primary_bc <- make_mrm_plot(
  df        = pd_primary$bc,
  ylimits   = ylim_bc,
  title_str = "MRM: Wild Samples (Bray-Curtis)",
  n_df      = NULL
)

fig_primary_jac <- make_mrm_plot(
  df        = pd_primary$jac,
  ylimits   = ylim_jac,
  title_str = "MRM: Wild Samples (Jaccard)",
  n_df      = n_labels_primary
)

fig_mw_bc <- make_mrm_plot(
  df        = pd_mw$bc,
  ylimits   = ylim_bc,
  title_str = "MRM: Wild + Managed Samples (Bray-Curtis)",
  n_df      = NULL
)

fig_mw_jac <- make_mrm_plot(
  df        = pd_mw$jac,
  ylimits   = ylim_jac,
  title_str = "MRM: Wild + Managed Samples (Jaccard)",
  n_df      = n_labels_mw
)

# ── Multipanel figures ────────────────────────────────────────────────────────

# Main text figure: Bray-Curtis, Wild (A) | Wild+Managed (B), shared legend
fig_multipanel_bc <- (
  (fig_primary_bc + theme(legend.position = "none") +
     labs(tag = "A)", title = "Wild only (Bray-Curtis)")) |
    (fig_mw_bc +
       labs(tag = "B)", title = "Wild + Managed (Bray-Curtis)"))
) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Supplementary figure: Jaccard, Wild (A) / Wild+Managed (B), shared legend
fig_multipanel_jac <- (
  (fig_primary_jac + theme(legend.position = "none") +
     labs(tag = "A)", title = "Wild only (Jaccard)")) /
    (fig_mw_jac +
       labs(tag = "B)", title = "Wild + Managed (Jaccard)"))
) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ── Save all figures ──────────────────────────────────────────────────────────
figs <- list(
  list(obj = fig_primary_bc,      file = "MRM_wild_BrayCurtis",        w = 9,  h = 6),
  list(obj = fig_primary_jac,     file = "MRM_wild_Jaccard",           w = 9,  h = 6),
  list(obj = fig_mw_bc,           file = "MRM_wildManaged_BrayCurtis", w = 9,  h = 6),
  list(obj = fig_mw_jac,          file = "MRM_wildManaged_Jaccard",    w = 9,  h = 6),
  list(obj = fig_multipanel_bc,   file = "MRM_multipanel_BrayCurtis",  w = 16, h = 6),
  list(obj = fig_multipanel_jac,  file = "MRM_multipanel_Jaccard",     w = 9,  h = 12)
)

for (fig in figs) {
  fname <- file.path(dir_figures, paste0(fig$file, ".pdf"))
  ggsave(fname, fig$obj,
         width = fig$w, height = fig$h,
         device = cairo_pdf, dpi = 400)
  cat("Saved:", fname, "\n")
}

# ── Supplementary table ───────────────────────────────────────────────────────
cat("\n=== SUPPLEMENTARY TABLE ===\n")

extract_table_rows <- function(results_list, analysis_label, analysis_type) {
  bind_rows(
    lapply(names(results_list), function(key) {
      parts        <- strsplit(key, "_")[[1]]
      marker       <- parts[1]
      clade        <- parts[2]
      dissim_label <- ifelse(grepl("Bray", key), "Bray-Curtis", "Jaccard")
      
      results_list[[key]]$summary %>%
        filter(predictor != "Int") %>%
        rename_predictors(analysis_type = analysis_type) %>%
        mutate(
          analysis        = analysis_label,
          dataset         = paste(marker, paste0(clade, "s")),
          distance_metric = dissim_label,
          n_species_mean_sd = paste0(round(mean_n_species, 1),
                                     " ± ", round(sd_n_species, 1))
        ) %>%
        select(analysis, dataset, distance_metric, predictor,
               mean_coef, sd_coef, mean_partial_r2, mean_pval,
               prop_sig_0.05, n_species_mean_sd, mean_r2_full)
    })
  )
}

pred_order <- c("Diet", "Habitat", "Geography", "Management status", "Phylogeny")

supp_table <- bind_rows(
  extract_table_rows(results_primary, "Wild only",      "primary"),
  extract_table_rows(results_managed, "Wild + Managed", "managed_wild")
) %>%
  mutate(
    analysis        = factor(analysis,
                             levels = c("Wild only", "Wild + Managed")),
    dataset         = factor(dataset,
                             levels = c("16S Amphibians", "16S Reptiles",
                                        "ITS Amphibians", "ITS Reptiles")),
    distance_metric = factor(distance_metric, levels = c("Bray-Curtis", "Jaccard")),
    predictor       = factor(predictor, levels = pred_order)
  ) %>%
  arrange(analysis, dataset, distance_metric, predictor) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  rename(
    Analysis          = analysis,
    Dataset           = dataset,
    Distance_metric   = distance_metric,
    Predictor         = predictor,
    Mean_coefficient  = mean_coef,
    SD_coefficient    = sd_coef,
    Mean_partial_R2   = mean_partial_r2,
    Mean_pval         = mean_pval,
    Prop_sig_0.05     = prop_sig_0.05,
    N_species_mean_SD = n_species_mean_sd,
    Mean_full_R2      = mean_r2_full
  )

supp_fname <- file.path(dir_mrm, "MRM_supplementary_table.csv")
write.csv(supp_table, supp_fname, row.names = FALSE)
cat("Saved supplementary table:", supp_fname, "\n")

cat("\n--- Table preview ---\n")
print(as.data.frame(supp_table), max = 200)


################################################################################
################################################################################
################################################################################


################################################################################
# SECTION 10: VALIDATION — Mantel test, old vs. new diet-phylogeny collinearity
#
# Checks whether the original coarse 4-category diet classification was
# collinear with host phylogeny (i.e., functioning partly as a taxonomy
# proxy), and whether the TetrapodTraits-based diet matrix reduces that
# collinearity. If the coarse diet variable was absorbing phylogenetic
# variance in the MRM models, we'd expect a HIGHER diet-phylogeny Mantel
# correlation under the old scheme than the new one — this would explain why
# diet's MRM coefficient/partial R2 dropped so consistently after switching
# to TetrapodTraits (see Section 2).
################################################################################
cat("\n\n################################################################\n")
cat("### SECTION 10: MANTEL VALIDATION — diet-phylogeny collinearity\n")
cat("################################################################\n")

mantel_clades <- list(
  list(label = "16S_Amphibian", meta = meta_16S_wild %>% filter(Clade_Order %in% amphibian_orders), tree = host_tree_16S_full),
  list(label = "16S_Reptile",   meta = meta_16S_wild %>% filter(Clade_Order %in% reptile_orders),   tree = host_tree_16S_full),
  list(label = "ITS_Amphibian", meta = meta_ITS_wild %>% filter(Clade_Order %in% amphibian_orders), tree = host_tree_ITS_full),
  list(label = "ITS_Reptile",   meta = meta_ITS_wild %>% filter(Clade_Order %in% reptile_orders),   tree = host_tree_ITS_full)
)

run_diet_phylo_mantel <- function(clade_entry) {
  
  # One row per species (diet/phylogeny are species-level, not sample-level)
  meta_sp <- clade_entry$meta %>% distinct(host_taxon, .keep_all = TRUE)
  spp_list <- meta_sp$host_taxon
  
  # Restrict to species with BOTH tree placement and complete new diet data,
  # so old and new Mantel tests are computed on the identical species set
  # (required for a fair comparison).
  spp_list <- intersect(spp_list, clade_entry$tree$tip.label)
  meta_sp <- meta_sp %>% filter(host_taxon %in% spp_list)
  meta_sp <- meta_sp %>% filter(host_taxon %in% diet_match$host_taxon[!is.na(diet_match$DietBreadth)])
  spp_list <- meta_sp$host_taxon
  
  if (length(spp_list) < 4) {
    cat(clade_entry$label, ": fewer than 4 species with complete data — skipping\n")
    return(NULL)
  }
  
  # Host phylogenetic distance (same approach as calc_all_distances())
  tree_sub  <- keep.tip(clade_entry$tree, unique(spp_list))
  phylo_mat <- cophenetic.phylo(tree_sub)[spp_list, spp_list]
  phylo_dist <- as.dist(phylo_mat)
  
  # OLD diet distance: coarse categorical Diet column, Gower distance
  # (mirrors the pre-revision calc_all_distances() diet block)
  old_diet_df <- meta_sp %>% select(Diet) %>% mutate(Diet = as.factor(as.character(Diet)))
  rownames(old_diet_df) <- spp_list
  old_diet_dist <- daisy(old_diet_df, metric = "gower")
  
  # NEW diet distance: TetrapodTraits 10-component binary matrix, Jaccard
  # (mirrors the current calc_all_distances() diet block; see Section 2)
  new_diet_df <- meta_sp %>% select(all_of(diet_binary_cols))
  rownames(new_diet_df) <- spp_list
  new_diet_dist <- vegan::vegdist(new_diet_df, method = "jaccard", binary = TRUE)
  
  # Mantel tests: diet distance vs. phylogenetic distance
  mantel_old <- vegan::mantel(old_diet_dist, phylo_dist, method = "pearson", permutations = 999)
  mantel_new <- vegan::mantel(new_diet_dist, phylo_dist, method = "pearson", permutations = 999)
  
  cat("\n---", clade_entry$label, "( n =", length(spp_list), "species ) ---\n")
  cat("  OLD (coarse Gower)      diet-phylogeny Mantel r =", round(mantel_old$statistic, 3),
      ", p =", mantel_old$signif, "\n")
  cat("  NEW (TetrapodTraits)    diet-phylogeny Mantel r =", round(mantel_new$statistic, 3),
      ", p =", mantel_new$signif, "\n")
  
  tibble(
    dataset        = clade_entry$label,
    n_species      = length(spp_list),
    mantel_r_old   = mantel_old$statistic,
    mantel_p_old   = mantel_old$signif,
    mantel_r_new   = mantel_new$statistic,
    mantel_p_new   = mantel_new$signif,
    r_reduction    = mantel_old$statistic - mantel_new$statistic
  )
}

mantel_results <- bind_rows(lapply(mantel_clades, run_diet_phylo_mantel))

cat("\n=== SUMMARY: diet-phylogeny collinearity, old vs. new diet scheme ===\n")
print(mantel_results)

f_mantel <- file.path(dir_mrm, "diet_phylogeny_mantel_validation.csv")
write_csv(mantel_results, f_mantel)
cat("\nSaved:", f_mantel, "\n")


#########

################################################################################
################################################################################
################################################################################


################################################################################
# SECTION 11: PHYLOGENETIC RANGE DIAGNOSTIC — amphibians vs. squamates
#             (Response to Reviewer 2, L263-4/268-9)
#
# REVIEWER CONCERN: There is more evolutionary distance spanned by the
# amphibian dataset than by the squamates in the wild dataset. Because
# MRM/Mantel-type tests derive power from variance in the predictor matrix,
# a compressed range of phylogenetic distances in squamates could reduce
# power to detect a phylogeny effect independent of whether the underlying
# phylosymbiosis relationship is equally strong per unit of divergence.
#
# WHAT THIS SECTION DOES:
#   1. Extracts species-level pairwise cophenetic (phylogenetic) distances
#      for each clade x marker combination, in RAW branch-length units
#      (i.e., NOT rescaled to [0,1] as in calc_all_distances() -- rescaling
#      would artificially equalize apparent range and defeat the purpose of
#      this check).
#   2. Reports summary statistics (min, max, mean, SD, CV) side by side for
#      amphibians vs. squamates -- becomes Table S_range.
#   3. Plots the pairwise distance distributions (violin + jitter) by clade,
#      faceted by marker -- becomes Fig S_range_a.
#   4. Recomputes the full (non-bootstrapped) wild-sample phylo vs. community
#      dissimilarity scatter per clade x marker x dissimilarity metric, so
#      the actual relationship each MRM model is fit to is visible directly
#      -- becomes Fig S_range_b.
#
# NOTE: This is descriptive/diagnostic, not a substitute for the subsampling
# robustness check (planned separately) -- it establishes whether restriction
# of phylogenetic range is a plausible contributor to the amphibian/squamate
# asymmetry in bacterial Jaccard MRM results (Section on bacterial
# phylosymbiosis, L263-269) before that more expensive analysis is run.
################################################################################
cat("\n\n################################################################\n")
cat("### SECTION 11: PHYLOGENETIC RANGE DIAGNOSTIC (Reviewer 2)\n")
cat("################################################################\n")

dir_range_diag <- file.path(dir_mrm, "phylo_range_diagnostic")
if (!dir.exists(dir_range_diag)) dir.create(dir_range_diag, recursive = TRUE)
cat("Range diagnostic outputs →", dir_range_diag, "\n")

# ── 11.1: Species-level pairwise phylogenetic distances (raw units) ──────────
# One row per unique species pair per clade x marker -- NOT per sample pair,
# to avoid pseudoreplicating the same species-level distance across many
# samples and inflating apparent n.
get_pairwise_phylo_dist <- function(tree, meta, orders, clade_label, marker_label) {
  spp <- meta %>%
    filter(Clade_Order %in% orders) %>%
    distinct(host_taxon) %>%
    pull(host_taxon)
  spp <- intersect(spp, tree$tip.label)
  
  if (length(spp) < 3) {
    cat(marker_label, clade_label, ": fewer than 3 tree-matched species — skipping\n")
    return(NULL)
  }
  
  tree_sub  <- keep.tip(tree, unique(spp))
  phylo_mat <- cophenetic.phylo(tree_sub)[spp, spp]
  d         <- as.dist(phylo_mat)
  
  tibble(
    marker  = marker_label,
    clade   = clade_label,
    n_spp   = length(spp),
    phylo_dist = as.vector(d)
  )
}

range_diag_combos <- list(
  list(tree = host_tree_16S_full, meta = meta_16S_wild, orders = amphibian_orders, clade = "Amphibian", marker = "16S"),
  list(tree = host_tree_16S_full, meta = meta_16S_wild, orders = reptile_orders,   clade = "Squamate",  marker = "16S"),
  list(tree = host_tree_ITS_full, meta = meta_ITS_wild, orders = amphibian_orders, clade = "Amphibian", marker = "ITS"),
  list(tree = host_tree_ITS_full, meta = meta_ITS_wild, orders = reptile_orders,   clade = "Squamate",  marker = "ITS")
)

phylo_range_long <- bind_rows(lapply(range_diag_combos, function(x) {
  get_pairwise_phylo_dist(x$tree, x$meta, x$orders, x$clade, x$marker)
}))

# ── 11.2: Summary table (min/max/mean/SD/CV) ──────────────────────────────────
phylo_range_summary <- phylo_range_long %>%
  group_by(marker, clade) %>%
  summarise(
    n_species        = first(n_spp),
    n_pairs          = n(),
    min_dist         = min(phylo_dist),
    max_dist         = max(phylo_dist),
    range_dist       = max_dist - min_dist,
    mean_dist        = mean(phylo_dist),
    sd_dist          = sd(phylo_dist),
    cv_dist          = sd_dist / mean_dist,
    .groups = "drop"
  ) %>%
  arrange(marker, clade)

cat("\n--- Pairwise phylogenetic distance summary (raw branch-length units) ---\n")
print(as.data.frame(phylo_range_summary))

f_range_summary <- file.path(dir_range_diag, "phylo_range_summary.csv")
write_csv(phylo_range_summary, f_range_summary)
cat("Saved:", f_range_summary, "\n")

# ── 11.3: Distribution plot — pairwise phylo distance, clade x marker ────────
clade_colors <- c("Amphibian" = "#0052B2", "Squamate" = "#009E73")

fig_range_dist <- ggplot(phylo_range_long, aes(x = clade, y = phylo_dist, fill = clade)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.35, color = "grey20") +
  facet_wrap(~marker, scales = "free_y") +
  scale_fill_manual(values = clade_colors) +
  labs(x = NULL, y = "Pairwise phylogenetic distance\n(raw cophenetic branch length)",
       title = "Phylogenetic distance range: amphibians vs. squamates",
       subtitle = "Species-level pairwise distances underlying each clade's MRM phylogeny predictor") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 11, face = "bold"),
    plot.subtitle   = element_text(size = 8.5, color = "grey30")
  )

f_range_fig <- file.path(dir_range_diag, "phylo_range_distribution.pdf")
ggsave(f_range_fig, fig_range_dist, width = 8, height = 5, device = cairo_pdf, dpi = 400)
cat("Saved:", f_range_fig, "\n")

# ── 11.4: Scatter — phylo distance vs. community dissimilarity (full wild data) ──
# Uses calc_all_distances() directly on the full wild sample set (no
# bootstrap subsampling) per clade x marker x dissimilarity metric, so the
# actual relationship underlying each MRM fit is visible. This is NOT a
# statistical test -- it's the visual companion to 11.2/11.3.
#
# NOTE: assumes sample IDs are accessible via phyloseq::sample_names(ps_wild)
# and align positionally with rows of meta_wild (true throughout the rest of
# this script, e.g. Section 6/7 bootstrap functions) -- adjust if your
# sample_id lives in a metadata column instead.
build_scatter_df <- function(ps_wild, meta_wild, tree, orders, clade_label, marker_label, dissim_name, dissim_method) {
  ids_sub <- sample_names(ps_wild)[meta_wild$Clade_Order %in% orders]
  if (length(ids_sub) < 3) return(NULL)
  
  d <- calc_all_distances(ps_wild, ids_sub, tree, dissim_method = dissim_method)
  
  tibble(
    marker  = marker_label,
    clade   = clade_label,
    dissim  = dissim_name,
    phylo   = as.vector(d$phylo),
    comm    = as.vector(d$comm)
  )
}

scatter_specs <- list(
  list(ps = ps_16S_wild, meta = meta_16S_wild, tree = host_tree_16S_full, orders = amphibian_orders, clade = "Amphibian", marker = "16S"),
  list(ps = ps_16S_wild, meta = meta_16S_wild, tree = host_tree_16S_full, orders = reptile_orders,   clade = "Squamate",  marker = "16S"),
  list(ps = ps_ITS_wild, meta = meta_ITS_wild, tree = host_tree_ITS_full, orders = amphibian_orders, clade = "Amphibian", marker = "ITS"),
  list(ps = ps_ITS_wild, meta = meta_ITS_wild, tree = host_tree_ITS_full, orders = reptile_orders,   clade = "Squamate",  marker = "ITS")
)

scatter_df <- bind_rows(lapply(scatter_specs, function(x) {
  bind_rows(
    build_scatter_df(x$ps, x$meta, x$tree, x$orders, x$clade, x$marker, "Bray-Curtis", "bray"),
    build_scatter_df(x$ps, x$meta, x$tree, x$orders, x$clade, x$marker, "Jaccard",     "jaccard")
  )
}))

fig_range_scatter <- ggplot(scatter_df, aes(x = phylo, y = comm, color = clade)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  facet_grid(dissim ~ marker, scales = "free") +
  scale_color_manual(values = clade_colors) +
  labs(x = "Phylogenetic distance (rescaled [0,1], as used in MRM)",
       y = "Community dissimilarity",
       color = "Clade",
       title = "Phylogenetic distance vs. community dissimilarity, full wild dataset",
       subtitle = "Non-bootstrapped; visualizes the relationship each clade's MRM phylogeny term is fit to") +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey30")
  )

f_scatter_fig <- file.path(dir_range_diag, "phylo_vs_community_scatter.pdf")
ggsave(f_scatter_fig, fig_range_scatter, width = 9, height = 7, device = cairo_pdf, dpi = 400)
cat("Saved:", f_scatter_fig, "\n")

cat("\n=== SECTION 11 COMPLETE ===\n")
cat("Range summary table:", f_range_summary, "\n")
cat("Distribution figure :", f_range_fig, "\n")
cat("Scatter figure      :", f_scatter_fig, "\n")


################################################################################
################################################################################
################################################################################


################################################################################
# SECTION 12: PHYLOGENETIC RANGE SUBSAMPLING — confirmatory check
#             (Response to Reviewer 2, L263-4/268-9)
#
# RATIONALE: Section 11 established that (a) squamates span less phylogenetic
# range/variance than amphibians even after [0,1] rescaling, but (b) this does
# not appear to explain the squamate Bray-Curtis-significant /
# Jaccard-non-significant asymmetry, because the phylogenetic predictor is
# identical between those two models and phylogeny IS detected under
# Bray-Curtis (98% of iterations; partial R2 = 0.109 -- higher than the
# amphibian Bray-Curtis partial R2 of 0.076).
#
# This section provides an empirical, confirmatory test of that logical
# argument: if restricted phylogenetic range were suppressing the squamate
# Jaccard phylogeny effect, then WIDENING that range (by subsampling toward
# more phylogenetically dispersed squamate species) should be more likely to
# reveal a phylogeny effect under Jaccard. If the Jaccard result remains
# non-significant even once range/variance is increased toward amphibian-like
# levels, that further supports the Section 11 conclusion that the pattern is
# not a range/power artifact.
#
# METHOD: Farthest-point (max-min) greedy selection on raw species-level
# cophenetic distances -- deterministic given a seed, and directly maximizes
# phylogenetic dispersion of the selected subset (as opposed to ad hoc
# "one species per family" selection). Run across several target subsample
# sizes as a sensitivity check, not just one.
#
# OPTIONAL MIRROR TEST (off by default): restricts amphibians to a tight,
# low-range phylogenetic cluster (opposite selection logic) to test whether
# the amphibian phylogeny effect weakens when range is compressed toward
# squamate-like levels.
################################################################################
cat("\n\n################################################################\n")
cat("### SECTION 12: PHYLOGENETIC RANGE SUBSAMPLING (Reviewer 2)\n")
cat("################################################################\n")

dir_subsample <- file.path(dir_mrm, "range_subsampling")
if (!dir.exists(dir_subsample)) dir.create(dir_subsample, recursive = TRUE)
cat("Subsampling outputs →", dir_subsample, "\n")

# ── 12.1: Selection algorithms ────────────────────────────────────────────────

# Greedy farthest-point (max-min) selection: maximizes phylogenetic dispersion
# of the selected subset. Starts from the single most divergent pair, then
# repeatedly adds whichever remaining species maximizes its minimum distance
# to the already-selected set.
farthest_point_subsample <- function(dist_mat, target_n, seed = 1) {
  set.seed(seed)
  sp <- rownames(dist_mat)
  if (target_n >= length(sp)) return(sp)
  
  # Seed the selection with the most divergent pair
  idx_max <- which(dist_mat == max(dist_mat), arr.ind = TRUE)[1, ]
  selected <- unique(c(sp[idx_max[1]], sp[idx_max[2]]))
  
  while (length(selected) < target_n) {
    remaining <- setdiff(sp, selected)
    min_dist_to_selected <- sapply(remaining, function(s) {
      min(dist_mat[s, selected])
    })
    # Ties broken randomly (seeded) rather than by row order
    best <- remaining[min_dist_to_selected == max(min_dist_to_selected)]
    selected <- c(selected, sample(best, 1))
  }
  return(selected)
}

# Greedy tightest-cluster selection (mirror test): minimizes phylogenetic
# dispersion. Starts from a random species, then repeatedly adds whichever
# remaining species produces the smallest increase in the cluster's maximum
# pairwise distance (i.e., grows the most compact possible cluster).
tightest_cluster_subsample <- function(dist_mat, target_n, seed = 1) {
  set.seed(seed)
  sp <- rownames(dist_mat)
  if (target_n >= length(sp)) return(sp)
  
  selected <- sample(sp, 1)
  while (length(selected) < target_n) {
    remaining <- setdiff(sp, selected)
    resulting_max <- sapply(remaining, function(s) {
      cand <- c(selected, s)
      max(dist_mat[cand, cand])
    })
    best <- remaining[resulting_max == min(resulting_max)]
    selected <- c(selected, sample(best, 1))
  }
  return(selected)
}

# ── 12.2: Build subsampled species sets + confirm range/variance shift ───────
# Reuses the same raw (non-rescaled) cophenetic-distance logic as Section 11.
build_cophenetic_dist_mat <- function(tree, species_pool) {
  spp <- intersect(unique(species_pool), tree$tip.label)
  tree_sub <- keep.tip(tree, spp)
  cophenetic.phylo(tree_sub)[spp, spp]
}

summarize_dist_vec <- function(d) {
  v <- as.vector(as.dist(d))
  tibble(
    n_species  = attr(as.dist(d), "Size"),
    n_pairs    = length(v),
    min_dist   = min(v), max_dist = max(v), range_dist = max(v) - min(v),
    mean_dist  = mean(v), sd_dist = sd(v), cv_dist = sd(v) / mean(v)
  )
}

target_sizes_squamate <- c(12, 15, 18)   # sensitivity range; adjust as needed
run_amphibian_mirror  <- FALSE           # set TRUE to also run the reverse test

squamate_specs <- list(
  list(marker = "16S", ps = ps_16S_wild, meta = meta_16S_wild, tree = host_tree_16S_full),
  list(marker = "ITS", ps = ps_ITS_wild, meta = meta_ITS_wild, tree = host_tree_ITS_full)
)

# ── 12.3: Run bootstrap MRM on each subsampled squamate species set ──────────
run_subsample_test <- function(spec, target_n, selection = "farthest") {
  
  meta_sq <- spec$meta %>% filter(Clade_Order %in% reptile_orders)
  spp_pool <- unique(meta_sq$host_taxon)
  spp_pool <- intersect(spp_pool, spec$tree$tip.label)
  
  if (length(spp_pool) < target_n) {
    cat(spec$marker, ": fewer species available (", length(spp_pool),
        ") than target_n =", target_n, "— skipping\n")
    return(NULL)
  }
  
  full_dist_mat <- build_cophenetic_dist_mat(spec$tree, spp_pool)
  
  sel_spp <- if (selection == "farthest") {
    farthest_point_subsample(full_dist_mat, target_n, seed = 52325)
  } else {
    tightest_cluster_subsample(full_dist_mat, target_n, seed = 52325)
  }
  
  # Confirm the subsample actually shifted range/variance in the intended
  # direction before spending bootstrap iterations on it
  diag_before <- summarize_dist_vec(full_dist_mat) %>% mutate(stage = "full_squamate")
  diag_after  <- summarize_dist_vec(full_dist_mat[sel_spp, sel_spp]) %>% mutate(stage = paste0("subsample_n", target_n))
  diag_tbl <- bind_rows(diag_before, diag_after) %>%
    mutate(marker = spec$marker, selection = selection, .before = 1)
  
  cat("\n---", spec$marker, "| target_n =", target_n, "| selection =", selection, "---\n")
  print(as.data.frame(diag_tbl))
  
  # Restrict metadata/phyloseq object to the selected species, then reuse the
  # existing bootstrap machinery unchanged
  meta_sub <- meta_sq %>% filter(host_taxon %in% sel_spp)
  ps_sub   <- prune_samples(meta_sub$sample_name, spec$ps)
  
  boot_results <- list()
  for (dissim in dissim_methods) {
    key <- paste(spec$marker, "SquamateSubsample", selection, target_n, dissim$name, sep = "_")
    cat("Running", n_bootstrap, "bootstrap iterations —", key, "\n")
    
    boot_df <- run_bootstrap(
      boot_fn       = bootstrap_mrm,
      label         = key,
      meta          = meta_sub,
      ps            = ps_sub,
      tree          = spec$tree,
      dissim_method = dissim$method
    )
    
    if (nrow(boot_df) == 0) {
      cat("  No successful iterations for", key, "— skipping summary\n")
      next
    }
    
    summary_df <- summarize_bootstrap(boot_df) %>%
      rename_predictors(analysis_type = "primary") %>%
      mutate(marker = spec$marker, selection = selection,
             target_n = target_n, dissim_metric = dissim$name)
    
    write.csv(boot_df,
              file.path(dir_subsample, paste0(key, "_bootstrap.csv")),
              row.names = FALSE)
    write.csv(summary_df,
              file.path(dir_subsample, paste0(key, "_summary.csv")),
              row.names = FALSE)
    
    boot_results[[dissim$name]] <- summary_df
  }
  
  list(diagnostics = diag_tbl, results = bind_rows(boot_results))
}

# ── 12.4: Execute across target sizes for both markers ───────────────────────
subsample_output <- list()
for (spec in squamate_specs) {
  for (n_target in target_sizes_squamate) {
    key <- paste(spec$marker, n_target, sep = "_n")
    subsample_output[[key]] <- run_subsample_test(spec, n_target, selection = "farthest")
  }
}

# ── 12.5: Combine diagnostics + results, compare to full-dataset baseline ────
all_diag <- bind_rows(lapply(subsample_output, function(x) x$diagnostics))
f_diag <- file.path(dir_subsample, "subsample_range_diagnostics.csv")
write_csv(all_diag, f_diag)
cat("\nSaved:", f_diag, "\n")

all_results <- bind_rows(lapply(subsample_output, function(x) x$results))

# Pull the ORIGINAL (full-dataset) squamate summaries from Section 6 results
# for direct side-by-side comparison. Assumes results_primary keys follow the
# "<marker>_Reptile_<dissim>" convention established in Section 6.
baseline_rows <- bind_rows(lapply(squamate_specs, function(spec) {
  bind_rows(lapply(dissim_methods, function(dissim) {
    key <- paste(spec$marker, "Reptile", dissim$name, sep = "_")
    if (is.null(results_primary[[key]])) return(NULL)
    results_primary[[key]]$summary %>%
      rename_predictors(analysis_type = "primary") %>%
      filter(predictor == "Phylogeny") %>%
      mutate(marker = spec$marker, selection = "full_dataset",
             target_n = length(unique(spec$meta %>%
                                        filter(Clade_Order %in% reptile_orders) %>%
                                        pull(host_taxon))),
             dissim_metric = dissim$name)
  }))
}))

comparison_tbl <- bind_rows(baseline_rows, all_results %>% filter(predictor == "Phylogeny")) %>%
  select(marker, selection, target_n, dissim_metric, predictor,
         mean_coef, mean_partial_r2, mean_pval, prop_sig_0.05, mean_r2_full) %>%
  arrange(marker, dissim_metric, target_n)

f_comparison <- file.path(dir_subsample, "phylogeny_subsample_vs_full_comparison.csv")
write_csv(comparison_tbl, f_comparison)
cat("Saved:", f_comparison, "\n")

cat("\n--- Phylogeny predictor: full dataset vs. range-widened subsamples ---\n")
print(as.data.frame(comparison_tbl))

###
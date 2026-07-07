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
#   Part 1  — Predictor label lookup
#   Part 2  — Primary MRM (wild samples only)
#   Part 3  — Wild + Managed MRM
#   Part 4  — Null phylogeny tests
#   Part 5  — Figures and supplementary table
#
# OUTPUT STRUCTURE:
#   output/mrm/                    <- main results (figures, summary table)
#   output/mrm/extended_outputs/   <- raw bootstrap CSVs/RDS per analysis
#   output/mrm/null_tests/         <- null phylogeny test outputs
#
# Alexander Rurik
################################################################################

# ========================== LOAD PACKAGES =====================================
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

# ========================== SET PARAMETERS ====================================

set.seed(52325)
n_bootstrap <- 100

# ========================== PROJECT DIRECTORIES ===============================
library(here)
dir_processed <- here("data", "processed")
dir_trees     <- here("data", "raw", "trees")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")
dir_mrm       <- here("output", "mrm")

# ========================== COLOR PALETTE =====================================
# 16S = darker shades; ITS = lighter/brighter shades
dataset_colors <- c(
  "16S Amphibians" = "#0052B2",
  "16S Reptiles"   = "#009E73",
  "ITS Amphibians" = "#11E4E9",
  "ITS Reptiles"   = "#69FDAA"
)

################################################################################
# SETUP: Load data and harmonize species names
################################################################################
cat("=== LOADING DATA ===\n")

ps_16S <- readRDS(file.path(dir_processed, "16S_abs_final_fecal.rds"))
ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_fecal.rds"))

cat("16S samples loaded:", nsamples(ps_16S), "\n")
cat("ITS samples loaded:", nsamples(ps_ITS), "\n")

# ======================== OUTPUT DIRECTORIES ==================================
dir_extended     <- file.path(dir_mrm, "extended_outputs")
dir_null         <- file.path(dir_mrm, "null_tests")

for (d in c(dir_mrm, dir_extended, dir_null)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("Output directories created:\n")
cat("  Extended outputs →", dir_extended,     "\n")
cat("  Null tests       →", dir_null,         "\n")

# == RECODE MANAGEMENT STATUS: "Captive" → "Managed" ===========================
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

# ==== SPECIES NAME HARMONIZATION ==============================================

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
# PRIMARY ANALYSIS SETUP: Wild samples only
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

# === FIX GPS COORDINATE ERRORS ================================================

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

# === LOAD HOST TREES ==========================================================

cat("\n=== Loading host trees ===\n")

host_tree_16S_full <- read.tree(file.path(dir_trees, "16S_species_for_timetree_clean.nwk"))
host_tree_ITS_full <- read.tree(file.path(dir_trees, "ITS_species_for_timetree_clean.nwk"))

for (tr in c("host_tree_16S_full", "host_tree_ITS_full")) {
  t <- get(tr)
  t$tip.label <- gsub("_", " ", t$tip.label)
  assign(tr, t)
}

cat("16S tree tips:", length(host_tree_16S_full$tip.label), "\n")
cat("ITS tree tips:", length(host_tree_ITS_full$tip.label), "\n")

# ======== DEFINE CLADES AND ANALYSIS COMBINATIONS =============================
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

# ===== VALIDATE SPECIES COVERAGE ===============================================

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
# HELPER FUNCTIONS
################################################################################

# == Rescale distance matrix to [0, 1] ─────────────────────────────────────────
rescale_dist <- function(d) {
  m <- as.matrix(d)
  as.dist((m - min(m)) / (max(m) - min(m)))
}

# == BOOTSTRAP SUMMARY =========================================================

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

# == Distance matrices: PRIMARY (wild-only) =====================================
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
  
  # Diet dissimilarity (Gower; categorical)
  diet_df      <- meta_sub %>% select(Diet) %>%
    mutate(Diet = as.factor(as.character(Diet)))
  rownames(diet_df) <- sample_ids
  diet_dist    <- daisy(diet_df, metric = "gower")
  
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

# == Distance matrices: MANAGED + WILD ==========================================
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
  
  # Diet dissimilarity
  diet_df   <- meta_sub %>% select(Diet) %>%
    mutate(Diet = as.factor(as.character(Diet)))
  rownames(diet_df) <- sample_ids
  diet_dist <- daisy(diet_df, metric = "gower")
  
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

# == Bootstrap function: PRIMARY (wild-only) ==================================
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

# == Bootstrap function: MANAGED + WILD =========================================
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

# == Generic bootstrap runner =================================================
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
# PART 1: Predictor label lookup
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
# PART 2: PRIMARY MRM ANALYSES — Wild samples only
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
saveRDS(results_primary, file.path(dir_extended, "primary_all_results.rds"))


################################################################################
# PART 3: WILD + MANAGED MRM
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
saveRDS(results_managed, file.path(dir_extended, "managed_wild_all_results.rds"))


################################################################################
# PART 4: NULL PHYLOGENY TESTS
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

# == Null test figure: 2×2 multipanel =========================================
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

cairo_pdf(file.path(dir_figures, "null_phylogeny_multipanel.pdf"), width = 14, height = 12)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
plot_null_panel(null_results_list[["16S_Amphibian"]], "A)", "16S Amphibians")
plot_null_panel(null_results_list[["16S_Reptile"]],   "B)", "16S Reptiles")
plot_null_panel(null_results_list[["ITS_Amphibian"]], "C)", "ITS Amphibians")
plot_null_panel(null_results_list[["ITS_Reptile"]],   "D)", "ITS Reptiles")
dev.off()
cat("Saved:", file.path(dir_figures, "null_phylogeny_multipanel.pdf"), "\n")


################################################################################
# PART 4b: NULL PHYLOGENY TESTS — JACCARD
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

# == Null test figure: 2×2 multipanel (Jaccard) ===============================
cairo_pdf(file.path(dir_figures, "null_phylogeny_multipanel_Jaccard.pdf"), width = 14, height = 12)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
plot_null_panel(null_results_list_jaccard[["16S_Amphibian"]], "A)", "16S Amphibians")
plot_null_panel(null_results_list_jaccard[["16S_Reptile"]],   "B)", "16S Reptiles")
plot_null_panel(null_results_list_jaccard[["ITS_Amphibian"]], "C)", "ITS Amphibians")
plot_null_panel(null_results_list_jaccard[["ITS_Reptile"]],   "D)", "ITS Reptiles")
dev.off()
cat("Saved:", file.path(dir_figures, "null_phylogeny_multipanel_Jaccard.pdf"), "\n")


################################################################################
# PART 5: FIGURES AND SUPPLEMENTARY TABLE
################################################################################
cat("\n=== PART 5: FIGURES AND SUPPLEMENTARY TABLE ===\n")

# == Helper: prepare summary for plotting =====================================
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
          prop_sig_label = paste0(round(prop_sig_0.05 * 100), "%")
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

# == Shared y-axis limits =======================================================
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

# == Sample size annotation =====================================================
build_n_label <- function(df) {
  df %>%
    group_by(dataset, dissim, analysis) %>%
    summarise(n_label = paste0("n=", round(unique(mean_n_species)), " spp"),
              .groups = "drop")
}

n_labels_primary <- build_n_label(plot_primary)
n_labels_mw      <- build_n_label(plot_managed_wild)

# ==Figure constructor ==========================================================
# Bars show mean MRM coefficient ± SD across 100 bootstrap iterations.
# Labels above bars show the proportion of iterations with BH-adjusted p < 0.05
# (e.g. "94%"), providing a measure of effect consistency across species
# subsamples. Predictor order: Diet, Habitat, Geography/Management status, Phylogeny.
make_mrm_plot <- function(df, ylimits, title_str, n_df = NULL, show_legend = TRUE) {
  
  pred_levels <- c("Diet", "Habitat", "Geography", "Management status", "Phylogeny")
  df <- df %>%
    mutate(predictor = factor(predictor,
                              levels = intersect(pred_levels, unique(predictor))))
  
  # Prop-sig label y-position: just above the top of the error bar
  sig_offset <- diff(ylimits) * 0.025
  
  p <- ggplot(df, aes(x = predictor, y = mean_coef, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_coef - sd_coef, ymax = mean_coef + sd_coef),
                  position = position_dodge(width = 0.8), width = 0.25,
                  linewidth = 0.5) +
    # Proportion of iterations significant, displayed instead of asterisks
    geom_text(aes(label = prop_sig_label,
                  y     = mean_coef + sd_coef + sig_offset),
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

# == Build individual figures =====================================================
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

# == Multipanel figures ==============================================================

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

# == Save all figures ==========================================================
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

# == Supplementary table =======================================================
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

supp_fname <- file.path(dir_tables, "MRM_supplementary_table.csv")
write.csv(supp_table, supp_fname, row.names = FALSE)
cat("Saved supplementary table:", supp_fname, "\n")

cat("\n--- Table preview ---\n")
print(as.data.frame(supp_table), max = 200)


################################################################################
################################################################################
################################################################################
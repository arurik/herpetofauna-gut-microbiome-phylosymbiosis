################################################################################
# DEMO: MRM Phylosymbiosis Workflow (Simulated Data)
#
# Purpose: minimal, self-contained demonstration of the bootstrap MRM
# phylosymbiosis testing approach used in this study (see scripts/10_mrm_
# phylosymbiosis.R for the full analysis on the published dataset). Runs
# entirely on simulated data generated in-memory. No external files or
# downloads required.
#
# NOTE: this is NOT a reduced version of the real dataset. Results are not
# biologically meaningful and will not resemble the manuscript's findings.
# This script exists solely to verify the core analysis code runs correctly.
#
# ANALYSIS ORDER:
#   Part 1 - Simulate host phylogeny, community, GPS, and predictor data
#   Part 2 - Bootstrap MRM (community ~ phylogeny + geography + diet + habitat)
#   Part 3 - Summarize and print results
#
################################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(ecodist)
  library(cluster)
  library(geosphere)
})

# ========================== SET PARAMETERS ====================================
set.seed(52325)
n_species   <- 10
n_per_sp    <- 3
n_otus      <- 15
n_bootstrap <- 5   # reduced from 100 in the real analysis, for demo speed

################################################################################
# PART 1: Simulate host phylogeny, community, GPS, and predictor data
################################################################################
cat("=== SIMULATING DEMO DATASET ===\n")

demo_tree <- rtree(n_species)
demo_tree$tip.label <- paste0("sp", seq_len(n_species))
cat("Simulated host phylogeny:", n_species, "tips\n")

demo_counts <- matrix(
  rpois(n_species * n_per_sp * n_otus, lambda = 20),
  nrow = n_species * n_per_sp, ncol = n_otus,
  dimnames = list(
    paste0(rep(demo_tree$tip.label, each = n_per_sp), "_", seq_len(n_per_sp)),
    paste0("OTU", seq_len(n_otus))
  )
)

# One set of "collection coordinates" per species (roughly continental US
# bounding box), same as our real pipeline uses per-sample GPS
species_gps <- data.frame(
  species = demo_tree$tip.label,
  gps_w   = runif(n_species, min = -100, max = -80),
  gps_n   = runif(n_species, min = 28, max = 40)
)

demo_meta <- data.frame(
  sample_id = rownames(demo_counts),
  species   = rep(demo_tree$tip.label, each = n_per_sp),
  Diet      = as.factor(sample(c("Carnivore", "Herbivore", "Omnivore"),
                               n_species, replace = TRUE)[
                                 rep(seq_len(n_species), each = n_per_sp)]),
  Habitat   = as.factor(sample(c("Aquatic", "Terrestrial", "Arboreal"),
                               n_species, replace = TRUE)[
                                 rep(seq_len(n_species), each = n_per_sp)])
) %>%
  left_join(species_gps, by = "species")

cat("Simulated community:", nrow(demo_counts), "samples,",
    ncol(demo_counts), "OTUs across", n_species, "host species\n\n")

################################################################################
# HELPER FUNCTIONS (identical to scripts/10_mrm_phylosymbiosis.R)
################################################################################

# == Rescale distance matrix to [0, 1] ─────────────────────────────────────────
rescale_dist <- function(d) {
  m <- as.matrix(d)
  as.dist((m - min(m)) / (max(m) - min(m)))
}

################################################################################
# PART 2: Bootstrap MRM (one individual per species per iteration)
################################################################################
cat("=== RUNNING BOOTSTRAP MRM (", n_bootstrap, "iterations) ===\n")

comm_dist_full  <- vegdist(demo_counts, method = "bray")
phylo_dist_full <- rescale_dist(cophenetic(demo_tree))

boot_results <- vector("list", n_bootstrap)

for (i in seq_len(n_bootstrap)) {
  
  # one random individual per species, mirroring the real bootstrap design
  boot_ids <- demo_meta %>%
    group_by(species) %>%
    slice_sample(n = 1) %>%
    pull(sample_id)
  
  meta_sub <- demo_meta %>%
    filter(sample_id %in% boot_ids) %>%
    arrange(match(sample_id, boot_ids))
  
  comm_sub <- as.matrix(comm_dist_full)[boot_ids, boot_ids]
  
  phylo_sub <- as.matrix(phylo_dist_full)[
    match(meta_sub$species, demo_tree$tip.label),
    match(meta_sub$species, demo_tree$tip.label)
  ]
  dimnames(phylo_sub) <- list(boot_ids, boot_ids)
  
  # Geographic distance (great-circle km, rescaled to [0,1]; same approach
  # as calc_all_distances() in scripts/10_mrm_phylosymbiosis.R)
  coord_mat <- as.matrix(meta_sub[, c("gps_w", "gps_n")])
  g_mat     <- geosphere::distm(coord_mat, fun = geosphere::distGeo) / 1000
  dimnames(g_mat) <- list(boot_ids, boot_ids)
  geo_sub   <- rescale_dist(as.dist(g_mat))
  
  diet_dist    <- daisy(meta_sub["Diet"],    metric = "gower")
  habitat_dist <- daisy(meta_sub["Habitat"], metric = "gower")
  
  fit <- MRM(as.dist(comm_sub) ~ as.dist(phylo_sub) + geo_sub +
               diet_dist + habitat_dist,
             nperm = 99)
  
  boot_results[[i]] <- fit$coef %>%
    as.data.frame() %>%
    rownames_to_column("predictor") %>%
    rename(coef = `as.dist(comm_sub)`) %>%
    mutate(iteration = i, r2_full = fit$r.squared[1])
}

################################################################################
# PART 3: Summarize results
################################################################################
cat("\n=== SUMMARY ACROSS BOOTSTRAP ITERATIONS ===\n")

boot_df <- bind_rows(boot_results) %>% filter(predictor != "Int")

summary_table <- boot_df %>%
  mutate(predictor = recode(predictor,
                            "as.dist(phylo_sub)" = "Phylogeny",
                            "geo_sub"            = "Geography",
                            "diet_dist"          = "Diet",
                            "habitat_dist"       = "Habitat"
  )) %>%
  group_by(predictor) %>%
  summarise(
    mean_coef = mean(coef, na.rm = TRUE),
    sd_coef   = sd(coef,   na.rm = TRUE),
    mean_pval = mean(pval, na.rm = TRUE),
    mean_r2   = mean(r2_full, na.rm = TRUE),
    .groups   = "drop"
  )

print(as.data.frame(summary_table))

cat("\nDemo complete. Runtime should be well under 1 minute.\n")
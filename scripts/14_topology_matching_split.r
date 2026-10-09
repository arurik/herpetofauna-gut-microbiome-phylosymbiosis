################################################################################
# Topology-Based Phylosymbiosis Testing: Matching Split Distance
#
# Companion analysis to script 10 (MRM phylosymbiosis test) — addresses Reviewer #2,
# who requested a second class of phylosymbiosis test per Lim & Bordenstein
# 2020 (PRSB). MRM/partial Mantel are matrix-correlation methods; this script
# adds the topological-congruence class, which is mathematically independent
# (uses dendrogram/tree topology only, not pairwise distance magnitudes).
#
# METRIC: Matching Split distance (Bogdanowicz & Giaro 2012), implemented in
# TreeDist (Smith 2020). Selected per Lim & Bordenstein 2020's identification
# of the matching-cluster/split family as more sensitive than Robinson-Foulds
# for detecting partial topological congruence. TreeDist implements only the
# unrooted variant of this metric (no maintained R package I can find implements 
# the rooted "Matching Cluster" version); host and microbiome trees are therefore
# compared as unrooted topologies throughout.
#
# Approach:  Bootstrap 100 iterations, 1 sample per species per iteration
#            (identical resampling logic to bootstrap_mrm() in MRM script) to
#            account for intraspecific variation and uneven sampling across species
# Analyses:  Amphibians & Reptiles | 16S & ITS | Bray-Curtis & Jaccard
#            UPGMA microbiome dendrogram vs. host phylogeny, each iteration
#
# NOTE ON CLADE DEFINITIONS (same as MRM script):
#   Reptilia: Squamata only. Testudines and Crocodilia sampled via
#   cloacal swab (excluded from fecal-only analyses) or are ex situ managed samples
#   Amphibia: Caudata + Anura.
#
# ANALYSIS ORDER:
#   Part 1 — Setup / load data (mirrors MRM script Part 1)
#   Part 2 — Bootstrap Matching Split distance (observed, wild samples only)
#   Part 3 — Null topology tests (random trees vs. observed host tree),
#            Bray-Curtis and Jaccard
#   Part 4 — Supplementary table (sole reported output)
#
# Alexander Rurik
################################################################################

# ========== Load packages ===================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(ape)
  library(phangorn)   # upgma()
  library(TreeDist)   # MatchingSplitDistance()
})

set.seed(52325)
n_bootstrap     <- 100    # matches MRM script
n_null          <- 1000   # random host-tree topologies per combo (Part 3)
n_boot_per_null <- 20     # matches MRM script's null-test asymmetry logic

# ========== Project directories =============================================
library(here)
dir_raw       <- here("data", "raw")
dir_processed <- here("data", "processed")
dir_trees     <- here("data", "raw", "trees")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")
dir_topology      <- file.path(dir_tables, "topology")
dir_topology_null <- file.path(dir_topology, "null_tests")
dir.create(dir_topology,      recursive = TRUE, showWarnings = FALSE)
dir.create(dir_topology_null, recursive = TRUE, showWarnings = FALSE)

################################################################################
# PART 1: SETUP — Load data and host trees
################################################################################
cat("=== LOADING DATA ===\n")

ps_16S <- readRDS(file.path(dir_processed, "16S_abs_final_fecal.rds"))
ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_fecal.rds"))
cat("16S samples loaded:", nsamples(ps_16S), "| ITS samples loaded:", nsamples(ps_ITS), "\n")

# ── Species name harmonization (identical to 260624_mrm_v5_AJR.r) ─────────────
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
  changes <- 0
  for (old in names(replacements)) {
    idx <- meta$host_taxon == old
    if (sum(idx, na.rm = TRUE) > 0) {
      meta$host_taxon[idx] <- replacements[old]
      changes <- changes + sum(idx)
    }
  }
  sample_data(ps) <- meta
  cat(dataset_name, "harmonization: updated", changes, "samples\n")
  return(ps)
}

ps_16S <- apply_replacements(ps_16S, replacements_16S, "16S")
ps_ITS <- apply_replacements(ps_ITS, replacements_ITS, "ITS")

# ── Filter to wild samples only ────────────────────────────────────────────────
filter_to_wild <- function(ps, dataset_name) {
  meta <- data.frame(sample_data(ps))
  wild_ids <- meta$sample_name[meta$env_broad_scale == "Wild"]
  ps_wild  <- prune_samples(wild_ids, ps)
  ps_wild  <- prune_taxa(taxa_sums(ps_wild) > 0, ps_wild)
  cat(dataset_name, "wild filter:", nsamples(ps_wild), "/", nrow(meta), "samples retained\n")
  return(ps_wild)
}

ps_16S_wild <- filter_to_wild(ps_16S, "16S")
ps_ITS_wild <- filter_to_wild(ps_ITS, "ITS")

# ── Load host trees ─────────────────────────────────────────────────────────
# NOTE: host and microbiome trees are compared UNROOTED inside
# bootstrap_topology() below (see METRIC note in header) -- no rooting is
# applied here, avoiding any arbitrary outgroup choice for clade subsets.
host_tree_16S_full <- read.tree(file.path(dir_trees, "host_tree_16S.nwk"))
host_tree_ITS_full <- read.tree(file.path(dir_trees, "host_tree_ITS.nwk"))
host_tree_16S_full$tip.label <- gsub("_", " ", host_tree_16S_full$tip.label)
host_tree_ITS_full$tip.label <- gsub("_", " ", host_tree_ITS_full$tip.label)
cat("16S tree tips:", length(host_tree_16S_full$tip.label),
    "| ITS tree tips:", length(host_tree_ITS_full$tip.label), "\n")

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

meta_16S_wild <- data.frame(sample_data(ps_16S_wild))
meta_ITS_wild <- data.frame(sample_data(ps_ITS_wild))

################################################################################
# HELPER FUNCTIONS
################################################################################

# ── One bootstrap iteration: Matching Split distance, host tree vs. UPGMA ────
bootstrap_topology <- function(meta, ps, tree, dissim_method = "bray") {
  
  sampled    <- meta %>% group_by(host_taxon) %>% slice_sample(n = 1) %>% ungroup()
  spp        <- unique(sampled$host_taxon)
  
  tree_sub <- keep.tip(tree, intersect(tree$tip.label, spp))
  keep_spp <- tree_sub$tip.label
  sampled  <- sampled %>% filter(host_taxon %in% keep_spp)
  
  ps_sub <- prune_samples(sampled$sample_name, ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  comm_dist <- phyloseq::distance(ps_sub, method = dissim_method,
                                  binary = (dissim_method == "jaccard"))
  comm_mat  <- as.matrix(comm_dist)
  
  meta_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
  rownames(comm_mat) <- meta_sub$host_taxon[match(rownames(comm_mat), meta_sub$sample_name)]
  colnames(comm_mat) <- rownames(comm_mat)
  comm_mat  <- comm_mat[keep_spp, keep_spp]
  
  microbe_tree <- upgma(as.dist(comm_mat))
  host_tree    <- unroot(tree_sub)
  microbe_tree <- unroot(microbe_tree)
  
  # Raw (unnormalized) Matching Split distance -- see header METRIC note.
  # normalize=TRUE is not supported for this metric (no analytic maximum);
  # raw units are standard practice (Brooks et al. 2016; Lim & Bordenstein
  # 2020) and sufficient for the empirical null-distribution p-value below,
  # since that comparison only requires observed vs. null trees of the SAME
  # size -- scale cancels out.
  mc_raw <- as.numeric(MatchingSplitDistance(host_tree, microbe_tree, normalize = FALSE))
  
  data.frame(mc_raw = mc_raw, n_species = length(host_tree$tip.label))
}

run_bootstrap_topology <- function(label, meta, ps, tree, dissim_method, n_iter = n_bootstrap) {
  cat("  Bootstrapping:", label, "(", n_iter, "iterations )\n")
  out <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    set.seed(52325 + i)
    out[[i]] <- tryCatch(
      bootstrap_topology(meta, ps, tree, dissim_method),
      error = function(e) { cat("    iter", i, "failed:", conditionMessage(e), "\n"); NULL }
    )
  }
  bind_rows(out)
}

summarize_topology_bootstrap <- function(boot_df) {
  boot_df %>%
    summarise(
      n_iterations   = n(),
      mean_mc_raw    = mean(mc_raw, na.rm = TRUE),
      sd_mc_raw      = sd(mc_raw,   na.rm = TRUE),
      mean_n_species = mean(n_species, na.rm = TRUE),
      sd_n_species   = sd(n_species,   na.rm = TRUE)
    )
}

################################################################################
# PART 2: BOOTSTRAP TOPOLOGY TESTS — WILD SAMPLES ONLY
################################################################################
cat("\n=== PART 2: BOOTSTRAP TOPOLOGY TESTS (WILD) ===\n")

datasets <- list(
  list(name = "16S", ps = ps_16S_wild, meta = meta_16S_wild, tree = host_tree_16S_full),
  list(name = "ITS", ps = ps_ITS_wild, meta = meta_ITS_wild, tree = host_tree_ITS_full)
)

results_topology <- list()

for (dataset in datasets) {
  for (clade in clades) {
    for (dissim in dissim_methods) {
      
      key <- paste(dataset$name, clade$name, dissim$name, sep = "_")
      cat("\n--- TOPOLOGY:", key, "---\n")
      
      meta_clade <- dataset$meta %>% filter(Clade_Order %in% clade$orders)
      cat("Pool:", nrow(meta_clade), "samples |",
          length(unique(meta_clade$host_taxon)), "species\n")
      
      if (nrow(meta_clade) < 10 || length(unique(meta_clade$host_taxon)) < 5) {
        cat("Skipping — too few species\n"); next
      }
      
      ps_clade <- prune_samples(meta_clade$sample_name, dataset$ps)
      
      boot_df <- run_bootstrap_topology(key, meta_clade, ps_clade, dataset$tree, dissim$method)
      write.csv(boot_df, file.path(dir_topology, paste0("topology_", key, "_bootstrap.csv")), row.names = FALSE)

      summary_df <- summarize_topology_bootstrap(boot_df) %>%
        mutate(dataset = dataset$name, clade = clade$name, dissim = dissim$name)
      write.csv(summary_df, file.path(dir_topology, paste0("topology_", key, "_summary.csv")), row.names = FALSE)
      
      results_topology[[key]] <- list(bootstrap = boot_df, summary = summary_df)
      cat("Saved:", key, "| mean Matching Split =", round(summary_df$mean_mc_raw, 3), "\n")
    }
  }
}

saveRDS(results_topology, file.path(dir_topology, "topology_all_results.rds"))
cat("\nTOPOLOGY BOOTSTRAP COMPLETE —", length(results_topology), "analyses\n")

################################################################################
# PART 3: NULL TOPOLOGY TESTS
# Empirical p = proportion of null Matching Split distances <= observed (i.e.,
# how often a random host tree is AS OR MORE congruent with the microbiome
# dendrogram than the real host phylogeny is)
################################################################################
cat("\n=== PART 3: NULL TOPOLOGY TESTS ===\n")

run_topology_null_test <- function(ps, meta, tree, dissim_method = "bray",
                                   n_null = 1000, n_boot_per_null = 20) {
  
  n_spp <- length(unique(meta$host_taxon))
  cat("  Species:", n_spp, "| Null trees:", n_null, "| Bootstrap per null:", n_boot_per_null, "\n")
  
  obs <- vector("list", n_boot_per_null)
  for (i in seq_len(n_boot_per_null)) {
    set.seed(52325 + i)
    obs[[i]] <- tryCatch(bootstrap_topology(meta, ps, tree, dissim_method), error = function(e) NULL)
  }
  mean_obs_mc <- mean(bind_rows(obs)$mc_raw, na.rm = TRUE)
  cat("  Observed mean Matching Split =", round(mean_obs_mc, 4), "\n")
  
  null_mc <- numeric(n_null)
  for (i in seq_len(n_null)) {
    set.seed(52325 + i)
    if (i %% 100 == 0) cat("    Null", i, "/", n_null, "\n")
    null_tree           <- rtree(n_spp, rooted = TRUE)
    null_tree$tip.label <- sample(unique(meta$host_taxon), n_spp, replace = FALSE)
    r <- tryCatch(bootstrap_topology(meta, ps, null_tree, dissim_method), error = function(e) NULL)
    null_mc[i] <- if (!is.null(r)) r$mc_raw[1] else NA
  }
  null_mc <- null_mc[!is.na(null_mc)]
  p_mc <- mean(null_mc <= mean_obs_mc)
  cat("  p (Matching Split):", round(p_mc, 4), "\n")
  
  list(mean_obs_mc = mean_obs_mc, null_mc = null_mc, p_mc = p_mc,
       n_null = length(null_mc), n_species = n_spp, dissim_method = dissim_method)
}

null_combos <- list(
  list(label = "16S_Amphibian_BrayCurtis", ps = ps_16S_wild,
       meta = meta_16S_wild %>% filter(Clade_Order %in% amphibian_orders),
       tree = host_tree_16S_full, dissim = "bray"),
  list(label = "16S_Amphibian_Jaccard", ps = ps_16S_wild,
       meta = meta_16S_wild %>% filter(Clade_Order %in% amphibian_orders),
       tree = host_tree_16S_full, dissim = "jaccard"),
  list(label = "16S_Reptile_BrayCurtis", ps = ps_16S_wild,
       meta = meta_16S_wild %>% filter(Clade_Order %in% reptile_orders),
       tree = host_tree_16S_full, dissim = "bray"),
  list(label = "16S_Reptile_Jaccard", ps = ps_16S_wild,
       meta = meta_16S_wild %>% filter(Clade_Order %in% reptile_orders),
       tree = host_tree_16S_full, dissim = "jaccard"),
  list(label = "ITS_Amphibian_BrayCurtis", ps = ps_ITS_wild,
       meta = meta_ITS_wild %>% filter(Clade_Order %in% amphibian_orders),
       tree = host_tree_ITS_full, dissim = "bray"),
  list(label = "ITS_Amphibian_Jaccard", ps = ps_ITS_wild,
       meta = meta_ITS_wild %>% filter(Clade_Order %in% amphibian_orders),
       tree = host_tree_ITS_full, dissim = "jaccard"),
  list(label = "ITS_Reptile_BrayCurtis", ps = ps_ITS_wild,
       meta = meta_ITS_wild %>% filter(Clade_Order %in% reptile_orders),
       tree = host_tree_ITS_full, dissim = "bray"),
  list(label = "ITS_Reptile_Jaccard", ps = ps_ITS_wild,
       meta = meta_ITS_wild %>% filter(Clade_Order %in% reptile_orders),
       tree = host_tree_ITS_full, dissim = "jaccard")
)

null_results_topology <- list()
for (combo in null_combos) {
  cat("\nNULL TEST:", combo$label, "\n")
  nr <- run_topology_null_test(
    ps = prune_samples(combo$meta$sample_name, combo$ps),
    meta = combo$meta, tree = combo$tree, dissim_method = combo$dissim,
    n_null = n_null, n_boot_per_null = n_boot_per_null
  )
  null_results_topology[[combo$label]] <- nr
  saveRDS(nr, file.path(dir_topology_null, paste0("topology_null_", combo$label, ".rds")))
}

################################################################################
# PART 4: SUPPLEMENTARY TABLE - reported output for this analysis
################################################################################
cat("\n=== PART 4: SUPPLEMENTARY TABLE ===\n")

supp_table <- bind_rows(lapply(names(results_topology), function(key) results_topology[[key]]$summary)) %>%
  mutate(
    Dataset = paste(dataset, paste0(clade, "s")),
    Distance_metric = factor(dissim, levels = c("BrayCurtis", "Jaccard"))
  ) %>%
  arrange(Dataset, Distance_metric) %>%
  transmute(
    Dataset, Distance_metric,
    Mean_MatchingSplit = round(mean_mc_raw, 4),
    SD_MatchingSplit    = round(sd_mc_raw, 4),
    N_species_mean      = round(mean_n_species, 4),
    N_species_SD        = round(sd_n_species, 4)
  )

null_p_table <- bind_rows(lapply(names(null_results_topology), function(k) {
  nr <- null_results_topology[[k]]
  parts <- strsplit(k, "_")[[1]]   # "{Marker}_{Clade}_{Dissim}"
  tibble(Dataset = paste(parts[1], paste0(parts[2], "s")),
         Distance_metric = parts[3],
         p_MatchingSplit_null = round(nr$p_mc, 4))
}))

supp_table <- left_join(supp_table, null_p_table, by = c("Dataset", "Distance_metric"))

supp_fname <- file.path(dir_topology, "topology_supplementary_table.csv")
write.csv(supp_table, supp_fname, row.names = FALSE)
cat("Saved supplementary table:", supp_fname, "\n")
print(as.data.frame(supp_table), max = 200)


################################################################################
################################################################################
################################################################################
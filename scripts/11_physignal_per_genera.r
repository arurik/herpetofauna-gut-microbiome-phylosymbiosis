################################################################################
# Pagel's Lambda Analysis: Phylogenetic Signal in Microbiome Composition
#
# Tests whether microbial lineage abundances show phylogenetic conservatism
# across the host phylogeny using Pagel's λ (range: 0 = no signal,
# 1 = Brownian motion expectation).
#
# Based on: Youngblut et al. 2019 (Nat Commun); Youngblut et al. 2021 (Nat Micro)
# Bootstrap approach: 100 iterations, 1 sample per species per iteration.
#
# Analyses run:
#   16S (Genus level):
#     - Amphibians (Caudata + Anura pooled)
#     - Reptiles   (Squamata + Testudines + Crocodilia pooled)
#     - Caudata    (order-level, for direct comparison)
#     - Anura      (order-level, for direct comparison)
#     - Squamata   (order-level, for direct comparison)
#   ITS (Genus level):
#     - Amphibians
#     - Reptiles
#     - Caudata    (feasibility-checked; run only if host-overlap criteria met)
#   ITS (Family level):
#     - Amphibians
#     - Reptiles
#     - Caudata    (feasibility-checked)
#
# Quality filters (applied within each bootstrap iteration):
#   - Minimum host prevalence: 15% of host species, with hard floor of 10 hosts
#   - Minimum iteration frequency: 50 of 100 bootstrap iterations (50%)
#
################################################################################

# ===== LOAD PACKAGES ===================================================== #
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(ape)
  library(phytools)
  library(phylosignal)
  library(phylobase)
  library(patchwork)
  library(ggtext)
})

# Mask common conflicts
filter  <- dplyr::filter
select  <- dplyr::select

# ===== PARAMETERS =================================================== 
set.seed(52325)

n_bootstrap       <- 100    # Bootstrap iterations
permutations      <- 999    # Permutations per phyloSignal call

MIN_HOST_PREVALENCE <- 0.15  # Taxon must occur in ≥15% of host species per iter
MIN_HOSTS_ABSOLUTE  <- 10    # Hard floor: require at least 10 host species (guards against
#   degenerate λ → 1 estimates from sparse, zero-inflated vectors)
MIN_ITERATIONS      <- 50    # Taxon must appear in ≥50 of 100 iterations

cat("  PAGEL'S LAMBDA ANALYSIS — PHYLOGENETIC SIGNAL IN MICROBIOMES\n")
cat("Bootstrap iterations:      ", n_bootstrap, "\n")
cat("Permutations per test:     ", permutations, "\n")
cat("Min host prevalence:       ", MIN_HOST_PREVALENCE * 100, "% (min", MIN_HOSTS_ABSOLUTE, "hosts absolute)\n")
cat("Min iteration frequency:   ", MIN_ITERATIONS, "of", n_bootstrap, "iterations\n\n")

# ===== COLOR PALETTE ================================================ 

# Hardcoded bacterial phylum colors — must match bact_colors_rel used in all
# other manuscript figures (host phylogeny, composition barplots, lollipop).
# NOTE: Synergistota and Planctomycetota share #6A3D9A; if both appear in the
# same figure panel, one will need to be swapped manually in Inkscape.
# Chloroflexi, Spirochaetota, and Deferribacterota are rare phyla not in the
# top-N composition figures; assigned distinct colorblind-friendly colors here.
bact_colors_rel <- c(
  "Actinobacteriota"  = "#A6CEE3",
  "Bacteroidota"      = "#1F78B4",
  "Campilobacterota"  = "#B2DF8A",
  "Chloroflexi"       = "#CAB2D6",
  "Deferribacterota"  = "#FFFF99",
  "Desulfobacterota"  = "#FB9A99",
  "Firmicutes"        = "#33A02C",
  "Fusobacteriota"    = "#E31A1C",
  "Planctomycetota"   = "#6A3D9A",
  "Proteobacteria"    = "#FDBF6F",
  "Spirochaetota"     = "#B15928",
  "Synergistota"      = "#6A3D9A",
  "Verrucomicrobiota" = "#FF7F00",
  "Other"             = "#222222"
)

# Fungal phylum colors for ITS figures.
its_colors_rel <- c(
  "Ascomycota"        = "#41AC66",
  "Basidiobolomycota" = "#66FFFF",
  "Basidiomycota"     = "#6366CC",
  "Mortierellomycota" = "#EFCCE5",
  "Mucoromycota"      = "#F266FF",
  "Rozellomycota"     = "#AAF266"
)

# Fallback for any phylum not listed above (should not occur in final dataset;
# a warning is printed if triggered so you can add the phylum explicitly).
get_phylum_color <- function(phylum_vec) {
  missing <- setdiff(phylum_vec, names(bact_colors_rel))
  if (length(missing) > 0) {
    warning("Phyla not in bact_colors_rel — assign manually: ",
            paste(missing, collapse = ", "))
  }
  bact_colors_rel[phylum_vec]
}

# ===== HELPER FUNCTIONS ============================================= #

# ===== PROJECT DIRECTORIES =================================================== #
library(here)
dir_processed <- here("data", "processed")
dir_trees     <- here("data", "raw", "trees")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")
output_base_dir <- here("output", "physignal_results")

# Create output directory: output_base_dir / phyloSignal_lambda_<tax_level>_level
create_output_dir <- function(base_dir, tax_level) {
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  tax_dir <- file.path(base_dir, paste0("phyloSignal_lambda_", tax_level, "_level"))
  if (!dir.exists(tax_dir)) dir.create(tax_dir)
  
  cat("Output directory:", tax_dir, "\n")
  return(tax_dir)
}

# Randomly draw 1 sample per host species (bootstrap unit)
bootstrap_one_sample_per_species <- function(ps, species_col = "host_taxon") {
  meta <- data.frame(sample_data(ps))
  
  sampled_ids <- meta %>%
    group_by(.data[[species_col]]) %>%
    slice_sample(n = 1) %>%
    pull(sample_name)
  
  ps_boot <- prune_samples(sampled_ids, ps)
  ps_boot <- prune_taxa(taxa_sums(ps_boot) > 0, ps_boot)
  return(ps_boot)
}

# Remap host_taxon values to match TimeTree tip labels
apply_replacements <- function(ps, replacements, dataset_name) {
  meta <- data.frame(sample_data(ps))
  meta$host_taxon <- as.character(meta$host_taxon)
  
  cat("\n--- Applying species name harmonization:", dataset_name, "---\n")
  changes <- 0
  for (old_name in names(replacements)) {
    new_name <- replacements[old_name]
    m <- meta$host_taxon == old_name
    if (sum(m, na.rm = TRUE) > 0) {
      cat("  ", old_name, "->", new_name,
          "(", sum(m, na.rm = TRUE), "samples)\n")
      meta$host_taxon[m] <- new_name
      changes <- changes + sum(m, na.rm = TRUE)
    }
  }
  sample_data(ps) <- meta
  cat("Total samples updated:", changes, "\n")
  return(ps)
}

# Retain only wild-collected samples (env_broad_scale == "Wild")
filter_to_wild <- function(ps, dataset_name) {
  meta <- data.frame(sample_data(ps))
  cat("\n--- Filtering to wild samples:", dataset_name, "---\n")
  cat("  Before:", nrow(meta), "samples\n")
  cat("  Wild:   ", sum(meta$env_broad_scale == "Wild", na.rm = TRUE), "\n")
  cat("  Captive:", sum(meta$env_broad_scale == "Captive", na.rm = TRUE), "\n")
  
  wild_ids <- meta$sample_name[meta$env_broad_scale == "Wild"]
  ps_wild  <- prune_samples(wild_ids, ps)
  cat("  After:", nsamples(ps_wild), "samples\n")
  return(ps_wild)
}

# ===== DATA LOADING ================================================= #
cat("===== LOADING DATA =====\n")

ps_16S <- readRDS(file.path(dir_processed, "16S_abs_final_fecal.rds"))
ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_fecal.rds"))

cat("16S samples loaded:", nsamples(ps_16S), "\n")
cat("ITS samples loaded:", nsamples(ps_ITS), "\n")

# Load host phylogeny
host_tree_wild <- read.tree(file.path(dir_trees, "tree_common_hosts_clean.nwk"))
host_tree_wild$tip.label <- gsub("_", " ", host_tree_wild$tip.label)
cat("Host tree loaded:", length(host_tree_wild$tip.label), "species\n\n")

# ===== SPECIES NAME HARMONIZATION =================================== #
cat("===== HARMONIZING SPECIES NAMES =====\n")

# Subspecies collapsed to species, or mapped to closest congener in TimeTree.
# All substitutions are phylogenetically informed (see Methods).
replacements_16S <- c(
  # Subspecies → species
  "Sceloporus occidentalis bocourtii" = "Sceloporus occidentalis",
  "Pituophis catenifer pumilus"        = "Pituophis catenifer",
  "Thamnophis elegans terrestris"      = "Thamnophis elegans",
  "Coluber constrictor mormon"         = "Coluber constrictor",
  "Thamnophis atratus atratus"         = "Thamnophis atratus",
  
  # Closest available congener in TimeTree
  "Agama picticauda"       = "Agama atra",
  "Terrapene carolina"     = "Terrapene ornata",
  "Desmognathus adatsihi"  = "Desmognathus ocoee",
  "Nerodia rhombifer"      = "Nerodia erythrogaster",
  "Thamnophis proximus"    = "Thamnophis sirtalis",
  "Anolis distichus"       = "Anolis cristatellus",
  
  # Species primarily present in captive samples
  "Atelopus balios"               = "Atelopus longirostris",
  "Dendrobates tinctorius azureus"= "Dendrobates tinctorius",
  "Ambystoma annulatum"           = "Ambystoma opacum",
  "Heloderma exasperatum"         = "Heloderma horridum",
  "Hyla cinerea"                  = "Dryophytes cinereus",
  "Hyla avivoca"                  = "Dryophytes avivoca",
  "Hyla chrysoscelis"             = "Dryophytes chrysoscelis",
  "Geochelone gigantea"           = "Aldabrachelys gigantea",
  "Litoria caerulea"              = "Ranoidea caerulea"
)

# ITS carries all 16S replacements plus ITS-specific subspecies
replacements_ITS <- c(
  replacements_16S,
  "Crotalus willardi silus"   = "Crotalus willardi",
  "Crotalus lepidus klauberi" = "Crotalus lepidus"
)

ps_16S <- apply_replacements(ps_16S, replacements_16S, "16S")
ps_ITS <- apply_replacements(ps_ITS, replacements_ITS, "ITS")

cat("\n16S unique species after harmonization:", length(unique(sample_data(ps_16S)$host_taxon)), "\n")
cat("ITS unique species after harmonization:", length(unique(sample_data(ps_ITS)$host_taxon)), "\n")

# ===== FILTER TO WILD SAMPLES ======================================= #
cat("===== FILTERING TO WILD SAMPLES =====\n")

ps_16S_wild <- filter_to_wild(ps_16S, "16S")
ps_ITS_wild <- filter_to_wild(ps_ITS, "ITS")

# ===== VALIDATE SPECIES COVERAGE ==================================== #
cat("===== VALIDATING SPECIES COVERAGE =====\n")

# Host order groupings used throughout
amphibian_orders <- c("Caudata", "Anura")
reptile_orders   <- c("Squamata", "Testudines", "Crocodilia")

meta_16S_wild <- data.frame(sample_data(ps_16S_wild))
meta_ITS_wild <- data.frame(sample_data(ps_ITS_wild))

species_16S_wild  <- unique(meta_16S_wild$host_taxon)
species_ITS_wild  <- unique(meta_ITS_wild$host_taxon)
species_tree_wild <- host_tree_wild$tip.label

cat("16S wild species:", length(species_16S_wild), "\n")
cat("ITS wild species:", length(species_ITS_wild), "\n")
cat("Tree species:    ", length(species_tree_wild), "\n\n")

mismatch_16S <- setdiff(species_16S_wild, species_tree_wild)
mismatch_ITS <- setdiff(species_ITS_wild, species_tree_wild)

if (length(mismatch_16S) > 0) {
  cat("WARNING: 16S species not in tree:\n"); print(mismatch_16S)
}
if (length(mismatch_ITS) > 0) {
  cat("WARNING: ITS species not in tree:\n"); print(mismatch_ITS)
}

if (length(mismatch_16S) > 0 | length(mismatch_ITS) > 0) {
  stop("STOPPING: Species in data not found in tree. ",
       "Update harmonization table or host phylogeny before proceeding.")
}

cat("✓ All species validated.\n\n")

# ===== PAGEL'S LAMBDA BOOTSTRAP FUNCTION =========================== #

# Runs 100 bootstrap iterations (1 sample per species per iteration), computes
# Pagel's lambda and permutation p-value for each qualifying taxon, and
# summarizes results across iterations with BH-corrected median p-values.
#
# Returns a list with:
#   $raw                — per-iteration results (all taxa)
#   $summary            — iteration-filtered summary (main results)
#   $summary_unfiltered — summary before iteration filter (for reference)
#   $excluded_taxa      — significant taxa excluded by iteration filter
#   $params             — analysis parameters

run_pagel_bootstrap <- function(
    ps,
    host_tree,
    tax_level   = "Genus",
    n_iter      = 100,
    species_col = "host_taxon",
    seed        = 52325,
    permutations = 999,
    output_dir  = NULL,
    analysis_name = ""
) {
  set.seed(seed)
  options(warn = -1)
  
  cat("\n  Running:", analysis_name, "| Tax level:", tax_level,
      "| Iter:", n_iter, "| Perm:", permutations, "\n")
  
  # Agglomerate to requested taxonomic level; NArm = FALSE retains unclassified
  ps_glom <- tax_glom(ps, taxrank = tax_level, NArm = FALSE)
  cat("  Taxa after agglomeration:", ntaxa(ps_glom), "\n\n")
  
  # Build taxonomy lookup for annotating results
  tax_tab_full <- as.data.frame(tax_table(ps_glom)) %>%
    rownames_to_column("OTU_ID") %>%
    mutate(
      Genus   = ifelse(is.na(Genus)   | Genus   == "", paste0("Unclassified_", Family), Genus),
      Family  = ifelse(is.na(Family)  | Family  == "", paste0("Unclassified_", Order),  Family),
      Order   = ifelse(is.na(Order)   | Order   == "", paste0("Unclassified_", Class),  Order),
      Class   = ifelse(is.na(Class)   | Class   == "", paste0("Unclassified_", Phylum), Class),
      Phylum  = ifelse(is.na(Phylum)  | Phylum  == "", paste0("Unclassified_", Kingdom), Phylum),
      Kingdom = ifelse(is.na(Kingdom) | Kingdom == "", "Unclassified_Kingdom", Kingdom)
    )
  
  # taxon_name = the focal rank column for labeling outputs
  tax_tab_full <- tax_tab_full %>%
    mutate(taxon_name = .data[[tax_level]])
  
  taxon_lookup <- tax_tab_full %>%
    select(OTU_ID, taxon_name, Kingdom, Phylum, Class, Order, Family, Genus)
  
  all_results <- vector("list", n_iter)
  start_time  <- Sys.time()
  
  for (i in seq_len(n_iter)) {
    if (i %% 10 == 0 || i == 1) {
      cat(sprintf("  Iteration %d/%d [%s]\n", i, n_iter,
                  format(Sys.time(), "%H:%M:%S")))
    }
    
    ps_boot  <- bootstrap_one_sample_per_species(ps_glom, species_col)
    meta     <- data.frame(sample_data(ps_boot))
    otu_mat  <- as(otu_table(ps_boot), "matrix")
    if (taxa_are_rows(ps_boot)) otu_mat <- t(otu_mat)
    
    species_abund           <- as.data.frame(otu_mat)
    rownames(species_abund) <- meta[[species_col]]
    
    shared_species <- intersect(rownames(species_abund), host_tree$tip.label)
    if (length(shared_species) < 3) next
    
    species_abund <- species_abund[shared_species, , drop = FALSE]
    tree_pruned   <- suppressMessages(
      drop.tip(host_tree, setdiff(host_tree$tip.label, shared_species))
    )
    
    # Reorder species_abund rows to match tree tip label order.
    # intersect() preserves the order of the first argument (species_abund
    # rownames), which may differ from tree_pruned$tip.label after drop.tip().
    # phylo4d requires rows to correspond to tips in tip order.
    species_abund <- species_abund[tree_pruned$tip.label, , drop = FALSE]
    
    # Prevalence filter: taxon must occur in ≥15% of hosts (min 10)
    host_prevalence     <- colSums(species_abund > 0)
    n_host_species      <- length(shared_species)
    min_hosts_required  <- max(MIN_HOSTS_ABSOLUTE,
                               ceiling(n_host_species * MIN_HOST_PREVALENCE))
    
    if (i == 1) {
      cat(sprintf("  Prevalence filter: %d host species, min %d hosts (%.0f%%)\n",
                  n_host_species, min_hosts_required,
                  (min_hosts_required / n_host_species) * 100))
      # Sanity check: after explicit reordering, rows must match tip labels.
      row_match <- all(rownames(species_abund) == tree_pruned$tip.label)
      if (!row_match) {
        stop("ALIGNMENT ERROR (iter 1): reordering failed — rownames still do not ",
             "match tree_pruned$tip.label. Check for duplicate or NA species names.")
      } else {
        cat("  ✓ Tip label alignment confirmed (iter 1)\n")
      }
    }
    
    keep_taxa <- host_prevalence >= min_hosts_required
    if (sum(keep_taxa) == 0) next
    
    species_abund <- species_abund[, keep_taxa, drop = FALSE]
    
    # suppressWarnings() silences the "Found more than one class 'phylo' in cache"
    # message from phylobase when phyloseq and RNeXML are both loaded.
    p4d <- tryCatch(
      suppressWarnings(
        phylobase::phylo4d(tree_pruned, tip.data = species_abund)
      ),
      error = function(e) NULL
    )
    if (is.null(p4d)) next
    
    physig_result <- tryCatch(
      suppressMessages(
        phylosignal::phyloSignal(p4d, methods = "Lambda", reps = permutations)
      ),
      error = function(e) NULL
    )
    if (is.null(physig_result)) next
    
    lambda_vals <- physig_result$stat$Lambda
    pvals       <- physig_result$pvalue$Lambda
    otu_ids     <- rownames(physig_result$stat)
    
    iter_results <- tibble(
      OTU_ID    = otu_ids,
      lambda    = lambda_vals,
      p_value   = pvals,
      iteration = i
    ) %>%
      left_join(taxon_lookup, by = "OTU_ID") %>%
      mutate(
        n_hosts        = host_prevalence[OTU_ID],
        mean_abundance = colMeans(species_abund[, OTU_ID, drop = FALSE])
      )
    
    all_results[[i]] <- iter_results
  }
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("\n  ✓ Done: %.1f minutes\n\n", elapsed))
  
  results_raw <- bind_rows(all_results)
  if (nrow(results_raw) == 0) stop("No results generated. Check input data.")
  
  # Summarize across iterations
  summary_all <- results_raw %>%
    group_by(OTU_ID, taxon_name, Kingdom, Phylum, Class, Order, Family, Genus) %>%
    summarise(
      n_iterations   = n(),
      mean_lambda    = mean(lambda,   na.rm = TRUE),
      median_lambda  = median(lambda, na.rm = TRUE),
      sd_lambda      = sd(lambda,     na.rm = TRUE),
      min_lambda     = min(lambda,    na.rm = TRUE),
      max_lambda     = max(lambda,    na.rm = TRUE),
      pct_sig        = mean(p_value < 0.05, na.rm = TRUE) * 100,
      median_p       = median(p_value, na.rm = TRUE),
      mean_n_hosts   = mean(n_hosts,  na.rm = TRUE),
      mean_abundance = mean(mean_abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      adj_p        = p.adjust(median_p, method = "BH"),
      sig_category = case_when(
        adj_p < 0.001 ~ "***",
        adj_p < 0.01  ~ "**",
        adj_p < 0.05  ~ "*",
        TRUE           ~ "ns"
      )
    ) %>%
    mutate(across(c(mean_lambda, median_lambda, sd_lambda,
                    min_lambda, max_lambda),         ~ round(.x, 3)),
           across(c(median_p, adj_p),               ~ round(.x, 4)),
           pct_sig        = round(pct_sig, 1),
           mean_n_hosts   = round(mean_n_hosts, 1),
           mean_abundance = round(mean_abundance, 1)) %>%
    arrange(adj_p, desc(mean_lambda))
  
  # Apply iteration frequency filter
  cat("  Before iteration filter:", nrow(summary_all), "taxa\n")
  summary_filtered <- summary_all %>% filter(n_iterations >= MIN_ITERATIONS)
  cat("  After iteration filter (≥", MIN_ITERATIONS, "iterations):",
      nrow(summary_filtered), "taxa\n")
  cat("  Significant (adj. p < 0.05):",
      sum(summary_filtered$adj_p < 0.05), "/", nrow(summary_filtered), "\n")
  
  excluded_taxa <- summary_all %>%
    filter(n_iterations < MIN_ITERATIONS, adj_p < 0.05) %>%
    arrange(n_iterations)
  
  if (nrow(excluded_taxa) > 0) {
    cat("  WARNING:", nrow(excluded_taxa),
        "significant taxa excluded by iteration filter:\n")
    print(excluded_taxa %>%
            select(taxon_name, n_iterations, mean_lambda, median_p, adj_p))
  } else {
    cat("  ✓ No significant taxa excluded by iteration filter.\n")
  }
  
  options(warn = 0)
  
  return(list(
    raw                = results_raw,
    summary            = summary_filtered,
    summary_unfiltered = summary_all,
    excluded_taxa      = excluded_taxa,
    params             = list(
      tax_level           = tax_level,
      n_iter              = n_iter,
      min_host_prevalence = MIN_HOST_PREVALENCE,
      min_hosts_absolute  = MIN_HOSTS_ABSOLUTE,
      min_iterations      = MIN_ITERATIONS,
      permutations        = permutations,
      seed                = seed
    )
  ))
}

# ===== FEASIBILITY CHECK FOR ORDER-LEVEL ITS ANALYSES ============= #
cat("===== FEASIBILITY CHECK — ORDER-LEVEL ITS =====\n")

# For ITS order-level analyses, check whether sufficient host-species overlap
# exists between data and tree. An order group is included only if ≥10 wild
# species are present in both the ITS dataset and the host tree — a conservative
# threshold to ensure the prevalence filter can operate meaningfully.

check_its_feasibility <- function(ps_wild, order_filter, host_tree,
                                  label, tax_level, min_species = 10) {
  meta_all   <- data.frame(sample_data(ps_wild))
  keep_ids   <- rownames(meta_all)[meta_all$Clade_Order %in% order_filter]
  ps_sub     <- prune_samples(keep_ids, ps_wild)
  meta   <- data.frame(sample_data(ps_sub))
  spp    <- unique(meta$host_taxon)
  overlap <- intersect(spp, host_tree$tip.label)
  
  cat(sprintf("  %s (%s): %d species in data, %d overlap with tree\n",
              label, tax_level, length(spp), length(overlap)))
  
  if (length(overlap) < min_species) {
    cat(sprintf("  → SKIPPED: fewer than %d species with tree overlap.\n\n",
                min_species))
    return(FALSE)
  } else {
    cat(sprintf("  → INCLUDED.\n\n"))
    return(TRUE)
  }
}

run_caudata_its_genus  <- check_its_feasibility(
  ps_ITS_wild, "Caudata",       host_tree_wild, "ITS Caudata Genus",  "Genus")
run_caudata_its_family <- check_its_feasibility(
  ps_ITS_wild, "Caudata",       host_tree_wild, "ITS Caudata Family", "Family")
run_anura_its_family   <- check_its_feasibility(
  ps_ITS_wild, "Anura",         host_tree_wild, "ITS Anura Family",   "Family")
run_reptile_its_family <- check_its_feasibility(
  ps_ITS_wild, reptile_orders,  host_tree_wild, "ITS Reptiles Family","Family")
run_amp_its_family     <- check_its_feasibility(
  ps_ITS_wild, amphibian_orders,host_tree_wild, "ITS Amphibians Family","Family")

# ===== DEFINE ANALYSES LIST ======================================== #
cat("===== DEFINING ANALYSES =====\n")

# --- 16S Genus-level analyses ---
# Core amphibian/reptile split (mirrors previous results for direct comparison)
# plus three focal orders for finer-scale comparison.
analyses_16S_genus <- list(
  list(name = "16S_Amphibians", tax_level = "Genus",
       ps = phyloseq::subset_samples(ps_16S_wild, Clade_Order %in% amphibian_orders)),
  list(name = "16S_Reptiles",   tax_level = "Genus",
       ps = phyloseq::subset_samples(ps_16S_wild, Clade_Order %in% reptile_orders)),
  list(name = "16S_Caudata",    tax_level = "Genus",
       ps = phyloseq::subset_samples(ps_16S_wild, Clade_Order == "Caudata")),
  list(name = "16S_Anura",      tax_level = "Genus",
       ps = phyloseq::subset_samples(ps_16S_wild, Clade_Order == "Anura")),
  list(name = "16S_Squamata",   tax_level = "Genus",
       ps = phyloseq::subset_samples(ps_16S_wild, Clade_Order == "Squamata"))
)

# --- ITS Genus-level analyses ---
# Amphibians and Reptiles run unconditionally; Caudata only if feasible.
analyses_ITS_genus <- list(
  list(name = "ITS_Amphibians", tax_level = "Genus",
       ps = phyloseq::subset_samples(ps_ITS_wild, Clade_Order %in% amphibian_orders)),
  list(name = "ITS_Reptiles",   tax_level = "Genus",
       ps = phyloseq::subset_samples(ps_ITS_wild, Clade_Order %in% reptile_orders))
)
if (run_caudata_its_genus) {
  analyses_ITS_genus <- c(analyses_ITS_genus, list(
    list(name = "ITS_Caudata_Genus", tax_level = "Genus",
         ps = phyloseq::subset_samples(ps_ITS_wild, Clade_Order == "Caudata"))
  ))
}

# --- ITS Family-level analyses ---
# Family level is added because fungal genera are often sparsely sampled;
# family-level agglomeration increases host prevalence and bootstrap stability.
analyses_ITS_family <- list()
if (run_amp_its_family) {
  analyses_ITS_family <- c(analyses_ITS_family, list(
    list(name = "ITS_Amphibians_Family", tax_level = "Family",
         ps = phyloseq::subset_samples(ps_ITS_wild, Clade_Order %in% amphibian_orders))
  ))
}
if (run_reptile_its_family) {
  analyses_ITS_family <- c(analyses_ITS_family, list(
    list(name = "ITS_Reptiles_Family", tax_level = "Family",
         ps = phyloseq::subset_samples(ps_ITS_wild, Clade_Order %in% reptile_orders))
  ))
}
if (run_caudata_its_family) {
  analyses_ITS_family <- c(analyses_ITS_family, list(
    list(name = "ITS_Caudata_Family", tax_level = "Family",
         ps = phyloseq::subset_samples(ps_ITS_wild, Clade_Order == "Caudata"))
  ))
}

# Combine all analyses into a single list for the main run loop
all_analyses <- c(analyses_16S_genus, analyses_ITS_genus, analyses_ITS_family)

cat("Total analyses to run:", length(all_analyses), "\n\n")

# ===== RUN ALL ANALYSES ============================================ #
cat("===== RUNNING PAGEL'S LAMBDA ANALYSES =====\n")

results_all <- list()

for (analysis in all_analyses) {
  
  cat("\n--- ANALYSIS:", analysis$name, "---\n")
  
  # Create tax-level-specific output directory
  output_dir <- create_output_dir(output_base_dir,
                                  tax_level = analysis$tax_level)
  
  result <- run_pagel_bootstrap(
    ps            = analysis$ps,
    host_tree     = host_tree_wild,
    tax_level     = analysis$tax_level,
    n_iter        = n_bootstrap,
    species_col   = "host_taxon",
    seed          = 52325,
    permutations  = permutations,
    output_dir    = output_dir,
    analysis_name = analysis$name
  )
  
  results_all[[analysis$name]] <- result
  
  # Save RDS for each analysis
  saveRDS(result, file.path(output_dir,
                            paste0(analysis$name, "_results.rds")))
  
  # The taxon_name column holds the focal rank label (Genus or Family).
  # Before renaming taxon_name -> rank_label, drop the existing taxonomy column
  # of the same name (e.g., the "Genus" column when tax_level == "Genus") to
  # prevent a duplicate-column error. The information is already in taxon_name.
  rank_label <- analysis$tax_level  # "Genus" or "Family"
  
  drop_and_rename <- function(df, rank) {
    df %>% select(-any_of(rank)) %>% rename(!!rank := taxon_name)
  }
  
  result_summary_out            <- drop_and_rename(result$summary,            rank_label)
  result_summary_unfiltered_out <- drop_and_rename(result$summary_unfiltered, rank_label)
  
  write_csv(result_summary_out,
            file.path(output_dir, paste0(analysis$name, "_summary.csv")))
  write_csv(result_summary_unfiltered_out,
            file.path(output_dir, paste0(analysis$name, "_summary_unfiltered.csv")))
  write_csv(result$raw,
            file.path(output_dir, paste0(analysis$name, "_raw_all_iterations.csv")))
  
  if (nrow(result$excluded_taxa) > 0) {
    excl_out <- drop_and_rename(result$excluded_taxa, rank_label)
    write_csv(excl_out,
              file.path(output_dir, paste0(analysis$name, "_excluded_taxa.csv")))
  }
  
  cat("\n  Top significant taxa:\n")
  print(result$summary %>%
          filter(adj_p < 0.05) %>%
          select(taxon_name, Phylum, n_iterations,
                 mean_lambda, median_p, adj_p, pct_sig) %>%
          head(10))
}

# Save combined RDS
saveRDS(results_all, file.path(output_base_dir,
                               "all_pagels_results_combined.rds"))

# (results_all is already in memory from the loop above — no need to reload it
# here. To resume from a saved run without re-running the bootstrap, use:
# results_all <- readRDS(file.path(output_base_dir, "all_pagels_results_combined.rds"))

# ===== SUMMARY CSV (SUPPLEMENTARY TABLE) =========================== #
cat("===== GENERATING SUMMARY CSV =====\n")

# Build one row per analysis for the supplementary overview table.
# taxon_name is renamed to the focal rank within each analysis block.
summary_rows <- lapply(names(results_all), function(nm) {
  res <- results_all[[nm]]
  rank_label <- res$params$tax_level
  results_all[[nm]]$summary %>%
    select(-any_of(rank_label)) %>%
    rename(!!rank_label := taxon_name) %>%
    mutate(Analysis = nm,
           Tax_level = rank_label,
           Dataset = ifelse(grepl("^16S", nm), "16S", "ITS"),
           Host_group = sub("^16S_|^ITS_", "", nm))
})

summary_combined <- bind_rows(summary_rows)

# High-level count table (for supp table header rows)
summary_table <- summary_combined %>%
  group_by(Analysis, Dataset, Host_group, Tax_level) %>%
  summarise(
    Total_taxa_tested = n(),
    Sig_taxa          = sum(adj_p < 0.05),
    Pct_sig           = sprintf("%.1f", 100 * Sig_taxa / Total_taxa_tested),
    Mean_lambda       = sprintf("%.3f", mean(mean_lambda, na.rm = TRUE)),
    Mean_iterations   = sprintf("%.1f", mean(n_iterations, na.rm = TRUE)),
    .groups = "drop"
  )

cat("\nOverall summary:\n")
print(summary_table)

# Output directory for genus and family level summaries (use genus-level dir
# as the home for the combined overview, since it's the primary analysis)
genus_output_dir <- create_output_dir(output_base_dir, "Genus")

write_csv(summary_table,   file.path(genus_output_dir, "overall_summary.csv"))
write_csv(summary_combined,file.path(genus_output_dir, "all_taxa_combined.csv"))

# Combine excluded taxa across all analyses
excluded_rows <- lapply(names(results_all), function(nm) {
  res <- results_all[[nm]]
  rank_label <- res$params$tax_level
  if (nrow(res$excluded_taxa) == 0) return(NULL)
  res$excluded_taxa %>%
    select(-any_of(rank_label)) %>%
    rename(!!rank_label := taxon_name) %>%
    mutate(Analysis = nm, Tax_level = rank_label)
})

excluded_combined <- bind_rows(excluded_rows)

if (!is.null(excluded_combined) && nrow(excluded_combined) > 0) {
  write_csv(excluded_combined,
            file.path(genus_output_dir, "all_excluded_taxa_combined.csv"))
  cat("\nExcluded taxa (significant but < MIN_ITERATIONS):\n")
  print(excluded_combined %>%
          group_by(Analysis) %>%
          summarise(n_excluded = n(), .groups = "drop"))
} else {
  cat("\n✓ No significant taxa excluded by iteration filter across any analysis.\n")
}

# ===== FIGURES ===================================================== #
cat("===== GENERATING FIGURES =====\n")

# Figures are generated for each analysis group separately.
# Point color = bacterial/fungal phylum; point size = mean number of host species.
# Panels are labeled in "A)" format per manuscript conventions.
# For ITS analyses, colors map to fungal phyla — see note below.
#
# ITS figures use its_colors_rel (defined in Part 3) for fungal phylum colors.

figures_dir <- file.path(output_base_dir, "figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Shared ggplot theme for all Pagel's lambda dot plots
theme_lambda <- function() {
  theme_classic() +
    theme(
      axis.text.y    = element_text(size = 8,  color = "black", lineheight = 0.9),
      axis.text.x    = element_text(size = 10, color = "black"),
      axis.title     = element_text(size = 11, face = "bold"),
      axis.line      = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      strip.text     = element_text(size = 10, face = "plain"),
      legend.position     = "right",
      legend.title        = element_text(size = 10, face = "bold"),
      legend.text         = element_text(size = 9),
      legend.key.size     = unit(1.2, "lines"),
      plot.margin    = margin(10, 10, 10, 10)
    )
}

# Prepare a data.frame for plotting from a single analysis result.
# tax_label_mode controls how the y-axis label is constructed:
#   "genus"  → "Order\nFamily\nGenus"
#   "family" → "Order\nFamily"
#
# null_summary: the null_summary_table produced in Part 13.5. Only taxa with
#   pct_null_below_obs >= null_pct_min are included in figures. This replaces
#   the former sd_lambda-based degenerate filter: rather than using a proxy
#   for estimation quality, we now use empirical null separation directly.
#   With a universal null median of 0 in these data, pct_null_below_obs >= 85
#   means the observed median lambda sits at or above the 85th percentile of
#   the null distribution for that taxon — a straightforward, defensible cut.
prepare_plot_data <- function(result, host_clade,
                              null_summary   = null_summary_table,
                              abundance_min  = 100,
                              tax_label_mode = "genus",
                              null_pct_min   = 85) {
  
  n_sig   <- sum(result$summary$adj_p < 0.05)
  
  # Join pct_null_below_obs onto observed summary by OTU_ID.
  # null_summary_table covers all analysis groups, so a distinct() ensures
  # no duplication when an OTU appears in multiple groups.
  null_info <- null_summary %>%
    select(OTU_ID, pct_null_below_obs) %>%
    distinct(OTU_ID, .keep_all = TRUE)
  
  n_filtered <- result$summary %>%
    filter(adj_p < 0.05) %>%
    left_join(null_info, by = "OTU_ID") %>%
    filter(is.na(pct_null_below_obs) | pct_null_below_obs < null_pct_min) %>%
    nrow()
  
  if (n_filtered > 0) {
    cat(sprintf("  [%s] Excluding %d taxa with pct_null_below_obs < %d%% from figure\n",
                host_clade, n_filtered, null_pct_min))
  }
  
  df <- result$summary %>%
    filter(adj_p < 0.05,
           ceiling(mean_abundance) >= abundance_min) %>%
    left_join(null_info, by = "OTU_ID") %>%
    filter(!is.na(pct_null_below_obs),
           pct_null_below_obs >= null_pct_min) %>%
    mutate(host_clade = host_clade)
  
  if (nrow(df) == 0) return(df)
  
  df <- df %>%
    arrange(desc(mean_lambda)) %>%
    mutate(rank = row_number(),
           mean_n_hosts_rounded = round(mean_n_hosts),
           Phylum = trimws(Phylum))
  
  if (tax_label_mode == "genus") {
    df <- df %>%
      mutate(tax_label = paste0(Order, "\n", Family, "\n", taxon_name))
  } else {
    df <- df %>%
      mutate(tax_label = paste0(Order, "\n", taxon_name))
  }
  
  df <- df %>%
    mutate(tax_label = fct_reorder(tax_label, rank))
  
  return(df)
}

# Single-panel dot plot for one analysis
create_lambda_panel <- function(data, panel_tag, panel_title,
                                color_palette = bact_colors_rel,
                                phylum_legend_title = "Bacterial Phylum",
                                x_limits = c(0.4, 1.05),
                                x_breaks = seq(0.4, 1, 0.2)) {
  
  if (nrow(data) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "No significant taxa", size = 5, color = "gray50") +
      labs(title = paste0(panel_tag, " ", panel_title)) +
      theme_void() +
      theme(plot.title = element_text(size = 11))
    return(p)
  }
  
  # Warn about any phyla not covered by the palette
  present_phyla <- unique(data$Phylum)
  missing       <- setdiff(present_phyla, names(color_palette))
  if (length(missing) > 0) {
    warning("Phyla missing from palette (assign in Inkscape): ",
            paste(missing, collapse = ", "))
    # Add placeholder grey so plot renders without error
    color_palette <- c(color_palette,
                       setNames(rep("#999999", length(missing)), missing))
  }
  
  data$Phylum <- factor(data$Phylum, levels = names(color_palette))
  
  ggplot(data, aes(x = mean_lambda, y = tax_label)) +
    geom_vline(xintercept = 1.0, linetype = "dashed",
               color = "gray40", linewidth = 0.5) +
    geom_point(aes(color = Phylum, size = mean_n_hosts_rounded), alpha = 0.8) +
    scale_color_manual(
      values = color_palette,
      name   = phylum_legend_title,
      guide  = guide_legend(override.aes = list(size = 5))
    ) +
    scale_size_continuous(
      name   = "Avg. host species\ntaxon detected in",
      range  = c(2, 6),
      breaks = function(x) pretty(x, n = 4)
    ) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks,
                       expand = c(0.02, 0.02)) +
    labs(x     = expression(paste("Pagel's ", lambda)),
         y     = NULL,
         title = paste0(panel_tag, " ", panel_title)) +
    theme_lambda()
}

# --- FIGURE: Amphibians vs Reptiles (16S Genus) ---
# Mirrors the previous Fig. 5C/5D layout for direct comparison.
amp_data <- prepare_plot_data(results_all[["16S_Amphibians"]], "Amphibians",
                              tax_label_mode = "genus")
rep_data <- prepare_plot_data(results_all[["16S_Reptiles"]],   "Reptiles",
                              tax_label_mode = "genus")

cat("16S Genus | Amphibians:", nrow(amp_data), "significant genera (≥100 reads)\n")
cat("16S Genus | Reptiles:  ", nrow(rep_data), "significant genera (≥100 reads)\n")

if (nrow(amp_data) > 0 | nrow(rep_data) > 0) {
  panel_A <- create_lambda_panel(amp_data, "A)", "Amphibians")
  panel_B <- create_lambda_panel(rep_data, "B)", "Reptiles")
  
  fig_amp_rep <- panel_A + panel_B +
    plot_layout(widths = c(1, 1.3), guides = "collect") &
    theme(legend.position = "right")
  
  cairo_pdf(
    filename = file.path(figures_dir,
                         "pagels_lambda_16S_AmpRep.pdf"),
    width  = 12,
    height = max(7, max(nrow(amp_data), nrow(rep_data)) * 0.35),
    bg     = "transparent"
  )
  print(fig_amp_rep)
  dev.off()
  cat("✓ Saved: pagels_lambda_16S_AmpRep.pdf\n")
}

# --- FIGURE: Order-level 16S (Caudata / Anura / Squamata) ---
order_analyses_16S <- c("16S_Caudata", "16S_Anura", "16S_Squamata")
order_labels       <- c("C) Caudata",  "D) Anura",  "E) Squamata")
panel_tags         <- c("C)", "D)", "E)")
panel_names        <- c("Caudata", "Anura", "Squamata")

order_data_list <- lapply(order_analyses_16S, function(nm) {
  if (!is.null(results_all[[nm]])) {
    prepare_plot_data(results_all[[nm]], sub("16S_", "", nm),
                      tax_label_mode = "genus")
  } else {
    data.frame()
  }
})

any_order_data <- any(sapply(order_data_list, nrow) > 0)
if (any_order_data) {
  panels_list <- mapply(create_lambda_panel,
                        data        = order_data_list,
                        panel_tag   = panel_tags,
                        panel_title = panel_names,
                        SIMPLIFY    = FALSE)
  
  fig_orders_16S <- wrap_plots(panels_list, nrow = 1,
                               guides = "collect") &
    theme(legend.position = "right")
  
  max_n <- max(sapply(order_data_list, nrow))
  cairo_pdf(
    filename = file.path(figures_dir,
                         "pagels_lambda_16S_orders.pdf"),
    width  = 18,
    height = max(7, max_n * 0.35),
    bg     = "transparent"
  )
  print(fig_orders_16S)
  dev.off()
  cat("✓ Saved: pagels_lambda_16S_orders.pdf\n")
} else {
  cat("No significant genera in order-level 16S analyses — skipping figure.\n")
}

# --- FIGURE: ITS Genus-level ---
its_genus_names <- grep("^ITS_.*(?<!Family)$", names(results_all),
                        perl = TRUE, value = TRUE)
its_genus_names <- its_genus_names[!grepl("Family", its_genus_names)]

if (length(its_genus_names) > 0) {
  its_genus_data <- lapply(seq_along(its_genus_names), function(i) {
    nm <- its_genus_names[i]
    prepare_plot_data(results_all[[nm]],
                      sub("^ITS_", "", nm),
                      tax_label_mode = "genus")
  })
  
  its_panel_tags   <- LETTERS[seq_along(its_genus_names)]
  its_panel_labels <- sub("^ITS_", "", its_genus_names)
  
  any_its_genus <- any(sapply(its_genus_data, nrow) > 0)
  if (any_its_genus) {
    its_panels <- mapply(create_lambda_panel,
                         data                = its_genus_data,
                         panel_tag           = paste0(its_panel_tags, ")"),
                         panel_title         = its_panel_labels,
                         color_palette       = list(its_colors_rel),
                         phylum_legend_title = "Fungal Phylum",
                         SIMPLIFY            = FALSE)
    
    fig_its_genus <- wrap_plots(its_panels, nrow = 1,
                                guides = "collect") &
      theme(legend.position = "right")
    
    max_n_its <- max(sapply(its_genus_data, nrow))
    cairo_pdf(
      filename = file.path(figures_dir,
                           "pagels_lambda_ITS_genus.pdf"),
      width  = 7 * length(its_genus_names),
      height = max(7, max_n_its * 0.35),
      bg     = "transparent"
    )
    print(fig_its_genus)
    dev.off()
    cat("✓ Saved: pagels_lambda_ITS_genus.pdf\n")
  } else {
    cat("No significant ITS genera — skipping genus-level ITS figure.\n")
  }
}

# --- FIGURE: ITS Family-level ---
its_family_names <- grep("Family", names(results_all), value = TRUE)
its_family_names <- its_family_names[grepl("^ITS", its_family_names)]

if (length(its_family_names) > 0) {
  its_fam_data <- lapply(seq_along(its_family_names), function(i) {
    nm <- its_family_names[i]
    prepare_plot_data(results_all[[nm]],
                      sub("^ITS_|_Family$", "", nm),
                      tax_label_mode = "family")
  })
  
  its_fam_tags   <- LETTERS[seq_along(its_family_names)]
  its_fam_labels <- paste0(sub("^ITS_|_Family$", "", its_family_names),
                           " (Family)")
  
  any_its_fam <- any(sapply(its_fam_data, nrow) > 0)
  if (any_its_fam) {
    its_fam_panels <- mapply(create_lambda_panel,
                             data                = its_fam_data,
                             panel_tag           = paste0(its_fam_tags, ")"),
                             panel_title         = its_fam_labels,
                             color_palette       = list(its_colors_rel),
                             phylum_legend_title = "Fungal Phylum",
                             SIMPLIFY            = FALSE)
    
    fig_its_fam <- wrap_plots(its_fam_panels, nrow = 1,
                              guides = "collect") &
      theme(legend.position = "right")
    
    max_n_fam <- max(sapply(its_fam_data, nrow))
    cairo_pdf(
      filename = file.path(figures_dir,
                           "pagels_lambda_ITS_family.pdf"),
      width  = 7 * length(its_family_names),
      height = max(7, max_n_fam * 0.35),
      bg     = "transparent"
    )
    print(fig_its_fam)
    dev.off()
    cat("✓ Saved: pagels_lambda_ITS_family.pdf\n")
  } else {
    cat("No significant ITS families — skipping family-level ITS figure.\n")
  }
}

# ===== FINAL CONSOLE SUMMARY ======================================= #
cat("===== FINAL SUMMARY =====\n\n")

cat("--- Observed results ---\n")
for (nm in names(results_all)) {
  res     <- results_all[[nm]]
  n_sig   <- sum(res$summary$adj_p < 0.05)
  n_total <- nrow(res$summary)
  n_excl  <- nrow(res$excluded_taxa)
  
  # Count how many significant taxa pass the null separation threshold
  n_null_pass <- if (!is.null(null_summary_table) && nrow(null_summary_table) > 0) {
    null_summary_table %>%
      filter(Analysis == nm, adj_p < 0.05, pct_null_below_obs >= 85) %>%
      nrow()
  } else NA_integer_
  
  cat(sprintf("  %-30s  %d / %d significant  (%d pass null filter, %d excluded)\n",
              nm, n_sig, n_total,
              ifelse(is.na(n_null_pass), 0, n_null_pass),
              n_excl))
}

cat("\n--- Null model results ---\n")
if (length(null_results) == 0) {
  cat("  No null model results (no analyses had significant taxa).\n")
} else {
  for (nm in names(null_results)) {
    obs_sig <- results_all[[nm]]$summary %>% filter(adj_p < 0.05) %>% pull(OTU_ID)
    null_df <- null_results[[nm]] %>% filter(OTU_ID %in% obs_sig)
    obs_df  <- results_all[[nm]]$raw %>% filter(OTU_ID %in% obs_sig)
    if (nrow(null_df) > 0 && nrow(obs_df) > 0) {
      cat(sprintf("  %-30s  obs median λ = %.3f  null median λ = %.3f\n",
                  nm,
                  median(obs_df$lambda,  na.rm = TRUE),
                  median(null_df$lambda, na.rm = TRUE)))
    }
  }
}

cat("\nAll outputs saved to:", output_base_dir, "\n")
cat("Note: Manually verify ITS phylum colors in Inkscape.\n")
cat("Note: paths are resolved via here() relative to the project root.\n\n")

# ===== FIGURE 4 — MAIN TEXT MULTIPANEL ============================= #
cat("===== GENERATING FIGURE 4 (MAIN TEXT) =====\n\n")

# Layout:
#   Left column  — Panel A: 16S Caudata (full height)
#   Right column — Panel B: 16S Squamata (top)
#                  Panel C: ITS combined Caudata + Reptiles (bottom)
#                  Panels B and C share the right column width and x-axis.
#
# Y-axis labels: three stacked lines — Order / Family / Genus — using plain
#   element_text(). Genus italics are applied manually in Inkscape after export.
#
# Legends:
#   Panels A + B share a bacterial phylum legend (bact_colors_rel).
#   Panel C has its own fungal phylum legend (its_colors_rel).
#   All three panels share a single circle size legend.
#   All legends are placed to the right of the figure.
#
# Panel C (ITS combined): Caudata and Reptile genera share a single y-axis,
#   separated by a horizontal divider line. Group labels ("Caudata" /
#   "Squamata") are placed at the right edge of the plot with annotate()
#   so they cannot overlap any data points.

# ---------------------------------------------------------------------------
# 16.1  Figure 4 theme
# ---------------------------------------------------------------------------

# theme_lambda_fig4 is kept as a self-contained alias of theme_lambda so
# Part 16 can be run independently. The two are currently identical; any
# Figure-4-specific tweaks (e.g. plot.margin adjustments for long labels)
# can be made here without affecting other figures in the script.
theme_lambda_fig4 <- function() {
  theme_classic() +
    theme(
      axis.text.y    = element_text(size = 8, color = "black", lineheight = 0.9),
      axis.text.x    = element_text(size = 10, color = "black"),
      axis.title     = element_text(size = 11, face = "bold"),
      axis.line      = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      strip.text     = element_text(size = 10, face = "plain"),
      legend.position     = "right",
      legend.title        = element_text(size = 10, face = "bold"),
      legend.text         = element_text(size = 9),
      legend.key.size     = unit(1.2, "lines"),
      plot.margin    = margin(10, 10, 10, 10)
    )
}

# ---------------------------------------------------------------------------
# 16.2  Y-axis label builder
# ---------------------------------------------------------------------------

# Produces stacked plain-text labels for the y-axis:
#   genus mode  → "Order\nFamily\nGenus"   (three lines)
#   family mode → "Order\nFamily"           (two lines)
# Genus/family italics are applied manually in Inkscape after PDF export,
# consistent with mycological convention for incertae sedis assignments.
make_stacked_label <- function(order, family, taxon, mode = "genus") {
  if (mode == "genus") {
    paste0(order, "\n", family, "\n", taxon)
  } else {
    paste0(order, "\n", taxon)
  }
}

# ---------------------------------------------------------------------------
# 16.3  prepare_plot_data_fig4 — wraps prepare_plot_data and rebuilds labels
#        using make_stacked_label for consistent three-line y-axis formatting
# ---------------------------------------------------------------------------

prepare_plot_data_fig4 <- function(result, host_clade,
                                   null_summary   = null_summary_table,
                                   abundance_min  = 100,
                                   tax_label_mode = "genus",
                                   null_pct_min   = 85) {
  
  df <- prepare_plot_data(
    result         = result,
    host_clade     = host_clade,
    null_summary   = null_summary,
    abundance_min  = abundance_min,
    tax_label_mode = tax_label_mode,
    null_pct_min   = null_pct_min
  )
  
  if (nrow(df) == 0) return(df)
  
  # Rebuild tax_label using the stacked plain-text format for Figure 4
  df <- df %>%
    mutate(tax_label = make_stacked_label(Order, Family, taxon_name,
                                          mode = tax_label_mode)) %>%
    mutate(tax_label = fct_reorder(tax_label, rank))
  
  return(df)
}

# ---------------------------------------------------------------------------
# 16.4  Unified size scale for all Figure 4 panels
# ---------------------------------------------------------------------------

# A shared scale_size_continuous object ensures circle sizes are visually
# comparable across panels A, B, and C regardless of the host count ranges
# in each dataset.  Limits are set to the observed data range across all
# three panels; breaks are chosen for a clean legend.
FIG4_SIZE_RANGE  <- c(2, 7)      # point radius range in mm
FIG4_SIZE_LIMITS <- c(5, 30)     # host count range across all panels
FIG4_SIZE_BREAKS <- c(5, 10, 15, 20, 25)

scale_size_fig4 <- function() {
  scale_size_continuous(
    name   = "Avg. host\nspecies",
    range  = FIG4_SIZE_RANGE,
    limits = FIG4_SIZE_LIMITS,
    breaks = FIG4_SIZE_BREAKS
  )
}

# ---------------------------------------------------------------------------
# 16.5  Shared x-axis limits across all panels
# ---------------------------------------------------------------------------

# All panels use the same x range so λ values are visually comparable.
# The ITS data have a slightly higher effective ceiling (λ ≈ 1.07) so
# we extend to 1.1 to avoid clipping.
FIG4_X_LIMITS <- c(0.55, 1.1)
FIG4_X_BREAKS <- c(0.6, 0.7, 0.8, 0.9, 1.0)

# ---------------------------------------------------------------------------
# 16.6  Panel builder for Figure 4
# ---------------------------------------------------------------------------

# Builds a single dot-plot panel using the Figure 4 theme and shared scales.
# show_size_legend / show_phylum_legend control which guide keys appear —
# used to suppress duplicate legends when panels are combined.
build_fig4_panel <- function(data, panel_tag, panel_title,
                             color_palette       = bact_colors_rel,
                             phylum_legend_title = "Bacterial phylum",
                             show_size_legend    = TRUE,
                             show_phylum_legend  = TRUE) {
  
  if (nrow(data) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.8, y = 0.5,
                 label = "No significant taxa", size = 4, color = "gray50") +
        labs(title = paste0(panel_tag, " ", panel_title)) +
        theme_void() +
        theme(plot.title = element_text(size = 11))
    )
  }
  
  # Ensure all phyla present in data are in palette; grey fallback for missing
  present_phyla <- unique(data$Phylum)
  missing_phyla <- setdiff(present_phyla, names(color_palette))
  if (length(missing_phyla) > 0) {
    warning("Phyla missing from palette (assign in Inkscape): ",
            paste(missing_phyla, collapse = ", "))
    color_palette <- c(color_palette,
                       setNames(rep("#999999", length(missing_phyla)),
                                missing_phyla))
  }
  
  data$Phylum <- factor(data$Phylum, levels = names(color_palette))
  
  size_guide   <- if (show_size_legend)   guide_legend(order = 2) else "none"
  phylum_guide <- if (show_phylum_legend) guide_legend(
    order = 1,
    override.aes = list(size = 4)
  ) else "none"
  
  ggplot(data, aes(x = mean_lambda, y = tax_label)) +
    geom_vline(xintercept = 1.0, linetype = "dashed",
               color = "gray40", linewidth = 0.5) +
    geom_point(aes(color = Phylum, size = mean_n_hosts_rounded), alpha = 0.85) +
    scale_color_manual(
      values = color_palette,
      name   = phylum_legend_title,
      guide  = phylum_guide
    ) +
    scale_size_fig4() +
    guides(size = size_guide) +
    scale_x_continuous(limits = FIG4_X_LIMITS, breaks = FIG4_X_BREAKS,
                       expand = c(0.02, 0.02)) +
    labs(
      x     = expression(paste("Pagel's ", lambda)),
      y     = NULL,
      title = paste0(panel_tag, " ", panel_title)
    ) +
    theme_lambda_fig4()
}

# ---------------------------------------------------------------------------
# 16.7  Prepare data for each panel
# ---------------------------------------------------------------------------

fig4_caudata_16S <- prepare_plot_data_fig4(
  results_all[["16S_Caudata"]], "Caudata"
)

fig4_squamata_16S <- prepare_plot_data_fig4(
  results_all[["16S_Squamata"]], "Squamata"
)

# ITS Caudata: genus level (Rozellomycota, Mortierellomycota, Umbelopsis)
fig4_caudata_ITS <- prepare_plot_data_fig4(
  results_all[["ITS_Caudata_Genus"]], "Caudata"
)

# ITS Reptiles: genus level. Nigrospora (80.3% null separation) is below
# the 85% threshold and is excluded by prepare_plot_data automatically.
# Only Basidiobolus (89%) will appear.
fig4_reptiles_ITS <- prepare_plot_data_fig4(
  results_all[["ITS_Reptiles"]], "Squamata"
)

cat("Figure 4 data:\n")
cat("  Panel A — 16S Caudata:  ", nrow(fig4_caudata_16S),  "genera\n")
cat("  Panel B — 16S Squamata: ", nrow(fig4_squamata_16S), "genera\n")
cat("  Panel C — ITS Caudata:  ", nrow(fig4_caudata_ITS),  "genera\n")
cat("  Panel C — ITS Squamata: ", nrow(fig4_reptiles_ITS), "genera\n\n")

# ---------------------------------------------------------------------------
# 16.8  Build Panel C: ITS combined with horizontal divider
# ---------------------------------------------------------------------------

# Combine Caudata and Reptile ITS genera into one data frame.
# A group column identifies which host order each taxon belongs to;
# the divider is drawn between the two groups using geom_hline at
# the boundary rank.  Group labels are added with annotate().
#
# Ordering: Caudata genera first (higher λ, primary result), then
# a gap/divider, then Reptile genera.

fig4_ITS_combined <- bind_rows(
  fig4_caudata_ITS  %>% mutate(host_group = "Caudata"),
  fig4_reptiles_ITS %>% mutate(host_group = "Squamata")
) %>%
  # Re-rank within the combined frame: Caudata first, Reptiles below,
  # each sub-group still ordered by descending mean_lambda
  arrange(host_group == "Squamata", desc(mean_lambda)) %>%
  mutate(combined_rank = row_number())

# Rebuild factor levels for the combined y-axis
fig4_ITS_combined <- fig4_ITS_combined %>%
  mutate(tax_label = fct_reorder(tax_label, combined_rank))

# Divider position: horizontal line between last Caudata and first Squamata
n_caudata_ITS   <- sum(fig4_ITS_combined$host_group == "Caudata")
n_squamata_ITS  <- sum(fig4_ITS_combined$host_group == "Squamata")
divider_y       <- n_squamata_ITS + 0.5   # in factor level coordinates

# Group label y positions (midpoint of each group's rows)
y_label_caudata  <- n_squamata_ITS + n_caudata_ITS / 2 + 0.5
y_label_squamata <- n_squamata_ITS / 2 + 0.5

# Build Panel C
if (nrow(fig4_ITS_combined) > 0) {
  
  # Build palette subset: only fungal phyla present in the combined ITS data
  present_its_phyla <- unique(fig4_ITS_combined$Phylum)
  its_palette_fig4  <- its_colors_rel[names(its_colors_rel) %in% present_its_phyla]
  missing_its <- setdiff(present_its_phyla, names(its_palette_fig4))
  if (length(missing_its) > 0) {
    its_palette_fig4 <- c(its_palette_fig4,
                          setNames(rep("#999999", length(missing_its)), missing_its))
  }
  
  fig4_ITS_combined$Phylum <- factor(fig4_ITS_combined$Phylum,
                                     levels = names(its_palette_fig4))
  
  panel_C <- ggplot(fig4_ITS_combined,
                    aes(x = mean_lambda, y = tax_label)) +
    geom_vline(xintercept = 1.0, linetype = "dashed",
               color = "gray40", linewidth = 0.5) +
    # Horizontal divider between Caudata and Squamata results
    geom_hline(yintercept = divider_y, color = "gray60",
               linewidth = 0.4, linetype = "solid") +
    geom_point(aes(color = Phylum, size = mean_n_hosts_rounded), alpha = 0.85) +
    # Group labels placed at the right edge of the plot, right-aligned,
    # so they cannot overlap with any data points (all points are left of x = 1.07).
    annotate("text", x = FIG4_X_LIMITS[2] - 0.002,
             y = y_label_caudata,
             label = "Caudata", hjust = 1, vjust = 0.5,
             size = 3, fontface = "italic", color = "gray35") +
    annotate("text", x = FIG4_X_LIMITS[2] - 0.002,
             y = y_label_squamata,
             label = "Squamata", hjust = 1, vjust = 0.5,
             size = 3, fontface = "italic", color = "gray35") +
    scale_color_manual(
      values = its_palette_fig4,
      name   = "Fungal phylum",
      guide  = guide_legend(order = 1, override.aes = list(size = 4))
    ) +
    scale_size_fig4() +
    guides(size = "none") +   # size legend collected from Panel A
    scale_x_continuous(limits = FIG4_X_LIMITS, breaks = FIG4_X_BREAKS,
                       expand = c(0.02, 0.02)) +
    labs(
      x     = expression(paste("Pagel's ", lambda)),
      y     = NULL,
      title = "C) Fungi"
    ) +
    theme_lambda_fig4()
  
} else {
  panel_C <- ggplot() +
    annotate("text", x = 0.8, y = 0.5,
             label = "No significant ITS taxa", size = 4, color = "gray50") +
    labs(title = "C) Fungi") +
    theme_void()
}

# ---------------------------------------------------------------------------
# 16.9  Build Panels A and B
# ---------------------------------------------------------------------------

panel_A <- build_fig4_panel(
  data                = fig4_caudata_16S,
  panel_tag           = "A)",
  panel_title         = "Bacteria — Caudata",
  color_palette       = bact_colors_rel,
  phylum_legend_title = "Bacterial phylum",
  show_size_legend    = TRUE,    # size legend lives here
  show_phylum_legend  = TRUE
)

panel_B <- build_fig4_panel(
  data                = fig4_squamata_16S,
  panel_tag           = "B)",
  panel_title         = "Bacteria — Squamata",
  color_palette       = bact_colors_rel,
  phylum_legend_title = "Bacterial phylum",
  show_size_legend    = FALSE,   # suppress: collected from Panel A
  show_phylum_legend  = TRUE
)

# ---------------------------------------------------------------------------
# 16.10  Assemble with patchwork
# ---------------------------------------------------------------------------

# Layout strategy:
#   - Left column:  Panel A alone, spans full height
#   - Right column: Panel B top, padding, Panel C bottom
#   - plot_layout(guides = "collect") merges legends per & operator scope
#
# The right-column stack is built first as its own patchwork unit, then
# combined with Panel A using a 2-column outer layout.
#
# Heights for right-column panels are proportional to the number of taxa
# so that point spacing is visually similar across panels.

n_A <- nrow(fig4_caudata_16S)
n_B <- nrow(fig4_squamata_16S)
n_C <- nrow(fig4_ITS_combined)

# Row height proportions for the right column
# Add 2 rows of padding between B and C via a plot_spacer()
spacer_rows    <- 2
right_heights  <- c(n_B, spacer_rows, n_C)

right_col <- (panel_B / plot_spacer() / panel_C) +
  plot_layout(heights = right_heights)

# Combine 16S legend across A+B panels and ITS legend for C separately.
# patchwork's guide collection works within a & call; we handle the two
# legend types by building separate patchwork units for 16S vs ITS,
# then joining them into the outer two-column layout.
#
# NOTE: guide collection for mixed-legend multipanels can behave
# unexpectedly — if the bacterial and fungal phylum legends merge
# incorrectly, separate them manually in Inkscape.

fig4 <- (panel_A | right_col) +
  plot_layout(
    widths  = c(1.6, 1),   # left column wider to accommodate more taxa
    guides  = "collect"
  ) &
  theme(legend.position = "right")

# ---------------------------------------------------------------------------
# 16.11  Save Figure 4
# ---------------------------------------------------------------------------

# Height: driven by the left column (Panel A, n_A genera).
# Allow 0.38 inches per genus plus base margins.
fig4_height <- max(8, n_A * 0.38 + 2)
fig4_width  <- 13   # two columns plus legend space

cairo_pdf(
  filename = file.path(figures_dir,
                       "Figure4_pagels_lambda_main.pdf"),
  width  = fig4_width,
  height = fig4_height,
  bg     = "transparent"
)
print(fig4)
dev.off()

cat("✓ Figure 4 saved:", file.path(dir_figures, "Figure4_pagels_lambda_main.pdf"), "\n")
cat("  Dimensions:", fig4_width, "×", sprintf("%.1f", fig4_height), "in\n")
cat("  Panels: A (", n_A, "genera), B (", n_B, "genera),",
    "C (", n_C, "total ITS genera)\n")
cat("  Note: if bacterial and fungal phylum legends merge, separate in Inkscape.\n\n")


###############################################################################
###############################################################################
###############################################################################
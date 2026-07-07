###############################################################################
# Alpha Diversity Analysis — Wild Samples Only
#
# Figures produced:
#   Fig S_a: 6-panel — Bacteria vs. Fungi overall + diet guild comparisons
#   Fig S_b: 2-panel — Alpha diversity by host family (Richness | Shannon),
#              Bact vs. Fungi within-family significance brackets
#   Fig S_c: 2-panel — Diet × Order (Richness | Shannon), no stats
#   Fig S_d: 4-panel — Diversity mapped onto host phylogeny (circular trees)
#              [to be finalized in Inkscape]
#   Fig S_e: 6-panel — Phylogenetic signal (all / amphibians / reptiles)
#              A–B: all species  C–D: Amphibia only  E–F: Reptilia only
#
# Alexander Rurik
###############################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(vegan)
  library(ape)
  library(phytools)
  library(ggpubr)
  library(patchwork)
  library(rstatix)
  library(ggtree)
  library(ggnewscale)
})

select <- dplyr::select
filter <- dplyr::filter

# ========================== BASIC PARAMETERS ==================================
set.seed(52325)

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_trees     <- here("data", "raw", "trees")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")

# ========================== LOAD DATA =========================================
ps_16S <- readRDS(file.path(dir_processed, "16S_abs_final_fecal.rds"))
ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_fecal.rds"))

# ========================== COLOR PALETTES ====================================
order_colors <- c(
  "Anura"      = "#CC79A7",
  "Caudata"    = "#E69F00",
  "Testudines" = "#009E73",
  "Squamata"   = "#18B1E9",
  "Crocodilia" = "#111189"
)

dataset_colors <- c(
  "Bacteria" = "#0072B2",
  "Fungi"    = "#E69F00"
)

amphibian_orders <- c("Anura", "Caudata")
reptile_orders   <- c("Squamata", "Testudines", "Crocodilia")

###############################################################################
# ========================== PART 1: CALCULATE ALPHA DIVERSITY ================
###############################################################################
cat("\n========== CALCULATING ALPHA DIVERSITY ==========\n")

calculate_diversity <- function(ps, dataset_name) {
  cat("\nCalculating diversity for", dataset_name, "...\n")
  
  otu <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) otu <- t(otu)
  
  metadata <- data.frame(sample_data(ps)) %>%
    select(-any_of(c("Observed", "Shannon", "Simpson", "InvSimpson", "Pielou")))
  
  diversity_df <- data.frame(
    sample_name = rownames(metadata),
    Observed    = colSums(otu > 0),
    Shannon     = diversity(t(otu), index = "shannon"),
    Simpson     = diversity(t(otu), index = "simpson"),
    InvSimpson  = diversity(t(otu), index = "invsimpson"),
    Pielou      = diversity(t(otu), index = "shannon") / log(colSums(otu > 0))
  )
  
  # Join metadata — handle both rowname and column sample_name formats
  if ("sample_name" %in% colnames(metadata)) {
    diversity_df <- left_join(diversity_df, metadata, by = "sample_name")
  } else {
    diversity_df <- left_join(diversity_df,
                              rownames_to_column(metadata, "sample_name"),
                              by = "sample_name")
  }
  
  diversity_df <- diversity_df %>%
    mutate(
      Clade = case_when(
        Clade_Order %in% amphibian_orders ~ "Amphibia",
        Clade_Order %in% reptile_orders   ~ "Reptilia",
        TRUE ~ NA_character_
      ),
      dataset = dataset_name
    )
  
  return(diversity_df)
}

div_16S <- calculate_diversity(ps_16S, "Bacteria")
div_ITS <- calculate_diversity(ps_ITS, "Fungi")

###############################################################################
# ========================== PART 2: FILTER TO WILD ===========================
###############################################################################
cat("\n========== FILTERING TO WILD SAMPLES ==========\n")

filter_to_wild <- function(div_df, label) {
  wild <- filter(div_df, env_broad_scale == "Wild")
  cat(label, "— total:", nrow(div_df), "| wild:", nrow(wild),
      "| species:", n_distinct(wild$host_taxon), "\n")
  return(wild)
}

div_16S_wild <- filter_to_wild(div_16S, "16S Bacteria")
div_ITS_wild <- filter_to_wild(div_ITS, "ITS Fungi")

cat("\n16S Wild Diversity Summary:\n")
print(summary(div_16S_wild[, c("Observed", "Shannon", "Simpson", "Pielou")]))
cat("\nITS Wild Diversity Summary:\n")
print(summary(div_ITS_wild[, c("Observed", "Shannon", "Simpson", "Pielou")]))

###############################################################################
# ========================== PART 3: STATISTICAL TESTS =======================
###############################################################################
cat("\n========== STATISTICAL TESTS (WILD ONLY) ==========\n")

# Kruskal-Wallis across host orders and diet guilds
run_kw_tests <- function(div_df, label) {
  cat("\n===", label, "===\n")
  for (grp in c("Clade_Order", "Diet")) {
    for (metric in c("Observed", "Shannon")) {
      kw <- kruskal.test(as.formula(paste(metric, "~", grp)), data = div_df)
      cat(metric, "~", grp, ": chi2 =", round(kw$statistic, 2),
          ", p =", format.pval(kw$p.value, digits = 3), "\n")
    }
  }
  # Wilcoxon: Amphibia vs Reptilia
  for (metric in c("Observed", "Shannon")) {
    wx <- wilcox.test(as.formula(paste(metric, "~ Clade")), data = div_df)
    cat(metric, "~ Clade (Amphibia vs Reptilia): W =", wx$statistic,
        ", p =", format.pval(wx$p.value, digits = 3), "\n")
  }
}

run_kw_tests(div_16S_wild, "16S Bacteria")
run_kw_tests(div_ITS_wild, "ITS Fungi")

#saveRDS(list(bacteria = div_16S_wild, fungi = div_ITS_wild), file.path(dir_processed, "alpha_diversity_WILD.rds"))

###############################################################################
# ========================== PART 4: SPECIES NAME HARMONIZATION ===============
###############################################################################
# Required before phylogenetic signal tests. Collapses subspecies and substitutes
# species absent from the TimeTree phylogeny with closest available congeners.
cat("\n========== HARMONIZING SPECIES NAMES ==========\n")

tree_16S <- read.tree(file.path(dir_trees, "16S_species_for_timetree_clean.nwk"))
tree_ITS <- read.tree(file.path(dir_trees, "ITS_species_for_timetree_clean.nwk"))
tree_16S$tip.label <- gsub("_", " ", tree_16S$tip.label)
tree_ITS$tip.label <- gsub("_", " ", tree_ITS$tip.label)
cat("16S tree tips:", length(tree_16S$tip.label),
    "| ITS tree tips:", length(tree_ITS$tip.label), "\n")

replacements_16S <- c(
  # Subspecies collapsed to species
  "Sceloporus occidentalis bocourtii" = "Sceloporus occidentalis",
  "Pituophis catenifer pumilus"        = "Pituophis catenifer",
  "Thamnophis elegans terrestris"      = "Thamnophis elegans",
  "Coluber constrictor mormon"         = "Coluber constrictor",
  "Thamnophis atratus atratus"         = "Thamnophis atratus",
  # Closest available congeners
  "Agama picticauda"      = "Agama atra",
  "Terrapene carolina"    = "Terrapene ornata",
  "Desmognathus adatsihi" = "Desmognathus ocoee",
  "Nerodia rhombifer"     = "Nerodia erythrogaster",
  "Thamnophis proximus"   = "Thamnophis sirtalis",
  "Anolis distichus"      = "Anolis cristatellus",
  # Species present only in captive samples
  "Atelopus balios"              = "Atelopus longirostris",
  "Dendrobates tinctorius azureus" = "Dendrobates tinctorius",
  "Ambystoma annulatum"          = "Ambystoma opacum",
  "Heloderma exasperatum"        = "Heloderma horridum",
  "Hyla cinerea"                 = "Dryophytes cinereus",
  "Hyla avivoca"                 = "Dryophytes avivoca",
  "Hyla chrysoscelis"            = "Dryophytes chrysoscelis",
  "Geochelone gigantea"          = "Aldabrachelys gigantea",
  "Litoria caerulea"             = "Ranoidea caerulea"
)

replacements_ITS <- c(
  replacements_16S,
  "Crotalus willardi silus"    = "Crotalus willardi",
  "Crotalus lepidus klauberi"  = "Crotalus lepidus"
)

harmonize_species <- function(div_df, replacements, label) {
  div_df$host_taxon <- as.character(div_df$host_taxon)
  n_changed <- 0
  for (old in names(replacements)) {
    idx <- div_df$host_taxon == old
    if (any(idx, na.rm = TRUE)) {
      cat("  ", old, "->", replacements[[old]], "(", sum(idx), "samples)\n")
      div_df$host_taxon[idx] <- replacements[[old]]
      n_changed <- n_changed + sum(idx)
    }
  }
  cat(label, ": total samples updated =", n_changed,
      "| unique species after =", n_distinct(div_df$host_taxon), "\n")
  return(div_df)
}

div_16S_wild <- harmonize_species(div_16S_wild, replacements_16S, "16S")
div_ITS_wild <- harmonize_species(div_ITS_wild, replacements_ITS, "ITS")

###############################################################################
# ========================== PART 5: PHYLOGENETIC SIGNAL ======================
###############################################################################
cat("\n========== PHYLOGENETIC SIGNAL (WILD ONLY) ==========\n")

# Returns a list with lambda + K results, pruned tree, species-mean div data,
# and the species list for use in downstream figures.
run_phylo_signal <- function(div_df, tree, label, clade_filter = NULL) {
  
  if (!is.null(clade_filter)) {
    div_df <- filter(div_df, Clade == clade_filter)
    label  <- paste0(label, " (", clade_filter, ")")
  }
  
  div_sp <- div_df %>%
    group_by(host_taxon) %>%
    summarise(Observed_mean = mean(Observed, na.rm = TRUE),
              Shannon_mean  = mean(Shannon,  na.rm = TRUE),
              n             = n(), .groups = "drop")
  
  shared <- intersect(tree$tip.label, div_sp$host_taxon)
  cat("\n===", label, "===\n")
  cat("Species in data:", nrow(div_sp), "| in tree:", length(tree$tip.label),
      "| in both:", length(shared), "\n")
  
  missing <- setdiff(div_sp$host_taxon, tree$tip.label)
  if (length(missing) > 0) {
    cat("WARNING: not in tree after harmonization:\n"); print(missing)
  }
  if (length(shared) < 10) {
    cat("Too few species for phylogenetic signal test (<10)\n"); return(NULL)
  }
  
  tree_p  <- keep.tip(tree, shared)
  div_p   <- filter(div_sp, host_taxon %in% shared) %>%
    arrange(match(host_taxon, tree_p$tip.label))
  
  obs_vec  <- setNames(div_p$Observed_mean, div_p$host_taxon)
  shan_vec <- setNames(div_p$Shannon_mean,  div_p$host_taxon)
  
  cat("Pagel's lambda...\n")
  lam_obs  <- phylosig(tree_p, obs_vec,  method = "lambda", test = TRUE)
  lam_shan <- phylosig(tree_p, shan_vec, method = "lambda", test = TRUE)
  
  cat("Blomberg's K (999 permutations)...\n")
  K_obs    <- phylosig(tree_p, obs_vec,  method = "K", test = TRUE, nsim = 999)
  K_shan   <- phylosig(tree_p, shan_vec, method = "K", test = TRUE, nsim = 999)
  
  cat("Observed — λ =", round(lam_obs$lambda, 3), "p =", format.pval(lam_obs$P, digits = 3),
      "| K =", round(K_obs$K, 3), "p =", format.pval(K_obs$P, digits = 3), "\n")
  cat("Shannon  — λ =", round(lam_shan$lambda, 3), "p =", format.pval(lam_shan$P, digits = 3),
      "| K =", round(K_shan$K, 3), "p =", format.pval(K_shan$P, digits = 3), "\n")
  
  list(lambda_obs = lam_obs, lambda_shan = lam_shan,
       K_obs = K_obs, K_shan = K_shan,
       n_species = length(shared), species_list = shared,
       div_species = div_p, tree_pruned = tree_p)
}

# All species
ps16_all  <- run_phylo_signal(div_16S_wild, tree_16S, "Bacteria")
psITS_all <- run_phylo_signal(div_ITS_wild, tree_ITS, "Fungi")

# Amphibia only
ps16_amph  <- run_phylo_signal(div_16S_wild, tree_16S, "Bacteria", "Amphibia")
psITS_amph <- run_phylo_signal(div_ITS_wild, tree_ITS, "Fungi",    "Amphibia")

# Reptilia only
ps16_rept  <- run_phylo_signal(div_16S_wild, tree_16S, "Bacteria", "Reptilia")
psITS_rept <- run_phylo_signal(div_ITS_wild, tree_ITS, "Fungi",    "Reptilia")

# Save all results
saveRDS(list(all_16S = ps16_all,  all_ITS = psITS_all,
             amph_16S = ps16_amph, amph_ITS = psITS_amph,
             rept_16S = ps16_rept, rept_ITS = psITS_rept),
        file.path(dir_processed, "phylo_signal_WILD.rds"))

# Summary CSV
extract_phylo_row <- function(res, label, metric) {
  if (is.null(res)) return(NULL)
  lam <- if (metric == "Observed") res$lambda_obs  else res$lambda_shan
  K   <- if (metric == "Observed") res$K_obs        else res$K_shan
  data.frame(
    Group = label, Metric = metric,
    Lambda = round(lam$lambda, 3), Lambda_p = lam$P,
    K      = round(K$K,        3), K_p      = K$P,
    N_species = res$n_species
  )
}

phylo_summary <- bind_rows(
  extract_phylo_row(ps16_all,  "All — Bacteria", "Observed"),
  extract_phylo_row(ps16_all,  "All — Bacteria", "Shannon"),
  extract_phylo_row(psITS_all, "All — Fungi",    "Observed"),
  extract_phylo_row(psITS_all, "All — Fungi",    "Shannon"),
  extract_phylo_row(ps16_amph,  "Amphibia — Bacteria", "Observed"),
  extract_phylo_row(ps16_amph,  "Amphibia — Bacteria", "Shannon"),
  extract_phylo_row(psITS_amph, "Amphibia — Fungi",    "Observed"),
  extract_phylo_row(psITS_amph, "Amphibia — Fungi",    "Shannon"),
  extract_phylo_row(ps16_rept,  "Reptilia — Bacteria", "Observed"),
  extract_phylo_row(ps16_rept,  "Reptilia — Bacteria", "Shannon"),
  extract_phylo_row(psITS_rept, "Reptilia — Fungi",    "Observed"),
  extract_phylo_row(psITS_rept, "Reptilia — Fungi",    "Shannon")
) %>%
  mutate(
    # BH-FDR correction applied separately within each statistic family
    # (12 lambda tests, 12 K tests) — they have different null distributions,
    # so they are not pooled into one correction.
    Lambda_p_BH = p.adjust(Lambda_p, method = "BH"),
    K_p_BH      = p.adjust(K_p,      method = "BH"),
    Lambda_sig  = case_when(Lambda_p_BH < 0.001 ~ "***", Lambda_p_BH < 0.01 ~ "**",
                            Lambda_p_BH < 0.05  ~ "*",   TRUE               ~ "ns"),
    K_sig       = case_when(K_p_BH < 0.001      ~ "***", K_p_BH < 0.01      ~ "**",
                            K_p_BH < 0.05       ~ "*",   TRUE               ~ "ns")
  )
write.csv(phylo_summary,
          file.path(dir_tables, "phylo_signal_summary_WILD.csv"), row.names = FALSE)
print(phylo_summary)

# Lookup BH-adjusted p-value from phylo_summary, for use in Fig S_e panels
get_bh_p <- function(group_label, dataset, test, metric) {
  grp_short <- recode(group_label,
                      "All species"   = "All",
                      "Amphibia only" = "Amphibia",
                      "Reptilia only" = "Reptilia")
  row <- phylo_summary %>%
    filter(Group == paste0(grp_short, " \u2014 ", dataset), Metric == metric)
  if (test == "Lambda") row$Lambda_p_BH else row$K_p_BH
}

###############################################################################
# ========================== PART 6: SAVE DIVERSITY DATA ======================
###############################################################################
write.csv(div_16S_wild, file.path(dir_tables, "16S_alpha_diversity_WILD.csv"), row.names = FALSE)
write.csv(div_ITS_wild, file.path(dir_tables, "ITS_alpha_diversity_WILD.csv"), row.names = FALSE)

# Summary tables
bind_rows(
  div_16S_wild %>% group_by(Clade_Order, Clade) %>%
    summarise(n = n(), n_spp = n_distinct(host_taxon),
              Observed_med = median(Observed), Shannon_med = median(Shannon),
              .groups = "drop") %>% mutate(dataset = "Bacteria"),
  div_ITS_wild %>% group_by(Clade_Order, Clade) %>%
    summarise(n = n(), n_spp = n_distinct(host_taxon),
              Observed_med = median(Observed), Shannon_med = median(Shannon),
              .groups = "drop") %>% mutate(dataset = "Fungi")
) %>% write.csv(file.path(dir_tables, "alpha_div_summary_by_order_WILD.csv"), row.names = FALSE)

cat("Data saved.\n")

###############################################################################
# ========================== FIGURES ==========================================
###############################################################################

##------------------------------------------------------------------------------
## Fig S_a: 6-panel — Bacteria vs. Fungi overall + diet guild comparisons ####
##------------------------------------------------------------------------------
# Panels A–B: overall Bacteria vs. Fungi richness and Shannon
# Panels C–D: richness by diet guild (Bacteria | Fungi)
# Panels E–F: Shannon by diet guild (Bacteria | Fungi)

# ---- Panel A: overall Bacteria vs. Fungi ----
div_bvf <- bind_rows(
  div_16S_wild %>% select(sample_name, Observed, Shannon) %>% mutate(Kingdom = "Bacteria"),
  div_ITS_wild %>% select(sample_name, Observed, Shannon) %>% mutate(Kingdom = "Fungi")
) %>% mutate(Kingdom = factor(Kingdom, levels = c("Bacteria", "Fungi")))

wilcox_obs_bvf  <- wilcox.test(Observed ~ Kingdom, data = div_bvf)
wilcox_shan_bvf <- wilcox.test(Shannon  ~ Kingdom, data = div_bvf)

make_bvf_stat <- function(wx, metric) {
  data.frame(
    group1    = "Bacteria", group2 = "Fungi",
    p         = wx$p.value,
    p.signif  = case_when(wx$p.value < 0.001 ~ "***", wx$p.value < 0.01 ~ "**",
                          wx$p.value < 0.05  ~ "*",   TRUE ~ "ns"),
    y.position = max(div_bvf[[metric]], na.rm = TRUE) * 1.08
  )
}
stat_obs_bvf  <- make_bvf_stat(wilcox_obs_bvf,  "Observed")
stat_shan_bvf <- make_bvf_stat(wilcox_shan_bvf, "Shannon")

n_bvf <- div_bvf %>%
  group_by(Kingdom) %>%
  summarise(n = n(), y_min_obs = min(Observed, na.rm = TRUE),
            y_min_shan = min(Shannon, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = paste0("n=", n))

pA_obs <- ggplot(div_bvf, aes(x = Kingdom, y = Observed, fill = Kingdom)) +
  geom_jitter(width = 0.2, alpha = 0.35, size = 1, color = "black") +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, color = "black", width = 0.5) +
  stat_pvalue_manual(stat_obs_bvf, label = "p.signif",
                     tip.length = 0.01, bracket.size = 0.5, size = 5) +
  geom_text(data = n_bvf, aes(x = Kingdom, y = y_min_obs, label = label),
            vjust = 1.8, size = 3, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = dataset_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = NULL, y = "Observed OTU Richness") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10))

pA_shan <- ggplot(div_bvf, aes(x = Kingdom, y = Shannon, fill = Kingdom)) +
  geom_jitter(width = 0.2, alpha = 0.35, size = 1, color = "black") +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, color = "black", width = 0.5) +
  stat_pvalue_manual(stat_shan_bvf, label = "p.signif",
                     tip.length = 0.01, bracket.size = 0.5, size = 5) +
  geom_text(data = n_bvf, aes(x = Kingdom, y = y_min_shan, label = label),
            vjust = 1.8, size = 3, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = dataset_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = NULL, y = "Shannon Diversity Index") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10))

# ---- Panels C–F: diet guild comparisons ----
diet_levels <- c("Herbivore", "Omnivore", "Insectivore", "Carnivore")

prepare_diet <- function(div_df, kingdom) {
  div_df %>%
    filter(!is.na(Diet), Diet %in% diet_levels) %>%
    filter(Diet %in% (count(., Diet) %>% filter(n >= 3) %>% pull(Diet))) %>%
    mutate(Diet = factor(Diet, levels = diet_levels), Kingdom = kingdom)
}

diet_16S <- prepare_diet(div_16S_wild, "Bacteria")
diet_ITS <- prepare_diet(div_ITS_wild, "Fungi")

get_diet_stats <- function(diet_df, metric) {
  kw <- kruskal.test(as.formula(paste(metric, "~ Diet")), data = diet_df)
  if (kw$p.value >= 0.05) return(NULL)
  diet_df %>%
    wilcox_test(as.formula(paste(metric, "~ Diet")), p.adjust.method = "fdr") %>%
    filter(p.adj < 0.05) %>%
    add_significance("p.adj") %>%
    rename(p.signif = p.adj.signif) %>%
    add_xy_position(x = "Diet", fun = "max", step.increase = 0.12)
}

make_n_diet <- function(diet_df, metric) {
  diet_df %>% group_by(Diet) %>%
    summarise(n = n(), y_pos = min(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    mutate(label = paste0("n=", n))
}

plot_diet_panel <- function(diet_df, metric, ylabel, kingdom, fill_color,
                            stat_df, n_df) {
  p <- ggplot(diet_df, aes(x = Diet, y = .data[[metric]])) +
    geom_jitter(width = 0.2, alpha = 0.35, size = 1, color = "black") +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, color = "black",
                 fill = fill_color, width = 0.55) +
    scale_y_continuous(expand = expansion(mult = c(0.14, 0.14))) +
    labs(x = NULL, y = ylabel, title = kingdom) +
    theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.x   = element_text(angle = 35, hjust = 1, size = 9),
          axis.title.y  = element_text(size = 10))
  if (!is.null(stat_df) && nrow(stat_df) > 0)
    p <- p + stat_pvalue_manual(stat_df, label = "p.signif",
                                tip.length = 0.01, bracket.size = 0.45,
                                size = 4.5, hide.ns = TRUE)
  p + geom_text(data = n_df, aes(x = Diet, y = y_pos, label = label),
                vjust = 2.0, size = 2.8, color = "black", inherit.aes = FALSE)
}

# Compute stats once and store — never call get_diet_stats inline
stats_16S_obs  <- get_diet_stats(diet_16S, "Observed")
stats_ITS_obs  <- get_diet_stats(diet_ITS, "Observed")
stats_16S_shan <- get_diet_stats(diet_16S, "Shannon")
stats_ITS_shan <- get_diet_stats(diet_ITS, "Shannon")

# Print to console so you can verify before plotting
cat("16S Observed diet stats:\n"); print(stats_16S_obs)
cat("ITS Observed diet stats:\n"); print(stats_ITS_obs)
cat("16S Shannon diet stats:\n");  print(stats_16S_shan)
cat("ITS Shannon diet stats:\n");  print(stats_ITS_shan)

pC_bact_obs  <- plot_diet_panel(diet_16S, "Observed", "Observed OTU Richness",
                                "Bacteria", dataset_colors["Bacteria"],
                                stats_16S_obs, make_n_diet(diet_16S, "Observed"))
pD_fung_obs  <- plot_diet_panel(diet_ITS, "Observed", "Observed OTU Richness",
                                "Fungi",    dataset_colors["Fungi"],
                                stats_ITS_obs, make_n_diet(diet_ITS, "Observed"))
pE_bact_shan <- plot_diet_panel(diet_16S, "Shannon", "Shannon Diversity Index",
                                "Bacteria", dataset_colors["Bacteria"],
                                stats_16S_shan, make_n_diet(diet_16S, "Shannon"))
pF_fung_shan <- plot_diet_panel(diet_ITS, "Shannon", "Shannon Diversity Index",
                                "Fungi",    dataset_colors["Fungi"],
                                stats_ITS_shan, make_n_diet(diet_ITS, "Shannon"))

fig_S_a <- (
  (pA_obs | pA_shan) /
    (pC_bact_obs | pD_fung_obs) /
    (pE_bact_shan | pF_fung_shan)
) +
  plot_annotation(
    tag_levels = list(c("A)", "B)", "C)", "D)", "E)", "F)")),
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  ) +
  plot_layout(heights = c(1, 1.1, 1.1))

ggsave(file.path(dir_figures, "FigS_alpha_bvf_diet_WILD.pdf"), fig_S_a, width = 10, height = 14, device = cairo_pdf, dpi = 400)

# Key stats for reporting
cat("\n--- Key stats ---\n")
cat("BvF Richness: W =", wilcox_obs_bvf$statistic,
    ", p =", format.pval(wilcox_obs_bvf$p.value, digits = 3),
    "| Bact median =", median(div_16S_wild$Observed, na.rm = TRUE),
    "| Fungi median =", median(div_ITS_wild$Observed, na.rm = TRUE), "\n")
cat("BvF Shannon:  W =", wilcox_shan_bvf$statistic,
    ", p =", format.pval(wilcox_shan_bvf$p.value, digits = 3), "\n")

##------------------------------------------------------------------------------
## Fig S_b: Alpha diversity by host family (Richness | Shannon) ####
## Bacteria vs. Fungi within-family significance brackets only
##------------------------------------------------------------------------------
cat("\n========== FIG S_b: DIVERSITY BY HOST FAMILY ==========\n")

build_family_panel <- function(div_16S, div_ITS, metric, ylabel, min_n = 3) {
  
  # Families present with >= min_n wild samples in BOTH datasets
  fam_16S <- div_16S %>% group_by(Family) %>%
    summarise(n = n(), Clade = first(Clade), .groups = "drop") %>% filter(n >= min_n)
  fam_ITS <- div_ITS %>% group_by(Family) %>%
    summarise(n = n(), .groups = "drop") %>% filter(n >= min_n)
  fam_keep <- intersect(fam_16S$Family, fam_ITS$Family)
  cat(ylabel, "— families in both datasets (n >=", min_n, "):", length(fam_keep), "\n")
  
  # Combine
  df <- bind_rows(
    div_16S %>% filter(Family %in% fam_keep) %>%
      select(sample_name, Family, Clade, all_of(metric)) %>%
      mutate(dataset = "Bacteria"),
    div_ITS %>% filter(Family %in% fam_keep) %>%
      select(sample_name, Family, Clade, all_of(metric)) %>%
      mutate(dataset = "Fungi")
  ) %>% mutate(dataset = factor(dataset, levels = c("Bacteria", "Fungi")))
  
  # Family order: Amphibia alphabetically, then Reptilia alphabetically
  fam_order <- df %>% select(Family, Clade) %>% distinct() %>%
    arrange(Clade, Family) %>% pull(Family)
  df$Family <- factor(df$Family, levels = fam_order)
  
  # Within-family Wilcoxon (Bacteria vs. Fungi), BH-corrected
  fam_stats <- df %>%
    group_by(Family) %>%
    filter(all(c("Bacteria", "Fungi") %in% dataset)) %>%
    group_modify(~ {
      wx <- wilcox.test(.x[[metric]] ~ .x$dataset)
      data.frame(p = wx$p.value,
                 y.position = max(.x[[metric]], na.rm = TRUE) * 1.12)
    }) %>%
    ungroup() %>%
    mutate(
      p.adj    = p.adjust(p, method = "BH"),
      p.signif = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
                           p.adj < 0.05  ~ "*",   TRUE          ~ "ns"),
      group1   = "Bacteria", group2 = "Fungi"
    ) %>%
    filter(p.adj < 0.05)
  
  # Plot
  # Sample sizes per family per dataset
  n_labels <- df %>%
    group_by(Family, dataset) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0("n=", n))
  
  ggplot(df, aes(x = dataset, y = .data[[metric]], fill = dataset)) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.7, color = "black") +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, color = "black", width = 0.6) +
    { if (nrow(fam_stats) > 0)
      stat_pvalue_manual(fam_stats, label = "p.signif",
                         tip.length = 0.01, bracket.size = 0.45, size = 3.5)
      else NULL } +
    geom_text(data = n_labels,
              aes(x = dataset, y = -Inf, label = label),
              vjust = -0.4, size = 2.5, color = "black", inherit.aes = FALSE) +
    facet_wrap(~ Family, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = dataset_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    labs(x = NULL, y = ylabel) +
    theme_classic() +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = 9, face = "plain"),   # <-- changed
      legend.position  = "bottom",
      legend.title     = element_blank(),
      axis.title.y     = element_text(size = 10)
    )
}

pFam_obs  <- build_family_panel(div_16S_wild, div_ITS_wild, "Observed",
                                "Observed OTU Richness")
pFam_shan <- build_family_panel(div_16S_wild, div_ITS_wild, "Shannon",
                                "Shannon Diversity Index")

fig_S_b <- (pFam_obs / pFam_shan) +
  plot_annotation(
    title      = paste0("Alpha diversity by host family (wild samples only; n \u2265 3)"),
    tag_levels = list(c("A)", "B)")),
    theme = theme(
      plot.title = element_text(face = "plain", size = 12, hjust = 0.5),
      plot.tag   = element_text(face = "bold",  size = 12)
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(dir_figures, "FigS_alpha_by_family_WILD.pdf"), fig_S_b, width = 12, height = 11, device = cairo_pdf, dpi = 400)

##------------------------------------------------------------------------------
## Fig S_c: Diet × Host Order (Richness | Shannon) ####
##------------------------------------------------------------------------------
cat("\n========== FIG S_c: DIET x ORDER ==========\n")

build_diet_order_panel <- function(metric, ylabel) {
  
  df <- bind_rows(
    div_16S_wild %>%
      filter(!is.na(Diet), !is.na(Clade_Order)) %>%
      select(Diet, Clade_Order, all_of(metric)) %>%
      mutate(dataset = "Bacteria"),
    div_ITS_wild %>%
      filter(!is.na(Diet), !is.na(Clade_Order)) %>%
      select(Diet, Clade_Order, all_of(metric)) %>%
      mutate(dataset = "Fungi")
  ) %>%
    mutate(Diet = factor(Diet, levels = diet_levels))
  
  # Keep only diet–order–dataset combos with >= 3 samples
  keep <- df %>%
    group_by(dataset, Diet, Clade_Order) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n >= 3)
  df <- semi_join(df, keep, by = c("dataset", "Diet", "Clade_Order"))
  
  dodge_width <- 0.85
  
  # Per-facet (dataset) y-range — used for bracket scaling
  facet_yrange <- df %>%
    group_by(dataset) %>%
    summarise(ymax   = max(.data[[metric]], na.rm = TRUE),
              yrange = diff(range(.data[[metric]], na.rm = TRUE)),
              .groups = "drop")
  
  stats_list <- list()
  
  for (ds in unique(df$dataset)) {
    
    y_info <- facet_yrange %>% filter(dataset == ds)
    y_step <- y_info$yrange * 0.10
    tick_h <- y_info$yrange * 0.02
    y_base <- y_info$ymax
    
    for (diet in levels(df$Diet)) {
      
      d <- df %>% filter(dataset == ds, Diet == diet)
      eligible_orders <- keep %>%
        filter(dataset == ds, Diet == diet) %>%
        pull(Clade_Order)
      
      if (length(eligible_orders) < 2) next
      d <- filter(d, Clade_Order %in% eligible_orders)
      
      if (length(eligible_orders) == 2) {
        # With 2 groups, wilcox_test returns 'p' not 'p.adj' — handle directly
        pw <- d %>%
          wilcox_test(as.formula(paste(metric, "~ Clade_Order"))) %>%
          filter(p < 0.05) %>%
          mutate(
            p.adj    = p,
            p.signif = case_when(
              p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"
            )
          )
      } else {
        # With 3+ groups, use KW as omnibus gate before post-hoc
        kw <- kruskal.test(as.formula(paste(metric, "~ Clade_Order")), data = d)
        if (kw$p.value >= 0.05) next
        pw <- d %>%
          wilcox_test(as.formula(paste(metric, "~ Clade_Order")),
                      p.adjust.method = "fdr") %>%
          filter(p.adj < 0.05) %>%
          add_significance("p.adj") %>%
          rename(p.signif = p.adj.signif)
      }
      
      if (nrow(pw) == 0) next
      
      # Dodge x-positions for orders within this Diet level
      diet_x    <- as.integer(factor(diet, levels = levels(df$Diet)))
      orders_in <- sort(unique(d$Clade_Order))
      n_orders  <- length(orders_in)
      offsets   <- seq(-dodge_width / 2 + dodge_width / (2 * n_orders),
                       dodge_width / 2 - dodge_width / (2 * n_orders),
                       length.out = n_orders)
      order_xpos <- setNames(diet_x + offsets, orders_in)
      
      # Stack brackets above the facet-specific data maximum
      for (i in seq_len(nrow(pw))) {
        g1    <- pw$group1[i]
        g2    <- pw$group2[i]
        y_top <- y_base + y_step * i
        
        stats_list[[length(stats_list) + 1]] <- data.frame(
          dataset    = ds,
          Diet       = diet,
          group1     = g1,
          group2     = g2,
          x          = order_xpos[g1],
          xend       = order_xpos[g2],
          y.position = y_top,
          tick_h     = tick_h,
          p.signif   = pw$p.signif[i],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # n labels
  n_labels <- keep %>% mutate(label = paste0("n=", n))
  
  p <- ggplot(df, aes(x = Diet, y = .data[[metric]], fill = Clade_Order)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width = dodge_width),
                alpha = 0.5, size = 0.8, color = "black") +
    geom_boxplot(position = position_dodge(dodge_width), outlier.shape = NA,
                 alpha = 0.7, width = 0.7, color = "black") +
    scale_fill_manual(values = order_colors, name = "Host Order") +
    facet_wrap(~ dataset, scales = "free_y") +
    geom_text(data = n_labels,
              aes(x = Diet, y = -Inf, label = label, group = Clade_Order),
              position = position_dodge(dodge_width), vjust = -0.4,
              size = 2.2, color = "black", inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.20))) +
    labs(x = "Diet", y = ylabel) +
    theme_classic() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold", size = 11),
      legend.position  = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1))
  
  # Draw brackets as geom_segment + geom_text, per facet
  if (length(stats_list) > 0) {
    brackets <- bind_rows(stats_list)
    
    p <- p +
      # Horizontal bracket line
      geom_segment(
        data = brackets,
        aes(x = x, xend = xend, y = y.position, yend = y.position),
        inherit.aes = FALSE, linewidth = 0.45, color = "black"
      ) +
      # Left tick
      geom_segment(
        data = brackets,
        aes(x = x, xend = x,
            y = y.position - tick_h, yend = y.position),
        inherit.aes = FALSE, linewidth = 0.45, color = "black"
      ) +
      # Right tick
      geom_segment(
        data = brackets,
        aes(x = xend, xend = xend,
            y = y.position - tick_h, yend = y.position),
        inherit.aes = FALSE, linewidth = 0.45, color = "black"
      ) +
      # Significance label
      geom_text(
        data = brackets,
        aes(x = (x + xend) / 2, y = y.position, label = p.signif),
        inherit.aes = FALSE, vjust = -0.3, size = 4.5
      )
  }
  
  p
}

pDO_obs  <- build_diet_order_panel("Observed", "Observed OTU Richness")
pDO_shan <- build_diet_order_panel("Shannon",  "Shannon Diversity Index")

fig_S_c <- (pDO_obs / pDO_shan) +
  plot_annotation(
    tag_levels = list(c("A)", "B)")),
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(dir_figures, "FigS_alpha_diet_order_WILD.pdf"), fig_S_c, width = 10, height = 10, device = cairo_pdf, dpi = 400)

##------------------------------------------------------------------------------
## Fig S_d: Diversity mapped onto host phylogeny (4-panel circular trees) ####
## Finalized in Inkscape; code kept for full reproducibility
##------------------------------------------------------------------------------

if (!is.null(ps16_all) && !is.null(psITS_all)) {
  
  cat("\n========== FIG S_d: PHYLOGENETIC TREE FIGURES ==========\n")
  
  plot_phylo_tree <- function(tree, div_sp, metric, label,
                              color_option = "viridis") {
    div_tree <- div_sp %>%
      select(host_taxon, !!sym(paste0(metric, "_mean"))) %>%
      rename(label = host_taxon, diversity = 2)
    tree_p <- keep.tip(tree, div_tree$label)
    cat(" ", label, metric, ":", length(tree_p$tip.label), "species\n")
    ggtree(tree_p, layout = "circular", size = 0.6) %<+% div_tree +
      geom_tippoint(aes(color = diversity), size = 4, alpha = 0.9) +
      scale_color_viridis_c(option = color_option, name = metric, direction = 1) +
      theme(legend.position  = "right",
            legend.title     = element_text(size = 14, face = "bold"),
            legend.text      = element_text(size = 12),
            legend.key.size  = unit(1.0, "cm"),
            plot.title       = element_text(hjust = 0.5, face = "bold", size = 16)) +
      labs(title = paste(label, "-", metric))
  }
  
  pT_16S_obs  <- plot_phylo_tree(tree_16S, ps16_all$div_species,
                                 "Observed", "Bacteria", "viridis")
  pT_16S_shan <- plot_phylo_tree(tree_16S, ps16_all$div_species,
                                 "Shannon",  "Bacteria", "plasma")
  pT_ITS_obs  <- plot_phylo_tree(tree_ITS, psITS_all$div_species,
                                 "Observed", "Fungi",    "mako")
  pT_ITS_shan <- plot_phylo_tree(tree_ITS, psITS_all$div_species,
                                 "Shannon",  "Fungi",    "rocket")
  
  fig_S_d <- (pT_16S_obs | pT_16S_shan) / (pT_ITS_obs | pT_ITS_shan) +
    plot_annotation(
      tag_levels = list(c("A)", "B)", "C)", "D)")),
      theme = theme(plot.tag  = element_text(size = 18, face = "bold"))
    )
  
  ggsave(file.path(dir_figures, "FigS_phylo_tree_diversity_WILD.pdf"),
         fig_S_d, width = 14, height = 14, device = cairo_pdf, dpi = 400)
  cat("Saved: Fig S_d (finalize background shading and dotted line in Inkscape)\n")
}

##------------------------------------------------------------------------------
## Fig S_e: Phylogenetic signal bar plot — 6 panels ####
## A–B: all species | C–D: Amphibia only | E–F: Reptilia only
## Each panel pair: Observed Richness (left) | Shannon (right)
##------------------------------------------------------------------------------
cat("\n========== FIG S_e: PHYLOGENETIC SIGNAL FIGURE ==========\n")

# panel_title: displayed in the grey strip (e.g. "All species — Bacteria (n = 87)")
# metric:      "Observed" or "Shannon"
# ylabel:      y-axis label
build_phylo_signal_panel <- function(res_bact, res_fung, metric, ylabel,
                                     group_label) {
  
  if (is.null(res_bact) || is.null(res_fung)) {
    message("Skipping panel — insufficient data"); return(NULL)
  }
  
  lam_b <- if (metric == "Observed") res_bact$lambda_obs  else res_bact$lambda_shan
  K_b   <- if (metric == "Observed") res_bact$K_obs        else res_bact$K_shan
  lam_f <- if (metric == "Observed") res_fung$lambda_obs   else res_fung$lambda_shan
  K_f   <- if (metric == "Observed") res_fung$K_obs         else res_fung$K_shan
  
  sig_label <- function(p)
    case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
  
  df <- data.frame(
    Dataset = rep(c("Bacteria", "Fungi"), each = 2),
    Test    = rep(c("Pagel's \u03bb", "Blomberg's K"), 2),
    Value   = c(lam_b$lambda, K_b$K, lam_f$lambda, K_f$K),
    P_value = c(
      get_bh_p(group_label, "Bacteria", "Lambda", metric),
      get_bh_p(group_label, "Bacteria", "K",      metric),
      get_bh_p(group_label, "Fungi",    "Lambda", metric),
      get_bh_p(group_label, "Fungi",    "K",      metric)
    ),    # n per dataset embedded in facet label
    Facet   = c(
      paste0("Bacteria (n = ", res_bact$n_species, ")"),
      paste0("Bacteria (n = ", res_bact$n_species, ")"),
      paste0("Fungi (n = ",    res_fung$n_species, ")"),
      paste0("Fungi (n = ",    res_fung$n_species, ")")
    )
  ) %>%
    mutate(
      Significance = sig_label(P_value),
      Dataset      = factor(Dataset, levels = c("Bacteria", "Fungi")),
      Facet        = factor(Facet,   levels = unique(Facet)),
      label_y      = Value + max(Value) * 0.07
    )
  
  ggplot(df, aes(x = Test, y = Value, fill = Dataset)) +
    geom_col(position = position_dodge(0.75), alpha = 0.85, width = 0.65,
             color = "black", linewidth = 0.3) +
    geom_text(aes(y = label_y, label = Significance),
              position = position_dodge(0.75), vjust = 0, size = 5) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "grey40", linewidth = 0.6) +
    facet_wrap(~ Facet, nrow = 1) +
    scale_fill_manual(values = dataset_colors, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(x = NULL, y = ylabel,
         title = group_label) +
    theme_classic() +
    theme(
      legend.position  = "bottom",
      axis.text.x      = element_text(size = 10),
      axis.title.y     = element_text(size = 10),
      plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = 9, face = "plain")
    )
}

# Build all 6 panels — group_label goes in the panel title
pSig_all_obs   <- build_phylo_signal_panel(ps16_all,  psITS_all,
                                           "Observed", "Observed OTU Richness",
                                           "All species")
pSig_all_shan  <- build_phylo_signal_panel(ps16_all,  psITS_all,
                                           "Shannon",  "Shannon Diversity Index",
                                           "All species")
pSig_amph_obs  <- build_phylo_signal_panel(ps16_amph, psITS_amph,
                                           "Observed", "Observed OTU Richness",
                                           "Amphibia only")
pSig_amph_shan <- build_phylo_signal_panel(ps16_amph, psITS_amph,
                                           "Shannon",  "Shannon Diversity Index",
                                           "Amphibia only")
pSig_rept_obs  <- build_phylo_signal_panel(ps16_rept, psITS_rept,
                                           "Observed", "Observed OTU Richness",
                                           "Reptilia only")
pSig_rept_shan <- build_phylo_signal_panel(ps16_rept, psITS_rept,
                                           "Shannon",  "Shannon Diversity Index",
                                           "Reptilia only")

# n per group for subtitle
n_bact_all  <- ps16_all$n_species
n_fung_all  <- psITS_all$n_species
n_bact_amph <- if (!is.null(ps16_amph))  ps16_amph$n_species  else 0
n_fung_amph <- if (!is.null(psITS_amph)) psITS_amph$n_species else 0
n_bact_rept <- if (!is.null(ps16_rept))  ps16_rept$n_species  else 0
n_fung_rept <- if (!is.null(psITS_rept)) psITS_rept$n_species else 0

fig_S_e <- (pSig_all_obs   | pSig_all_shan) /
  (pSig_amph_obs  | pSig_amph_shan) /
  (pSig_rept_obs  | pSig_rept_shan) +
  plot_annotation(
    title      = "Phylogenetic signal in alpha diversity (wild samples only)",
    tag_levels = list(c("A)", "B)", "C)", "D)", "E)", "F)")),
    theme = theme(
      plot.title = element_text(face = "plain", size = 12, hjust = 0.5),
      plot.tag   = element_text(face = "bold",  size = 12)
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(dir_figures, "FigS_phylo_signal_WILD.pdf"), fig_S_e, width = 10, height = 11, device = cairo_pdf, dpi = 400)

cat("\n========== ALL FIGURES COMPLETE ==========\n")


###############################################################################
###############################################################################
###############################################################################
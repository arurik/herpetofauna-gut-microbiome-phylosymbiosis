###############################################################################
# Sample type analysis: beta-dispersion, PERMANOVA, taxonomy barplots,
# phyloseq splitting.
#
# Workflow:
#   1.  Beta-dispersion (PERMDISP) — full dataset, Bray-Curtis & Jaccard
#   2.  Subsampled beta-dispersion — balanced group sizes to check size bias
#   3.  PERMANOVA (adonis2) — full dataset
#   4.  Sensitivity analysis — Fecal vs. Cloacal Swab restricted to Squamata
#   5.  Taxonomy barplots — species with paired fecal + cloacal swab samples
#   6.  Split phyloseq objects by sample type and save RDS
#   7.  ITS phyloseq splitting and save RDS
#
# Alexander Rurik
###############################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(vegan)
  library(ape)
  library(reshape2)
  library(RColorBrewer)
  library(patchwork)
  library(kableExtra)
  library(webshot2)
  library(magrittr)
  library(scales)
})

select <- dplyr::select
filter <- dplyr::filter

# ========================== HELPER FUNCTIONS ===================================

# Pairwise PERMANOVA across levels of a grouping factor
#
# Runs adonis2 on every pairwise subset of group_col levels (rather than
# relying on the omnibus multi-level test), with BH correction across the
# resulting family of pairwise p-values. This is the relevant test for
# decisions about collapsing specific factor levels (e.g., "are fecal and
# lower-GI similar enough to combine?"), which an omnibus 3+ level PERMANOVA
# cannot address on its own.
#
# dist_obj   a dist object (e.g., from phyloseq::distance)
# metadata   data.frame with rownames matching labels(dist_obj)
# group_col  character; column name in metadata giving group levels
# strata_col character or NULL; column name in metadata for strata
# permutations integer; number of permutations passed to adonis2
# returns - tibble with one row per pairwise comparison, BH-adjusted p-values

pairwise_permanova <- function(dist_obj, metadata, group_col,
                               strata_col = NULL, permutations = 999) {
  
  ids      <- labels(dist_obj)
  meta_sub <- metadata[ids, , drop = FALSE]   # enforce same order as dist_obj
  dist_mat <- as.matrix(dist_obj)
  
  levs  <- sort(unique(as.character(meta_sub[[group_col]])))
  pairs <- combn(levs, 2, simplify = FALSE)
  
  results <- purrr::map_dfr(pairs, function(pr) {
    
    keep_ids <- ids[meta_sub[[group_col]] %in% pr]
    sub_dist <- as.dist(dist_mat[keep_ids, keep_ids])
    sub_meta <- meta_sub[keep_ids, , drop = FALSE]
    sub_meta[[group_col]] <- factor(sub_meta[[group_col]], levels = pr)
    
    form <- as.formula(paste0("sub_dist ~ ", group_col))
    
    mod <- if (!is.null(strata_col)) {
      adonis2(form, data = sub_meta, permutations = permutations,
              strata = sub_meta[[strata_col]])
    } else {
      adonis2(form, data = sub_meta, permutations = permutations)
    }
    
    tibble(
      group1   = pr[1],
      group2   = pr[2],
      n1       = sum(sub_meta[[group_col]] == pr[1]),
      n2       = sum(sub_meta[[group_col]] == pr[2]),
      Df       = mod$Df[1],
      SumOfSqs = mod$SumOfSqs[1],
      R2       = mod$R2[1],
      F        = mod$F[1],
      p_raw    = mod$`Pr(>F)`[1]
    )
  })
  
  results$p_adj <- p.adjust(results$p_raw, method = "BH")
  results$sig   <- symnum(results$p_adj, corr = FALSE, na = FALSE,
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                          symbols   = c("***", "**", "*", "ns"))
  results
}

# ========================== BASIC PARAMETERS ==================================
set.seed(52325)

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")

# ========================== SHARED COLOR PALETTES =============================
clade_colors <- c(
  "Anura"      = "#CC79A7",
  "Caudata"    = "#E69F00",
  "Testudines" = "#009E73",
  "Squamata"   = "#18B1E9",
  "Crocodilia" = "#111189"
)

sampletype_colors <- c(
  "Cloacal swab" = "#1b9e77",
  "Lower GI"     = "#d95f02",
  "Fecal"        = "#7570b3"
)

###############################################################################
# ========================== 16S DATASET ======================================
###############################################################################

# ========================== LOAD DATA =========================================
ps_16S <- readRDS(file.path(dir_processed, "16S_absolute_all_final.rds"))

# ========================== DISTANCE MATRICES & METADATA (16S) ===============
full_ds_bray_16S <- phyloseq::distance(ps_16S, method = "bray")
# binary = TRUE: true presence/absence Jaccard. Without this, vegdist's default
# (binary = FALSE) computes the abundance-weighted quantitative Jaccard
# (2*Bray/(1+Bray) — a monotonic transform of Bray-Curtis), not presence/absence.
full_ds_jacc_16S <- phyloseq::distance(ps_16S, method = "jaccard", binary = TRUE)
full_meta_16S    <- data.frame(sample_data(ps_16S))

##------------------------------------------------------------------------------
## Beta-dispersion — full 16S dataset ####
##------------------------------------------------------------------------------
# PERMANOVA assumes equal within-group dispersion. Significant PERMDISP means
# compositional differences may reflect dispersion as well as centroid shifts.
# Expected in diverse microbiome datasets; downstream Mantel/MRM analyses are
# robust to heterogeneous dispersion.

# ---- Bray-Curtis ----
full_bd_order_bc  <- betadisper(full_ds_bray_16S, full_meta_16S$Clade_Order,     type = "centroid")
full_bd_wc_bc     <- betadisper(full_ds_bray_16S, full_meta_16S$env_broad_scale, type = "centroid")
full_bd_diet_bc   <- betadisper(full_ds_bray_16S, full_meta_16S$Diet,            type = "centroid")
full_bd_sample_bc <- betadisper(full_ds_bray_16S, full_meta_16S$env_medium,      type = "centroid")

perm_order_bc  <- permutest(full_bd_order_bc,  pairwise = FALSE, permutations = 999)
perm_wc_bc     <- permutest(full_bd_wc_bc,     pairwise = FALSE, permutations = 999)
perm_diet_bc   <- permutest(full_bd_diet_bc,   pairwise = FALSE, permutations = 999)
perm_sample_bc <- permutest(full_bd_sample_bc, pairwise = FALSE, permutations = 999)

# Convenience labels with group sizes for boxplot axes
order_labels  <- paste0(names(table(full_meta_16S$Clade_Order)),     "\n(n=", table(full_meta_16S$Clade_Order),     ")")
wc_labels     <- paste0(names(table(full_meta_16S$env_broad_scale)), "\n(n=", table(full_meta_16S$env_broad_scale), ")")
diet_labels   <- paste0(names(table(full_meta_16S$Diet)),            "\n(n=", table(full_meta_16S$Diet),            ")")
sample_labels <- paste0(names(table(full_meta_16S$env_medium)),      "\n(n=", table(full_meta_16S$env_medium),      ")")

boxplot(full_bd_order_bc,  names = order_labels,  las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Host Order (Bray-Curtis; full dataset)",
        xlab = "Host Order",    ylab = "Distance to centroid")
boxplot(full_bd_wc_bc,     names = wc_labels,     las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Management (Bray-Curtis; full dataset)",
        xlab = "Management",    ylab = "Distance to centroid")
boxplot(full_bd_diet_bc,   names = diet_labels,   las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Diet (Bray-Curtis; full dataset)",
        xlab = "Diet",          ylab = "Distance to centroid")
boxplot(full_bd_sample_bc, names = sample_labels, las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Sample Type (Bray-Curtis; full dataset)",
        xlab = "Sample Type",   ylab = "Distance to centroid")

# Annotate sample type plot with permutest F and p
f_val_bc <- perm_sample_bc$tab$F[1]
p_val_bc <- perm_sample_bc$tab$`Pr(>F)`[1]
label_bc <- paste0("permutest:\nF = ", round(f_val_bc, 2), ", p = ", signif(p_val_bc, 3))
usr <- par("usr")
text(usr[1] + 0.02*(usr[2]-usr[1]), usr[3] + 0.02*(usr[4]-usr[3]),
     labels = label_bc, adj = c(0, 0), cex = 0.9)

# ---- Jaccard ----
full_bd_order_jc  <- betadisper(full_ds_jacc_16S, full_meta_16S$Clade_Order,     type = "centroid")
full_bd_wc_jc     <- betadisper(full_ds_jacc_16S, full_meta_16S$env_broad_scale, type = "centroid")
full_bd_diet_jc   <- betadisper(full_ds_jacc_16S, full_meta_16S$Diet,            type = "centroid")
full_bd_sample_jc <- betadisper(full_ds_jacc_16S, full_meta_16S$env_medium,      type = "centroid")

perm_order_jc  <- permutest(full_bd_order_jc,  pairwise = FALSE, permutations = 999)
perm_wc_jc     <- permutest(full_bd_wc_jc,     pairwise = FALSE, permutations = 999)
perm_diet_jc   <- permutest(full_bd_diet_jc,   pairwise = FALSE, permutations = 999)
perm_sample_jc <- permutest(full_bd_sample_jc, pairwise = FALSE, permutations = 999)

boxplot(full_bd_order_jc,  names = order_labels,  las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Host Order (Jaccard; full dataset)",
        xlab = "Host Order",    ylab = "Distance to centroid")
boxplot(full_bd_wc_jc,     names = wc_labels,     las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Management (Jaccard; full dataset)",
        xlab = "Management",    ylab = "Distance to centroid")
boxplot(full_bd_diet_jc,   names = diet_labels,   las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Diet (Jaccard; full dataset)",
        xlab = "Diet",          ylab = "Distance to centroid")
boxplot(full_bd_sample_jc, names = sample_labels, las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Sample Type (Jaccard; full dataset)",
        xlab = "Sample Type",   ylab = "Distance to centroid")

f_val_jc <- perm_sample_jc$tab$F[1]
p_val_jc <- perm_sample_jc$tab$`Pr(>F)`[1]
label_jc <- paste0("permutest:\nF = ", round(f_val_jc, 2), ", p = ", signif(p_val_jc, 3))
usr <- par("usr")
text(usr[1] + 0.02*(usr[2]-usr[1]), usr[3] + 0.02*(usr[4]-usr[3]),
     labels = label_jc, adj = c(0, 0), cex = 0.9)

##------------------------------------------------------------------------------
## Subsampled beta-dispersion (16S) — check for sample-size bias ####
##------------------------------------------------------------------------------
# Subsample each sample type to the size of the smallest group (Lower GI).
# Confirms that dispersion differences are biological, not driven by unequal n.

sample_type_counts <- table(sample_data(ps_16S)$env_medium)
min_size <- min(sample_type_counts)
cat("Subsampling to n =", min_size, "per sample type group\n")

meta_sub_samp <- full_meta_16S %>%
  tibble::rownames_to_column("SampleID") %>%
  group_by(env_medium) %>%
  slice_sample(n = min_size) %>%
  pull(SampleID)

ps_16S_sub    <- prune_samples(meta_sub_samp, ps_16S)
full_meta_sub <- data.frame(sample_data(ps_16S_sub))

full_ds_bray_sub <- phyloseq::distance(ps_16S_sub, method = "bray")
full_ds_jacc_sub <- phyloseq::distance(ps_16S_sub, method = "jaccard", binary = TRUE)

sample_labels_sub <- paste0(names(table(full_meta_sub$env_medium)),
                            "\n(n=", table(full_meta_sub$env_medium), ")")

# Bray-Curtis
bd_sample_sub_bc <- betadisper(full_ds_bray_sub, full_meta_sub$env_medium, type = "centroid")
perm_sub_bc      <- permutest(bd_sample_sub_bc, pairwise = FALSE, permutations = 999)

boxplot(bd_sample_sub_bc, names = sample_labels_sub, las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Sample Type (Bray-Curtis; subsampled)",
        xlab = "Sample Type", ylab = "Distance to centroid")
f_sub_bc <- perm_sub_bc$tab$F[1]
p_sub_bc <- perm_sub_bc$tab$`Pr(>F)`[1]
label_sub_bc <- paste0("permutest:\nF = ", round(f_sub_bc, 2), ", p = ", signif(p_sub_bc, 3))
usr <- par("usr")
text(usr[1] + 0.38*(usr[2]-usr[1]), usr[3] + 0.03*(usr[4]-usr[3]),
     labels = label_sub_bc, adj = c(0, 0), cex = 0.9)

# Jaccard
bd_sample_sub_jc <- betadisper(full_ds_jacc_sub, full_meta_sub$env_medium, type = "centroid")
perm_sub_jc      <- permutest(bd_sample_sub_jc, pairwise = FALSE, permutations = 999)

boxplot(bd_sample_sub_jc, names = sample_labels_sub, las = 1, cex.axis = 0.7,
        main = "Beta-dispersion - Sample Type (Jaccard; subsampled)",
        xlab = "Sample Type", ylab = "Distance to centroid")
f_sub_jc <- perm_sub_jc$tab$F[1]
p_sub_jc <- perm_sub_jc$tab$`Pr(>F)`[1]
label_sub_jc <- paste0("permutest:\nF = ", round(f_sub_jc, 2), ", p = ", signif(p_sub_jc, 3))
usr <- par("usr")
text(usr[1] + 0.38*(usr[2]-usr[1]), usr[3] + 0.03*(usr[4]-usr[3]),
     labels = label_sub_jc, adj = c(0, 0), cex = 0.9)

# Note: If subsampled results match full-dataset results, dispersion patterns
# are biological rather than driven by unequal group sizes.

##------------------------------------------------------------------------------
## Save beta-dispersion boxplot figures ####
##------------------------------------------------------------------------------
boxplots_dir <- file.path(dir_figures, "boxplots")
dir.create(boxplots_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: save a betadisper boxplot as PDF
save_betadisper_plot <- function(bd_obj, labels, main_title, label_text,
                                 label_x_frac, fname_base) {
  pdf(paste0(fname_base, ".pdf"), width = 6, height = 5)
  boxplot(bd_obj, names = labels, las = 1, cex.axis = 0.7, cex.lab = 0.9,
          main = main_title, cex.main = 0.9,
          xlab = "Sample Type", ylab = "Distance to centroid")
  usr <- par("usr")
  text(usr[1] + label_x_frac*(usr[2]-usr[1]),
       usr[3] + 0.03*(usr[4]-usr[3]),
       labels = label_text, adj = c(0, 0), cex = 0.7)
  dev.off()
}

save_betadisper_plot(full_bd_sample_bc, sample_labels,
                     "Beta-dispersion for Sample Type (Bray-Curtis; full dataset)",
                     label_bc, 0.02,
                     file.path(boxplots_dir, "16S_boxplot_full_bray"))

save_betadisper_plot(full_bd_sample_jc, sample_labels,
                     "Beta-dispersion for Sample Type (Jaccard; full dataset)",
                     label_jc, 0.02,
                     file.path(boxplots_dir, "16S_boxplot_full_jaccard"))

save_betadisper_plot(bd_sample_sub_bc, sample_labels_sub,
                     "Beta-dispersion for Sample Type (Bray-Curtis; subsampled)",
                     label_sub_bc, 0.38,
                     file.path(boxplots_dir, "16S_boxplot_subsample_bray"))

save_betadisper_plot(bd_sample_sub_jc, sample_labels_sub,
                     "Beta-dispersion for Sample Type (Jaccard; subsampled)",
                     label_sub_jc, 0.38,
                     file.path(boxplots_dir, "16S_boxplot_subsample_jaccard"))

# 2x2 multipanel: full Bray | full Jaccard / subsampled Bray | subsampled Jaccard
pdf(file.path(boxplots_dir, "16S_boxplot_multipanel.pdf"), width = 8, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(2, 2, 2, 0))

boxplot(full_bd_sample_bc, names = sample_labels, las = 1, cex.axis = 0.7,
        xlab = "", ylab = "Distance to centroid",
        main = "A) Beta-dispersion (Bray-Curtis; full dataset)",
        cex.main = 0.9, cex.lab = 0.9)
usr <- par("usr")
text(usr[1]+0.02*(usr[2]-usr[1]), usr[3]+0.02*(usr[4]-usr[3]),
     labels = label_bc, adj = c(0,0), cex = 0.7)

boxplot(full_bd_sample_jc, names = sample_labels, las = 1, cex.axis = 0.7,
        xlab = "", ylab = "",
        main = "B) Beta-dispersion (Jaccard; full dataset)",
        cex.main = 0.9, cex.lab = 0.9)
usr <- par("usr")
text(usr[1]+0.02*(usr[2]-usr[1]), usr[3]+0.02*(usr[4]-usr[3]),
     labels = label_jc, adj = c(0,0), cex = 0.7)

boxplot(bd_sample_sub_bc, names = sample_labels_sub, las = 1, cex.axis = 0.7,
        xlab = "Sample Type", ylab = "Distance to centroid",
        main = "C) Beta-dispersion (Bray-Curtis; subsampled)",
        cex.main = 0.9, cex.lab = 0.9)
usr <- par("usr")
text(usr[1]+0.38*(usr[2]-usr[1]), usr[3]+0.03*(usr[4]-usr[3]),
     labels = label_sub_bc, adj = c(0,0), cex = 0.7)

boxplot(bd_sample_sub_jc, names = sample_labels_sub, las = 1, cex.axis = 0.7,
        xlab = "Sample Type", ylab = "",
        main = "D) Beta-dispersion (Jaccard; subsampled)",
        cex.main = 0.9, cex.lab = 0.9)
usr <- par("usr")
text(usr[1]+0.38*(usr[2]-usr[1]), usr[3]+0.03*(usr[4]-usr[3]),
     labels = label_sub_jc, adj = c(0,0), cex = 0.7)

dev.off()

##------------------------------------------------------------------------------
## PERMANOVA (16S) — full dataset ####
##------------------------------------------------------------------------------
full_bc_permanova_16S <- adonis2(
  full_ds_bray_16S ~ Clade_Order * Diet * env_broad_scale * env_medium,
  data         = full_meta_16S,
  strata       = full_meta_16S$site,
  permutations = 999,
  by           = "terms"
)

full_jacc_permanova_16S <- adonis2(
  full_ds_jacc_16S ~ Clade_Order * Diet * env_broad_scale * env_medium,
  data         = full_meta_16S,
  strata       = full_meta_16S$site,
  permutations = 999,
  by           = "terms"
)

saveRDS(full_bc_permanova_16S,   file.path(dir_processed, "16S_full_bc_permanova_sampletype.rds"))
saveRDS(full_jacc_permanova_16S, file.path(dir_processed, "16S_full_jacc_permanova_sampletype.rds"))

cat("=== 16S PERMANOVA (Bray-Curtis) ===\n"); print(full_bc_permanova_16S)
cat("=== 16S PERMANOVA (Jaccard) ===\n");     print(full_jacc_permanova_16S)

##------------------------------------------------------------------------------
## Pairwise PERMANOVA (16S) — Sample Type, full dataset ####
##------------------------------------------------------------------------------
# The omnibus 3-level Sample Type term above (Cloacal swab / Fecal / Lower GI)
# tests only whether the three groups differ overall — it does not address
# whether fecal and lower-GI specifically are similar enough to justify
# combining them downstream. Pairwise comparisons (BH-corrected across the
# 3 pairs) give the effect size that decision actually rests on.

pw_bc_sampletype_16S <- pairwise_permanova(
  dist_obj    = full_ds_bray_16S,
  metadata    = full_meta_16S,
  group_col   = "env_medium",
  strata_col  = "site",
  permutations = 999
)

pw_jc_sampletype_16S <- pairwise_permanova(
  dist_obj    = full_ds_jacc_16S,
  metadata    = full_meta_16S,
  group_col   = "env_medium",
  strata_col  = "site",
  permutations = 999
)

cat("\n=== 16S Pairwise PERMANOVA — Sample Type (Bray-Curtis) ===\n")
print(pw_bc_sampletype_16S)
cat("\n=== 16S Pairwise PERMANOVA — Sample Type (Jaccard) ===\n")
print(pw_jc_sampletype_16S)

# Save results
#saveRDS(pw_bc_sampletype_16S, file.path(dir_processed, "16S_pairwise_permanova_sampletype_bray.rds"))
#saveRDS(pw_jc_sampletype_16S, file.path(dir_processed, "16S_pairwise_permanova_sampletype_jaccard.rds"))
#write.csv(pw_bc_sampletype_16S, file.path(dir_tables, "16S_pairwise_permanova_sampletype_bray.csv"), row.names = FALSE)
#write.csv(pw_jc_sampletype_16S, file.path(dir_tables, "16S_pairwise_permanova_sampletype_jaccard.csv"), row.names = FALSE)

# Note: the Squamata-only sensitivity analysis below has only two Sample Type
# levels (Cloacal swab, Fecal), so its omnibus adonis2 term already *is* the
# pairwise comparison — no separate pairwise call is needed there.

##------------------------------------------------------------------------------
## Diagnose p = 1 in Lower GI pairwise comparisons ####
##------------------------------------------------------------------------------
# Both pairs involving Lower GI returned p_raw = 1.000 exactly, despite R2
# values as large as or larger than the (clearly significant) Cloacal vs
# Fecal pair. p = 1 exactly is a signature of a degenerate permutation space
# under `strata`, not evidence of no difference: Lower GI samples derive
# largely from dissected, euthanized individuals concentrated at a small
# number of sites (Methods). If a site stratum contains only Lower GI
# samples, within-stratum permutation can never reassign group labels there,
# pinning the observed F at the edge of an effectively fixed null
# distribution. Check the cross-tab below before interpreting these results.

check_strata_confound <- function(metadata, group_col, strata_col, pair) {
  sub <- metadata[metadata[[group_col]] %in% pair, ]
  tab <- table(sub[[strata_col]], sub[[group_col]])
  cat("\nSite x Sample Type cross-tab —", paste(pair, collapse = " vs "), "\n")
  print(tab)
  mono_strata <- rownames(tab)[rowSums(tab > 0) == 1]
  if (length(mono_strata) > 0) {
    cat("WARNING:", length(mono_strata),
        "site(s) contain only one group — degenerate for within-strata permutation:\n")
    print(mono_strata)
  } else {
    cat("No mono-group strata detected.\n")
  }
  invisible(tab)
}

check_strata_confound(full_meta_16S, "env_medium", "site", c("Cloacal swab", "Fecal"))
check_strata_confound(full_meta_16S, "env_medium", "site", c("Cloacal swab", "Lower GI"))
check_strata_confound(full_meta_16S, "env_medium", "site", c("Fecal", "Lower GI"))

##------------------------------------------------------------------------------
## Pairwise PERMANOVA (16S) — Sample Type, UNSTRATIFIED sensitivity check ####
##------------------------------------------------------------------------------
# If the diagnostic above confirms mono-group strata for Lower GI pairs,
# the stratified test is uninterpretable for those pairs specifically.
# Re-run without strata as a sensitivity check; report alongside the
# stratified result with a note on which is appropriate for which pair,
# rather than silently dropping strata for the whole analysis.

pw_bc_sampletype_16S_unstrat <- pairwise_permanova(
  dist_obj     = full_ds_bray_16S,
  metadata     = full_meta_16S,
  group_col    = "env_medium",
  strata_col   = NULL,
  permutations = 999
)

pw_jc_sampletype_16S_unstrat <- pairwise_permanova(
  dist_obj     = full_ds_jacc_16S,
  metadata     = full_meta_16S,
  group_col    = "env_medium",
  strata_col   = NULL,
  permutations = 999
)

cat("\n=== 16S Pairwise PERMANOVA — Sample Type, UNSTRATIFIED (Bray-Curtis) ===\n")
print(pw_bc_sampletype_16S_unstrat)
cat("\n=== 16S Pairwise PERMANOVA — Sample Type, UNSTRATIFIED (Jaccard) ===\n")
print(pw_jc_sampletype_16S_unstrat)

# Save results
#saveRDS(pw_bc_sampletype_16S_unstrat, file.path(dir_processed, "16S_pairwise_permanova_sampletype_bray_unstrat.rds"))
#saveRDS(pw_jc_sampletype_16S_unstrat, file.path(dir_processed, "16S_pairwise_permanova_sampletype_jaccard_unstrat.rds"))
#write.csv(pw_bc_sampletype_16S_unstrat, file.path(dir_tables, "16S_pairwise_permanova_sampletype_bray_unstrat.csv"), row.names = FALSE)
#write.csv(pw_jc_sampletype_16S_unstrat, file.path(dir_tables, "16S_pairwise_permanova_sampletype_jaccard_unstrat.csv"), row.names = FALSE)

bind_rows(
  pw_bc_sampletype_16S_unstrat %>% mutate(distance_metric = "bray-curtis", .before = 1),
  pw_jc_sampletype_16S_unstrat %>% mutate(distance_metric = "jaccard",     .before = 1)
) %>%
  write.csv(file.path(dir_tables, "16S_pairwise_permanova_sampletype_unstrat_combined.csv"), row.names = FALSE)

##------------------------------------------------------------------------------
## NMDS ordinations (16S) — colored by sample type ####
##------------------------------------------------------------------------------
nmds_dir <- file.path(dir_figures, "nmds_ordinations")
dir.create(nmds_dir, showWarnings = FALSE, recursive = TRUE)

ds.nmds.bray.16S.st <- metaMDS(full_ds_bray_16S, k = 5)
ds.nmds.jacc.16S.st <- metaMDS(full_ds_jacc_16S, k = 5)

cat("16S Bray-Curtis NMDS stress:", round(ds.nmds.bray.16S.st$stress, 4), "\n")
cat("16S Jaccard NMDS stress:",     round(ds.nmds.jacc.16S.st$stress, 4), "\n\n")

# Build NMDS dataframes
df.nmds.bray.16S.st <- as.data.frame(ds.nmds.bray.16S.st[["points"]]) %>%
  rownames_to_column("sample_name") %>%
  left_join(full_meta_16S, by = "sample_name")

df.nmds.jacc.16S.st <- as.data.frame(ds.nmds.jacc.16S.st[["points"]]) %>%
  rownames_to_column("sample_name") %>%
  left_join(full_meta_16S, by = "sample_name")

# Extract PERMANOVA R² and p for Sample Type (env_medium) — full dataset.
# Pulled from the omnibus model fit above (full_bc_permanova_16S /
# full_jacc_permanova_16S); same main-effect term reported in the main text.
r2_bc_full <- round(full_bc_permanova_16S["env_medium", "R2"], 3)
p_bc_full  <- full_bc_permanova_16S["env_medium", "Pr(>F)"]
p_bc_full_label <- ifelse(p_bc_full < 0.001, "p < 0.001",
                          paste0("p = ", signif(p_bc_full, 2)))

r2_jc_full <- round(full_jacc_permanova_16S["env_medium", "R2"], 3)
p_jc_full  <- full_jacc_permanova_16S["env_medium", "Pr(>F)"]
p_jc_full_label <- ifelse(p_jc_full < 0.001, "p < 0.001",
                          paste0("p = ", signif(p_jc_full, 2)))

# ---- Bray-Curtis NMDS (Panel A) ----
bc.nmds.st <- df.nmds.bray.16S.st %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = env_medium), linewidth = 1, linetype = "longdash") +
  geom_point(aes(color = env_medium, shape = Clade_Order), size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = env_medium), size = 7, shape = 21) +
  annotate("text",
           x = min(df.nmds.bray.16S.st$MDS1), y = min(df.nmds.bray.16S.st$MDS2),
           hjust = 0, vjust = 0.5, size = 3,
           label = paste0("Stress = ", round(ds.nmds.bray.16S.st$stress, 3),
                          "\nPERMANOVA Sample type: R\u00b2 = ", r2_bc_full,
                          ", ", p_bc_full_label)) +
  labs(title    = "A) Bray-Curtis NMDS",
       subtitle = "Bacteria; full dataset",
       fill = "Sample type", color = "Sample type", shape = "Host Order") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title    = element_text(face = "plain"),
        plot.subtitle = element_text(size = 10)) +
  guides(shape = guide_legend(order = 2),
         fill  = guide_legend(order = 1),
         color = guide_legend(order = 1)) +
  scale_color_manual(values = sampletype_colors) +
  scale_fill_manual(values  = sampletype_colors)

# ---- Jaccard NMDS (Panel B) ----
jacc.nmds.st <- df.nmds.jacc.16S.st %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = env_medium), linewidth = 1, linetype = "longdash") +
  geom_point(aes(color = env_medium, shape = Clade_Order), size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = env_medium), size = 7, shape = 21) +
  annotate("text",
           x = min(df.nmds.jacc.16S.st$MDS1), y = min(df.nmds.jacc.16S.st$MDS2),
           hjust = 0, vjust = 0.5, size = 3,
           label = paste0("Stress = ", round(ds.nmds.jacc.16S.st$stress, 3),
                          "\nPERMANOVA Sample type: R\u00b2 = ", r2_jc_full,
                          ", ", p_jc_full_label)) +
  labs(title    = "B) Jaccard NMDS",
       subtitle = "Bacteria; full dataset",
       fill = "Sample type", color = "Sample type", shape = "Host Order") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title    = element_text(face = "plain"),
        plot.subtitle = element_text(size = 10)) +
  guides(shape = guide_legend(order = 2),
         fill  = guide_legend(order = 1),
         color = guide_legend(order = 1)) +
  scale_color_manual(values = sampletype_colors) +
  scale_fill_manual(values  = sampletype_colors)

plot(bc.nmds.st)
plot(jacc.nmds.st)

ggsave(file.path(nmds_dir, "16S_bc_nmds_sampletype.pdf"), bc.nmds.st, device = cairo_pdf, width = 8, height = 5, units = "in")
ggsave(file.path(nmds_dir, "16S_jacc_nmds_sampletype.pdf"), jacc.nmds.st, device = cairo_pdf, width = 8, height = 5, units = "in")

# Stacked multipanel (Bray-Curtis top, Jaccard bottom)
combined_st <- (bc.nmds.st / jacc.nmds.st) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(nmds_dir, "16S_nmds_sampletype_multi.pdf"), combined_st, device = cairo_pdf, width = 10, height = 9, units = "in")

##------------------------------------------------------------------------------
## Taxonomy barplots — species with paired fecal + cloacal swab samples ####
##------------------------------------------------------------------------------
# Generates phylum-level relative abundance barplots for any species with
# >= 4 samples of both fecal and cloacal swab types, using consistent colors.


# Identify qualifying species
species_sample_counts <- full_meta_16S %>%
  dplyr::filter(env_medium %in% c("Fecal", "Cloacal swab")) %>%
  group_by(host_taxon, env_medium) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = env_medium, values_from = n, values_fill = 0) %>%
  dplyr::filter(Fecal > 3 & `Cloacal swab` > 3) %>%
  arrange(desc(Fecal + `Cloacal swab`))

species_vec   <- species_sample_counts$host_taxon
rank_to_plot  <- "Phylum"
top_n         <- 8
width_pdf     <- 10
height_pdf    <- 5
outdir_bars   <- file.path(dir_figures, "taxa_barplots")
dir.create(outdir_bars, showWarnings = FALSE, recursive = TRUE)

# Identify all top taxa across qualifying species for consistent color mapping
all_top_taxa <- c()
for (sp in species_vec) {
  ps_sp <- prune_samples(sample_data(ps_16S)$host_taxon == sp, ps_16S)
  if (nsamples(ps_sp) == 0) next
  ps_gl  <- tax_glom(ps_sp, taxrank = rank_to_plot, NArm = FALSE)
  ps_rel <- transform_sample_counts(ps_gl, function(x) if (sum(x) == 0) x else x / sum(x))
  df_sp  <- psmelt(ps_rel)
  if (!rank_to_plot %in% colnames(df_sp)) df_sp[[rank_to_plot]] <- "Unassigned"
  df_sp[[rank_to_plot]][is.na(df_sp[[rank_to_plot]]) | df_sp[[rank_to_plot]] == ""] <- "Unassigned"
  top_taxa <- df_sp %>%
    group_by(across(all_of(rank_to_plot))) %>%
    summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_ab)) %>%
    slice_head(n = top_n) %>%
    pull(1)
  all_top_taxa <- unique(c(all_top_taxa, top_taxa))
}

# Assign consistent colors
n_taxa     <- length(all_top_taxa)
max_colors <- brewer.pal.info["Paired", "maxcolors"]
if (n_taxa <= max_colors) {
  global_colors <- setNames(brewer.pal(n_taxa, "Paired"), all_top_taxa)
} else {
  global_colors <- setNames(colorRampPalette(brewer.pal(max_colors, "Paired"))(n_taxa), all_top_taxa)
}
global_colors <- c(global_colors, Other = "grey80")

# Generate and save per-species barplots
for (sp in species_vec) {
  message("Processing: ", sp)
  ps_sp <- prune_samples(sample_data(ps_16S)$host_taxon == sp, ps_16S)
  if (nsamples(ps_sp) == 0) { message("  No samples — skipping"); next }
  
  ps_gl  <- tax_glom(ps_sp, taxrank = rank_to_plot, NArm = FALSE)
  ps_rel <- transform_sample_counts(ps_gl, function(x) if (sum(x) == 0) x else x / sum(x))
  df_sp  <- psmelt(ps_rel)
  if (!rank_to_plot %in% colnames(df_sp)) df_sp[[rank_to_plot]] <- "Unassigned"
  df_sp[[rank_to_plot]][is.na(df_sp[[rank_to_plot]]) | df_sp[[rank_to_plot]] == ""] <- "Unassigned"
  
  top_taxa_sp   <- intersect(all_top_taxa, df_sp[[rank_to_plot]])
  df_sp$tax_group <- ifelse(df_sp[[rank_to_plot]] %in% top_taxa_sp, df_sp[[rank_to_plot]], "Other")
  df_sp$Sample     <- factor(df_sp$Sample, levels = unique(df_sp$Sample))
  df_sp$env_medium <- factor(df_sp$env_medium)
  
  p <- ggplot(df_sp, aes(x = Sample, y = Abundance, fill = tax_group)) +
    geom_col(width = 0.8) +
    facet_wrap(~ env_medium, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = global_colors, name = "Top Bacterial Phyla") +
    theme_light(base_size = 14) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(size = 12, face = "bold", color = "black"),
      legend.title     = element_text(size = 12),
      legend.text      = element_text(size = 10)
    ) +
    labs(title = sp, x = "Host ID", y = "Relative Abundance")
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", sp)
  ggsave(file.path(outdir_bars, paste0("16S_", safe_name, "_taxa_barplot.pdf")),
         p, width = width_pdf, height = height_pdf)
  ggsave(file.path(outdir_bars, paste0("16S_", safe_name, "_taxa_barplot.png")),
         p, width = width_pdf, height = height_pdf, dpi = 300)
  message("  Saved: ", sp)
}
message("All barplots saved to: ", normalizePath(outdir_bars))

# Multipanel: Coluber constrictor (A) / Crotalus viridis (B)
make_species_plot <- function(sp, panel_label) {
  ps_sp <- prune_samples(sample_data(ps_16S)$host_taxon == sp, ps_16S)
  ps_gl  <- tax_glom(ps_sp, taxrank = rank_to_plot, NArm = FALSE)
  ps_rel <- transform_sample_counts(ps_gl, function(x) if (sum(x) == 0) x else x / sum(x))
  df_sp  <- psmelt(ps_rel)
  if (!rank_to_plot %in% colnames(df_sp)) df_sp[[rank_to_plot]] <- "Unassigned"
  df_sp[[rank_to_plot]][is.na(df_sp[[rank_to_plot]]) | df_sp[[rank_to_plot]] == ""] <- "Unassigned"
  top_taxa_sp     <- intersect(all_top_taxa, df_sp[[rank_to_plot]])
  df_sp$tax_group <- ifelse(df_sp[[rank_to_plot]] %in% top_taxa_sp, df_sp[[rank_to_plot]], "Other")
  df_sp$Sample     <- factor(df_sp$Sample, levels = unique(df_sp$Sample))
  df_sp$env_medium <- factor(df_sp$env_medium)
  
  ggplot(df_sp, aes(x = Sample, y = Abundance, fill = tax_group)) +
    geom_col(width = 0.8) +
    facet_wrap(~ env_medium, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = global_colors, name = "Top Bacterial Phyla") +
    theme_light(base_size = 14) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(size = 12, face = "bold", color = "black"),
      legend.title     = element_text(size = 12),
      legend.text      = element_text(size = 10)
    ) +
    labs(title = paste0(panel_label, ") ", sp), x = "Host ID", y = "Relative Abundance")
}

p_cc <- make_species_plot("Coluber constrictor mormon", "A")
p_cv <- make_species_plot("Crotalus viridis",    "B")

multi_species <- p_cc / p_cv +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(dir_figures, "16S_Coluber_Crotalus_multipanel.pdf"),
       multi_species, device = cairo_pdf, width = width_pdf, height = height_pdf * 2, units = "in")

##------------------------------------------------------------------------------
## Snake-only sensitivity analysis: PERMDISP + PERMANOVA + NMDS ####
##------------------------------------------------------------------------------
# Rationale: Lizards in this dataset carry only fecal samples; snakes are the
# only Squamata with both fecal and cloacal swab samples. Restricting to snakes
# provides the cleanest within-lineage test that sample type drives community
# separation independently of host clade or sample-type confounding with
# animal_type. Ellipse color = sample type (env_medium); point shape = sample
# type (env_medium) — no animal_type shape mapping needed since all points
# are snakes.


# ---- Subset to snakes only ----
ps_16S_snake <- subset_samples(ps_16S, animal_type == "Snake")
ps_16S_snake <- prune_taxa(taxa_sums(ps_16S_snake) > 0, ps_16S_snake)
meta_snake   <- data.frame(sample_data(ps_16S_snake))

cat("\nSnake subset — sample type counts:\n")
print(table(meta_snake$env_medium))
cat("\nSnake subset — animal type check (should be Snake only):\n")
print(table(meta_snake$animal_type))

# ---- Distance matrices ----
ds_bray_snake <- phyloseq::distance(ps_16S_snake, method = "bray")
ds_jacc_snake <- phyloseq::distance(ps_16S_snake, method = "jaccard", binary = TRUE)

# ---- Beta-dispersion (PERMDISP) — env_medium only ----
# Tests whether within-group dispersion differs between sample types.
# Significant result means communities differ in spread as well as centroid.
# NOTE: retained as a console-reported robustness check only (for a single
# Methods-text sentence confirming PERMANOVA results aren't a dispersion
# artifact) — no longer rendered as its own supplementary panel. Dispersion
# direction reversed between Bray-Curtis and Jaccard in this snake-only
# subset, which would be confusing to present visually without equivalent
# checks elsewhere in the paper; see Fig. S8 (NMDS-only, both metrics).
bd_snake_bc <- betadisper(ds_bray_snake, meta_snake$env_medium, type = "centroid")
bd_snake_jc <- betadisper(ds_jacc_snake, meta_snake$env_medium, type = "centroid")

perm_snake_bc <- permutest(bd_snake_bc, pairwise = FALSE, permutations = 999)
perm_snake_jc <- permutest(bd_snake_jc, pairwise = FALSE, permutations = 999)

cat("\n=== Snake PERMDISP — Bray-Curtis ===\n"); print(perm_snake_bc)
cat("\n=== Snake PERMDISP — Jaccard ===\n");     print(perm_snake_jc)

# Convenience labels with group sizes for betadisper boxplots
snake_sample_labels <- paste0(names(table(meta_snake$env_medium)),
                              "\n(n=", table(meta_snake$env_medium), ")")

# Extract PERMDISP stats for annotation
f_bd_bc <- perm_snake_bc$tab$F[1]
p_bd_bc <- perm_snake_bc$tab$`Pr(>F)`[1]
label_bd_bc <- paste0("permutest:\nF = ", round(f_bd_bc, 2),
                      ", p = ", signif(p_bd_bc, 3))

f_bd_jc <- perm_snake_jc$tab$F[1]
p_bd_jc <- perm_snake_jc$tab$`Pr(>F)`[1]
label_bd_jc <- paste0("permutest:\nF = ", round(f_bd_jc, 2),
                      ", p = ", signif(p_bd_jc, 3))

# ---- PERMANOVA (adonis2) ----
# env_medium + Diet; env_broad_scale dropped — confounded with site strata
# (wild/captive maps ~1:1 with collection site in the snake subset, producing
# p = 1.000 in the stratified model; see cross-tab above). Family may be added
# here pending the diagnostic output above.
snake_bc_permanova <- adonis2(
  ds_bray_snake ~ env_medium + Diet,
  data         = meta_snake,
  strata       = meta_snake$site,
  permutations = 999,
  by           = "terms"
)

snake_jacc_permanova <- adonis2(
  ds_jacc_snake ~ env_medium + Diet,
  data         = meta_snake,
  strata       = meta_snake$site,
  permutations = 999,
  by           = "terms"
)

cat("\n=== Snake PERMANOVA — Bray-Curtis ===\n"); print(snake_bc_permanova)
cat("\n=== Snake PERMANOVA — Jaccard ===\n");     print(snake_jacc_permanova)

#saveRDS(snake_bc_permanova,   file.path(dir_processed, "16S_snake_bc_permanova_sampletype.rds"))
#saveRDS(snake_jacc_permanova, file.path(dir_processed, "16S_snake_jacc_permanova_sampletype.rds"))

# ---- NMDS ordinations ----
ds.nmds.bray.snake <- metaMDS(ds_bray_snake, k = 5)
ds.nmds.jacc.snake <- metaMDS(ds_jacc_snake, k = 5)

cat("Snake Bray-Curtis NMDS stress:", round(ds.nmds.bray.snake$stress, 4), "\n")
cat("Snake Jaccard NMDS stress:",     round(ds.nmds.jacc.snake$stress, 4), "\n\n")

# Build ordination data frames
df.nmds.bray.snake <- as.data.frame(ds.nmds.bray.snake[["points"]]) %>%
  rownames_to_column("sample_name") %>%
  left_join(meta_snake, by = "sample_name") %>%
  mutate(env_medium = factor(env_medium))

df.nmds.jacc.snake <- as.data.frame(ds.nmds.jacc.snake[["points"]]) %>%
  rownames_to_column("sample_name") %>%
  left_join(meta_snake, by = "sample_name") %>%
  mutate(env_medium = factor(env_medium))

# Extract PERMANOVA R² and p for env_medium for plot annotations
r2_bc_snake <- round(snake_bc_permanova["env_medium", "R2"], 3)
p_bc_snake  <- snake_bc_permanova["env_medium", "Pr(>F)"]
p_bc_snake_label <- ifelse(p_bc_snake < 0.001, "p < 0.001",
                           paste0("p = ", signif(p_bc_snake, 2)))

r2_jc_snake <- round(snake_jacc_permanova["env_medium", "R2"], 3)
p_jc_snake  <- snake_jacc_permanova["env_medium", "Pr(>F)"]
p_jc_snake_label <- ifelse(p_jc_snake < 0.001, "p < 0.001",
                           paste0("p = ", signif(p_jc_snake, 2)))

# ---- Bray-Curtis NMDS — snakes ----
snake.bc.nmds.plot <- df.nmds.bray.snake %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = env_medium), linewidth = 1, linetype = "longdash") +
  geom_point(aes(color = env_medium, shape = env_medium), size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = env_medium), size = 7, shape = 21) +
  annotate("text",
           x = min(df.nmds.bray.snake$MDS1), y = min(df.nmds.bray.snake$MDS2),
           hjust = 0, vjust = 0.5, size = 3,
           label = paste0("Stress = ", round(ds.nmds.bray.snake$stress, 3),
                          "\nPERMANOVA env_medium: R\u00b2 = ", r2_bc_snake,
                          ", ", p_bc_snake_label)) +
  labs(title    = "C) Bray-Curtis NMDS",
       subtitle = "Bacteria; snakes-only",
       fill     = "Sample type",
       color    = "Sample type",
       shape    = "Sample type") +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "plain"),
    plot.subtitle   = element_text(size = 10)
  ) +
  guides(fill  = guide_legend(order = 1),
         color = guide_legend(order = 1),
         shape = guide_legend(order = 1)) +
  scale_color_manual(values = sampletype_colors) +
  scale_fill_manual(values  = sampletype_colors) +
  scale_shape_manual(values = c("Cloacal swab" = 16, "Fecal" = 17))

# ---- Jaccard NMDS — snakes ----
snake.jacc.nmds.plot <- df.nmds.jacc.snake %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = env_medium), linewidth = 1, linetype = "longdash") +
  geom_point(aes(color = env_medium, shape = env_medium), size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = env_medium), size = 7, shape = 21) +
  annotate("text",
           x = min(df.nmds.jacc.snake$MDS1), y = min(df.nmds.jacc.snake$MDS2),
           hjust = 0, vjust = 0.5, size = 3,
           label = paste0("Stress = ", round(ds.nmds.jacc.snake$stress, 3),
                          "\nPERMANOVA env_medium: R\u00b2 = ", r2_jc_snake,
                          ", ", p_jc_snake_label)) +
  labs(title    = "D) Jaccard NMDS",
       subtitle = "Bacteria; snakes-only",
       fill     = "Sample type",
       color    = "Sample type",
       shape    = "Sample type") +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "plain"),
    plot.subtitle   = element_text(size = 10)
  ) +
  guides(fill  = guide_legend(order = 1),
         color = guide_legend(order = 1),
         shape = guide_legend(order = 1)) +
  scale_color_manual(values = sampletype_colors) +
  scale_fill_manual(values  = sampletype_colors) +
  scale_shape_manual(values = c("Cloacal swab" = 16, "Fecal" = 17))

# ---- Save individual plots ----
plot(snake.bc.nmds.plot)
plot(snake.jacc.nmds.plot)

ggsave(file.path(nmds_dir, "16S_snake_bc_nmds_sampletype.pdf"),
       snake.bc.nmds.plot, device = cairo_pdf, width = 8, height = 5, units = "in")

ggsave(file.path(nmds_dir, "16S_snake_jacc_nmds_sampletype.pdf"),
       snake.jacc.nmds.plot, device = cairo_pdf, width = 8, height = 5, units = "in")

# ---- 2x2 multipanel: A=full Bray | B=full Jaccard / C=snake Bray | D=snake Jaccard ----
# Betadisper panels (previously A/B here) have been dropped: the dispersion
# direction was inconsistent between Bray-Curtis and Jaccard in the snake-only
# subset (see snake-only betadisper output above), and betadisper was not run
# as a standard check elsewhere in this analysis. PERMANOVA + NMDS already
# establish the centroid (location) differences this figure is meant to show;
# all four panels below are ggplot NMDS objects with PERMANOVA stats inset.

sampletype_nmds_4panel <- (bc.nmds.st | jacc.nmds.st) / (snake.bc.nmds.plot | snake.jacc.nmds.plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(dir_figures, "16S_sampletype_nmds_4panel_FigS8.pdf"),
       sampletype_nmds_4panel, device = cairo_pdf, width = 12, height = 9, units = "in")

cat("4-panel sample type NMDS figure (Fig. S8) saved: full dataset (A,B) + snakes-only (C,D).\n")

##------------------------------------------------------------------------------
## Split 16S phyloseq by sample type and save ####
##------------------------------------------------------------------------------
cat("\n16S sample type counts:\n")
print(table(sample_data(ps_16S)$env_medium))

# Fecal + Lower GI combined (treated as equivalent GI samples for downstream analyses)
ps_16S_fecal_lowerGI <- ps_16S %>%
  subset_samples(env_medium %in% c("Fecal", "Lower GI")) %>%
  prune_taxa(taxa_sums(.) > 0, .)

# Cloacal swab samples (analyzed separately)
ps_16S_cloacal <- ps_16S %>%
  subset_samples(env_medium == "Cloacal swab") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# Sanity checks
cat("\n--- 16S Fecal + Lower GI ---\n"); print(ps_16S_fecal_lowerGI)
print(table(sample_data(ps_16S_fecal_lowerGI)$env_medium))

cat("\n--- 16S Cloacal ---\n"); print(ps_16S_cloacal)
print(table(sample_data(ps_16S_cloacal)$env_medium))

saveRDS(ps_16S_fecal_lowerGI, file.path(dir_processed, "16S_abs_final_fecal.rds"))
saveRDS(ps_16S_cloacal,       file.path(dir_processed, "16S_abs_final_cloacal.rds"))
cat("16S split phyloseq objects saved.\n")

# Write as .csv
meta_16S_fecal <- full_meta_16S[full_meta_16S$env_medium %in% c("Fecal", "Lower GI"), ]
write.csv(meta_16S_fecal, file.path(dir_tables, "16S_abs_final_fecal_metadata.csv"))

###############################################################################
# ========================== ITS DATASET ======================================
###############################################################################

# ========================== LOAD DATA =========================================
ps_ITS <- readRDS(file.path(dir_processed, "ITS_absolute_all_final.rds"))

full_meta_ITS <- data.frame(sample_data(ps_ITS))
cat("\nITS sample type counts:\n")
print(table(full_meta_ITS$env_medium))
cat("\nITS Kingdom check:\n")
print(table(tax_table(ps_ITS)[, "Kingdom"], useNA = "ifany"))

##------------------------------------------------------------------------------
## Split ITS phyloseq by sample type and save ####
##------------------------------------------------------------------------------

cat("\nITS sample type counts:\n")
print(table(sample_data(ps_ITS)$env_medium))

ps_ITS_fecal_lowerGI <- ps_ITS %>%
  subset_samples(env_medium %in% c("Fecal", "Lower GI")) %>%
  prune_taxa(taxa_sums(.) > 0, .)

ps_ITS_cloacal <- ps_ITS %>%
  subset_samples(env_medium == "Cloacal swab") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# Sanity checks
cat("\n--- ITS Fecal + Lower GI ---\n"); print(ps_ITS_fecal_lowerGI)
print(table(tax_table(ps_ITS_fecal_lowerGI)[, "Kingdom"], useNA = "ifany"))
print(table(sample_data(ps_ITS_fecal_lowerGI)$env_medium))

cat("\n--- ITS Cloacal ---\n"); print(ps_ITS_cloacal)
print(table(sample_data(ps_ITS_cloacal)$env_medium))

saveRDS(ps_ITS_fecal_lowerGI, file.path(dir_processed, "ITS_abs_final_fecal.rds"))
saveRDS(ps_ITS_cloacal,       file.path(dir_processed, "ITS_abs_final_cloacal.rds"))
cat("ITS split phyloseq objects saved.\n")

# Write as .csv
meta_ITS_fecal <- full_meta_ITS[full_meta_ITS$env_medium %in% c("Fecal", "Lower GI"), ]
write.csv(meta_ITS_fecal, file.path(dir_tables, "ITS_abs_final_fecal_metadata.csv"))

###############################################################################
# ========================== DATASET SUMMARIES ================================
###############################################################################
# Concise summary of each split phyloseq object for reporting and QC.

summarize_ps <- function(ps, label) {
  meta  <- data.frame(sample_data(ps))
  reads <- sample_sums(ps)
  taxa  <- tax_table(ps)
  
  # Top 5 taxa at phylum level by mean relative abundance
  ps_rel   <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_phyl  <- tax_glom(ps_rel, taxrank = "Phylum", NArm = TRUE)
  top_phyla <- tapply(taxa_sums(ps_phyl), tax_table(ps_phyl)[, "Phylum"], sum)
  top_phyla <- sort(top_phyla, decreasing = TRUE)
  top5      <- head(names(top_phyla), 5)
  
  cat("\n", strrep("=", 60), "\n", label, "\n", strrep("=", 60), "\n", sep = "")
  cat("Samples:        ", nsamples(ps), "\n")
  cat("OTUs:           ", ntaxa(ps), "\n")
  cat("Total reads:    ", format(sum(reads), big.mark = ","), "\n")
  cat("Mean reads:     ", format(round(mean(reads)), big.mark = ","), "\n")
  cat("Median reads:   ", format(round(median(reads)), big.mark = ","), "\n")
  cat("Range reads:    ", format(min(reads), big.mark = ","), "-",
      format(max(reads), big.mark = ","), "\n")
  cat("Host species:   ", n_distinct(meta$host_taxon), "\n")
  cat("Host families:  ", n_distinct(meta$Family), "\n")
  cat("Host orders:    ", n_distinct(meta$Clade_Order),
      paste0("(", paste(sort(unique(meta$Clade_Order)), collapse = ", "), ")"), "\n")
  cat("Wild samples:   ", sum(meta$env_broad_scale == "Wild"),
      paste0("(", round(100 * mean(meta$env_broad_scale == "Wild"), 1), "%)"), "\n")
  cat("Sample types:   "); print(table(meta$env_medium))
  cat("Top 5 phyla:    ", paste(top5, collapse = ", "), "\n")
}

summarize_ps(ps_16S_fecal_lowerGI, "16S — Fecal + Lower GI")
summarize_ps(ps_16S_cloacal,       "16S — Cloacal Swab")
summarize_ps(ps_ITS_fecal_lowerGI, "ITS — Fecal + Lower GI")
summarize_ps(ps_ITS_cloacal,       "ITS — Cloacal Swab")


###############################################################################
###############################################################################
###############################################################################
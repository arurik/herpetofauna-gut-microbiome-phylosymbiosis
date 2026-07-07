###############################################################################
# Gut microbiome community structure (β-diversity) — 16S & ITS full datasets
#
# For each dataset:
#   1. Remove outlier samples (if applicable)
#   2. Compute Bray-Curtis and Jaccard distance matrices
#   3. Beta-dispersion tests (PERMDISP) for key metadata factors
#   4. PERMANOVA (adonis2) — Bray-Curtis and Jaccard
#   5. NMDS ordinations — Bray-Curtis and Jaccard
#   6. Multipanel figure assembly (patchwork)
#   7. Export combined PERMANOVA results table
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
  library(ggtext)       
  library(kableExtra)
  library(webshot2)
  library(magrittr)
})

select <- dplyr::select
filter <- dplyr::filter

set.seed(52325)

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")

# ========================== SHARED COLOR PALETTE ==============================
clade_colors <- c(
  "Anura"      = "#CC79A7",
  "Caudata"    = "#E69F00",
  "Testudines" = "#009E73",
  "Squamata"   = "#18B1E9",
  "Crocodilia" = "#111189"
)

# ========================== HELPER: PERMANOVA ANNOTATION LABEL ================
# Extracts main effects from an adonis2 object and formats as HTML richtext
# for use in ggplot annotate(geom = "richtext"). Interactions excluded.
# R² rendered with superscript via R<sup>2</sup>.

make_permanova_label <- function(permanova_obj, distance_name) {
  res <- as.data.frame(permanova_obj)
  res$term <- rownames(res)
  
  # Rename Pr(>F) immediately to avoid special-character issues
  # in downstream dplyr/base operations
  names(res)[names(res) == "Pr(>F)"] <- "pval"
  
  main_effects <- res %>%
    dplyr::filter(!term %in% c("Residual", "Total")) %>%
    dplyr::filter(!grepl(":", term))
  
  label_map <- c(
    "Clade_Order"     = "Host Order",
    "Diet"            = "Diet",
    "env_broad_scale" = "Management",
    "env_medium"      = "Sample type"
  )
  
  paste0(
    "<b>PERMANOVA (", distance_name, ")</b><br>",
    paste(sapply(seq_len(nrow(main_effects)), function(i) {
      term     <- main_effects$term[i]
      r2_fmt   <- sprintf("%.3f", main_effects$R2[i])
      pval     <- main_effects$pval[i]                   # use renamed column
      pval_fmt <- ifelse(pval < 0.001, "< 0.001", sprintf("= %.3f", pval))
      display  <- ifelse(term %in% names(label_map), label_map[term], term)
      paste0(display, ": R<sup>2</sup> = ", r2_fmt, ", p ", pval_fmt)
    }), collapse = "<br>")
  )
}

# ========================== HELPER: BETADISPER + SAVE =========================
# Runs betadisper (centroid method, matching original scripts) and permutest,
# saves boxplot to subfolder as PDF, and prints to active device.

run_betadisper <- function(dist_mat, group_vec, label, distance_name,
                           dataset_name, outdir) {
  # Added type = "centroid" to match original scripts
  bd  <- betadisper(dist_mat, group_vec, type = "centroid")
  pt  <- permutest(bd, pairwise = FALSE, permutations = 999)
  lbs <- paste0(names(table(group_vec)), "\n(n=", table(group_vec), ")")
  
  plot_title    <- paste0("Beta-dispersion — ", label, " (", distance_name, ")")
  safe_label    <- gsub(" ", "_", tolower(label))
  safe_distance <- gsub("-", "", tolower(distance_name))
  fbase         <- file.path(outdir,
                             paste0(dataset_name, "_betadisper_",
                                    safe_label, "_", safe_distance))
  
  pdf(paste0(fbase, ".pdf"), width = 7, height = 5)
  boxplot(bd, names = lbs, las = 1, cex.axis = 0.7,
          main = plot_title, xlab = label, ylab = "Distance to centroid")
  dev.off()
  
  boxplot(bd, names = lbs, las = 1, cex.axis = 0.7,
          main = plot_title, xlab = label, ylab = "Distance to centroid")
  
  cat("PERMDISP —", label, "(", distance_name, "):\n")
  print(pt)
  return(list(bd = bd, permutest = pt))
}

# ========================== HELPER: FORMAT PERMANOVA TABLE ====================
# Converts an adonis2 object to a tidy data frame with dataset/distance labels.
# [FIX Bug 3 & 4] Renames F and Pr(>F) immediately after coercion to avoid
# conflicts with base R's F (= FALSE) and special-character column names.

format_permanova <- function(permanova_obj, dataset_label, distance_label) {
  
  factor_map <- c(
    "Clade_Order"     = "Host Order",
    "Diet"            = "Diet",
    "env_broad_scale" = "Management Status",
    "env_medium"      = "Sample Type"
  )
  
  res <- as.data.frame(permanova_obj)
  res$term <- rownames(res)
  
  # Rename problematic columns immediately
  names(res)[names(res) == "F"]       <- "Fstat"
  names(res)[names(res) == "Pr(>F)"]  <- "pval"
  
  res %>%
    mutate(
      factor = case_when(
        term == "Residual" ~ "Residual",
        term == "Total"    ~ "Total",
        TRUE ~ {
          t <- term
          for (raw in names(factor_map)) {
            t <- gsub(raw, factor_map[raw], t, fixed = TRUE)
          }
          t
        }
      ),
      dataset         = dataset_label,
      distance_matrix = distance_label,
      df              = Df,
      `sum of squares` = round(SumOfSqs, 3),
      `R²`            = round(R2, 5),
      F               = ifelse(is.na(Fstat), "-", sprintf("%.4f", Fstat)),  
      `p-value`       = ifelse(is.na(pval),  "-", sprintf("%.3f", pval)),   
      significance    = case_when(
        is.na(pval)    ~ "-",
        pval < 0.001   ~ "***",
        pval < 0.01    ~ "**",
        pval < 0.05    ~ "*",
        pval < 0.1     ~ ".",
        TRUE           ~ "ns"
      )
    ) %>%
    select(factor, dataset, distance_matrix, df,
           `sum of squares`, `R²`, F, `p-value`, significance)
}

###############################################################################
# ========================== 16S DATASET ======================================
###############################################################################

# ========================== LOAD DATA =========================================
ps_16S <- readRDS(file.path(dir_processed, "16S_absolute_all_final.rds"))

# ========================== HABITAT CLEANUP (16S) =============================
# "Fossorial" is a very small group with lower dispersion than other habitat
# categories. Collapse into "Fossorial-Terrestrial" to avoid instability.
full_meta_16S <- data.frame(sample_data(ps_16S))

# sample_data already contains a "sample_name" column — no rownames manipulation needed.

full_meta_16S$animal_ecomode <- as.character(full_meta_16S$animal_ecomode)

if ("Fossorial" %in% full_meta_16S$animal_ecomode) {
  n_fossorial <- sum(full_meta_16S$animal_ecomode == "Fossorial")
  cat("16S: Collapsing", n_fossorial,
      "'Fossorial' samples into 'Fossorial-Terrestrial'\n")
  full_meta_16S$animal_ecomode[full_meta_16S$animal_ecomode == "Fossorial"] <-
    "Fossorial-Terrestrial"
} else {
  cat("16S: No 'Fossorial' entries found — no collapse needed\n")
}
full_meta_16S$animal_ecomode <- factor(full_meta_16S$animal_ecomode)

# Propagate corrected ecomode back into the phyloseq sample_data
# sample_name is already a column AND the rownames — assign rownames directly
meta_for_ps_16S <- full_meta_16S
rownames(meta_for_ps_16S) <- meta_for_ps_16S$sample_name
sample_data(ps_16S) <- sample_data(meta_for_ps_16S)

# ========================== DISTANCE MATRICES (16S) ==========================
full_ds_bray_16S <- phyloseq::distance(ps_16S, method = "bray")
full_ds_jacc_16S <- phyloseq::distance(ps_16S, method = "jaccard", binary = TRUE)

cat("16S distance matrices computed:",
    nrow(as.matrix(full_ds_bray_16S)), "samples\n\n")

 ##------------------------------------------------------------------------------
## Beta-dispersion (16S) ####
##------------------------------------------------------------------------------
# PERMANOVA assumes equal dispersion. Significant PERMDISP indicates unequal
# dispersion — interpret PERMANOVA results as centroid shifts + dispersion
# differences. Downstream phylogenetically explicit analyses (partial Mantel,
# MRM) are robust to dispersion differences.

betadisper_dir_16S <- file.path(dir_figures, "betadisper_plots", "16S")
dir.create(betadisper_dir_16S, recursive = TRUE, showWarnings = FALSE)
dir.create(betadisper_dir_16S, recursive = TRUE, showWarnings = FALSE)

# Make sure R Studio's "Plots" tab is large enough for figure viewing

# Bray-Curtis
bd_16S_bc_order   <- run_betadisper(full_ds_bray_16S, full_meta_16S$Clade_Order,     "Host Order",  "Bray-Curtis", "16S", betadisper_dir_16S)
bd_16S_bc_wc      <- run_betadisper(full_ds_bray_16S, full_meta_16S$env_broad_scale, "Management",  "Bray-Curtis", "16S", betadisper_dir_16S)
bd_16S_bc_diet    <- run_betadisper(full_ds_bray_16S, full_meta_16S$Diet,             "Diet",        "Bray-Curtis", "16S", betadisper_dir_16S)
bd_16S_bc_sample  <- run_betadisper(full_ds_bray_16S, full_meta_16S$env_medium,       "Sample Type", "Bray-Curtis", "16S", betadisper_dir_16S)
bd_16S_bc_habitat <- run_betadisper(full_ds_bray_16S, full_meta_16S$animal_ecomode,   "Habitat",     "Bray-Curtis", "16S", betadisper_dir_16S)

# Jaccard
bd_16S_jc_order   <- run_betadisper(full_ds_jacc_16S, full_meta_16S$Clade_Order,     "Host Order",  "Jaccard", "16S", betadisper_dir_16S)
bd_16S_jc_wc      <- run_betadisper(full_ds_jacc_16S, full_meta_16S$env_broad_scale, "Management",  "Jaccard", "16S", betadisper_dir_16S)
bd_16S_jc_diet    <- run_betadisper(full_ds_jacc_16S, full_meta_16S$Diet,             "Diet",        "Jaccard", "16S", betadisper_dir_16S)
bd_16S_jc_sample  <- run_betadisper(full_ds_jacc_16S, full_meta_16S$env_medium,       "Sample Type", "Jaccard", "16S", betadisper_dir_16S)
bd_16S_jc_habitat <- run_betadisper(full_ds_jacc_16S, full_meta_16S$animal_ecomode,   "Habitat",     "Jaccard", "16S", betadisper_dir_16S)

##------------------------------------------------------------------------------
## PERMANOVA (16S) ####
##------------------------------------------------------------------------------
# Main effects + interactions, stratified by site to account for sampling design

# Bray-Curtis
full_bc_permanova_16S <- adonis2(
  full_ds_bray_16S ~ Clade_Order * Diet * env_broad_scale * env_medium,
  data         = full_meta_16S,
  strata       = full_meta_16S$site,
  permutations = 999,
  by           = "terms"
)

# Jaccard
full_jacc_permanova_16S <- adonis2(
  full_ds_jacc_16S ~ Clade_Order * Diet * env_broad_scale * env_medium,
  data         = full_meta_16S,
  strata       = full_meta_16S$site,
  permutations = 999,
  by           = "terms"
)

# Save results
#saveRDS(full_bc_permanova_16S,   file = file.path(dir_processed, "16S_full_bc_permanova.rds"))
#saveRDS(full_jacc_permanova_16S, file = file.path(dir_processed, "16S_full_jacc_permanova.rds"))

cat("=== 16S PERMANOVA (Bray-Curtis) ===\n"); print(full_bc_permanova_16S)
cat("=== 16S PERMANOVA (Jaccard) ===\n");     print(full_jacc_permanova_16S)

##------------------------------------------------------------------------------
## NMDS Ordinations (16S) ####
##------------------------------------------------------------------------------
# Non-metric multidimensional scaling in k=5 dimensions.
# Stress value quantifies how well the ordination preserves pairwise distances.
# previous.best = FALSE on Jaccard prevents it from inheriting the Bray-Curtis
# solution as its starting configuration, which would produce near-identical plots.

ds.nmds.bray.16S <- metaMDS(full_ds_bray_16S, k = 5)
ds.nmds.jacc.16S <- metaMDS(full_ds_jacc_16S, k = 5)

cat("16S Bray-Curtis NMDS stress:", round(ds.nmds.bray.16S$stress, 4), "\n")
cat("16S Jaccard NMDS stress:",     round(ds.nmds.jacc.16S$stress, 4), "\n\n")

# Save NMDS results
#saveRDS(ds.nmds.bray.16S, file = file.path(dir_processed, "16S_nmds_bray.rds"))
#saveRDS(ds.nmds.jacc.16S, file = file.path(dir_processed, "16S_nmds_jacc.rds"))

# Build NMDS dataframes with metadata
# metaMDS points have sample names as rownames — rownames_to_column is correct here.
# full_meta_16S already has a "sample_name" column, so left_join matches correctly.
df.nmds.bray.16S <- as.data.frame(ds.nmds.bray.16S[["points"]]) %>%
  rownames_to_column("sample_name") %>%
  left_join(full_meta_16S, by = "sample_name")

df.nmds.jacc.16S <- as.data.frame(ds.nmds.jacc.16S[["points"]]) %>%
  rownames_to_column("sample_name") %>%
  left_join(full_meta_16S, by = "sample_name")

# Build PERMANOVA annotation strings from adonis2 results
label_bc_16S   <- make_permanova_label(full_bc_permanova_16S,   "Bray-Curtis")
label_jacc_16S <- make_permanova_label(full_jacc_permanova_16S, "Jaccard")

# ---- 16S Bray-Curtis NMDS ----
bc.nmds.16S <- df.nmds.bray.16S %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = Clade_Order), linewidth = 1, linetype = "longdash") +
  geom_point(aes(color = Clade_Order, shape = env_medium), size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = Clade_Order),
               size = 7, shape = 21) +
  annotate(geom = "richtext",
           x = max(df.nmds.bray.16S$MDS1), y = min(df.nmds.bray.16S$MDS2),
           hjust = 1, vjust = 0.5, size = 2.75,
           label = label_bc_16S, fill = NA, label.color = NA) +
  annotate(geom = "text",
           x = min(df.nmds.bray.16S$MDS1), y = min(df.nmds.bray.16S$MDS2),
           hjust = 0, vjust = 1, size = 3,
           label = paste0("Stress = ", round(ds.nmds.bray.16S$stress, 3))) +
  labs(title = "Bacteria", subtitle = "Bray-Curtis NMDS",
       fill = "Host Order", color = "Host Order", shape = "Sample type") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10)) +
  guides(shape = guide_legend(order = 2),
         fill  = guide_legend(order = 1),
         color = guide_legend(order = 1)) +
  scale_color_manual(values = clade_colors) +
  scale_fill_manual(values  = clade_colors)

# ---- 16S Jaccard NMDS ----
jacc.nmds.16S <- df.nmds.jacc.16S %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = Clade_Order), linewidth = 1, linetype = "longdash") +
  geom_point(aes(color = Clade_Order, shape = env_medium), size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = Clade_Order),
               size = 7, shape = 21) +
  annotate(geom = "richtext",
           x = max(df.nmds.jacc.16S$MDS1), y = min(df.nmds.jacc.16S$MDS2),
           hjust = 1, vjust = 0.2, size = 2.75,
           label = label_jacc_16S, fill = NA, label.color = NA) +
  annotate(geom = "text",
           x = min(df.nmds.jacc.16S$MDS1), y = min(df.nmds.jacc.16S$MDS2),
           hjust = 0, vjust = 1, size = 3,
           label = paste0("Stress = ", round(ds.nmds.jacc.16S$stress, 3))) +
  labs(title = "Bacteria", subtitle = "Jaccard NMDS",
       fill = "Host Order", color = "Host Order", shape = "Sample type") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10)) +
  guides(shape = guide_legend(order = 2),
         fill  = guide_legend(order = 1),
         color = guide_legend(order = 1)) +
  scale_color_manual(values = clade_colors) +
  scale_fill_manual(values  = clade_colors)

plot(bc.nmds.16S)
plot(jacc.nmds.16S)

# Save 16S NMDS plots
#ggsave(file.path(dir_figures, "16S_bray_nmds_full.pdf"), bc.nmds.16S, width = 7, height = 6, units = "in", device = cairo_pdf)
#ggsave(file.path(dir_figures, "16S_jaccard_nmds_full.pdf"), jacc.nmds.16S, width = 7, height = 6, units = "in", device = cairo_pdf)

###############################################################################
# ========================== ITS DATASET ======================================
###############################################################################

# ========================== LOAD DATA =========================================
ps_ITS <- readRDS(file.path(dir_processed, "ITS_absolute_all_final.rds"))

# ========================== OUTLIER REMOVAL (ITS) =============================
# UHM1483-34786 identified as an outlier in the ITS ordination.
# Removed before all analyses for consistency.
outlier_ITS <- "UHM1483-34786"
ps_ITS      <- prune_samples(sample_names(ps_ITS) != outlier_ITS, ps_ITS)
cat("ITS: Removed outlier sample:", outlier_ITS, "\n")
cat("ITS: Samples remaining:", nsamples(ps_ITS), "\n\n")

# ========================== HABITAT CLEANUP (ITS) =============================
full_meta_ITS <- data.frame(sample_data(ps_ITS))

# sample_data already contains a "sample_name" column — no rownames manipulation needed.

full_meta_ITS$animal_ecomode <- as.character(full_meta_ITS$animal_ecomode)

if ("Fossorial" %in% full_meta_ITS$animal_ecomode) {
  n_fossorial <- sum(full_meta_ITS$animal_ecomode == "Fossorial")
  cat("ITS: Collapsing", n_fossorial,
      "'Fossorial' samples into 'Fossorial-Terrestrial'\n")
  full_meta_ITS$animal_ecomode[full_meta_ITS$animal_ecomode == "Fossorial"] <-
    "Fossorial-Terrestrial"
} else {
  cat("ITS: No 'Fossorial' entries found — no collapse needed\n")
}
full_meta_ITS$animal_ecomode <- factor(full_meta_ITS$animal_ecomode)

# Propagate corrected ecomode back into the phyloseq sample_data
# sample_name is already a column AND the rownames — assign rownames directly
meta_for_ps_ITS <- full_meta_ITS
rownames(meta_for_ps_ITS) <- meta_for_ps_ITS$sample_name
sample_data(ps_ITS) <- sample_data(meta_for_ps_ITS)

# ========================== DISTANCE MATRICES (ITS) ==========================
full_ds_bray_ITS <- phyloseq::distance(ps_ITS, method = "bray")
full_ds_jacc_ITS <- phyloseq::distance(ps_ITS, method = "jaccard", binary = TRUE)

cat("ITS distance matrices computed:",
    nrow(as.matrix(full_ds_bray_ITS)), "samples\n\n")

##------------------------------------------------------------------------------
## Beta-dispersion (ITS) ####
##------------------------------------------------------------------------------
betadisper_dir_ITS <- file.path(dir_figures, "betadisper_plots", "ITS")
dir.create(betadisper_dir_ITS, recursive = TRUE, showWarnings = FALSE)
dir.create(betadisper_dir_ITS, recursive = TRUE, showWarnings = FALSE)

# Make sure R Studio's "Plots" tab is large enough for figure viewing

# Bray-Curtis
bd_ITS_bc_order   <- run_betadisper(full_ds_bray_ITS, full_meta_ITS$Clade_Order,     "Host Order",  "Bray-Curtis", "ITS", betadisper_dir_ITS)
bd_ITS_bc_wc      <- run_betadisper(full_ds_bray_ITS, full_meta_ITS$env_broad_scale, "Management",  "Bray-Curtis", "ITS", betadisper_dir_ITS)
bd_ITS_bc_diet    <- run_betadisper(full_ds_bray_ITS, full_meta_ITS$Diet,             "Diet",        "Bray-Curtis", "ITS", betadisper_dir_ITS)
bd_ITS_bc_sample  <- run_betadisper(full_ds_bray_ITS, full_meta_ITS$env_medium,       "Sample Type", "Bray-Curtis", "ITS", betadisper_dir_ITS)
bd_ITS_bc_habitat <- run_betadisper(full_ds_bray_ITS, full_meta_ITS$animal_ecomode,   "Habitat",     "Bray-Curtis", "ITS", betadisper_dir_ITS)

# Jaccard
bd_ITS_jc_order   <- run_betadisper(full_ds_jacc_ITS, full_meta_ITS$Clade_Order,     "Host Order",  "Jaccard", "ITS", betadisper_dir_ITS)
bd_ITS_jc_wc      <- run_betadisper(full_ds_jacc_ITS, full_meta_ITS$env_broad_scale, "Management",  "Jaccard", "ITS", betadisper_dir_ITS)
bd_ITS_jc_diet    <- run_betadisper(full_ds_jacc_ITS, full_meta_ITS$Diet,             "Diet",        "Jaccard", "ITS", betadisper_dir_ITS)
bd_ITS_jc_sample  <- run_betadisper(full_ds_jacc_ITS, full_meta_ITS$env_medium,       "Sample Type", "Jaccard", "ITS", betadisper_dir_ITS)
bd_ITS_jc_habitat <- run_betadisper(full_ds_jacc_ITS, full_meta_ITS$animal_ecomode,   "Habitat",     "Jaccard", "ITS", betadisper_dir_ITS)

##------------------------------------------------------------------------------
## PERMANOVA (ITS) ####
##------------------------------------------------------------------------------
# Bray-Curtis
full_bc_permanova_ITS <- adonis2(
  full_ds_bray_ITS ~ Clade_Order * Diet * env_broad_scale * env_medium,
  data         = full_meta_ITS,
  strata       = full_meta_ITS$site,
  permutations = 999,
  by           = "terms"
)

# Jaccard
full_jacc_permanova_ITS <- adonis2(
  full_ds_jacc_ITS ~ Clade_Order * Diet * env_broad_scale * env_medium,
  data         = full_meta_ITS,
  strata       = full_meta_ITS$site,
  permutations = 999,
  by           = "terms"
)

# Save results
#saveRDS(full_bc_permanova_ITS,   file = file.path(dir_processed, "ITS_full_bc_permanova.rds"))
#saveRDS(full_jacc_permanova_ITS, file = file.path(dir_processed, "ITS_full_jacc_permanova.rds"))

cat("=== ITS PERMANOVA (Bray-Curtis) ===\n"); print(full_bc_permanova_ITS)
cat("=== ITS PERMANOVA (Jaccard) ===\n");     print(full_jacc_permanova_ITS)

##------------------------------------------------------------------------------
## NMDS Ordinations (ITS) ####
##------------------------------------------------------------------------------
ds.nmds.bray.ITS <- metaMDS(full_ds_bray_ITS, k = 5)
ds.nmds.jacc.ITS <- metaMDS(full_ds_jacc_ITS, k = 5)

cat("ITS Bray-Curtis NMDS stress:", round(ds.nmds.bray.ITS$stress, 4), "\n")
cat("ITS Jaccard NMDS stress:",     round(ds.nmds.jacc.ITS$stress, 4), "\n\n")

# Save NMDS results
#saveRDS(ds.nmds.bray.ITS, file = file.path(dir_processed, "ITS_nmds_bray.rds"))
#saveRDS(ds.nmds.jacc.ITS, file = file.path(dir_processed, "ITS_nmds_jacc.rds"))

# Build NMDS dataframes with metadata
# metaMDS points have sample names as rownames — rownames_to_column is correct here.
# full_meta_ITS already has a "sample_name" column, so left_join matches correctly.
df.nmds.bray.ITS <- as.data.frame(ds.nmds.bray.ITS[["points"]]) %>%
  rownames_to_column("sample_name") %>%
  left_join(full_meta_ITS, by = "sample_name")

df.nmds.jacc.ITS <- as.data.frame(ds.nmds.jacc.ITS[["points"]]) %>%
  rownames_to_column("sample_name") %>%
  left_join(full_meta_ITS, by = "sample_name")

# Build PERMANOVA annotation strings
label_bc_ITS   <- make_permanova_label(full_bc_permanova_ITS,   "Bray-Curtis")
label_jacc_ITS <- make_permanova_label(full_jacc_permanova_ITS, "Jaccard")

# ---- ITS Bray-Curtis NMDS ----
bc.nmds.ITS <- df.nmds.bray.ITS %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = Clade_Order), linewidth = 1, linetype = "longdash") +
  geom_point(aes(color = Clade_Order, shape = env_medium), size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = Clade_Order),
               size = 7, shape = 21) +
  annotate(geom = "richtext",
           x = max(df.nmds.bray.ITS$MDS1), y = min(df.nmds.bray.ITS$MDS2),
           hjust = 1, vjust = 0.75, size = 2.75,
           label = label_bc_ITS, fill = NA, label.color = NA) +
  annotate(geom = "text",
           x = min(df.nmds.bray.ITS$MDS1), y = min(df.nmds.bray.ITS$MDS2),
           hjust = 0, vjust = 1, size = 3,
           label = paste0("Stress = ", round(ds.nmds.bray.ITS$stress, 3))) +
  labs(title = "Fungi", subtitle = "Bray-Curtis NMDS",
       fill = "Host Order", color = "Host Order", shape = "Sample type") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10)) +
  guides(shape = guide_legend(order = 2),
         fill  = guide_legend(order = 1),
         color = guide_legend(order = 1)) +
  scale_color_manual(values = clade_colors) +
  scale_fill_manual(values  = clade_colors)

# ---- ITS Jaccard NMDS ----
jacc.nmds.ITS <- df.nmds.jacc.ITS %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = Clade_Order), linewidth = 1, linetype = "longdash") +
  geom_point(aes(color = Clade_Order, shape = env_medium), size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = Clade_Order),
               size = 7, shape = 21) +
  annotate(geom = "richtext",
           x = max(df.nmds.jacc.ITS$MDS1), y = min(df.nmds.jacc.ITS$MDS2),
           hjust = 1, vjust = 1, size = 2.75,
           label = label_jacc_ITS, fill = NA, label.color = NA) +
  annotate(geom = "text",
           x = min(df.nmds.jacc.ITS$MDS1), y = min(df.nmds.jacc.ITS$MDS2),
           hjust = 0, vjust = 1, size = 3,
           label = paste0("Stress = ", round(ds.nmds.jacc.ITS$stress, 3))) +
  labs(title = "Fungi", subtitle = "Jaccard NMDS",
       fill = "Host Order", color = "Host Order", shape = "Sample type") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10)) +
  guides(shape = guide_legend(order = 2),
         fill  = guide_legend(order = 1),
         color = guide_legend(order = 1)) +
  scale_color_manual(values = clade_colors) +
  scale_fill_manual(values  = clade_colors)

plot(bc.nmds.ITS)
plot(jacc.nmds.ITS)

# Save ITS NMDS plots
#ggsave(file.path(dir_figures, "ITS_bray_nmds_full.pdf"), bc.nmds.ITS, width = 7, height = 6, units = "in", device = cairo_pdf)
#ggsave(file.path(dir_figures, "ITS_jaccard_nmds_full.pdf"), jacc.nmds.ITS, width = 7, height = 6, units = "in", device = cairo_pdf)

###############################################################################
# ========================== MULTIPANEL FIGURES ================================
###############################################################################
# Bray-Curtis: Bacteria (left) | Fungi (right) — shared right-side legend
# Jaccard:     Bacteria (left) | Fungi (right) — shared right-side legend

multi_bray <- (bc.nmds.16S + bc.nmds.ITS) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(multi_bray)
ggsave(file.path(dir_figures, "nmds_bray_multipanel.pdf"), multi_bray,
       width = 14, height = 6, units = "in", device = cairo_pdf, dpi = 750)

multi_jacc <- (jacc.nmds.16S + jacc.nmds.ITS) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(multi_jacc)
ggsave(file.path(dir_figures, "nmds_jaccard_multipanel.pdf"), multi_jacc,
       width = 14, height = 6, units = "in", device = cairo_pdf, dpi = 750)

# Stacked Jaccard: Bacteria (top) / Fungi (bottom) — shared right-side legend
multi_jacc_stacked <- (jacc.nmds.16S / jacc.nmds.ITS) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(multi_jacc_stacked)
ggsave(file.path(dir_figures, "nmds_jaccard_multipanel_stacked.pdf"), multi_jacc_stacked,
       width = 7, height = 12, units = "in", device = cairo_pdf, dpi = 750)

###############################################################################
# ========================== EXPORT PERMANOVA RESULTS TABLE ===================
###############################################################################
# Combines all four PERMANOVA results (16S Bray, 16S Jaccard, ITS Bray,
# ITS Jaccard) into a single CSV matching the supplementary table format.
# Significance codes: *** <0.001, ** <0.01, * <0.05, . <0.1, ns >=0.1

# Build combined table
permanova_table <- bind_rows(
  format_permanova(full_bc_permanova_16S,   "bacteria_16S", "bray-curtis"),
  format_permanova(full_jacc_permanova_16S, "bacteria_16S", "jaccard"),
  format_permanova(full_bc_permanova_ITS,   "fungi_ITS",    "bray-curtis"),
  format_permanova(full_jacc_permanova_ITS, "fungi_ITS",    "jaccard")
)

write.csv(permanova_table, file = file.path(dir_tables, "PERMANOVA_results_full.csv"), row.names = FALSE)

cat("PERMANOVA results table saved:",
    file.path(dir_tables, "PERMANOVA_results_full.csv"), "\n")
cat("Rows:", nrow(permanova_table), "\n")

###############################################################################
###############################################################################
###############################################################################
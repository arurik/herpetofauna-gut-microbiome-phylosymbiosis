###############################################################################
# Basidiobolus sensitivity analysis — NMDS + PERMANOVA
# Distance metrics:
#
# Goal: Test whether Basidiobolus drives the dominant ordination separation
#       between Anura and Caudata, and test whether the signal persists
#       under presence/absence (Jaccard) — isolating abundance effects
#       from occupancy effects.
# Outlier removed: UHM1483-34786 (Telmatobius culeus)
#   MDS1 = 910.9 in Bray-Curtis NMDS — extreme outlier collapsing ordination
#   space. Aquatic obligate specialist; removed prior to all analyses.
#
# Alexander Rurik
###############################################################################

# ========== Load packages ===================================================
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(ggtext)
})

select <- dplyr::select
filter <- dplyr::filter

set.seed(52325)

# ========== Project directories =============================================
library(here)
dir_processed <- here("data", "processed")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")

# ===========================================================================
# 1) Load data
# ===========================================================================
ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_fecal.rds"))

# ===========================================================================
# 2) Subset to amphibians only and remove outlier
# ===========================================================================
ps_amph <- subset_samples(ps_ITS, Clade_Order %in% c("Anura", "Caudata"))

# Remove outlier: UHM1483-34786 (Telmatobius culeus)
# Identified as extreme NMDS outlier (MDS1 = 910.9)
ps_amph <- prune_samples(sample_names(ps_amph) != "UHM1483-34786", ps_amph)
ps_amph <- prune_taxa(taxa_sums(ps_amph) > 0, ps_amph)

cat("Amphibian-only phyloseq object (outlier removed):\n")
print(ps_amph)
cat("Sample counts by order:\n")
print(table(sample_data(ps_amph)$Clade_Order))

# ===========================================================================
# 3) Extract OTU matrix and identify Basidiobolus OTUs
# ===========================================================================
tax <- as.data.frame(tax_table(ps_amph), stringsAsFactors = FALSE)
otu <- as.data.frame(otu_table(ps_amph))
if (taxa_are_rows(ps_amph)) otu <- t(otu)
# otu is now samples x OTUs

basid_otus <- rownames(tax)[tax$Genus == "Basidiobolus"]
cat("\nNumber of OTUs assigned to Basidiobolus:", length(basid_otus), "\n")

if (length(basid_otus) == 0) {
  stop("No OTUs found with Genus == 'Basidiobolus'. Check taxonomy column name/spelling.")
}

# Basidiobolus-removed OTU matrix
otu_no_basid <- otu[, !colnames(otu) %in% basid_otus, drop = FALSE]
cat("OTUs removed:", length(basid_otus), "\n")
cat("OTU matrix dimensions after removal:", dim(otu_no_basid), "\n")

# Presence/absence matrix for Jaccard
otu_pa       <- decostand(otu,          method = "pa")
otu_pa_nb    <- decostand(otu_no_basid, method = "pa")

# ===========================================================================
# 4) Extract metadata
# ===========================================================================
env <- data.frame(sample_data(ps_amph), check.names = FALSE)
env$SampleID <- rownames(env)
stopifnot(all(rownames(env) == rownames(otu)))


# ===========================================================================
# 5) Shared aesthetics
# ===========================================================================
order_colors <- c(
  "Caudata" = "#E69F00",
  "Anura"   = "#CC79A7"
)

order_shapes <- c(
  "Anura"   = 16,
  "Caudata" = 17
)

# ===========================================================================
# 6) Helper: PERMANOVA annotation label (HTML richtext for ggtext)
# ===========================================================================
make_permanova_label <- function(permanova_obj, distance_name) {
  res      <- as.data.frame(permanova_obj)
  res$term <- rownames(res)
  names(res)[names(res) == "Pr(>F)"] <- "pval"
  names(res)[names(res) == "F"]      <- "Fstat"
  
  row      <- res %>% dplyr::filter(term == "Clade_Order")
  r2_fmt   <- sprintf("%.3f", row$R2)
  pval     <- row$pval
  pval_fmt <- ifelse(pval < 0.001, "< 0.001", sprintf("= %.3f", pval))
  
  paste0(
    "<b>PERMANOVA (", distance_name, ")</b><br>",
    "Host Order: R<sup>2</sup> = ", r2_fmt, ", p ", pval_fmt
  )
}

# ===========================================================================
# 7) Helper: build standard NMDS plot
# ===========================================================================
build_nmds_plot <- function(df,
                            perm_label,
                            stress_val,
                            title_str,
                            sub_str,
                            suppress_legend = FALSE) {
  p <- df %>%
    ggplot(aes(x = MDS1, y = MDS2)) +
    stat_ellipse(
      aes(color = Clade_Order),
      linewidth = 1,
      linetype  = "longdash"
    ) +
    geom_point(
      aes(color = Clade_Order, shape = Clade_Order),
      size  = 2.5,
      alpha = 0.75
    ) +
    stat_ellipse(
      level = 1e-10,
      geom  = "point",
      aes(fill = Clade_Order),
      size  = 7,
      shape = 21
    ) +
    annotate(
      geom        = "richtext",
      x           = max(df$MDS1, na.rm = TRUE),
      y           = min(df$MDS2, na.rm = TRUE),
      hjust       = 1,
      vjust       = 0.75,
      size        = 2.75,
      label       = perm_label,
      fill        = NA,
      label.color = NA
    ) +
    annotate(
      geom  = "text",
      x     = min(df$MDS1, na.rm = TRUE),
      y     = min(df$MDS2, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size  = 3,
      label = paste0("Stress = ", round(stress_val, 3))
    ) +
    scale_color_manual(values = order_colors, name = "Host Order") +
    scale_fill_manual(values  = order_colors, name = "Host Order") +
    scale_shape_manual(values = order_shapes, name = "Host Order") +
    labs(title = title_str, subtitle = sub_str) +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold", size = 10),
      plot.subtitle   = element_text(size = 9)
    ) +
    guides(
      fill  = guide_legend(order = 1),
      color = guide_legend(order = 1),
      shape = guide_legend(order = 1)
    )
  
  if (suppress_legend) p <- p + theme(legend.position = "none")
  return(p)
}

# ===========================================================================
# 8) Distance matrices
# ===========================================================================
cat("\n=== Computing distance matrices ===\n")

# Bray-Curtis (raw absolute abundance)
dist_bc     <- vegdist(otu,          method = "bray")
dist_bc_nb  <- vegdist(otu_no_basid, method = "bray")

# Jaccard (presence/absence)
dist_jac    <- vegdist(otu_pa,    method = "jaccard", binary = TRUE)
dist_jac_nb <- vegdist(otu_pa_nb, method = "jaccard", binary = TRUE)

# ===========================================================================
# 9) PERMANOVA — both distance matrices
#     adonis2, strata = site, 9999 permutations
# ===========================================================================
cat("\n=== PERMANOVA (9999 permutations, strata = site) ===\n")

env_fac <- env %>% mutate(Clade_Order = as.factor(Clade_Order))

# Bray-Curtis
perm_bc_full <- adonis2(dist_bc    ~ Clade_Order, data = env_fac,
                        strata = env$site, permutations = 9999, by = "terms")
perm_bc_nb   <- adonis2(dist_bc_nb ~ Clade_Order, data = env_fac,
                        strata = env$site, permutations = 9999, by = "terms")

cat("\nBray-Curtis — full community:\n");    print(perm_bc_full)
cat("\nBray-Curtis — Basidiobolus removed:\n"); print(perm_bc_nb)

# Jaccard
perm_jac_full <- adonis2(dist_jac    ~ Clade_Order, data = env_fac,
                         strata = env$site, permutations = 9999, by = "terms")
perm_jac_nb   <- adonis2(dist_jac_nb ~ Clade_Order, data = env_fac,
                         strata = env$site, permutations = 9999, by = "terms")

cat("\nJaccard — full community:\n");    print(perm_jac_full)
cat("\nJaccard — Basidiobolus removed:\n"); print(perm_jac_nb)

# Combined summary table
r2_bc_full  <- perm_bc_full["Clade_Order",  "R2"]
r2_bc_nb    <- perm_bc_nb["Clade_Order",    "R2"]
r2_jac_full <- perm_jac_full["Clade_Order", "R2"]
r2_jac_nb   <- perm_jac_nb["Clade_Order",   "R2"]

p_bc_full   <- perm_bc_full["Clade_Order",  "Pr(>F)"]
p_bc_nb     <- perm_bc_nb["Clade_Order",    "Pr(>F)"]
p_jac_full  <- perm_jac_full["Clade_Order", "Pr(>F)"]
p_jac_nb    <- perm_jac_nb["Clade_Order",   "Pr(>F)"]

perm_summary <- data.frame(
  Distance     = c(rep("Bray-Curtis (raw absolute abundance)", 2),
                   rep("Jaccard (presence/absence)", 2)),
  Basidiobolus = rep(c("Included", "Removed"), 2),
  R2           = round(c(r2_bc_full, r2_bc_nb, r2_jac_full, r2_jac_nb), 5),
  p            = c(p_bc_full, p_bc_nb, p_jac_full, p_jac_nb),
  R2_reduction_pct = c(
    NA, round((1 - r2_bc_nb  / r2_bc_full)  * 100, 1),
    NA, round((1 - r2_jac_nb / r2_jac_full) * 100, 1)
  )
)

cat("\n=== PERMANOVA combined summary ===\n")
print(perm_summary)

write.csv(perm_summary,
          file      = file.path(dir_tables, "ITS_amph_PERMANOVA_summary.csv"),
          row.names = FALSE)
cat("PERMANOVA summary CSV saved.\n")

# Annotation labels for plots
label_bc_full  <- make_permanova_label(perm_bc_full,  "Bray-Curtis")
label_bc_nb    <- make_permanova_label(perm_bc_nb,    "Bray-Curtis")
label_jac_full <- make_permanova_label(perm_jac_full, "Jaccard")
label_jac_nb   <- make_permanova_label(perm_jac_nb,   "Jaccard")

# ===========================================================================
# 10) NMDS
# ===========================================================================
cat("\n=== Running NMDS ===\n")

# Bray-Curtis NMDS
cat("Bray-Curtis NMDS...\n")
nmds_bc_full <- metaMDS(dist_bc,    k = 4, trymax = 100,
                        autotransform = FALSE, trace = FALSE)
nmds_bc_nb   <- metaMDS(dist_bc_nb, k = 4, trymax = 100,
                        autotransform = FALSE, trace = FALSE)

# Jaccard NMDS
cat("Jaccard NMDS...\n")
nmds_jac_full <- metaMDS(dist_jac,    k = 4, trymax = 100,
                         autotransform = FALSE, trace = FALSE)
nmds_jac_nb   <- metaMDS(dist_jac_nb, k = 4, trymax = 100,
                         autotransform = FALSE, trace = FALSE)

cat("\nStress values:\n")
cat("  Bray-Curtis full:       ", round(nmds_bc_full$stress,  4), "\n")
cat("  Bray-Curtis no Basid:   ", round(nmds_bc_nb$stress,    4), "\n")
cat("  Jaccard full:           ", round(nmds_jac_full$stress, 4), "\n")
cat("  Jaccard no Basid:       ", round(nmds_jac_nb$stress,   4), "\n")

df_bc_full <- as.data.frame(nmds_bc_full$points) %>%
  rownames_to_column("SampleID") %>%
  left_join(env, by = "SampleID")

df_bc_nb <- as.data.frame(nmds_bc_nb$points) %>%
  rownames_to_column("SampleID") %>%
  left_join(env, by = "SampleID")

df_jac_full <- as.data.frame(nmds_jac_full$points) %>%
  rownames_to_column("SampleID") %>%
  left_join(env, by = "SampleID")

df_jac_nb <- as.data.frame(nmds_jac_nb$points) %>%
  rownames_to_column("SampleID") %>%
  left_join(env, by = "SampleID")

# ===========================================================================
# 11) Build plots
# ===========================================================================

# ---- Bray-Curtis panels ----
plot_A <- build_nmds_plot(
  df         = df_bc_full,
  perm_label = label_bc_full,
  stress_val = nmds_bc_full$stress,
  title_str  = "Fungi (Amphibians)",
  sub_str    = "Bray-Curtis NMDS \u2014 Basidiobolus included"
)
plot_A

plot_B <- build_nmds_plot(
  df              = df_bc_nb,
  perm_label      = label_bc_nb,
  stress_val      = nmds_bc_nb$stress,
  title_str       = "Fungi (Amphibians)",
  sub_str         = "Bray-Curtis NMDS \u2014 Basidiobolus removed",
  suppress_legend = TRUE
)
plot_B

# ---- Jaccard panels ----
plot_D <- build_nmds_plot(
  df         = df_jac_full,
  perm_label = label_jac_full,
  stress_val = nmds_jac_full$stress,
  title_str  = "Fungi (Amphibians)",
  sub_str    = "Jaccard NMDS (axes 1\u20132) \u2014 Basidiobolus included"
)
plot_D

plot_E <- build_nmds_plot(
  df              = df_jac_nb,
  perm_label      = label_jac_nb,
  stress_val      = nmds_jac_nb$stress,
  title_str       = "Fungi (Amphibians)",
  sub_str         = "Jaccard NMDS (axes 1\u20132) \u2014 Basidiobolus removed",
  suppress_legend = TRUE
)
plot_E

# ===========================================================================
# 12) Combined panel figures
# ===========================================================================

# Bray-Curtis two-panel (manuscript candidate)
twopanel_bc <- (plot_A + plot_B) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 11, face = "bold"))
  )
twopanel_bc

# Jaccard two-panel
twopanel_jac <- (plot_D + plot_E) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 11, face = "bold"))
  )
twopanel_jac

# ===========================================================================
# 13) Save figures (PDF only)
# ===========================================================================

# Bray-Curtis
ggsave(file.path(dir_figures, "ITS_amph_bc_nmds_full.pdf"),
       plot = plot_A, width = 7, height = 6, units = "in", device = cairo_pdf)
ggsave(file.path(dir_figures, "ITS_amph_bc_nmds_sensitivity.pdf"),
       plot = plot_B, width = 7, height = 6, units = "in", device = cairo_pdf)
ggsave(file.path(dir_figures, "ITS_amph_bc_nmds_twopanel.pdf"),
       plot = twopanel_bc, width = 13, height = 6,
       units = "in", device = cairo_pdf)

# Jaccard
ggsave(file.path(dir_figures, "ITS_amph_jac_nmds_full.pdf"),
       plot = plot_D, width = 7, height = 6, units = "in", device = cairo_pdf)
ggsave(file.path(dir_figures, "ITS_amph_jac_nmds_sensitivity.pdf"),
       plot = plot_E, width = 7, height = 6, units = "in", device = cairo_pdf)
ggsave(file.path(dir_figures, "ITS_amph_jac_nmds_twopanel.pdf"),
       plot = twopanel_jac, width = 13, height = 6,
       units = "in", device = cairo_pdf)

cat("\nAll figures saved.\n")
cat("\nScript complete.\n")


###############################################################################
###############################################################################
###############################################################################
################################################################################
# Core Microbiome Analysis — Wild Samples Only
#
# Goal: describe the core gut bacterial and fungal microbiomes for reptiles & amphibians
# Outputs:
#   CSV  : <date>_core_microbiome_WILD.csv                      — Amphibia vs. Reptilia (bacteria + fungi)
#   CSV  : <date>_core_microbiome_amphibian_orders_WILD.csv     — Caudata vs. Anura (bacteria + fungi)
#   CSV  : <date>_core_microbiome_squamata_WILD.csv             — Lizards vs. Snakes (bacteria + fungi)
#   PDF  : <date>_core_microbiome_GENUS_WILD_fungi.pdf          — 4-panel (fungi lollipop)
#   PDF  : <date>_core_microbiome_GENUS_WILD_bacteria.pdf       — 4-panel (bacteria lollipop)
#   PDF  : <date>_core_bacteria_lollipop_amphibian_orders.pdf   — Caudata vs. Anura lollipop
#   PDF  : <date>_core_bacteria_lollipop_squamata.pdf           — Lizards vs. Snakes lollipop
#
# "Common core" threshold decision: ≥50% prevalence across individual samples within each group
#
# Alexander Rurik 
################################################################################

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
  library(VennDiagram)
  library(gridExtra)
  library(grid)
  library(cowplot)
})

select <- dplyr::select
filter <- dplyr::filter

# ========================== BASIC PARAMETERS ==================================
set.seed(52325)
PREV_THRESHOLD <- 0.5 # Chosen core prevalence threshold (proportion of samples)

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")

# ========================== LOAD DATA =========================================
ps_16S <- readRDS(file.path(dir_processed, "16S_abs_final_fecal.rds"))
ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_fecal.rds"))

# ========================== COLOR PALETTES ====================================
clade_colors <- c(
  "Amphibia" = "#0072B2",
  "Reptilia" = "#009E73"
)

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
# ==================== PART 1: PREPARE PHYLOSEQ OBJECTS =======================
###############################################################################
cat("\n====== CORE MICROBIOME ANALYSIS (GENUS LEVEL, WILD ONLY) ======\n")
cat("Core threshold: ≥", PREV_THRESHOLD * 100, "% prevalence across individual samples\n\n")

# ---- Subset to wild samples and add Clade column ----
add_clade_column <- function(ps) {
  sdata <- data.frame(sample_data(ps))
  sdata$Clade <- case_when(
    sdata$Clade_Order %in% amphibian_orders ~ "Amphibia",
    sdata$Clade_Order %in% reptile_orders   ~ "Reptilia",
    TRUE                                    ~ NA_character_
  )
  sample_data(ps) <- sdata
  return(ps)
}

ps_16S_wild <- subset_samples(ps_16S, env_broad_scale == "Wild") |> add_clade_column()
ps_ITS_wild <- subset_samples(ps_ITS, env_broad_scale == "Wild") |> add_clade_column()

# ---- Agglomerate to genus level ----
# NArm = FALSE retains taxa assigned at higher ranks but not to genus
cat("Agglomerating to genus level...\n")
ps_16S_genus <- tax_glom(ps_16S_wild, taxrank = "Genus", NArm = FALSE)
ps_ITS_genus <- tax_glom(ps_ITS_wild, taxrank = "Genus", NArm = FALSE)

ps_16S_genus <- prune_taxa(taxa_sums(ps_16S_genus) > 0, ps_16S_genus)
ps_ITS_genus <- prune_taxa(taxa_sums(ps_ITS_genus) > 0, ps_ITS_genus)

cat("16S:", nsamples(ps_16S_genus), "samples |", ntaxa(ps_16S_genus), "genera\n")
cat("ITS:", nsamples(ps_ITS_genus), "samples |", ntaxa(ps_ITS_genus), "genera\n")

# ---- Filter unassigned / incertae sedis taxa ----
# Removes: (1) Kingdom-level unassigned; (2) Fungi_phy_incertae_sedis genera
# Retains:  Rozellomycota_gen_incertae_sedis (has valid phylum-level assignment)
filter_unassigned <- function(ps) {
  tax <- as.data.frame(tax_table(ps))
  remove_idx <- (is.na(tax$Kingdom) | tax$Kingdom == "Unassigned") |
    (tax$Phylum == "Fungi_phy_incertae_sedis" & grepl("incertae_sedis", tax$Genus))
  taxa_to_remove <- rownames(tax)[remove_idx]
  ps_out <- prune_taxa(setdiff(rownames(tax), taxa_to_remove), ps)
  cat("  Removed", length(taxa_to_remove), "unassigned taxa;",
      ntaxa(ps_out), "genera retained\n")
  return(ps_out)
}

cat("\nFiltering unassigned taxa...\n")
ps_16S_genus <- filter_unassigned(ps_16S_genus)
ps_ITS_genus <- filter_unassigned(ps_ITS_genus)

###############################################################################
# ========================== PART 2: IDENTIFY CORE TAXA =======================
###############################################################################

# Returns OTU names present in >= prevalence_threshold proportion of samples
get_core_taxa <- function(ps, prevalence_threshold = PREV_THRESHOLD, min_abundance = 0) {
  otu <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) otu <- t(otu)
  prevalence <- rowSums(otu > min_abundance) / ncol(otu)
  return(names(prevalence[prevalence >= prevalence_threshold]))
}

# Subsets phyloseq to a group, then calls get_core_taxa
get_core_by_group <- function(ps, group_var, group_value,
                              prevalence_threshold = PREV_THRESHOLD) {
  sdata        <- data.frame(sample_data(ps))
  samples_keep <- rownames(sdata)[sdata[[group_var]] == group_value]
  if (length(samples_keep) == 0) {
    warning("No samples found for ", group_var, " = ", group_value)
    return(character(0))
  }
  cat("  ", group_value, ": n =", length(samples_keep), "samples\n")
  ps_sub <- prune_samples(samples_keep, ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(get_core_taxa(ps_sub, prevalence_threshold))
}

# ---- Overall core (all wild samples) ----
core_16S_all <- get_core_taxa(ps_16S_genus)
core_ITS_all <- get_core_taxa(ps_ITS_genus)

cat("\nCore genera (≥50% prevalence, all wild samples):\n")
cat("  Bacteria:", length(core_16S_all), "genera\n")
cat("  Fungi:   ", length(core_ITS_all), "genera\n")

# ---- Core by clade ----
cat("\nCalculating core by clade...\n")
core_16S_amph <- get_core_by_group(ps_16S_genus, "Clade", "Amphibia")
core_16S_rept <- get_core_by_group(ps_16S_genus, "Clade", "Reptilia")
core_ITS_amph <- get_core_by_group(ps_ITS_genus, "Clade", "Amphibia")
core_ITS_rept <- get_core_by_group(ps_ITS_genus, "Clade", "Reptilia")

cat("\nCore genera by clade:\n")
cat("Bacteria — Amphibia:", length(core_16S_amph),
    "| Reptilia:", length(core_16S_rept),
    "| Shared:", length(intersect(core_16S_amph, core_16S_rept)), "\n")
cat("Fungi    — Amphibia:", length(core_ITS_amph),
    "| Reptilia:", length(core_ITS_rept),
    "| Shared:", length(intersect(core_ITS_amph, core_ITS_rept)), "\n")

# ---- Clade-specific core (not shared between clades) ----
core_16S_amph_only <- setdiff(core_16S_amph, core_16S_rept)
core_16S_rept_only <- setdiff(core_16S_rept, core_16S_amph)
core_ITS_amph_only <- setdiff(core_ITS_amph, core_ITS_rept)
core_ITS_rept_only <- setdiff(core_ITS_rept, core_ITS_amph)

cat("\nClade-specific core (not shared):\n")
cat("Bacteria — Amphibia only:", length(core_16S_amph_only),
    "| Reptilia only:", length(core_16S_rept_only), "\n")
cat("Fungi    — Amphibia only:", length(core_ITS_amph_only),
    "| Reptilia only:", length(core_ITS_rept_only), "\n")

###############################################################################
# ========================== PART 3: BUILD SUPPLEMENTARY TABLE ================
#
# Matches format of reference table:
#   Kingdom | Phylum | Class | Order | Family | Genus |
#   Sample_prevalence | Species_prevalence |
#   mean_absolute_abundance | median_absolute_abundance | total_absolute_abundance |
#   Taxon_ID | Group
###############################################################################

# Clade-level phyloseq subsets (needed for correct prevalence denominators)
ps_16S_amph <- subset_samples(ps_16S_genus, Clade == "Amphibia")
ps_16S_rept <- subset_samples(ps_16S_genus, Clade == "Reptilia")
ps_ITS_amph <- subset_samples(ps_ITS_genus, Clade == "Amphibia")
ps_ITS_rept <- subset_samples(ps_ITS_genus, Clade == "Reptilia")

# Returns a tidy taxonomy + prevalence + abundance table for a list of core taxa.
# Uses absolute abundance values directly from the (spike-in normalized) OTU table.
build_core_table <- function(ps, core_taxa_list, group_label) {
  
  if (length(core_taxa_list) == 0) return(data.frame())
  
  tax   <- as.data.frame(tax_table(ps))
  otu   <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) otu <- t(otu)
  sdata <- data.frame(sample_data(ps))
  
  core_otu <- otu[core_taxa_list, , drop = FALSE]
  
  total_samples <- ncol(core_otu)
  total_species <- length(unique(sdata$host_taxon))
  
  # Sample-level prevalence (vectorized)
  n_samples_present <- rowSums(core_otu > 0)
  
  # Species-level prevalence
  n_species_present <- sapply(core_taxa_list, function(tid) {
    samples_with_taxon <- colnames(core_otu)[core_otu[tid, ] > 0]
    length(unique(sdata$host_taxon[rownames(sdata) %in% samples_with_taxon]))
  })
  
  # Absolute abundance metrics (DspikeIn-normalized counts)
  mean_abs   <- rowMeans(core_otu)
  median_abs <- apply(core_otu, 1, median)
  total_abs  <- rowSums(core_otu)
  
  # Assemble table
  tax[core_taxa_list, ] %>%
    rownames_to_column("Taxon_ID") %>%
    mutate(
      Sample_prevalence          = paste0(n_samples_present, "/", total_samples),
      Species_prevalence         = paste0(n_species_present, "/", total_species),
      mean_absolute_abundance    = round(mean_abs,   2),
      median_absolute_abundance  = round(median_abs, 0),
      total_absolute_abundance   = round(total_abs,  0),
      Group                      = group_label
    ) %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus,
           Sample_prevalence, Species_prevalence,
           mean_absolute_abundance, median_absolute_abundance,
           total_absolute_abundance,
           Taxon_ID, Group) %>%
    arrange(desc(
      as.numeric(sub("/.*", "", Sample_prevalence)) /
        as.numeric(sub(".*/", "", Sample_prevalence))
    ))
}

cat("\nBuilding supplementary table...\n")

supp_table <- bind_rows(
  build_core_table(ps_16S_amph, core_16S_amph, "Amphibia"),
  build_core_table(ps_16S_rept, core_16S_rept, "Reptilia"),
  build_core_table(ps_ITS_amph, core_ITS_amph, "Amphibia"),
  build_core_table(ps_ITS_rept, core_ITS_rept, "Reptilia")
)

write.csv(supp_table, file.path(dir_tables, "core_microbiome_WILD.csv"), row.names = FALSE)

cat("Saved:", file.path(dir_tables, "core_microbiome_WILD.csv"), "\n")
cat("  Rows:", nrow(supp_table),
    "| Bacteria:", sum(supp_table$Kingdom == "Bacteria"),
    "| Fungi:", sum(supp_table$Kingdom == "Fungi"), "\n")

##------------------------------------------------------------------------------
## Helper: build and save a core microbiome supplementary table for any
## two-group split, combining bacteria (16S) and fungi (ITS).
##
## Arguments:
##   ps_16S       — genus-level 16S phyloseq object (already subsetted as needed)
##   ps_ITS       — genus-level ITS phyloseq object (already subsetted as needed)
##   group_var    — metadata column to split on (must exist in both objects)
##   group_a / b  — group values (must match entries in group_var column)
##   label_a / b  — display labels for the Group column in the output CSV
##   csv_filename — output CSV filename
##------------------------------------------------------------------------------
make_core_supp_table <- function(ps_16S, ps_ITS, group_var,
                                 group_a, group_b,
                                 label_a, label_b,
                                 csv_filename) {
  
  build_kingdom_tables <- function(ps, kingdom_label) {
    sdata <- data.frame(sample_data(ps))
    ps_a  <- prune_samples(rownames(sdata)[sdata[[group_var]] == group_a], ps)
    ps_b  <- prune_samples(rownames(sdata)[sdata[[group_var]] == group_b], ps)
    ps_a  <- prune_taxa(taxa_sums(ps_a) > 0, ps_a)
    ps_b  <- prune_taxa(taxa_sums(ps_b) > 0, ps_b)
    
    core_a <- get_core_taxa(ps_a)
    core_b <- get_core_taxa(ps_b)
    
    cat(" ", kingdom_label, "— core genera:",
        label_a, ":", length(core_a),
        "| ", label_b, ":", length(core_b),
        "| Shared:", length(intersect(core_a, core_b)), "\n")
    
    bind_rows(
      build_core_table(ps_a, core_a, label_a),
      build_core_table(ps_b, core_b, label_b)
    )
  }
  
  tbl <- bind_rows(
    build_kingdom_tables(ps_16S, "Bacteria"),
    build_kingdom_tables(ps_ITS,  "Fungi")
  )
  
  write.csv(tbl, csv_filename, row.names = FALSE)
  cat("  Saved:", csv_filename,
      "| Rows:", nrow(tbl),
      "| Bacteria:", sum(tbl$Kingdom == "Bacteria"),
      "| Fungi:", sum(tbl$Kingdom == "Fungi"), "\n")
  
  invisible(tbl)
}

###############################################################################
# ========================== PART 4: BASIDIOBOLUS CHECK =======================
###############################################################################
cat("\n========== CHECKING FOR BASIDIOBOLUS ==========\n")

check_basidiobolus <- function(df, label) {
  hit <- filter(df, grepl("Basidiobolus", Genus, ignore.case = TRUE))
  if (nrow(hit) > 0) {
    cat("✓ Basidiobolus is core in:", label, "\n")
    cat("  Sample prevalence:", hit$Sample_prevalence, "\n")
    cat("  Species prevalence:", hit$Species_prevalence, "\n")
  } else {
    cat("✗ Basidiobolus NOT core in:", label, "\n")
  }
  invisible(hit)
}

check_basidiobolus(filter(supp_table, Kingdom == "Fungi"),           "All groups (combined table)")
check_basidiobolus(filter(supp_table, Kingdom == "Fungi", Group == "Amphibia"), "Amphibia")
check_basidiobolus(filter(supp_table, Kingdom == "Fungi", Group == "Reptilia"), "Reptilia")

###############################################################################
# ========================== PART 5: PREVALENCE DISTRIBUTIONS =================
###############################################################################

# Returns genus-level prevalence and mean abundance across all samples
calc_prevalence_df <- function(ps, dataset_label) {
  otu        <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) otu <- t(otu)
  prevalence <- rowSums(otu > 0) / ncol(otu)
  data.frame(
    Taxon_ID       = names(prevalence),
    prevalence     = prevalence,
    mean_abundance = rowMeans(otu),
    dataset        = dataset_label,
    stringsAsFactors = FALSE
  )
}

cat("\nCalculating prevalence distributions...\n")
prev_16S_all <- calc_prevalence_df(ps_16S_genus, "Bacteria")
prev_ITS_all <- calc_prevalence_df(ps_ITS_genus, "Fungi")

###############################################################################
# ========================== PART 6: FIGURES ==================================
###############################################################################
cat("\nBuilding multi-panel figures...\n")

# ---- Panel A: Genus-level prevalence histogram ----
prev_all <- bind_rows(prev_16S_all, prev_ITS_all) %>%
  mutate(dataset = factor(dataset, levels = c("Bacteria", "Fungi")))

p_prev_all <- ggplot(prev_all, aes(x = prevalence, fill = dataset)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = PREV_THRESHOLD, linetype = "dashed",
             color = "red", linewidth = 0.5) +
  annotate("text", x = PREV_THRESHOLD + 0.02, y = Inf,
           label = "Core threshold\n(50%)", hjust = 0, vjust = 1.5,
           size = 3, color = "red") +
  scale_fill_manual(values = dataset_colors) +
  theme_classic() +
  labs(
    x    = "Prevalence (proportion of samples)",
    y    = "Number of genera",
    title = "Genus-level prevalence across all wild hosts",
    fill  = "Dataset"
  ) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# ---- Panel B: Prevalence vs. mean absolute abundance ----
prev_all_labeled <- prev_all %>%
  mutate(core_status = ifelse(prevalence >= PREV_THRESHOLD,
                              "Core (\u226550%)", "Non-core (<50%)"))

p_prev_abund <- ggplot(prev_all_labeled,
                       aes(x = prevalence,
                           y = log10(mean_abundance + 1),
                           color = dataset,
                           shape = core_status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_vline(xintercept = PREV_THRESHOLD, linetype = "dashed",
             color = "red", linewidth = 0.5) +
  scale_color_manual(values = dataset_colors) +
  scale_shape_manual(values = c("Core (\u226550%)" = 16, "Non-core (<50%)" = 1)) +
  facet_wrap(~ dataset, scales = "free_y") +
  theme_classic() +
  labs(
    x     = "Prevalence (proportion of samples)",
    y     = expression(log[10]~"(mean absolute abundance + 1)"),
    title = "Prevalence vs. abundance across all wild hosts",
    color = "Dataset",
    shape = "Status"
  ) +
  theme(legend.position   = "bottom",
        strip.background  = element_rect(fill = "grey92", color = NA),
        strip.text        = element_text(face = "bold"),
        plot.title        = element_text(hjust = 0.5, face = "bold"))

# ---- Panel C: Venn diagram — bacterial core by clade ----
venn.plot.16S <- suppressMessages(
  venn.diagram(
    x              = list(Amphibia = core_16S_amph, Reptilia = core_16S_rept),
    category.names = c("Amphibia", "Reptilia"),
    filename       = NULL,
    fill           = c(clade_colors["Amphibia"], clade_colors["Reptilia"]),
    alpha          = 0.5,
    cex            = 1.5,
    cat.cex        = 1.2,
    cat.pos        = c(-20, 20),
    cat.dist       = c(0.05, 0.05),
    main           = "Core bacterial genera by host clade",
    main.cex       = 1.2,
    main.fontface  = "bold"
  )
)
p_venn_16S <- ggdraw() + draw_grob(grobTree(venn.plot.16S))

# ---- Helper: parse prevalence column from supp_table into numeric % ----
# Genus_clean rules (applied in order):
#   1. Rozellomycota incertae sedis  → "Rozellomycota sp."
#   2. Any genus matching "uncultured" (case-insensitive) → "{Family}_uncultured"
#   3. All other genera kept as-is.
parse_prev_pct <- function(df, group) {
  df %>%
    filter(Group == group) %>%
    mutate(
      n_samples      = as.numeric(sub("/.*", "", Sample_prevalence)),
      total_samples  = as.numeric(sub(".*/", "", Sample_prevalence)),
      prevalence_pct = n_samples / total_samples * 100,
      Genus_clean    = case_when(
        grepl("Rozellomycota.*incertae", Genus, ignore.case = TRUE) ~
          "Rozellomycota sp.",
        grepl("uncultured", Genus, ignore.case = TRUE) ~
          paste0(Family, "_uncultured"),
        TRUE ~ Genus
      ),
      Clade = group
    ) %>%
    select(Phylum, Family, Genus_clean, prevalence_pct, Clade, total_samples)
}

# ---- Panel D (fungi version): Core fungal genera lollipop ----
core_fungi_prev_df <- bind_rows(
  supp_table %>% filter(Kingdom == "Fungi") %>% parse_prev_pct("Amphibia"),
  supp_table %>% filter(Kingdom == "Fungi") %>% parse_prev_pct("Reptilia")
) %>% mutate(Clade = factor(Clade, levels = c("Amphibia", "Reptilia")))

p_fungi_prev <- ggplot(core_fungi_prev_df,
                       aes(x     = reorder(Genus_clean, prevalence_pct),
                           y     = prevalence_pct,
                           color = Clade)) +
  geom_segment(aes(xend = Genus_clean, y = 50, yend = prevalence_pct),
               linewidth = 1.5, alpha = 0.7) +
  geom_point(size = 4, alpha = 0.9) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_color_manual(values = clade_colors) +
  coord_flip() +
  ylim(50, 100) +
  theme_classic() +
  labs(
    x     = NULL,
    y     = "Prevalence (% of samples)",
    title = "Core fungal genera by host clade",
    color = "Clade"
  ) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.y     = element_text(face = "italic")
  )

# ---- Assemble and save 4-panel figure (Panels A–D, fungi lollipop) ----

fig_core <- (p_prev_all | p_prev_abund) /
  (p_venn_16S  | p_fungi_prev) +
  plot_annotation(
    title      = "Core Microbiome — Genus Level (Wild Samples Only)",
    subtitle   = "Core defined as \u226550% prevalence across individual samples within each group",
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

ggsave(file.path(dir_figures, "core_microbiome_GENUS_WILD.pdf"), fig_core, width = 13, height = 11, device = cairo_pdf, dpi = 400)
cat("Saved: core microbiome figure\n")

##------------------------------------------------------------------------------
## Standalone bacteria lollipop — Amphibia | Reptilia side-by-side
## Phylum-colored; shared-core genera flagged with *
## Genera grouped by family: families with ≥3 core genera get their own facet;
## all others go into a "Singletons" facet at the bottom.
##------------------------------------------------------------------------------

core_bact_prev_df <- bind_rows(
  supp_table %>% filter(Kingdom == "Bacteria") %>% parse_prev_pct("Amphibia"),
  supp_table %>% filter(Kingdom == "Bacteria") %>% parse_prev_pct("Reptilia")
) %>%
  mutate(
    Clade      = factor(Clade, levels = c("Amphibia", "Reptilia")),
    shared     = Genus_clean %in% intersect(
      filter(supp_table, Kingdom == "Bacteria", Group == "Amphibia")$Genus,
      filter(supp_table, Kingdom == "Bacteria", Group == "Reptilia")$Genus
    ),
    label_text = if_else(shared, paste0(Genus_clean, "*"), Genus_clean)
  )

# Phylum palette — matches host phylogeny figure colors; any phyla not listed
# there fall back to additional colorblind-friendly values
phylum_palette <- c(
  "Actinobacteriota"  = "#A6CEE3",
  "Bacteroidota"      = "#1F78B4",
  "Campilobacterota"  = "#B2DF8A",
  "Desulfobacterota"  = "#FB9A99",
  "Firmicutes"        = "#33A02C",
  "Fusobacteriota"    = "#E31A1C",
  "Planctomycetota"   = "#6A3D9A",
  "Proteobacteria"    = "#FDBF6F",
  "Synergistota"      = "#6A3D9A",
  "Verrucomicrobiota" = "#FF7F00",
  # Additional phyla not in host phylogeny figure
  "Cyanobacteria"     = "#CAB2D6",
  "Deferribacterota"  = "#FFFF99",
  "Spirochaetota"     = "#B15928",
  "Chloroflexota"     = "#009E73",
  "Other"             = "#222222"
)

##------------------------------------------------------------------------------
## Family-grouping logic
##
## Within each clade panel, genera are grouped by family:
##   - Families with ≥ 3 core genera in that clade  → own named facet strip
##   - All other families (< 3 genera)              → "Singletons" facet at bottom
##
## Ordering within each facet: ascending prevalence (lowest at bottom) so the
## most-prevalent genus sits closest to the facet boundary — mirrors the
## original lollipop sort direction.
##
## Facet order: named family groups (alpha) on top, "Singletons" always last.
##------------------------------------------------------------------------------

FAMILY_GROUP_MIN <- 3   # minimum genera per family to earn its own facet

assign_family_group <- function(df) {
  # Count distinct Genus_clean per Family within this clade slice
  fam_counts <- df %>%
    group_by(Family) %>%
    summarise(n_genera = n_distinct(Genus_clean), .groups = "drop")
  
  df %>%
    left_join(fam_counts, by = "Family") %>%
    mutate(
      family_group = if_else(n_genera >= FAMILY_GROUP_MIN,
                             Family,
                             "Singletons")
    )
}

# Build one panel per clade with family-grouped facets
make_bact_panel <- function(clade) {
  df <- filter(core_bact_prev_df, Clade == clade) %>%
    assign_family_group()
  
  # Determine facet order: named families (alphabetical) then "Singletons" last
  named_fams <- sort(unique(df$family_group[df$family_group != "Singletons"]))
  facet_order <- c(named_fams, "Singletons")
  df <- df %>%
    mutate(family_group = factor(family_group, levels = facet_order))
  
  # Within each facet sort by ascending prevalence, then set factor levels
  # so coord_flip puts lowest prevalence at the bottom of the strip
  df <- df %>%
    arrange(family_group, prevalence_pct) %>%
    mutate(label_text = factor(label_text, levels = unique(label_text)))
  
  ggplot(df, aes(x = label_text, y = prevalence_pct, color = Phylum)) +
    geom_segment(aes(xend = label_text, y = 50, yend = prevalence_pct),
                 linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.9) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red",
               linewidth = 0.5) +
    scale_color_manual(values = phylum_palette, drop = FALSE) +
    # Free y-scales so each facet strip is sized to its own genera count
    facet_grid(family_group ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    coord_flip() +
    ylim(50, 100) +
    theme_classic() +
    labs(
      x     = NULL,
      y     = "Prevalence (% of samples)",
      title = paste0(clade, " (n = ",unique(filter(core_bact_prev_df, Clade == clade)$total_samples)," samples)"),
      color = "Phylum"
    ) +
    theme(
      legend.position   = "right",
      legend.key.size   = unit(0.4, "cm"),
      legend.text       = element_text(size = 8),
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 12),
      # Genus labels: italic (family name is in strip, not in axis label)
      axis.text.y       = element_text(face = "italic", size = 10),
      # Facet strip styling: plain (non-italic) family names, left-side placement
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "grey92", color = NA),
      strip.text.y.left = element_text(face = "plain", size = 8,
                                       angle = 0, hjust = 1),
      panel.spacing     = unit(0.3, "lines")
    )
}

p_bact_amph <- make_bact_panel("Amphibia")
p_bact_rept <- make_bact_panel("Reptilia")

# Height scales with total genus count across both panels (facets add overhead)
bact_fig_height <- max(8, max(
  nrow(filter(core_bact_prev_df, Clade == "Amphibia")),
  nrow(filter(core_bact_prev_df, Clade == "Reptilia"))
) * 0.28 + 4)

fig_bact_lollipop <- (p_bact_amph | p_bact_rept) +
  plot_annotation(
    title      = "Core bacterial genera by host clade (wild samples only)",
    subtitle   = paste0(
      "* = shared core (present in both Amphibia and Reptilia); color = phylum\n",
      "Families with \u2265", FAMILY_GROUP_MIN,
      " core genera shown as named groups; remaining genera in \u2018Singletons\u2019"
    ),
    tag_levels = list(c("A)", "B)")),
    theme = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      plot.tag      = element_text(face = "bold", size = 12)
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(dir_figures, "core_bacteria_lollipop_WILD.pdf"), fig_bact_lollipop, width = 14, height = bact_fig_height, device = cairo_pdf, dpi = 400)

cat("Saved: bacteria lollipop figure\n")

##------------------------------------------------------------------------------
## Helper: core lollipop figure for any two-group split
##
## Arguments:
##   ps            — genus-level phyloseq object (wild samples, already subsetted)
##   group_var     — metadata column to split on (e.g. "Clade_Order", "animal_type")
##   group_a       — value for panel A
##   group_b       — value for panel B
##   label_a       — display title for panel A
##   label_b       — display title for panel B
##   fig_filename  — output PDF filename (without path)
##   fig_title     — overall figure title
##------------------------------------------------------------------------------

make_bact_lollipop_figure <- function(ps, group_var,
                                      group_a, group_b,
                                      label_a, label_b,
                                      fig_filename, fig_title) {
  
  # ---- Compute core taxa for each group ----
  core_a <- get_core_by_group(ps, group_var, group_a)
  core_b <- get_core_by_group(ps, group_var, group_b)
  
  cat("  Core genera —", label_a, ":", length(core_a),
      "| ", label_b, ":", length(core_b),
      "| Shared:", length(intersect(core_a, core_b)), "\n")
  
  # ---- Subset phyloseq objects for correct prevalence denominators ----
  sdata  <- data.frame(sample_data(ps))
  ps_a   <- prune_samples(rownames(sdata)[sdata[[group_var]] == group_a], ps)
  ps_b   <- prune_samples(rownames(sdata)[sdata[[group_var]] == group_b], ps)
  
  # ---- Build prevalence data frames ----
  # Reuses build_core_table + parse_prev_pct already defined above.
  # parse_prev_pct expects a Group column, so we label each subset accordingly.
  tbl_a <- build_core_table(ps_a, core_a, label_a)
  tbl_b <- build_core_table(ps_b, core_b, label_b)
  
  prev_df <- bind_rows(
    tbl_a %>% parse_prev_pct(label_a),
    tbl_b %>% parse_prev_pct(label_b)
  ) %>%
    mutate(
      Clade      = factor(Clade, levels = c(label_a, label_b)),
      shared     = Genus_clean %in% intersect(tbl_a$Genus, tbl_b$Genus),
      label_text = if_else(shared, paste0(Genus_clean, "*"), Genus_clean)
    )
  
  # ---- Panel builder (same family-grouping logic as main figure) ----
  make_panel <- function(grp_label) {
    df <- filter(prev_df, Clade == grp_label) %>%
      assign_family_group()
    
    named_fams <- sort(unique(df$family_group[df$family_group != "Singletons"]))
    df <- df %>%
      mutate(family_group = factor(family_group,
                                   levels = c(named_fams, "Singletons"))) %>%
      arrange(family_group, prevalence_pct) %>%
      mutate(label_text = factor(label_text, levels = unique(label_text)))
    
    ggplot(df, aes(x = label_text, y = prevalence_pct, color = Phylum)) +
      geom_segment(aes(xend = label_text, y = 50, yend = prevalence_pct),
                   linewidth = 1.2, alpha = 0.7) +
      geom_point(size = 3, alpha = 0.9) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "red",
                 linewidth = 0.5) +
      scale_color_manual(values = phylum_palette, drop = FALSE) +
      facet_grid(family_group ~ ., scales = "free_y", space = "free_y",
                 switch = "y") +
      coord_flip() +
      ylim(50, 100) +
      theme_classic() +
      labs(
        x     = NULL,
        y     = "Prevalence (% of samples)",
        title = paste0(grp_label, " (n = ",unique(filter(prev_df, Clade == grp_label)$total_samples)," samples)"),
        color = "Phylum"
      ) +
      theme(
        legend.position   = "right",
        legend.key.size   = unit(0.4, "cm"),
        legend.text       = element_text(size = 8),
        plot.title        = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.y       = element_text(face = "italic", size = 10),
        strip.placement   = "outside",
        strip.background  = element_rect(fill = "grey92", color = NA),
        strip.text.y.left = element_text(face = "plain", size = 8,
                                         angle = 0, hjust = 1),
        panel.spacing     = unit(0.3, "lines")
      )
  }
  
  p_a <- make_panel(label_a)
  p_b <- make_panel(label_b)
  
  fig_height <- max(8, max(
    nrow(filter(prev_df, Clade == label_a)),
    nrow(filter(prev_df, Clade == label_b))
  ) * 0.28 + 4)
  
  fig <- (p_a | p_b) +
    plot_annotation(
      title      = fig_title,
      subtitle   = paste0(
        "* = shared core (present in both groups); color = phylum\n",
        "Families with \u2265", FAMILY_GROUP_MIN,
        " core genera shown as named groups; remaining genera in \u2018Singletons\u2019"
      ),
      tag_levels = list(c("A)", "B)")),
      theme = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        plot.tag      = element_text(face = "bold", size = 12)
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
  
  ggsave(fig_filename, fig, width = 14, height = fig_height,
         device = cairo_pdf, dpi = 400)
  cat("Saved:", fig_filename, "\n")
}

##------------------------------------------------------------------------------
## Figure: Amphibian orders — Caudata (A) vs. Anura (B)
##------------------------------------------------------------------------------
cat("\nBuilding amphibian orders lollipop (Caudata vs. Anura)...\n")

make_bact_lollipop_figure(
  ps           = ps_16S_genus,
  group_var    = "Clade_Order",
  group_a      = "Caudata",
  group_b      = "Anura",
  label_a      = "Caudata",
  label_b      = "Anura",
  fig_filename = file.path(dir_figures, "core_bacteria_lollipop_amphibian_orders.pdf"),
  fig_title    = "Core bacterial genera — amphibian orders (wild samples only)"
)

cat("\nBuilding amphibian orders supplementary table (Caudata vs. Anura)...\n")
make_core_supp_table(
  ps_16S       = ps_16S_genus,
  ps_ITS       = ps_ITS_genus,
  group_var    = "Clade_Order",
  group_a      = "Caudata",
  group_b      = "Anura",
  label_a      = "Caudata",
  label_b      = "Anura",
  csv_filename = file.path(dir_tables, "core_microbiome_amphibian_orders_WILD.csv")
)

##------------------------------------------------------------------------------
## Figure: Squamata — Lizards (A) vs. Snakes (B)
##------------------------------------------------------------------------------
cat("\nBuilding Squamata lollipop (Lizards vs. Snakes)...\n")

# Subset ps_16S_genus to Squamata only before passing in, so prevalence
# denominators are correct (i.e. % of lizard samples, not all reptile samples)
ps_16S_squamata <- subset_samples(ps_16S_genus, Clade_Order == "Squamata")
ps_16S_squamata <- prune_taxa(taxa_sums(ps_16S_squamata) > 0, ps_16S_squamata)

make_bact_lollipop_figure(
  ps           = ps_16S_squamata,
  group_var    = "animal_type",
  group_a      = "Lizard",
  group_b      = "Snake",
  label_a      = "Lizards",
  label_b      = "Snakes",
  fig_filename = file.path(dir_figures, "core_bacteria_lollipop_squamata.pdf"),
  fig_title    = "Core bacterial genera — Squamata (wild samples only)"
)

cat("\nBuilding Squamata supplementary table (Lizards vs. Snakes)...\n")

# Subset ITS to Squamata only — same rationale as for 16S above
ps_ITS_squamata <- subset_samples(ps_ITS_genus, Clade_Order == "Squamata")
ps_ITS_squamata <- prune_taxa(taxa_sums(ps_ITS_squamata) > 0, ps_ITS_squamata)

make_core_supp_table(
  ps_16S       = ps_16S_squamata,
  ps_ITS       = ps_ITS_squamata,
  group_var    = "animal_type",
  group_a      = "Lizard",
  group_b      = "Snake",
  label_a      = "Lizards",
  label_b      = "Snakes",
  csv_filename = file.path(dir_tables, "core_microbiome_squamata_WILD.csv")
)


################################################################################
################################################################################
################################################################################
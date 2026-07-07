###############################################################################
# Figure 2 — Host phylogeny with bacterial (16S) and fungal (ITS) community
# composition (absolute abundance data -> relative abundance) and 
# host ecological metadata
#
# Loads post-decontam + metadata-cleaned phyloseq objects, applies species
# name substitutions to align with TimeTree tip labels, aggregates to phylum-
# level relative abundance per host species, and plots alongside the host
# phylogeny using ggtree + ggtreeExtra.
#
# Also produces supplementary figure (wild/fecal+LGI only) at the bottom.
#
# Alexander Rurik
###############################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(vegan)
  library(RColorBrewer)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(cowplot)
  library(ape)
})

# ========================== BASIC PARAMETERS ==================================
set.seed(52325)

# Explicitly prioritize dplyr functions
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
first  <- dplyr::first

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_trees     <- here("data", "raw", "trees")
dir_figures   <- here("output", "figures")

# ========================== LOAD DATA =========================================
ps_bact_abs_rel <- readRDS(file.path(dir_processed, "16S_absolute_all_final_ALL_SPECIES.rds"))
ps_fung_abs_rel <- readRDS(file.path(dir_processed, "ITS_absolute_all_final_ALL_SPECIES.rds"))

host_tree <- read.tree(file.path(dir_trees, "tree_common_hosts_clean.nwk"))
host_tree$tip.label <- gsub("_", " ", host_tree$tip.label)
tree_hosts <- host_tree$tip.label

# ========================== SPECIES FILTERING & NAME ALIGNMENT ================

# Remove Desmognathus adatsihi — unresolved species not represented in TimeTree
ps_bact_abs_rel <- prune_samples(
  sample_data(ps_bact_abs_rel)$host_taxon != "Desmognathus adatsihi", ps_bact_abs_rel)
ps_fung_abs_rel <- prune_samples(
  sample_data(ps_fung_abs_rel)$host_taxon != "Desmognathus adatsihi", ps_fung_abs_rel)

# Species name substitutions to align microbiome metadata with TimeTree tip labels.
# (a) Subspecies collapsed to species — TimeTree does not resolve subspecies
# (b) Species/genus substitutions — taxonomic synonyms or updates used in TimeTree
# NOTE: species_subs is also used by the supplementary figure section below
species_subs <- c(
  # --- Subspecies collapsed to species ---
  "Sceloporus occidentalis bocourtii" = "Sceloporus occidentalis",
  "Pituophis catenifer pumilus"       = "Pituophis catenifer",
  "Pituophis catenifer catenifer"     = "Pituophis catenifer",
  "Thamnophis elegans terrestris"     = "Thamnophis elegans",
  "Thamnophis atratus atratus"        = "Thamnophis atratus",
  "Dendrobates tinctorius azureus"    = "Dendrobates tinctorius",
  "Coluber constrictor mormon"        = "Coluber constrictor",
  # --- Species substitutions (TimeTree synonyms) ---
  "Agama picticauda"                  = "Agama atra",
  "Terrapene carolina"                = "Terrapene ornata",
  "Ambystoma annulatum"               = "Ambystoma opacum",
  "Heloderma exasperatum"             = "Heloderma horridum",
  "Atelopus balios"                   = "Atelopus longirostris",
  # --- Genus updates (modern taxonomy / TimeTree naming) ---
  "Hyla cinerea"                      = "Dryophytes cinereus",
  "Hyla avivoca"                      = "Dryophytes avivoca",
  "Hyla chrysoscelis"                 = "Dryophytes chrysoscelis",
  "Geochelone gigantea"               = "Aldabrachelys gigantea"
)

sample_data(ps_bact_abs_rel)$host_taxon <- as.character(sample_data(ps_bact_abs_rel)$host_taxon)
sample_data(ps_fung_abs_rel)$host_taxon <- as.character(sample_data(ps_fung_abs_rel)$host_taxon)

sample_data(ps_bact_abs_rel)$host_taxon <- recode(
  sample_data(ps_bact_abs_rel)$host_taxon, !!!species_subs,
  .default = sample_data(ps_bact_abs_rel)$host_taxon)
sample_data(ps_fung_abs_rel)$host_taxon <- recode(
  sample_data(ps_fung_abs_rel)$host_taxon, !!!species_subs,
  .default = sample_data(ps_fung_abs_rel)$host_taxon)

# Retain only samples whose host_taxon matches a tree tip label
ps_bact_abs_rel <- prune_samples(
  sample_data(ps_bact_abs_rel)$host_taxon %in% host_tree$tip.label, ps_bact_abs_rel)
ps_fung_abs_rel <- prune_samples(
  sample_data(ps_fung_abs_rel)$host_taxon %in% host_tree$tip.label, ps_fung_abs_rel)

# Verify alignment
cat("Tree tips with no 16S data:\n")
print(setdiff(host_tree$tip.label, unique(sample_data(ps_bact_abs_rel)$host_taxon)))
cat("Tree tips with no ITS data:\n")
print(setdiff(host_tree$tip.label, unique(sample_data(ps_fung_abs_rel)$host_taxon)))

# Prune to species present in tree + 16S + ITS
species_with_both    <- intersect(
  unique(sample_data(ps_bact_abs_rel)$host_taxon),
  unique(sample_data(ps_fung_abs_rel)$host_taxon))
species_in_all_three <- intersect(species_with_both, host_tree$tip.label)

cat("Species in tree:", length(host_tree$tip.label), "\n")
cat("Species with 16S data:", length(unique(sample_data(ps_bact_abs_rel)$host_taxon)), "\n")
cat("Species with ITS data:", length(unique(sample_data(ps_fung_abs_rel)$host_taxon)), "\n")
cat("Species in all three (final figure):", length(species_in_all_three), "\n\n")

ps_bact_abs_rel <- prune_samples(
  sample_data(ps_bact_abs_rel)$host_taxon %in% species_in_all_three, ps_bact_abs_rel)
ps_fung_abs_rel <- prune_samples(
  sample_data(ps_fung_abs_rel)$host_taxon %in% species_in_all_three, ps_fung_abs_rel)
host_tree <- keep.tip(host_tree, species_in_all_three)

cat("Final verification — should both be 0:\n")
cat("  Tree tips with no 16S data:",
    length(setdiff(host_tree$tip.label, unique(sample_data(ps_bact_abs_rel)$host_taxon))), "\n")
cat("  Tree tips with no ITS data:",
    length(setdiff(host_tree$tip.label, unique(sample_data(ps_fung_abs_rel)$host_taxon))), "\n\n")

# ========================== RELATIVE ABUNDANCE AGGREGATION ===================
# Aggregate to phylum-level relative abundance per host species.
# Relative abundance calculated within each host species (not globally).
# Top 8 phyla by total abundance retained; all others collapsed to "Other".

# --- Bacteria ---
df_bact_rel <- psmelt(ps_bact_abs_rel) %>%
  group_by(Sample, host_taxon, Phylum) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(host_taxon) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  mutate(Phylum = ifelse(
    is.na(Phylum) | Phylum == "" | tolower(Phylum) %in% c("unassigned", "unknown", "na"),
    "Unknown", Phylum))

top_bact_rel <- df_bact_rel %>%
  filter(Phylum != "Unknown") %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%
  pull(Phylum)

cat("Top 8 bacterial phyla (main figure):\n"); print(top_bact_rel)

df_bact_rel <- df_bact_rel %>%
  mutate(Phylum = ifelse(Phylum %in% c("Unknown", "Other") | !(Phylum %in% top_bact_rel),
                         "Other", Phylum)) %>%
  group_by(host_taxon, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_bact_rel$Phylum <- factor(df_bact_rel$Phylum,
                             levels = c(sort(setdiff(unique(df_bact_rel$Phylum), "Other")), "Other"))

# --- Fungi ---
df_fung_rel <- psmelt(ps_fung_abs_rel) %>%
  group_by(Sample, host_taxon, Phylum) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(host_taxon) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  mutate(Phylum = ifelse(
    is.na(Phylum) | Phylum == "" | tolower(Phylum) %in% c("unassigned", "unknown", "na"),
    "Unknown", Phylum))

top_fung_rel <- df_fung_rel %>%
  filter(Phylum != "Unknown") %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%
  pull(Phylum)

df_fung_rel <- df_fung_rel %>%
  mutate(Phylum = ifelse(Phylum %in% c("Unknown", "Other") | !(Phylum %in% top_fung_rel),
                         "Other", Phylum)) %>%
  group_by(host_taxon, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_fung_rel$Phylum <- factor(df_fung_rel$Phylum,
                             levels = c(sort(setdiff(unique(df_fung_rel$Phylum), "Other")), "Other"))

# ---- Optional: save/load aggregated data to avoid re-running psmelt() ----
#saveRDS(df_bact_rel, file = file.path(dir_processed, "df_bact_ra_psmelt_result.rds"))
#saveRDS(df_fung_rel, file = file.path(dir_processed, "df_fungi_ra_psmelt_result.rds"))
#df_bact_rel <- readRDS(file.path(dir_processed, "df_bact_ra_psmelt_result.rds"))
#df_fung_rel <- readRDS(file.path(dir_processed, "df_fungi_ra_psmelt_result.rds"))

# ========================== ALIGN HOSTS TO TREE TIP ORDER ====================
common_species_rel <- intersect(host_tree$tip.label, unique(df_bact_rel$host_taxon))

df_bact_rel <- df_bact_rel %>% filter(host_taxon %in% common_species_rel)
df_fung_rel <- df_fung_rel %>% filter(host_taxon %in% common_species_rel)

df_bact_rel$host_taxon <- factor(df_bact_rel$host_taxon,
                                 levels = host_tree$tip.label[host_tree$tip.label %in% common_species_rel])
df_fung_rel$host_taxon <- factor(df_fung_rel$host_taxon,
                                 levels = host_tree$tip.label[host_tree$tip.label %in% common_species_rel])

# ========================== PREPARE HOST METADATA ============================
metadata_rel <- data.frame(sample_data(ps_bact_abs_rel)) %>%
  filter(host_taxon %in% host_tree$tip.label) %>%
  group_by(host_taxon) %>%
  summarise(
    Diet = first(Diet),
    env_broad_scale = case_when(
      any(env_broad_scale == "Captive") & any(env_broad_scale == "Wild") ~ "Both",
      any(env_broad_scale == "Captive")                                   ~ "Captive",
      any(env_broad_scale == "Wild")                                      ~ "Wild",
      TRUE                                                                ~ NA_character_),
    env_medium = case_when(
      any(env_medium == "Fecal") & any(env_medium == "Cloacal swab")     ~ "Fecal/Cloacal swab",
      any(env_medium == "Fecal")                                          ~ "Fecal",
      any(env_medium == "Cloacal swab")                                   ~ "Cloacal swab",
      any(env_medium == "Lower GI")                                       ~ "Lower GI",
      TRUE                                                                ~ NA_character_),
    .groups = "drop")

metadata_rel$env_broad_scale <- factor(metadata_rel$env_broad_scale,
                                       levels = c("Wild", "Captive", "Both"))
metadata_rel$env_medium      <- factor(metadata_rel$env_medium,
                                       levels = c("Fecal", "Cloacal swab", "Lower GI", "Fecal/Cloacal swab"))
metadata_rel$Diet            <- factor(metadata_rel$Diet,
                                       levels = c("Insectivore", "Omnivore", "Carnivore", "Herbivore"))

df_blank_rel <- metadata_rel %>% select(host_taxon) %>% mutate(blank = 1)

# ========================== COLOR PALETTES ====================================

# Bacteria: hardcoded per-phylum assignments to prevent color shuffling when
# the top-8 composition changes between dataset subsets or reruns.
# Colors follow RColorBrewer "Paired" order applied alphabetically to the
# confirmed top-8 phyla for this dataset, with Firmicutes/Desulfobacterota
# swapped for visual clarity. "Other" is always near-black (#222222).
#
# Verified top 8 for full dataset (psmelt abundance totals, April 2026):
#   Firmicutes, Proteobacteria, Bacteroidota, Verrucomicrobiota, Fusobacteriota,
#   Actinobacteriota, Desulfobacterota, Campilobacterota
#
# Verified top 8 for wild/fecal subset:
#   Firmicutes, Proteobacteria, Bacteroidota, Fusobacteriota, Verrucomicrobiota,
#   Actinobacteriota, Desulfobacterota, Planctomycetota
#
# NOTE: "Campilobacterota" is the spelling used in the SILVA taxonomy as
# classified by DADA2 in this pipeline — do not correct to "Campylobacterota"
#
# All phyla that could appear in either figure's top 8 are listed here.
# The subset line retains only those present in the current dataset run.
bact_colors_rel <- c(
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
  "Other"             = "#222222"
)
bact_colors_rel <- bact_colors_rel[names(bact_colors_rel) %in% levels(df_bact_rel$Phylum)]

# Confirm mapping — useful for checking which phyla are present this run
print(data.frame(Phylum = names(bact_colors_rel), Color = bact_colors_rel))

# Fungi: named custom palette; positional assignment is stable here because
# the same 8 phyla consistently dominate across both dataset subsets
fung_colors_rel <- c(
  "#41AC66", "#66FFFF", "#6366CC", "#FFD966", "#F6735C",
  "#EFCCE5", "#F266FF", "#AAF266", "#222222"
)
names(fung_colors_rel) <- levels(df_fung_rel$Phylum)

# Metadata annotation strips
captive_wild_colors_rel <- c("Captive" = "#E69F00", "Wild" = "#86B4E9", "Both" = "#444444")

sample_type_colors_rel <- c(
  "Fecal"              = "#1a1a1a",
  "Cloacal swab"       = "#666666",
  "Lower GI"           = "#b3b3b3",
  "Fecal/Cloacal swab" = "#e6e6e6"
)

diet_colors_rel <- c(
  "Insectivore" = "#FF8ED6",
  "Omnivore"    = "#6B5EFF",
  "Carnivore"   = "#FFB500",
  "Herbivore"   = "#00C060"
)

# ========================== BUILD MAIN FIGURE =================================
p_tree_rel <- ggtree(host_tree) +
  geom_tiplab(aes(label = label), size = 2, offset = 3) +
  
  # A. Management status (Wild / Captive / Both)
  new_scale_fill() +
  geom_fruit(data = metadata_rel, geom = geom_col,
             mapping = aes(y = host_taxon, x = 1, fill = env_broad_scale),
             width = 0.9, pwidth = 0.1, offset = 0.78) +
  scale_fill_manual(name = "A) Management Status",
                    values = captive_wild_colors_rel,
                    guide = guide_legend(order = 1)) +
  
  # B. Sample type
  new_scale_fill() +
  geom_fruit(data = metadata_rel, geom = geom_col,
             mapping = aes(y = host_taxon, x = 1, fill = env_medium),
             width = 0.9, pwidth = 0.1, offset = 0.01) +
  scale_fill_manual(name = "B) Sample Type",
                    values = sample_type_colors_rel,
                    guide = guide_legend(order = 2)) +
  
  # C. Host diet
  new_scale_fill() +
  geom_fruit(data = metadata_rel, geom = geom_col,
             mapping = aes(y = host_taxon, x = 1, fill = Diet),
             width = 0.9, pwidth = 0.1, offset = 0.01) +
  scale_fill_manual(name = "C) Host Diet",
                    values = diet_colors_rel,
                    guide = guide_legend(order = 3)) +
  
  # D. Bacterial phyla
  new_scale_fill() +
  geom_fruit(data = df_bact_rel, geom = geom_col,
             mapping = aes(y = host_taxon, x = Abundance, fill = Phylum),
             orientation = "y", width = 0.9, pwidth = 1.1, offset = 0.01) +
  scale_fill_manual(name = "D) Bacterial Phyla",
                    values = bact_colors_rel,
                    guide = guide_legend(order = 4)) +
  
  # E. Fungal phyla
  new_scale_fill() +
  geom_fruit(data = df_fung_rel, geom = geom_col,
             mapping = aes(y = host_taxon, x = Abundance, fill = Phylum),
             orientation = "y", width = 0.9, pwidth = 1.1, offset = 0.01) +
  scale_fill_manual(name = "E) Fungal Phyla",
                    values = fung_colors_rel,
                    guide = guide_legend(order = 5)) +
  
  # Blank spacer
  new_scale_fill() +
  geom_fruit(data = df_blank_rel, geom = geom_col,
             mapping = aes(y = host_taxon, x = blank),
             orientation = "y", width = 0.9, pwidth = 0.19, offset = 0.25,
             fill = "white", color = NA) +
  scale_fill_identity() +
  
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.4, "cm"),
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 14,
                                     margin = margin(b = 18)),
    plot.margin       = margin(t = 10, r = 10, b = 25, l = 10)
  ) +
  ggtitle("Top Bacterial and Fungal Phyla in Herpetofauna Gut Microbiomes")

print(p_tree_rel)

#ggsave(filename = file.path(dir_figures, "host_tree_microbiome_rel_abundance.pdf"), plot = p_tree_rel, width = 10, height = 12, units = "in", device = cairo_pdf, dpi = 750)

###############################################################################
###############################################################################
###############################################################################

###############################################################################
# SUPPLEMENTARY FIGURE — Wild samples only, Fecal and Lower GI only
#
# Supplementary version of Figure 2 restricted to wild-caught individuals
# with fecal or lower GI samples only. Removes captive husbandry effects
# and cloacal swab mucosal signal to show natural gut community composition.
#
# Reloads fresh phyloseq objects and host tree because the main figure
# section above modifies them in place. species_subs is reused from above.
###############################################################################

# ========================== RELOAD DATA (fresh copies) ========================
ps_bact_wild <- readRDS(file.path(dir_processed, "16S_absolute_all_final.rds"))
ps_fung_wild <- readRDS(file.path(dir_processed, "ITS_absolute_all_final.rds"))

host_tree_wild <- read.tree(file.path(dir_trees, "tree_common_hosts_clean.nwk"))
host_tree_wild$tip.label <- gsub("_", " ", host_tree_wild$tip.label)

# ========================== SPECIES FILTERING & NAME ALIGNMENT ================

ps_bact_wild <- prune_samples(
  sample_data(ps_bact_wild)$host_taxon != "Desmognathus adatsihi", ps_bact_wild)
ps_fung_wild <- prune_samples(
  sample_data(ps_fung_wild)$host_taxon != "Desmognathus adatsihi", ps_fung_wild)

# Filter to wild samples with fecal or lower GI sample types only
ps_bact_wild <- prune_samples(
  sample_data(ps_bact_wild)$env_broad_scale == "Wild" &
    sample_data(ps_bact_wild)$env_medium %in% c("Fecal", "Lower GI"), ps_bact_wild)
ps_fung_wild <- prune_samples(
  sample_data(ps_fung_wild)$env_broad_scale == "Wild" &
    sample_data(ps_fung_wild)$env_medium %in% c("Fecal", "Lower GI"), ps_fung_wild)

cat("Wild/fecal+LGI filter — samples retained:\n")
cat("  16S:", nsamples(ps_bact_wild), "\n")
cat("  ITS:", nsamples(ps_fung_wild), "\n\n")

# Apply same species substitutions as main figure (species_subs defined above)
sample_data(ps_bact_wild)$host_taxon <- as.character(sample_data(ps_bact_wild)$host_taxon)
sample_data(ps_fung_wild)$host_taxon <- as.character(sample_data(ps_fung_wild)$host_taxon)

sample_data(ps_bact_wild)$host_taxon <- recode(
  sample_data(ps_bact_wild)$host_taxon, !!!species_subs,
  .default = sample_data(ps_bact_wild)$host_taxon)
sample_data(ps_fung_wild)$host_taxon <- recode(
  sample_data(ps_fung_wild)$host_taxon, !!!species_subs,
  .default = sample_data(ps_fung_wild)$host_taxon)

ps_bact_wild <- prune_samples(
  sample_data(ps_bact_wild)$host_taxon %in% host_tree_wild$tip.label, ps_bact_wild)
ps_fung_wild <- prune_samples(
  sample_data(ps_fung_wild)$host_taxon %in% host_tree_wild$tip.label, ps_fung_wild)

# Prune to species with both 16S and ITS after wild/fecal filter
species_with_both_wild    <- intersect(
  unique(sample_data(ps_bact_wild)$host_taxon),
  unique(sample_data(ps_fung_wild)$host_taxon))
species_in_all_three_wild <- intersect(species_with_both_wild, host_tree_wild$tip.label)

cat("Species in tree:", length(host_tree_wild$tip.label), "\n")
cat("Species with 16S (wild/fecal):", length(unique(sample_data(ps_bact_wild)$host_taxon)), "\n")
cat("Species with ITS (wild/fecal):", length(unique(sample_data(ps_fung_wild)$host_taxon)), "\n")
cat("Species in supplementary figure:", length(species_in_all_three_wild), "\n\n")

ps_bact_wild <- prune_samples(
  sample_data(ps_bact_wild)$host_taxon %in% species_in_all_three_wild, ps_bact_wild)
ps_fung_wild <- prune_samples(
  sample_data(ps_fung_wild)$host_taxon %in% species_in_all_three_wild, ps_fung_wild)
host_tree_wild <- keep.tip(host_tree_wild, species_in_all_three_wild)

cat("Final verification — should both be 0:\n")
cat("  Tree tips with no 16S data:",
    length(setdiff(host_tree_wild$tip.label, unique(sample_data(ps_bact_wild)$host_taxon))), "\n")
cat("  Tree tips with no ITS data:",
    length(setdiff(host_tree_wild$tip.label, unique(sample_data(ps_fung_wild)$host_taxon))), "\n\n")

# ========================== RELATIVE ABUNDANCE AGGREGATION ====================

# --- Bacteria ---
df_bact_wild <- psmelt(ps_bact_wild) %>%
  group_by(Sample, host_taxon, Phylum) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(host_taxon) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  mutate(Phylum = ifelse(
    is.na(Phylum) | Phylum == "" | tolower(Phylum) %in% c("unassigned", "unknown", "na"),
    "Unknown", Phylum))

top_bact_wild <- df_bact_wild %>%
  filter(Phylum != "Unknown") %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%
  pull(Phylum)

cat("Top 8 bacterial phyla (wild/fecal figure):\n"); print(top_bact_wild)

df_bact_wild <- df_bact_wild %>%
  mutate(Phylum = ifelse(Phylum %in% c("Unknown", "Other") | !(Phylum %in% top_bact_wild),
                         "Other", Phylum)) %>%
  group_by(host_taxon, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_bact_wild$Phylum <- factor(df_bact_wild$Phylum,
                              levels = c(sort(setdiff(unique(df_bact_wild$Phylum), "Other")), "Other"))

# --- Fungi ---
df_fung_wild <- psmelt(ps_fung_wild) %>%
  group_by(Sample, host_taxon, Phylum) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(host_taxon) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  mutate(Phylum = ifelse(
    is.na(Phylum) | Phylum == "" | tolower(Phylum) %in% c("unassigned", "unknown", "na"),
    "Unknown", Phylum))

top_fung_wild <- df_fung_wild %>%
  filter(Phylum != "Unknown") %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%
  pull(Phylum)

df_fung_wild <- df_fung_wild %>%
  mutate(Phylum = ifelse(Phylum %in% c("Unknown", "Other") | !(Phylum %in% top_fung_wild),
                         "Other", Phylum)) %>%
  group_by(host_taxon, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_fung_wild$Phylum <- factor(df_fung_wild$Phylum,
                              levels = c(sort(setdiff(unique(df_fung_wild$Phylum), "Other")), "Other"))

# ---- Optional: save/load aggregated data to avoid re-running psmelt() ----
# saveRDS(df_bact_wild, file = file.path(dir_processed, "df_bact_ra_psmelt_WILD_result.rds"))
# saveRDS(df_fung_wild, file = file.path(dir_processed, "df_fungi_ra_psmelt_WILD_result.rds"))
# df_bact_wild <- readRDS(file.path(dir_processed, "df_bact_ra_psmelt_WILD_result.rds"))
# df_fung_wild <- readRDS(file.path(dir_processed, "df_fungi_ra_psmelt_WILD_result.rds"))

# ========================== ALIGN HOSTS TO TREE TIP ORDER ====================
common_species_wild <- intersect(host_tree_wild$tip.label, unique(df_bact_wild$host_taxon))

df_bact_wild <- df_bact_wild %>% filter(host_taxon %in% common_species_wild)
df_fung_wild <- df_fung_wild %>% filter(host_taxon %in% common_species_wild)

df_bact_wild$host_taxon <- factor(df_bact_wild$host_taxon,
                                  levels = host_tree_wild$tip.label[host_tree_wild$tip.label %in% common_species_wild])
df_fung_wild$host_taxon <- factor(df_fung_wild$host_taxon,
                                  levels = host_tree_wild$tip.label[host_tree_wild$tip.label %in% common_species_wild])

# ========================== PREPARE HOST METADATA ============================
# Wild-only: env_broad_scale always "Wild"; sample types are Fecal and/or Lower GI only
metadata_wild <- data.frame(sample_data(ps_bact_wild)) %>%
  filter(host_taxon %in% host_tree_wild$tip.label) %>%
  group_by(host_taxon) %>%
  summarise(
    Diet            = first(Diet),
    env_broad_scale = first(env_broad_scale),
    env_medium      = case_when(
      any(env_medium == "Fecal") & any(env_medium == "Lower GI") ~ "Fecal/Lower GI",
      any(env_medium == "Fecal")                                  ~ "Fecal",
      any(env_medium == "Lower GI")                               ~ "Lower GI",
      TRUE                                                        ~ NA_character_),
    .groups = "drop")

metadata_wild$env_broad_scale <- factor(metadata_wild$env_broad_scale, levels = "Wild")
metadata_wild$env_medium      <- factor(metadata_wild$env_medium,
                                        levels = c("Fecal", "Lower GI", "Fecal/Lower GI"))
metadata_wild$Diet            <- factor(metadata_wild$Diet,
                                        levels = c("Insectivore", "Omnivore", "Carnivore", "Herbivore"))

df_blank_wild <- metadata_wild %>% select(host_taxon) %>% mutate(blank = 1)

# ========================== COLOR PALETTES ====================================
# Bacteria: same master palette as main figure — subset line handles which
# phyla are actually present in the wild/fecal dataset.
# Campilobacterota drops out (14th in wild/fecal); Cyanobacteria enters (8th).
# NOTE: these are hard-coded and may need to be changed based on the dataset.
bact_colors_wild <- c(
  "Actinobacteriota"  = "#A6CEE3",
  "Bacteroidota"      = "#1F78B4",
  "Campilobacterota"  = "#B2DF8A",
  "Desulfobacterota"  = "#FB9A99",  
  "Firmicutes"        = "#33A02C", 
  "Fusobacteriota"    = "#E31A1C",
  "Cyanobacteria"     = "#6A3D9A",
  "Proteobacteria"    = "#FDBF6F",
  "Synergistota"      = "#6A3D9A",
  "Verrucomicrobiota" = "#FF7F00",
  "Other"             = "#222222"
)
bact_colors_wild <- bact_colors_wild[names(bact_colors_wild) %in% levels(df_bact_wild$Phylum)]

print(data.frame(Phylum = names(bact_colors_wild), Color = bact_colors_wild))

# Fungi: same pool as main figure; subset to phyla present
fung_color_pool <- c(
  "#41AC66", "#66FFFF", "#6366CC", "#FFD966", "#F6735C",
  "#EFCCE5", "#F266FF", "#AAF266", "#222222"
)
fung_colors_wild <- setNames(
  fung_color_pool[seq_along(levels(df_fung_wild$Phylum))],
  levels(df_fung_wild$Phylum))
if ("Other" %in% names(fung_colors_wild)) fung_colors_wild["Other"] <- "#222222"

captive_wild_colors_wild <- c("Wild" = "#86B4E9")

sample_type_colors_wild <- c(
  "Fecal"          = "#1a1a1a",
  "Lower GI"       = "#b3b3b3",
  "Fecal/Lower GI" = "#e6e6e6"
)

diet_colors_wild <- diet_colors_rel  # reuse from main figure

# ========================== BUILD SUPPLEMENTARY FIGURE =======================
p_tree_wild <- ggtree(host_tree_wild) +
  geom_tiplab(aes(label = label), size = 2, offset = 3) +
  
  # A. Management status (Wild only)
  new_scale_fill() +
  geom_fruit(data = metadata_wild, geom = geom_col,
             mapping = aes(y = host_taxon, x = 1, fill = env_broad_scale),
             width = 0.9, pwidth = 0.1, offset = 0.78) +
  scale_fill_manual(name = "A) Management Status",
                    values = captive_wild_colors_wild,
                    guide = guide_legend(order = 1)) +
  
  # B. Sample type (Fecal / Lower GI only)
  new_scale_fill() +
  geom_fruit(data = metadata_wild, geom = geom_col,
             mapping = aes(y = host_taxon, x = 1, fill = env_medium),
             width = 0.9, pwidth = 0.1, offset = 0.01) +
  scale_fill_manual(name = "B) Sample Type",
                    values = sample_type_colors_wild,
                    guide = guide_legend(order = 2)) +
  
  # C. Host diet
  new_scale_fill() +
  geom_fruit(data = metadata_wild, geom = geom_col,
             mapping = aes(y = host_taxon, x = 1, fill = Diet),
             width = 0.9, pwidth = 0.1, offset = 0.01) +
  scale_fill_manual(name = "C) Host Diet",
                    values = diet_colors_wild,
                    guide = guide_legend(order = 3)) +
  
  # D. Bacterial phyla
  new_scale_fill() +
  geom_fruit(data = df_bact_wild, geom = geom_col,
             mapping = aes(y = host_taxon, x = Abundance, fill = Phylum),
             orientation = "y", width = 0.9, pwidth = 1.1, offset = 0.01) +
  scale_fill_manual(name = "D) Bacterial Phyla",
                    values = bact_colors_wild,
                    guide = guide_legend(order = 4)) +
  
  # E. Fungal phyla
  new_scale_fill() +
  geom_fruit(data = df_fung_wild, geom = geom_col,
             mapping = aes(y = host_taxon, x = Abundance, fill = Phylum),
             orientation = "y", width = 0.9, pwidth = 1.1, offset = 0.01) +
  scale_fill_manual(name = "E) Fungal Phyla",
                    values = fung_colors_wild,
                    guide = guide_legend(order = 5)) +
  
  # Blank spacer
  new_scale_fill() +
  geom_fruit(data = df_blank_wild, geom = geom_col,
             mapping = aes(y = host_taxon, x = blank),
             orientation = "y", width = 0.9, pwidth = 0.19, offset = 0.25,
             fill = "white", color = NA) +
  scale_fill_identity() +
  
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.4, "cm"),
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 14,
                                     margin = margin(b = 18)),
    plot.margin       = margin(t = 10, r = 10, b = 25, l = 10)
  ) +
  ggtitle("Wild Herpetofauna Gut Microbiomes (Fecal and Lower GI Samples Only)")

print(p_tree_wild)

ggsave(filename = file.path(dir_figures, "host_tree_microbiome_rel_abundance_wild_fecal_lowerGI.pdf"), plot = p_tree_wild, width = 10, height = 10, units = "in", device = cairo_pdf, dpi = 750)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
# DSpikeIn Script - 16S
# Converting relative abundance microbiome data to absolute abundance 
# Alexander Rurik
###############################################################################

# Load required packages #

# Function to check and load required packages
load_pkgs <- function(pkgs, type = "CRAN") {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing) > 0) {
    warning(
      paste0(
        "The following ", type, " packages are not installed:\n",
        paste(missing, collapse = ", "), 
        "\nInstall them before running this script."
      )
    )
  }
  invisible(lapply(pkgs, require, character.only = TRUE))
}

# ----------------------------- CRAN PACKAGES ---------------------------------

cran_pkgs <- c(
  "stats", "dplyr", "ggplot2", "flextable", "ggpubr", "randomForest",
  "ggridges", "ggalluvial", "tibble", "matrixStats", "RColorBrewer",
  "ape", "rlang", "scales", "magrittr", "phangorn", "igraph", "tidyr",
  "xml2", "data.table", "reshape2", "vegan", "patchwork", "officer",
  "stringr"
)

load_pkgs(cran_pkgs, type = "CRAN")

# -------------------------- BIOCONDUCTOR PACKAGES ----------------------------

bioc_pkgs <- c(
  "phyloseq", "msa", "DESeq2", "edgeR", "Biostrings", "ggtree",
  "DECIPHER", "microbiome", "limma", "S4Vectors",
  "SummarizedExperiment", "TreeSummarizedExperiment"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  warning("BiocManager not installed. Install using install.packages('BiocManager')")
} else {
  load_pkgs(bioc_pkgs, type = "Bioconductor")
}

# ---------------------------- GITHUB PACKAGES --------------------------------

github_pkgs <- c("speedyseq", "microbiomeutilities", "DspikeIn")

# Warn if missing
load_pkgs(github_pkgs, type = "GitHub")

# ---------------------------------- DONE -------------------------------------
message("All available packages loaded.")


##------------------------------------------------------------------------------
## Begin DspikeIn Run ####
##------------------------------------------------------------------------------


# ----------------------------- Load additional packages ----------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggpubr)
  library(data.table)
  library(decontam)
  library(vegan)
  library(ggtext)
  library(phyloseq)
  library(stringr)
})

# ----------------------------- Prioritize dplyr functions --------------------
# Ensures select, filter, rename use dplyr versions
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# ----------------------------- Set basic parameters --------------------------
set.seed(52325)                        # Set seed for reproducibility

# ----------------------------- Project directories ----------------------------
# here() anchors all paths to the project root, so this script runs identically
# regardless of whose machine or OS it's run on, as long as the repo's folder
# structure (data/raw, data/processed, output/figures, output/tables) is intact.
library(here)
dir_raw       <- here("data", "raw", "16S")
dir_processed <- here("data", "processed")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")

# ----------------------------- Load tables -----------------------------------
# Read in OTU table, taxonomy, and sample metadata
otu.table <- fread(file.path(dir_raw, "16S_count_table.csv"))
tax.table <- fread(file.path(dir_raw, "16S_taxonomy.csv"))
metadata  <- fread(file.path(dir_raw, "16S_metadata_final.csv"))

# ----------------------------- Clean these files -------------------------------
# OTU table: move OTU ID column to rownames; rows = OTUs, columns = samples
otu.table <- otu.table %>%
  column_to_rownames("OTU")

# Taxonomy: standardize taxonomy strings and set OTU IDs as rownames
tax.table <- tax.table %>%
  mutate(across(Kingdom:Species, 
                ~ ifelse(is.na(.) | . == "", "", str_remove(., "^[a-z]__") %>% str_to_title()))) %>%
  column_to_rownames('OTU')
colnames(tax.table) # Inspect cleaned taxonomy columns

# Metdata: arrange samples and convert relevant columns to factors
metadata.2 <- metadata %>% 
  arrange(sample_name) %>% 
  mutate(across(all_of(c('plate_ID','host_taxon', 'host_genus', 'host_species', 'state', 'site')), 
                as.factor)) %>%
  column_to_rownames('sample_name')

# ----------------------------- Build phyloseq object -------------------------
# Create OTU, taxonomy, and metadata components
OTU <- otu_table(as.matrix(otu.table), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax.table))
metadata.phyloseq <- sample_data(metadata.2)

# Merge components into a single phyloseq object
dspikein.PS <- phyloseq(OTU, TAX, metadata.phyloseq)
dspikein.PS  

# Remove samples from unspiked STP plate
dspikein.PS <- subset_samples(dspikein.PS, plate_ID != "220220_STP_P1")
dspikein.PS

# Remove OTUs that are absent across all samples
dspikein.PS <- prune_taxa(taxa_sums(dspikein.PS) > 0, dspikein.PS)
dspikein.PS

# Check spike-in counts
unique(sample_data(dspikein.PS)$Tetra_cell_count)

# ---- Check read counts after removing STP plate ----

# Total reads across entire dataset
total_reads <- sum(otu_table(dspikein.PS))
cat("Total reads in dataset:", format(total_reads, big.mark = ","), "\n")

# Number of samples
n_samples <- nsamples(dspikein.PS)
cat("Number of samples:", n_samples, "\n")

# Number of OTUs
n_otus <- ntaxa(dspikein.PS)
cat("Number of OTUs:", n_otus, "\n")

##------------------------------------------------------------------------------
## DspikeIn Pre-processing ####
##------------------------------------------------------------------------------
# Rename column in metadata to match DspikeIn expectations
sample.data <- microbiome::meta(dspikein.PS)
sample.data$spiked.volume <- sample.data$spiked_volume 
sample_data(dspikein.PS) <- sample_data(sample.data)

# Convert phyloseq object to tidy format compatible with DspikeIn
physeq_16SOTU <- tidy_phyloseq_tse(dspikein.PS)
physeq_16SOTU@sam_data$spiked.volume # Quick check: ensure metadata contains spiked volumes

# ----------------------------- Define spiked species -------------------------
species_name <- spiked_species <- c("Tetragenococcus")  # Spiked species genus
merged_spiked_species <- "Tetragenococcus"

# Select OTUs corresponding to the spiked species
Tetra <- subset_taxa(physeq_16SOTU, Genus == "Tetragenococcus")
hashcodes <- row.names(phyloseq::tax_table(Tetra))

# Keep all samples (including blanks) for decontamination purposes
spiked_16S_OTU <- tidy_phyloseq_tse(physeq_16SOTU)

# ------------------ Preprocessing: One Species Scaling Factor ----------------
# Merge OTUs/ASVs of the spiked species and calculate scaling factor
# merge_method = "max": retains the most abundant ASV
# merge_method = "sum": sums abundances across all OTUs/ASVs

# If rerunning analyses, can load processed object directly instead of re-merging:
# Spiked_16S_sum_scaled <- readRDS(file.path(dir_processed, "16S_spike_merged_physeq_sum.rds"))

# Merge hashcodes (specific for QIIME2/DADA2 outputs)
Spiked_16S_sum_scaled <- Pre_processing_hashcodes(
  spiked_16S_OTU, 
  hashcodes = hashcodes, 
  merge_method = "sum", 
  output_prefix = file.path(dir_processed, "16S_spike_merged_physeq_sum")
)

# ----------------------------- Identify merged spike taxon -------------------
# Confirm that the merged spiked species is present
ta_df <- as.data.frame(tax_table(Spiked_16S_sum_scaled))
tetra_idx <- which(ta_df$Genus == "Tetragenococcus")
if (length(tetra_idx) == 0) stop("No Tetragenococcus taxa found in Spiked_16S_sum_scaled")

# Extract taxa names and sums for the spiked species
tetra_names <- taxa_names(Spiked_16S_sum_scaled)[tetra_idx]
tetra_sums <- taxa_sums(Spiked_16S_sum_scaled)[tetra_names]

# Pick the top carrier (most abundant merged taxon)
merged_spiked_hashcodes <- names(sort(tetra_sums, decreasing = TRUE))[1]
message("Using merged spike taxon: ", merged_spiked_hashcodes)

##------------------------------------------------------------------------------
## Calculate Spiked Species Retrieval Percentage ####
##------------------------------------------------------------------------------
## This section summarizes OTU counts and calculates the spiked species recovery
## percentage for the merged spiked taxon (Tetragenococcus). 

# Convert merged phyloseq object to tidy format for DspikeIn functions
Spiked_16S_OTU_scaled <- tidy_phyloseq_tse(Spiked_16S_sum_scaled)

# ----------------------------- Calculate spike percentage --------------------
# Customize the threshold (passed_range) and spiked species/hashcodes as needed
# passed_range = c(x, y): acceptable range (%) for spike recovery
result <- calculate_spike_percentage(
  obj = Spiked_16S_sum_scaled,
  merged_spiked_hashcodes = merged_spiked_hashcodes,  # use the merged top spike hashcode
  passed_range = c(0, 100)                            # adjust to desired range
)

##-----------------------------------------------------------------------------------
## Alpha & Beta Diversity Analysis / Examining System-specific Spike-in Cutoffs ####
##-----------------------------------------------------------------------------------
## This section calculates alpha diversity metrics, merges them with metadata,
## computes beta diversity (distance to global centroid), and generates regression
## plots to evaluate spike-in recovery and community diversity.

# ----------------------------- Alpha Diversity Metrics -----------------------
# Calculate Observed richness and Shannon diversity
alphab <- estimate_richness(Spiked_16S_sum_scaled, measures = c("Observed", "Shannon"))

# Fix sample names to match metadata / OTU table (replace "." with "-")
rownames(alphab) <- str_replace_all(rownames(alphab), "\\.", "-")
alphab$Sample <- rownames(alphab)

# Calculate Pielou's Evenness
alphab$Pielou_evenness <- alphab$Shannon / alphab$Observed

# Subset OTU matrix to match alpha diversity samples
otu_mat <- as.matrix(otu_table(Spiked_16S_sum_scaled))
if (taxa_are_rows(Spiked_16S_sum_scaled)) {
  otu_mat <- t(otu_mat)
}
otu_mat_subset <- otu_mat[rownames(alphab), , drop = FALSE]

# Calculate Hill number q = 1 (exp of Shannon) on matching samples
alphab$Hill_q1 <- exp(vegan::diversity(otu_mat_subset, index = "shannon"))

# --------------------- Merge Alpha Diversity with Metadata ------------------
metadata.dspikein <- as.data.frame(microbiome::meta(Spiked_16S_sum_scaled))
metadata.dspikein$Sample <- rownames(metadata.dspikein)

# Merge alpha metrics and spike-in statistics
metadata.4 <- metadata.dspikein %>%
  left_join(alphab[, c("Sample", "Observed", "Shannon", "Pielou_evenness", "Hill_q1")], by = "Sample") %>%
  left_join(result[, c("Sample", "Spiked_Reads", "Percentage")], by = "Sample") %>% 
  filter(sample_or_blank == 'sample') %>%
  column_to_rownames(var = "Sample")

# ----------------------------- Beta Diversity: Distance to Centroid ---------
# Convert OTU matrix to relative abundances
otu_mat_rel <- vegan::decostand(otu_mat_subset, method = "total") 
# Compute the centroid profile across samples
centroid_profile <- colMeans(otu_mat_rel) 
# Compute Bray-Curtis distance from each sample to the centroid
dist_to_centroid <- apply(otu_mat_rel, 1, function(x) { 
  vegan::vegdist(rbind(x, centroid_profile), method = "bray")[1]
})

# Add distance metric to metadata
metadata.4$Dist_to_Centroid <- dist_to_centroid[rownames(metadata.4)]
# Update phyloseq object with new metadata
sample_data(Spiked_16S_sum_scaled) <- sample_data(metadata.4)

# Quick check: ensure spike reads column exists
if (!"Spiked_Reads" %in% colnames(metadata.4)) {
  stop("Column 'Spiked_Reads' not found in metadata.")
}

# ----------------------------- Merge metadata with spike-in results ----------
# Move first column of result to rownames (or ensure rownames alignment)
result_aligned <- result %>%
  column_to_rownames(var = "Sample")

# Move rownames of metadata.4 to a column called "Sample"
metadata_with_sample <- metadata.4 %>%
  rownames_to_column(var = "Sample")

# Identify common samples between metadata.4 and result
common_samples <- intersect(rownames(metadata.4), rownames(result_aligned))

# Create combined dataframe with metadata and result columns
combined_data <- metadata_with_sample %>%
  left_join(
    result_aligned %>% rownames_to_column(var = "Sample"),
    by = "Sample"
  ) %>%
  filter(Sample %in% common_samples)

# Save combined metadata + spike-in results
write.csv(combined_data, file = file.path(dir_tables, "16S_metadata_spike_results_combined.csv"), row.names = FALSE)

message("Combined metadata and spike-in data saved to: ",
        file.path(dir_tables, "16S_metadata_spike_results_combined.csv"))

# ----------------------------- Generate Regression Plots ---------------------
# Pielou's Evenness
p1 <- regression_plot(
  data = metadata.4,
  y_var = "Pielou_evenness",
  x_var = "Spiked_Reads",
  custom_range = c(0.01, 10, 20, 30, 40, 50, 100),
  plot_title = "Pielou_evenness vs Spike-in Reads"
)
# Save
pdf_filename <- file.path(dir_figures, "16S_PielouEvenness_vs_SpikeReads.pdf")
pdf(pdf_filename, width = 12, height = 7)  
print(p1)
dev.off()

# ----------------------------- Hill number q = 1 -----------------------------
p2 <- regression_plot(
  data = metadata.4,
  y_var = "Hill_q1",
  x_var = "Spiked_Reads",
  custom_range = c(0.01, 20, 30, 40, 100),
  plot_title = "Hill number vs Spike-in Reads"
)
# Save
pdf_filename <- file.path(dir_figures, "16S_HillNumber_vs_SpikeReads.pdf")
pdf(pdf_filename, width = 12, height = 7)
print(p2)
dev.off()

# ----------------------------- Distance to Global Centroid --------------------
p3 <- regression_plot(
  data = metadata.4,
  y_var = "Dist_to_Centroid",
  x_var = "Spiked_Reads",
  custom_range = c(0.01, 20, 30, 40, 100),
  plot_title = "Distance to Global Centroid vs Spike-in Reads"
)
# Save
pdf_filename <- file.path(dir_figures, "16S_DistToCentroid_vs_SpikeReads.pdf")
pdf(pdf_filename, width = 12, height = 7)
print(p3)
dev.off()

##-----------------------------------------------------------------------------------
## Create Multipanel Figure for Spike-in Analysis (with smaller text) ####
##-----------------------------------------------------------------------------------

library(cowplot)  # For plot_grid

# Generate the three plots with panel labels in titles and smaller text
p1 <- regression_plot(
  data = metadata.4,
  y_var = "Pielou_evenness",
  x_var = "Spiked_Reads",
  custom_range = c(0.01, 10, 20, 30, 40, 50, 100),
  plot_title = "A) Pielou's Evenness vs Spike-in Reads"
) + 
  theme(
    text = element_text(size = 9),           # Overall text size
    axis.title = element_text(size = 10),    # Axis labels
    axis.text = element_text(size = 9),      # Axis tick labels
    plot.title = element_text(size = 11),    # Plot title
    legend.text = element_text(size = 9),    # Legend text
    legend.title = element_text(size = 9),   # Legend title
    strip.text = element_text(size = 9)      # Facet labels
  )

p2 <- regression_plot(
  data = metadata.4,
  y_var = "Hill_q1",
  x_var = "Spiked_Reads",
  custom_range = c(0.01, 20, 30, 40, 100),
  plot_title = "B) Hill Number vs Spike-in Reads"
) + 
  theme(
    text = element_text(size = 9),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    strip.text = element_text(size = 9)
  )

p3 <- regression_plot(
  data = metadata.4,
  y_var = "Dist_to_Centroid",
  x_var = "Spiked_Reads",
  custom_range = c(0.01, 20, 30, 40, 100),
  plot_title = "C) Distance to Global Centroid vs Spike-in Reads"
) + 
  theme(
    text = element_text(size = 9),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    strip.text = element_text(size = 8)
  )

# Extract the shared legend from p2 (or p3, they should be the same)
shared_legend <- get_legend(p2)

# Remove legends from p2 and p3
p2_no_legend <- p2 + theme(legend.position = "none")
p3_no_legend <- p3 + theme(legend.position = "none")

# Combine p2 and p3 horizontally with the shared legend
bottom_row <- plot_grid(
  p2_no_legend, 
  p3_no_legend, 
  shared_legend,
  ncol = 3,
  rel_widths = c(1, 1, 0.3)  # Adjust legend width as needed
)

# Combine top panel (p1) with bottom row
multipanel <- plot_grid(
  p1,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 0.8)  # Adjust relative heights as needed
)

# Save the multipanel figure
pdf_filename <- file.path(dir_figures, "16S_SpikeinAnalysis_Multipanel.pdf")
pdf(pdf_filename, width = 12, height = 13)
print(multipanel)
dev.off()

###########################################################
### Based on these plots, selecting a 0.01%-40% cutoff    #
### Failed reads will be removed during the decontam step #
###########################################################

##------------------------------------------------------------------------------
## Selecting Spike-in Cutoffs and Estimating Scaling Factors ####
##------------------------------------------------------------------------------

# ----------------------------- Split Samples by Spike Level ------------------
meta_16S <- data.frame(sample_data(Spiked_16S_OTU_scaled))

meta_16S <- meta_16S %>%
  mutate(spike_group = case_when(
    Tetra_cell_count %in% c(1847, 924, 0) ~ "full_half_spike",
    Tetra_cell_count == 154 ~ "one_twelfth_spike"
  ))

# Subset phyloseq objects by spike group
Spiked_16S_fullhalf <- prune_samples(
  rownames(meta_16S[meta_16S$spike_group == "full_half_spike", ]),
  Spiked_16S_OTU_scaled
)
Spiked_16S_one12th <- prune_samples(
  rownames(meta_16S[meta_16S$spike_group == "one_twelfth_spike", ]),
  Spiked_16S_OTU_scaled
)

# ----------------------------- Assign Spiked Cells per Well -----------------
# Define the number of spiked cells corresponding to each spike level
spiked_cells_fullhalf <- 1847 # for full and half spikes (1847 and 924)
spiked_cells_one12th <- 308   # for 1/12 spikes (154 cells in 1 µl spike volume. So use 308 cells for 1/12 spike, since listed as 1ul spike volume)

# ----------------------------- Calculate Spike-in Scaling Factors -----------
result_fullhalf_16S <- calculate_spikeIn_factors(
  obj = Spiked_16S_fullhalf,
  spiked_cells = spiked_cells_fullhalf,
  merged_spiked_species = merged_spiked_species
)
result_one12th_16S <- calculate_spikeIn_factors(
  obj = Spiked_16S_one12th,
  spiked_cells = spiked_cells_one12th,
  merged_spiked_species = merged_spiked_species
)

# Extract scaling factors
scaling_factors_fullhalf_16S <- result_fullhalf_16S$scaling_factors
scaling_factors_one12th_16S  <- result_one12th_16S$scaling_factors

# Optional: view scaling factors
view(scaling_factors_fullhalf_16S)
view(scaling_factors_one12th_16S)

# ------------ Convert Relative Counts to Absolute Counts ------------------
# Convert relative abundance to absolute counts using scaling factors
abs_fullhalf_16S <- convert_to_absolute_counts(
  obj = Spiked_16S_fullhalf,
  scaling_factors = scaling_factors_fullhalf_16S
)

abs_one12th_16S <- convert_to_absolute_counts(
  obj = Spiked_16S_one12th,
  scaling_factors = scaling_factors_one12th_16S
)

# Store the adjusted phyloseq objects with absolute counts
physeq_absolute_fullhalf_16S <- abs_fullhalf_16S$obj_adj
physeq_absolute_one12th_16S  <- abs_one12th_16S$obj_adj

# ----------------------------- Quick Verification ---------------------------
# Check that the sum of a specific Tetragenococcus OTU matches expectations
taxa_sums(physeq_absolute_fullhalf_16S)["2ddb215ff668b6a24a9ccf2bb076a453"] +
  taxa_sums(physeq_absolute_one12th_16S)["2ddb215ff668b6a24a9ccf2bb076a453"]
# The absolute abundance estimate after applying the spike-in scaling factors

##------------------------------------------------------------------------------
## Merge spike-level phyloseq objects and export ####
##------------------------------------------------------------------------------

# Merge full/half and one-twelfth spike phyloseq objects
physeq_absolute_16S <- merge_phyloseq(physeq_absolute_fullhalf_16S,
                                      physeq_absolute_one12th_16S)

# Extract metadata and add alpha/beta diversity and spike info
metadata_final_16S <- as(sample_data(physeq_absolute_16S), "data.frame") %>%
  rownames_to_column("sample_name") %>%
  left_join(alphab[, c("Sample", "Observed", "Shannon", "Pielou_evenness", "Hill_q1")],
            by = c("sample_name" = "Sample")) %>%
  left_join(result[, c("Sample", "Spiked_Reads", "Percentage")],
            by = c("sample_name" = "Sample"))

# Put enriched metadata back into the phyloseq object
sample_data(physeq_absolute_16S) <- sample_data(metadata_final_16S %>%
                                                  column_to_rownames("sample_name"))

# Export count and taxonomy tables
count_table_16S <- as.data.frame(as(otu_table(physeq_absolute_16S), "matrix")) %>%
  rownames_to_column("OTU")
taxonomy_16S <- as.data.frame(as(tax_table(physeq_absolute_16S), "matrix")) %>%
  rownames_to_column("OTU")

# Save RDS and CSVs
saveRDS(physeq_absolute_16S, file = file.path(dir_processed, "physeq_absolute_16S.rds"))
write.csv(count_table_16S, file = file.path(dir_tables, "16S_abs_abund_otu_table.csv"), row.names = FALSE)
write.csv(taxonomy_16S,   file = file.path(dir_tables, "16S_abs_abund_taxonomy.csv"), row.names = FALSE)
write.csv(metadata_final_16S, file = file.path(dir_tables, "16S_abs_abund_metadata.csv"), row.names = FALSE)

# Done. Can take the phyloseq object and move to decontam 
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
###############################################################################
# Supplementary tables: top microbial taxa by host taxonomic group (16S & ITS)
#
# Produces a CSV table showing the top N bacterial and fungal taxa
# for each host group at a specified taxonomic level, for a specified 
# dataset type (in this case, all samples, or fecal or cloacal swabs only.
#
# Parameters to modify before running:
#   DATASET_TYPE     — "all", "fecal", or "cloacal"
#   HOST_LEVEL       — "Order", "Family", "Genus", or "Species"
#   MICROBE_TAX_LEVEL — "Phylum", "Class", "Order", "Family", or "Genus"
#   N_TOP_TAXA       — number of top taxa to report per host group
#
# Alexander Rurik
###############################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
})

# ========================== BASIC PARAMETERS ==================================
set.seed(52325)

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_tables    <- here("output", "tables")

# ========================== USER-DEFINED PARAMETERS ==========================
DATASET_TYPE      <- "fecal"  # "all", "fecal", or "cloacal"
HOST_LEVEL        <- "Order"  # "Order", "Family", "Genus", or "Species"
MICROBE_TAX_LEVEL <- "Phylum"  # "Phylum", "Class", "Order", "Family", or "Genus"
N_TOP_TAXA        <- 5

# ========================== LOAD DATA =========================================
if (DATASET_TYPE == "all") {
  ps_16S <- readRDS(file.path(dir_processed, "16S_absolute_all_final.rds"))
  ps_ITS <- readRDS(file.path(dir_processed, "ITS_absolute_all_final.rds"))
} else if (DATASET_TYPE == "fecal") {
  ps_16S <- readRDS(file.path(dir_processed, "16S_abs_final_fecal.rds"))
  ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_fecal.rds"))
} else if (DATASET_TYPE == "cloacal") {
  ps_16S <- readRDS(file.path(dir_processed, "16S_abs_final_cloacal.rds"))
  ps_ITS <- readRDS(file.path(dir_processed, "ITS_abs_final_cloacal.rds"))
} else {
  stop("Invalid DATASET_TYPE. Choose 'all', 'fecal', or 'cloacal'.")
}

cat("Loaded:", DATASET_TYPE, "dataset\n")
cat("Host level:", HOST_LEVEL, "\n")
cat("Microbial level:", MICROBE_TAX_LEVEL, "\n\n")

# ========================== STANDARDIZE METADATA COLUMNS =====================
# Renames phyloseq sample_data columns to standardized names used throughout
# this script. Preserves rownames so sample_data assignment doesn't break.

standardize_host_columns <- function(ps_obj) {
  sdata <- data.frame(sample_data(ps_obj))
  rn    <- rownames(sdata)  # preserve sample name rownames
  
  rename_map <- c(
    "Clade_Order" = "Order",
    "host_genus"  = "Genus",
    "host_taxon"  = "Species"
  )
  for (old in names(rename_map)) {
    if (old %in% colnames(sdata)) {
      colnames(sdata)[colnames(sdata) == old] <- rename_map[[old]]
    }
  }
  # Family column is expected to already be named "host_family" or "Family";
  # add a rename here if your metadata uses a different name.
  
  rownames(sdata) <- rn
  sample_data(ps_obj) <- sample_data(sdata)
  return(ps_obj)
}

ps_16S <- standardize_host_columns(ps_16S)
ps_ITS <- standardize_host_columns(ps_ITS)

# ========================== HELPER FUNCTIONS ==================================

# Returns the plural form of a microbial taxonomic rank name for column headers
pluralize_taxon <- function(tax_level) {
  plurals <- c(
    Phylum  = "Phyla",
    Class   = "Classes",
    Order   = "Orders",
    Family  = "Families",
    Genus   = "Genera",
    Species = "Species"
  )
  if (tax_level %in% names(plurals)) plurals[[tax_level]] else paste0(tax_level, "s")
}

# Returns top N microbial taxa by summed absolute abundance for each host group.
# Abundances are converted to percentages relative to the total within each
# host group (not across the full dataset).
get_top_taxa_by_host_level <- function(ps_obj, microbe_tax_level, host_level, n_top = 5) {
  
  if (is.null(ps_obj) || nsamples(ps_obj) == 0) return(NULL)
  
  # Agglomerate to target microbial rank; NArm = FALSE retains unclassified taxa
  ps_glom <- tax_glom(ps_obj, taxrank = microbe_tax_level, NArm = FALSE)
  
  sample_df        <- data.frame(sample_data(ps_glom))
  sample_df$SampleID <- rownames(sample_df)
  
  otu_df <- data.frame(otu_table(ps_glom))
  tax_df <- data.frame(tax_table(ps_glom), stringsAsFactors = FALSE)
  
  # Ensure samples are rows, taxa are columns
  if (taxa_are_rows(ps_glom)) otu_df <- t(otu_df)
  
  otu_df$SampleID <- rownames(otu_df)
  merged_df       <- left_join(otu_df, sample_df, by = "SampleID")
  
  taxon_cols  <- seq_len(ncol(otu_df) - 1)  # all columns except SampleID
  host_groups <- na.omit(unique(merged_df[[host_level]]))
  
  results_list <- list()
  
  for (host_group in host_groups) {
    host_samples <- merged_df[merged_df[[host_level]] == host_group, ]
    total_reads  <- colSums(host_samples[, taxon_cols, drop = FALSE], na.rm = TRUE)
    total_sum    <- sum(total_reads)
    percentages  <- if (total_sum > 0) (total_reads / total_sum) * 100 else total_reads * 0
    
    top_indices   <- order(total_reads, decreasing = TRUE)[seq_len(min(n_top, length(total_reads)))]
    top_taxa_ids  <- names(total_reads)[top_indices]
    top_taxa_names <- as.character(tax_df[top_taxa_ids, microbe_tax_level])
    top_taxa_perc  <- percentages[top_taxa_ids]
    
    # Replace missing or blank taxonomy with "Unknown"
    top_taxa_names[is.na(top_taxa_names) | str_trim(top_taxa_names) == ""] <- "Unknown"
    
    formatted_taxa <- paste0(top_taxa_names, " (", round(top_taxa_perc, 1), "%)")
    results_list[[host_group]] <- paste(formatted_taxa, collapse = "\n")
  }
  
  return(results_list)
}

# ========================== MAIN ANALYSIS =====================================
cat("Calculating top", N_TOP_TAXA, "bacterial taxa at", MICROBE_TAX_LEVEL, "level...\n")
top_16S <- get_top_taxa_by_host_level(ps_16S, MICROBE_TAX_LEVEL, HOST_LEVEL, N_TOP_TAXA)

cat("Calculating top", N_TOP_TAXA, "fungal taxa at", MICROBE_TAX_LEVEL, "level...\n")
top_ITS <- get_top_taxa_by_host_level(ps_ITS, MICROBE_TAX_LEVEL, HOST_LEVEL, N_TOP_TAXA)

# Build host group metadata table, including Order for sorting (unless HOST_LEVEL == "Order")
host_info_16S <- data.frame(sample_data(ps_16S)) %>%
  { if (HOST_LEVEL == "Order") select(., all_of(HOST_LEVEL))
    else select(., all_of(c(HOST_LEVEL, "Order"))) } %>%
  distinct()

host_info_ITS <- data.frame(sample_data(ps_ITS)) %>%
  { if (HOST_LEVEL == "Order") select(., all_of(HOST_LEVEL))
    else select(., all_of(c(HOST_LEVEL, "Order"))) } %>%
  distinct()

# Sample counts per host group for each dataset
sample_counts_16S <- data.frame(sample_data(ps_16S)) %>%
  group_by(across(all_of(HOST_LEVEL))) %>%
  summarise(n_16S = n(), .groups = "drop")

sample_counts_ITS <- data.frame(sample_data(ps_ITS)) %>%
  group_by(across(all_of(HOST_LEVEL))) %>%
  summarise(n_ITS = n(), .groups = "drop")

# Combine host info from both datasets, attach sample counts
host_info <- bind_rows(host_info_16S, host_info_ITS) %>%
  distinct() %>%
  filter(!is.na(.data[[HOST_LEVEL]])) %>%
  left_join(sample_counts_16S, by = HOST_LEVEL) %>%
  left_join(sample_counts_ITS, by = HOST_LEVEL) %>%
  mutate(
    n_16S = replace_na(n_16S, 0),
    n_ITS = replace_na(n_ITS, 0)
  )

# Sort: amphibians first (Anura, Caudata), then reptiles (Crocodilia, Squamata, Testudines)
if ("Order" %in% colnames(host_info)) {
  host_info$Order <- factor(host_info$Order,
                            levels = c("Anura", "Caudata", "Crocodilia", "Squamata", "Testudines"))
  host_info <- arrange(host_info, Order, .data[[HOST_LEVEL]])
} else {
  host_info <- arrange(host_info, .data[[HOST_LEVEL]])
}

# ========================== BUILD FINAL TABLE =================================
microbe_tax_plural <- pluralize_taxon(MICROBE_TAX_LEVEL)
col_name_host      <- paste0("Host ", HOST_LEVEL, " (# samples: 16S/ITS)")
col_name_bacteria  <- paste0("Top ", N_TOP_TAXA, " Bacterial ", microbe_tax_plural)
col_name_fungi     <- paste0("Top ", N_TOP_TAXA, " Fungal ",    microbe_tax_plural)

# Safe lookup: uses match() instead of [[ ]] to avoid "subscript out of bounds"
# when genus names contain special characters or whitespace
lookup_taxa <- function(host_groups, taxa_list) {
  vapply(host_groups, function(hg) {
    idx <- match(hg, names(taxa_list))
    if (!is.na(idx)) taxa_list[[idx]] else "NA"
  }, character(1))
}

final_table <- host_info %>%
  mutate(
    !!col_name_host     := paste0(.data[[HOST_LEVEL]], " (", n_16S, "/", n_ITS, ")"),
    !!col_name_bacteria := lookup_taxa(.data[[HOST_LEVEL]], top_16S),
    !!col_name_fungi    := lookup_taxa(.data[[HOST_LEVEL]], top_ITS)
  ) %>%
  select(all_of(c(col_name_host, col_name_bacteria, col_name_fungi)))

# ========================== SAVE OUTPUT =======================================
output_filename <- file.path(dir_tables, paste0(DATASET_TYPE, "_top_taxa_by_host_",
                          HOST_LEVEL, "_microbe_", MICROBE_TAX_LEVEL, ".csv"))
write.csv(final_table, output_filename, row.names = FALSE)

cat("\nAnalysis complete!\n")
cat("Dataset:          ", DATASET_TYPE, "\n")
cat("Host level:       ", HOST_LEVEL, "\n")
cat("Microbial level:  ", MICROBE_TAX_LEVEL, "\n")
cat("Output:           ", output_filename, "\n")
cat("Host groups:      ", nrow(final_table), "\n")
cat("Groups with 16S:  ", sum(final_table[[col_name_bacteria]] != "NA"), "\n")
cat("Groups with ITS:  ", sum(final_table[[col_name_fungi]]    != "NA"), "\n")

###############################################################################
###############################################################################
###############################################################################
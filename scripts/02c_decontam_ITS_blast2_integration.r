###############################################################################
# ITS BLAST Integration & Final Dataset Export — Script 2 of 2
#
# Picks up from Script 1 (decontam_ITS_blast1.r) after the BLAST search
# has been run in terminal against the UNITE+INSD all-eukaryotes database.
#
# Pipeline overview:
#   Part 1:  Load pre-BLAST phyloseq + BLAST results
#   Part 2:  Parse BLAST taxonomy from sseqid field
#   Part 3:  Apply keep/remove logic and track filter reason per OTU
#   Part 4:  Validate filtering — confirm no non-fungal OTUs remain
#   Part 5:  Update taxonomy table with BLAST-rescued assignments
#   Part 6:  Build and export final fungi-only phyloseq object
#   Part 7:  Full validation of final phyloseq object
#
# Required inputs:
#   ITS_physeq_absolute_preblast.rds   (output of Script 1)
#   blastn_all_results.tsv             (output of terminal BLAST run; place in data/processed)
#
# Alexander Rurik
###############################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(Biostrings)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# ========================== BASIC PARAMETERS ==================================
set.seed(52325)

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_tables    <- here("output", "tables")

##------------------------------------------------------------------------------
## Part 1: Load data ####
##------------------------------------------------------------------------------

# --- Load pre-BLAST phyloseq object from Script 1 ---
ps_preblast <- readRDS(file.path(dir_processed, "ITS_physeq_absolute_preblast.rds"))

cat("========== PRE-BLAST DATASET (from Script 1) ==========\n")
cat("Samples:", nsamples(ps_preblast), "\n")
cat("OTUs:", ntaxa(ps_preblast), "\n")
cat("Total reads:", format(sum(sample_sums(ps_preblast)), big.mark = ","), "\n\n")

# --- Load BLAST results ---
# This .tsv is the output of the terminal BLAST run (see note at end of script)
blast_results <- read.delim(file.path(dir_processed, "ITS_blastn_all_results.tsv"),
                            header           = FALSE,
                            stringsAsFactors = FALSE)

colnames(blast_results) <- c("OTU_ID", "sseqid", "pident", "qcov", "evalue", "bitscore")

# Keep only the top hit per OTU (highest bitscore; BLAST -max_target_seqs 5 returns up to 5)
blast_results <- blast_results %>%
  group_by(OTU_ID) %>%
  slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("========== BLAST RESULTS ==========\n")
cat("OTUs with BLAST results:", nrow(blast_results), "\n")
cat("OTUs with a hit:", sum(!is.na(blast_results$sseqid)), "\n")
cat("OTUs without a hit:", sum(is.na(blast_results$sseqid)), "\n\n")

##------------------------------------------------------------------------------
## Part 2: Parse BLAST taxonomy ####
##------------------------------------------------------------------------------
# UNITE+INSD sseqid format: k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
# Extract Kingdom, Phylum, and Class for use in filtering and taxonomy update

blast_results <- blast_results %>%
  mutate(
    BLAST_Kingdom = ifelse(!is.na(sseqid),
                           str_extract(sseqid, "k__([^;]+)") %>% str_remove("k__"), NA),
    BLAST_Phylum  = ifelse(!is.na(sseqid),
                           str_extract(sseqid, "p__([^;]+)") %>% str_remove("p__"), NA),
    BLAST_Class   = ifelse(!is.na(sseqid),
                           str_extract(sseqid, "c__([^;]+)") %>% str_remove("c__"), NA)
  )

cat("========== BLAST KINGDOM ASSIGNMENTS ==========\n")
print(table(blast_results$BLAST_Kingdom, useNA = "always"))

##------------------------------------------------------------------------------
## Part 3: Apply keep/remove logic ####
##------------------------------------------------------------------------------
# Decision rules:
#
# KEEP if:
#   (a) Original UNITE Kingdom = "Fungi" AND BLAST does not contradict
#       (i.e., BLAST is also Fungi, or no BLAST hit)
#   (b) BLAST Kingdom = "Fungi" with high-quality alignment
#       (pident >= 80 AND qcov >= 70) — BLAST-rescued OTU
#
# REMOVE if:
#   (a) BLAST explicitly assigns a non-fungal kingdom (regardless of UNITE assignment)
#   (b) OTU is unassigned or incertae sedis AND BLAST either fails quality thresholds
#       or returns no hit
#
# Note: The quality thresholds (pident >= 80, qcov >= 70) are consistent with
# those used in the original analysis. Adjust here if revisiting.

# Extract original taxonomy from phyloseq
tax_orig        <- as.data.frame(tax_table(ps_preblast))
tax_orig$OTU_ID <- rownames(tax_orig)

# Merge original taxonomy with BLAST results
tax_combined <- tax_orig %>%
  left_join(
    blast_results %>% select(OTU_ID, BLAST_Kingdom, BLAST_Phylum, BLAST_Class,
                             pident, qcov, sseqid),
    by = "OTU_ID"
  )

# Apply filtering logic
tax_combined <- tax_combined %>%
  mutate(
    Keep = case_when(
      # REMOVE: BLAST explicitly says non-fungal (overrides UNITE)
      !is.na(BLAST_Kingdom) & BLAST_Kingdom != "Fungi"                        ~ FALSE,
      
      # KEEP: UNITE assigned as Fungi and BLAST did not contradict
      Kingdom == "Fungi" & Phylum != "Fungi_phy_incertae_sedis"               ~ TRUE,
      
      # KEEP: BLAST rescued — Fungi_phy_incertae_sedis or unassigned,
      #       confirmed fungal by BLAST with good alignment quality
      BLAST_Kingdom == "Fungi" & pident >= 80 & qcov >= 70                    ~ TRUE,
      
      # REMOVE: Everything else (unassigned/incertae sedis that failed or had no BLAST hit)
      TRUE                                                                      ~ FALSE
    ),
    
    Filter_Reason = case_when(
      # UNITE said Fungi but BLAST contradicts
      Kingdom == "Fungi" & !is.na(BLAST_Kingdom) & BLAST_Kingdom != "Fungi"  ~
        paste0("UNITE_Fungi_BLAST_says_", BLAST_Kingdom),
      
      # Straightforward UNITE Fungi, no BLAST hit
      Kingdom == "Fungi" & Phylum != "Fungi_phy_incertae_sedis" &
        is.na(BLAST_Kingdom)                                                   ~ "UNITE_Fungi_no_BLAST",
      
      # UNITE Fungi confirmed by BLAST
      Kingdom == "Fungi" & Phylum != "Fungi_phy_incertae_sedis" &
        BLAST_Kingdom == "Fungi"                                               ~ "UNITE_and_BLAST_Fungi",
      
      # BLAST rescued from incertae sedis
      Kingdom == "Fungi" & Phylum == "Fungi_phy_incertae_sedis" &
        BLAST_Kingdom == "Fungi" & pident >= 80 & qcov >= 70                  ~ "BLAST_rescued_incertae_sedis",
      
      # BLAST rescued from unassigned
      (is.na(Kingdom) | Kingdom == "" | Kingdom == "Unassigned") &
        BLAST_Kingdom == "Fungi" & pident >= 80 & qcov >= 70                  ~ "BLAST_rescued_unassigned",
      
      # Non-fungal confirmed by BLAST
      !is.na(BLAST_Kingdom) & BLAST_Kingdom != "Fungi"                        ~
        paste0("Non_fungal_", BLAST_Kingdom),
      
      # Low-quality BLAST hit — could not confirm
      BLAST_Kingdom == "Fungi" & (pident < 80 | qcov < 70)                   ~ "Low_quality_BLAST_hit",
      
      # No BLAST hit at all
      is.na(BLAST_Kingdom)                                                     ~ "No_BLAST_hit",
      
      TRUE                                                                     ~ "Other"
    )
  )

cat("========== FILTERING SUMMARY ==========\n")
cat("\nAll OTUs by filter reason and keep decision:\n")
print(table(tax_combined$Filter_Reason, tax_combined$Keep))

cat("\nOTUs KEPT by reason:\n")
print(table(tax_combined$Filter_Reason[tax_combined$Keep == TRUE]))

cat("\nOTUs REMOVED by reason:\n")
print(table(tax_combined$Filter_Reason[tax_combined$Keep == FALSE]))

##------------------------------------------------------------------------------
## Part 4: Validate filtering ####
##------------------------------------------------------------------------------
# Hard stop if any confirmed non-fungal OTUs are being kept — indicates a
# logic error in the filtering above that must be resolved before continuing

non_fungal_kept <- tax_combined %>%
  filter(Keep == TRUE & grepl("Non_fungal", Filter_Reason))

if (nrow(non_fungal_kept) > 0) {
  cat("\nERROR: The following non-fungal OTUs are marked to keep — fix filtering logic:\n")
  print(non_fungal_kept %>% select(OTU_ID, Kingdom, BLAST_Kingdom, Filter_Reason))
  stop("Filtering validation failed. Resolve before continuing.")
} else {
  cat("\n✓ Validation passed: no non-fungal OTUs in kept set\n")
}

# Check for zero-read OTUs in the keep set
current_abundances   <- taxa_sums(ps_preblast)
otus_to_keep_ids     <- tax_combined %>% filter(Keep == TRUE) %>% pull(OTU_ID)
kept_abundances      <- current_abundances[otus_to_keep_ids]
zero_abundance_kept  <- sum(kept_abundances == 0, na.rm = TRUE)

cat("\nZero-read check for kept OTUs:\n")
cat("  OTUs marked to keep:", length(otus_to_keep_ids), "\n")
cat("  Of these, zero reads:", zero_abundance_kept, "\n")
cat("  Of these, >0 reads:", sum(kept_abundances > 0, na.rm = TRUE), "\n")
if (zero_abundance_kept > 0)
  cat("  Note: zero-read OTUs will be pruned in Part 6\n")

##----------------------------------------------------------------------------------------------
## Part 5: Update taxonomy table with BLAST-rescued assignments and filter low-read samples ####
##----------------------------------------------------------------------------------------------
# For OTUs rescued by BLAST, update Kingdom, Phylum, and Class with
# BLAST-derived assignments where the original UNITE taxonomy was missing.
# Standardize "incertae sedis" capitalization across all ranks for consistency.

tax_combined <- tax_combined %>%
  mutate(
    # Update Kingdom for BLAST-rescued OTUs
    Kingdom_Final = case_when(
      grepl("BLAST_rescued", Filter_Reason) ~ BLAST_Kingdom,
      TRUE                                  ~ Kingdom
    ),
    
    # Update Phylum for rescued OTUs where original was empty
    Phylum_Final = case_when(
      grepl("BLAST_rescued", Filter_Reason) & (is.na(Phylum) | Phylum == "") ~ BLAST_Phylum,
      TRUE                                                                    ~ Phylum
    ),
    
    # Update Class for rescued OTUs where original was empty
    Class_Final = case_when(
      grepl("BLAST_rescued", Filter_Reason) & (is.na(Class) | Class == "")   ~ BLAST_Class,
      TRUE                                                                    ~ Class
    )
  )

# Standardize "incertae sedis" capitalization — normalize to lowercase convention
# (UNITE inconsistently uses both "_Incertae_sedis" and "_incertae_sedis")
for (rank in c("Phylum_Final", "Class_Final", "Order", "Family")) {
  if (rank %in% colnames(tax_combined)) {
    tax_combined[[rank]] <- gsub("_Incertae_sedis", "_incertae_sedis",
                                 tax_combined[[rank]], ignore.case = FALSE)
  }
}

n_rescued <- sum(grepl("BLAST_rescued", tax_combined$Filter_Reason) & tax_combined$Keep)
cat("BLAST-rescued OTUs with updated taxonomy:", n_rescued, "\n\n")

##------------------------------------------------------------------------------
## Part 6: Build and export final fungi-only phyloseq object ####
##------------------------------------------------------------------------------

# Filter phyloseq to retained OTUs
ps_fungi_only <- prune_taxa(otus_to_keep_ids, ps_preblast)

# Build updated taxonomy matrix with final (potentially BLAST-corrected) assignments
tax_new <- tax_combined %>%
  filter(Keep == TRUE) %>%
  select(
    Kingdom = Kingdom_Final,
    Phylum  = Phylum_Final,
    Class   = Class_Final,
    Order, Family, Genus, Species
  ) %>%
  as.matrix()
rownames(tax_new) <- otus_to_keep_ids

# Apply updated taxonomy to phyloseq object
tax_table(ps_fungi_only) <- tax_table(tax_new)

cat("Samples with zero reads after BLAST filter:", 
    sum(sample_sums(ps_fungi_only) == 0), "\n")
cat("Sample names of zero-read samples:\n")
print(names(sample_sums(ps_fungi_only)[sample_sums(ps_fungi_only) == 0]))

# Prune samples and OTUs with zero reads after filtering
ps_fungi_only <- prune_samples(sample_sums(ps_fungi_only) > 0, ps_fungi_only)

zero_read_otus <- sum(taxa_sums(ps_fungi_only) == 0)
if (zero_read_otus > 0) {
  cat("Removing", zero_read_otus, "zero-read OTUs after BLAST filtering\n")
  ps_fungi_only <- prune_taxa(taxa_sums(ps_fungi_only) > 0, ps_fungi_only)
}

cat("========== POST-BLAST FILTERING SUMMARY ==========\n")
cat("OTUs before BLAST filtering:", ntaxa(ps_preblast), "\n")
cat("OTUs after BLAST filtering:", ntaxa(ps_fungi_only), "\n")
cat("OTUs removed:", ntaxa(ps_preblast) - ntaxa(ps_fungi_only), "\n")
cat("Samples:", nsamples(ps_fungi_only), "\n")
cat("Total reads:", format(sum(sample_sums(ps_fungi_only)), big.mark = ","), "\n")
cat("Reads retained vs. pre-BLAST:",
    round(100 * sum(sample_sums(ps_fungi_only)) / sum(sample_sums(ps_preblast)), 2), "%\n\n")

##----------------------------------------------------------------------------
## Save pre-read-depth-filter phyloseq object (for Figure 2 — all species) ##
##----------------------------------------------------------------------------
# Mirrors the equivalent save in the 16S decontam script. This object retains
# all biological samples that passed BLAST filtering (Parts 1–5) but have NOT
# yet been filtered by read depth. Used for Figure 2 where all host species
# are required regardless of sequencing depth.

ps_ITS_all_species <- ps_fungi_only  # capture state before depth filter

cat("--- Pre-read-depth-filter phyloseq object (all species, ITS) ---\n")
cat("  Samples:", nsamples(ps_ITS_all_species),
    "| OTUs:", ntaxa(ps_ITS_all_species), "\n")
cat("  Host species:",
    n_distinct(as(sample_data(ps_ITS_all_species), "data.frame")$host_taxon),
    "\n\n")

saveRDS(ps_ITS_all_species, file = file.path(dir_processed, "ITS_physeq_absolute_PRE_depth_filter.rds"))
cat("Saved:", file.path(dir_processed, "ITS_physeq_absolute_PRE_depth_filter.rds"), "\n\n")

# --- Post-BLAST read depth filter ---
# Some samples may drop below the read depth threshold after BLAST filtering
# removes non-fungal/unassigned OTUs. Apply the same 500-read floor here.
post_blast_depth_min <- 300

post_blast_depths     <- sample_sums(ps_fungi_only)
low_post_blast        <- names(post_blast_depths[post_blast_depths < post_blast_depth_min])

if (length(low_post_blast) > 0) {
  cat("Post-BLAST read depth filter (>=", post_blast_depth_min, "reads):\n")
  cat("  Samples removed:", length(low_post_blast), "\n")
  
  # Pull metadata from phyloseq — use as() to get a plain data.frame
  its_meta_df <- as(sample_data(ps_fungi_only), "data.frame")
  
  its_meta_lookup <- its_meta_df %>%
    tibble::rownames_to_column("sample_name") %>%
    filter(sample_name %in% low_post_blast) %>%
    select(sample_name, host_taxon, sample_type, env_broad_scale)
  
  removed_depth_df <- data.frame(
    sample_name = low_post_blast,
    read_depth  = post_blast_depths[low_post_blast]
  ) %>%
    left_join(its_meta_lookup, by = "sample_name") %>%
    arrange(read_depth)
  
  print(removed_depth_df)
  
  # Export removed samples to CSV before pruning
  write.csv(removed_depth_df,
            file      = file.path(dir_tables, "ITS_samples_removed_read_depth_filter.csv"),
            row.names = FALSE)
  cat("  Removed samples exported to:", file.path(dir_tables, "ITS_samples_removed_read_depth_filter.csv"), "\n")
  
  ps_fungi_only <- prune_samples(
    setdiff(sample_names(ps_fungi_only), low_post_blast),
    ps_fungi_only
  )
  
  # Remove any OTUs now at zero reads after sample removal
  zero_after_depth <- sum(taxa_sums(ps_fungi_only) == 0)
  if (zero_after_depth > 0) {
    cat("  OTUs removed (now zero-read):", zero_after_depth, "\n")
    ps_fungi_only <- prune_taxa(taxa_sums(ps_fungi_only) > 0, ps_fungi_only)
  }
  
  cat("  Samples retained:", nsamples(ps_fungi_only), "\n")
  cat("  OTUs retained:", ntaxa(ps_fungi_only), "\n")
  cat("  Total reads:", format(sum(sample_sums(ps_fungi_only)), big.mark = ","), "\n\n")
} else {
  cat("No samples below post-BLAST read depth threshold — no CSV written.\n\n")
}

# Save final phyloseq object
saveRDS(ps_fungi_only, file = file.path(dir_processed, "ITS_physeq_absolute_decontam_output.rds"))
cat("Final phyloseq saved:", file.path(dir_processed, "ITS_physeq_absolute_decontam_output.rds"), "\n\n")

# Export component tables as CSVs
# Note: use as(sample_data(ps_fungi_only), "data.frame") instead of 
# as.data.frame() to avoid validObject() trigger on sample_data slot
otu_final_df  <- as.data.frame(otu_table(ps_fungi_only))
tax_final_df  <- as.data.frame(tax_table(ps_fungi_only))
meta_final_df <- as(sample_data(ps_fungi_only), "data.frame")

write.csv(otu_final_df,  file = file.path(dir_tables, "ITS_count_absolute_dc.csv"))
write.csv(tax_final_df,  file = file.path(dir_tables, "ITS_tax_absolute_dc.csv"))
write.csv(meta_final_df, file = file.path(dir_tables, "ITS_metadata_absolute_dc.csv"), row.names = FALSE)

##------------------------------------------------------------------------------
## Part 7: Full validation of final phyloseq object ####
##------------------------------------------------------------------------------
cat("========== FINAL DATASET VALIDATION ==========\n\n")

tax_final_val  <- as.data.frame(tax_table(ps_fungi_only))
otu_abundances <- taxa_sums(ps_fungi_only)
sample_depths  <- sample_sums(ps_fungi_only)

# --- Sample depth summary ---
cat("--- Sample depth ---\n")
cat("  Samples:", nsamples(ps_fungi_only), "\n")
cat("  Min:", format(min(sample_depths), big.mark = ","), "\n")
cat("  25th pct:", format(quantile(sample_depths, 0.25), big.mark = ","), "\n")
cat("  Median:", format(median(sample_depths), big.mark = ","), "\n")
cat("  Mean:", format(round(mean(sample_depths)), big.mark = ","), "\n")
cat("  75th pct:", format(quantile(sample_depths, 0.75), big.mark = ","), "\n")
cat("  Max:", format(max(sample_depths), big.mark = ","), "\n")
low_read_samples <- sum(sample_depths < 300)
if (low_read_samples > 0)
  cat("  Note:", low_read_samples, "samples have < 500 reads — consider reviewing\n")

# Examine the low-read samples before deciding to remove
low_depth_df <- data.frame(
  sample     = names(sample_depths[sample_depths < 300]),
  read_depth = sample_depths[sample_depths < 300],
  host_taxon = as.data.frame(sample_data(ps_fungi_only))$host_taxon[
    match(names(sample_depths[sample_depths < 300]),
          sample_names(ps_fungi_only))]
) %>% arrange(read_depth)

cat("Read depth distribution of low-read samples:\n")
print(summary(low_depth_df$read_depth))
cat("\nFull sample list:\n")
print(low_depth_df)

# --- OTU abundance summary ---
cat("--- OTU abundance ---\n")
cat("  OTUs:", ntaxa(ps_fungi_only), "\n")
cat("  Total reads:", format(sum(otu_abundances), big.mark = ","), "\n")
cat("  Min per OTU:", format(min(otu_abundances), big.mark = ","), "\n")
cat("  Median per OTU:", format(median(otu_abundances), big.mark = ","), "\n")
cat("  Mean per OTU:", format(round(mean(otu_abundances)), big.mark = ","), "\n")
cat("  Max per OTU:", format(max(otu_abundances), big.mark = ","), "\n")
cat("  Host species:", n_distinct(as.data.frame(sample_data(ps_fungi_only))$host_taxon), "\n\n")

# --- Taxonomy completeness ---
cat("--- Taxonomy completeness ---\n")
for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  if (rank %in% colnames(tax_final_val)) {
    assigned <- sum(!is.na(tax_final_val[[rank]]) &
                      tax_final_val[[rank]] != "" &
                      tax_final_val[[rank]] != "Unassigned")
    pct <- round(100 * assigned / ntaxa(ps_fungi_only), 1)
    cat(sprintf("  %-10s: %4d / %4d assigned (%5.1f%%)\n",
                rank, assigned, ntaxa(ps_fungi_only), pct))
  }
}

# --- Kingdom verification (must be 100% Fungi) ---
cat("--- Kingdom verification ---\n")
kingdom_table <- table(tax_final_val$Kingdom, useNA = "always")
print(kingdom_table)
if (all(names(kingdom_table)[!is.na(names(kingdom_table))] == "Fungi")) {
  cat("✓ PASS: All assigned OTUs are Kingdom Fungi\n\n")
} else {
  cat("✗ FAIL: Non-fungal kingdoms detected — investigate before proceeding\n\n")
}

# --- Phylum distribution ---
cat("--- Top 10 phyla ---\n")
phylum_counts <- sort(table(tax_final_val$Phylum), decreasing = TRUE)
for (i in 1:min(10, length(phylum_counts))) {
  ph_name  <- names(phylum_counts)[i]
  ph_count <- phylum_counts[i]
  ph_otus  <- rownames(tax_final_val)[tax_final_val$Phylum == ph_name]
  ph_reads <- sum(otu_abundances[ph_otus])
  cat(sprintf("  %2d. %-35s %4d OTUs (%5.1f%%)  %s reads (%5.1f%%)\n",
              i, ph_name, ph_count,
              round(100 * ph_count / ntaxa(ps_fungi_only), 1),
              format(ph_reads, big.mark = ","),
              round(100 * ph_reads / sum(otu_abundances), 1)))
}

# --- Top 10 most abundant OTUs ---
cat("--- Top 10 most abundant OTUs ---\n")
top10_ids  <- names(sort(otu_abundances, decreasing = TRUE)[1:10])
top10_info <- tax_final_val[top10_ids, c("Kingdom", "Phylum", "Class", "Genus")]
top10_info$Total_Reads  <- otu_abundances[top10_ids]
top10_info$Pct_Total    <- round(100 * otu_abundances[top10_ids] / sum(otu_abundances), 2)
print(top10_info)

# --- Automated validation checks ---
cat("--- Automated checks ---\n")
checks_passed <- 0
checks_total  <- 5

# 1. All samples have reads
if (all(sample_depths > 0)) {
  cat("✓ All samples have reads\n"); checks_passed <- checks_passed + 1
} else { cat("✗ Some samples have zero reads\n") }

# 2. All OTUs have reads
if (all(otu_abundances > 0)) {
  cat("✓ All OTUs have reads\n"); checks_passed <- checks_passed + 1
} else { cat("✗ Some OTUs have zero reads\n") }

# 3. All assigned OTUs are Kingdom Fungi
all_fungi <- all(tax_final_val$Kingdom == "Fungi" | is.na(tax_final_val$Kingdom))
if (all_fungi) {
  cat("✓ All assigned OTUs are Kingdom Fungi\n"); checks_passed <- checks_passed + 1
} else { cat("✗ Non-fungal kingdoms present\n") }

# 4. Sample names match across components
meta_df <- as.data.frame(sample_data(ps_fungi_only))
if (all(sample_names(ps_fungi_only) == rownames(meta_df))) {
  cat("✓ Sample names match across components\n"); checks_passed <- checks_passed + 1
} else { cat("✗ Sample name mismatch between OTU table and metadata\n") }

# 5. OTU names match across components
if (all(taxa_names(ps_fungi_only) == rownames(tax_final_val))) {
  cat("✓ OTU names match across components\n"); checks_passed <- checks_passed + 1
} else { cat("✗ OTU name mismatch between OTU table and taxonomy\n") }

cat(sprintf("\n✓ Passed %d / %d validation checks\n", checks_passed, checks_total))
cat("\nAll outputs saved. ITS pipeline complete — proceed to downstream analyses.\n")

###############################################################################
###############################################################################
###############################################################################
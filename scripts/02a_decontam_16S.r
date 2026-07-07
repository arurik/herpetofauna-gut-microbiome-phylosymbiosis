###############################################################################
# Decontam Script - 16S rRNA
# Contaminant identification and removal workflow
# Includes: spike removal, pre-decontam abundance filtering, decontamination,
#           post-decontam abundance/prevalence filtering, read depth filtering,
#           and final taxonomic cleaning
# Alexander Rurik
###############################################################################

# ========================== LOAD PACKAGES =====================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(decontam)
  library(phyloseq)
  library(vegan)
  library(ggpubr)
  library(data.table)
  library(ggtext)
})

# Explicitly prioritize dplyr versions of commonly masked functions
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# ========================== BASIC PARAMETERS ==================================
set.seed(52325)                       # Reproducibility

# ========================== PROJECT DIRECTORIES ================================
library(here)
dir_processed <- here("data", "processed")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")

# ========================== LOAD DATA =========================================
# Load phyloseq object output from Script 01 (DspikeIn absolute abundance conversion)
ps <- readRDS(file.path(dir_processed, "physeq_absolute_16S.rds"))

# ========================== INITIAL DATASET CHECK =============================
cat("\n========== DATASET CHECK: START OF DECONTAM SCRIPT ==========\n")
cat("Phyloseq object loaded from DspikeIn output:\n")
print(ps)
cat("\nTotal reads (absolute abundance):", format(sum(otu_table(ps)), big.mark = ","), "\n")
cat("Number of samples:", nsamples(ps), "\n")
cat("Number of OTUs:", ntaxa(ps), "\n")

# ========================== EXTRACT COMPONENTS ================================
# Extract OTU table (samples as rows)
otu <- as.data.frame(phyloseq::otu_table(ps))
if (taxa_are_rows(ps)) otu <- t(otu)
otu <- as.data.frame(otu)

# Extract taxonomy table, move OTU IDs from rownames to a column
tax <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
tax <- tibble::rownames_to_column(tax, var = "OTU")

# Extract sample metadata, move sample names from rownames to a column
metadata <- data.frame(sample_data(ps))
metadata <- tibble::rownames_to_column(metadata, var = "sample_name")

# ========================== SPIKE-IN REMOVAL ==================================
# Identify Tetragenococcus OTUs (spike-in internal standard) in the taxonomy table
# and confirm they are present in the OTU table before removing
tetra.tax <- tax %>%
  filter(Genus == "Tetragenococcus")
tetra_intersect <- intersect(tetra.tax$OTU, colnames(otu))

cat("Spike-in (Tetragenococcus) OTUs found:", length(tetra_intersect), "\n")

# Remove spike-in OTUs from taxonomy table; rename OTU column and add "Otu_" prefix
tax_clean_filtered <- tax %>%
  filter(OTU %in% colnames(otu)) %>%
  filter(Genus != "Tetragenococcus") %>%
  rename(otu = OTU) %>%
  mutate(otu = paste0("Otu_", otu))

# Remove spike-in OTUs from OTU table and add "Otu_" prefix to all OTU column names
otu.spikeout <- otu %>%
  select(-all_of(tetra_intersect)) %>%
  rename_with(~ paste0("Otu_", .))

cat("OTUs remaining after spike removal:", ncol(otu.spikeout), "\n\n")

# ========================== PRE-DECONTAM: LOW-ABUNDANCE OTU FILTERING =========
# Context: Because absolute abundance scaling (DspikeIn) inflates raw read counts
# relative to conventional relative-abundance datasets, total reads per OTU are
# substantially higher. A threshold of 250 total reads is therefore appropriate
# to remove rare OTUs likely representing noise or sequencing artifacts, while
# preserving the vast majority of true biological signal.
# NOTE: This filter is applied BEFORE decontam intentionally — decontam's frequency
# and prevalence models perform better when not fed very low-abundance OTUs.
# A second (post-decontam) filter is applied later to catch marginal OTUs
# that may survive decontam.

# --- Retention curve: explore effect of different cutoffs on OTU and read retention ---
cutoffs <- c(1, 10, 30, 50, 100, 150, 250, 350, 450)
selected_threshold <- 250  # Final selected threshold (see retention plots below)

total_reads_pre <- sum(otu.spikeout)

results <- tibble(
  cutoff        = numeric(),
  OTUs_retained = numeric(),
  reads_retained = numeric()
)

for (c in cutoffs) {
  otu_filtered  <- otu.spikeout %>% select(where(~ is.numeric(.x) && sum(.x) > c))
  results <- results %>%
    add_row(
      cutoff         = c,
      OTUs_retained  = ncol(otu_filtered),
      reads_retained = sum(otu_filtered) / total_reads_pre
    )
}

# Plot OTUs retained vs. cutoff threshold
p1 <- ggplot(results, aes(x = cutoff, y = OTUs_retained)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = selected_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_x_continuous(breaks = cutoffs) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(
    x     = "OTU read abundance threshold",
    y     = "Number of OTUs retained",
    title = "A) Effect of abundance threshold on OTU retention - 16S"
  ) +
  theme_minimal()

# Plot % reads retained vs. cutoff threshold
p2 <- ggplot(results, aes(x = cutoff, y = reads_retained)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = selected_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_x_continuous(breaks = cutoffs) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    x     = "OTU read abundance threshold",
    y     = "% of total reads retained",
    title = "B) Effect of abundance threshold on read retention — 16S"
  ) +
  theme_minimal()

combined_plot <- ggarrange(p1, p2, ncol = 2, nrow = 1)
print(combined_plot)

# Save retention plots (uncomment to save)
# ggsave(file.path(dir_figures, "16S_prevalence_cutoff.pdf"), combined_plot, width = 11, height = 5, units = "in")
# ggsave(file.path(dir_figures, "16S_prevalence_cutoff.png"), combined_plot, width = 11, height = 5, units = "in", dpi = 600)

# --- Apply selected abundance threshold (250 total reads) ---
otu.mat.abund <- otu.spikeout %>%
  select(where(~ is.numeric(.x) && sum(.x) > 250))

# Remove any samples that now have zero reads after OTU filtering
zero_reads <- rowSums(otu.mat.abund) == 0
if (any(zero_reads)) {
  cat("Removing", sum(zero_reads), "samples with zero reads after abundance filtering\n")
  otu.mat.abund <- otu.mat.abund[!zero_reads, ]
  metadata      <- metadata %>% filter(sample_name %in% rownames(otu.mat.abund))
}

# Report read retention after pre-decontam filter
reads_before <- sum(otu.spikeout)
reads_after  <- sum(otu.mat.abund)
cat("\n--- Pre-decontam abundance filter summary ---\n")
cat("Threshold applied: >250 total reads per OTU\n")
cat("OTUs before filter:", ncol(otu.spikeout), "| After:", ncol(otu.mat.abund),
    "| Removed:", ncol(otu.spikeout) - ncol(otu.mat.abund), "\n")
cat("Reads before:", format(reads_before, big.mark = ","),
    "| After:", format(reads_after, big.mark = ","),
    "| Retained:", round(reads_after / reads_before * 100, 2), "%\n\n")

# ========================== DECONTAMINATION ====================================
# Uses decontam's "combined" method, which integrates:
#   - Frequency method: true taxa should not correlate with DNA concentration
#   - Prevalence method: contaminants are more prevalent in negative controls
# The workflow iterates over a range of probability thresholds and uses an
# analytical criterion to select the optimal threshold (see below).

# Control variable: set to "Y" to run decontam fresh, "N" to load a saved run
calculate <- "Y"

# File naming parameters for save/load logic
environ.desig <- "decontam.out"                                  # Name of object in R environment
decontam_cache_file <- file.path(dir_processed, "16S_decontam_out.rds")  # Path for cached output

if (calculate == "N") {
  # --- Load previously computed decontam results ---
  if (environ.desig %in% ls()) {
    print("Decontam output already exists in environment — skipping load")

  } else if (file.exists(decontam_cache_file)) {
    assign(environ.desig, readRDS(decontam_cache_file))
    print(paste0("Loaded cached decontam output from ", decontam_cache_file))

  } else {
    print("No saved decontam output found — set calculate = 'Y' to run")
  }
  
} else if (calculate == "Y") {
  
  print("Running decontamination — this may take several minutes...")
  
  # Reorder metadata rows to match the order of samples in the OTU matrix
  metadata <- otu.mat.abund %>%
    rownames_to_column("sample_name") %>%
    select(sample_name) %>%
    left_join(metadata, by = "sample_name")
  
  # Join OTU table with metadata and prepare for decontam
  # Samples are flagged as negative controls (TRUE) or true samples (FALSE)
  # based on sample name patterns (blanks, NTCs, pools, standards, etc.)
  otus.seqdat <- otu.mat.abund %>%
    rownames_to_column("sample_name") %>%
    right_join(metadata, ., by = "sample_name") %>%
    filter(!is.na(ampliconlibrary_quantification_ng.ul)) %>%  # Require known DNA concentration
    mutate(
      sample_type.1 = (
        grepl("blank",        sample_name, ignore.case = TRUE) |
          grepl("NTC",          sample_name, ignore.case = TRUE) |
          grepl("Pool",         sample_name, ignore.case = TRUE) |
          grepl("swab",         sample_name, ignore.case = TRUE) |
          grepl("Std",          sample_name, ignore.case = TRUE) |
          grepl("spiked",       sample_name, ignore.case = TRUE) |
          grepl("UHMspikedblank", sample_name, ignore.case = TRUE) |
          grepl("ziplock",      sample_name, ignore.case = TRUE) |
          grepl("ctrl",         sample_name, ignore.case = TRUE) |
          grepl("UHMNTC",       sample_name, ignore.case = TRUE)
      ),
      .after = sample_name
    )
  row.names(otus.seqdat) <- NULL  # Remove rownames artifact from OTU matrix join
  
  cat("Samples flagged as negative controls:", sum(otus.seqdat$sample_type.1), "\n")
  cat("Samples flagged as true samples:", sum(!otus.seqdat$sample_type.1), "\n\n")
  
  # Define probability thresholds to iterate over
  threshold <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
  
  # Initialize storage lists for loop outputs
  sapply(
    c("decontam.results.list", "decontam.summary", "plots.sampletype", "plots.persample"),
    function(k) assign(k, list(), envir = .GlobalEnv)
  )
  
  # Extract OTU matrix for decontam input (samples x OTUs, no metadata columns)
  decontam.seqs <- otus.seqdat %>%
    column_to_rownames("sample_name") %>%
    select(starts_with("Otu")) %>%
    as.matrix()
  
  # --- Main decontam loop: iterate over threshold values ---
  for (j in seq_along(threshold)) {
    
    print(paste("Running threshold:", threshold[j], paste0("(", j, "/", length(threshold), ")")))
    t.start <- Sys.time()
    
    # Run decontam: combined method uses both frequency (DNA conc.) and prevalence (neg controls)
    decontam.results <- isContaminant(
      seqtab    = decontam.seqs,
      neg       = otus.seqdat$sample_type.1,
      conc      = otus.seqdat$ampliconlibrary_quantification_ng.ul,
      method    = "combined",
      threshold = threshold[j],
      normalize = TRUE,
      detailed  = TRUE
    ) %>%
      filter(contaminant == TRUE) %>%
      rownames_to_column("otu") %>%
      left_join(., tax_clean_filtered, by = "otu")
    
    decontam.results.list[[j]] <- decontam.results
    names(decontam.results.list)[j] <- paste("Threshold:", threshold[j])
    
    contaminants <- decontam.results.list[[j]]$otu
    
    # Build long-format dataframe for visualization
    # Classifies each OTU per sample as "contamination" or "sample taxa"
    df.otu.long <- otus.seqdat %>%
      pivot_longer(starts_with("Otu"), names_to = "otu", values_to = "reads") %>%
      mutate(
        otu.type    = ifelse(otu %in% contaminants, "contamination", "sample taxa"),
        reads       = as.numeric(reads)
      ) %>%
      group_by(sample_type.1) %>%
      mutate(total.reads = sum(reads)) %>%
      ungroup() %>%
      group_by(sample_name) %>%
      mutate(sample.reads = sum(reads)) %>%
      ungroup() %>%
      mutate(otu.type = as.factor(str_to_title(otu.type)))
    
    # Ensure legend always shows both levels even if no contaminants found at this threshold
    if (!("Contamination" %in% levels(df.otu.long$otu.type))) {
      levels(df.otu.long$otu.type) <- c(levels(df.otu.long$otu.type), "Contamination")
    }
    
    # --- Plot: % contamination by sample type (blanks vs. true samples) ---
    # sample_n(1) is used here to reduce from OTU-level rows to one row per
    # sample x otu.type combination after computing per-group summary stats
    df.plot.sampletype <- df.otu.long %>%
      group_by(sample_name, otu.type) %>%
      mutate(
        otu.id.reads      = sum(reads),
        otu.id.proportion = otu.id.reads / sample.reads
      ) %>%
      sample_n(1) %>%
      group_by(sample_type.1, otu.type) %>%
      mutate(
        mean.otu.id.pro = mean(otu.id.proportion, na.rm = TRUE),
        se.otu.id.pro   = sd(otu.id.proportion, na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(sample_type.1 = ifelse(sample_type.1, "Extraction Blank", "Sample")) %>%
      ungroup()
    
    plot.sampletype <- df.plot.sampletype %>%
      ggplot(aes(x = sample_type.1, y = -mean.otu.id.pro, fill = otu.type)) +
      geom_errorbar(
        aes(ymax = -(mean.otu.id.pro + se.otu.id.pro),
            ymin = -(mean.otu.id.pro - se.otu.id.pro)),
        position = position_dodge(width = 0.9), width = 0.25
      ) +
      geom_col(color = "black", position = position_dodge(width = 0.9)) +
      scale_y_reverse(limits = c(0, -1), labels = scales::percent) +
      scale_fill_manual(
        drop   = FALSE,
        values = c("Contamination" = "#F8766D", "Sample Taxa" = "#00BFC4")
      ) +
      labs(
        title = paste("Threshold =", threshold[j]),
        x     = "Sample Type",
        y     = "Sample Composition"
      ) +
      theme_classic() +
      theme(plot.title = element_text(size = 12), legend.title = element_blank())
    
    decontam.summary[[j]] <- df.plot.sampletype %>%
      mutate(threshold = threshold[j]) %>%
      as.data.frame()
    
    plots.sampletype[[j]] <- plot.sampletype
    names(plots.sampletype)[j] <- paste("Threshold:", threshold[j])
    
    # --- Plot: reads per sample colored by contamination status ---
    # sample_n(1) used to collapse to one row per sample x otu.type after
    # computing read count sums — avoids duplicate rows in stacked bar chart
    #
    # NOTE: The y = read.counts - 100000 offset in geom_col below creates a
    # visual offset for the stacked bar outlines. This value is hardcoded and
    # may need adjustment if absolute abundance scales differ substantially
    # across dataset versions. Consider replacing with a dynamic offset:
    #   offset <- max(rowSums(otu.mat.abund)) * 0.01
    df.persample <- df.otu.long %>%
      group_by(sample_name) %>%
      mutate(sample.reads = sum(reads)) %>%
      ungroup() %>%
      mutate(
        sample_name  = fct_reorder(as.factor(sample_name), sample.reads, .desc = TRUE),
        sample.reads = as.numeric(sample.reads)
      ) %>%
      group_by(sample_name, otu.type) %>%
      mutate(read.counts = sum(reads)) %>%
      sample_n(1) %>%
      select(-otu)
    
    plot.persample <- df.persample %>%
      filter(read.counts > 0) %>%
      ggplot(aes(x = sample_name, y = read.counts, fill = otu.type)) +
      geom_col(aes(y = read.counts - 100000), color = "black", alpha = 1, linewidth = 0.5) +
      geom_col(width = 1, show.legend = FALSE) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, as.numeric(quantile(df.persample$read.counts, 0.999)))) +
      labs(title = paste("Threshold =", threshold[j])) +
      scale_fill_manual(values = c("Contamination" = "#F8766D", "Sample Taxa" = "#00BFC4")) +
      theme_classic() +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title   = element_text(size = 12)
      )
    
    plots.persample[[j]] <- plot.persample
    names(plots.persample)[j] <- paste("Threshold:", threshold[j])
    
    print(paste0("Loop ", j, "/", length(threshold),
                 "  |  Elapsed: ", round(abs(t.start - Sys.time()), 2), "s"))
  }
  
  # --- Summary plot: delta contamination removed at each threshold step ---
  # Used to analytically select the optimal threshold (see selection logic below)
  summary.df <- decontam.summary %>%
    bind_rows() %>%
    select(-c(otu, reads)) %>%
    filter(otu.type == "Contamination") %>%
    group_by(threshold, sample_type.1) %>%
    mutate(
      sample.type.reads    = total.reads,
      thresh.contam.reads  = sum(otu.id.reads)
    ) %>%
    sample_n(1) %>%
    ungroup() %>%
    mutate(
      thresh.contam.pro   = thresh.contam.reads / sample.type.reads,
      delta.thresh.contam = -(thresh.contam.pro - data.table::shift(thresh.contam.pro, n = 2))
    ) %>%
    mutate(delta.thresh.contam = ifelse(is.na(delta.thresh.contam), thresh.contam.pro, delta.thresh.contam))
  
  decontam.summary.plot <- summary.df %>%
    ggplot(aes(x = as.factor(threshold), y = delta.thresh.contam, fill = sample_type.1)) +
    geom_col(position = position_dodge(), color = "black") +
    scale_y_reverse(labels = scales::percent) +
    labs(y = "\u0394 Total Reads", x = "Threshold Value") +
    scale_fill_brewer(palette = "Paired") +
    theme_classic() +
    theme(
      legend.title       = element_blank(),
      legend.position    = "bottom",
      text               = element_text(size = 12),
      legend.text        = element_text(size = 14),
      legend.spacing.x   = unit(0.5, "cm"),
      axis.title.x       = element_text(margin = margin(10, 0, 0, 0)),
      axis.title.y       = ggtext::element_markdown(margin = margin(0, 15, 0, 10))
    )
  
  # Bundle all decontam output into a single named list for easy saving/loading
  decontam.out <- list(
    decontam.results.list = decontam.results.list,
    plots.persample       = plots.persample,
    plots.sampletype      = plots.sampletype,
    decontam.summary      = decontam.summary,
    summary.df            = summary.df,
    decontam.summary.plot = decontam.summary.plot
  )
  
  # Save decontam output — avoids re-running the computationally expensive loop
  saveRDS(decontam.out, file = decontam_cache_file)
}

# ========================== LOAD SAVED DECONTAM OUTPUT ========================
# Load from a specific saved run (update filename as needed)
#decontam.out <- readRDS(decontam_cache_file)

# ========================== THRESHOLD SELECTION & REVIEW ======================
# Visualize decontam results across thresholds before selecting final value
(decontam.summary.plot <- decontam.out$decontam.summary.plot)
decontam.out$plots.persample
decontam.out$plots.sampletype

# Analytical threshold selection criterion:
# Select the lowest threshold at which proportionally more reads are removed
# from negative controls/blanks than from true samples, after the cumulative
# reads removed from controls exceeds 10% (treating the first 10% as background noise)
threshold.decontam <- decontam.out$summary.df %>%
  group_by(sample_type.1) %>%
  mutate(reads.remove.cumulative = cumsum(thresh.contam.pro)) %>%
  select(threshold, sample_type.1, thresh.contam.pro, reads.remove.cumulative) %>%
  mutate(ntc.reads.remove.cumulative = ifelse(sample_type.1 == "Extraction Blank",
                                              reads.remove.cumulative, NA)) %>%
  group_by(threshold) %>%
  fill(ntc.reads.remove.cumulative) %>%
  ungroup() %>%
  filter(ntc.reads.remove.cumulative > 0.1) %>%
  select(threshold, sample_type.1, thresh.contam.pro) %>%
  mutate(sample_type.1 = recode(sample_type.1,
                                "Extraction Blank" = "ntc",
                                "Sample"           = "sample")) %>%
  pivot_wider(names_from = "sample_type.1", values_from = thresh.contam.pro) %>%
  filter(ntc > sample) %>%
  arrange(threshold) %>%
  slice_head(n = 1) %>%
  pull(threshold)

threshold.decontam.index <- paste("Threshold:", threshold.decontam)
cat("Analytically selected decontam threshold:", threshold.decontam, "\n")

# ========================== APPLY DECONTAM THRESHOLD ==========================
# Extract contaminant OTU list at selected threshold and remove from OTU table
decontam.final           <- decontam.out[["decontam.results.list"]][[threshold.decontam.index]]
contaminants             <- decontam.final$otu
otu.mat.abund.decontam   <- select(otu.mat.abund, -any_of(contaminants))

# Report decontam summary statistics
reads_before_decontam  <- sum(otu.mat.abund)
reads_after_decontam   <- sum(otu.mat.abund.decontam)
otus_before_decontam   <- ncol(otu.mat.abund)
otus_after_decontam    <- ncol(otu.mat.abund.decontam)

cat("\n========== DECONTAMINATION SUMMARY ==========\n")
cat("Threshold applied:", threshold.decontam, "\n")
cat("Reads before:", format(reads_before_decontam, big.mark = ","),
    "| After:", format(reads_after_decontam, big.mark = ","),
    "| Removed:", format(reads_before_decontam - reads_after_decontam, big.mark = ","),
    "| Retained:", round(reads_after_decontam / reads_before_decontam * 100, 2), "%\n")
cat("OTUs before:", otus_before_decontam,
    "| After:", otus_after_decontam,
    "| Removed:", otus_before_decontam - otus_after_decontam,
    "| Retained:", round(otus_after_decontam / otus_before_decontam * 100, 2), "%\n")

# ========================== GENERATE DECONTAM FIGURES =========================
# Finalize plots at the selected threshold for supplementary materials

decontam.out$plots.persample[[threshold.decontam.index]] +
  labs(
    x = "Samples<br><span style='font-size:11pt'>(arranged by sequencing depth)</span>",
    y = "16S Reads"
  ) +
  theme(
    legend.title     = element_blank(),
    plot.title       = element_blank(),
    legend.position  = c(0.6, 0.9),
    text             = element_text(size = 12),
    legend.text      = element_text(size = 14),
    legend.spacing.x = unit(0.5, "cm"),
    axis.title.x     = element_markdown(margin = margin(10, 0, 0, 0)),
    axis.title.y     = element_text(margin = margin(0, 15, 0, 10))
  ) -> persample.plot

decontam.out$plots.sampletype[[threshold.decontam.index]] +
  labs(y = "% Reads<br><span style='font-size:12pt'>(Mean &plusmn; SE)</span>") +
  theme(
    legend.title     = element_blank(),
    plot.title       = element_blank(),
    legend.position  = "none",
    text             = element_text(size = 12),
    legend.text      = element_text(size = 14),
    legend.spacing.x = unit(0.5, "cm"),
    axis.title.x     = element_blank(),
    axis.title.y     = element_markdown(margin = margin(0, 15, 0, 10))
  ) -> sampletype.plot

decontam.summary.plot <- decontam.out$decontam.summary.plot +
  annotate(
    "text",
    x     = 6,
    y     = max(decontam.out$summary.df$delta.thresh.contam) * -12,
    label = paste("Analytically selected threshold:", threshold.decontam),
    size  = 4,
    color = "black"
  )

# Assemble multi-panel supplementary figure
multi_panel <- ggpubr::ggarrange(
  decontam.summary.plot,
  ggpubr::ggarrange(
    sampletype.plot, persample.plot,
    ncol        = 2,
    align       = "hv",
    labels      = c("B)", "C)"),
    font.label  = list(size = 12, color = "black", face = "plain", family = "sans")
  ),
  labels     = c("A)", ""),
  font.label = list(size = 12, color = "black", face = "plain", family = "sans"),
  nrow       = 2
)
multi_panel <- ggpubr::annotate_figure(
  multi_panel,
  top = text_grob(
    "Analytical selection of decontamination threshold — 16S",
    size = 16,
    face = "plain"  
  )
)

print(multi_panel)

# Save decontam figures (uncomment to save)
# ggsave(file.path(dir_figures, "16S_decontam_threshold.pdf"), multi_panel, width = 10, height = 7, units = "in", dpi = 600)
# ggsave(file.path(dir_figures, "16S_decontam_threshold.png"), multi_panel, width = 10, height = 7, units = "in", dpi = 600)

##------------------------------------------------------------------------------
## Post-decontam filtering & dataset finalization ####
##------------------------------------------------------------------------------
# Filtering is applied in the following order to ensure read depth statistics
# reflect only true biological samples with clean OTU tables:
#
#   1. Remove controls/blanks from metadata & OTU table
#   2. Remove non-target OTUs (Archaea, Mitochondria, Chloroplast, non-bacterial)
#   3. Second abundance filter (remove low-abundance OTUs post-decontam)
#   4. Read depth filter (remove samples with unreliably low sequencing depth)
#
# Counts are reported before and after each step.
##------------------------------------------------------------------------------

# Snapshot dimensions entering this section (immediately after decontam removal)
cat("\n========== POST-DECONTAM FILTERING ==========\n")
cat("Entering post-decontam section:\n")
cat("  Samples:", nrow(otu.mat.abund.decontam),
    "| OTUs:", ncol(otu.mat.abund.decontam),
    "| Total reads:", format(sum(otu.mat.abund.decontam), big.mark = ","), "\n\n")

##------------------------------------------------------------------------------
## Step 1: Remove controls, blanks, and samples in the spike range ####
##------------------------------------------------------------------------------
# Keep only true biological samples by filtering on metadata fields.
# This must happen first so that subsequent read depth stats reflect
# only animal samples, not extraction blanks or NTCs.

metadata.final <- metadata %>%
  filter(sample_or_blank == "sample") %>%       # Remove blanks/controls/NTCs
  filter(Spiked_Reads > 0) %>%                  # Must have spike-in reads
  filter(Percentage >= 0.01, Percentage <= 40)  # Spike-in % within DspikeIn QC range

# Subset OTU table to match filtered metadata
otu.working <- otu.mat.abund.decontam %>%
  rownames_to_column("sample_name") %>%
  filter(sample_name %in% metadata.final$sample_name) %>%
  column_to_rownames("sample_name")

# Sync metadata to samples present in OTU table
metadata.final <- metadata.final %>%
  filter(sample_name %in% rownames(otu.working))

cat("--- Step 1: Remove controls/blanks/samples outside the spike range ---\n")
cat("  Samples retained:", nrow(otu.working),
    "| Removed:", nrow(otu.mat.abund.decontam) - nrow(otu.working), "\n")
cat("  OTUs:", ncol(otu.working),
    "| Total reads:", format(sum(otu.working), big.mark = ","), "\n")
cat("  Host species represented:", n_distinct(metadata.final$host_taxon), "\n\n")

##------------------------------------------------------------------------------
## Step 2: Remove non-target OTUs ####
##------------------------------------------------------------------------------
# Remove Archaea, Mitochondria, and Chloroplasts.
# These classifications can appear at different taxonomic ranks depending on
# the reference database, so we check multiple ranks to ensure complete removal.

# Subset taxonomy to OTUs present in current working table
tax.working <- tax_clean_filtered %>%
  filter(otu %in% colnames(otu.working))

otus_to_remove <- tax.working %>%
  filter(
    Kingdom == "Archaea"          |
      Family  == "Mitochondria"     |
      Order   == "Mitochondria"     |
      Order   == "Chloroplast"      |
      Class   == "Chloroplast"
  ) %>%
  pull(otu)

if (length(otus_to_remove) > 0) {
  removed_tax <- tax.working %>%
    filter(otu %in% otus_to_remove) %>%
    mutate(category = case_when(
      Kingdom == "Archaea"                                  ~ "Archaea",
      Family  == "Mitochondria" | Order == "Mitochondria"  ~ "Mitochondria",
      Order   == "Chloroplast"  | Class == "Chloroplast"   ~ "Chloroplast",
      TRUE                                                  ~ "Other"
    ))
  cat("--- Step 2: Remove non-target OTUs ---\n")
  cat("  OTUs removed by category:\n")
  print(table(removed_tax$category))
  otu.working  <- otu.working  %>% select(-all_of(otus_to_remove))
  tax.working  <- tax.working  %>% filter(!otu %in% otus_to_remove)
} else {
  cat("--- Step 2: Remove non-target OTUs ---\n")
  cat("  No Archaea, Mitochondria, or Chloroplast OTUs detected\n")
}

# Remove non-bacterial kingdoms (Eukaryota, unassigned, etc.)
non_bacterial <- tax.working %>%
  filter(!Kingdom %in% "Bacteria" | is.na(Kingdom))

if (nrow(non_bacterial) > 0) {
  cat("  Removing", nrow(non_bacterial), "non-bacterial/unassigned OTUs\n")
  cat("  Kingdoms:\n"); print(table(non_bacterial$Kingdom))
  otu.working  <- otu.working  %>% select(-all_of(non_bacterial$otu))
  tax.working  <- tax.working  %>% filter(!otu %in% non_bacterial$otu)
} else {
  cat("  No non-bacterial kingdoms detected\n")
}

# Remove any OTUs now at zero reads across all samples
zero_read_otus <- colSums(otu.working) == 0
if (any(zero_read_otus)) {
  cat("  Removing", sum(zero_read_otus), "zero-read OTUs after taxonomic cleaning\n")
  otu.working  <- otu.working[, !zero_read_otus]
  tax.working  <- tax.working %>% filter(!otu %in% names(which(zero_read_otus)))
}

cat("  After Step 2:\n")
cat("  Samples:", nrow(otu.working),
    "| OTUs:", ncol(otu.working),
    "| Total reads:", format(sum(otu.working), big.mark = ","), "\n")
cat("  Host species represented:", n_distinct(metadata.final$host_taxon), "\n\n")

##------------------------------------------------------------------------------
## Step 3: Second abundance filter ####
##------------------------------------------------------------------------------
# Re-apply abundance filter after decontam. Rationale: removing contaminant OTUs
# can reduce total reads for some OTUs below the original threshold, and new
# zero-abundance OTUs may have emerged after taxonomic cleaning above.
# Threshold matches the pre-decontam filter (250 reads) for consistency.

post_abund_threshold <- 250  # Min. total reads per OTU across all retained samples

otus_before_abund2 <- ncol(otu.working)

otu.working <- otu.working %>%
  select(where(~ is.numeric(.x) && sum(.x) > post_abund_threshold))

cat("--- Step 3: Second abundance filter (>", post_abund_threshold, "total reads) ---\n")
cat("  OTUs retained:", ncol(otu.working),
    "| Removed:", otus_before_abund2 - ncol(otu.working), "\n")
cat("  Samples:", nrow(otu.working),
    "| Total reads:", format(sum(otu.working), big.mark = ","), "\n")
cat("  Host species represented:", n_distinct(metadata.final$host_taxon), "\n\n")

# Sync taxonomy to OTUs remaining after abundance filter
tax.working <- tax.working %>%
  filter(otu %in% colnames(otu.working))

##------------------------------------------------------------------------------
## Save pre-read-depth-filter phyloseq object (for Figure 2 — all species) ####
##------------------------------------------------------------------------------
# This object retains all biological samples that passed Steps 1–3
# (blank removal, taxonomic cleaning, abundance filtering) but have NOT yet
# been filtered by read depth. This is intentional: Figure 2 requires all
# host species regardless of sequencing depth, and the read depth filter
# exists primarily to improve reliability of beta-diversity estimates rather
# than to exclude species from occurrence-based analyses.

physeq_all_species <- phyloseq(
  otu_table(otu.working, taxa_are_rows = FALSE),
  tax_table(as.matrix(tax.working %>% column_to_rownames("otu"))),
  sample_data(metadata.final %>% column_to_rownames("sample_name"))
)

cat("--- Pre-read-depth-filter phyloseq object (all species) ---\n")
cat("  Samples:", nsamples(physeq_all_species),
    "| OTUs:", ntaxa(physeq_all_species), "\n")
cat("  Host species:", n_distinct(metadata.final$host_taxon), "\n\n")

saveRDS(physeq_all_species,
        file = file.path(dir_processed, "16S_physeq_absolute_PRE_depth_filter.rds"))
cat("Saved:", file.path(dir_processed, "16S_physeq_absolute_PRE_depth_filter.rds"), "\n")

##------------------------------------------------------------------------------
## Step 4: Read depth filter ####
##------------------------------------------------------------------------------
# Remove samples with very low total reads. This filter is applied last —
# after blanks and non-target reads are removed — so the depth calculation
# reflects clean biological signal only.
#
# Rationale for 500-read floor:
# In absolute abundance data, low-read samples have DspikeIn scaling factors
# estimated from very few spike-in reads, making absolute abundance estimates
# unreliable regardless of apparent read count. Review the printed sample list
# and spike-in recovery values before finalizing this threshold.

post_read_depth_min <- 500  # Adjust after reviewing distribution below

sample_read_depths <- rowSums(otu.working)

cat("--- Step 4: Read depth filter (>=", post_read_depth_min, "reads per sample) ---\n")
cat("  Read depth distribution (biological samples, clean OTUs):\n")
cat("    Min:", format(min(sample_read_depths), big.mark = ","), "\n")
cat("    25th pct:", format(quantile(sample_read_depths, 0.25), big.mark = ","), "\n")
cat("    Median:", format(median(sample_read_depths), big.mark = ","), "\n")
cat("    Mean:", format(round(mean(sample_read_depths)), big.mark = ","), "\n")
cat("    75th pct:", format(quantile(sample_read_depths, 0.75), big.mark = ","), "\n")
cat("    Max:", format(max(sample_read_depths), big.mark = ","), "\n")
cat("  Samples below threshold:", sum(sample_read_depths < post_read_depth_min), "\n\n")

# Print specific samples falling below threshold for manual review before removal
low_depth_samples <- names(sample_read_depths[sample_read_depths < post_read_depth_min])
if (length(low_depth_samples) > 0) {
  cat("  Samples flagged for removal (read depth <", post_read_depth_min, "):\n")
  meta_lookup <- metadata.final %>%
    filter(sample_name %in% low_depth_samples) %>%
    select(sample_name, host_taxon, sample_type, env_broad_scale)
  
  low_depth_df <- data.frame(
    sample_name = low_depth_samples,
    read_depth  = sample_read_depths[low_depth_samples]
  ) %>%
    left_join(meta_lookup, by = "sample_name") %>%
    arrange(read_depth)
  print(low_depth_df)
  cat("\n")
}

# Apply depth filter
samples_pass_depth <- names(sample_read_depths[sample_read_depths >= post_read_depth_min])
otu.working        <- otu.working[samples_pass_depth, ]
metadata.final     <- metadata.final %>% filter(sample_name %in% samples_pass_depth)

# Export removed samples to CSV for records
if (length(low_depth_samples) > 0) {
  write.csv(low_depth_df,
            file      = file.path(dir_tables, "16S_samples_removed_read_depth_filter.csv"),
            row.names = FALSE)
  cat("  Removed samples exported to:", file.path(dir_tables, "16S_samples_removed_read_depth_filter.csv"), "\n")
} else {
  cat("  No samples removed by read depth filter — no CSV written.\n")
}

# Remove any OTUs now at zero reads after sample removal
zero_after_depth <- colSums(otu.working) == 0
if (any(zero_after_depth)) {
  cat("  Removing", sum(zero_after_depth),
      "OTUs with zero reads after low-depth sample removal\n")
  otu.working <- otu.working[, !zero_after_depth]
  tax.working <- tax.working %>% filter(!otu %in% names(which(zero_after_depth)))
}

cat("  After Step 4:\n")
cat("  Samples:", nrow(otu.working),
    "| OTUs:", ncol(otu.working),
    "| Total reads:", format(sum(otu.working), big.mark = ","), "\n")
cat("  Host species represented:", n_distinct(metadata.final$host_taxon), "\n\n")

##------------------------------------------------------------------------------
## Assign final objects & print summary ####
##------------------------------------------------------------------------------
otu.final  <- otu.working
tax.final  <- tax.working

cat("========== FINAL DATASET SUMMARY ==========\n")
cat("Samples:", nrow(otu.final), "\n")
cat("OTUs:", ncol(otu.final), "\n")
cat("Total reads:", format(sum(otu.final), big.mark = ","), "\n")
cat("Host species (host_taxon):", n_distinct(metadata.final$host_taxon), "\n")
cat("\nReads per sample:\n")
cat("  Min:", format(min(rowSums(otu.final)), big.mark = ","), "\n")
cat("  25th pct:", format(quantile(rowSums(otu.final), 0.25), big.mark = ","), "\n")
cat("  Median:", format(median(rowSums(otu.final)), big.mark = ","), "\n")
cat("  Mean:", format(round(mean(rowSums(otu.final))), big.mark = ","), "\n")
cat("  75th pct:", format(quantile(rowSums(otu.final), 0.75), big.mark = ","), "\n")
cat("  Max:", format(max(rowSums(otu.final)), big.mark = ","), "\n")
cat("============================================\n\n")

##------------------------------------------------------------------------------
## Build & export final phyloseq object ####
##------------------------------------------------------------------------------
physeq_final <- phyloseq(
  otu_table(otu.final, taxa_are_rows = FALSE),
  tax_table(as.matrix(tax.final %>% column_to_rownames("otu"))),
  sample_data(metadata.final %>% column_to_rownames("sample_name"))
)

cat("Final phyloseq object:", nsamples(physeq_final), "samples,",
    ntaxa(physeq_final), "OTUs\n")

# Export phyloseq object
saveRDS(physeq_final, file = file.path(dir_processed, "16S_physeq_absolute_decontam_output.rds"))

# Export component tables as CSVs
write.csv(otu.final,      file = file.path(dir_tables, "16S_count_absolute_dc.csv"))
write.csv(metadata.final, file = file.path(dir_tables, "16S_metadata_absolute_dc.csv"), row.names = FALSE)
write.csv(tax.final,      file = file.path(dir_tables, "16S_tax_absolute_dc.csv"),      row.names = FALSE)

cat("All outputs saved. Proceed to downstream analyses.\n")


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
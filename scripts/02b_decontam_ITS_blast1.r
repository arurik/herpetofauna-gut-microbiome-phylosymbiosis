###############################################################################
# Decontam Script - ITS1 rDNA
# Contaminant identification and removal workflow
# Includes: spike removal, pre-decontam abundance filtering, decontamination,
#           post-decontam taxonomic cleaning, abundance filtering, read depth
#           filtering, and extracts unknown/uncertain OTU sequences for BLAST
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
dir_raw       <- here("data", "raw", "ITS")
dir_processed <- here("data", "processed")
dir_figures   <- here("output", "figures")
dir_tables    <- here("output", "tables")

# ========================== LOAD DATA =========================================
# Load phyloseq object output from Script 01 (DspikeIn absolute abundance conversion)
ps <- readRDS(file.path(dir_processed, "physeq_absolute_ITS.rds"))

###############################################################################
# *** DATABASE COMPARISON TEST ***
# Uncomment this block to load the UNITE all-eukaryotes test dataset instead
# of the standard UNITE fungi-only dataset. Used to evaluate the effect of
# reference database choice on taxonomy assignment and OTU composition.
# Comment out the standard load above before running this block.
#
# ps <- readRDS(file.path(dir_raw, "physeq_absolute_ITS_allEUK.rds"))
###############################################################################

# ========================== INITIAL DATASET CHECK =============================
cat("\n========== DATASET CHECK: START OF DECONTAM SCRIPT ==========\n")
cat("Phyloseq object loaded from DspikeIn output:\n")
print(ps)
cat("\nTotal reads (absolute abundance):", format(sum(otu_table(ps)), big.mark = ","), "\n")
cat("Number of samples:", nsamples(ps), "\n")
cat("Number of OTUs:", ntaxa(ps), "\n")
cat("==============================================================\n\n")

# ========================== EXTRACT COMPONENTS ================================
# Extract OTU table (samples as rows)
otu <- as.data.frame(phyloseq::otu_table(ps))
if (taxa_are_rows(ps)) otu <- t(otu)
otu <- as.data.frame(otu)

# Extract taxonomy table; move OTU IDs from rownames to a column
tax <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
tax <- tibble::rownames_to_column(tax, var = "OTU")

# Extract sample metadata; move sample names from rownames to a column
metadata <- data.frame(sample_data(ps))
metadata <- tibble::rownames_to_column(metadata, var = "sample_name")

# ========================== SPIKE-IN REMOVAL ==================================
# ITS spike-in is Dekkera (a yeast genus used as an internal standard in DspikeIn)
# Identify Dekkera OTUs in the taxonomy table and confirm presence in OTU table
dekkera.tax       <- tax %>% filter(Genus == "Dekkera")
dekkera_intersect <- intersect(dekkera.tax$OTU, colnames(otu))

cat("Spike-in (Dekkera) OTUs found:", length(dekkera_intersect), "\n")

# Remove spike-in OTUs from taxonomy table; rename OTU column and add "Otu_" prefix
tax_clean_filtered <- tax %>%
  filter(OTU %in% colnames(otu)) %>%
  filter(Genus != "Dekkera") %>%
  rename(otu = OTU) %>%
  mutate(otu = paste0("Otu_", otu))

# Remove spike-in OTUs from OTU table and add "Otu_" prefix to all OTU column names
otu.spikeout <- otu %>%
  select(-all_of(dekkera_intersect)) %>%
  rename_with(~ paste0("Otu_", .))

cat("OTUs remaining after spike removal:", ncol(otu.spikeout), "\n\n")

# ========================== PRE-DECONTAM: LOW-ABUNDANCE OTU FILTERING =========
# Context: Because absolute abundance scaling (DspikeIn) inflates raw read counts
# relative to conventional relative-abundance datasets, total reads per OTU are
# substantially higher. A threshold of 150 total reads is appropriate for the
# ITS dataset (vs. 250 for 16S) given its smaller overall size and read depth.
# NOTE: Applied BEFORE decontam intentionally — decontam's frequency and
# prevalence models perform better without very low-abundance noise OTUs.
# A second (post-decontam) filter is applied later to catch marginal OTUs
# that may survive decontam.

# --- Retention curve: explore effect of different cutoffs on OTU and read retention ---
cutoffs            <- c(1, 10, 30, 50, 100, 150, 250, 350, 450)
selected_threshold <- 150  # Final selected threshold (see retention plots below)

total_reads_pre <- sum(otu.spikeout)

results <- tibble(
  cutoff         = numeric(),
  OTUs_retained  = numeric(),
  reads_retained = numeric()
)

for (c in cutoffs) {
  otu_filtered <- otu.spikeout %>% select(where(~ is.numeric(.x) && sum(.x) > c))
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
    title = "C) Effect of abundance threshold on OTU retention — ITS"
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
    title = "D) Effect of abundance threshold on read retention — ITS"
  ) +
  theme_minimal()

combined_plot <- ggarrange(p1, p2, ncol = 2, nrow = 1)
print(combined_plot)

# Save retention plots (uncomment to save)
# ggsave(file.path(dir_figures, "ITS_prevalence_cutoff.pdf"), combined_plot, width = 11, height = 5, units = "in")

# --- Apply selected abundance threshold (150 total reads) ---
otu.mat.abund <- otu.spikeout %>%
  select(where(~ is.numeric(.x) && sum(.x) > 150))

# Remove any samples that now have zero reads after OTU filtering
zero_reads <- rowSums(otu.mat.abund) == 0
if (any(zero_reads)) {
  cat("Removing", sum(zero_reads), "samples with zero reads after abundance filtering\n")
  cat("  Samples removed:", paste(rownames(otu.mat.abund)[zero_reads], collapse = ", "), "\n")
  otu.mat.abund <- otu.mat.abund[!zero_reads, ]
  metadata      <- metadata %>% filter(sample_name %in% rownames(otu.mat.abund))
} else {
  cat("No samples with zero reads after abundance filtering\n")
}

# Report read retention after pre-decontam filter
reads_before <- sum(otu.spikeout)
reads_after  <- sum(otu.mat.abund)
cat("\n--- Pre-decontam abundance filter summary ---\n")
cat("Threshold applied: >150 total reads per OTU\n")
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
decontam_cache_file <- file.path(dir_processed, "ITS_decontam_out.rds")  # Path for cached output

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
        grepl("blank",           sample_name, ignore.case = TRUE) |
          grepl("NTC",             sample_name, ignore.case = TRUE) |
          grepl("Pool",            sample_name, ignore.case = TRUE) |
          grepl("swab",            sample_name, ignore.case = TRUE) |
          grepl("Std",             sample_name, ignore.case = TRUE) |
          grepl("spiked",          sample_name, ignore.case = TRUE) |
          grepl("UHMspikedblank",  sample_name, ignore.case = TRUE) |
          grepl("ziplock",         sample_name, ignore.case = TRUE) |
          grepl("ctrl",            sample_name, ignore.case = TRUE) |
          grepl("UHMNTC",          sample_name, ignore.case = TRUE)
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
        otu.type = ifelse(otu %in% contaminants, "contamination", "sample taxa"),
        reads    = as.numeric(reads)
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
    # sample_n(1) used to reduce from OTU-level rows to one row per
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
    # visual separation for bar outlines. This hardcoded value may need
    # adjustment if ITS absolute abundance scales differ across dataset versions.
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
      sample.type.reads   = total.reads,
      thresh.contam.reads = sum(otu.id.reads)
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
      legend.title     = element_blank(),
      legend.position  = "bottom",
      text             = element_text(size = 12),
      legend.text      = element_text(size = 14),
      legend.spacing.x = unit(0.5, "cm"),
      axis.title.x     = element_text(margin = margin(10, 0, 0, 0)),
      axis.title.y     = ggtext::element_markdown(margin = margin(0, 15, 0, 10))
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
decontam.final         <- decontam.out[["decontam.results.list"]][[threshold.decontam.index]]
contaminants           <- decontam.final$otu
otu.mat.abund.decontam <- select(otu.mat.abund, -any_of(contaminants))

# Report decontam summary statistics
reads_before_decontam <- sum(otu.mat.abund)
reads_after_decontam  <- sum(otu.mat.abund.decontam)
otus_before_decontam  <- ncol(otu.mat.abund)
otus_after_decontam   <- ncol(otu.mat.abund.decontam)

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
    y = "ITS Reads"
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

# Explicitly set axis limits before passing to ggarrange —
# scale_y_reverse() can interact badly with ggarrange/annotate_figure
# causing axis expansion. Locking limits and breaks here prevents that.
decontam.summary.plot <- decontam.out$decontam.summary.plot +
  annotate(
    "text",
    x     = 6,
    y     = max(decontam.out$summary.df$delta.thresh.contam) * 0.5,
    label = paste("Analytically selected threshold:", threshold.decontam),
    size  = 4,
    color = "black"
  ) +
  scale_y_reverse(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0.50, -0.20),   # positive = bottom of reversed axis, negative = top
    breaks = seq(0.50, -0.20, by = -0.10)
  )

# Assemble multi-panel supplementary figure
multi_panel <- ggpubr::ggarrange(
  decontam.summary.plot,
  ggpubr::ggarrange(
    sampletype.plot, persample.plot,
    ncol       = 2,
    align      = "hv",
    labels     = c("B)", "C)"),
    font.label = list(size = 12, color = "black", face = "plain", family = "sans")
  ),
  labels     = c("A)", ""),
  font.label = list(size = 12, color = "black", face = "plain", family = "sans"),
  nrow       = 2
)
multi_panel <- ggpubr::annotate_figure(
  multi_panel,
  top = text_grob(
    "Analytical selection of decontamination threshold — ITS",
    size = 16,
    face = "plain"
  )
)
print(multi_panel)

# Save decontam figures (uncomment to save)
# ggsave(file.path(dir_figures, "ITS_decontam_threshold.pdf"), multi_panel, width = 10, height = 7, units = "in", dpi = 600)

##------------------------------------------------------------------------------
## Post-decontam filtering & dataset finalization ####
##------------------------------------------------------------------------------
# Filtering is applied in the following order to ensure read depth statistics
# reflect only true biological samples with clean OTU tables:
#
#   1. Remove controls/blanks from metadata & OTU table
#   2. Remove non-fungal OTUs (confirmed non-fungal kingdoms only;
#        unassigned OTUs are RETAINED for downstream BLAST-based reclassification)
#   3. Second abundance filter (remove low-abundance OTUs post-decontam)
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
## Step 1: Remove controls, blanks, and samples outside the spike-in QC range ####
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

cat("--- Step 1: Remove controls/blanks/samples outside spike-in QC range ---\n")
cat("  Samples retained:", nrow(otu.working),
    "| Removed:", nrow(otu.mat.abund.decontam) - nrow(otu.working), "\n")
cat("  OTUs:", ncol(otu.working),
    "| Total reads:", format(sum(otu.working), big.mark = ","), "\n")
cat("  Host species represented:", n_distinct(metadata.final$host_taxon), "\n\n")

##------------------------------------------------------------------------------
## Step 2: Remove confirmed non-fungal OTUs ####
##------------------------------------------------------------------------------
# IMPORTANT: Two categories of OTUs are intentionally RETAINED at this step
# for downstream BLAST-based reclassification (handled in Script 2):
#
#   1. Kingdom-unassigned OTUs (NA, blank, or "Unassigned") — may include
#      true fungi that the UNITE classifier could not confidently assign
#   2. Kingdom == "Fungi", Phylum == "Fungi_phy_incertae_sedis" — likely fungi
#      but with uncertain higher-level placement; BLAST can sometimes resolve these
#
# Only OTUs with a CONFIRMED non-fungal kingdom assignment are removed here.

tax.working <- tax_clean_filtered %>%
  filter(otu %in% colnames(otu.working))

# Confirmed non-fungal: has a kingdom assignment that is neither Fungi nor unassigned
non_fungal_confirmed <- tax.working %>%
  filter(!is.na(Kingdom) &
           Kingdom != "" &
           Kingdom != "Unassigned" &
           Kingdom != "Fungi")

cat("--- Step 2: Remove confirmed non-fungal OTUs ---\n")
if (nrow(non_fungal_confirmed) > 0) {
  cat("  Confirmed non-fungal OTUs removed:", nrow(non_fungal_confirmed), "\n")
  cat("  Kingdoms:\n")
  print(table(non_fungal_confirmed$Kingdom))
  otu.working <- otu.working %>% select(-all_of(non_fungal_confirmed$otu))
  tax.working <- tax.working %>% filter(!otu %in% non_fungal_confirmed$otu)
} else {
  cat("  No confirmed non-fungal OTUs detected\n")
}

# Report what is being carried forward to BLAST
unassigned_otus     <- tax.working %>%
  filter(is.na(Kingdom) | Kingdom == "" | Kingdom == "Unassigned")
incertae_sedis_otus <- tax.working %>%
  filter(Kingdom == "Fungi" & Phylum == "Fungi_phy_incertae_sedis")

cat("  OTUs retained for BLAST reclassification:\n")
cat("    Kingdom-unassigned:", nrow(unassigned_otus), "\n")
cat("    Fungi_phy_incertae_sedis:", nrow(incertae_sedis_otus), "\n")
cat("    Total queued for BLAST:", nrow(unassigned_otus) + nrow(incertae_sedis_otus), "\n")

# Remove any OTUs now at zero reads
zero_read_otus <- colSums(otu.working) == 0
if (any(zero_read_otus)) {
  cat("  Removing", sum(zero_read_otus), "zero-read OTUs after taxonomic cleaning\n")
  otu.working <- otu.working[, !zero_read_otus]
  tax.working <- tax.working %>% filter(!otu %in% names(which(zero_read_otus)))
}

cat("  After Step 2:\n")
cat("  Samples:", nrow(otu.working),
    "| OTUs:", ncol(otu.working),
    "| Total reads:", format(sum(otu.working), big.mark = ","), "\n")
cat("  Host species represented:", n_distinct(metadata.final$host_taxon), "\n\n")

##------------------------------------------------------------------------------
## Step 3: Second abundance filter ####
##------------------------------------------------------------------------------
post_abund_threshold <- 150

otus_before_abund2 <- ncol(otu.working)
otu.working <- otu.working %>%
  select(where(~ is.numeric(.x) && sum(.x) > post_abund_threshold))

cat("--- Step 3: Second abundance filter (>", post_abund_threshold, "total reads) ---\n")
cat("  OTUs retained:", ncol(otu.working),
    "| Removed:", otus_before_abund2 - ncol(otu.working), "\n")
cat("  Samples:", nrow(otu.working),
    "| Total reads:", format(sum(otu.working), big.mark = ","), "\n")
cat("  Host species represented:", n_distinct(metadata.final$host_taxon), "\n\n")

tax.working <- tax.working %>% filter(otu %in% colnames(otu.working))

##------------------------------------------------------------------------------
## Assign working objects & interim summary ####
##------------------------------------------------------------------------------
otu.final  <- otu.working
tax.final  <- tax.working

cat("========== INTERIM DATASET SUMMARY (pre-BLAST) ==========\n")
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
cat("==========================================================\n\n")

##------------------------------------------------------------------------------
## Build interim phyloseq object ####
##------------------------------------------------------------------------------
# This object is passed to Script 2 for BLAST integration.
# It retains unassigned and Fungi_phy_incertae_sedis OTUs intentionally.
physeq_preblast <- phyloseq(
  otu_table(otu.final, taxa_are_rows = FALSE),
  tax_table(as.matrix(tax.final %>% column_to_rownames("otu"))),
  sample_data(metadata.final %>% column_to_rownames("sample_name"))
)

saveRDS(physeq_preblast, file = file.path(dir_processed, "ITS_physeq_absolute_preblast.rds"))
cat("Interim phyloseq saved:", file.path(dir_processed, "ITS_physeq_absolute_preblast.rds"), "\n\n")

##------------------------------------------------------------------------------
## Part 4: Extract unknown/uncertain OTU sequences for BLAST ####
##------------------------------------------------------------------------------
# Identifies OTUs that require secondary BLAST-based reclassification:
#   - Kingdom-unassigned OTUs (NA, blank, "Unassigned")
#   - Fungi_phy_incertae_sedis OTUs (likely fungi but uncertain placement)
# Extracts their representative sequences from the QIIME2 FASTA output and
# writes sorted FASTA files for use as BLAST query input.

library(Biostrings)  # for reading/writing FASTA files

cat("========== EXTRACTING SEQUENCES FOR BLAST ==========\n")

# Path to QIIME2 representative sequences FASTA
fasta_file <- file.path(dir_raw, "clustered_sequences_dna-sequences.fasta")

if (!file.exists(fasta_file)) {
  stop("ERROR: QIIME2 FASTA file not found at:\n", fasta_file,
       "\nUpdate the fasta_file path before proceeding.")
}

# Read all representative sequences
all_seqs <- readDNAStringSet(fasta_file)
cat("Loaded", length(all_seqs), "sequences from QIIME2 FASTA\n")

# Add "Otu_" prefix to FASTA names to match phyloseq OTU ID format
fasta_names_original    <- names(all_seqs)
fasta_names_with_prefix <- paste0("Otu_", fasta_names_original)

# Identify OTUs to send to BLAST from the final taxonomy table
# Category 1: kingdom-unassigned
blast_targets_unassigned <- tax.final %>%
  filter(is.na(Kingdom) | Kingdom == "" | Kingdom == "Unassigned") %>%
  pull(otu)

# Category 2: Fungi_phy_incertae_sedis
blast_targets_incertae <- tax.final %>%
  filter(Kingdom == "Fungi" & Phylum == "Fungi_phy_incertae_sedis") %>%
  pull(otu)

blast_targets_all <- unique(c(blast_targets_unassigned, blast_targets_incertae))

cat("OTUs queued for BLAST:\n")
cat("  Kingdom-unassigned:", length(blast_targets_unassigned), "\n")
cat("  Fungi_phy_incertae_sedis:", length(blast_targets_incertae), "\n")
cat("  Total:", length(blast_targets_all), "\n\n")

# Match BLAST target OTU IDs to FASTA sequences
matched_indices <- which(fasta_names_with_prefix %in% blast_targets_all)

if (length(matched_indices) == 0) {
  stop("ERROR: No sequences matched between phyloseq OTU IDs and FASTA names.\n",
       "Check that OTU ID formats are consistent (both should have 'Otu_' prefix).")
}

matched_seqs <- all_seqs[matched_indices]
names(matched_seqs) <- fasta_names_with_prefix[matched_indices]

cat(sprintf("Matched %d / %d OTU sequences (%.1f%%)\n",
            length(matched_seqs), length(blast_targets_all),
            100 * length(matched_seqs) / length(blast_targets_all)))

if (length(matched_seqs) < length(blast_targets_all)) {
  cat("Warning:", length(blast_targets_all) - length(matched_seqs),
      "OTUs could not be matched — may have been filtered post-clustering\n")
}

# Sort sequences by total abundance (descending) so most important OTUs are first
otu_abundance_all <- colSums(otu.final)
abundance_targets <- otu_abundance_all[names(otu_abundance_all) %in% names(matched_seqs)]
abundance_sorted  <- sort(abundance_targets, decreasing = TRUE)
blast_seqs_sorted <- matched_seqs[names(abundance_sorted)]

# Write FASTA outputs
# A. Full set — use this for the actual BLAST run
writeXStringSet(blast_seqs_sorted,
                filepath = file.path(dir_tables, "ITS_all_blast_targets.fasta"),
                format = "fasta")
cat(sprintf("Saved full BLAST query FASTA: ITS_all_blast_targets.fasta (%d sequences)\n",
            length(blast_seqs_sorted)))

# Write metadata table for BLAST targets
blast_target_metadata <- data.frame(
  OTU_ID             = names(abundance_sorted),
  Total_Reads        = as.vector(abundance_sorted),
  Pct_of_Targets     = round(100 * abundance_sorted / sum(abundance_sorted), 3),
  Cumulative_Pct     = round(cumsum(100 * abundance_sorted / sum(abundance_sorted)), 3),
  Sequence_Length_bp = width(blast_seqs_sorted),
  BLAST_Category     = ifelse(names(abundance_sorted) %in% blast_targets_incertae,
                              "Fungi_phy_incertae_sedis", "Kingdom_unassigned")
) %>%
  left_join(tax.final %>% rename(OTU_ID = otu), by = "OTU_ID")

write.csv(blast_target_metadata, file = file.path(dir_tables, "ITS_blast_target_metadata.csv"), row.names = FALSE)
cat("Saved BLAST target metadata: ITS_blast_target_metadata.csv\n\n")

# Sequence length summary
cat("Sequence length distribution for BLAST targets:\n")
cat(sprintf("  Range: %d – %d bp\n", min(width(blast_seqs_sorted)), max(width(blast_seqs_sorted))))
cat(sprintf("  Mean: %.1f bp | Median: %.1f bp\n",
            mean(width(blast_seqs_sorted)), median(width(blast_seqs_sorted))))

short_seqs <- sum(width(blast_seqs_sorted) < 100)
if (short_seqs > 0)
  cat(sprintf("  Warning: %d sequences < 100 bp (may be low quality)\n", short_seqs))
long_seqs <- sum(width(blast_seqs_sorted) > 1000)
if (long_seqs > 0)
  cat(sprintf("  Note: %d sequences > 1000 bp\n", long_seqs))

###############################################################################
# ======================= BLAST BREAK POINT ===================================
#
# Script 1 is complete. Before running Script 2, perform the BLAST search
# in terminal using the FASTA output above.
#
# Required inputs:
#   Query:    [date]_all_blast_targets.fasta  (full run)
#             [date]_top100_blast_targets.fasta  (test run)
#   Database: UNITE+INSD all-eukaryotes (build with makeblastdb if not already done)
#
# Build database (one-time setup):
#   makeblastdb -in UNITE_public_all_19.02.2025.fasta -dbtype nucl -out UNITE_db
#
# Run BLAST:
#   blastn \
#     -query [date]_all_blast_targets.fasta \
#     -db UNITE_db \
#     -outfmt "6 qseqid sseqid pident qcovs evalue bitscore" \
#     -evalue 1e-6 \
#     -max_target_seqs 5 \
#     -num_threads 4 \
#     > [date]_blastn_all_results.tsv
#
# Then open Script 2: ITS_BLAST_integration.R
# Load: [date]_ITS_physeq_absolute_preblast.rds + [date]_blastn_all_results.tsv
###############################################################################

###############################################################################
###############################################################################
###############################################################################
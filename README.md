# Herpetofauna Gut Microbiome Phylosymbiosis

Code accompanying:
> Rurik et al., in prep

This repository contains the R scripts used to process, clean, and analyze 16S rRNA
(bacterial) and ITS1 rDNA (fungal) amplicon sequencing data characterizing gut
microbiomes across a phylogenetically diverse set of reptile and amphibian hosts
(snakes, lizards, turtles, tortoises, crocodilians, frogs, toads, and salamanders).
Absolute microbial abundances were quantified using the DspikeIn spike-in
normalization protocol (Ghotbi et al. 2025, *The ISME Journal*).

## Data availability
- Raw sequencing reads: deposited in NCBI SRA under accession PRJNA932855
- Raw data: (OTU/count tables, taxonomy tables, sample metadata,
  phylogenetic trees, and other files used as direct input to these scripts)
  archived on FigShare with DOI [https://doi.org/10.6084/m9.figshare.32927786]
- Code (this repository): archived at the time of publication; see [Zenodo DOI, in progress]

## System requirements

**Note:** the versions listed below reflect what was used to generate the results
in this paper and are provided for exact reproducibility; they are not strict
requirements. These analyses do not rely on version-specific behavior and are
expected to run correctly on reasonably current versions of R and these
packages, subject to normal package-to-package compatibility.

**Bioinformatic processing (upstream of this repository):**
QIIME2 v2025.7, Cutadapt v5.1, DADA2 v1.26.0, FIGARO v1.1.2, VSEARCH v2.22.1,
ITSXpress v2.1.4, SILVA v138, UNITE v10 / UNITE+INSDC v10, BLAST+, decontam v1.28.0

**R environment (this repository):** R v4.5.1, using:
- `phyloseq` v1.52.0
- `vegan` v2.7.2
- `ecodist` v2.1.3
- `ape` v5.8-1
- `geosphere` v1.5-20
- `cluster` v2.1.8.1
- `phytools` v2.5-2
- `phylosignal` v1.3.1
- `SpiecEasi` v1.1.1
- `igraph`
- `rstatix` v0.7.2
- `DspikeIn` (see Ghotbi et al. 2025)
- `tidyverse`, `ggtree`, `here`

Network visualizations were additionally rendered in Cytoscape (not required to
reproduce underlying statistical results, only figure layout).

Tested on: macOS 26.5.2 (Tahoe)
No non-standard or specialized hardware is required.

## Installation guide
1. Clone this repository: `git clone https://github.com/arurik/herpetofauna-gut-microbiome-phylosymbiosis`
2. Install R (v4.5.1) from CRAN
3. Install required packages from CRAN (`vegan`, `ecodist`, `ape`, `geosphere`,
   `cluster`, `phytools`, `phylosignal`, `igraph`, `rstatix`, `phyloseq`,
   `tidyverse`, `ggtree`, `here`)
4. Install `SpiecEasi` and `DspikeIn` via GitHub (`devtools::install_github(...)`)
5. Typical install time: ~15 minutes on a normal desktop

## Demo
A separate demo dataset was not created for this repository. These scripts are
standard analysis code (not standalone software) designed to run directly on
this study's own published dataset, which is fully public (Figshare DOI:
10.6084/m9.figshare.32927786; NCBI SRA: PRJNA932855). Running the pipeline on
this published dataset, as described below, serves the same verification
purpose a separate demo dataset would.

## Repository contents
All analysis code is in `scripts/`, numbered in the order they are meant to be run:

| Script | Purpose |
|---|---|
| `01a_dspikein_16S` / `01b_dspikein_ITS` | Spike-in based absolute abundance quantification (DspikeIn) |
| `02a_decontam_16S` | Contaminant removal (16S), read-depth filtering |
| `02b_decontam_ITS_blast1` / `02c_decontam_ITS_blast2_integration` | Contaminant removal (ITS), BLAST-based taxonomic verification |
| `03_phyloseq_cleaning` | Metadata cleaning and correction; phyloseq object finalization |
| `04a_host_metadata_table` / `04b_top_microbial_taxa_table` | Supplementary host metadata and top-taxa summary tables |
| `05_host_microbe_tree_figs` | Host phylogeny + microbiome composition figure |
| `06_full_nmds_permanova` | Beta-diversity ordinations and PERMANOVA, full dataset |
| `07_sample_type_test_split` | Sample-type (fecal/cloacal/lower-GI) comparison and dataset splitting |
| `08_alpha_div` | Alpha diversity analyses |
| `09_core_mb` | Core microbiome analysis |
| `10_mrm_phylosymbiosis` | Multiple regression on distance matrices (MRM), phylosymbiosis testing |
| `11_physignal_per_genera` | Phylogenetic signal (Pagel's λ) in microbial lineage abundances |
| `12_basid_nmds_permanova_test` | Sensitivity analysis: influence of *Basidiobolus* on fungal ordination |
| `13_network_analyses` | Cross-kingdom bacterial-fungal association network inference (SPIEC-EASI), topological role classification, and network-phylosymbiosis (Mantel) tests |

## Instructions for use
1. Download the published dataset from Figshare (DOI: 10.6084/m9.figshare.32927786)
2. Place files into the folder structure below
3. Run scripts in `scripts/` in numerical order — each script's inputs are
   produced by an earlier script in the sequence; see each file's header
   comment for a description
4. Running the full pipeline (scripts 01–13) on the complete published dataset
   reproduces all quantitative results, figures, and supplementary tables
   reported in the manuscript. Detailed parameters (filtering thresholds,
   model specifications, statistical test settings) for each analysis step
   are described in the manuscript's Methods section.

These scripts expect a local project structure (not included in this repository)
populated from the data above:
- data/raw/16S/(count table, taxonomy table, metadata CSV)
- data/raw/ITS/(count table, taxonomy table, metadata CSV, BLAST result files)
- data/raw/trees/(.nwk phylogenies used for figures and phylogenetic analyses)
- data/processed/(intermediate/final phyloseq objects; generated by the pipeline)
- output/figures/
- output/tables/

Scripts use the `here` R package to reference these folders relative to the
project root, so paths will resolve correctly as long as this structure exists
alongside the `scripts/` folder and R is run from the project root (e.g., via an
`.Rproj` file, or by calling `here::i_am()` manually).

## Questions
For questions about the code, data, or requests for unmasked GPS coordinates
(withheld from public data/metadata for a subset of samples to protect
collection sites of poaching-sensitive species), contact the corresponding
author, Donald Walker, at Donald.Walker@mtsu.edu.

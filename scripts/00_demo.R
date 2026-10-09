################################################################################
# Demo: MRM Phylosymbiosis Workflow (simulated data, no external downloads)
#
# Runs a reduced version of the MRM test in 10_mrm on an in-memory simulated
# dataset (~15 species) so the method can be inspected in under a minute
# without the full study data. Mirrors 10_mrm's predictor set: phylogeny,
# diet (10-column binary, asymmetric-binary Gower -- see 10_mrm Section 2 for
# the TetrapodTraits-based version used on real data), habitat, and
# geography, regressed against simulated community dissimilarity via
# ecodist::MRM().
#
# Alexander Rurik
################################################################################

suppressPackageStartupMessages({
  library(ape)
  library(vegan)
  library(ecodist)
  library(cluster)
})

set.seed(52325)
n_spp <- 15

# ── Simulate a host phylogeny and host-trait data ────────────────────────────
host_tree <- rtree(n_spp)
host_tree$tip.label <- paste0("Species_", seq_len(n_spp))

# 10-column binary diet matrix (EltonTraits-style, as used for the real
# TetrapodTraits-derived predictor in 10_mrm Section 2)
diet_cols <- c("Inv", "Vend", "Vect", "Vfish", "Vunk",
               "Scav", "Fruit", "Nect", "Seed", "Plant")
diet_mat  <- matrix(rbinom(n_spp * length(diet_cols), 1, 0.3),
                     nrow = n_spp, dimnames = list(host_tree$tip.label, diet_cols))

habitat   <- sample(c("Terrestrial", "Aquatic", "Arboreal"), n_spp, replace = TRUE)
geo_coords <- data.frame(lat = runif(n_spp, 25, 45), lon = runif(n_spp, -100, -80))
rownames(geo_coords) <- host_tree$tip.label

# ── Simulate a community dissimilarity matrix with a phylogenetic signal ────
# (toy relationship only -- for demonstrating the workflow, not a real result)
phylo_dist <- cophenetic.phylo(host_tree)
comm_dist  <- phylo_dist / max(phylo_dist) + matrix(runif(n_spp^2, 0, 0.3), n_spp)
comm_dist  <- as.dist((comm_dist + t(comm_dist)) / 2)

# ── Build predictor distance matrices ────────────────────────────────────────
dist_phylo   <- as.dist(phylo_dist / max(phylo_dist))
dist_diet    <- daisy(diet_mat, metric = "gower", type = list(asymm = seq_along(diet_cols)))
dist_habitat <- daisy(data.frame(habitat), metric = "gower")
dist_geo     <- as.dist(geosphere_dist <- as.matrix(dist(geo_coords)))  # Euclidean proxy; real script uses geosphere::distm()

# ── MRM ───────────────────────────────────────────────────────────────────────
cat("=== Demo MRM: simulated community dissimilarity ~ phylogeny + diet + habitat + geography ===\n")
fit <- MRM(comm_dist ~ dist_phylo + dist_diet + dist_habitat + dist_geo, nperm = 999)
print(fit)

cat("\nDemo complete. For the full analysis on the study's real data, see 10_mrm\n",
    "(and the data availability section of the README for obtaining the input files).\n")

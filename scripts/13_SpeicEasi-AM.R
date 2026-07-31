# ============================================================
# SPIEC-EASI host-genus networks
# Figure: host phylogeny, network turnover, and
# host-specific SPIEC-EASI networks
# Mitra Ghotbi 2026
# ============================================================

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

packages <- c(
  "ape",
  "phyloseq",
  "vegan",
  "dplyr",
  "tibble",
  "readr",
  "purrr",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "geosphere",
  "ggtree",
  "igraph",
  "ggraph",
  "tidygraph",
  "grid"
)

missing_packages <- packages[!packages %in% rownames(installed.packages())]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

project_dir <- "~/Library/CloudStorage/GoogleDrive-mitra.ghotbi@gmail.com/My Drive/Alex-NetPhylosim"
setwd(project_dir)

out_dir <- "figure5_complete_spiec_all_genera"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

host_tree_path <- "260506_host_genera_tree.nwk"

spiec_results_path <- file.path(
  "spieceasi_host_genus_networks_high_stringency",
  "spiec_host_genus_results_all_high_stringency.rds")

hosts_to_plot <- c(
  "Plethodon",
  "Desmognathus",
  "Eurycea",
  "Ambystoma",
  "Lithobates",
  "Anolis")

if (!exists("ps_network")) {
  stop("`ps_network` must exist in the R environment.")
}

# ------------------------------------------------------------
# 2. Load networks
# ------------------------------------------------------------

spiec_results <- readRDS(spiec_results_path)

stopifnot(all(hosts_to_plot %in% names(spiec_results)))

edge_counts <- purrr::map_dbl(
  spiec_results,
  \(x) igraph::ecount(x$graph)
)

node_counts <- purrr::map_dbl(
  spiec_results,
  \(x) igraph::vcount(x$graph)
)

print(edge_counts)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

get_sample_df <- function(ps) {
  data.frame(phyloseq::sample_data(ps), check.names = FALSE)
}

get_edge_table <- function(graph) {
  igraph::as_data_frame(graph, what = "edges") |>
    dplyr::mutate(
      edge_id = ifelse(
        from < to,
        paste(from, to, sep = "__"),
        paste(to, from, sep = "__")
      ),
      value = 1
    ) |>
    dplyr::select(edge_id, value)
}

theme_nature <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 4),
      legend.position = "none"
    )
}

prepare_complete_graph <- function(graph) {
  genus_attr <- igraph::vertex_attr(graph, "Genus")
  kingdom_attr <- igraph::vertex_attr(graph, "Kingdom")
  
  if (is.null(genus_attr)) {
    genus_attr <- rep(NA_character_, igraph::vcount(graph))
  }
  
  if (is.null(kingdom_attr)) {
    kingdom_attr <- rep(NA_character_, igraph::vcount(graph))
  }
  
  igraph::V(graph)$is_basidiobolus <-
    grepl("Basidiobolus", genus_attr, ignore.case = TRUE) |
    grepl("Basidiobolus", igraph::V(graph)$name, ignore.case = TRUE)
  
  igraph::V(graph)$kingdom_group <- factor(
    dplyr::case_when(
      igraph::V(graph)$is_basidiobolus ~ "Basidiobolus",
      tolower(as.character(kingdom_attr)) %in% c("bacteria", "bacterium") ~ "Bacteria",
      tolower(as.character(kingdom_attr)) == "fungi" ~ "Other fungi",
      TRUE ~ "Other/unknown"
    ),
    levels = c("Bacteria", "Basidiobolus", "Other fungi", "Other/unknown")
  )
  
  igraph::V(graph)$degree_plot <- igraph::degree(graph)
  
  if (!"weight" %in% igraph::edge_attr_names(graph)) {
    igraph::E(graph)$weight <- 1
  }
  
  igraph::E(graph)$original_weight <- igraph::E(graph)$weight
  
  igraph::E(graph)$edge_sign <- factor(
    ifelse(igraph::E(graph)$original_weight >= 0, "Positive", "Negative"),
    levels = c("Positive", "Negative")
  )
  
  graph
}

plot_complete_network <- function(
    graph,
    host_genus_name,
    degree_limits,
    kingdom_cols,
    edge_cols,
    seed = 10010
) {
  set.seed(seed)
  
  graph <- prepare_complete_graph(graph)
  
  # Use a separate layout graph so negative SPIEC-EASI weights do not break
  # shortest-path-based layout calculations.
  graph_layout <- graph
  igraph::E(graph_layout)$weight <- abs(igraph::E(graph_layout)$original_weight)
  
  ggraph::ggraph(graph_layout, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_colour = edge_sign),
      linewidth = 0.22,
      alpha = 0.28,
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill = kingdom_group, size = degree_plot),
      shape = 21,
      color = "white",
      stroke = 0.15,
      alpha = 0.88,
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      data = function(x) dplyr::filter(x, is_basidiobolus),
      ggplot2::aes(size = degree_plot),
      shape = 21,
      fill = "#D73027",
      color = "black",
      stroke = 0.45,
      show.legend = FALSE
    ) +
    ggraph::scale_edge_colour_manual(values = edge_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = kingdom_cols, drop = FALSE) +
    ggplot2::scale_size_continuous(
      range = c(1.2, 5.2),
      limits = degree_limits
    ) +
    ggraph::theme_graph(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 10
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 8
      ),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = host_genus_name,
      subtitle = paste0(
        igraph::vcount(graph),
        " nodes, ",
        igraph::ecount(graph),
        " edges"
      )
    )
}

make_manual_legend <- function(
    kingdom_cols,
    edge_cols,
    degree_values
) {
  node_legend <- tibble::tibble(
    x = c(1.1, 2.55, 4.0),
    y = 2,
    label = c("Bacteria", "Basidiobolus", "Other fungi"),
    fill = c("Bacteria", "Basidiobolus", "Other fungi"))
  
  edge_legend <- tibble::tibble(
    x = c(5.7, 7.25),
    xend = c(6.35, 7.9),
    y = c(2, 2),
    yend = c(2, 2),
    label = c("Positive", "Negative"),
    edge_type = c("Positive", "Negative"))
  
  degree_legend <- tibble::tibble(
    x = c(9.15, 10.05, 10.95),
    y = 2,
    degree = degree_values )
  
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.15,
      y = 2,
      label = "Node identity",
      fontface = "bold",
      hjust = 0,
      size = 3
    ) +
    ggplot2::geom_point(
      data = node_legend,
      ggplot2::aes(x, y, fill = fill),
      shape = 21,
      size = 4,
      color = "white",
      stroke = 0.25
    ) +
    ggplot2::geom_text(
      data = node_legend,
      ggplot2::aes(x + 0.18, y, label = label),
      hjust = 0,
      size = 2.8 ) +
    ggplot2::annotate(
      "text",
      x = 4.95,
      y = 2,
      label = "Association",
      fontface = "bold",
      hjust = 0,
      size = 3
    ) +
    ggplot2::geom_segment(
      data = edge_legend,
      ggplot2::aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend,
        color = edge_type
      ),
      linewidth = 0.8 ) +
    ggplot2::geom_text(
      data = edge_legend,
      ggplot2::aes(x = xend + 0.1, y = y, label = label),
      hjust = 0,
      size = 2.8
    ) +
    ggplot2::annotate(
      "text",
      x = 8.45,
      y = 2,
      label = "Degree",
      fontface = "bold",
      hjust = 0,
      size = 3) +
    ggplot2::geom_point(
      data = degree_legend,
      ggplot2::aes(x, y, size = degree),
      shape = 21,
      fill = "grey75",
      color = "black",
      stroke = 0.25
    ) +
    ggplot2::geom_text(
      data = degree_legend,
      ggplot2::aes(x, y - 0.45, label = degree),
      size = 2.6 ) +
    ggplot2::scale_fill_manual(values = kingdom_cols, guide = "none") +
    ggplot2::scale_color_manual(values = edge_cols, guide = "none") +
    ggplot2::scale_size_continuous(
      range = c(2.2, 6),
      limits = range(degree_values),
      guide = "none"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 12),
      ylim = c(1.25, 2.6),
      clip = "off"
    ) +
    ggplot2::theme_void()
}

# ------------------------------------------------------------
# 4. Colors
# ------------------------------------------------------------

host_cols <- c(
  "Plethodon" = "#5E4FA2",
  "Desmognathus" = "#3288BD",
  "Lithobates" = "#66C2A5",
  "Eurycea" = "#E6AB02",
  "Anolis" = "#F46D43",
  "Ambystoma" = "#9E0142"
)

kingdom_cols <- c(
  "Bacteria" = "#4575B4",
  "Basidiobolus" = "#D73027",
  "Other fungi" = "#E6AB02",
  "Other/unknown" = "grey75"
)

edge_cols <- c(
  "Positive" = "#2C7BB6",
  "Negative" = "#D7191C")

# ------------------------------------------------------------
# 5. SPIEC-EASI edge-set distance
# ------------------------------------------------------------

spiec_edge_presence <- purrr::imap(
  spiec_results,
  function(x, host_name) {
    get_edge_table(x$graph) |>
      dplyr::rename(!!host_name := value)
  }
) |>
  purrr::reduce(dplyr::full_join, by = "edge_id")

spiec_edge_presence[is.na(spiec_edge_presence)] <- 0

spiec_edge_matrix <- spiec_edge_presence |>
  tibble::column_to_rownames("edge_id") |>
  as.matrix()

spiec_jaccard_dist <- vegan::vegdist(
  t(spiec_edge_matrix),
  method = "jaccard",
  binary = TRUE)

network_hosts <- labels(spiec_jaccard_dist)
network_dist_mat <- as.matrix(spiec_jaccard_dist)

readr::write_csv(
  as.data.frame(network_dist_mat) |>
    tibble::rownames_to_column("host_genus"),
  file.path(out_dir, "complete_spiec_jaccard_distance_matrix.csv"))

# ------------------------------------------------------------
# 6. Host phylogeny
# ------------------------------------------------------------

host_tree <- ape::read.tree(host_tree_path)

if (!setequal(network_hosts, host_tree$tip.label)) {
  stop("Network genera and host-tree tips do not match.")
}

host_tree <- ape::keep.tip(host_tree, network_hosts)

phylo_dist_mat <- ape::cophenetic.phylo(host_tree)
phylo_dist_mat <- phylo_dist_mat[network_hosts, network_hosts]
phylo_dist <- as.dist(phylo_dist_mat)

stopifnot(identical(labels(phylo_dist), labels(spiec_jaccard_dist)))

readr::write_csv(
  as.data.frame(phylo_dist_mat) |>
    tibble::rownames_to_column("host_genus"),
  file.path(out_dir, "host_phylogenetic_distance_matrix.csv"))

# ------------------------------------------------------------
# 7. Geographic distance
# ------------------------------------------------------------

meta_df <- get_sample_df(ps_network)

host_n <- meta_df |>
  dplyr::mutate(host_genus = as.character(host_genus)) |>
  dplyr::filter(host_genus %in% network_hosts) |>
  dplyr::count(host_genus, name = "n_samples")

geo_centroids <- meta_df |>
  dplyr::mutate(
    host_genus = as.character(host_genus),
    gps_n = as.numeric(gps_n),
    gps_w = as.numeric(gps_w)
  ) |>
  dplyr::filter(
    host_genus %in% network_hosts,
    is.finite(gps_n),
    is.finite(gps_w)
  ) |>
  dplyr::group_by(host_genus) |>
  dplyr::summarise(
    mean_lat = mean(gps_n, na.rm = TRUE),
    mean_lon = mean(gps_w, na.rm = TRUE),
    n_geo_samples = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::slice(match(network_hosts, host_genus))

geo_dist_mat <- geosphere::distm(
  geo_centroids[, c("mean_lon", "mean_lat")],
  fun = geosphere::distHaversine) / 1000

rownames(geo_dist_mat) <- geo_centroids$host_genus
colnames(geo_dist_mat) <- geo_centroids$host_genus

geo_dist <- as.dist(geo_dist_mat)

stopifnot(identical(labels(geo_dist), labels(phylo_dist)))

readr::write_csv(
  geo_centroids,
  file.path(out_dir, "host_genus_geographic_centroids.csv"))

readr::write_csv(
  as.data.frame(geo_dist_mat) |>
    tibble::rownames_to_column("host_genus"),
  file.path(out_dir, "host_geographic_distance_matrix_km.csv"))

# ------------------------------------------------------------
# 8. Mantel tests
# ------------------------------------------------------------

set.seed(10010)

mantel_res <- vegan::mantel(
  xdis = phylo_dist,
  ydis = spiec_jaccard_dist,
  method = "pearson",
  permutations = 9999)

partial_mantel_res <- vegan::mantel.partial(
  xdis = spiec_jaccard_dist,
  ydis = phylo_dist,
  zdis = geo_dist,
  method = "pearson",
  permutations = 9999)

mantel_summary <- tibble::tibble(
  test = c("Mantel", "Partial Mantel"),
  network_method = "SPIEC-EASI",
  network_scope = "complete network",
  statistic_r = c(
    unname(mantel_res$statistic),
    unname(partial_mantel_res$statistic)),
  p_value = c(
    mantel_res$signif,
    partial_mantel_res$signif),
  permutations = c(
    mantel_res$permutations,
    partial_mantel_res$permutations))

readr::write_csv(
  mantel_summary,
  file.path(out_dir, "mantel_complete_spiec_network_summary.csv"))

# ------------------------------------------------------------
# 9. Pairwise data for Panel B
# ------------------------------------------------------------

pairwise_df <- expand.grid(
  host_1 = network_hosts,
  host_2 = network_hosts,
  stringsAsFactors = FALSE
) |>
  dplyr::filter(host_1 < host_2) |>
  dplyr::mutate(
    phylogenetic_distance = phylo_dist_mat[cbind(host_1, host_2)],
    network_dissimilarity = network_dist_mat[cbind(host_1, host_2)],
    pair = paste(host_1, host_2, sep = " – ") )

readr::write_csv(
  pairwise_df,
  file.path(out_dir, "pairwise_complete_spiec_phylogeny_distances.csv"))

mantel_label <- paste0(
  "Mantel r = ",
  sprintf("%.2f", mantel_summary$statistic_r[1]),
  ", P = ",
  format.pval(mantel_summary$p_value[1], digits = 2, eps = 1e-4),
  "\nPartial Mantel r = ",
  sprintf("%.2f", mantel_summary$statistic_r[2]),
  ", P = ",
  format.pval(mantel_summary$p_value[2], digits = 2, eps = 1e-4))

y_range <- range(pairwise_df$network_dissimilarity, na.rm = TRUE)
y_pad <- diff(y_range) * 0.18

if (!is.finite(y_pad) || y_pad == 0) {
  y_pad <- 0.005
}

y_limits <- c(
  max(0, y_range[1] - y_pad),
  min(1, y_range[2] + y_pad)
)

# ------------------------------------------------------------
# 10. Panel A: host phylogeny
# ------------------------------------------------------------

tree_label_df <- tibble::tibble(label = host_tree$tip.label) |>
  dplyr::left_join(host_n, by = c("label" = "host_genus")) |>
  dplyr::mutate(label_n = paste0(label, "  n=", n_samples))

tree_xmax <- max(ape::node.depth.edgelength(host_tree), na.rm = TRUE)

p_tree <- ggtree::ggtree(host_tree, linewidth = 0.7) %<+% tree_label_df +
  ggtree::geom_tiplab(
    ggplot2::aes(label = label_n, color = label),
    size = 3.2,
    fontface = "italic",
    align = TRUE,
    linesize = 0.25
  ) +
  ggplot2::scale_color_manual(values = host_cols, guide = "none") +
  ggplot2::xlim(0, tree_xmax * 1.55) +
  ggtree::theme_tree2(base_size = 9) +
  ggplot2::labs(
    title = "Host phylogeny",
    x = "Root-to-tip distance in host tree")

# ------------------------------------------------------------
# 11. Panel B: complete-network turnover
# ------------------------------------------------------------

p_mantel <- ggplot2::ggplot(
  pairwise_df,
  ggplot2::aes(phylogenetic_distance, network_dissimilarity)
) +
  ggplot2::geom_point(size = 3.0, color = "black", alpha = 0.9) +
  ggplot2::geom_smooth(
    method = "loess",
    formula = y ~ x,
    span = 0.9,
    se = TRUE,
    color = "#9E0142",
    fill = "#9E0142",
    alpha = 0.16,
    linewidth = 0.9
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = pair),
    size = 2.1,
    color = "grey25",
    segment.size = 0.2,
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.35,
    seed = 10010
  ) +
  ggplot2::annotate(
    "text",
    x = min(pairwise_df$phylogenetic_distance),
    y = y_limits[2],
    label = mantel_label,
    hjust = 0,
    vjust = 1,
    size = 2.8
  ) +
  ggplot2::scale_y_continuous(limits = y_limits) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_nature(base_size = 9) +
  ggplot2::theme(
    plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5)
  ) +
  ggplot2::labs(
    title = "Phylogenetic signal in complete-network turnover",
    x = "Pairwise host phylogenetic distance",
    y = "SPIEC-EASI complete-network dissimilarity\n(Jaccard edge-set distance)")

# ------------------------------------------------------------
# 12. Panels C-H: complete SPIEC-EASI networks
# ------------------------------------------------------------

complete_graphs <- purrr::map(
  hosts_to_plot,
  \(host_name) spiec_results[[host_name]]$graph)

names(complete_graphs) <- hosts_to_plot

global_degree_max <- max(
  purrr::map_dbl(
    complete_graphs,
    \(graph) max(igraph::degree(graph), na.rm = TRUE)
  ),
  na.rm = TRUE)

degree_limits <- c(1, global_degree_max)

degree_values <- pretty(c(1, global_degree_max), n = 3)
degree_values <- degree_values[degree_values > 0]
degree_values <- degree_values[seq_len(min(3, length(degree_values)))]

network_panels <- purrr::imap(
  complete_graphs,
  \(graph, host_name) {
    plot_complete_network(
      graph = graph,
      host_genus_name = host_name,
      degree_limits = degree_limits,
      kingdom_cols = kingdom_cols,
      edge_cols = edge_cols
    )
  }
)

p_legend <- make_manual_legend(
  kingdom_cols = kingdom_cols,
  edge_cols = edge_cols,
  degree_values = degree_values)

# ------------------------------------------------------------
# 13. ov figure
# ------------------------------------------------------------

fig5_complete <- (p_tree | p_mantel) /
  (
    network_panels[["Plethodon"]] |
      network_panels[["Desmognathus"]] |
      network_panels[["Eurycea"]]
  ) /
  (
    network_panels[["Ambystoma"]] |
      network_panels[["Lithobates"]] |
      network_panels[["Anolis"]]
  ) /
  p_legend +
  patchwork::plot_layout(
    heights = c(1.15, 1, 1, 0.16)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Host evolutionary divergence is associated with microbial network turnover",
    subtitle = "Complete host-specific SPIEC-EASI networks",
    theme = ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 10, hjust = 0)
      ))

ggplot2::ggsave(
  file.path(out_dir, "figure5_complete_spiec_all_genera_networks.pdf"),
  fig5_complete,
  width = 13.5,
  height = 14.0,
  device = cairo_pdf)

ggplot2::ggsave(
  file.path(out_dir, "figure5_complete_spiec_all_genera_networks.png"),
  fig5_complete,
  width = 13.5,
  height = 14.0,
  dpi = 1000,
  limitsize = FALSE)

# ------------------------------------------------------------
# 14. summary
# ------------------------------------------------------------

cat("\nComplete SPIEC-EASI network figure complete.\n")
cat("Output directory: ", out_dir, "\n", sep = "")

 
cat(
  "\nFigure legend clarification:\n",
  "Panel A x-axis shows root-to-tip distance in the genus-level host tree. ",
  "The y-axis represents graphical ordering of taxa and has no quantitative ",
  "interpretation. Panel B x-axis shows pairwise cophenetic phylogenetic ",
  "distance between host genera. Panel B y-axis shows complete SPIEC-EASI ",
  "network dissimilarity quantified as Jaccard edge-set distance. Panels C-H ",
  "show complete host-specific SPIEC-EASI networks. Nodes represent microbial ",
  "taxa, Basidiobolus nodes are highlighted in red, node size is proportional ",
  "to degree, and edge color denotes positive or negative conditional ",
  "associations.\n",
  sep = ""
)

  
sessionInfo()

 

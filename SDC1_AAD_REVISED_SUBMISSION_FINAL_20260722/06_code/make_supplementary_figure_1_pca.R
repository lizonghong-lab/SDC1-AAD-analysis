options(stringsAsFactors = FALSE, scipen = 999)

# Figure contract
# Core conclusion: an outcome-blind high-completeness PCA summarizes global
# proteomic covariance without serving as a biomarker-selection step.
# Evidence: one PC1-versus-PC2 score plot; all technical-QC panels are omitted.
# Archetype: single quantitative scatter plot.
# Review risk: post-hoc disease colouring must not be read as supervised fitting.

package_root <- Sys.getenv("SDC1_R_LIBRARY", unset = "")
if (dir.exists(package_root)) .libPaths(c(package_root, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

root <- normalizePath(
  Sys.getenv(
    "SDC1_REVISION_PACKAGE",
    unset = file.path(getwd(), "SDC1_AAD_REVISED_SUBMISSION_FINAL_20260722")
  ),
  winslash = "/",
  mustWork = TRUE
)

source_dir <- file.path(root, "04_figures", "supplementary", "source_data")
output_dir <- file.path(root, "04_figures", "supplementary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(
  file.path(source_dir, "proteomics_pca_scores.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
variance <- read.csv(
  file.path(source_dir, "proteomics_pca_variance.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pc1 <- variance$percent_variance[variance$principal_component == "PC1"]
pc2 <- variance$percent_variance[variance$principal_component == "PC2"]

scores$group <- factor(
  scores$group,
  levels = c("Healthy control", "AMI", "AAD")
)

palette <- c(
  "Healthy control" = "#6F8F88",
  "AMI" = "#D69032",
  "AAD" = "#B6404C"
)

plot_pca <- ggplot(scores, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 1.65, alpha = 0.72) +
  stat_ellipse(level = 0.80, linewidth = 0.48, show.legend = FALSE) +
  scale_colour_manual(values = palette, drop = FALSE) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pc1),
    y = sprintf("PC2 (%.1f%%)", pc2),
    colour = "Diagnostic group"
  ) +
  theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(size = 7.0, colour = "black"),
    axis.title = element_text(size = 7.5, colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7.0),
    legend.text = element_text(size = 6.8),
    legend.key.width = grid::unit(4.5, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 6, 4, 5)
  )

stem <- file.path(output_dir, "Supplementary_Figure_1_outcome_blind_PCA")
width_mm <- 120
height_mm <- 92
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
print(plot_pca)
dev.off()

grDevices::cairo_pdf(
  paste0(stem, ".pdf"),
  width = width_in,
  height = height_in,
  family = "Arial"
)
print(plot_pca)
dev.off()

ragg::agg_tiff(
  paste0(stem, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(plot_pca)
dev.off()

ragg::agg_png(
  paste0(stem, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300,
  background = "white"
)
print(plot_pca)
dev.off()

message("Created: ", stem)

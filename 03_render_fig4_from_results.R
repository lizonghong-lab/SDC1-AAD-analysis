options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
legacy_rlib <- Sys.getenv(
  "SDC1_RLIB",
  unset = ""
)
if (dir.exists(legacy_rlib)) .libPaths(c(legacy_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
figure_dir <- file.path(root, "04_results", "figures", "main")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(name) read.csv(file.path(result_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
read_source <- function(name) read.csv(file.path(source_dir, name), stringsAsFactors = FALSE, check.names = FALSE)

format_probability <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.0001) return("<0.0001")
  if (value < 0.001) return(sprintf("%.4f", value))
  sprintf("%.3f", value)
}

a_data <- read_source("fig4_panel_a_source_data.csv")
b_data <- read_source("fig4_panel_b_source_data.csv")
htn_data <- read_source("fig4_panels_c_to_f_source_data.csv")
a_means <- read_result("fig4_adjusted_total_sdc1_group_means.csv")
b_means <- read_result("fig4_adjusted_d_dimer_group_means.csv")
htn_means <- read_result("fig4_adjusted_hypertension_means.csv")
a_contrasts <- read_result("fig4_adjusted_total_sdc1_group_contrasts.csv")
b_contrasts <- read_result("fig4_adjusted_d_dimer_group_contrasts.csv")
htn_effects <- read_result("fig4_adjusted_hypertension_effects.csv")
performance <- read_result("fig4_clinically_adjusted_validated_performance.csv")
calibration <- read_source("fig4_clinically_adjusted_calibration.csv")
roc_curves <- read_source("fig4_clinically_adjusted_roc_curves.csv")
roc_auc <- read_result("fig4_clinically_adjusted_roc_auc.csv")

all_group_levels <- c("Healthy control", "AAD", "AMI", "Angina", "PE")
symptomatic_levels <- c("AAD", "AMI", "Angina", "PE")
a_data$group <- factor(a_data$group, levels = all_group_levels)
b_data$group <- factor(b_data$group, levels = symptomatic_levels)
a_means$group <- factor(a_means$group, levels = all_group_levels)
b_means$group <- factor(b_means$group, levels = symptomatic_levels)
htn_data$group <- factor(htn_data$group, levels = symptomatic_levels)
htn_means$group <- factor(htn_means$group, levels = symptomatic_levels)
htn_data$hypertension_label <- factor(htn_data$hypertension_label, levels = c("No hypertension", "Hypertension"))
htn_means$hypertension_label <- factor(htn_means$hypertension_label, levels = c("No hypertension", "Hypertension"))

palette <- c("Healthy control" = "#4B8F8A", "AAD" = "#BF3F4C", "AMI" = "#3F78A8", "Angina" = "#D78B2B", "PE" = "#7B5AA6")
model_palette <- c("Clinical" = "#92999F", "+ D-dimer" = "#397AAE", "+ D-dimer + SDC1" = "#BE4050")
htn_palettes <- list(
  "AAD" = c("No hypertension" = "#E8B7BD", "Hypertension" = "#B83C4A"),
  "AMI" = c("No hypertension" = "#BAD3E7", "Hypertension" = "#397AAE"),
  "Angina" = c("No hypertension" = "#F0D2A4", "Hypertension" = "#D78B2B"),
  "PE" = c("No hypertension" = "#D8C9E7", "Hypertension" = "#7B5AA6")
)

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6.1, colour = "black"),
      plot.title = element_text(size = 7, face = "bold", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 5.8, colour = "#555B60", margin = margin(b = 3)),
      strip.text = element_text(size = 5.9, face = "bold", margin = margin(2.5, 1.5, 2.5, 1.5)),
      strip.background = element_rect(fill = "white", colour = "#222222", linewidth = 0.35),
      legend.text = element_text(size = 5.7),
      plot.margin = margin(5, 5, 4, 5)
    )
}

a_labels <- setNames(
  paste0(c("HC", "AAD", "AMI", "Angina", "PE"), "\n(n=", a_means$n_group[match(all_group_levels, a_means$group)], ")"),
  all_group_levels
)
a_annotation <- merge(
  data.frame(reference = c("Healthy control", "AMI", "Angina", "PE"), stringsAsFactors = FALSE),
  a_contrasts[, c("reference", "bh_q")],
  by = "reference"
)
a_annotation$x <- match(a_annotation$reference, all_group_levels)
a_annotation$y <- max(a_data$total_sdc1) * 1.34
a_annotation$label <- paste0("vs AAD\nq", vapply(a_annotation$bh_q, format_probability, character(1)))

p_a <- ggplot(a_data, aes(group, total_sdc1, colour = group)) +
  geom_jitter(position = position_jitter(width = 0.11, height = 0, seed = 4101), size = 0.60, alpha = 0.42) +
  geom_errorbar(data = a_means, aes(x = group, ymin = lower_display_95, ymax = upper_display_95), inherit.aes = FALSE, width = 0.14, linewidth = 0.55, colour = "#202428") +
  geom_point(data = a_means, aes(x = group, y = estimate_display), inherit.aes = FALSE, shape = 21, size = 2.2, stroke = 0.6, fill = "white", colour = "#202428") +
  geom_text(data = a_annotation, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 1.9, lineheight = 0.9) +
  scale_colour_manual(values = palette) +
  scale_x_discrete(labels = a_labels) +
  scale_y_log10(expand = expansion(mult = c(0.04, 0.23))) +
  labs(
    x = NULL,
    y = "Total-plasma SDC1 (pg/mL, log scale)",
    title = "Total-plasma SDC1",
    subtitle = "Common-covariate adjusted estimates"
  ) +
  theme_pub() +
  theme(legend.position = "none")

b_labels <- setNames(
  paste0(symptomatic_levels, "\n(n=", b_means$n_group[match(symptomatic_levels, b_means$group)], ")"),
  symptomatic_levels
)
b_annotation <- merge(
  data.frame(reference = c("AMI", "Angina", "PE"), stringsAsFactors = FALSE),
  b_contrasts[, c("reference", "bh_q")],
  by = "reference"
)
b_annotation$x <- match(b_annotation$reference, symptomatic_levels)
b_annotation$y <- max(b_data$log2_d_dimer) + 1.1
b_annotation$label <- paste0("vs AAD\nq", vapply(b_annotation$bh_q, format_probability, character(1)))

p_b <- ggplot(b_data, aes(group, log2_d_dimer, colour = group)) +
  geom_jitter(position = position_jitter(width = 0.11, height = 0, seed = 4102), size = 0.62, alpha = 0.44) +
  geom_errorbar(data = b_means, aes(x = group, ymin = lower_log2_95, ymax = upper_log2_95), inherit.aes = FALSE, width = 0.14, linewidth = 0.55, colour = "#202428") +
  geom_point(data = b_means, aes(x = group, y = estimate_log2), inherit.aes = FALSE, shape = 21, size = 2.2, stroke = 0.6, fill = "white", colour = "#202428") +
  geom_text(data = b_annotation, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 1.9, lineheight = 0.9) +
  scale_colour_manual(values = palette) +
  scale_x_discrete(labels = b_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(
    x = NULL,
    y = "log2 D-dimer (mg/L)",
    title = "D-dimer",
    subtitle = "Symptomatic-cohort adjusted estimates"
  ) +
  theme_pub() +
  theme(legend.position = "none")

add_log_bracket <- function(plot, y, label) {
  plot +
    annotate("segment", x = 1, xend = 2, y = y, yend = y, linewidth = 0.35) +
    annotate("segment", x = 1, xend = 1, y = y, yend = y / 1.08, linewidth = 0.35) +
    annotate("segment", x = 2, xend = 2, y = y, yend = y / 1.08, linewidth = 0.35) +
    annotate("text", x = 1.5, y = y * 1.04, label = label, size = 2.0, vjust = 0)
}

make_htn_panel <- function(group_name, show_y_title = FALSE, jitter_seed = 1L) {
  raw <- htn_data[htn_data$group == group_name, , drop = FALSE]
  mean_data <- htn_means[htn_means$group == group_name, , drop = FALSE]
  effect <- htn_effects[htn_effects$group == group_name, , drop = FALSE]
  labels <- setNames(
    paste0(c("No HTN", "HTN"), "\n(n=", mean_data$n_group[match(c("No hypertension", "Hypertension"), mean_data$hypertension_label)], ")"),
    c("No hypertension", "Hypertension")
  )
  bracket_y <- max(raw$total_sdc1) * 1.20
  plot <- ggplot(raw, aes(hypertension_label, total_sdc1, colour = hypertension_label)) +
    geom_jitter(position = position_jitter(width = 0.10, height = 0, seed = jitter_seed), size = 0.63, alpha = 0.46) +
    geom_errorbar(data = mean_data, aes(x = hypertension_label, ymin = lower_display_95, ymax = upper_display_95), inherit.aes = FALSE, width = 0.14, linewidth = 0.55, colour = "#202428") +
    geom_point(data = mean_data, aes(x = hypertension_label, y = estimate_display), inherit.aes = FALSE, shape = 21, size = 2.1, stroke = 0.6, fill = "white", colour = "#202428") +
    scale_colour_manual(values = htn_palettes[[group_name]]) +
    scale_x_discrete(labels = labels) +
    scale_y_log10(expand = expansion(mult = c(0.04, 0.22))) +
    labs(
      x = NULL,
      y = if (show_y_title) "Total-plasma SDC1 (pg/mL, log scale)" else NULL,
      title = group_name
    ) +
    theme_pub() +
    theme(legend.position = "none", plot.margin = margin(5, 3, 4, 3))
  add_log_bracket(plot, bracket_y, paste0("BH q", format_probability(effect$bh_q)))
}

p_c <- make_htn_panel("AAD", show_y_title = TRUE, jitter_seed = 4201)
p_d <- make_htn_panel("AMI", jitter_seed = 4202)
p_e <- make_htn_panel("Angina", jitter_seed = 4203)
p_f <- make_htn_panel("PE", jitter_seed = 4204)

model_levels <- c("Clinical_context", "Clinical_plus_D_dimer", "Clinical_plus_D_dimer_plus_SDC1")
model_display <- c("Clinical", "+ D-dimer", "+ D-dimer + SDC1")

make_roc_panel <- function(population_name, title) {
  curve_data <- roc_curves[roc_curves$population == population_name, , drop = FALSE]
  auc_data <- roc_auc[roc_auc$population == population_name, , drop = FALSE]
  curve_data$model_label <- factor(curve_data$model, levels = model_levels, labels = model_display)
  auc_data <- auc_data[match(model_levels, auc_data$model), , drop = FALSE]
  legend_labels <- setNames(
    sprintf("%s: %.3f (%.3f-%.3f)", model_display, auc_data$estimate, auc_data$lower_95, auc_data$upper_95),
    model_display
  )
  subtitle <- sprintf(
    "n=%d; %d AAD / %d non-AAD",
    auc_data$n_total[1], auc_data$aad_cases[1], auc_data$non_aad_controls[1]
  )
  ggplot(curve_data, aes(false_positive_rate, true_positive_rate, colour = model_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#8F969C", linewidth = 0.4) +
    geom_step(direction = "vh", linewidth = 0.72) +
    scale_colour_manual(values = model_palette, breaks = model_display, labels = legend_labels) +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = expansion(mult = 0)) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = expansion(mult = 0)) +
    labs(
      x = "False-positive rate",
      y = "True-positive rate",
      colour = "AUROC (95% full-pipeline bootstrap CI)",
      title = title,
      subtitle = subtitle
    ) +
    theme_pub() +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.title = element_text(size = 5.7, face = "bold"),
      legend.text = element_text(size = 5.4),
      legend.key.height = grid::unit(2.8, "mm"),
      legend.margin = margin(t = 0),
      plot.title = element_text(size = 7.5),
      aspect.ratio = 0.68,
      plot.margin = margin(5, 8, 4, 5)
    ) +
    guides(colour = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 0.9)))
}

p_g <- make_roc_panel("All symptomatic comparators", "Adjusted ROC: all comparators")
p_h <- make_roc_panel("AAD versus PE", "Adjusted ROC: AAD vs PE")

performance$model_label <- factor(
  performance$model,
  levels = model_levels,
  labels = model_display
)
performance$metric_label <- factor(
  performance$metric,
  levels = c("brier", "log_loss", "c_index"),
  labels = c("Brier score\n(lower)", "Log loss\n(lower)", "C-index / AUROC\n(higher)")
)

p_i <- ggplot(performance, aes(estimate, model_label, colour = model_label)) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0.14, linewidth = 0.5) +
  geom_point(size = 2.2) +
  facet_grid(. ~ metric_label, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = model_palette) +
  scale_x_continuous(n.breaks = 3, expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    x = "Cross-validated estimate (95% full-pipeline bootstrap interval)",
    y = NULL,
    title = "Validated probability performance",
    subtitle = paste0("n=", unique(performance$n_complete), "; lower Brier/log loss and higher C-index are favourable")
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 5.8),
    axis.text.x = element_text(size = 5.4, angle = 35, hjust = 1),
    panel.spacing.x = grid::unit(2.2, "mm"),
    plot.margin = margin(5, 6, 4, 5)
  )

calibration$model_label <- factor(
  calibration$model,
  levels = model_levels,
  labels = model_display
)
for (column in c("cross_validated_smooth_observed", "lower_95", "upper_95")) {
  calibration[[column]] <- pmin(pmax(calibration[[column]], 0), 1)
}

p_j <- ggplot(calibration, aes(predicted_probability, cross_validated_smooth_observed, colour = model_label, fill = model_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#90969C", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.72) +
  scale_colour_manual(values = model_palette) +
  scale_fill_manual(values = model_palette) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  coord_equal() +
  labs(
    x = "Cross-validated predicted AAD probability",
    y = "Smooth observed AAD proportion",
    colour = NULL,
    fill = NULL,
    title = "Cross-validated calibration",
    subtitle = "30 x 10-fold CV; 500 full-pipeline bootstrap bands"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.margin = margin(t = 0), legend.box.margin = margin(0, 0, 0, 0))

top_row <- (p_a | p_b) + plot_layout(widths = c(1.18, 1))
middle_row <- p_c | p_d | p_e | p_f
roc_row <- p_g | p_h
performance_row <- (p_i | p_j) + plot_layout(widths = c(1.50, 1))
figure <- top_row / middle_row / roc_row / performance_row +
  plot_layout(heights = c(0.86, 0.73, 0.94, 1.08)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 8, family = "Arial"))
  )

save_figure <- function(plot, stem, width_mm = 183, height_mm = 242, dpi = 600) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = cairo_pdf, family = "Arial")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, device = svglite::svglite)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300, device = ragg::agg_png)
  ggsave(paste0(stem, ".tiff"), plot, width = width, height = height, dpi = dpi, compression = "lzw", device = ragg::agg_tiff)
}

save_figure(figure, file.path(figure_dir, "Figure_4_clinically_adjusted_performance_and_calibration"))
cat("Figure 4 rendered from locked V9 statistical and source-data tables.\n")

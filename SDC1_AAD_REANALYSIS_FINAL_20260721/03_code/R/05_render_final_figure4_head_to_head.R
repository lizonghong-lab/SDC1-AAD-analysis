options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source(file.path(root, "03_code", "R", "03_render_fig4_from_results.R"), local = FALSE)

standalone_curves <- read_source("fig4_standalone_biomarker_roc_curves.csv")
standalone_auc <- read_result("fig4_standalone_biomarker_roc_auc.csv")
standalone_difference <- read_result("fig4_standalone_biomarker_roc_paired_difference.csv")

biomarker_levels <- c("D_dimer", "SDC1")
biomarker_display <- c("D-dimer", "SDC1")
biomarker_palette <- c("D-dimer" = "#3F78A8", "SDC1" = "#D55E00")

make_standalone_roc_panel <- function(population_name, title) {
  curve_data <- standalone_curves[standalone_curves$population == population_name, , drop = FALSE]
  auc_data <- standalone_auc[standalone_auc$population == population_name, , drop = FALSE]
  difference_data <- standalone_difference[standalone_difference$population == population_name, , drop = FALSE]
  curve_data$biomarker_label <- factor(curve_data$biomarker, levels = biomarker_levels, labels = biomarker_display)
  auc_data <- auc_data[match(biomarker_levels, auc_data$biomarker), , drop = FALSE]
  legend_labels <- setNames(
    sprintf("%s: %.3f (%.3f-%.3f)", biomarker_display, auc_data$estimate, auc_data$lower_95, auc_data$upper_95),
    biomarker_display
  )
  subtitle <- sprintf(
    "n=%d; %d AAD / %d non-AAD; Delta AUROC=%.3f",
    auc_data$n_total[1], auc_data$aad_cases[1], auc_data$non_aad_controls[1], difference_data$difference[1]
  )
  ggplot(curve_data, aes(false_positive_rate, true_positive_rate, colour = biomarker_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#8F969C", linewidth = 0.4) +
    geom_step(direction = "vh", linewidth = 0.82) +
    scale_colour_manual(values = biomarker_palette, breaks = biomarker_display, labels = legend_labels) +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = expansion(mult = 0)) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = expansion(mult = 0)) +
    labs(
      x = "False-positive rate",
      y = "True-positive rate",
      colour = "AUROC (95% paired cluster-bootstrap CI)",
      title = title,
      subtitle = subtitle
    ) +
    theme_pub() +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.title = element_text(size = 5.7, face = "bold"),
      legend.text = element_text(size = 5.5),
      legend.key.height = grid::unit(3.0, "mm"),
      legend.margin = margin(t = 0),
      plot.title = element_text(size = 7),
      aspect.ratio = 0.68,
      plot.margin = margin(5, 8, 4, 5)
    ) +
    guides(colour = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 1.0)))
}

p_g_head <- make_standalone_roc_panel("All symptomatic comparators", "Standalone biomarkers: all comparators")
p_h_head <- make_standalone_roc_panel("AAD versus PE", "Standalone biomarkers: AAD vs PE")

figure4_head_to_head <- top_row / middle_row / (p_g_head | p_h_head) +
  plot_layout(heights = c(0.90, 0.80, 1.10)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold", size = 8, family = "Arial"))
  )

save_figure(
  figure4_head_to_head,
  file.path(figure_dir, "Figure_4_covariate_adjusted_results_and_standalone_biomarker_ROC"),
  width_mm = 183,
  height_mm = 170,
  dpi = 600
)

panel_dir <- file.path(root, "04_results", "figures", "panels", "figure_4")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
panel_list <- list(A = p_a, B = p_b, C = p_c, D = p_d, E = p_e, F = p_f, G = p_g_head, H = p_h_head)
for (panel_name in names(panel_list)) {
  width_mm <- if (panel_name %in% c("G", "H")) 92 else if (panel_name %in% c("A", "B")) 92 else 46
  height_mm <- if (panel_name %in% c("G", "H")) 92 else 64
  save_figure(
    panel_list[[panel_name]],
    file.path(panel_dir, paste0("Figure_4", panel_name)),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = 600
  )
}
cat("Figure 4 A-H rendered with standalone SDC1 versus D-dimer ROC panels; combined models are excluded.\n")

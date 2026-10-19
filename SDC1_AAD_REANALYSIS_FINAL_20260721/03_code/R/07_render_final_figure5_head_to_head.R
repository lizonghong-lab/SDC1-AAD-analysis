options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

# Reuse the locked FINAL_20260721 analysis objects and R-only rendering functions.
source(file.path(root, "03_code", "R", "05_render_final_figure4_head_to_head.R"), local = FALSE)

log_dir_5 <- file.path(root, "04_results", "logs")
dir.create(log_dir_5, recursive = TRUE, showWarnings = FALSE)

p_c_5 <- p_c +
  labs(y = "Total-plasma SDC1\n(pg/mL, log scale)") +
  theme(axis.title.y = element_text(lineheight = 0.9))
middle_row_5 <- p_c_5 | p_d | p_e | p_f

figure5_head_to_head <- top_row / middle_row_5 / (p_g_head | p_h_head) +
  plot_layout(heights = c(0.90, 0.80, 1.10)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 8, family = "Arial", margin = margin(r = 3, b = 2))
    )
  )

save_figure(
  figure5_head_to_head,
  file.path(figure_dir, "Figure_5_clinical_biomarker_head_to_head"),
  width_mm = 183,
  height_mm = 170,
  dpi = 600
)

panel_dir_5 <- file.path(root, "04_results", "figures", "panels", "figure_5")
dir.create(panel_dir_5, recursive = TRUE, showWarnings = FALSE)
panel_list_5 <- list(A = p_a, B = p_b, C = p_c_5, D = p_d, E = p_e, F = p_f, G = p_g_head, H = p_h_head)

for (panel_name in names(panel_list_5)) {
  width_mm <- if (panel_name %in% c("G", "H")) 92 else if (panel_name %in% c("A", "B")) 92 else 46
  height_mm <- if (panel_name %in% c("G", "H")) 92 else 64
  save_figure(
    panel_list_5[[panel_name]],
    file.path(panel_dir_5, paste0("Figure_5", panel_name)),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = 600
  )
}

# Figure-number aliases keep the submission package traceable without changing
# any locked statistical estimate or re-running a different analysis.
source_aliases <- c(
  "fig4_panel_a_source_data.csv" = "fig5_panel_a_source_data.csv",
  "fig4_panel_b_source_data.csv" = "fig5_panel_b_source_data.csv",
  "fig4_panels_c_to_f_source_data.csv" = "fig5_panels_c_to_f_source_data.csv",
  "fig4_standalone_biomarker_roc_curves.csv" = "fig5_standalone_biomarker_roc_curves.csv",
  "fig4_standalone_biomarker_roc_source_data.csv" = "fig5_standalone_biomarker_roc_source_data.csv",
  "fig4_standalone_biomarker_roc_bootstrap_draws.csv" = "fig5_standalone_biomarker_roc_bootstrap_draws.csv"
)
result_aliases <- c(
  "fig4_adjusted_total_sdc1_group_means.csv" = "fig5_adjusted_total_sdc1_group_means.csv",
  "fig4_adjusted_total_sdc1_group_contrasts.csv" = "fig5_adjusted_total_sdc1_group_contrasts.csv",
  "fig4_adjusted_d_dimer_group_means.csv" = "fig5_adjusted_d_dimer_group_means.csv",
  "fig4_adjusted_d_dimer_group_contrasts.csv" = "fig5_adjusted_d_dimer_group_contrasts.csv",
  "fig4_adjusted_hypertension_means.csv" = "fig5_adjusted_hypertension_means.csv",
  "fig4_adjusted_hypertension_effects.csv" = "fig5_adjusted_hypertension_effects.csv",
  "fig4_standalone_biomarker_roc_auc.csv" = "fig5_standalone_biomarker_roc_auc.csv",
  "fig4_standalone_biomarker_roc_paired_difference.csv" = "fig5_standalone_biomarker_roc_paired_difference.csv"
)

copy_aliases <- function(directory, aliases) {
  for (source_name in names(aliases)) {
    source_path <- file.path(directory, source_name)
    target_path <- file.path(directory, unname(aliases[[source_name]]))
    if (!file.exists(source_path)) stop("Missing locked source: ", source_path)
    if (!file.copy(source_path, target_path, overwrite = TRUE, copy.date = TRUE)) {
      stop("Failed to create Figure 5 alias: ", target_path)
    }
  }
}

copy_aliases(source_dir, source_aliases)
copy_aliases(result_dir, result_aliases)

log_lines <- c(
  "Figure 5 final head-to-head rendering",
  "Data lock: FINAL_20260721",
  "Panels A-F: prespecified covariate-adjusted biomarker results",
  "Panels G-H: same-patient standalone SDC1 versus D-dimer ROC comparisons",
  "Combined clinical models, calibration and proper scores are excluded and reserved for Figure 6",
  "Main size: 183 x 170 mm",
  "Exports: editable PDF/SVG, 600 dpi TIFF, PNG preview"
)
writeLines(log_lines, file.path(log_dir_5, "fig5_head_to_head_render_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir_5, "fig5_sessionInfo.txt"), useBytes = TRUE)

cat("Figure 5 A-H rendered as the final standalone biomarker head-to-head figure.\n")

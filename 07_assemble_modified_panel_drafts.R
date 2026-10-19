options(stringsAsFactors = FALSE, scipen = 999)

root_dir <- normalizePath(getwd(), mustWork = TRUE)
project_dir <- file.path(root_dir, "SDC1_AAD_HARRELL_BOOTCAL_20260728")
.libPaths(c(file.path(project_dir, "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

# Rebuild the locked individual plots in the same R session so that the
# grouped drafts use the original vector plot objects rather than raster tiles.
render_script <- Sys.getenv(
  "SDC1_RENDER_SCRIPT",
  unset = "06_render_all_modified_panels.R"
)
source(file.path(project_dir, "03_code", "R", render_script), local = FALSE)

group_run <- Sys.getenv(
  "SDC1_GROUP_RUN",
  unset = "GROUPED_DRAFTS_FOR_AUTHOR_ADJUSTMENT_20260728"
)
group_dir <- file.path(
  project_dir, "04_results", "figures",
  group_run
)
dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

save_group <- function(plot, stem, width_mm, height_mm, dpi = 600L) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "Arial"
  )
  print(plot)
  dev.off()
  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw", background = "white"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  dev.off()
}

compact_roc_panel <- function(plot, scalar_data, population_label) {
  scalar_c <- scalar_data[
    scalar_data$metric == "C-index" & scalar_data$model %in% biomarker_models,
    , drop = FALSE
  ]
  scalar_c <- scalar_c[match(biomarker_models, scalar_c$model), ]
  wrapped_labels <- sprintf(
    "%s:\n%.3f (%.3f-%.3f)",
    model_labels[scalar_c$model], scalar_c$corrected, scalar_c$lower, scalar_c$upper
  )
  plot +
    labs(title = paste0("Joint models:\n", population_label)) +
    scale_colour_manual(
      values = model_colours, breaks = biomarker_models,
      labels = wrapped_labels, name = "Validated AUROC\n(ABCLOC 95% CI)"
    ) +
    theme(
      plot.title = element_text(size = 6.1, lineheight = 0.92),
      legend.position = "inside",
      legend.position.inside = c(0.03, 0.03),
      legend.justification = c(0, 0),
      legend.key.width = grid::unit(3.6, "mm"),
      legend.title = element_text(size = 4.3, lineheight = 0.88),
      legend.text = element_text(size = 3.9, lineheight = 0.86),
      legend.spacing.y = grid::unit(0.15, "mm")
    )
}

prepare_nine_panel_set <- function(panel_list, scalar_data, population_label) {
  panel_list$A <- panel_list$A +
    labs(title = "Smooth optimism-corrected\ncalibration") +
    theme(
      plot.title = element_text(size = 6.25, lineheight = 0.92),
      plot.subtitle = element_text(size = 5.25),
      legend.text = element_text(size = 4.9),
      legend.key.width = grid::unit(5.0, "mm")
    )
  panel_list$B <- compact_roc_panel(panel_list$B, scalar_data, population_label)
  panel_list
}

assemble_nine_panel <- function(panel_list) {
  top <- (panel_list$A | panel_list$B | panel_list$C) +
    plot_layout(widths = c(1.00, 1.00, 1.08))
  middle <- (panel_list$D | panel_list$E | panel_list$F) +
    plot_layout(widths = c(1.00, 1.02, 0.98))
  bottom <- (panel_list$G | panel_list$H | panel_list$I) +
    plot_layout(widths = c(1.00, 1.00, 1.00))
  top / middle / bottom +
    plot_layout(heights = c(1.06, 1.03, 0.97)) &
    theme(plot.margin = margin(3.5, 4.0, 3.5, 4.0))
}

figure5_grouped <- assemble_nine_panel(prepare_nine_panel_set(
  figure5$panels, figure5$performance, "AAD vs symptomatic non-AAD"
))
figure6_grouped <- assemble_nine_panel(prepare_nine_panel_set(
  figure6$panels, figure6$performance, "AAD vs PE"
))

save_group(
  figure5_grouped,
  file.path(group_dir, "Figure_5_GROUPED_DRAFT_HARRELL"),
  width_mm = 183, height_mm = 177
)
save_group(
  figure6_grouped,
  file.path(group_dir, "Figure_6_GROUPED_DRAFT_HARRELL"),
  width_mm = 183, height_mm = 177
)

# Rebuild unchanged Supplementary Figure 2A-C from their locked Source Data,
# then combine them with the strict fully nested replacement panel D.
supp_source <- file.path(
  root_dir, "SDC1_AAD_SUPPLEMENT_REVIEW_20260722",
  "04_figures", "supplementary", "source_data"
)
detection <- read.csv(file.path(supp_source, "sdc1_detection_by_group_v7.csv"), check.names = FALSE)
draws <- read.csv(file.path(supp_source, "sdc1_multiple_imputation_draws_v7.csv"), check.names = FALSE)
sensitivity <- read.csv(file.path(supp_source, "sdc1_imputation_parameter_sensitivity_v7.csv"), check.names = FALSE)

supp_palette <- c(
  "Healthy control" = "#6F8F88",
  "AMI" = "#D69032",
  "AAD" = "#B6404C",
  "teal" = "#2F8A7C"
)

detection$group <- factor(detection$group, levels = c("Healthy control", "AMI", "AAD"))
p2a <- ggplot(detection, aes(group, detection_fraction, fill = group)) +
  geom_col(width = 0.64, colour = "black", linewidth = 0.28) +
  geom_text(
    aes(label = paste0(detected, "/", total, "\n", scales::percent(detection_fraction, accuracy = 0.1))),
    vjust = -0.28, size = 2.35, lineheight = 0.88
  ) +
  scale_fill_manual(values = supp_palette[c("Healthy control", "AMI", "AAD")], guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.17), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "SDC1 detection", title = "Observed detection pattern") +
  theme_pub(6.5)
p2a <- tag_panel(p2a, "A")

draws$contrast <- factor(
  draws$contrast,
  levels = c("AAD versus AMI", "AAD versus Healthy control"),
  labels = c("AAD vs AMI", "AAD vs healthy control")
)
p2b <- ggplot(draws, aes(contrast, log2FC, fill = contrast)) +
  geom_violin(width = 0.82, alpha = 0.35, colour = NA, trim = TRUE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_jitter(width = 0.08, size = 0.55, alpha = 0.45) +
  scale_fill_manual(
    values = c("AAD vs AMI" = unname(supp_palette["AMI"]), "AAD vs healthy control" = unname(supp_palette["Healthy control"])),
    guide = "none"
  ) +
  labs(
    x = NULL, y = "Adjusted log2 fold change",
    title = "Between-imputation effects",
    subtitle = "50 primary imputations"
  ) +
  theme_pub(6.5) +
  theme(axis.text.x = element_text(angle = 12, hjust = 0.5))
p2b <- tag_panel(p2b, "B")

sensitivity_long <- rbind(
  data.frame(
    shift = sensitivity$shift,
    scale = sensitivity$scale,
    contrast = "AAD vs AMI",
    log2FC = sensitivity$sdc1_log2FC_aad_vs_ami,
    rank = sensitivity$sdc1_joint_rank
  ),
  data.frame(
    shift = sensitivity$shift,
    scale = sensitivity$scale,
    contrast = "AAD vs healthy control",
    log2FC = sensitivity$sdc1_log2FC_aad_vs_healthy,
    rank = sensitivity$sdc1_joint_rank
  )
)
sensitivity_long$shift_label <- factor(
  paste0("shift ", sensitivity_long$shift),
  levels = paste0("shift ", sort(unique(sensitivity_long$shift), decreasing = TRUE))
)
sensitivity_long$scale_label <- factor(
  paste0("scale ", sensitivity_long$scale),
  levels = paste0("scale ", sort(unique(sensitivity_long$scale)))
)
p2c <- ggplot(sensitivity_long, aes(scale_label, shift_label, fill = log2FC)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(aes(label = paste0("rank ", rank, "\n", sprintf("%.2f", log2FC))), size = 2.05, lineheight = 0.88) +
  facet_wrap(~contrast, nrow = 1) +
  scale_fill_gradient(low = "#DDE9E7", high = unname(supp_palette["AAD"]), name = "log2FC") +
  labs(x = NULL, y = NULL, title = "Pre-fixed imputation sensitivity grid") +
  theme_pub(6.3) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
p2c <- tag_panel(p2c, "C")

s2_top <- (p2a | p2b) + plot_layout(widths = c(0.82, 1.18))
s2_bottom <- (p2c | p_s2d) + plot_layout(widths = c(1.18, 0.82))
supplementary_figure2_grouped <- s2_top / s2_bottom +
  plot_layout(heights = c(0.92, 1.08)) &
  theme(plot.margin = margin(4, 5, 4, 5))

save_group(
  supplementary_figure2_grouped,
  file.path(group_dir, "Supplementary_Figure_2_GROUPED_DRAFT_NESTED_SELECTION"),
  width_mm = 183, height_mm = 138
)

supplementary_lm_grouped <- (p_l | p_m) +
  plot_layout(widths = c(1, 1)) &
  theme(plot.margin = margin(4, 5, 4, 5))
save_group(
  supplementary_lm_grouped,
  file.path(group_dir, "Supplementary_Figure_LM_GROUPED_DRAFT_HARRELL"),
  width_mm = 183, height_mm = 112
)

stems <- file.path(
  group_dir,
  c(
    "Figure_5_GROUPED_DRAFT_HARRELL",
    "Figure_6_GROUPED_DRAFT_HARRELL",
    "Supplementary_Figure_2_GROUPED_DRAFT_NESTED_SELECTION",
    "Supplementary_Figure_LM_GROUPED_DRAFT_HARRELL"
  )
)
exports <- unlist(lapply(stems, function(stem) paste0(stem, c(".svg", ".pdf", ".tiff", ".png"))))
qa <- data.frame(
  check = c(
    "four_grouped_drafts", "all_four_export_formats",
    "figure5_has_nine_panels", "figure6_has_nine_panels",
    "supplementary2_uses_strict_panel_D"
  ),
  pass = c(
    length(stems) == 4L,
    all(file.exists(exports)),
    length(figure5$panels) == 9L,
    length(figure6$panels) == 9L,
    sum(selection_draws$sdc1_retained & selection_draws$sdc1_rank == 1, na.rm = TRUE) == 295L
  ),
  stringsAsFactors = FALSE
)
if (calibration_split) {
  qa <- rbind(
    qa,
    data.frame(
      check = c(
        "figure5_D_E_include_addrs",
        "figure6_D_E_include_addrs",
        "paired_improvement_panels_replaced"
      ),
      pass = c(
        all(c("ADD_RS", "ADD_RS_Ddimer", "ADD_RS_SDC1", "Full") %in%
          as.character(figure5$panels$D$data$model)) &&
          all(c("ADD_RS", "ADD_RS_Ddimer", "ADD_RS_SDC1", "Full") %in%
            as.character(figure5$panels$E$data$model)),
        all(c("ADD_RS", "ADD_RS_Ddimer", "ADD_RS_SDC1", "Full") %in%
          as.character(figure6$panels$D$data$model)) &&
          all(c("ADD_RS", "ADD_RS_Ddimer", "ADD_RS_SDC1", "Full") %in%
            as.character(figure6$panels$E$data$model)),
        identical(unique(as.character(figure5$panels$D$data$metric)), "Calibration intercept") &&
          identical(unique(as.character(figure5$panels$E$data$metric)), "Calibration slope") &&
          identical(unique(as.character(figure6$panels$D$data$metric)), "Calibration intercept") &&
          identical(unique(as.character(figure6$panels$E$data$metric)), "Calibration slope")
      ),
      stringsAsFactors = FALSE
    )
  )
}
write.csv(qa, file.path(group_dir, "GROUPED_DRAFTS_QA.csv"), row.names = FALSE)
writeLines(
  c(
    "Grouped draft figures for author adjustment",
    "Figure 5 and Figure 6: 3 x 3 A-I layouts, 183 x 177 mm.",
    if (calibration_split) {
      "Panels D/E: calibration intercept/slope shown separately for all four models, including ADD-RS alone; former paired-improvement panels removed."
    } else {
      "Panels D/E: original combined-calibration and paired-improvement layout."
    },
    "Supplementary Figure 2: locked A-C plus strict nested-selection D, 183 x 138 mm.",
    "Supplementary L/M: strict standalone Harrell panels in a two-column layout, 183 x 112 mm.",
    "All grouped drafts exported as SVG, PDF, 600 dpi TIFF and PNG.",
    paste0("QA: ", if (all(qa$pass)) "PASS" else "FAIL")
  ),
  file.path(group_dir, "README.txt")
)

cat("Grouped drafts generated. QA:", if (all(qa$pass)) "PASS" else "FAIL", "\n")
cat("Output:", group_dir, "\n")
if (!all(qa$pass)) quit(status = 1)

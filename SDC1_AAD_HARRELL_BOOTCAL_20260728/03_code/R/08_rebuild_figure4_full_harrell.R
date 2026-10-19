options(stringsAsFactors = FALSE, scipen = 999)

# Figure contract: see 00_project_docs/FIGURE4_FULL_HARRELL_CONTRACT_20260729.md
# This script intentionally rebuilds all panels in R and never reuses the old
# repeated-CV probability-score intervals in panels G-M.

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(
  file.path(root_dir, "SDC1_AAD_HARRELL_BOOTCAL_20260728"),
  winslash = "/", mustWork = TRUE
)
base_dir <- normalizePath(
  file.path(root_dir, "SDC1_AAD_REANALYSIS_FINAL_20260721"),
  winslash = "/", mustWork = TRUE
)
.libPaths(c(file.path(project_dir, "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(rms)
  library(Hmisc)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

standalone_run <- "FINAL_STANDALONE_B1000_LOWESS_ABCLOC_LOCKED"
output_run <- "FIGURE4_FULL_HARRELL_ABCLOC_LOCKED_20260729"

standalone_rds_file <- file.path(
  project_dir, "04_results", "source_data", standalone_run,
  "standalone_harrell_bootcal_complete_results.rds"
)
clinical_file <- file.path(base_dir, "02_analysis_data", "elisa_evaluation_primary.csv")
base_source_dir <- file.path(base_dir, "04_results", "source_data")
base_stat_dir <- file.path(base_dir, "04_results", "statistical_results")

required_files <- c(
  standalone_rds_file, clinical_file,
  file.path(base_source_dir, "fig4_panel_a_source_data.csv"),
  file.path(base_source_dir, "fig4_panel_b_source_data.csv"),
  file.path(base_source_dir, "fig4_panels_c_to_f_source_data.csv"),
  file.path(base_stat_dir, "fig4_adjusted_total_sdc1_group_means.csv"),
  file.path(base_stat_dir, "fig4_adjusted_d_dimer_group_means.csv"),
  file.path(base_stat_dir, "fig4_adjusted_hypertension_means.csv"),
  file.path(base_stat_dir, "fig4_adjusted_total_sdc1_group_contrasts.csv"),
  file.path(base_stat_dir, "fig4_adjusted_d_dimer_group_contrasts.csv"),
  file.path(base_stat_dir, "fig4_adjusted_hypertension_effects.csv")
)
if (!all(file.exists(required_files))) {
  stop("Missing locked input(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

figure_root <- file.path(project_dir, "04_results", "figures", output_run)
panel_dir <- file.path(figure_root, "panels")
source_dir <- file.path(project_dir, "04_results", "source_data", output_run)
stat_dir <- file.path(project_dir, "04_results", "statistical_results", output_run)
qc_dir <- file.path(project_dir, "07_quality_control", output_run)
log_dir <- file.path(project_dir, "04_results", "logs", output_run)
for (d in c(figure_root, panel_dir, source_dir, stat_dir, qc_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

read_source <- function(name) {
  read.csv(file.path(base_source_dir, name), check.names = FALSE)
}
read_result <- function(name) {
  read.csv(file.path(base_stat_dir, name), check.names = FALSE)
}

clip01 <- function(x) pmin(pmax(x, 0), 1)
percent_label <- function(x) paste0(round(100 * x), "%")
format_probability <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.0001) return("<0.0001")
  if (value < 0.001) return(sprintf("%.4f", value))
  sprintf("%.3f", value)
}

diagnosis_palette <- c(
  "Healthy control" = "#4B8F8A", "AAD" = "#BF3F4C",
  "AMI" = "#3F78A8", "Angina" = "#D78B2B", "PE" = "#7B5AA6"
)
standalone_colours <- c(D_dimer = "#3F78A8", SDC1 = "#D86613")
standalone_labels <- c(D_dimer = "D-dimer", SDC1 = "SDC1")
population_labels <- c(
  symptomatic = "AAD vs symptomatic non-AAD",
  aad_vs_pe = "AAD vs PE"
)

theme_pub <- function(base_size = 6.6) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.36, colour = "#111111"),
      axis.ticks = element_line(linewidth = 0.34, colour = "#111111"),
      axis.title = element_text(size = base_size + 0.1, colour = "#111111"),
      axis.text = element_text(size = base_size - 0.4, colour = "#222222"),
      legend.title = element_text(size = base_size - 0.5, face = "bold"),
      legend.text = element_text(size = base_size - 0.8),
      strip.text = element_text(size = base_size - 0.3, face = "bold"),
      plot.title = element_text(size = base_size + 0.7, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 0.8, colour = "#5F666C", hjust = 0.5),
      plot.tag = element_text(size = 8.5, face = "bold", family = "Arial"),
      plot.tag.position = c(0.01, 1.025),
      panel.grid = element_blank(),
      plot.margin = margin(7.0, 5, 4, 5)
    )
}

tag_panel <- function(plot, tag) {
  plot + labs(tag = tag) +
    theme(plot.tag = element_text(size = 8.5, face = "bold", family = "Arial"))
}

save_plot_bundle <- function(plot, stem, width_mm, height_mm, dpi = 600L) {
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

# -----------------------------------------------------------------------------
# Panels A-F: locked covariate-adjusted descriptive evidence
# -----------------------------------------------------------------------------
a_data <- read_source("fig4_panel_a_source_data.csv")
b_data <- read_source("fig4_panel_b_source_data.csv")
htn_data <- read_source("fig4_panels_c_to_f_source_data.csv")
a_means <- read_result("fig4_adjusted_total_sdc1_group_means.csv")
b_means <- read_result("fig4_adjusted_d_dimer_group_means.csv")
htn_means <- read_result("fig4_adjusted_hypertension_means.csv")
a_contrasts <- read_result("fig4_adjusted_total_sdc1_group_contrasts.csv")
b_contrasts <- read_result("fig4_adjusted_d_dimer_group_contrasts.csv")
htn_effects <- read_result("fig4_adjusted_hypertension_effects.csv")

all_group_levels <- c("Healthy control", "AAD", "AMI", "Angina", "PE")
symptomatic_levels <- c("AAD", "AMI", "Angina", "PE")
a_data$group <- factor(a_data$group, levels = all_group_levels)
b_data$group <- factor(b_data$group, levels = symptomatic_levels)
a_means$group <- factor(a_means$group, levels = all_group_levels)
b_means$group <- factor(b_means$group, levels = symptomatic_levels)
htn_data$group <- factor(htn_data$group, levels = symptomatic_levels)
htn_means$group <- factor(htn_means$group, levels = symptomatic_levels)
htn_data$hypertension_label <- factor(
  htn_data$hypertension_label,
  levels = c("No hypertension", "Hypertension")
)
htn_means$hypertension_label <- factor(
  htn_means$hypertension_label,
  levels = c("No hypertension", "Hypertension")
)

a_labels <- setNames(
  paste0(c("HC", "AAD", "AMI", "Angina", "PE"), "\n(n=",
         a_means$n_group[match(all_group_levels, a_means$group)], ")"),
  all_group_levels
)
a_annotation <- merge(
  data.frame(reference = c("Healthy control", "AMI", "Angina", "PE")),
  a_contrasts[, c("reference", "bh_q")], by = "reference"
)
a_annotation$x <- match(a_annotation$reference, all_group_levels)
a_annotation$y <- max(a_data$total_sdc1) * 1.34
a_annotation$label <- paste0(
  "vs AAD\nq", vapply(a_annotation$bh_q, format_probability, character(1))
)

p_a <- ggplot(a_data, aes(group, total_sdc1, colour = group)) +
  geom_jitter(
    position = position_jitter(width = 0.11, height = 0, seed = 4101),
    size = 0.58, alpha = 0.42
  ) +
  geom_errorbar(
    data = a_means,
    aes(x = group, ymin = lower_display_95, ymax = upper_display_95),
    inherit.aes = FALSE, width = 0.14, linewidth = 0.52, colour = "#202428"
  ) +
  geom_point(
    data = a_means, aes(x = group, y = estimate_display), inherit.aes = FALSE,
    shape = 21, size = 2.05, stroke = 0.58, fill = "white", colour = "#202428"
  ) +
  geom_text(
    data = a_annotation, aes(x = x, y = y, label = label), inherit.aes = FALSE,
    size = 1.8, lineheight = 0.9
  ) +
  scale_colour_manual(values = diagnosis_palette, guide = "none") +
  scale_x_discrete(labels = a_labels) +
  scale_y_log10(expand = expansion(mult = c(0.04, 0.23))) +
  labs(
    x = NULL, y = "Total-plasma SDC1 (pg/mL, log scale)",
    title = "Total-plasma SDC1", subtitle = "Common-covariate adjusted estimates"
  ) + theme_pub(6.5)
p_a <- tag_panel(p_a, "A")
p_a <- p_a + theme(plot.tag.position = c(0.10, 1.025))

b_labels <- setNames(
  paste0(symptomatic_levels, "\n(n=",
         b_means$n_group[match(symptomatic_levels, b_means$group)], ")"),
  symptomatic_levels
)
b_annotation <- merge(
  data.frame(reference = c("AMI", "Angina", "PE")),
  b_contrasts[, c("reference", "bh_q")], by = "reference"
)
b_annotation$x <- match(b_annotation$reference, symptomatic_levels)
b_annotation$y <- max(b_data$log2_d_dimer) + 1.1
b_annotation$label <- paste0(
  "vs AAD\nq", vapply(b_annotation$bh_q, format_probability, character(1))
)

p_b <- ggplot(b_data, aes(group, log2_d_dimer, colour = group)) +
  geom_jitter(
    position = position_jitter(width = 0.11, height = 0, seed = 4102),
    size = 0.60, alpha = 0.44
  ) +
  geom_errorbar(
    data = b_means, aes(x = group, ymin = lower_log2_95, ymax = upper_log2_95),
    inherit.aes = FALSE, width = 0.14, linewidth = 0.52, colour = "#202428"
  ) +
  geom_point(
    data = b_means, aes(x = group, y = estimate_log2), inherit.aes = FALSE,
    shape = 21, size = 2.05, stroke = 0.58, fill = "white", colour = "#202428"
  ) +
  geom_text(
    data = b_annotation, aes(x = x, y = y, label = label), inherit.aes = FALSE,
    size = 1.8, lineheight = 0.9
  ) +
  scale_colour_manual(values = diagnosis_palette, guide = "none") +
  scale_x_discrete(labels = b_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(
    x = NULL, y = "log2 D-dimer (mg/L)", title = "D-dimer",
    subtitle = "Symptomatic-cohort adjusted estimates"
  ) + theme_pub(6.5)
p_b <- tag_panel(p_b, "B")

htn_palettes <- list(
  AAD = c("No hypertension" = "#E8B7BD", "Hypertension" = "#B83C4A"),
  AMI = c("No hypertension" = "#BAD3E7", "Hypertension" = "#397AAE"),
  Angina = c("No hypertension" = "#F0D2A4", "Hypertension" = "#D78B2B"),
  PE = c("No hypertension" = "#D8C9E7", "Hypertension" = "#7B5AA6")
)

add_log_bracket <- function(plot, y, label) {
  plot +
    annotate("segment", x = 1, xend = 2, y = y, yend = y, linewidth = 0.35) +
    annotate("segment", x = 1, xend = 1, y = y, yend = y / 1.08, linewidth = 0.35) +
    annotate("segment", x = 2, xend = 2, y = y, yend = y / 1.08, linewidth = 0.35) +
    annotate("text", x = 1.5, y = y * 1.04, label = label, size = 1.85, vjust = 0)
}

make_htn_panel <- function(group_name, tag, show_y_title = FALSE, jitter_seed = 1L) {
  raw_group <- htn_data[htn_data$group == group_name, , drop = FALSE]
  means_group <- htn_means[htn_means$group == group_name, , drop = FALSE]
  effect <- htn_effects[htn_effects$group == group_name, , drop = FALSE]
  labels <- setNames(
    paste0(c("No HTN", "HTN"), "\n(n=",
           means_group$n_group[match(
             c("No hypertension", "Hypertension"), means_group$hypertension_label
           )], ")"),
    c("No hypertension", "Hypertension")
  )
  bracket_y <- max(raw_group$total_sdc1) * 1.20
  plot <- ggplot(raw_group, aes(hypertension_label, total_sdc1, colour = hypertension_label)) +
    geom_jitter(
      position = position_jitter(width = 0.10, height = 0, seed = jitter_seed),
      size = 0.58, alpha = 0.46
    ) +
    geom_errorbar(
      data = means_group,
      aes(x = hypertension_label, ymin = lower_display_95, ymax = upper_display_95),
      inherit.aes = FALSE, width = 0.14, linewidth = 0.50, colour = "#202428"
    ) +
    geom_point(
      data = means_group, aes(x = hypertension_label, y = estimate_display),
      inherit.aes = FALSE, shape = 21, size = 1.95, stroke = 0.56,
      fill = "white", colour = "#202428"
    ) +
    scale_colour_manual(values = htn_palettes[[group_name]], guide = "none") +
    scale_x_discrete(labels = labels) +
    scale_y_log10(expand = expansion(mult = c(0.04, 0.22))) +
    labs(
      x = NULL,
      y = if (show_y_title) "Total-plasma SDC1 (pg/mL)" else NULL,
      title = group_name
    ) +
    theme_pub(6.1) +
    theme(plot.margin = margin(4.5, 3, 4, 3))
  tag_panel(
    add_log_bracket(plot, bracket_y, paste0("BH q", format_probability(effect$bh_q))),
    tag
  )
}

p_c <- make_htn_panel("AAD", "C", show_y_title = TRUE, jitter_seed = 4201)
p_c <- p_c + theme(plot.tag.position = c(0.12, 1.025))
p_d <- make_htn_panel("AMI", "D", jitter_seed = 4202)
p_e <- make_htn_panel("Angina", "E", jitter_seed = 4203)
p_f <- make_htn_panel("PE", "F", jitter_seed = 4204)

# -----------------------------------------------------------------------------
# Panels G-M: one locked standalone Harrell/ABCLOC evidence chain
# -----------------------------------------------------------------------------
standalone <- readRDS(standalone_rds_file)
clinical <- read.csv(clinical_file, check.names = FALSE)
required_clinical <- c(
  "row_uid", "cluster_id", "group", "aad", "symptomatic",
  "d_dimer_mg_l", "sdc1_source_value"
)
missing_clinical <- setdiff(required_clinical, names(clinical))
if (length(missing_clinical)) {
  stop("Missing clinical columns: ", paste(missing_clinical, collapse = ", "))
}

eligible <- with(
  clinical,
  symptomatic == 1 & !is.na(aad) & !is.na(d_dimer_mg_l) &
    !is.na(sdc1_source_value) & d_dimer_mg_l > 0 & sdc1_source_value > 0
)
d_sym <- clinical[eligible, , drop = FALSE]
d_sym$aad <- as.integer(d_sym$aad)
d_pe <- d_sym[d_sym$group %in% c("AAD", "PE"), , drop = FALSE]
populations <- list(symptomatic = d_sym, aad_vs_pe = d_pe)

stopifnot(
  nrow(d_sym) == 329L, sum(d_sym$aad) == 99L,
  nrow(d_pe) == 142L, sum(d_pe$aad) == 99L,
  length(unique(d_sym$cluster_id)) == 329L
)

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  folds <- integer(length(y))
  for (outcome in c(0L, 1L)) {
    idx <- which(y == outcome)
    idx <- sample(idx, length(idx), replace = FALSE)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

fit_predict_fold <- function(train, test, raw_variable) {
  train_x <- log2(train[[raw_variable]])
  test_x <- log2(test[[raw_variable]])
  train_basis <- Hmisc::rcspline.eval(train_x, nk = 3, inclx = TRUE)
  knots <- attr(train_basis, "knots")
  test_basis <- Hmisc::rcspline.eval(test_x, knots = knots, inclx = TRUE)
  train_model <- data.frame(
    aad = train$aad,
    rcs_1 = train_basis[, 1],
    rcs_2 = train_basis[, 2]
  )
  test_model <- data.frame(
    rcs_1 = test_basis[, 1],
    rcs_2 = test_basis[, 2]
  )
  fit <- suppressWarnings(glm(
    aad ~ rcs_1 + rcs_2, data = train_model, family = binomial(),
    control = glm.control(maxit = 100, epsilon = 1e-9)
  ))
  if (!isTRUE(fit$converged)) stop("OOF fold fit did not converge: ", raw_variable)
  probability <- predict(fit, newdata = test_model, type = "response")
  if (any(!is.finite(probability))) stop("Non-finite OOF prediction: ", raw_variable)
  probability
}

make_repeated_oof <- function(data, seed_base, repeats = 50L, k = 10L) {
  marker_variables <- c(D_dimer = "d_dimer_mg_l", SDC1 = "sdc1_source_value")
  probability_sum <- matrix(
    0, nrow = nrow(data), ncol = length(marker_variables),
    dimnames = list(NULL, names(marker_variables))
  )
  for (repeat_id in seq_len(repeats)) {
    folds <- make_stratified_folds(data$aad, k, seed_base + repeat_id)
    for (fold in seq_len(k)) {
      test_idx <- which(folds == fold)
      train <- data[folds != fold, , drop = FALSE]
      test <- data[test_idx, , drop = FALSE]
      for (model in names(marker_variables)) {
        probability_sum[test_idx, model] <- probability_sum[test_idx, model] +
          fit_predict_fold(train, test, marker_variables[[model]])
      }
    }
  }
  probability_sum <- probability_sum / repeats
  data.frame(
    row_uid = data$row_uid,
    cluster_id = data$cluster_id,
    group = data$group,
    aad = data$aad,
    probability_D_dimer = probability_sum[, "D_dimer"],
    probability_SDC1 = probability_sum[, "SDC1"],
    stringsAsFactors = FALSE
  )
}

auc_rank <- function(y, score) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_step <- function(y, score) {
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  score <- score[ord]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  ends <- c(which(diff(score) != 0), length(score))
  unique(data.frame(
    false_positive_rate = c(0, fp[ends] / sum(y == 0), 1),
    true_positive_rate = c(0, tp[ends] / sum(y == 1), 1)
  ))
}

oof <- list(
  symptomatic = make_repeated_oof(d_sym, seed_base = 202607290L),
  aad_vs_pe = make_repeated_oof(d_pe, seed_base = 202607390L)
)

oof_long <- do.call(rbind, lapply(names(oof), function(population) {
  do.call(rbind, lapply(names(standalone_labels), function(model) {
    data.frame(
      population = population,
      row_uid = oof[[population]]$row_uid,
      cluster_id = oof[[population]]$cluster_id,
      group = oof[[population]]$group,
      aad = oof[[population]]$aad,
      model = model,
      model_label = standalone_labels[[model]],
      repeated_oof_probability = oof[[population]][[paste0("probability_", model)]],
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(oof_long) <- NULL

roc_data <- do.call(rbind, lapply(names(populations), function(population) {
  do.call(rbind, lapply(names(standalone_labels), function(model) {
    curve <- roc_step(
      populations[[population]]$aad,
      oof[[population]][[paste0("probability_", model)]]
    )
    curve$population <- population
    curve$model <- model
    curve$model_label <- standalone_labels[[model]]
    curve
  }))
}))
rownames(roc_data) <- NULL

scalar_all <- rbind(standalone$symptomatic$scalar, standalone$aad_vs_pe$scalar)
calibration_all <- rbind(
  standalone$symptomatic$calibration,
  standalone$aad_vs_pe$calibration
)
difference_all <- rbind(
  standalone$symptomatic$differences,
  standalone$aad_vs_pe$differences
)

oof_auc <- do.call(rbind, lapply(names(populations), function(population) {
  do.call(rbind, lapply(names(standalone_labels), function(model) {
    scalar_row <- scalar_all[
      scalar_all$population == population & scalar_all$model == model &
        scalar_all$metric == "C-index", , drop = FALSE
    ]
    curve_auc <- auc_rank(
      populations[[population]]$aad,
      oof[[population]][[paste0("probability_", model)]]
    )
    data.frame(
      population = population,
      model = model,
      repeated_oof_curve_auc = curve_auc,
      harrell_corrected_c_index = scalar_row$corrected,
      absolute_difference = abs(curve_auc - scalar_row$corrected),
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(oof_auc) <- NULL

make_roc_panel <- function(population, tag) {
  curve <- roc_data[roc_data$population == population, , drop = FALSE]
  curve$model <- factor(curve$model, levels = names(standalone_labels))
  scalar_c <- scalar_all[
    scalar_all$population == population & scalar_all$metric == "C-index", , drop = FALSE
  ]
  scalar_c <- scalar_c[match(names(standalone_labels), scalar_c$model), ]
  legend_labels <- sprintf(
    "%s: %.3f (%.3f-%.3f)",
    standalone_labels[scalar_c$model], scalar_c$corrected, scalar_c$lower, scalar_c$upper
  )
  data <- populations[[population]]
  subtitle <- sprintf(
    "n=%d (%d AAD/%d non-AAD); 50 x 10-fold OOF curve\nSpline knots and model refitted within each fold",
    nrow(data), sum(data$aad), sum(data$aad == 0)
  )
  plot <- ggplot(curve, aes(false_positive_rate, true_positive_rate, colour = model, group = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#929AA1", linewidth = 0.50) +
    geom_step(direction = "hv", linewidth = 0.90) +
    scale_colour_manual(
      values = standalone_colours, breaks = names(standalone_labels),
      labels = legend_labels, name = "Harrell-validated AUROC (ABCLOC 95% CI)"
    ) +
    scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label,
      expand = c(0, 0)
    ) +
    coord_equal() +
    labs(
      x = "False-positive rate", y = "True-positive rate",
      title = paste0("Standalone biomarkers: ", population_labels[[population]]),
      subtitle = subtitle
    ) +
    theme_pub(6.2) +
    theme(
      legend.position = "bottom", legend.direction = "vertical",
      legend.key.width = grid::unit(6.5, "mm"),
      legend.title = element_text(size = 5.1),
      legend.text = element_text(size = 4.9),
      plot.title = element_text(size = 6.9),
      plot.subtitle = element_text(size = 5.0, lineheight = 0.95)
    )
  tag_panel(plot, tag)
}

make_performance_panel <- function(population, tag) {
  plot_data <- scalar_all[
    scalar_all$population == population &
      scalar_all$metric %in% c("Brier score", "Log loss", "C-index"), , drop = FALSE
  ]
  plot_data$model <- factor(plot_data$model, levels = c("SDC1", "D_dimer"))
  plot_data$model_label <- factor(
    standalone_labels[as.character(plot_data$model)],
    levels = c("SDC1", "D-dimer")
  )
  plot_data$metric <- factor(
    plot_data$metric,
    levels = c("Brier score", "Log loss", "C-index")
  )
  plot <- ggplot(plot_data, aes(corrected, model_label, colour = model)) +
    geom_errorbar(
      aes(xmin = lower, xmax = upper), orientation = "y",
      width = 0, linewidth = 0.56
    ) +
    geom_point(size = 2.05) +
    facet_wrap(~metric, scales = "free_x", ncol = 1) +
    scale_colour_manual(values = standalone_colours, guide = "none") +
    labs(
      x = "Optimism-corrected estimate (ABCLOC 95% CI)", y = NULL,
      title = paste0("Validated probability performance: ", population_labels[[population]]),
      subtitle = "1,000 bootstraps; lower Brier/log loss and higher C-index are better"
    ) +
    theme_pub(6.2) +
    theme(
      strip.background = element_rect(fill = "#F1F3F4", colour = NA),
      panel.spacing.y = grid::unit(0.9, "mm"),
      axis.text.y = element_text(size = 5.5),
      axis.text.x = element_text(size = 5.3),
      plot.title = element_text(size = 6.7),
      plot.subtitle = element_text(size = 5.0)
    )
  tag_panel(plot, tag)
}

make_calibration_panel <- function(population, tag) {
  cal <- calibration_all[calibration_all$population == population, , drop = FALSE]
  cal$model <- factor(cal$model, levels = names(standalone_labels))
  cal$corrected_plot <- clip01(cal$corrected_observed)
  cal$apparent_plot <- clip01(cal$apparent_observed)
  cal$lower_plot <- clip01(cal$lower_observed)
  cal$upper_plot <- clip01(cal$upper_observed)
  cal$band_estimable <- is.finite(cal$lower_observed) & is.finite(cal$upper_observed)
  plot <- ggplot(cal, aes(predicted_probability, group = model)) +
    geom_ribbon(
      data = cal[cal$band_estimable, , drop = FALSE],
      aes(ymin = lower_plot, ymax = upper_plot, fill = model),
      alpha = 0.14, colour = NA
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#8A9299", linewidth = 0.50) +
    geom_line(aes(y = apparent_plot, colour = model), linetype = "dotted", linewidth = 0.45, alpha = 0.75) +
    geom_line(aes(y = corrected_plot, colour = model), linewidth = 0.82) +
    scale_colour_manual(
      values = standalone_colours, labels = standalone_labels,
      breaks = names(standalone_labels), name = NULL
    ) +
    scale_fill_manual(values = standalone_colours, guide = "none") +
    scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label,
      expand = c(0, 0)
    ) +
    coord_equal() +
    labs(
      x = "Predicted probability of AAD", y = "Observed AAD proportion",
      title = population_labels[[population]],
      subtitle = "Harrell LOWESS: solid corrected; dotted apparent; dashed ideal"
    ) +
    theme_pub(5.9) +
    theme(
      legend.position = "top", legend.direction = "horizontal",
      legend.key.width = grid::unit(6.5, "mm"),
      plot.title = element_text(size = 6.6),
      plot.subtitle = element_text(size = 4.6)
    )
  tag_panel(plot, tag)
}

make_calibration_parameter_panel <- function(tag) {
  scalar <- scalar_all[
    scalar_all$metric %in% c("Calibration intercept", "Calibration slope"),
    , drop = FALSE
  ]
  scalar$population_label <- factor(
    ifelse(scalar$population == "symptomatic", "Overall", "AAD vs PE"),
    levels = c("Overall", "AAD vs PE")
  )
  scalar$model <- factor(scalar$model, levels = c("SDC1", "D_dimer"))
  scalar$model_label <- factor(
    standalone_labels[as.character(scalar$model)],
    levels = c("SDC1", "D-dimer")
  )
  scalar$metric_label <- factor(
    scalar$metric,
    levels = c("Calibration intercept", "Calibration slope"),
    labels = c("Intercept (target 0)", "Slope (target 1)")
  )
  targets <- expand.grid(
    population_label = levels(scalar$population_label),
    metric_label = levels(scalar$metric_label),
    stringsAsFactors = FALSE
  )
  targets$target <- ifelse(targets$metric_label == "Intercept (target 0)", 0, 1)
  plot <- ggplot(scalar, aes(corrected, model_label, colour = model)) +
    geom_vline(
      data = targets, aes(xintercept = target),
      linetype = "22", colour = "#8A9299", linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    geom_errorbar(
      aes(xmin = lower, xmax = upper), orientation = "y",
      width = 0, linewidth = 0.52
    ) +
    geom_point(size = 1.8) +
    facet_grid(population_label ~ metric_label, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = standalone_colours, guide = "none") +
    labs(
      x = "Optimism-corrected estimate (ABCLOC 95% CI)", y = NULL,
      title = "Calibration parameters"
    ) +
    theme_pub(5.7) +
    theme(
      strip.background = element_rect(fill = "#F1F3F4", colour = NA),
      strip.text = element_text(size = 4.8, face = "bold"),
      axis.text = element_text(size = 4.7),
      axis.title.x = element_text(size = 5.0),
      panel.spacing = grid::unit(0.8, "mm"),
      plot.title = element_text(size = 6.5)
    )
  tag_panel(plot, tag)
}

p_g <- make_roc_panel("symptomatic", "G")
p_h <- make_performance_panel("symptomatic", "H")
p_i <- make_roc_panel("aad_vs_pe", "I")
p_j <- make_performance_panel("aad_vs_pe", "J")
p_k <- make_calibration_panel("symptomatic", "K")
p_l <- make_calibration_panel("aad_vs_pe", "L")
p_m <- make_calibration_parameter_panel("M")

top_row <- (p_a | p_b) + plot_layout(widths = c(1.15, 1.00))
middle_row <- p_c | p_d | p_e | p_f
symptomatic_row <- (p_g | p_h) + plot_layout(widths = c(1.03, 0.97))
pe_row <- (p_i | p_j) + plot_layout(widths = c(1.03, 0.97))
calibration_row <- (p_k | p_l | p_m) + plot_layout(widths = c(0.93, 0.93, 1.14))

figure4 <- top_row / middle_row / symptomatic_row / pe_row / calibration_row +
  plot_layout(heights = c(0.78, 0.70, 0.97, 0.97, 0.88)) &
  theme(plot.margin = margin(3.5, 4, 3.5, 4))

figure_stem <- file.path(figure_root, "Figure_4_FULL_HARRELL_ABCLOC")
save_plot_bundle(figure4, figure_stem, width_mm = 183, height_mm = 260)

panels <- list(
  A = p_a, B = p_b, C = p_c, D = p_d, E = p_e, F = p_f,
  G = p_g, H = p_h, I = p_i, J = p_j, K = p_k, L = p_l, M = p_m
)
panel_sizes <- data.frame(
  panel = names(panels),
  width_mm = c(92, 84, 46, 46, 46, 46, 91, 86, 91, 86, 60, 60, 70),
  height_mm = c(66, 66, 61, 61, 61, 61, 88, 88, 88, 88, 74, 74, 74)
)
for (panel_name in names(panels)) {
  size_row <- panel_sizes[panel_sizes$panel == panel_name, ]
  save_plot_bundle(
    panels[[panel_name]],
    file.path(panel_dir, paste0("Figure_4", panel_name)),
    width_mm = size_row$width_mm,
    height_mm = size_row$height_mm
  )
}

write.csv(
  oof_long,
  file.path(source_dir, "figure4_repeated_oof_patient_predictions.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  roc_data,
  file.path(source_dir, "figure4_roc_curve_source_data.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  calibration_all,
  file.path(source_dir, "figure4_harrell_calibration_source_data.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  scalar_all,
  file.path(stat_dir, "figure4_harrell_scalar_performance.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  difference_all,
  file.path(stat_dir, "figure4_harrell_paired_differences.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  oof_auc,
  file.path(qc_dir, "figure4_oof_curve_vs_harrell_auc_alignment.csv"),
  row.names = FALSE, na = ""
)

panel_manifest <- data.frame(
  panel = LETTERS[1:13],
  content = c(
    "Adjusted total-plasma SDC1 by diagnosis",
    "Adjusted D-dimer by symptomatic diagnosis",
    "Adjusted hypertension contrast within AAD",
    "Adjusted hypertension contrast within AMI",
    "Adjusted hypertension contrast within angina",
    "Adjusted hypertension contrast within PE",
    "Repeated OOF ROC: AAD vs symptomatic non-AAD",
    "Harrell scalar performance: AAD vs symptomatic non-AAD",
    "Repeated OOF ROC: AAD vs PE",
    "Harrell scalar performance: AAD vs PE",
    "Harrell LOWESS calibration: AAD vs symptomatic non-AAD",
    "Harrell LOWESS calibration: AAD vs PE",
    "Harrell calibration intercept and slope"
  ),
  validation_target = c(
    rep("Covariate-adjusted descriptive estimate", 6),
    "50 x 10-fold OOF display; B=1,000 Harrell scalar AUROC",
    "B=1,000 Harrell optimism correction and ABCLOC CI",
    "50 x 10-fold OOF display; B=1,000 Harrell scalar AUROC",
    "B=1,000 Harrell optimism correction and ABCLOC CI",
    rep("B=1,000 Harrell optimism correction and ABCLOC CI", 3)
  ),
  stringsAsFactors = FALSE
)
write.csv(panel_manifest, file.path(figure_root, "PANEL_MANIFEST.csv"), row.names = FALSE)

expected_stems <- c(
  figure_stem,
  file.path(panel_dir, paste0("Figure_4", LETTERS[1:13]))
)
expected_exports <- unlist(lapply(
  expected_stems,
  function(stem) paste0(stem, c(".svg", ".pdf", ".tiff", ".png"))
))

qa <- data.frame(
  check = c(
    "symptomatic_population_329_99",
    "aad_vs_pe_population_142_99",
    "all_oof_predictions_finite",
    "all_scalar_resamples_1000",
    "scalar_intervals_contain_estimates",
    "calibration_corrected_values_present",
    "oof_curve_auc_close_to_harrell_corrected_auc",
    "all_A_to_M_exports_present",
    "panel_manifest_has_13_rows"
  ),
  pass = c(
    nrow(d_sym) == 329L && sum(d_sym$aad) == 99L,
    nrow(d_pe) == 142L && sum(d_pe$aad) == 99L,
    all(is.finite(oof_long$repeated_oof_probability)),
    all(scalar_all$successful_resamples == 1000L),
    all(scalar_all$lower <= scalar_all$corrected & scalar_all$corrected <= scalar_all$upper),
    all(is.finite(calibration_all$corrected_observed)),
    max(oof_auc$absolute_difference) < 0.03,
    all(file.exists(expected_exports)),
    nrow(panel_manifest) == 13L
  ),
  stringsAsFactors = FALSE
)
write.csv(qa, file.path(qc_dir, "FIGURE4_FULL_HARRELL_QA.csv"), row.names = FALSE)

summary_lines <- c(
  "Figure 4 full Harrell/ABCLOC rebuild",
  "A-F: locked covariate-adjusted descriptive panels.",
  "G/I: 50 repeats of stratified 10-fold OOF prediction; log2 transformation, RCS knots, model fit and prediction rebuilt within each fold.",
  "H/J/K/L/M: locked standalone B=1,000 Harrell optimism bootstrap and ABCLOC intervals.",
  "No historical repeated-CV percentile intervals are used in G-M.",
  paste0("Maximum absolute difference between displayed OOF-curve AUROC and Harrell-corrected AUROC: ", sprintf("%.4f", max(oof_auc$absolute_difference))),
  paste0("Computational QA: ", if (all(qa$pass)) "PASS" else "FAIL"),
  paste0("Figure: ", figure_stem)
)
writeLines(summary_lines, file.path(log_dir, "figure4_full_harrell_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
if (!all(qa$pass)) quit(status = 1)

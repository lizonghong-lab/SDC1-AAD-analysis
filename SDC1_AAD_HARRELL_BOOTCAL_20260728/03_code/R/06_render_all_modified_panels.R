options(stringsAsFactors = FALSE, scipen = 999)

root_dir <- normalizePath(getwd(), mustWork = TRUE)
project_dir <- normalizePath(
  file.path(root_dir, "SDC1_AAD_HARRELL_BOOTCAL_20260728"),
  mustWork = TRUE
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

joint_run <- "FINAL_B1000_LOWESS_ABCLOC_LOCKED"
standalone_run <- "FINAL_STANDALONE_B1000_LOWESS_ABCLOC_LOCKED"
output_run <- Sys.getenv(
  "SDC1_OUTPUT_RUN",
  unset = "MODIFIED_PANELS_FOR_MANUAL_ASSEMBLY_20260728"
)
calibration_split <- identical(
  Sys.getenv("SDC1_CALIBRATION_SPLIT", unset = "0"),
  "1"
)

joint_rds_file <- file.path(
  project_dir, "04_results", "source_data", joint_run,
  "harrell_bootcal_complete_results.rds"
)
standalone_rds_file <- file.path(
  project_dir, "04_results", "source_data", standalone_run,
  "standalone_harrell_bootcal_complete_results.rds"
)
clinical_file <- file.path(
  root_dir, "SDC1_AAD_REANALYSIS_FINAL_20260721",
  "02_analysis_data", "elisa_evaluation_primary.csv"
)
selection_draw_file <- file.path(
  project_dir, "04_results", "source_data", "discovery_nested_selection",
  "FINAL_B300_M50_NESTED_LOCKED", "nested_mi_selection_bootstrap_draws.csv"
)
selection_summary_file <- file.path(
  project_dir, "04_results", "statistical_results", "discovery_nested_selection",
  "FINAL_B300_M50_NESTED_LOCKED", "nested_mi_sdc1_selection_stability.csv"
)

required_files <- c(
  joint_rds_file, standalone_rds_file, clinical_file,
  selection_draw_file, selection_summary_file
)
if (!all(file.exists(required_files))) {
  stop("Missing required locked inputs: ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

figure_root <- file.path(project_dir, "04_results", "figures", output_run)
figure5_dir <- file.path(figure_root, "Figure_5")
figure6_dir <- file.path(figure_root, "Figure_6")
supp_dir <- file.path(figure_root, "Supplementary")
source_dir <- file.path(project_dir, "04_results", "source_data", output_run)
stat_dir <- file.path(project_dir, "04_results", "statistical_results", output_run)
qc_dir <- file.path(project_dir, "07_quality_control", output_run)
log_dir <- file.path(project_dir, "04_results", "logs", output_run)
for (d in c(figure_root, figure5_dir, figure6_dir, supp_dir, source_dir, stat_dir, qc_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

joint <- readRDS(joint_rds_file)
standalone <- readRDS(standalone_rds_file)
raw <- read.csv(clinical_file, check.names = FALSE)
selection_draws <- read.csv(selection_draw_file, check.names = FALSE)
selection_summary <- read.csv(selection_summary_file, check.names = FALSE)

eligible <- with(
  raw,
  symptomatic == 1 & !is.na(aad) & !is.na(add_rs_corrected) &
    !is.na(log2_d_dimer) & !is.na(log2_total_sdc1)
)
d_sym <- raw[eligible, , drop = FALSE]
d_sym$aad <- as.integer(d_sym$aad)
d_pe <- d_sym[d_sym$group %in% c("AAD", "PE"), , drop = FALSE]
populations <- list(symptomatic = d_sym, aad_vs_pe = d_pe)
population_results <- list(symptomatic = joint$symptomatic, aad_vs_pe = joint$aad_vs_pe)

stopifnot(
  nrow(d_sym) == 329L, sum(d_sym$aad) == 99L,
  nrow(d_pe) == 142L, sum(d_pe$aad) == 99L,
  nrow(selection_draws) == 300L
)

model_order <- c("ADD_RS", "ADD_RS_Ddimer", "ADD_RS_SDC1", "Full")
biomarker_models <- c("ADD_RS_Ddimer", "ADD_RS_SDC1", "Full")
model_labels <- c(
  ADD_RS = "ADD-RS",
  ADD_RS_Ddimer = "ADD-RS + D-dimer",
  ADD_RS_SDC1 = "ADD-RS + SDC1",
  Full = "ADD-RS + D-dimer + SDC1"
)
model_short <- c(
  ADD_RS = "ADD-RS",
  ADD_RS_Ddimer = "+ D-dimer",
  ADD_RS_SDC1 = "+ SDC1",
  Full = "+ D-dimer + SDC1"
)
model_colours <- c(
  ADD_RS = "#9AA0A6",
  ADD_RS_Ddimer = "#3F78A8",
  ADD_RS_SDC1 = "#D86613",
  Full = "#2A9D8F"
)
outcome_colours <- c(
  AAD = "#BE4652",
  `Symptomatic non-AAD` = "#3F78A8",
  PE = "#3F78A8"
)
population_labels <- c(
  symptomatic = "AAD vs symptomatic non-AAD",
  aad_vs_pe = "AAD vs PE"
)
panel_prefix <- c(symptomatic = "Figure_5", aad_vs_pe = "Figure_6")
panel_directories <- list(symptomatic = figure5_dir, aad_vs_pe = figure6_dir)

theme_pub <- function(base_size = 6.6) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.38, colour = "#111111"),
      axis.ticks = element_line(linewidth = 0.36, colour = "#111111"),
      axis.title = element_text(size = base_size + 0.2, colour = "#111111"),
      axis.text = element_text(size = base_size - 0.3, colour = "#222222"),
      legend.title = element_text(size = base_size - 0.4, face = "bold"),
      legend.text = element_text(size = base_size - 0.7),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 1.0, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 0.6, colour = "#626970", hjust = 0.5),
      plot.tag = element_text(size = 9.0, face = "bold"),
      plot.tag.position = c(0.01, 0.995),
      panel.grid = element_blank(),
      plot.margin = margin(6, 6, 5, 6)
    )
}

percent_label <- function(x) paste0(round(100 * x), "%")
clip01 <- function(x) pmin(pmax(x, 0), 1)

tag_panel <- function(plot, tag) {
  plot + labs(tag = tag) + theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"))
}

save_panel <- function(plot, stem, width_mm = 64, height_mm = 61, dpi = 600L) {
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

add_rcs_basis <- function(data, fit) {
  d_knots <- unname(fit$Design$parms$log2_d_dimer)
  s_knots <- unname(fit$Design$parms$log2_total_sdc1)
  d_basis <- Hmisc::rcspline.eval(data$log2_d_dimer, knots = d_knots, inclx = TRUE)
  s_basis <- Hmisc::rcspline.eval(data$log2_total_sdc1, knots = s_knots, inclx = TRUE)
  data$d_rcs_1 <- d_basis[, 1]
  data$d_rcs_2 <- d_basis[, 2]
  data$s_rcs_1 <- s_basis[, 1]
  data$s_rcs_2 <- s_basis[, 2]
  data
}

oof_formulas <- list(
  ADD_RS = aad ~ add_rs_corrected,
  ADD_RS_Ddimer = aad ~ add_rs_corrected + d_rcs_1 + d_rcs_2,
  ADD_RS_SDC1 = aad ~ add_rs_corrected + s_rcs_1 + s_rcs_2,
  Full = aad ~ add_rs_corrected + d_rcs_1 + d_rcs_2 + s_rcs_1 + s_rcs_2
)

make_repeated_oof <- function(data, population_result, seed_base, repeats = 50L, k = 10L) {
  d <- add_rcs_basis(data, population_result$fits$Full)
  prediction_sum <- matrix(0, nrow = nrow(d), ncol = length(oof_formulas))
  colnames(prediction_sum) <- names(oof_formulas)
  for (repeat_id in seq_len(repeats)) {
    folds <- make_stratified_folds(d$aad, k = k, seed = seed_base + repeat_id)
    for (fold in seq_len(k)) {
      train <- d[folds != fold, , drop = FALSE]
      test_idx <- which(folds == fold)
      test <- d[test_idx, , drop = FALSE]
      for (model in names(oof_formulas)) {
        fit <- suppressWarnings(glm(
          oof_formulas[[model]], data = train, family = binomial(),
          control = glm.control(maxit = 100, epsilon = 1e-9)
        ))
        if (!isTRUE(fit$converged)) stop("OOF fit failed: ", model)
        prediction_sum[test_idx, model] <- prediction_sum[test_idx, model] +
          predict(fit, newdata = test, type = "response")
      }
    }
  }
  prediction_sum <- prediction_sum / repeats
  out <- data.frame(
    row_uid = d$row_uid,
    group = d$group,
    aad = d$aad,
    stringsAsFactors = FALSE
  )
  for (model in colnames(prediction_sum)) out[[paste0("prob_", model)]] <- prediction_sum[, model]
  out
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
  out <- data.frame(
    false_positive_rate = c(0, fp[ends] / sum(y == 0)),
    true_positive_rate = c(0, tp[ends] / sum(y == 1))
  )
  out <- rbind(out, data.frame(false_positive_rate = 1, true_positive_rate = 1))
  unique(out)
}

raw_metric <- function(raw_object, metric, component) {
  raw_columns <- names(raw_object$orig)
  if (metric == "C-index") {
    return((raw_object[[component]][, match("Dxy", raw_columns)] + 1) / 2)
  }
  column <- switch(metric, "Brier score" = "Brier", "Log loss" = "LogLoss")
  raw_object[[component]][, match(column, raw_columns)]
}

orig_metric <- function(raw_object, metric) {
  if (metric == "C-index") return((unname(raw_object$orig["Dxy"]) + 1) / 2)
  column <- switch(metric, "Brier score" = "Brier", "Log loss" = "LogLoss")
  unname(raw_object$orig[column])
}

abcloc_summary <- function(orig, btrain, btest, conf = 0.95) {
  ok <- is.finite(btrain) & is.finite(btest)
  btrain <- btrain[ok]
  btest <- btest[ok]
  optimism_draw <- btrain - btest
  corrected <- orig - mean(optimism_draw)
  directional_sd <- Hmisc::dualSD(btrain - 1.25 * btest)
  z <- qnorm(1 - (1 - conf) / 2)
  c(
    corrected = corrected,
    lower = corrected - z * unname(directional_sd["top"]),
    upper = corrected + z * unname(directional_sd["bottom"]),
    successful_resamples = length(btrain)
  )
}

improvement_vs_addrs <- function(result, population) {
  augmented <- c("ADD_RS_Ddimer", "ADD_RS_SDC1", "Full")
  rows <- list()
  row_id <- 0L
  for (model in augmented) {
    for (metric in c("Brier score", "Log loss", "C-index")) {
      row_id <- row_id + 1L
      a <- result$validations[[model]]$raw
      b <- result$validations$ADD_RS$raw
      if (metric == "C-index") {
        ans <- abcloc_summary(
          orig_metric(a, metric) - orig_metric(b, metric),
          raw_metric(a, metric, "btrain") - raw_metric(b, metric, "btrain"),
          raw_metric(a, metric, "btest") - raw_metric(b, metric, "btest")
        )
      } else {
        ans <- abcloc_summary(
          orig_metric(b, metric) - orig_metric(a, metric),
          raw_metric(b, metric, "btrain") - raw_metric(a, metric, "btrain"),
          raw_metric(b, metric, "btest") - raw_metric(a, metric, "btest")
        )
      }
      rows[[row_id]] <- data.frame(
        population = population,
        model = model,
        model_label = unname(model_short[model]),
        metric = metric,
        corrected_improvement = ans["corrected"],
        lower = ans["lower"],
        upper = ans["upper"],
        successful_resamples = ans["successful_resamples"],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  do.call(rbind, rows)
}

direction_summary <- function(update_data, seed, B = 2000L) {
  update_data$correct <- ifelse(
    update_data$aad == 1,
    update_data$probability_change > 0,
    update_data$probability_change < 0
  )
  labels <- unique(update_data$outcome)
  point <- aggregate(correct ~ outcome, data = update_data, FUN = mean)
  draw_rows <- vector("list", B)
  set.seed(seed)
  event_idx <- which(update_data$aad == 1)
  control_idx <- which(update_data$aad == 0)
  for (b in seq_len(B)) {
    idx <- c(
      sample(event_idx, length(event_idx), replace = TRUE),
      sample(control_idx, length(control_idx), replace = TRUE)
    )
    x <- update_data[idx, , drop = FALSE]
    z <- aggregate(correct ~ outcome, data = x, FUN = mean)
    z$bootstrap <- b
    draw_rows[[b]] <- z
  }
  draws <- do.call(rbind, draw_rows)
  interval <- do.call(rbind, lapply(labels, function(label) {
    x <- draws$correct[draws$outcome == label]
    data.frame(
      outcome = label,
      lower = unname(quantile(x, 0.025)),
      upper = unname(quantile(x, 0.975)),
      stringsAsFactors = FALSE
    )
  }))
  names(point)[2] <- "estimate"
  list(summary = merge(point, interval, by = "outcome", sort = FALSE), draws = draws)
}

make_panel_set <- function(population, data, result, oof) {
  prefix <- unname(panel_prefix[population])
  out_dir <- panel_directories[[population]]
  comparator_label <- if (population == "symptomatic") "Symptomatic non-AAD" else "PE"

  calibration <- result$calibration[result$calibration$model %in% biomarker_models, , drop = FALSE]
  calibration$model <- factor(calibration$model, levels = biomarker_models)
  calibration$model_label <- factor(model_labels[as.character(calibration$model)], levels = model_labels[biomarker_models])
  calibration$corrected_plot <- clip01(calibration$corrected_observed)
  calibration$lower_plot <- clip01(calibration$lower_observed)
  calibration$upper_plot <- clip01(calibration$upper_observed)
  calibration$band_estimable <- is.finite(calibration$lower_observed) & is.finite(calibration$upper_observed)

  p_a <- ggplot(calibration, aes(predicted_probability, corrected_plot, colour = model, group = model)) +
    geom_ribbon(
      data = calibration[calibration$band_estimable, , drop = FALSE],
      aes(ymin = lower_plot, ymax = upper_plot, fill = model),
      alpha = 0.14, colour = NA, inherit.aes = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#8A9299", linewidth = 0.55) +
    geom_line(linewidth = 0.85) +
    scale_colour_manual(values = model_colours, labels = model_labels[biomarker_models], breaks = biomarker_models) +
    scale_fill_manual(values = model_colours, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    coord_equal() +
    labs(
      x = "Predicted probability of AAD",
      y = "Observed AAD proportion",
      colour = NULL,
      title = "Smooth optimism-corrected calibration",
      subtitle = "LOWESS; ABCLOC 95% CI; 1,000 bootstraps"
    ) +
    theme_pub(6.2) +
    theme(
      legend.position = "top", legend.direction = "horizontal",
      legend.key.width = grid::unit(6.5, "mm"),
      legend.margin = margin(0, 0, 1, 0),
      plot.title = element_text(size = 7.0)
    ) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
  p_a <- tag_panel(p_a, "A")

  roc_rows <- list()
  for (model in biomarker_models) {
    curve <- roc_step(data$aad, oof[[paste0("prob_", model)]])
    curve$model <- model
    roc_rows[[model]] <- curve
  }
  roc_data <- do.call(rbind, roc_rows)
  scalar_c <- result$scalar[
    result$scalar$metric == "C-index" & result$scalar$model %in% biomarker_models,
    , drop = FALSE
  ]
  scalar_c <- scalar_c[match(biomarker_models, scalar_c$model), ]
  legend_labels <- sprintf(
    "%s: %.3f (%.3f-%.3f)",
    model_labels[scalar_c$model], scalar_c$corrected, scalar_c$lower, scalar_c$upper
  )
  p_b <- ggplot(roc_data, aes(false_positive_rate, true_positive_rate, colour = model, group = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#9AA1A7", linewidth = 0.55) +
    geom_step(linewidth = 0.92, direction = "hv") +
    scale_colour_manual(
      values = model_colours, breaks = biomarker_models,
      labels = legend_labels, name = "Validated AUROC (ABCLOC 95% CI)"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    coord_equal() +
    labs(
      x = "False-positive rate", y = "True-positive rate",
      title = paste0("Joint models: ", unname(population_labels[population]))
    ) +
    theme_pub(6.2) +
    theme(
      legend.position = "inside", legend.position.inside = c(0.98, 0.03),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.88), colour = NA),
      legend.key.width = grid::unit(5.5, "mm"),
      legend.title = element_text(size = 5.2),
      legend.text = element_text(size = 4.9),
      plot.title = element_text(size = 6.8)
    )
  p_b <- tag_panel(p_b, "B")

  performance <- result$scalar[result$scalar$metric %in% c("Brier score", "Log loss", "C-index"), , drop = FALSE]
  performance$model <- factor(performance$model, levels = rev(model_order))
  performance$model_short <- factor(model_short[as.character(performance$model)], levels = model_short[rev(model_order)])
  performance$metric <- factor(performance$metric, levels = c("Brier score", "Log loss", "C-index"))
  p_c <- ggplot(performance, aes(corrected, model_short, colour = model)) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0, linewidth = 0.55) +
    geom_point(size = 2.0) +
    facet_wrap(~metric, scales = "free_x", ncol = 1) +
    scale_colour_manual(values = model_colours, guide = "none") +
    labs(x = "Optimism-corrected estimate (95% CI)", y = NULL) +
    theme_pub(6.2) +
    theme(
      strip.background = element_rect(fill = "#F1F3F4", colour = NA),
      panel.spacing.y = grid::unit(1.0, "mm"),
      axis.text.y = element_text(size = 5.4),
      axis.text.x = element_text(size = 5.4)
    )
  p_c <- tag_panel(p_c, "C")

  cal_scalar <- result$scalar[
    result$scalar$metric %in% c("Calibration intercept", "Calibration slope") &
      result$scalar$model %in% biomarker_models,
    , drop = FALSE
  ]
  cal_scalar$model <- factor(cal_scalar$model, levels = rev(biomarker_models))
  cal_scalar$model_short <- factor(model_short[as.character(cal_scalar$model)], levels = model_short[rev(biomarker_models)])
  cal_scalar$metric_label <- factor(
    cal_scalar$metric,
    levels = c("Calibration intercept", "Calibration slope"),
    labels = c("Intercept (target 0)", "Slope (target 1)")
  )
  targets <- data.frame(
    metric_label = factor(c("Intercept (target 0)", "Slope (target 1)"), levels = levels(cal_scalar$metric_label)),
    target = c(0, 1)
  )
  p_d <- ggplot(cal_scalar, aes(corrected, model_short, colour = model)) +
    geom_vline(data = targets, aes(xintercept = target), linetype = "22", colour = "#8A9299", linewidth = 0.5) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0, linewidth = 0.58) +
    geom_point(size = 2.0) +
    facet_wrap(~metric_label, scales = "free_x", ncol = 1) +
    scale_colour_manual(values = model_colours, guide = "none") +
    labs(x = "Optimism-corrected estimate (95% CI)", y = NULL, title = "Calibration intercept and slope") +
    theme_pub(6.2) +
    theme(
      strip.background = element_rect(fill = "#F1F3F4", colour = NA),
      panel.spacing.y = grid::unit(1.2, "mm"),
      axis.text.y = element_text(size = 5.5),
      plot.title = element_text(size = 7.0)
    )
  p_d <- tag_panel(p_d, "D")

  improvements <- improvement_vs_addrs(result, population)
  improvements$model <- factor(improvements$model, levels = rev(biomarker_models))
  improvements$model_label <- factor(model_short[as.character(improvements$model)], levels = model_short[rev(biomarker_models)])
  improvements$metric <- factor(improvements$metric, levels = c("Brier score", "Log loss", "C-index"))
  p_e <- ggplot(improvements, aes(corrected_improvement, model_label, colour = model)) +
    geom_vline(xintercept = 0, linetype = "22", colour = "#8A9299", linewidth = 0.5) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0, linewidth = 0.58) +
    geom_point(size = 2.0) +
    facet_wrap(~metric, scales = "free_x", ncol = 1) +
    scale_colour_manual(values = model_colours, guide = "none") +
    labs(
      x = "Improvement vs ADD-RS (ABCLOC 95% CI)", y = NULL,
      title = unname(population_labels[population])
    ) +
    theme_pub(6.1) +
    theme(
      strip.background = element_rect(fill = "#F1F3F4", colour = NA),
      panel.spacing.y = grid::unit(1.0, "mm"),
      axis.text.y = element_text(size = 5.4),
      axis.text.x = element_text(size = 5.2),
      plot.title = element_text(size = 7.0)
    )
  p_e <- tag_panel(p_e, "E")

  # Optional final-layout revision: replace the paired-improvement panel with
  # a standalone calibration-slope panel, split calibration intercept/slope
  # across D/E, and include the ADD-RS-only model in both panels.
  if (calibration_split) {
    calibration_all <- result$scalar[
      result$scalar$metric %in% c("Calibration intercept", "Calibration slope") &
        result$scalar$model %in% model_order,
      , drop = FALSE
    ]
    calibration_all$model <- factor(calibration_all$model, levels = rev(model_order))
    calibration_all$model_short <- factor(
      model_short[as.character(calibration_all$model)],
      levels = model_short[rev(model_order)]
    )

    calibration_intercept <- calibration_all[
      calibration_all$metric == "Calibration intercept",
      , drop = FALSE
    ]
    p_d <- ggplot(calibration_intercept, aes(corrected, model_short, colour = model)) +
      geom_vline(xintercept = 0, linetype = "22", colour = "#8A9299", linewidth = 0.5) +
      geom_errorbar(
        aes(xmin = lower, xmax = upper),
        orientation = "y", width = 0, linewidth = 0.58
      ) +
      geom_point(size = 2.0) +
      scale_colour_manual(values = model_colours, guide = "none") +
      labs(
        x = "Optimism-corrected estimate (95% CI)", y = NULL,
        title = "Calibration intercept (target 0)"
      ) +
      theme_pub(6.2) +
      theme(
        axis.text.y = element_text(size = 5.5),
        axis.text.x = element_text(size = 5.3),
        plot.title = element_text(size = 7.0)
      )
    p_d <- tag_panel(p_d, "D")

    calibration_slope <- calibration_all[
      calibration_all$metric == "Calibration slope",
      , drop = FALSE
    ]
    p_e <- ggplot(calibration_slope, aes(corrected, model_short, colour = model)) +
      geom_vline(xintercept = 1, linetype = "22", colour = "#8A9299", linewidth = 0.5) +
      geom_errorbar(
        aes(xmin = lower, xmax = upper),
        orientation = "y", width = 0, linewidth = 0.58
      ) +
      geom_point(size = 2.0) +
      scale_colour_manual(values = model_colours, guide = "none") +
      labs(
        x = "Optimism-corrected estimate (95% CI)", y = NULL,
        title = "Calibration slope (target 1)"
      ) +
      theme_pub(6.2) +
      theme(
        axis.text.y = element_text(size = 5.5),
        axis.text.x = element_text(size = 5.3),
        plot.title = element_text(size = 7.0)
      )
    p_e <- tag_panel(p_e, "E")
  }

  full_lr <- unname(result$fits$Full$stats["Model L.R."])
  benchmark_lr <- unname(result$fits$ADD_RS_Ddimer$stats["Model L.R."])
  lr_information <- data.frame(
    component = c("ADD-RS + D-dimer", "Additional SDC1 information"),
    lr_chisq = c(benchmark_lr, full_lr - benchmark_lr),
    fraction = c(benchmark_lr / full_lr, (full_lr - benchmark_lr) / full_lr),
    stringsAsFactors = FALSE
  )
  lr_information$component <- factor(lr_information$component, levels = rev(lr_information$component))
  p_f <- ggplot(lr_information, aes(fraction, component, fill = component)) +
    geom_col(width = 0.58) +
    geom_text(aes(label = paste0(round(100 * fraction, 1), "%")), hjust = 1.06, colour = "white", fontface = "bold", size = 2.35) +
    scale_fill_manual(values = c("ADD-RS + D-dimer" = "#A8AFB5", "Additional SDC1 information" = "#BE4652"), guide = "none") +
    scale_x_continuous(labels = percent_label, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(x = "Fraction of full-model LR information", y = NULL, title = "Relative new information") +
    theme_pub(6.2) +
    theme(axis.text.y = element_text(size = 5.7), plot.title = element_text(size = 7.2))
  p_f <- tag_panel(p_f, "F")

  update <- data.frame(
    row_uid = oof$row_uid,
    group = oof$group,
    aad = oof$aad,
    outcome = ifelse(oof$aad == 1, "AAD", comparator_label),
    benchmark_probability = oof$prob_ADD_RS_Ddimer,
    full_probability = oof$prob_Full,
    stringsAsFactors = FALSE
  )
  update$probability_change <- update$full_probability - update$benchmark_probability
  p_g <- ggplot(update, aes(benchmark_probability, full_probability, colour = outcome)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#8A9299", linewidth = 0.5) +
    geom_point(size = 1.1, alpha = 0.63) +
    scale_colour_manual(values = outcome_colours[c("AAD", comparator_label)], breaks = c("AAD", comparator_label)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    coord_equal() +
    labs(
      x = "ADD-RS + D-dimer probability", y = "Full-model probability",
      colour = NULL, title = "Patient probability updating"
    ) +
    theme_pub(6.2) +
    theme(
      legend.position = "top", legend.direction = "horizontal",
      legend.key.width = grid::unit(4.0, "mm"),
      plot.title = element_text(size = 7.0)
    )
  p_g <- tag_panel(p_g, "G")

  p_h <- ggplot(update, aes(probability_change, colour = outcome, fill = outcome)) +
    geom_vline(xintercept = 0, linetype = "22", colour = "#8A9299", linewidth = 0.5) +
    geom_density(alpha = 0.18, linewidth = 0.75, adjust = 0.9) +
    scale_colour_manual(values = outcome_colours[c("AAD", comparator_label)]) +
    scale_fill_manual(values = outcome_colours[c("AAD", comparator_label)]) +
    scale_x_continuous(labels = percent_label) +
    labs(x = "Probability change after adding SDC1", y = "Density") +
    theme_pub(6.2) +
    theme(legend.position = "none")
  p_h <- tag_panel(p_h, "H")

  direction <- direction_summary(update, seed = if (population == "symptomatic") 202607289L else 202607290L)
  direction$summary$outcome <- factor(direction$summary$outcome, levels = c("AAD", comparator_label))
  p_i <- ggplot(direction$summary, aes(outcome, estimate, fill = outcome)) +
    geom_col(width = 0.58) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.16, linewidth = 0.55) +
    geom_text(aes(y = pmin(upper + 0.045, 1.045), label = percent_label(estimate)), size = 2.45) +
    scale_fill_manual(values = outcome_colours[c("AAD", comparator_label)], guide = "none") +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    labs(x = NULL, y = "Correct-direction updates", title = "Directionally coherent updating") +
    theme_pub(6.2) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1), plot.title = element_text(size = 7.0))
  p_i <- tag_panel(p_i, "I")

  panels <- list(A = p_a, B = p_b, C = p_c, D = p_d, E = p_e, F = p_f, G = p_g, H = p_h, I = p_i)
  panel_heights <- c(A = 64, B = 64, C = 63, D = 62, E = 63, F = 58, G = 61, H = 58, I = 58)
  for (letter in names(panels)) {
    save_panel(
      panels[[letter]],
      file.path(out_dir, paste0(prefix, letter)),
      width_mm = 64,
      height_mm = unname(panel_heights[letter])
    )
  }

  oof_long <- do.call(rbind, lapply(model_order, function(model) {
    data.frame(
      population = population,
      row_uid = oof$row_uid,
      group = oof$group,
      aad = oof$aad,
      model = model,
      repeated_oof_probability = oof[[paste0("prob_", model)]],
      stringsAsFactors = FALSE
    )
  }))
  roc_data$population <- population
  update$population <- population
  direction$summary$population <- population
  direction$draws$population <- population
  improvements$population <- population
  lr_information$population <- population

  list(
    panels = panels,
    calibration = calibration,
    performance = result$scalar,
    improvements = improvements,
    lr_information = lr_information,
    oof = oof_long,
    roc = roc_data,
    update = update,
    direction = direction$summary,
    direction_draws = direction$draws,
    oof_auc = data.frame(
      population = population,
      model = biomarker_models,
      oof_auc = vapply(biomarker_models, function(model) auc_rank(data$aad, oof[[paste0("prob_", model)]]), numeric(1)),
      corrected_auc = scalar_c$corrected,
      lower = scalar_c$lower,
      upper = scalar_c$upper,
      stringsAsFactors = FALSE
    )
  )
}

cat("Generating repeated OOF predictions with locked RCS bases...\n")
oof_sym <- make_repeated_oof(d_sym, joint$symptomatic, seed_base = 202607300L)
oof_pe <- make_repeated_oof(d_pe, joint$aad_vs_pe, seed_base = 202607400L)

figure5 <- make_panel_set("symptomatic", d_sym, joint$symptomatic, oof_sym)
figure6 <- make_panel_set("aad_vs_pe", d_pe, joint$aad_vs_pe, oof_pe)

all_oof <- rbind(figure5$oof, figure6$oof)
all_roc <- rbind(figure5$roc, figure6$roc)
all_update <- rbind(figure5$update, figure6$update)
all_direction <- rbind(figure5$direction, figure6$direction)
all_direction_draws <- rbind(figure5$direction_draws, figure6$direction_draws)
all_improvements <- rbind(figure5$improvements, figure6$improvements)
all_lr <- rbind(figure5$lr_information, figure6$lr_information)
all_oof_auc <- rbind(figure5$oof_auc, figure6$oof_auc)
all_scalar <- rbind(joint$symptomatic$scalar, joint$aad_vs_pe$scalar)
all_calibration <- rbind(joint$symptomatic$calibration, joint$aad_vs_pe$calibration)

write.csv(all_oof, file.path(source_dir, "figure5_6_repeated_oof_patient_predictions.csv"), row.names = FALSE, na = "")
write.csv(all_roc, file.path(source_dir, "figure5_6_roc_curve_source_data.csv"), row.names = FALSE, na = "")
write.csv(all_update, file.path(source_dir, "figure5_6_probability_update_source_data.csv"), row.names = FALSE, na = "")
write.csv(all_direction_draws, file.path(source_dir, "figure5_6_direction_participant_bootstrap_draws.csv"), row.names = FALSE, na = "")
write.csv(all_calibration, file.path(source_dir, "figure5_6_harrell_calibration_source_data.csv"), row.names = FALSE, na = "")
write.csv(all_scalar, file.path(stat_dir, "figure5_6_harrell_scalar_performance.csv"), row.names = FALSE, na = "")
write.csv(all_improvements, file.path(stat_dir, "figure5_6_improvement_vs_addrs.csv"), row.names = FALSE, na = "")
write.csv(all_lr, file.path(stat_dir, "figure5_6_likelihood_ratio_information.csv"), row.names = FALSE, na = "")
write.csv(all_direction, file.path(stat_dir, "figure5_6_directionally_coherent_updates.csv"), row.names = FALSE, na = "")
write.csv(all_oof_auc, file.path(qc_dir, "figure5_6_oof_curve_auc_alignment.csv"), row.names = FALSE, na = "")

# Supplementary Figure 2D: strict nested selection status.
selection_draws$status <- ifelse(
  selection_draws$sdc1_retained & selection_draws$sdc1_rank == 1,
  "Retained; rank 1",
  ifelse(selection_draws$sdc1_retained, "Retained; rank >1", "Below detection filter")
)
selection_counts <- as.data.frame(table(selection_draws$status), stringsAsFactors = FALSE)
names(selection_counts) <- c("status", "count")
selection_counts$frequency <- selection_counts$count / nrow(selection_draws)
selection_counts$status <- factor(
  selection_counts$status,
  levels = c("Retained; rank 1", "Retained; rank >1", "Below detection filter")
)
selection_counts <- selection_counts[order(selection_counts$status), ]
p_s2d <- ggplot(selection_counts, aes(status, count, fill = status)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.32) +
  geom_text(
    aes(y = count + 18, label = paste0(count, "/300\n", scales::percent(frequency, accuracy = 0.1))),
    vjust = 0.5, size = 2.55, lineheight = 0.90
  ) +
  scale_fill_manual(
    values = c(
      "Retained; rank 1" = "#2F8A7C",
      "Retained; rank >1" = "#D69032",
      "Below detection filter" = "#C9CED3"
    ),
    guide = "none"
  ) +
  scale_x_discrete(labels = c(
    "Retained; rank 1" = "Retained;\nrank 1",
    "Retained; rank >1" = "Retained;\nrank >1",
    "Below detection filter" = "Below detection\nfilter"
  )) +
  scale_y_continuous(limits = c(0, 335), breaks = c(0, 100, 200, 300), expand = c(0, 0)) +
  labs(
    x = NULL, y = "Bootstrap resamples",
    title = "Fully nested selection stability",
    subtitle = "300 resamples; 50 fresh imputations per resample"
  ) +
  theme_pub(6.6) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title = element_text(size = 7.4))
p_s2d <- tag_panel(p_s2d, "D")
save_panel(p_s2d, file.path(supp_dir, "Supplementary_Figure_2D"), width_mm = 82, height_mm = 64)
write.csv(selection_counts, file.path(source_dir, "supplementary_figure_2D_nested_selection_counts.csv"), row.names = FALSE, na = "")

# Strict standalone calibration panels replacing the former L/M panels.
standalone_colours <- c(D_dimer = "#3F78A8", SDC1 = "#D86613")
standalone_labels <- c(D_dimer = "D-dimer", SDC1 = "SDC1")

make_standalone_calibration_panel <- function(population, tag) {
  result <- standalone[[population]]
  cal <- result$calibration
  cal$model <- factor(cal$model, levels = c("D_dimer", "SDC1"))
  cal$corrected_plot <- clip01(cal$corrected_observed)
  cal$apparent_plot <- clip01(cal$apparent_observed)
  cal$lower_plot <- clip01(cal$lower_observed)
  cal$upper_plot <- clip01(cal$upper_observed)
  cal$band_estimable <- is.finite(cal$lower_observed) & is.finite(cal$upper_observed)

  top <- ggplot(cal, aes(predicted_probability, group = model)) +
    geom_ribbon(
      data = cal[cal$band_estimable, , drop = FALSE],
      aes(ymin = lower_plot, ymax = upper_plot, fill = model),
      alpha = 0.14, colour = NA
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#8A9299", linewidth = 0.55) +
    geom_line(aes(y = apparent_plot, colour = model), linetype = "dotted", linewidth = 0.48, alpha = 0.72) +
    geom_line(aes(y = corrected_plot, colour = model), linewidth = 0.88) +
    scale_colour_manual(values = standalone_colours, labels = standalone_labels, breaks = names(standalone_labels), name = NULL) +
    scale_fill_manual(values = standalone_colours, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_label, expand = c(0, 0)) +
    coord_equal() +
    labs(
      x = "Predicted probability of AAD", y = "Observed AAD proportion",
      title = unname(population_labels[population]),
      subtitle = "Solid, optimism-corrected; dotted, apparent; dashed, ideal"
    ) +
    theme_pub(6.4) +
    theme(
      legend.position = "top", legend.direction = "horizontal",
      legend.key.width = grid::unit(8, "mm"),
      plot.title = element_text(size = 7.3)
    )

  scalar <- result$scalar[result$scalar$metric %in% c("Calibration intercept", "Calibration slope"), , drop = FALSE]
  scalar$model <- factor(scalar$model, levels = rev(c("D_dimer", "SDC1")))
  scalar$model_label <- factor(standalone_labels[as.character(scalar$model)], levels = standalone_labels[rev(c("D_dimer", "SDC1"))])
  scalar$metric_label <- factor(
    scalar$metric,
    levels = c("Calibration intercept", "Calibration slope"),
    labels = c("Intercept (target 0)", "Slope (target 1)")
  )
  target <- data.frame(
    metric_label = factor(c("Intercept (target 0)", "Slope (target 1)"), levels = levels(scalar$metric_label)),
    target = c(0, 1)
  )
  bottom <- ggplot(scalar, aes(corrected, model_label, colour = model)) +
    geom_vline(data = target, aes(xintercept = target), linetype = "22", colour = "#8A9299", linewidth = 0.5) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0, linewidth = 0.58) +
    geom_point(size = 2.0) +
    facet_wrap(~metric_label, scales = "free_x", ncol = 1) +
    scale_colour_manual(values = standalone_colours, guide = "none") +
    labs(x = "Optimism-corrected estimate (ABCLOC 95% CI)", y = NULL) +
    theme_pub(6.3) +
    theme(
      strip.background = element_rect(fill = "#F1F3F4", colour = NA),
      panel.spacing.y = grid::unit(1.0, "mm")
    )

  top <- tag_panel(top, tag)
  top / bottom + plot_layout(heights = c(1.45, 0.75))
}

p_l <- make_standalone_calibration_panel("symptomatic", "L")
p_m <- make_standalone_calibration_panel("aad_vs_pe", "M")
save_panel(p_l, file.path(supp_dir, "Supplementary_Figure_L_strict_standalone_calibration"), width_mm = 86, height_mm = 112)
save_panel(p_m, file.path(supp_dir, "Supplementary_Figure_M_strict_standalone_calibration"), width_mm = 86, height_mm = 112)

standalone_calibration <- rbind(standalone$symptomatic$calibration, standalone$aad_vs_pe$calibration)
standalone_scalar <- rbind(standalone$symptomatic$scalar, standalone$aad_vs_pe$scalar)
standalone_difference <- rbind(standalone$symptomatic$differences, standalone$aad_vs_pe$differences)
write.csv(standalone_calibration, file.path(source_dir, "supplementary_LM_standalone_calibration_source_data.csv"), row.names = FALSE, na = "")
write.csv(standalone_scalar, file.path(stat_dir, "supplementary_LM_standalone_scalar_performance.csv"), row.names = FALSE, na = "")
write.csv(standalone_difference, file.path(stat_dir, "supplementary_LM_standalone_paired_differences.csv"), row.names = FALSE, na = "")

expected_main_stems <- c(
  file.path(figure5_dir, paste0("Figure_5", LETTERS[1:9])),
  file.path(figure6_dir, paste0("Figure_6", LETTERS[1:9]))
)
expected_supp_stems <- c(
  file.path(supp_dir, "Supplementary_Figure_2D"),
  file.path(supp_dir, "Supplementary_Figure_L_strict_standalone_calibration"),
  file.path(supp_dir, "Supplementary_Figure_M_strict_standalone_calibration")
)
all_stems <- c(expected_main_stems, expected_supp_stems)
all_exports <- unlist(lapply(all_stems, function(stem) paste0(stem, c(".svg", ".pdf", ".tiff", ".png"))))

qa <- data.frame(
  check = c(
    "figure5_n_329", "figure6_n_142",
    "all_joint_scalar_resamples_1000", "all_joint_improvement_resamples_1000",
    "oof_predictions_complete", "oof_auc_close_to_abcloc_auc",
    "selection_rank1_295", "selection_retained_296",
    "standalone_scalar_resamples_1000", "all_panel_exports_present"
  ),
  pass = c(
    nrow(d_sym) == 329L,
    nrow(d_pe) == 142L,
    all(all_scalar$successful_resamples == 1000L),
    all(all_improvements$successful_resamples == 1000L),
    all(is.finite(all_oof$repeated_oof_probability)),
    max(abs(all_oof_auc$oof_auc - all_oof_auc$corrected_auc)) < 0.025,
    sum(selection_draws$sdc1_retained & selection_draws$sdc1_rank == 1, na.rm = TRUE) == 295L,
    sum(selection_draws$sdc1_retained, na.rm = TRUE) == 296L,
    all(standalone_scalar$successful_resamples == 1000L),
    all(file.exists(all_exports))
  ),
  stringsAsFactors = FALSE
)
write.csv(qa, file.path(qc_dir, "modified_panels_computational_qa.csv"), row.names = FALSE)

joint_panel_d_method <- if (calibration_split) {
  "Harrell optimism-corrected calibration intercept across all four models with ABCLOC 95% CI"
} else {
  "Harrell optimism-corrected intercept/slope with ABCLOC 95% CI"
}
joint_panel_e_method <- if (calibration_split) {
  "Harrell optimism-corrected calibration slope across all four models with ABCLOC 95% CI"
} else {
  "Paired Harrell optimism-corrected improvement vs ADD-RS with ABCLOC 95% CI"
}

manifest <- data.frame(
  figure = c(
    paste0("Figure 5", LETTERS[1:9]),
    paste0("Figure 6", LETTERS[1:9]),
    "Supplementary Figure 2D", "Supplementary Figure L", "Supplementary Figure M"
  ),
  stem = all_stems,
  method = c(
    "Harrell LOWESS optimism-corrected calibration with ABCLOC 95% CI",
    "Repeated stratified 10-fold OOF ROC; numeric AUROC from Harrell ABCLOC",
    "Harrell optimism-corrected scalar performance with ABCLOC 95% CI",
    joint_panel_d_method,
    joint_panel_e_method,
    "Original-sample likelihood-ratio information under locked RCS specification",
    "Patient-level mean repeated OOF probability update",
    "Distribution of patient-level mean repeated OOF probability changes",
    "Correct-direction update with stratified participant-bootstrap 95% CI",
    "Harrell LOWESS optimism-corrected calibration with ABCLOC 95% CI",
    "Repeated stratified 10-fold OOF ROC; numeric AUROC from Harrell ABCLOC",
    "Harrell optimism-corrected scalar performance with ABCLOC 95% CI",
    joint_panel_d_method,
    joint_panel_e_method,
    "Original-sample likelihood-ratio information under locked RCS specification",
    "Patient-level mean repeated OOF probability update",
    "Distribution of patient-level mean repeated OOF probability changes",
    "Correct-direction update with stratified participant-bootstrap 95% CI",
    "Fully nested selection bootstrap: filter, 50 fresh imputations, two contrasts, Rubin pooling and ranking",
    "Standalone D-dimer vs SDC1 Harrell calibration in symptomatic cohort",
    "Standalone D-dimer vs SDC1 Harrell calibration in AAD-vs-PE cohort"
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(figure_root, "PANEL_MANIFEST.csv"), row.names = FALSE)

summary_lines <- c(
  "Modified figure panels for manual assembly",
  paste0("Computational QA: ", if (all(qa$pass)) "PASS" else "FAIL"),
  "Figure 5: A-I; symptomatic cohort n=329 (99 AAD, 230 symptomatic non-AAD)",
  "Figure 6: A-I; AAD versus PE n=142 (99 AAD, 43 PE)",
  "Joint scalar estimates and intervals: 1,000 outcome-stratified Harrell optimism bootstraps with ABCLOC intervals",
  if (calibration_split) {
    "Panels D/E: calibration intercept/slope separated and ADD-RS alone included; former paired-improvement panel removed"
  } else {
    "Panels D/E: original combined-calibration and paired-improvement layout"
  },
  "Calibration: LOWESS smooth, optimism-corrected; pointwise ABCLOC confidence limits only where estimable",
  "ROC and probability updating: 50 repeats of stratified 10-fold OOF prediction using the locked RCS bases",
  "Supplementary Figure 2D: 295/300 retained and rank 1; 1/300 retained and rank >1; 4/300 below detection filter",
  "Supplementary L/M: strict standalone D-dimer versus SDC1 Harrell calibration panels",
  paste0("Output directory: ", figure_root)
)
writeLines(summary_lines, file.path(log_dir, "modified_panels_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
if (!all(qa$pass)) quit(status = 1)

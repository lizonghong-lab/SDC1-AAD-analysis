options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

analysis_file <- file.path(root, "02_analysis_data", "elisa_evaluation_primary.csv")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
figure_dir <- file.path(root, "04_results", "figures", "main")
panel_dir <- file.path(root, "04_results", "figures", "panels", "figure_4_revised")
log_dir <- file.path(root, "04_results", "logs")
qa_dir <- file.path(root, "07_quality_control")
for (directory in c(result_dir, source_dir, figure_dir, panel_dir, log_dir, qa_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

set.seed(20260722)

data <- read.csv(analysis_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("row_uid", "source_id", "cluster_id", "cluster_rule", "group", "aad", "sdc1_source_value", "d_dimer_mg_l")
missing_columns <- setdiff(required, names(data))
if (length(missing_columns)) stop("Missing required columns: ", paste(missing_columns, collapse = ", "))

data <- data[data$group %in% c("AAD", "AMI", "Angina", "PE"), required, drop = FALSE]
data$aad <- as.integer(data$aad)
data$sdc1_source_value <- as.numeric(data$sdc1_source_value)
data$d_dimer_mg_l <- as.numeric(data$d_dimer_mg_l)
analysis_required <- c("row_uid", "cluster_id", "group", "aad", "sdc1_source_value", "d_dimer_mg_l")
if (any(!complete.cases(data[, analysis_required]))) {
  stop("Standalone biomarker variables or resolved patient keys are incomplete")
}
if (any(data$sdc1_source_value <= 0) || any(data$d_dimer_mg_l <= 0)) stop("Biomarker values must be positive")
if (nrow(data) != 329L || sum(data$aad == 1) != 99L || sum(data$aad == 0) != 230L) {
  stop("Unexpected symptomatic-cohort counts")
}
if (length(unique(data$cluster_id)) != 329L) stop("Expected 329 author-resolved patient keys")

data$log2_d_dimer <- log2(data$d_dimer_mg_l)
data$log2_sdc1 <- log2(data$sdc1_source_value)

model_formulas <- list(
  D_dimer = aad ~ splines::ns(log2_d_dimer, df = 2),
  SDC1 = aad ~ splines::ns(log2_sdc1, df = 2)
)
model_order <- names(model_formulas)
model_labels <- c(D_dimer = "D-dimer", SDC1 = "SDC1")
model_palette <- c("D-dimer" = "#3F78A8", "SDC1" = "#D55E00")

clip_probability <- function(p, eps = 1e-8) pmin(pmax(p, eps), 1 - eps)

auc_rank <- function(y, score) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calculate_metrics <- function(y, p) {
  p <- clip_probability(p)
  c(
    Brier = mean((y - p)^2),
    Log_loss = -mean(y * log(p) + (1 - y) * log(1 - p)),
    C_index = auc_rank(y, p)
  )
}

fit_model <- function(formula, training_data) {
  fit <- suppressWarnings(glm(
    formula,
    data = training_data,
    family = binomial(),
    control = glm.control(maxit = 100, epsilon = 1e-9)
  ))
  if (!isTRUE(fit$converged)) stop("Logistic model did not converge")
  fit
}

make_stratified_cluster_folds <- function(dataset, k = 10L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fold_map <- integer(0)
  for (outcome_value in c(0, 1)) {
    clusters <- unique(dataset$cluster_id[dataset$aad == outcome_value])
    clusters <- sample(clusters, length(clusters), replace = FALSE)
    assignments <- rep(seq_len(k), length.out = length(clusters))
    fold_map <- c(fold_map, setNames(assignments, clusters))
  }
  unname(fold_map[dataset$cluster_id])
}

cross_validated_predictions <- function(dataset, formulas, k = 10L, seed = NULL) {
  folds <- make_stratified_cluster_folds(dataset, k = k, seed = seed)
  predictions <- matrix(NA_real_, nrow = nrow(dataset), ncol = length(formulas))
  colnames(predictions) <- names(formulas)
  for (fold in seq_len(k)) {
    test_index <- which(folds == fold)
    train_index <- which(folds != fold)
    training_data <- dataset[train_index, , drop = FALSE]
    test_data <- dataset[test_index, , drop = FALSE]
    for (model_name in names(formulas)) {
      fit <- fit_model(formulas[[model_name]], training_data)
      predictions[test_index, model_name] <- predict(fit, newdata = test_data, type = "response")
    }
  }
  if (any(!is.finite(predictions))) stop("Non-finite held-out predictions")
  predictions
}

quantile_interval <- function(x) {
  c(
    lower_95 = unname(quantile(x, 0.025, na.rm = TRUE)),
    upper_95 = unname(quantile(x, 0.975, na.rm = TRUE))
  )
}

analyse_population <- function(dataset, population_name, seed_offset, n_repeats = 50L, n_bootstrap = 1000L) {
  message("Analysing ", population_name, ": ", nrow(dataset), " participants")
  repeated_rows <- vector("list", n_repeats)
  prediction_sum <- matrix(0, nrow = nrow(dataset), ncol = length(model_order), dimnames = list(NULL, model_order))

  for (repeat_index in seq_len(n_repeats)) {
    predictions <- cross_validated_predictions(
      dataset,
      model_formulas,
      k = 10L,
      seed = seed_offset + repeat_index
    )
    prediction_sum <- prediction_sum + predictions
    metric_rows <- lapply(model_order, function(model_name) {
      metrics <- calculate_metrics(dataset$aad, predictions[, model_name])
      data.frame(
        population = population_name,
        repeat_id = repeat_index,
        model = model_name,
        metric = names(metrics),
        value = as.numeric(metrics),
        stringsAsFactors = FALSE
      )
    })
    repeated_rows[[repeat_index]] <- do.call(rbind, metric_rows)
  }

  repeated_metrics <- do.call(rbind, repeated_rows)
  point <- aggregate(value ~ population + model + metric, repeated_metrics, mean)
  names(point)[names(point) == "value"] <- "estimate"

  case_indices <- which(dataset$aad == 1)
  control_indices <- which(dataset$aad == 0)
  bootstrap_rows <- vector("list", n_bootstrap)
  for (bootstrap_index in seq_len(n_bootstrap)) {
    set.seed(seed_offset + 100000L + bootstrap_index)
    sampled_indices <- c(
      sample(case_indices, length(case_indices), replace = TRUE),
      sample(control_indices, length(control_indices), replace = TRUE)
    )
    bootstrap_data <- dataset[sampled_indices, , drop = FALSE]
    predictions <- cross_validated_predictions(
      bootstrap_data,
      model_formulas,
      k = 10L,
      seed = seed_offset + 200000L + bootstrap_index
    )
    metric_rows <- lapply(model_order, function(model_name) {
      metrics <- calculate_metrics(bootstrap_data$aad, predictions[, model_name])
      data.frame(
        population = population_name,
        bootstrap_replicate = bootstrap_index,
        model = model_name,
        metric = names(metrics),
        value = as.numeric(metrics),
        stringsAsFactors = FALSE
      )
    })
    bootstrap_rows[[bootstrap_index]] <- do.call(rbind, metric_rows)
    if (bootstrap_index %% 100L == 0L) {
      message(population_name, ": completed bootstrap ", bootstrap_index, "/", n_bootstrap)
    }
  }
  bootstrap_metrics <- do.call(rbind, bootstrap_rows)

  summary_rows <- vector("list", nrow(point))
  for (row_index in seq_len(nrow(point))) {
    point_row <- point[row_index, , drop = FALSE]
    draws <- bootstrap_metrics$value[
      bootstrap_metrics$model == point_row$model &
        bootstrap_metrics$metric == point_row$metric
    ]
    interval <- quantile_interval(draws)
    summary_rows[[row_index]] <- data.frame(
      population = population_name,
      model = point_row$model,
      metric = point_row$metric,
      estimate = point_row$estimate,
      lower_95 = interval[["lower_95"]],
      upper_95 = interval[["upper_95"]],
      n_total = nrow(dataset),
      aad_cases = sum(dataset$aad == 1),
      non_aad_controls = sum(dataset$aad == 0),
      stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, summary_rows)

  difference_rows <- lapply(unique(summary$metric), function(metric_name) {
    sdc1_point <- summary$estimate[summary$model == "SDC1" & summary$metric == metric_name]
    d_dimer_point <- summary$estimate[summary$model == "D_dimer" & summary$metric == metric_name]
    sdc1_draws <- bootstrap_metrics$value[bootstrap_metrics$model == "SDC1" & bootstrap_metrics$metric == metric_name]
    d_dimer_draws <- bootstrap_metrics$value[bootstrap_metrics$model == "D_dimer" & bootstrap_metrics$metric == metric_name]
    difference_draws <- sdc1_draws - d_dimer_draws
    interval <- quantile_interval(difference_draws)
    data.frame(
      population = population_name,
      metric = metric_name,
      comparison = "SDC1 minus D-dimer",
      difference = sdc1_point - d_dimer_point,
      lower_95 = interval[["lower_95"]],
      upper_95 = interval[["upper_95"]],
      stringsAsFactors = FALSE
    )
  })

  mean_predictions <- prediction_sum / n_repeats
  patient_predictions <- data.frame(
    population = population_name,
    row_uid = dataset$row_uid,
    source_id = dataset$source_id,
    cluster_id = dataset$cluster_id,
    group = dataset$group,
    aad = dataset$aad,
    probability_D_dimer = mean_predictions[, "D_dimer"],
    probability_SDC1 = mean_predictions[, "SDC1"],
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    differences = do.call(rbind, difference_rows),
    repeated_metrics = repeated_metrics,
    bootstrap_metrics = bootstrap_metrics,
    patient_predictions = patient_predictions
  )
}

population_data <- list(
  "All symptomatic comparators" = data,
  "AAD versus PE" = data[data$aad == 1 | data$group == "PE", , drop = FALSE]
)
if (nrow(population_data[["AAD versus PE"]]) != 142L) stop("Expected 142 AAD-versus-PE participants")

analysis_all <- analyse_population(population_data[[1]], names(population_data)[1], 20264000L)
analysis_pe <- analyse_population(population_data[[2]], names(population_data)[2], 20265000L)

validated_summary <- rbind(analysis_all$summary, analysis_pe$summary)
paired_differences <- rbind(analysis_all$differences, analysis_pe$differences)
repeated_metrics <- rbind(analysis_all$repeated_metrics, analysis_pe$repeated_metrics)
bootstrap_metrics <- rbind(analysis_all$bootstrap_metrics, analysis_pe$bootstrap_metrics)
patient_predictions <- rbind(analysis_all$patient_predictions, analysis_pe$patient_predictions)

write.csv(validated_summary, file.path(result_dir, "fig4_standalone_validated_probability_performance.csv"), row.names = FALSE, na = "")
write.csv(paired_differences, file.path(result_dir, "fig4_standalone_validated_probability_differences.csv"), row.names = FALSE, na = "")
write.csv(repeated_metrics, file.path(source_dir, "fig4_standalone_repeated_cv_metric_draws.csv"), row.names = FALSE, na = "")
write.csv(bootstrap_metrics, file.path(source_dir, "fig4_standalone_full_pipeline_bootstrap_metric_draws.csv"), row.names = FALSE, na = "")
write.csv(patient_predictions, file.path(source_dir, "fig4_standalone_patient_level_repeated_oof_predictions.csv"), row.names = FALSE, na = "")

# Load the locked A-H panels. This also recreates the previous standalone figure,
# but no existing final PDF is overwritten by the revised export below.
source(file.path(root, "03_code", "R", "05_render_final_figure4_head_to_head.R"), local = FALSE)

validated_summary$model_label <- factor(
  model_labels[validated_summary$model],
  levels = rev(unname(model_labels[model_order]))
)
validated_summary$metric_label <- factor(
  validated_summary$metric,
  levels = c("Brier", "Log_loss", "C_index"),
  labels = c("Brier score", "Log loss", "C-index")
)

make_probability_metric_panel <- function(population_name, title) {
  plot_data <- validated_summary[validated_summary$population == population_name, , drop = FALSE]
  subtitle <- sprintf(
    "n=%d; %d AAD / %d non-AAD; 50 x 10-fold CV",
    plot_data$n_total[1], plot_data$aad_cases[1], plot_data$non_aad_controls[1]
  )
  ggplot(plot_data, aes(estimate, model_label, colour = model_label)) +
    geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0, linewidth = 0.55) +
    geom_point(size = 2.0) +
    facet_wrap(~metric_label, scales = "free_x", ncol = 1) +
    scale_colour_manual(values = biomarker_palette, guide = "none") +
    labs(
      x = "Full-pipeline validated estimate (95% CI)",
      y = NULL,
      title = title,
      subtitle = paste0(subtitle, "\nLower Brier/log loss and higher C-index are better")
    ) +
    theme_pub() +
    theme(
      strip.background = element_rect(fill = "#F1F3F4", colour = NA),
      strip.text = element_text(size = 6.0, face = "bold"),
      panel.spacing.y = grid::unit(1.2, "mm"),
      axis.text.y = element_text(size = 5.8),
      plot.subtitle = element_text(size = 5.4)
    )
}

p_i <- make_probability_metric_panel("All symptomatic comparators", "Validated single-marker probability performance")
p_j <- make_probability_metric_panel("AAD versus PE", "Validated single-marker probability performance: AAD vs PE")

all_comparator_row <- (p_g_head | p_i) + plot_layout(widths = c(1.05, 0.95))
pe_row <- (p_h_head | p_j) + plot_layout(widths = c(1.05, 0.95))
figure4_revised <- top_row / middle_row / all_comparator_row / pe_row +
  plot_layout(heights = c(0.78, 0.70, 1.02, 1.02)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 8, family = "Arial", margin = margin(r = 3, b = 2)))
  )

figure_stem <- file.path(figure_dir, "Figure_4_standalone_biomarker_comparison_revised")
save_figure(figure4_revised, figure_stem, width_mm = 183, height_mm = 232, dpi = 600)

revised_panels <- list(I = p_i, J = p_j)
for (panel_name in names(revised_panels)) {
  save_figure(
    revised_panels[[panel_name]],
    file.path(panel_dir, paste0("Figure_4", panel_name)),
    width_mm = 88,
    height_mm = 86,
    dpi = 600
  )
}

interval_gap <- pmax(validated_summary$lower_95 - validated_summary$estimate, validated_summary$estimate - validated_summary$upper_95, 0)
qa_checks <- c(
  symptomatic_n_329 = nrow(population_data[[1]]) == 329L,
  aad_vs_pe_n_142 = nrow(population_data[[2]]) == 142L,
  unique_patient_keys_329 = length(unique(population_data[[1]]$cluster_id)) == 329L,
  all_point_estimates_inside_intervals = max(interval_gap) < 1e-12,
  expected_bootstrap_draw_count = nrow(bootstrap_metrics) == 2L * 1000L * 2L * 3L,
  no_nonfinite_estimates = all(is.finite(unlist(validated_summary[, c("estimate", "lower_95", "upper_95")]))),
  revised_exports_present = all(file.exists(paste0(figure_stem, c(".pdf", ".svg", ".tiff", ".png"))))
)
qa_table <- data.frame(check = names(qa_checks), pass = as.logical(qa_checks), stringsAsFactors = FALSE)
write.csv(qa_table, file.path(qa_dir, "FINAL_FIGURE4_REVISED_QA.csv"), row.names = FALSE)

summary_lines <- c(
  "Revised Figure 4: standalone SDC1 versus D-dimer",
  paste0("Computational QA: ", if (all(qa_checks)) "PASS" else "FAIL"),
  "ROC panels retain prespecified continuous raw-marker AUROCs from the locked same-patient analysis.",
  "Probability metrics use separate single-marker logistic models with 2-df natural cubic splines.",
  "Point estimates use 50 repeats of stratified 10-fold cross-validation.",
  "Intervals use 1,000 full-pipeline outcome-stratified patient bootstraps with fold reconstruction.",
  "No clinical covariates, combined biomarkers, optimized cutoffs, Youden statistics or forced biomarker dichotomization were used.",
  paste0("Figure: ", figure_stem)
)
writeLines(summary_lines, file.path(log_dir, "fig4_revised_standalone_model_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "fig4_revised_sessionInfo.txt"), useBytes = TRUE)

cat(paste(summary_lines, collapse = "\n"), "\n")
if (!all(qa_checks)) quit(status = 1)

options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source_dir <- file.path(root, "04_results", "source_data")
result_dir <- file.path(root, "04_results", "statistical_results")
qa_dir <- file.path(root, "07_quality_control")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

pred <- read.csv(file.path(source_dir, "fig6_repeated_cv_prediction_draws.csv"), check.names = FALSE)
patient <- read.csv(file.path(source_dir, "fig6_patient_level_repeated_oof_predictions.csv"), check.names = FALSE)
boot_metric <- read.csv(file.path(source_dir, "fig6_full_pipeline_bootstrap_metric_draws.csv"), check.names = FALSE)
boot_delta <- read.csv(file.path(source_dir, "fig6_full_pipeline_bootstrap_paired_delta_draws.csv"), check.names = FALSE)
reported <- read.csv(file.path(result_dir, "fig6_validated_performance.csv"), check.names = FALSE)
reported_delta <- read.csv(file.path(result_dir, "fig6_paired_incremental_performance.csv"), check.names = FALSE)
reported_lr <- read.csv(file.path(result_dir, "fig6_sdc1_added_value_lr_test.csv"), check.names = FALSE)
analysis_data <- read.csv(file.path(root, "02_analysis_data", "elisa_evaluation_primary.csv"), check.names = FALSE)

clip <- function(x) pmin(pmax(x, 1e-8), 1 - 1e-8)

auc_manual <- function(y, score) {
  n_case <- sum(y == 1)
  n_control <- sum(y == 0)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n_case * (n_case + 1) / 2) / (n_case * n_control)
}

metric_manual <- function(y, probability) {
  probability <- clip(probability)
  linear_predictor <- qlogis(probability)
  calibration_intercept <- unname(coef(glm(y ~ 1, family = binomial(), offset = linear_predictor))[1])
  calibration_slope <- unname(coef(glm(y ~ linear_predictor, family = binomial()))[2])
  c(
    Brier = sum((y - probability)^2) / length(y),
    Log_loss = -sum(y * log(probability) + (1 - y) * log(1 - probability)) / length(y),
    C_index = auc_manual(y, probability),
    Calibration_intercept = calibration_intercept,
    Calibration_slope = calibration_slope
  )
}

models <- unique(pred$model)
repeats <- sort(unique(pred$repeat_id))
manual_metric_rows <- list()
manual_delta_rows <- list()

for (repeat_id in repeats) {
  metrics_by_model <- list()
  for (model_name in models) {
    x <- pred[pred$repeat_id == repeat_id & pred$model == model_name, , drop = FALSE]
    metrics <- metric_manual(x$aad, x$prediction)
    metrics_by_model[[model_name]] <- metrics
    manual_metric_rows[[length(manual_metric_rows) + 1L]] <- data.frame(
      repeat_id = repeat_id,
      model = model_name,
      metric = names(metrics),
      value = as.numeric(metrics)
    )
  }
  benchmark <- metrics_by_model[["ADD_RS_plus_D_dimer"]]
  full <- metrics_by_model[["ADD_RS_plus_D_dimer_plus_SDC1"]]
  manual_delta_rows[[length(manual_delta_rows) + 1L]] <- data.frame(
    repeat_id = repeat_id,
    metric = c("Brier improvement", "Log-loss improvement", "C-index improvement"),
    value = c(
      benchmark[["Brier"]] - full[["Brier"]],
      benchmark[["Log_loss"]] - full[["Log_loss"]],
      full[["C_index"]] - benchmark[["C_index"]]
    )
  )
}

manual_metrics <- do.call(rbind, manual_metric_rows)
manual_deltas <- do.call(rbind, manual_delta_rows)
manual_point <- aggregate(value ~ model + metric, manual_metrics, mean)
names(manual_point)[3] <- "manual_estimate"
manual_delta_point <- aggregate(value ~ metric, manual_deltas, mean)
names(manual_delta_point)[2] <- "manual_estimate"

metric_compare <- merge(reported, manual_point, by = c("model", "metric"), all = TRUE)
metric_compare$estimate_abs_difference <- abs(metric_compare$estimate - metric_compare$manual_estimate)

bootstrap_keys <- unique(boot_metric[, c("model", "metric")])
bootstrap_quantiles <- do.call(rbind, lapply(seq_len(nrow(bootstrap_keys)), function(index) {
  model_name <- bootstrap_keys$model[index]
  metric_name <- bootstrap_keys$metric[index]
  draws <- boot_metric$value[
    boot_metric$model == model_name & boot_metric$metric == metric_name
  ]
  interval <- unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE))
  data.frame(
    model = model_name,
    metric = metric_name,
    manual_lower_95 = interval[1],
    manual_upper_95 = interval[2]
  )
}))
metric_compare <- merge(metric_compare, bootstrap_quantiles, by = c("model", "metric"), all = TRUE)
metric_compare$lower_abs_difference <- abs(metric_compare$lower_95 - metric_compare$manual_lower_95)
metric_compare$upper_abs_difference <- abs(metric_compare$upper_95 - metric_compare$manual_upper_95)

delta_compare <- merge(reported_delta, manual_delta_point, by = "metric", all = TRUE)
delta_compare$estimate_abs_difference <- abs(delta_compare$estimate - delta_compare$manual_estimate)
manual_delta_interval <- do.call(rbind, lapply(unique(boot_delta$metric), function(metric_name) {
  draws <- boot_delta$value[boot_delta$metric == metric_name]
  interval <- unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE))
  data.frame(metric = metric_name, manual_lower_95 = interval[1], manual_upper_95 = interval[2])
}))
delta_compare <- merge(delta_compare, manual_delta_interval, by = "metric", all = TRUE)
delta_compare$lower_abs_difference <- abs(delta_compare$lower_95 - delta_compare$manual_lower_95)
delta_compare$upper_abs_difference <- abs(delta_compare$upper_95 - delta_compare$manual_upper_95)

prediction_means <- aggregate(prediction ~ row_uid + model, pred, mean)
prediction_mean_wide <- reshape(prediction_means, idvar = "row_uid", timevar = "model", direction = "wide")
prediction_columns <- paste0("prob_", models)
patient_compare <- merge(
  patient[, c("row_uid", prediction_columns)],
  prediction_mean_wide,
  by = "row_uid",
  all = TRUE
)
prediction_differences <- unlist(lapply(models, function(model_name) {
  abs(
    patient_compare[[paste0("prob_", model_name)]] -
      patient_compare[[paste0("prediction.", model_name)]]
  )
}))

analysis_data <- analysis_data[analysis_data$group %in% c("AAD", "AMI", "Angina", "PE"), , drop = FALSE]
analysis_data$add_rs <- as.numeric(analysis_data$add_rs_corrected)
analysis_data$log2_d_dimer <- log2(as.numeric(analysis_data$d_dimer_mg_l))
analysis_data$log2_sdc1 <- log2(as.numeric(analysis_data$sdc1_source_value))
benchmark_fit <- glm(
  aad ~ add_rs + splines::ns(log2_d_dimer, df = 2),
  data = analysis_data,
  family = binomial()
)
full_fit <- glm(
  aad ~ add_rs + splines::ns(log2_d_dimer, df = 2) + splines::ns(log2_sdc1, df = 2),
  data = analysis_data,
  family = binomial()
)
manual_lr <- 2 * (as.numeric(logLik(full_fit)) - as.numeric(logLik(benchmark_fit)))
manual_lr_df <- attr(logLik(full_fit), "df") - attr(logLik(benchmark_fit), "df")

checks <- c(
  repeated_prediction_rows_65800 = nrow(pred) == 50L * 4L * 329L,
  repeated_cv_point_estimates_match = max(metric_compare$estimate_abs_difference) < 2e-6,
  bootstrap_metric_intervals_match = max(c(
    metric_compare$lower_abs_difference,
    metric_compare$upper_abs_difference
  )) < 1e-12,
  paired_delta_point_estimates_match = max(delta_compare$estimate_abs_difference) < 1e-12,
  paired_delta_intervals_match = max(c(
    delta_compare$lower_abs_difference,
    delta_compare$upper_abs_difference
  )) < 1e-12,
  patient_probabilities_equal_50_repeat_means = max(prediction_differences) < 1e-12,
  lr_increment_matches = abs(manual_lr - reported_lr$lr_increment[1]) < 1e-10,
  lr_degrees_of_freedom_match = manual_lr_df == reported_lr$df[1]
)

write.csv(metric_compare, file.path(qa_dir, "FIGURE6_INDEPENDENT_METRIC_RECALCULATION.csv"), row.names = FALSE)
write.csv(delta_compare, file.path(qa_dir, "FIGURE6_INDEPENDENT_DELTA_RECALCULATION.csv"), row.names = FALSE)
write.csv(data.frame(check = names(checks), pass = unname(checks)), file.path(qa_dir, "FIGURE6_INDEPENDENT_NUMERICAL_AUDIT.csv"), row.names = FALSE)

report <- c(
  "# Figure 6 independent numerical audit",
  "",
  paste0("Overall status: **", if (all(checks)) "PASS" else "FAIL", "**"),
  "",
  "The audit recalculated all repeated-CV point estimates and bootstrap percentile intervals directly from the saved prediction and resampling draws, without calling the Figure 6 analysis functions.",
  "",
  paste0("- Repeated held-out prediction rows: ", format(nrow(pred), big.mark = ",")),
  paste0("- Maximum point-estimate discrepancy: ", format(max(metric_compare$estimate_abs_difference), scientific = TRUE)),
  paste0("- Maximum metric-interval discrepancy: ", format(max(c(metric_compare$lower_abs_difference, metric_compare$upper_abs_difference)), scientific = TRUE)),
  paste0("- Maximum paired-delta discrepancy: ", format(max(c(delta_compare$estimate_abs_difference, delta_compare$lower_abs_difference, delta_compare$upper_abs_difference)), scientific = TRUE)),
  paste0("- Maximum patient-level mean-prediction discrepancy: ", format(max(prediction_differences), scientific = TRUE)),
  paste0("- Recalculated SDC1 LR increment: ", sprintf("%.10f", manual_lr), " on ", manual_lr_df, " df"),
  "",
  "The audit checks numerical reproduction only; cohort independence and transportability remain study-design limitations."
)
writeLines(report, file.path(qa_dir, "FIGURE6_INDEPENDENT_NUMERICAL_AUDIT.md"), useBytes = TRUE)

cat("Figure 6 independent numerical audit: ", if (all(checks)) "PASS" else "FAIL", "\n", sep = "")
if (!all(checks)) quit(status = 1)

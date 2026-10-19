options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

data_file <- file.path(root, "02_analysis_data", "elisa_evaluation_primary.csv")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
figure_dir <- file.path(root, "04_results", "figures", "main")
panel_dir <- file.path(root, "04_results", "figures", "panels", "figure_5_revised")
log_dir <- file.path(root, "04_results", "logs")
qa_dir <- file.path(root, "07_quality_control")
for (directory in c(result_dir, source_dir, figure_dir, panel_dir, log_dir, qa_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

set.seed(20260721)

data <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "row_uid", "source_id", "cluster_id", "cluster_rule", "group", "aad",
  "add_rs_source", "add_rs_corrected", "sdc1_source_value", "d_dimer_mg_l"
)
missing_columns <- setdiff(required, names(data))
if (length(missing_columns)) stop("Missing required columns: ", paste(missing_columns, collapse = ", "))

data <- data[data$group %in% c("AAD", "AMI", "Angina", "PE"), required, drop = FALSE]
numeric_columns <- c("aad", "add_rs_source", "add_rs_corrected", "sdc1_source_value", "d_dimer_mg_l")
data[numeric_columns] <- lapply(data[numeric_columns], as.numeric)

if (nrow(data) != 329L) stop("Expected 329 symptomatic records; found ", nrow(data))
if (length(unique(data$cluster_id)) != 329L) stop("Expected 329 unique author-resolved patient keys")
if (sum(data$aad == 1) != 99L || sum(data$aad == 0) != 230L) stop("Unexpected outcome counts")
analysis_required <- c(
  "row_uid", "cluster_id", "cluster_rule", "group", "aad",
  "add_rs_source", "add_rs_corrected", "sdc1_source_value", "d_dimer_mg_l"
)
if (any(!complete.cases(data[, analysis_required]))) {
  stop("Figure 6 analysis variables or derived patient keys are not complete")
}
if (any(data$sdc1_source_value <= 0) || any(data$d_dimer_mg_l <= 0)) stop("Biomarkers must be positive")
if (!identical(as.integer(table(factor(data$add_rs_corrected, levels = 1:3))), c(239L, 75L, 15L))) {
  stop("Unexpected corrected ADD-RS distribution")
}
if (sum(data$add_rs_source == 4 & data$add_rs_corrected == 3) != 1L) {
  stop("Expected exactly one author-confirmed ADD-RS 4-to-3 correction")
}

data$add_rs <- data$add_rs_corrected
data$log2_d_dimer <- log2(data$d_dimer_mg_l)
data$log2_sdc1 <- log2(data$sdc1_source_value)
data$outcome_label <- ifelse(data$aad == 1, "AAD", "Symptomatic non-AAD")

model_formulas <- list(
  ADD_RS = aad ~ add_rs,
  ADD_RS_plus_D_dimer = aad ~ add_rs + splines::ns(log2_d_dimer, df = 2),
  ADD_RS_plus_SDC1 = aad ~ add_rs + splines::ns(log2_sdc1, df = 2),
  ADD_RS_plus_D_dimer_plus_SDC1 = aad ~ add_rs + splines::ns(log2_d_dimer, df = 2) + splines::ns(log2_sdc1, df = 2)
)

model_order <- names(model_formulas)
model_labels <- c(
  ADD_RS = "ADD-RS",
  ADD_RS_plus_D_dimer = "ADD-RS + D-dimer",
  ADD_RS_plus_SDC1 = "ADD-RS + SDC1",
  ADD_RS_plus_D_dimer_plus_SDC1 = "ADD-RS + D-dimer + SDC1"
)
model_palette <- c(
  "ADD-RS" = "#9AA0A6",
  "ADD-RS + D-dimer" = "#3F78A8",
  "ADD-RS + SDC1" = "#D55E00",
  "ADD-RS + D-dimer + SDC1" = "#2A9D8F"
)
outcome_palette <- c("AAD" = "#BE4652", "Symptomatic non-AAD" = "#3F78A8")

clip_probability <- function(p, eps = 1e-8) pmin(pmax(p, eps), 1 - eps)

auc_rank <- function(y, score) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calibration_intercept_slope <- function(y, p) {
  lp <- qlogis(clip_probability(p))
  intercept_fit <- suppressWarnings(glm(y ~ 1, family = binomial(), offset = lp))
  slope_fit <- suppressWarnings(glm(y ~ lp, family = binomial()))
  c(intercept = unname(coef(intercept_fit)[1]), slope = unname(coef(slope_fit)[2]))
}

calculate_metrics <- function(y, p) {
  p <- clip_probability(p)
  calibration <- calibration_intercept_slope(y, p)
  c(
    Brier = mean((y - p)^2),
    Log_loss = -mean(y * log(p) + (1 - y) * log(1 - p)),
    C_index = auc_rank(y, p),
    Calibration_intercept = calibration[["intercept"]],
    Calibration_slope = calibration[["slope"]]
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
    if (!length(test_index) || !length(train_index)) stop("Empty CV fold")
    training_data <- dataset[train_index, , drop = FALSE]
    test_data <- dataset[test_index, , drop = FALSE]
    for (model_name in names(formulas)) {
      fit <- fit_model(formulas[[model_name]], training_data)
      predictions[test_index, model_name] <- predict(fit, newdata = test_data, type = "response")
    }
  }
  if (any(!is.finite(predictions))) stop("Non-finite held-out predictions")
  list(predictions = predictions, folds = folds)
}

smooth_calibration <- function(y, p, grid) {
  smoother <- lowess(p, y, f = 2 / 3, iter = 0)
  keep <- !duplicated(smoother$x)
  calibrated <- approx(
    x = smoother$x[keep],
    y = smoother$y[keep],
    xout = grid,
    rule = 2,
    ties = mean
  )$y
  pmin(pmax(calibrated, 0), 1)
}

message("Running 50 repeats of stratified 10-fold cross-validation...")
n_repeats <- 50L
prediction_sum <- matrix(0, nrow = nrow(data), ncol = length(model_formulas))
colnames(prediction_sum) <- model_order
oof_long <- vector("list", n_repeats)

for (repeat_index in seq_len(n_repeats)) {
  cv <- cross_validated_predictions(
    data,
    model_formulas,
    k = 10L,
    seed = 20260721L + repeat_index
  )
  prediction_sum <- prediction_sum + cv$predictions
  oof_long[[repeat_index]] <- data.frame(
    repeat_id = repeat_index,
    row_uid = rep(data$row_uid, times = length(model_order)),
    cluster_id = rep(data$cluster_id, times = length(model_order)),
    group = rep(data$group, times = length(model_order)),
    aad = rep(data$aad, times = length(model_order)),
    model = rep(model_order, each = nrow(data)),
    prediction = as.vector(cv$predictions),
    fold = rep(cv$folds, times = length(model_order)),
    stringsAsFactors = FALSE
  )
}

oof_predictions <- prediction_sum / n_repeats
oof_long <- do.call(rbind, oof_long)

# Performance is calculated within each repeat before averaging. Averaging
# patient probabilities first can introduce meaningless within-score ordering
# for the three-level ADD-RS model when fold-specific intercepts differ.
repeated_cv_metric_rows <- list()
repeated_cv_delta_rows <- list()
for (repeat_index in seq_len(n_repeats)) {
  repeat_metrics <- list()
  for (model_name in model_order) {
    repeat_data <- oof_long[
      oof_long$repeat_id == repeat_index & oof_long$model == model_name,
      ,
      drop = FALSE
    ]
    metrics <- calculate_metrics(repeat_data$aad, repeat_data$prediction)
    repeat_metrics[[model_name]] <- metrics
    repeated_cv_metric_rows[[length(repeated_cv_metric_rows) + 1L]] <- data.frame(
      repeat_id = repeat_index,
      model = model_name,
      metric = names(metrics),
      value = as.numeric(metrics),
      stringsAsFactors = FALSE
    )
  }
  benchmark_metrics <- repeat_metrics[["ADD_RS_plus_D_dimer"]]
  full_metrics <- repeat_metrics[["ADD_RS_plus_D_dimer_plus_SDC1"]]
  repeated_cv_delta_rows[[repeat_index]] <- data.frame(
    repeat_id = repeat_index,
    metric = c("Brier improvement", "Log-loss improvement", "C-index improvement"),
    value = c(
      benchmark_metrics[["Brier"]] - full_metrics[["Brier"]],
      benchmark_metrics[["Log_loss"]] - full_metrics[["Log_loss"]],
      full_metrics[["C_index"]] - benchmark_metrics[["C_index"]]
    ),
    stringsAsFactors = FALSE
  )
}
repeated_cv_metrics <- do.call(rbind, repeated_cv_metric_rows)
repeated_cv_deltas <- do.call(rbind, repeated_cv_delta_rows)

patient_predictions <- data.frame(
  row_uid = data$row_uid,
  source_id = data$source_id,
  cluster_id = data$cluster_id,
  cluster_rule = data$cluster_rule,
  group = data$group,
  aad = data$aad,
  add_rs_source = data$add_rs_source,
  add_rs_corrected = data$add_rs_corrected,
  d_dimer_mg_l = data$d_dimer_mg_l,
  total_plasma_sdc1 = data$sdc1_source_value,
  stringsAsFactors = FALSE
)
for (model_name in model_order) {
  patient_predictions[[paste0("prob_", model_name)]] <- oof_predictions[, model_name]
}

performance_point <- aggregate(
  value ~ model + metric,
  data = repeated_cv_metrics,
  FUN = mean
)
names(performance_point)[names(performance_point) == "value"] <- "estimate"

calibration_models <- c(
  "ADD_RS_plus_D_dimer",
  "ADD_RS_plus_SDC1",
  "ADD_RS_plus_D_dimer_plus_SDC1"
)
calibration_grids <- lapply(calibration_models, function(model_name) {
  p <- oof_predictions[, model_name]
  seq(quantile(p, 0.02), quantile(p, 0.98), length.out = 50)
})
names(calibration_grids) <- calibration_models

message("Running 1,000 full-pipeline outcome-stratified patient bootstraps...")
n_bootstrap <- 1000L
bootstrap_metric_rows <- vector("list", n_bootstrap)
bootstrap_delta_rows <- vector("list", n_bootstrap)
bootstrap_direction_rows <- vector("list", n_bootstrap)
bootstrap_calibration_rows <- vector("list", n_bootstrap)

case_indices <- which(data$aad == 1)
control_indices <- which(data$aad == 0)

for (bootstrap_index in seq_len(n_bootstrap)) {
  set.seed(20261721L + bootstrap_index)
  sampled_indices <- c(
    sample(case_indices, length(case_indices), replace = TRUE),
    sample(control_indices, length(control_indices), replace = TRUE)
  )
  bootstrap_data <- data[sampled_indices, , drop = FALSE]
  bootstrap_data$bootstrap_row <- seq_len(nrow(bootstrap_data))
  bootstrap_data$cluster_id <- as.character(bootstrap_data$cluster_id)

  cv_bootstrap <- cross_validated_predictions(
    bootstrap_data,
    model_formulas,
    k = 10L,
    seed = 20262721L + bootstrap_index
  )
  boot_predictions <- cv_bootstrap$predictions
  boot_y <- bootstrap_data$aad

  metric_frame <- list()
  for (model_name in model_order) {
    metrics <- calculate_metrics(boot_y, boot_predictions[, model_name])
    metric_frame[[model_name]] <- data.frame(
      bootstrap_replicate = bootstrap_index,
      model = model_name,
      metric = names(metrics),
      value = as.numeric(metrics),
      stringsAsFactors = FALSE
    )
  }
  bootstrap_metric_rows[[bootstrap_index]] <- do.call(rbind, metric_frame)

  benchmark_metrics <- calculate_metrics(boot_y, boot_predictions[, "ADD_RS_plus_D_dimer"])
  full_metrics <- calculate_metrics(boot_y, boot_predictions[, "ADD_RS_plus_D_dimer_plus_SDC1"])
  bootstrap_delta_rows[[bootstrap_index]] <- data.frame(
    bootstrap_replicate = bootstrap_index,
    metric = c("Brier improvement", "Log-loss improvement", "C-index improvement"),
    value = c(
      benchmark_metrics[["Brier"]] - full_metrics[["Brier"]],
      benchmark_metrics[["Log_loss"]] - full_metrics[["Log_loss"]],
      full_metrics[["C_index"]] - benchmark_metrics[["C_index"]]
    ),
    stringsAsFactors = FALSE
  )

  probability_change <-
    boot_predictions[, "ADD_RS_plus_D_dimer_plus_SDC1"] -
    boot_predictions[, "ADD_RS_plus_D_dimer"]
  correct_direction <- ifelse(boot_y == 1, probability_change > 0, probability_change < 0)
  bootstrap_direction_rows[[bootstrap_index]] <- data.frame(
    bootstrap_replicate = bootstrap_index,
    outcome = c("AAD", "Symptomatic non-AAD"),
    proportion = c(
      mean(correct_direction[boot_y == 1]),
      mean(correct_direction[boot_y == 0])
    ),
    stringsAsFactors = FALSE
  )

  calibration_frame <- list()
  for (model_name in calibration_models) {
    grid <- calibration_grids[[model_name]]
    calibrated <- smooth_calibration(boot_y, boot_predictions[, model_name], grid)
    calibration_frame[[model_name]] <- data.frame(
      bootstrap_replicate = bootstrap_index,
      model = model_name,
      predicted_probability = grid,
      calibrated_probability = calibrated,
      stringsAsFactors = FALSE
    )
  }
  bootstrap_calibration_rows[[bootstrap_index]] <- do.call(rbind, calibration_frame)

  if (bootstrap_index %% 50L == 0L) message("Completed bootstrap ", bootstrap_index, "/", n_bootstrap)
}

bootstrap_metrics <- do.call(rbind, bootstrap_metric_rows)
bootstrap_deltas <- do.call(rbind, bootstrap_delta_rows)
bootstrap_direction <- do.call(rbind, bootstrap_direction_rows)
bootstrap_calibration <- do.call(rbind, bootstrap_calibration_rows)

quantile_interval <- function(x) {
  c(lower_95 = unname(quantile(x, 0.025, na.rm = TRUE)), upper_95 = unname(quantile(x, 0.975, na.rm = TRUE)))
}

performance_summary_rows <- list()
for (model_name in model_order) {
  for (metric_name in unique(performance_point$metric)) {
    estimate <- performance_point$estimate[
      performance_point$model == model_name & performance_point$metric == metric_name
    ]
    draws <- bootstrap_metrics$value[
      bootstrap_metrics$model == model_name & bootstrap_metrics$metric == metric_name
    ]
    interval <- quantile_interval(draws)
    performance_summary_rows[[length(performance_summary_rows) + 1]] <- data.frame(
      model = model_name,
      metric = metric_name,
      estimate = estimate,
      lower_95 = interval[["lower_95"]],
      upper_95 = interval[["upper_95"]],
      stringsAsFactors = FALSE
    )
  }
}
performance_summary <- do.call(rbind, performance_summary_rows)

delta_point <- aggregate(
  value ~ metric,
  data = repeated_cv_deltas,
  FUN = mean
)
names(delta_point)[names(delta_point) == "value"] <- "estimate"
delta_summary_rows <- list()
for (metric_name in delta_point$metric) {
  interval <- quantile_interval(bootstrap_deltas$value[bootstrap_deltas$metric == metric_name])
  delta_summary_rows[[metric_name]] <- data.frame(
    metric = metric_name,
    estimate = delta_point$estimate[delta_point$metric == metric_name],
    lower_95 = interval[["lower_95"]],
    upper_95 = interval[["upper_95"]],
    stringsAsFactors = FALSE
  )
}
delta_summary <- do.call(rbind, delta_summary_rows)
row.names(delta_summary) <- NULL

calibration_point_rows <- list()
for (model_name in calibration_models) {
  grid <- calibration_grids[[model_name]]
  calibration_point_rows[[model_name]] <- data.frame(
    model = model_name,
    predicted_probability = grid,
    calibrated_probability = smooth_calibration(data$aad, oof_predictions[, model_name], grid),
    stringsAsFactors = FALSE
  )
}
calibration_point <- do.call(rbind, calibration_point_rows)
row.names(calibration_point) <- NULL

calibration_summary_rows <- list()
for (model_name in calibration_models) {
  grid <- calibration_grids[[model_name]]
  point <- calibration_point[calibration_point$model == model_name, , drop = FALSE]
  for (grid_index in seq_along(grid)) {
    draws <- bootstrap_calibration$calibrated_probability[
      bootstrap_calibration$model == model_name &
        abs(bootstrap_calibration$predicted_probability - grid[grid_index]) < 1e-12
    ]
    draws <- draws[is.finite(draws)]
    interval <- if (length(draws) >= 800L) quantile_interval(draws) else c(lower_95 = NA_real_, upper_95 = NA_real_)
    calibration_summary_rows[[length(calibration_summary_rows) + 1]] <- data.frame(
      model = model_name,
      predicted_probability = grid[grid_index],
      calibrated_probability = point$calibrated_probability[grid_index],
      lower_95 = interval[["lower_95"]],
      upper_95 = interval[["upper_95"]],
      successful_bootstrap_curves = length(draws),
      stringsAsFactors = FALSE
    )
  }
}
calibration_summary <- do.call(rbind, calibration_summary_rows)

probability_update <- data.frame(
  row_uid = data$row_uid,
  source_id = data$source_id,
  cluster_id = data$cluster_id,
  group = data$group,
  aad = data$aad,
  outcome = data$outcome_label,
  benchmark_probability = oof_predictions[, "ADD_RS_plus_D_dimer"],
  full_model_probability = oof_predictions[, "ADD_RS_plus_D_dimer_plus_SDC1"],
  probability_change =
    oof_predictions[, "ADD_RS_plus_D_dimer_plus_SDC1"] -
    oof_predictions[, "ADD_RS_plus_D_dimer"],
  stringsAsFactors = FALSE
)
probability_update$correct_direction <- ifelse(
  probability_update$aad == 1,
  probability_update$probability_change > 0,
  probability_update$probability_change < 0
)

direction_point <- aggregate(correct_direction ~ outcome, data = probability_update, FUN = mean)
names(direction_point)[2] <- "estimate"
direction_summary_rows <- list()
for (outcome_name in direction_point$outcome) {
  interval <- quantile_interval(
    bootstrap_direction$proportion[bootstrap_direction$outcome == outcome_name]
  )
  direction_summary_rows[[outcome_name]] <- data.frame(
    outcome = outcome_name,
    estimate = direction_point$estimate[direction_point$outcome == outcome_name],
    lower_95 = interval[["lower_95"]],
    upper_95 = interval[["upper_95"]],
    stringsAsFactors = FALSE
  )
}
direction_summary <- do.call(rbind, direction_summary_rows)
row.names(direction_summary) <- NULL

apparent_fits <- lapply(model_formulas, fit_model, training_data = data)
null_fit <- fit_model(aad ~ 1, data)
benchmark_fit <- apparent_fits[["ADD_RS_plus_D_dimer"]]
full_fit <- apparent_fits[["ADD_RS_plus_D_dimer_plus_SDC1"]]
lr_total <- 2 * (as.numeric(logLik(full_fit)) - as.numeric(logLik(null_fit)))
lr_benchmark <- 2 * (as.numeric(logLik(benchmark_fit)) - as.numeric(logLik(null_fit)))
lr_increment <- 2 * (as.numeric(logLik(full_fit)) - as.numeric(logLik(benchmark_fit)))
lr_df <- attr(logLik(full_fit), "df") - attr(logLik(benchmark_fit), "df")
lr_p <- pchisq(lr_increment, df = lr_df, lower.tail = FALSE)
lr_information <- data.frame(
  component = c("ADD-RS + D-dimer information", "New information from SDC1"),
  lr_chisq = c(lr_benchmark, lr_increment),
  fraction_of_full_lr = c(lr_benchmark / lr_total, lr_increment / lr_total),
  stringsAsFactors = FALSE
)

analyse_joint_population <- function(
    dataset,
    formulas,
    population_label,
    seed_offset,
    n_repeats = 50L,
    n_bootstrap = 1000L) {
  selected_models <- names(formulas)
  repeated_rows <- vector("list", n_repeats)

  message("Running joint-model analysis for ", population_label, "...")
  for (repeat_index in seq_len(n_repeats)) {
    predictions <- cross_validated_predictions(
      dataset,
      formulas,
      k = 10L,
      seed = seed_offset + repeat_index
    )$predictions
    model_rows <- lapply(selected_models, function(model_name) {
      metrics <- calculate_metrics(dataset$aad, predictions[, model_name])
      data.frame(
        population = population_label,
        repeat_id = repeat_index,
        model = model_name,
        metric = names(metrics),
        value = as.numeric(metrics),
        stringsAsFactors = FALSE
      )
    })
    repeated_rows[[repeat_index]] <- do.call(rbind, model_rows)
  }
  repeated <- do.call(rbind, repeated_rows)
  point <- aggregate(value ~ population + model + metric, repeated, mean)
  names(point)[names(point) == "value"] <- "estimate"

  case_indices_local <- which(dataset$aad == 1)
  control_indices_local <- which(dataset$aad == 0)
  bootstrap_rows <- vector("list", n_bootstrap)
  bootstrap_delta_rows <- vector("list", n_bootstrap)
  comparator_models <- setdiff(selected_models, "ADD_RS_plus_D_dimer_plus_SDC1")

  for (bootstrap_index in seq_len(n_bootstrap)) {
    set.seed(seed_offset + 100000L + bootstrap_index)
    sampled_indices <- c(
      sample(case_indices_local, length(case_indices_local), replace = TRUE),
      sample(control_indices_local, length(control_indices_local), replace = TRUE)
    )
    bootstrap_data <- dataset[sampled_indices, , drop = FALSE]
    predictions <- cross_validated_predictions(
      bootstrap_data,
      formulas,
      k = 10L,
      seed = seed_offset + 200000L + bootstrap_index
    )$predictions

    metric_values <- lapply(selected_models, function(model_name) {
      calculate_metrics(bootstrap_data$aad, predictions[, model_name])
    })
    names(metric_values) <- selected_models
    bootstrap_rows[[bootstrap_index]] <- do.call(rbind, lapply(selected_models, function(model_name) {
      metrics <- metric_values[[model_name]]
      data.frame(
        population = population_label,
        bootstrap_replicate = bootstrap_index,
        model = model_name,
        metric = names(metrics),
        value = as.numeric(metrics),
        stringsAsFactors = FALSE
      )
    }))

    full_metrics_local <- metric_values[["ADD_RS_plus_D_dimer_plus_SDC1"]]
    bootstrap_delta_rows[[bootstrap_index]] <- do.call(rbind, lapply(comparator_models, function(comparator_name) {
      comparator_metrics <- metric_values[[comparator_name]]
      data.frame(
        population = population_label,
        bootstrap_replicate = bootstrap_index,
        comparator = comparator_name,
        metric = c("Brier improvement", "Log-loss improvement", "C-index improvement"),
        value = c(
          comparator_metrics[["Brier"]] - full_metrics_local[["Brier"]],
          comparator_metrics[["Log_loss"]] - full_metrics_local[["Log_loss"]],
          full_metrics_local[["C_index"]] - comparator_metrics[["C_index"]]
        ),
        stringsAsFactors = FALSE
      )
    }))

    if (bootstrap_index %% 100L == 0L) {
      message(population_label, ": completed bootstrap ", bootstrap_index, "/", n_bootstrap)
    }
  }

  bootstrap_metrics_local <- do.call(rbind, bootstrap_rows)
  bootstrap_deltas_local <- do.call(rbind, bootstrap_delta_rows)
  summary_rows <- vector("list", nrow(point))
  for (row_index in seq_len(nrow(point))) {
    point_row <- point[row_index, , drop = FALSE]
    draws <- bootstrap_metrics_local$value[
      bootstrap_metrics_local$model == point_row$model &
        bootstrap_metrics_local$metric == point_row$metric
    ]
    interval <- quantile_interval(draws)
    summary_rows[[row_index]] <- data.frame(
      population = population_label,
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

  repeated_delta_rows <- list()
  for (repeat_index in seq_len(n_repeats)) {
    repeat_data <- repeated[repeated$repeat_id == repeat_index, , drop = FALSE]
    full_metrics_local <- setNames(
      repeat_data$value[repeat_data$model == "ADD_RS_plus_D_dimer_plus_SDC1"],
      repeat_data$metric[repeat_data$model == "ADD_RS_plus_D_dimer_plus_SDC1"]
    )
    for (comparator_name in comparator_models) {
      comparator_metrics <- setNames(
        repeat_data$value[repeat_data$model == comparator_name],
        repeat_data$metric[repeat_data$model == comparator_name]
      )
      repeated_delta_rows[[length(repeated_delta_rows) + 1L]] <- data.frame(
        repeat_id = repeat_index,
        comparator = comparator_name,
        metric = c("Brier improvement", "Log-loss improvement", "C-index improvement"),
        value = c(
          comparator_metrics[["Brier"]] - full_metrics_local[["Brier"]],
          comparator_metrics[["Log_loss"]] - full_metrics_local[["Log_loss"]],
          full_metrics_local[["C_index"]] - comparator_metrics[["C_index"]]
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  repeated_deltas_local <- do.call(rbind, repeated_delta_rows)
  delta_point_local <- aggregate(value ~ comparator + metric, repeated_deltas_local, mean)
  names(delta_point_local)[names(delta_point_local) == "value"] <- "estimate"
  delta_summary_local <- do.call(rbind, lapply(seq_len(nrow(delta_point_local)), function(row_index) {
    point_row <- delta_point_local[row_index, , drop = FALSE]
    draws <- bootstrap_deltas_local$value[
      bootstrap_deltas_local$comparator == point_row$comparator &
        bootstrap_deltas_local$metric == point_row$metric
    ]
    interval <- quantile_interval(draws)
    data.frame(
      population = population_label,
      comparator = point_row$comparator,
      metric = point_row$metric,
      estimate = point_row$estimate,
      lower_95 = interval[["lower_95"]],
      upper_95 = interval[["upper_95"]],
      stringsAsFactors = FALSE
    )
  }))

  list(
    performance = do.call(rbind, summary_rows),
    paired_differences = delta_summary_local,
    repeated_metrics = repeated,
    bootstrap_metrics = bootstrap_metrics_local,
    bootstrap_differences = bootstrap_deltas_local
  )
}

aad_vs_pe_data <- data[data$aad == 1 | data$group == "PE", , drop = FALSE]
if (nrow(aad_vs_pe_data) != 142L || sum(aad_vs_pe_data$aad == 1) != 99L || sum(aad_vs_pe_data$aad == 0) != 43L) {
  stop("Unexpected AAD-versus-PE counts")
}
aad_vs_pe_models <- model_formulas[c(
  "ADD_RS_plus_D_dimer",
  "ADD_RS_plus_SDC1",
  "ADD_RS_plus_D_dimer_plus_SDC1"
)]
aad_vs_pe_analysis <- analyse_joint_population(
  aad_vs_pe_data,
  aad_vs_pe_models,
  population_label = "AAD versus PE",
  seed_offset = 20267000L,
  n_repeats = 50L,
  n_bootstrap = 1000L
)
aad_vs_pe_performance <- aad_vs_pe_analysis$performance
aad_vs_pe_differences <- aad_vs_pe_analysis$paired_differences

write.csv(patient_predictions, file.path(source_dir, "fig5_revised_patient_level_repeated_oof_predictions.csv"), row.names = FALSE, na = "")
write.csv(oof_long, file.path(source_dir, "fig5_revised_repeated_cv_prediction_draws.csv"), row.names = FALSE, na = "")
write.csv(repeated_cv_metrics, file.path(source_dir, "fig5_revised_repeated_cv_metric_draws.csv"), row.names = FALSE, na = "")
write.csv(repeated_cv_deltas, file.path(source_dir, "fig5_revised_repeated_cv_paired_delta_draws.csv"), row.names = FALSE, na = "")
write.csv(bootstrap_metrics, file.path(source_dir, "fig5_revised_full_pipeline_bootstrap_metric_draws.csv"), row.names = FALSE, na = "")
write.csv(bootstrap_deltas, file.path(source_dir, "fig5_revised_full_pipeline_bootstrap_paired_delta_draws.csv"), row.names = FALSE, na = "")
write.csv(bootstrap_direction, file.path(source_dir, "fig5_revised_full_pipeline_bootstrap_direction_draws.csv"), row.names = FALSE, na = "")
write.csv(bootstrap_calibration, file.path(source_dir, "fig5_revised_full_pipeline_bootstrap_calibration_draws.csv"), row.names = FALSE, na = "")
write.csv(calibration_summary, file.path(source_dir, "fig5_revised_calibration_curve_source_data.csv"), row.names = FALSE, na = "")
write.csv(probability_update, file.path(source_dir, "fig5_revised_probability_updating_source_data.csv"), row.names = FALSE, na = "")
write.csv(aad_vs_pe_analysis$repeated_metrics, file.path(source_dir, "fig5_revised_aad_vs_pe_repeated_cv_metric_draws.csv"), row.names = FALSE, na = "")
write.csv(aad_vs_pe_analysis$bootstrap_metrics, file.path(source_dir, "fig5_revised_aad_vs_pe_full_pipeline_bootstrap_metric_draws.csv"), row.names = FALSE, na = "")
write.csv(aad_vs_pe_analysis$bootstrap_differences, file.path(source_dir, "fig5_revised_aad_vs_pe_full_pipeline_bootstrap_difference_draws.csv"), row.names = FALSE, na = "")
write.csv(performance_summary, file.path(result_dir, "fig5_revised_validated_performance.csv"), row.names = FALSE, na = "")
write.csv(delta_summary, file.path(result_dir, "fig5_revised_paired_incremental_performance.csv"), row.names = FALSE, na = "")
write.csv(direction_summary, file.path(result_dir, "fig5_revised_correct_direction_updating.csv"), row.names = FALSE, na = "")
write.csv(lr_information, file.path(result_dir, "fig5_revised_likelihood_ratio_information.csv"), row.names = FALSE, na = "")
write.csv(data.frame(lr_increment = lr_increment, df = lr_df, p_value = lr_p), file.path(result_dir, "fig5_revised_sdc1_added_value_lr_test.csv"), row.names = FALSE, na = "")
write.csv(aad_vs_pe_performance, file.path(result_dir, "fig5_revised_aad_vs_pe_validated_performance.csv"), row.names = FALSE, na = "")
write.csv(aad_vs_pe_differences, file.path(result_dir, "fig5_revised_aad_vs_pe_paired_differences.csv"), row.names = FALSE, na = "")

interval_alignment <- transform(
  performance_summary,
  endpoint_gap = pmax(lower_95 - estimate, estimate - upper_95, 0)
)
write.csv(
  interval_alignment,
  file.path(qa_dir, "FIGURE5_REVISED_POINT_INTERVAL_ALIGNMENT_AUDIT.csv"),
  row.names = FALSE,
  na = ""
)

model_factor <- factor(model_labels[performance_summary$model], levels = rev(unname(model_labels[model_order])))
performance_summary$model_label <- model_factor
model_short_labels <- c(
  ADD_RS = "ADD-RS",
  ADD_RS_plus_D_dimer = "+ D-dimer",
  ADD_RS_plus_SDC1 = "+ SDC1",
  ADD_RS_plus_D_dimer_plus_SDC1 = "+ D-dimer + SDC1"
)
performance_summary$model_short <- factor(
  model_short_labels[performance_summary$model],
  levels = rev(unname(model_short_labels[model_order]))
)
performance_summary$metric_label <- factor(
  performance_summary$metric,
  levels = c("Brier", "Log_loss", "C_index"),
  labels = c("Brier score", "Log loss", "C-index")
)

theme_pub <- function() {
  theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 6.5),
      axis.text = element_text(size = 6.0, colour = "#222222"),
      legend.title = element_text(size = 5.8, face = "bold"),
      legend.text = element_text(size = 5.6),
      strip.text = element_text(size = 6.2, face = "bold"),
      plot.title = element_text(size = 7, face = "bold"),
      plot.subtitle = element_text(size = 5.8, colour = "#687078"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

percent_label <- function(x) paste0(round(100 * x), "%")

calibration_summary$model_label <- factor(
  model_labels[calibration_summary$model],
  levels = unname(model_labels[calibration_models])
)

p_a <- ggplot(
  calibration_summary,
  aes(predicted_probability, calibrated_probability, colour = model_label, fill = model_label)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, colour = "#8F969C") +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.13, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = model_palette, drop = FALSE) +
  scale_fill_manual(values = model_palette, drop = FALSE) +
  scale_x_continuous(labels = percent_label, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = 0)) +
  scale_y_continuous(labels = percent_label, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = 0)) +
  labs(
    x = "Internally validated predicted probability",
    y = "Observed AAD proportion",
    colour = NULL,
    fill = NULL,
    title = "Smooth internal calibration",
    subtitle = "Repeated OOF curves with 95% full-pipeline bootstrap intervals"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), fill = "none") +
  coord_equal()

performance_plot_data <- performance_summary[
  performance_summary$metric %in% c("Brier", "Log_loss", "C_index"),
  ,
  drop = FALSE
]
p_b <- ggplot(performance_plot_data, aes(estimate, model_short, colour = model_label)) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.0) +
  facet_wrap(~metric_label, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = model_palette, guide = "none") +
  labs(
    x = "Full-pipeline validated estimate (95% CI)",
    y = NULL,
    title = "Validated predictive performance",
    subtitle = "Lower Brier/log loss and higher C-index are better"
  ) +
  theme_pub() +
  theme(strip.background = element_rect(fill = "#F1F3F4", colour = NA), panel.spacing.y = grid::unit(1.4, "mm"))

calibration_scalar <- performance_summary[
  performance_summary$metric %in% c("Calibration_intercept", "Calibration_slope") &
    performance_summary$model %in% calibration_models,
  ,
  drop = FALSE
]
calibration_scalar$metric_label <- factor(
  calibration_scalar$metric,
  levels = c("Calibration_intercept", "Calibration_slope"),
  labels = c("Intercept (target 0)", "Slope (target 1)")
)
calibration_targets <- data.frame(
  metric_label = factor(c("Intercept (target 0)", "Slope (target 1)"), levels = levels(calibration_scalar$metric_label)),
  target = c(0, 1)
)
p_c <- ggplot(calibration_scalar, aes(estimate, model_short, colour = model_label)) +
  geom_vline(data = calibration_targets, aes(xintercept = target), linetype = 2, colour = "#8F969C", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 1.9) +
  facet_wrap(~metric_label, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = model_palette, guide = "none") +
  labs(x = "Validated estimate (95% CI)", y = NULL, title = "Calibration intercept and slope") +
  theme_pub() +
  theme(strip.background = element_rect(fill = "#F1F3F4", colour = NA), panel.spacing.y = grid::unit(1.5, "mm"))

delta_summary$metric_label <- factor(
  delta_summary$metric,
  levels = c("Brier improvement", "Log-loss improvement", "C-index improvement")
)
p_d <- ggplot(delta_summary, aes(estimate, metric_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#8F969C", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0, linewidth = 0.6, colour = "#2A9D8F") +
  geom_point(size = 2.2, colour = "#2A9D8F") +
  facet_wrap(~metric_label, scales = "free_x", ncol = 1) +
  labs(
    x = "Improvement after adding SDC1 (95% CI)",
    y = NULL,
    title = "Paired incremental performance",
    subtitle = "Full model versus ADD-RS + D-dimer; positive favours full model"
  ) +
  theme_pub() +
  theme(
    strip.background = element_rect(fill = "#F1F3F4", colour = NA),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.y = grid::unit(1.3, "mm")
  )

lr_information$component <- factor(lr_information$component, levels = rev(lr_information$component))
p_e <- ggplot(lr_information, aes(fraction_of_full_lr, component, fill = component)) +
  geom_col(width = 0.56) +
  geom_text(
    aes(label = paste0(round(100 * fraction_of_full_lr, 1), "%")),
    hjust = 1.08,
    colour = "white",
    fontface = "bold",
    size = 2.15
  ) +
  scale_fill_manual(values = c(
    "ADD-RS + D-dimer information" = "#A8AFB5",
    "New information from SDC1" = "#BE4652"
  ), guide = "none") +
  scale_x_continuous(labels = percent_label, limits = c(0, 0.80), expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "Fraction of full-model LR information",
    y = NULL,
    title = "Relative new information",
    subtitle = sprintf("SDC1 increment: LR chi-square %.1f (%d df)", lr_increment, lr_df)
  ) +
  theme_pub()

p_f <- ggplot(probability_update, aes(benchmark_probability, full_model_probability, colour = outcome)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#8F969C", linewidth = 0.4) +
  geom_point(size = 1.05, alpha = 0.62) +
  scale_colour_manual(values = outcome_palette) +
  scale_x_continuous(labels = percent_label, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = 0)) +
  scale_y_continuous(labels = percent_label, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = 0)) +
  labs(
    x = "ADD-RS + D-dimer probability",
    y = "Full-model probability",
    colour = NULL,
    title = "Patient probability updating",
    subtitle = "Continuous estimates; no diagnostic cutoff imposed"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  coord_equal()

p_g <- ggplot(probability_update, aes(probability_change, colour = outcome, fill = outcome)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#8F969C", linewidth = 0.4) +
  geom_density(alpha = 0.18, linewidth = 0.7, adjust = 0.9) +
  scale_colour_manual(values = outcome_palette) +
  scale_fill_manual(values = outcome_palette) +
  scale_x_continuous(labels = percent_label) +
  labs(
    x = "Probability change after adding SDC1",
    y = "Density",
    colour = NULL,
    fill = NULL,
    title = "Change after adding SDC1"
  ) +
  theme_pub() +
  guides(colour = "none", fill = "none") +
  theme(legend.position = "none")

direction_summary$outcome <- factor(direction_summary$outcome, levels = c("AAD", "Symptomatic non-AAD"))
p_h <- ggplot(direction_summary, aes(outcome, estimate, fill = outcome)) +
  geom_col(width = 0.58) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.16, linewidth = 0.5) +
  geom_text(aes(y = pmin(upper_95 + 0.04, 1.04), label = percent_label(estimate)), size = 2.2) +
  scale_fill_manual(values = outcome_palette, guide = "none") +
  scale_y_continuous(labels = percent_label, limits = c(0, 1.08), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL,
    y = "Correct-direction updates",
    title = "Directionally coherent updating",
    subtitle = "Rise for AAD; fall for symptomatic non-AAD"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

aad_vs_pe_plot_data <- aad_vs_pe_performance[
  aad_vs_pe_performance$metric %in% c("Brier", "Log_loss", "C_index"),
  ,
  drop = FALSE
]
aad_vs_pe_plot_data$model_label <- factor(
  model_labels[aad_vs_pe_plot_data$model],
  levels = rev(unname(model_labels[names(aad_vs_pe_models)]))
)
aad_vs_pe_plot_data$model_short <- factor(
  model_short_labels[aad_vs_pe_plot_data$model],
  levels = rev(unname(model_short_labels[names(aad_vs_pe_models)]))
)
aad_vs_pe_plot_data$metric_label <- factor(
  aad_vs_pe_plot_data$metric,
  levels = c("Brier", "Log_loss", "C_index"),
  labels = c("Brier score", "Log loss", "C-index")
)
p_i <- ggplot(aad_vs_pe_plot_data, aes(estimate, model_short, colour = model_label)) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0, linewidth = 0.58) +
  geom_point(size = 2.1) +
  facet_grid(. ~ metric_label, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = model_palette, guide = "none") +
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    x = "Full-pipeline validated estimate (95% CI)",
    y = NULL,
    title = "Joint-model comparison: AAD versus PE",
    subtitle = "n=142 (99 AAD; 43 PE); lower Brier/log loss and higher C-index are better"
  ) +
  theme_pub() +
  theme(
    strip.background = element_rect(fill = "#F1F3F4", colour = NA),
    panel.spacing.x = grid::unit(3.0, "mm"),
    axis.text.y = element_text(size = 5.8),
    axis.text.x = element_text(size = 5.6),
    plot.subtitle = element_text(size = 5.6)
  )

top_row <- (p_a | p_b) + plot_layout(widths = c(0.88, 1.12))
middle_row <- p_c | p_d | p_e
bottom_row <- (p_f | p_g | p_h) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")
figure5_revised <- top_row / middle_row / bottom_row / p_i +
  plot_layout(heights = c(1.08, 0.84, 0.94, 0.72)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 8, family = "Arial", margin = margin(r = 3, b = 2)))
  )

save_figure <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 300)
  print(plot)
  dev.off()
}

figure_stem <- file.path(figure_dir, "Figure_5_joint_models_revised")
save_figure(figure5_revised, figure_stem, width_mm = 183, height_mm = 238, dpi = 600)

panel_list <- list(A = p_a, B = p_b, C = p_c, D = p_d, E = p_e, F = p_f, G = p_g, H = p_h, I = p_i)
for (panel_name in names(panel_list)) {
  width_mm <- if (panel_name == "I") 183 else if (panel_name %in% c("A", "B")) 92 else 61
  height_mm <- if (panel_name == "I") 66 else if (panel_name %in% c("A", "B")) 80 else 65
  save_figure(
    panel_list[[panel_name]],
    file.path(panel_dir, paste0("Figure_5", panel_name)),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = 600
  )
}

qa_checks <- c(
  n_records_329 = nrow(data) == 329L,
  aad_vs_pe_records_142 = nrow(aad_vs_pe_data) == 142L,
  unique_patient_keys_329 = length(unique(data$cluster_id)) == 329L,
  aad_99 = sum(data$aad == 1) == 99L,
  non_aad_230 = sum(data$aad == 0) == 230L,
  no_model_variable_missingness = all(complete.cases(data[, c("aad", "add_rs", "log2_d_dimer", "log2_sdc1")])),
  add_rs_distribution_locked = identical(as.integer(table(factor(data$add_rs, levels = 1:3))), c(239L, 75L, 15L)),
  one_four_to_three_correction = sum(data$add_rs_source == 4 & data$add_rs == 3) == 1L,
  repeated_cv_predictions_finite = all(is.finite(oof_predictions)),
  maximum_point_interval_endpoint_gap_below_0_005 = max(interval_alignment$endpoint_gap) < 0.005,
  paired_delta_points_inside_intervals = all(
    delta_summary$estimate >= delta_summary$lower_95 &
      delta_summary$estimate <= delta_summary$upper_95
  ),
  aad_vs_pe_plotted_points_inside_intervals = all(
    aad_vs_pe_plot_data$estimate >= aad_vs_pe_plot_data$lower_95 &
      aad_vs_pe_plot_data$estimate <= aad_vs_pe_plot_data$upper_95
  ),
  aad_vs_pe_paired_differences_inside_intervals = all(
    aad_vs_pe_differences$estimate >= aad_vs_pe_differences$lower_95 &
      aad_vs_pe_differences$estimate <= aad_vs_pe_differences$upper_95
  ),
  bootstrap_draw_count = nrow(bootstrap_metrics) == n_bootstrap * length(model_order) * 5L,
  aad_vs_pe_bootstrap_draw_count = nrow(aad_vs_pe_analysis$bootstrap_metrics) == 1000L * length(aad_vs_pe_models) * 5L,
  calibration_curve_draws_present = nrow(bootstrap_calibration) == n_bootstrap * length(calibration_models) * 50L,
  calibration_curve_summary_finite = all(is.finite(unlist(calibration_summary[, c(
    "calibrated_probability", "lower_95", "upper_95"
  )]))),
  main_exports_present = all(file.exists(paste0(figure_stem, c(".pdf", ".svg", ".tiff", ".png")))),
  panel_exports_present = all(vapply(
    LETTERS[1:9],
    function(panel_name) all(file.exists(file.path(panel_dir, paste0("Figure_5", panel_name, c(".pdf", ".svg", ".tiff", ".png"))))),
    logical(1)
  ))
)

qa_table <- data.frame(check = names(qa_checks), pass = as.logical(qa_checks), stringsAsFactors = FALSE)
write.csv(qa_table, file.path(qa_dir, "FINAL_FIGURE5_REVISED_QA.csv"), row.names = FALSE)

summary_lines <- c(
  "Revised Figure 5: ADD-RS, D-dimer and SDC1 joint-model analysis",
  paste0("Overall computational QA: ", if (all(qa_checks)) "PASS" else "FAIL"),
  "Population: 329 symptomatic participants (99 AAD; 230 symptomatic non-AAD)",
  "AAD-versus-PE comparison: 142 participants (99 AAD; 43 PE)",
  "Model variables complete in all participants",
  "Point estimates: 50 repeats of stratified 10-fold cross-validation",
  "Intervals: 1,000 full-pipeline outcome-stratified patient bootstraps with fold reconstruction",
  paste0("SDC1 added-value LR chi-square: ", sprintf("%.4f", lr_increment), " on ", lr_df, " df; P=", format.pval(lr_p, digits = 4, eps = 1e-16)),
  "No stepwise selection, optimized threshold, Youden statistic or forced biomarker dichotomization",
  paste0("Main figure: ", figure_stem)
)
writeLines(summary_lines, file.path(log_dir, "fig5_revised_analysis_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "fig5_revised_sessionInfo.txt"), useBytes = TRUE)

cat(paste(summary_lines, collapse = "\n"), "\n")
if (!all(qa_checks)) quit(status = 1)

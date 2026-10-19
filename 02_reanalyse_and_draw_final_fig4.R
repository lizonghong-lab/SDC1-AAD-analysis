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

set.seed(20260721)

analysis_dir <- file.path(root, "02_analysis_data")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
figure_dir <- file.path(root, "04_results", "figures", "main")
log_dir <- file.path(root, "04_results", "logs")
for (directory in c(result_dir, source_dir, figure_dir, log_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
clip_probability <- function(p) pmin(pmax(as.numeric(p), 1e-8), 1 - 1e-8)

hc3_vcov <- function(fit) {
  x <- model.matrix(fit)
  residual <- residuals(fit)
  leverage <- pmin(hatvalues(fit), 1 - 1e-8)
  omega <- (residual / (1 - leverage))^2
  bread <- solve(crossprod(x))
  meat <- crossprod(x, sweep(x, 1, omega, "*"))
  bread %*% meat %*% bread
}

contrast_result <- function(fit, vcov_hc3, contrast) {
  estimate <- as.numeric(crossprod(contrast, coef(fit)))
  variance <- as.numeric(t(contrast) %*% vcov_hc3 %*% contrast)
  standard_error <- sqrt(max(variance, 0))
  statistic <- estimate / standard_error
  data.frame(
    estimate_log2 = estimate,
    robust_se = standard_error,
    lower_log2_95 = estimate - qnorm(0.975) * standard_error,
    upper_log2_95 = estimate + qnorm(0.975) * standard_error,
    p_value = 2 * pnorm(abs(statistic), lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

adjusted_predictions <- function(fit, vcov_hc3, newdata) {
  x_new <- model.matrix(delete.response(terms(fit)), newdata, contrasts.arg = fit$contrasts)
  estimate <- as.numeric(x_new %*% coef(fit))
  standard_error <- sqrt(pmax(diag(x_new %*% vcov_hc3 %*% t(x_new)), 0))
  data.frame(
    estimate_log2 = estimate,
    robust_se = standard_error,
    lower_log2_95 = estimate - qnorm(0.975) * standard_error,
    upper_log2_95 = estimate + qnorm(0.975) * standard_error,
    stringsAsFactors = FALSE
  )
}

format_probability <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.0001) return("<0.0001")
  if (value < 0.001) return(sprintf("%.4f", value))
  sprintf("%.3f", value)
}

scale_covariates <- function(data, include_wbc = FALSE, include_crp = FALSE) {
  data$z_age <- as.numeric(scale(data$age_years))
  data$z_creatinine <- as.numeric(scale(data$creatinine_umol_l))
  if (include_wbc) data$z_wbc <- as.numeric(scale(data$wbc_10e9_l))
  if (include_crp) data$z_crp <- as.numeric(scale(data$crp_mg_l))
  data
}

validation <- read.csv(file.path(analysis_dir, "clinical_evaluation_all.csv"), stringsAsFactors = FALSE, check.names = FALSE)
numeric_columns <- c(
  "total_sdc1", "d_dimer_mg_l", "age_years", "female", "hypertension",
  "diabetes", "smoking", "creatinine_umol_l", "crp_mg_l", "wbc_10e9_l"
)
validation[numeric_columns] <- lapply(validation[numeric_columns], as.numeric)

all_group_levels <- c("Healthy control", "AAD", "AMI", "Angina", "PE")
symptomatic_levels <- c("AAD", "AMI", "Angina", "PE")
validation$group <- factor(validation$group, levels = all_group_levels)

# Panel A: five-group total-plasma SDC1 with a common covariate set.
a_required <- c("total_sdc1", "group", "age_years", "female", "hypertension", "diabetes", "smoking", "creatinine_umol_l", "wbc_10e9_l")
a_keep <- complete.cases(validation[, a_required]) & validation$total_sdc1 > 0
a_data <- droplevels(validation[a_keep, , drop = FALSE])
a_data <- scale_covariates(a_data, include_wbc = TRUE)
a_data$log2_total_sdc1 <- log2(a_data$total_sdc1)
a_fit <- lm(log2_total_sdc1 ~ group + z_age + female + hypertension + diabetes + smoking + z_creatinine + z_wbc, data = a_data)
a_vcov <- hc3_vcov(a_fit)
a_common <- data.frame(
  group = factor(all_group_levels, levels = all_group_levels),
  z_age = 0,
  female = mean(a_data$female),
  hypertension = mean(a_data$hypertension),
  diabetes = mean(a_data$diabetes),
  smoking = mean(a_data$smoking),
  z_creatinine = 0,
  z_wbc = 0
)
a_means <- cbind(a_common["group"], adjusted_predictions(a_fit, a_vcov, a_common))
a_means$n_group <- as.integer(table(a_data$group)[as.character(a_means$group)])
a_means$estimate_display <- 2^a_means$estimate_log2
a_means$lower_display_95 <- 2^a_means$lower_log2_95
a_means$upper_display_95 <- 2^a_means$upper_log2_95
a_x <- model.matrix(delete.response(terms(a_fit)), a_common, contrasts.arg = a_fit$contrasts)
rownames(a_x) <- as.character(a_common$group)
a_contrasts <- do.call(rbind, lapply(c("Healthy control", "AMI", "Angina", "PE"), function(reference) {
  result <- contrast_result(a_fit, a_vcov, a_x["AAD", ] - a_x[reference, ])
  result$contrast <- paste("AAD versus", reference)
  result$reference <- reference
  result$geometric_mean_ratio <- 2^result$estimate_log2
  result$ratio_lower_95 <- 2^result$lower_log2_95
  result$ratio_upper_95 <- 2^result$upper_log2_95
  result
}))
a_contrasts$bh_q <- p.adjust(a_contrasts$p_value, method = "BH")
a_contrasts$n_complete <- nrow(a_data)
a_contrasts$adjustment_set <- "Age, sex, hypertension, diabetes, smoking, creatinine and WBC"

# Panel B: symptomatic D-dimer with the fixed clinical adjustment set.
b_required <- c("d_dimer_mg_l", "group", "age_years", "female", "hypertension", "diabetes", "smoking", "creatinine_umol_l", "crp_mg_l")
b_keep <- validation$group %in% symptomatic_levels & complete.cases(validation[, b_required]) & validation$d_dimer_mg_l > 0
b_data <- droplevels(validation[b_keep, , drop = FALSE])
b_data$group <- factor(b_data$group, levels = symptomatic_levels)
b_data <- scale_covariates(b_data, include_crp = TRUE)
b_data$log2_d_dimer <- log2(b_data$d_dimer_mg_l)
b_fit <- lm(log2_d_dimer ~ group + z_age + female + hypertension + diabetes + smoking + z_creatinine + z_crp, data = b_data)
b_vcov <- hc3_vcov(b_fit)
b_common <- data.frame(
  group = factor(symptomatic_levels, levels = symptomatic_levels),
  z_age = 0,
  female = mean(b_data$female),
  hypertension = mean(b_data$hypertension),
  diabetes = mean(b_data$diabetes),
  smoking = mean(b_data$smoking),
  z_creatinine = 0,
  z_crp = 0
)
b_means <- cbind(b_common["group"], adjusted_predictions(b_fit, b_vcov, b_common))
b_means$n_group <- as.integer(table(b_data$group)[as.character(b_means$group)])
b_x <- model.matrix(delete.response(terms(b_fit)), b_common, contrasts.arg = b_fit$contrasts)
rownames(b_x) <- as.character(b_common$group)
b_contrasts <- do.call(rbind, lapply(c("AMI", "Angina", "PE"), function(reference) {
  result <- contrast_result(b_fit, b_vcov, b_x["AAD", ] - b_x[reference, ])
  result$contrast <- paste("AAD versus", reference)
  result$reference <- reference
  result$geometric_mean_ratio <- 2^result$estimate_log2
  result$ratio_lower_95 <- 2^result$lower_log2_95
  result$ratio_upper_95 <- 2^result$upper_log2_95
  result
}))
b_contrasts$bh_q <- p.adjust(b_contrasts$p_value, method = "BH")
b_contrasts$n_complete <- nrow(b_data)
b_contrasts$adjustment_set <- "Age, sex, hypertension, diabetes, smoking, creatinine and CRP"

# Panels C-F: within-group hypertension effects on total-plasma SDC1.
htn_means <- list()
htn_effects <- list()
htn_plot_data <- list()
for (group_name in symptomatic_levels) {
  subset <- validation[validation$group == group_name, , drop = FALSE]
  required <- c("total_sdc1", "hypertension", "age_years", "female", "diabetes", "smoking", "creatinine_umol_l", "crp_mg_l")
  keep <- complete.cases(subset[, required]) & subset$total_sdc1 > 0
  data <- subset[keep, , drop = FALSE]
  data <- scale_covariates(data, include_crp = TRUE)
  data$log2_total_sdc1 <- log2(data$total_sdc1)
  fit <- lm(log2_total_sdc1 ~ hypertension + z_age + female + diabetes + smoking + z_creatinine + z_crp, data = data)
  vcov_hc3 <- hc3_vcov(fit)
  common <- data.frame(
    hypertension = c(0, 1),
    z_age = 0,
    female = mean(data$female),
    diabetes = mean(data$diabetes),
    smoking = mean(data$smoking),
    z_creatinine = 0,
    z_crp = 0
  )
  prediction <- adjusted_predictions(fit, vcov_hc3, common)
  prediction$hypertension <- common$hypertension
  prediction$hypertension_label <- factor(ifelse(common$hypertension == 1, "Hypertension", "No hypertension"), levels = c("No hypertension", "Hypertension"))
  prediction$group <- group_name
  prediction$n_group <- c(sum(data$hypertension == 0), sum(data$hypertension == 1))
  prediction$estimate_display <- 2^prediction$estimate_log2
  prediction$lower_display_95 <- 2^prediction$lower_log2_95
  prediction$upper_display_95 <- 2^prediction$upper_log2_95
  coefficient <- rep(0, length(coef(fit)))
  names(coefficient) <- names(coef(fit))
  coefficient["hypertension"] <- 1
  effect <- contrast_result(fit, vcov_hc3, coefficient)
  effect$group <- group_name
  effect$n_complete <- nrow(data)
  effect$n_no_hypertension <- sum(data$hypertension == 0)
  effect$n_hypertension <- sum(data$hypertension == 1)
  effect$geometric_mean_ratio <- 2^effect$estimate_log2
  effect$ratio_lower_95 <- 2^effect$lower_log2_95
  effect$ratio_upper_95 <- 2^effect$upper_log2_95
  effect$adjustment_set <- "Age, sex, diabetes, smoking, creatinine and CRP"
  data$hypertension_label <- factor(ifelse(data$hypertension == 1, "Hypertension", "No hypertension"), levels = c("No hypertension", "Hypertension"))
  htn_means[[group_name]] <- prediction
  htn_effects[[group_name]] <- effect
  htn_plot_data[[group_name]] <- data[, c("row_uid", "source_id", "group", "hypertension", "hypertension_label", "total_sdc1")]
}
htn_means <- do.call(rbind, htn_means)
htn_effects <- do.call(rbind, htn_effects)
htn_effects$bh_q <- p.adjust(htn_effects$p_value, method = "BH")
htn_plot_data <- do.call(rbind, htn_plot_data)

# Reviewer-aligned clinical prediction models for panels G-H.
auc_rank <- function(y, p) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

score_predictions <- function(y, p) {
  p <- clip_probability(p)
  c(
    brier = mean((y - p)^2),
    log_loss = -mean(y * log(p) + (1 - y) * log1p(-p)),
    c_index = auc_rank(y, p)
  )
}

cluster_outcome <- function(data) {
  split_y <- split(data$aad, data$cluster_id)
  inconsistent <- names(split_y)[vapply(split_y, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(inconsistent)) stop("Cluster IDs span outcomes")
  vapply(split_y, function(x) x[1], numeric(1))
}

make_cluster_folds <- function(data, k = 10L, seed = 1L) {
  y_by_cluster <- cluster_outcome(data)
  fold_by_cluster <- integer(length(y_by_cluster))
  names(fold_by_cluster) <- names(y_by_cluster)
  set.seed(seed)
  for (level in sort(unique(y_by_cluster))) {
    ids <- sample(names(y_by_cluster)[y_by_cluster == level])
    fold_by_cluster[ids] <- rep(seq_len(k), length.out = length(ids))
  }
  unname(fold_by_cluster[data$cluster_id])
}

safe_glm_predict <- function(formula_value, train, test) {
  fit <- suppressWarnings(glm(formula_value, data = train, family = binomial()))
  clip_probability(predict(fit, newdata = test, type = "response"))
}

run_repeated_cv <- function(data, formulas, repeats = 30L, folds = 10L, seed = 20279000L) {
  prediction_sum <- matrix(0, nrow(data), length(formulas), dimnames = list(NULL, names(formulas)))
  prediction_count <- matrix(0L, nrow(data), length(formulas), dimnames = list(NULL, names(formulas)))
  for (repetition in seq_len(repeats)) {
    fold_id <- make_cluster_folds(data, folds, seed + repetition)
    for (fold in seq_len(folds)) {
      train <- fold_id != fold
      test <- !train
      for (model_name in names(formulas)) {
        prediction <- safe_glm_predict(formulas[[model_name]], data[train, ], data[test, ])
        prediction_sum[test, model_name] <- prediction_sum[test, model_name] + prediction
        prediction_count[test, model_name] <- prediction_count[test, model_name] + 1L
      }
    }
  }
  prediction_sum / pmax(prediction_count, 1L)
}

bootstrap_cluster_indices <- function(data) {
  y_by_cluster <- cluster_outcome(data)
  selected <- unlist(lapply(sort(unique(y_by_cluster)), function(level) {
    ids <- names(y_by_cluster)[y_by_cluster == level]
    sample(ids, length(ids), replace = TRUE)
  }), use.names = FALSE)
  unlist(lapply(selected, function(id) which(data$cluster_id == id)), use.names = FALSE)
}

run_single_cv <- function(data, formulas, folds = 10L, seed = 1L) {
  fold_id <- make_cluster_folds(data, folds, seed)
  predictions <- matrix(NA_real_, nrow(data), length(formulas), dimnames = list(NULL, names(formulas)))
  for (fold in seq_len(folds)) {
    train <- fold_id != fold
    test <- !train
    for (model_name in names(formulas)) predictions[test, model_name] <- safe_glm_predict(formulas[[model_name]], data[train, ], data[test, ])
  }
  predictions
}

smooth_calibration <- function(y, p, grid) {
  fit <- lowess(p, y, f = 2 / 3, iter = 0)
  approx(fit$x, fit$y, xout = grid, rule = 2, ties = mean)$y
}

full_pipeline_bootstrap <- function(data, formulas, point_predictions, B = 500L, folds = 10L, seed = 20280000L) {
  model_names <- names(formulas)
  metrics <- c("brier", "log_loss", "c_index")
  raw <- array(NA_real_, dim = c(B, length(model_names), length(metrics)), dimnames = list(NULL, model_names, metrics))
  grids <- setNames(lapply(model_names, function(model_name) {
    p <- point_predictions[, model_name]
    seq(max(0.01, quantile(p, 0.01)), min(0.99, quantile(p, 0.99)), length.out = 100)
  }), model_names)
  calibration_draws <- setNames(lapply(model_names, function(model_name) matrix(NA_real_, B, length(grids[[model_name]]))), model_names)
  set.seed(seed)
  bootstrap_seeds <- sample.int(.Machine$integer.max, B)
  for (b in seq_len(B)) {
    set.seed(bootstrap_seeds[b])
    index <- bootstrap_cluster_indices(data)
    bootstrap_data <- data[index, , drop = FALSE]
    bootstrap_predictions <- run_single_cv(bootstrap_data, formulas, folds, bootstrap_seeds[b] + 1L)
    for (model_name in model_names) {
      raw[b, model_name, ] <- score_predictions(bootstrap_data$aad, bootstrap_predictions[, model_name])
      calibration_draws[[model_name]][b, ] <- smooth_calibration(bootstrap_data$aad, bootstrap_predictions[, model_name], grids[[model_name]])
    }
    if (b %% 50L == 0L) cat("Adjusted full-pipeline bootstrap", b, "of", B, "\n")
  }
  performance <- do.call(rbind, lapply(model_names, function(model_name) {
    observed <- score_predictions(data$aad, point_predictions[, model_name])
    do.call(rbind, lapply(metrics, function(metric) data.frame(
      model = model_name,
      metric = metric,
      estimate = observed[metric],
      lower_95 = quantile(raw[, model_name, metric], 0.025, na.rm = TRUE),
      upper_95 = quantile(raw[, model_name, metric], 0.975, na.rm = TRUE),
      validation = "30 repeated patient-key-separated 10-fold CV; 500 full-pipeline stratified author-resolved patient-key bootstrap resamples",
      stringsAsFactors = FALSE
    )))
  }))
  calibration <- do.call(rbind, lapply(model_names, function(model_name) data.frame(
    model = model_name,
    predicted_probability = grids[[model_name]],
    cross_validated_smooth_observed = smooth_calibration(data$aad, point_predictions[, model_name], grids[[model_name]]),
    lower_95 = apply(calibration_draws[[model_name]], 2, quantile, 0.025, na.rm = TRUE),
    upper_95 = apply(calibration_draws[[model_name]], 2, quantile, 0.975, na.rm = TRUE),
    method = "LOWESS of repeated out-of-fold predictions; full-pipeline bootstrap bands",
    stringsAsFactors = FALSE
  )))
  list(raw = raw, performance = performance, calibration = calibration)
}

symptomatic <- read.csv(file.path(analysis_dir, "elisa_evaluation_primary.csv"), stringsAsFactors = FALSE, check.names = FALSE)
symptomatic$aad <- as.integer(symptomatic$aad)
symptomatic$log2_sdc1 <- log2(symptomatic$sdc1_source_value)
symptomatic$log2_d_dimer <- log2(symptomatic$d_dimer_mg_l)
symptomatic$cluster_id <- ifelse(is.na(symptomatic$cluster_id) | !nzchar(symptomatic$cluster_id), symptomatic$row_uid, symptomatic$cluster_id)
adjustment_vars <- c("aad", "log2_sdc1", "log2_d_dimer", "age_years", "female", "hypertension", "diabetes", "smoking", "creatinine_umol_l", "crp_mg_l", "cluster_id", "row_uid")
adjusted_data <- symptomatic[complete.cases(symptomatic[, adjustment_vars]), , drop = FALSE]
if (nrow(adjusted_data) != 321L || sum(adjusted_data$aad) != 98L) stop("Unexpected adjusted modelling set")
adjusted_formulas <- list(
  Clinical_context = aad ~ age_years + female + hypertension + diabetes + smoking + creatinine_umol_l + crp_mg_l,
  Clinical_plus_D_dimer = aad ~ age_years + female + hypertension + diabetes + smoking + creatinine_umol_l + crp_mg_l + log2_d_dimer,
  Clinical_plus_D_dimer_plus_SDC1 = aad ~ age_years + female + hypertension + diabetes + smoking + creatinine_umol_l + crp_mg_l + log2_d_dimer + log2_sdc1
)
adjusted_predictions <- run_repeated_cv(adjusted_data, adjusted_formulas, repeats = 30L, folds = 10L, seed = 20279000L)
adjusted_boot <- full_pipeline_bootstrap(adjusted_data, adjusted_formulas, adjusted_predictions, B = 500L, folds = 10L, seed = 20280000L)
adjusted_boot$performance$n_complete <- nrow(adjusted_data)
adjusted_boot$performance$aad_cases <- sum(adjusted_data$aad)
adjusted_boot$performance$adjustment_set <- "Age, sex, hypertension, diabetes, smoking, creatinine and CRP; fixed a priori"

paired_difference <- do.call(rbind, lapply(c("brier", "log_loss", "c_index"), function(metric) {
  point <- score_predictions(adjusted_data$aad, adjusted_predictions[, "Clinical_plus_D_dimer_plus_SDC1"])[metric] -
    score_predictions(adjusted_data$aad, adjusted_predictions[, "Clinical_plus_D_dimer"])[metric]
  raw_difference <- adjusted_boot$raw[, "Clinical_plus_D_dimer_plus_SDC1", metric] - adjusted_boot$raw[, "Clinical_plus_D_dimer", metric]
  data.frame(
    comparison = "Clinical + D-dimer + SDC1 minus Clinical + D-dimer",
    metric = metric,
    difference = point,
    lower_95 = quantile(raw_difference, 0.025, na.rm = TRUE),
    upper_95 = quantile(raw_difference, 0.975, na.rm = TRUE),
    favourable_direction = ifelse(metric == "c_index", "positive", "negative"),
    stringsAsFactors = FALSE
  )
}))

prediction_export <- data.frame(
  row_uid = adjusted_data$row_uid,
  source_id = adjusted_data$source_id,
  cluster_id = adjusted_data$cluster_id,
  group = adjusted_data$group,
  aad = adjusted_data$aad,
  adjusted_predictions,
  stringsAsFactors = FALSE
)

metric_draws <- data.frame(
  bootstrap_replicate = rep(seq_len(dim(adjusted_boot$raw)[1]), times = dim(adjusted_boot$raw)[2] * dim(adjusted_boot$raw)[3]),
  model = rep(rep(dimnames(adjusted_boot$raw)[[2]], each = dim(adjusted_boot$raw)[1]), times = dim(adjusted_boot$raw)[3]),
  metric = rep(dimnames(adjusted_boot$raw)[[3]], each = dim(adjusted_boot$raw)[1] * dim(adjusted_boot$raw)[2]),
  value = as.vector(adjusted_boot$raw),
  stringsAsFactors = FALSE
)

write_csv(a_means, file.path(result_dir, "fig4_adjusted_total_sdc1_group_means.csv"))
write_csv(a_contrasts, file.path(result_dir, "fig4_adjusted_total_sdc1_group_contrasts.csv"))
write_csv(b_means, file.path(result_dir, "fig4_adjusted_d_dimer_group_means.csv"))
write_csv(b_contrasts, file.path(result_dir, "fig4_adjusted_d_dimer_group_contrasts.csv"))
write_csv(htn_means, file.path(result_dir, "fig4_adjusted_hypertension_means.csv"))
write_csv(htn_effects, file.path(result_dir, "fig4_adjusted_hypertension_effects.csv"))
write_csv(adjusted_boot$performance, file.path(result_dir, "fig4_clinically_adjusted_validated_performance.csv"))
write_csv(paired_difference, file.path(result_dir, "fig4_clinically_adjusted_paired_differences.csv"))
write_csv(adjusted_boot$calibration, file.path(source_dir, "fig4_clinically_adjusted_calibration.csv"))
write_csv(prediction_export, file.path(source_dir, "fig4_clinically_adjusted_oof_predictions.csv"))
write_csv(metric_draws, file.path(source_dir, "fig4_clinically_adjusted_bootstrap_metric_draws.csv"))
write_csv(a_data[, c("row_uid", "source_id", "group", "total_sdc1")], file.path(source_dir, "fig4_panel_a_source_data.csv"))
write_csv(b_data[, c("row_uid", "source_id", "group", "d_dimer_mg_l", "log2_d_dimer")], file.path(source_dir, "fig4_panel_b_source_data.csv"))
write_csv(htn_plot_data, file.path(source_dir, "fig4_panels_c_to_f_source_data.csv"))

palette <- c("Healthy control" = "#4B8F8A", "AAD" = "#BF3F4C", "AMI" = "#3F78A8", "Angina" = "#D78B2B", "PE" = "#7B5AA6")
model_palette <- c("Clinical" = "#9AA0A6", "+ D-dimer" = "#3F78A8", "+ D-dimer + SDC1" = "#BF3F4C")
htn_palette <- c("No hypertension" = "#B8BDC2", "Hypertension" = "#5C6470")

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6.1, colour = "black"),
      plot.title = element_text(size = 8, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 5.9, colour = "#555555", margin = margin(b = 3)),
      strip.text = element_text(size = 6.2, face = "bold"),
      legend.text = element_text(size = 5.7),
      plot.margin = margin(6, 5, 5, 5)
    )
}

# A: raw values and covariate-adjusted geometric means.
a_labels <- setNames(paste0(c("HC", "AAD", "AMI", "Angina", "PE"), "\n(n=", a_means$n_group, ")"), all_group_levels)
a_annotation <- merge(data.frame(reference = c("Healthy control", "AMI", "Angina", "PE"), stringsAsFactors = FALSE), a_contrasts[, c("reference", "bh_q")], by = "reference")
a_annotation$x <- match(a_annotation$reference, all_group_levels)
a_annotation$y <- max(a_data$total_sdc1) * 1.36
a_annotation$label <- paste0("vs AAD\nq", vapply(a_annotation$bh_q, format_probability, character(1)))
p_a <- ggplot(a_data, aes(group, total_sdc1, colour = group)) +
  geom_jitter(width = 0.11, size = 0.60, alpha = 0.42) +
  geom_errorbar(data = a_means, aes(x = group, ymin = lower_display_95, ymax = upper_display_95), inherit.aes = FALSE, width = 0.14, linewidth = 0.55, colour = "#202428") +
  geom_point(data = a_means, aes(x = group, y = estimate_display), inherit.aes = FALSE, shape = 21, size = 2.2, stroke = 0.6, fill = "white", colour = "#202428") +
  geom_text(data = a_annotation, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 1.9, lineheight = 0.9) +
  scale_colour_manual(values = palette) +
  scale_x_discrete(labels = a_labels) +
  scale_y_log10(expand = expansion(mult = c(0.04, 0.23))) +
  labs(x = NULL, y = "Total-plasma SDC1 (pg/mL, log scale)", title = "Total-plasma SDC1", subtitle = "Adjusted for age, sex, hypertension, diabetes, smoking, creatinine and WBC") +
  theme_pub() + theme(legend.position = "none")

# B: adjusted D-dimer in symptomatic groups.
b_labels <- setNames(paste0(symptomatic_levels, "\n(n=", b_means$n_group, ")"), symptomatic_levels)
b_annotation <- merge(data.frame(reference = c("AMI", "Angina", "PE"), stringsAsFactors = FALSE), b_contrasts[, c("reference", "bh_q")], by = "reference")
b_annotation$x <- match(b_annotation$reference, symptomatic_levels)
b_annotation$y <- max(b_data$log2_d_dimer) + 1.1
b_annotation$label <- paste0("vs AAD\nq", vapply(b_annotation$bh_q, format_probability, character(1)))
p_b <- ggplot(b_data, aes(group, log2_d_dimer, colour = group)) +
  geom_jitter(width = 0.11, size = 0.62, alpha = 0.44) +
  geom_errorbar(data = b_means, aes(x = group, ymin = lower_log2_95, ymax = upper_log2_95), inherit.aes = FALSE, width = 0.14, linewidth = 0.55, colour = "#202428") +
  geom_point(data = b_means, aes(x = group, y = estimate_log2), inherit.aes = FALSE, shape = 21, size = 2.2, stroke = 0.6, fill = "white", colour = "#202428") +
  geom_text(data = b_annotation, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 1.9, lineheight = 0.9) +
  scale_colour_manual(values = palette) +
  scale_x_discrete(labels = b_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = NULL, y = "log2 D-dimer (mg/L)", title = "D-dimer", subtitle = "Adjusted for age, sex, hypertension, diabetes, smoking, creatinine and CRP") +
  theme_pub() + theme(legend.position = "none")

add_log_bracket <- function(plot, y, label) {
  plot +
    annotate("segment", x = 1, xend = 2, y = y, yend = y, linewidth = 0.35) +
    annotate("segment", x = 1, xend = 1, y = y, yend = y / 1.08, linewidth = 0.35) +
    annotate("segment", x = 2, xend = 2, y = y, yend = y / 1.08, linewidth = 0.35) +
    annotate("text", x = 1.5, y = y * 1.04, label = label, size = 2.0, vjust = 0)
}

make_htn_panel <- function(group_name) {
  raw <- htn_plot_data[htn_plot_data$group == group_name, , drop = FALSE]
  mean_data <- htn_means[htn_means$group == group_name, , drop = FALSE]
  effect <- htn_effects[htn_effects$group == group_name, , drop = FALSE]
  labels <- setNames(paste0(c("No HTN", "HTN"), "\n(n=", mean_data$n_group, ")"), levels(raw$hypertension_label))
  bracket_y <- max(raw$total_sdc1) * 1.20
  plot <- ggplot(raw, aes(hypertension_label, total_sdc1, colour = hypertension_label)) +
    geom_jitter(width = 0.10, size = 0.63, alpha = 0.46) +
    geom_errorbar(data = mean_data, aes(x = hypertension_label, ymin = lower_display_95, ymax = upper_display_95), inherit.aes = FALSE, width = 0.14, linewidth = 0.55, colour = "#202428") +
    geom_point(data = mean_data, aes(x = hypertension_label, y = estimate_display), inherit.aes = FALSE, shape = 21, size = 2.1, stroke = 0.6, fill = "white", colour = "#202428") +
    scale_colour_manual(values = htn_palette) +
    scale_x_discrete(labels = labels) +
    scale_y_log10(expand = expansion(mult = c(0.04, 0.22))) +
    labs(x = NULL, y = "Total-plasma SDC1 (pg/mL, log scale)", title = group_name, subtitle = "Adjusted for age, sex, diabetes, smoking, creatinine and CRP") +
    theme_pub() + theme(legend.position = "none")
  add_log_bracket(plot, bracket_y, paste0("BH q", format_probability(effect$bh_q)))
}

p_c <- make_htn_panel("AAD")
p_d <- make_htn_panel("AMI")
p_e <- make_htn_panel("Angina")
p_f <- make_htn_panel("PE")

# G: validated performance beyond AUROC alone.
performance <- adjusted_boot$performance
performance$model_label <- factor(
  performance$model,
  levels = c("Clinical_context", "Clinical_plus_D_dimer", "Clinical_plus_D_dimer_plus_SDC1"),
  labels = c("Clinical", "+ D-dimer", "+ D-dimer + SDC1")
)
performance$metric_label <- factor(performance$metric, levels = c("brier", "log_loss", "c_index"), labels = c("Brier (lower)", "Log loss (lower)", "C-index / AUROC (higher)"))
p_g <- ggplot(performance, aes(estimate, model_label, colour = model_label)) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0.14, linewidth = 0.5) +
  geom_point(size = 2.2) +
  facet_grid(. ~ metric_label, scales = "free_x") +
  scale_colour_manual(values = model_palette) +
  labs(x = "Cross-validated estimate (95% full-pipeline bootstrap interval)", y = NULL, title = "Validated probability performance", subtitle = paste0("n=", nrow(adjusted_data), "; lower Brier/log loss and higher C-index are favourable")) +
  theme_pub() + theme(legend.position = "none", axis.text.y = element_text(size = 5.8))

# H: full-pipeline cross-validated calibration.
calibration <- adjusted_boot$calibration
calibration$model_label <- factor(
  calibration$model,
  levels = c("Clinical_context", "Clinical_plus_D_dimer", "Clinical_plus_D_dimer_plus_SDC1"),
  labels = c("Clinical", "+ D-dimer", "+ D-dimer + SDC1")
)
for (column in c("cross_validated_smooth_observed", "lower_95", "upper_95")) calibration[[column]] <- pmin(pmax(calibration[[column]], 0), 1)
p_h <- ggplot(calibration, aes(predicted_probability, cross_validated_smooth_observed, colour = model_label, fill = model_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#90969C", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.72) +
  scale_colour_manual(values = model_palette) +
  scale_fill_manual(values = model_palette) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  coord_equal() +
  labs(x = "Cross-validated predicted AAD probability", y = "Smooth observed AAD proportion", colour = NULL, fill = NULL, title = "Cross-validated calibration", subtitle = "30x10-fold CV; 500 full-pipeline bootstrap bands") +
  theme_pub() + theme(legend.position = "bottom")

top_row <- (p_a | p_b) + plot_layout(widths = c(1.18, 1))
middle_row <- p_c | p_d | p_e | p_f
bottom_row <- (p_g | p_h) + plot_layout(widths = c(1.25, 1))
figure <- top_row / middle_row / bottom_row +
  plot_layout(heights = c(0.92, 0.90, 1.18)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 9, family = "Arial")))

save_figure <- function(plot, stem, width_mm = 183, height_mm = 195, dpi = 600) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = cairo_pdf, family = "Arial")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, device = svglite::svglite)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300, device = ragg::agg_png)
  ggsave(paste0(stem, ".tiff"), plot, width = width, height = height, dpi = dpi, compression = "lzw", device = ragg::agg_tiff)
}

save_figure(figure, file.path(figure_dir, "Figure_4_clinically_adjusted_performance_and_calibration"))

log_lines <- c(
  paste("Panel A complete records:", nrow(a_data)),
  paste("Panel B complete records:", nrow(b_data)),
  paste("Adjusted modelling set:", nrow(adjusted_data), "records;", sum(adjusted_data$aad), "AAD"),
  paste("Full model Brier:", round(performance$estimate[performance$model == "Clinical_plus_D_dimer_plus_SDC1" & performance$metric == "brier"], 4)),
  paste("Full model log loss:", round(performance$estimate[performance$model == "Clinical_plus_D_dimer_plus_SDC1" & performance$metric == "log_loss"], 4)),
  paste("Full model C-index:", round(performance$estimate[performance$model == "Clinical_plus_D_dimer_plus_SDC1" & performance$metric == "c_index"], 4)),
  "Panels G-H include Brier, log loss, C-index/AUROC and smooth calibration; no diagnostic threshold or Youden J is used."
)
writeLines(log_lines, file.path(log_dir, "fig4_adjusted_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "fig4_sessionInfo.txt"), useBytes = TRUE)
cat(paste(log_lines, collapse = "\n"), "\n")

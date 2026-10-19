options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
legacy_rlib <- Sys.getenv(
  "SDC1_RLIB",
  unset = ""
)
if (dir.exists(legacy_rlib)) .libPaths(c(legacy_rlib, .libPaths()))

analysis_dir <- file.path(root, "02_analysis_data")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
log_dir <- file.path(root, "04_results", "logs")
for (directory in c(result_dir, source_dir, log_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
clip_probability <- function(p) pmin(pmax(as.numeric(p), 1e-8), 1 - 1e-8)

auc_rank <- function(y, p) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
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

run_single_cv <- function(data, formulas, folds = 10L, seed = 1L) {
  fold_id <- make_cluster_folds(data, folds, seed)
  predictions <- matrix(NA_real_, nrow(data), length(formulas), dimnames = list(NULL, names(formulas)))
  for (fold in seq_len(folds)) {
    train <- fold_id != fold
    test <- !train
    for (model_name in names(formulas)) {
      predictions[test, model_name] <- safe_glm_predict(formulas[[model_name]], data[train, ], data[test, ])
    }
  }
  predictions
}

bootstrap_cluster_indices <- function(data) {
  y_by_cluster <- cluster_outcome(data)
  selected <- unlist(lapply(sort(unique(y_by_cluster)), function(level) {
    ids <- names(y_by_cluster)[y_by_cluster == level]
    sample(ids, length(ids), replace = TRUE)
  }), use.names = FALSE)
  unlist(lapply(selected, function(id) which(data$cluster_id == id)), use.names = FALSE)
}

roc_curve <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(data.frame(false_positive_rate = numeric(), true_positive_rate = numeric(), threshold = numeric()))
  order_index <- order(p, decreasing = TRUE)
  y <- y[order_index]
  p <- p[order_index]
  endpoint <- c(which(diff(p) != 0), length(p))
  data.frame(
    false_positive_rate = c(0, cumsum(y == 0)[endpoint] / n0),
    true_positive_rate = c(0, cumsum(y == 1)[endpoint] / n1),
    threshold = c(Inf, p[endpoint]),
    stringsAsFactors = FALSE
  )
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

prediction_export <- read.csv(file.path(source_dir, "fig4_clinically_adjusted_oof_predictions.csv"), stringsAsFactors = FALSE, check.names = FALSE)
prediction_export <- prediction_export[match(adjusted_data$row_uid, prediction_export$row_uid), , drop = FALSE]
if (any(prediction_export$row_uid != adjusted_data$row_uid)) stop("OOF prediction rows do not match adjusted modelling data")
point_predictions <- as.matrix(prediction_export[, names(adjusted_formulas), drop = FALSE])

population_index <- list(
  "All symptomatic comparators" = rep(TRUE, nrow(adjusted_data)),
  "AAD versus PE" = adjusted_data$aad == 1 | adjusted_data$group == "PE"
)

roc_source <- do.call(rbind, lapply(names(population_index), function(population_name) {
  index <- population_index[[population_name]]
  do.call(rbind, lapply(names(adjusted_formulas), function(model_name) {
    curve <- roc_curve(adjusted_data$aad[index], point_predictions[index, model_name])
    curve$population <- population_name
    curve$model <- model_name
    curve
  }))
}))

point_auc <- do.call(rbind, lapply(names(population_index), function(population_name) {
  index <- population_index[[population_name]]
  do.call(rbind, lapply(names(adjusted_formulas), function(model_name) data.frame(
    population = population_name,
    model = model_name,
    estimate = auc_rank(adjusted_data$aad[index], point_predictions[index, model_name]),
    n_total = sum(index),
    aad_cases = sum(adjusted_data$aad[index] == 1),
    non_aad_controls = sum(adjusted_data$aad[index] == 0),
    stringsAsFactors = FALSE
  )))
}))

B <- 500L
bootstrap_draws <- vector("list", B)
set.seed(20280721)
bootstrap_seeds <- sample.int(.Machine$integer.max, B)
for (b in seq_len(B)) {
  set.seed(bootstrap_seeds[b])
  index <- bootstrap_cluster_indices(adjusted_data)
  bootstrap_data <- adjusted_data[index, , drop = FALSE]
  bootstrap_predictions <- run_single_cv(bootstrap_data, adjusted_formulas, folds = 10L, seed = bootstrap_seeds[b] + 1L)
  bootstrap_populations <- list(
    "All symptomatic comparators" = rep(TRUE, nrow(bootstrap_data)),
    "AAD versus PE" = bootstrap_data$aad == 1 | bootstrap_data$group == "PE"
  )
  bootstrap_draws[[b]] <- do.call(rbind, lapply(names(bootstrap_populations), function(population_name) {
    population_rows <- bootstrap_populations[[population_name]]
    do.call(rbind, lapply(names(adjusted_formulas), function(model_name) data.frame(
      bootstrap_replicate = b,
      population = population_name,
      model = model_name,
      auc = auc_rank(bootstrap_data$aad[population_rows], bootstrap_predictions[population_rows, model_name]),
      stringsAsFactors = FALSE
    )))
  }))
  if (b %% 50L == 0L) cat("Adjusted ROC full-pipeline bootstrap", b, "of", B, "\n")
}
bootstrap_draws <- do.call(rbind, bootstrap_draws)

auc_intervals <- do.call(rbind, lapply(seq_len(nrow(point_auc)), function(i) {
  row <- point_auc[i, , drop = FALSE]
  draws <- bootstrap_draws$auc[bootstrap_draws$population == row$population & bootstrap_draws$model == row$model]
  row$lower_95 <- quantile(draws, 0.025, na.rm = TRUE)
  row$upper_95 <- quantile(draws, 0.975, na.rm = TRUE)
  row$validation <- "30 repeated patient-key-separated 10-fold CV point estimate; 500 full-pipeline stratified author-resolved patient-key bootstrap intervals"
  row
}))

write_csv(roc_source, file.path(source_dir, "fig4_clinically_adjusted_roc_curves.csv"))
write_csv(bootstrap_draws, file.path(source_dir, "fig4_clinically_adjusted_roc_bootstrap_draws.csv"))
write_csv(auc_intervals, file.path(result_dir, "fig4_clinically_adjusted_roc_auc.csv"))

log_lines <- c(
  paste("Adjusted ROC modelling set:", nrow(adjusted_data), "records;", sum(adjusted_data$aad), "AAD"),
  paste("AAD versus PE subgroup:", sum(population_index[["AAD versus PE"]]), "records;", sum(adjusted_data$aad[population_index[["AAD versus PE"]]]), "AAD"),
  "ROC curves use repeated out-of-fold predictions; intervals use 500 full-pipeline bootstrap resamples.",
  "No ROC-derived diagnostic cut-point or Youden J statistic was selected."
)
writeLines(log_lines, file.path(log_dir, "fig4_adjusted_roc_summary.txt"), useBytes = TRUE)
cat(paste(log_lines, collapse = "\n"), "\n")

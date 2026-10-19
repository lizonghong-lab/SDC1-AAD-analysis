options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
suppressPackageStartupMessages({
  library(limma)
  library(glmnet)
})
set.seed(20260714)

analysis_dir <- file.path(root, "02_analysis_data")
results_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
logs_dir <- file.path(root, "04_results", "logs")
for (directory in c(results_dir, source_dir, logs_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")

raw <- readRDS(file.path(analysis_dir, "proteomics_raw_deidentified.rds"))
qc <- readRDS(file.path(analysis_dir, "proteomics_qc_reduction.rds"))
meta <- raw$sample_metadata
log2_raw <- log2(raw$intensity)
rownames(log2_raw) <- raw$protein_group_id
sdc1_id <- raw$protein_group_id[toupper(raw$gene_symbol) == "SDC1"]
if (length(sdc1_id) != 1L) stop("Expected one SDC1 protein group")

auc_rank <- function(y, p) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

score_predictions <- function(y, p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  c(
    brier = mean((y - p)^2),
    log_loss = -mean(y * log(p) + (1 - y) * log1p(-p)),
    c_index = auc_rank(y, p)
  )
}

make_folds <- function(y, k = 10L, seed = 1L) {
  set.seed(seed)
  fold <- integer(length(y))
  for (level in sort(unique(y))) {
    index <- sample(which(y == level))
    fold[index] <- rep(seq_len(k), length.out = length(index))
  }
  fold
}

# Joint ridge hurdle model: detection and conditional abundance enter together.
ridge_hurdle <- function(reference_group, seed) {
  keep <- meta$group %in% c("AAD", reference_group)
  d <- meta[keep, , drop = FALSE]
  values <- log2_raw[sdc1_id, keep]
  detected <- as.integer(is.finite(values))
  centre <- median(values[is.finite(values)])
  abundance_component <- ifelse(is.finite(values), values - centre, 0)
  y <- as.integer(d$group == "AAD")
  clinical <- model.matrix(~ age_years + female + hypertension, d)[, -1, drop = FALSE]
  hurdle <- cbind(clinical, sdc1_detected = detected, conditional_abundance = abundance_component)
  outer_fold <- make_folds(y, 10L, seed)
  predictions <- matrix(NA_real_, nrow = length(y), ncol = 2, dimnames = list(NULL, c("Clinical", "Clinical_plus_two_part_SDC1")))
  matrices <- list(Clinical = clinical, Clinical_plus_two_part_SDC1 = hurdle)

  for (fold in seq_len(10L)) {
    train <- outer_fold != fold
    test <- !train
    for (model_name in names(matrices)) {
      x <- matrices[[model_name]]
      inner_fold <- make_folds(y[train], 5L, seed + 100L * fold + match(model_name, names(matrices)))
      fitted <- cv.glmnet(
        x[train, , drop = FALSE], y[train], family = "binomial", alpha = 0,
        foldid = inner_fold, type.measure = "deviance", standardize = TRUE
      )
      predictions[test, model_name] <- as.numeric(predict(fitted, x[test, , drop = FALSE], s = "lambda.1se", type = "response"))
    }
  }

  scores <- lapply(colnames(predictions), function(model_name) {
    metric <- score_predictions(y, predictions[, model_name])
    data.frame(
      contrast = paste("AAD versus", reference_group),
      model = model_name,
      n = length(y),
      aad_records = sum(y),
      reference_records = sum(!y),
      brier = metric["brier"],
      log_loss = metric["log_loss"],
      c_index = metric["c_index"],
      validation = "Nested 10-fold cross-validation; ridge penalty selected inside each training fold",
      stringsAsFactors = FALSE
    )
  })
  score_table <- do.call(rbind, scores)
  delta <- score_table[2, ]
  delta$model <- "Increment from two-part SDC1"
  delta$brier <- score_table$brier[2] - score_table$brier[1]
  delta$log_loss <- score_table$log_loss[2] - score_table$log_loss[1]
  delta$c_index <- score_table$c_index[2] - score_table$c_index[1]
  rbind(score_table, delta)
}

hurdle_performance <- rbind(
  ridge_hurdle("AMI", 20264001),
  ridge_hurdle("Healthy control", 20264002)
)
write_csv(hurdle_performance, file.path(results_dir, "proteomics_two_part_ridge_hurdle_v4.csv"))

prepare_design <- function(metadata, positive_group, reference_group) {
  columns <- which(metadata$group %in% c(positive_group, reference_group))
  d <- metadata[columns, , drop = FALSE]
  complete <- complete.cases(d[, c("group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l")])
  columns <- columns[complete]
  d <- d[complete, , drop = FALSE]
  d$aad_indicator <- as.integer(d$group == positive_group)
  d$age_scaled <- as.numeric(scale(d$age_years))
  d$wbc_scaled <- as.numeric(scale(log1p(d$wbc_10e9_l)))
  d$crp_scaled <- as.numeric(scale(log1p(d$crp_mg_l)))
  design <- model.matrix(~ aad_indicator + age_scaled + female + hypertension + wbc_scaled + crp_scaled, d)
  list(columns = columns, metadata = d, design = design)
}

fit_sdc1 <- function(matrix_value, protein_ids, reference_group) {
  if (!sdc1_id %in% protein_ids) return(c(logFC = NA, t = NA, p = NA, rank = NA))
  prepared <- prepare_design(meta, "AAD", reference_group)
  matrix_pair <- matrix_value[, prepared$columns, drop = FALSE]
  fit <- eBayes(lmFit(matrix_pair, prepared$design), trend = TRUE, robust = TRUE)
  coef_index <- which(colnames(prepared$design) == "aad_indicator")
  table <- topTable(fit, coef = coef_index, number = Inf, sort.by = "none")
  row <- match(sdc1_id, rownames(table))
  c(
    logFC = table$logFC[row],
    t = table$t[row],
    p = table$P.Value[row],
    rank = rank(-abs(table$t), ties.method = "average")[row]
  )
}

normalise_matrix <- function(matrix_value, method) {
  if (method == "none") return(matrix_value)
  if (method == "sample_median") {
    medians <- apply(matrix_value, 2, median, na.rm = TRUE)
    return(sweep(matrix_value, 2, medians - median(medians, na.rm = TRUE), "-"))
  }
  if (method == "quantile") return(normalizeBetweenArrays(matrix_value, method = "quantile"))
  if (method == "cyclic_loess") return(normalizeCyclicLoess(matrix_value, method = "fast", iterations = 2))
  stop("Unknown normalisation method")
}

grid <- expand.grid(
  detection_filter = c(0.40, 0.50, 0.60, 0.70),
  normalisation = c("none", "sample_median", "quantile", "cyclic_loess"),
  reference_group = c("AMI", "Healthy control"),
  stringsAsFactors = FALSE
)
grid_results <- lapply(seq_len(nrow(grid)), function(i) {
  eligible <- rowMeans(is.finite(log2_raw)) >= grid$detection_filter[i]
  ids <- raw$protein_group_id[eligible]
  matrix_value <- normalise_matrix(log2_raw[eligible, , drop = FALSE], grid$normalisation[i])
  rownames(matrix_value) <- ids
  estimate <- fit_sdc1(matrix_value, ids, grid$reference_group[i])
  data.frame(
    grid[i, ],
    retained_protein_groups = sum(eligible),
    sdc1_eligible = sdc1_id %in% ids,
    sdc1_overall_detection_fraction = mean(is.finite(log2_raw[sdc1_id, ])),
    log2_fold_change = estimate["logFC"],
    moderated_t = estimate["t"],
    p_value = estimate["p"],
    absolute_t_rank = estimate["rank"],
    stringsAsFactors = FALSE
  )
})
preprocessing_grid <- do.call(rbind, grid_results)
write_csv(preprocessing_grid, file.path(results_dir, "proteomics_preprocessing_stability_grid_v4.csv"))

# Left-censoring imputation is deliberately a sensitivity analysis, not a primary estimator.
left_censor_impute <- function(matrix_value, seed = 20265000) {
  set.seed(seed)
  result <- matrix_value
  for (i in seq_len(nrow(result))) {
    observed <- result[i, is.finite(result[i, ])]
    if (!length(observed)) next
    sigma <- sd(observed)
    if (!is.finite(sigma) || sigma == 0) sigma <- 0.25
    missing <- which(!is.finite(result[i, ]))
    if (length(missing)) {
      result[i, missing] <- rnorm(length(missing), min(observed) - 1.8 * sigma, 0.3 * sigma)
    }
  }
  result
}

primary_ids <- qc$primary_protein_id
primary_raw <- log2_raw[match(primary_ids, rownames(log2_raw)), , drop = FALSE]
rownames(primary_raw) <- primary_ids
imputed <- left_censor_impute(primary_raw)
imputed <- normalise_matrix(imputed, "sample_median")
conditional <- qc$log2_median_normalised
rownames(conditional) <- primary_ids

missingness_rows <- lapply(c("AMI", "Healthy control"), function(reference_group) {
  conditional_result <- fit_sdc1(conditional, primary_ids, reference_group)
  imputed_result <- fit_sdc1(imputed, primary_ids, reference_group)
  detection_model <- read.csv(file.path(results_dir, "sdc1_proteomic_detection_models.csv"), stringsAsFactors = FALSE)
  detection_row <- detection_model[grepl(reference_group, detection_model$contrast, fixed = TRUE), ][1, ]
  rbind(
    data.frame(
      contrast = paste("AAD versus", reference_group),
      approach = "Detection component",
      estimate = detection_row$adjusted_detection_odds_ratio,
      estimate_scale = "Adjusted detection odds ratio",
      p_value = detection_row$p_value,
      direction = "Higher detection in AAD",
      stringsAsFactors = FALSE
    ),
    data.frame(
      contrast = paste("AAD versus", reference_group),
      approach = "Conditional abundance without imputation",
      estimate = conditional_result["logFC"],
      estimate_scale = "Adjusted log2 fold change among detected values",
      p_value = conditional_result["p"],
      direction = ifelse(conditional_result["logFC"] > 0, "Higher in AAD", "Lower in AAD"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      contrast = paste("AAD versus", reference_group),
      approach = "Left-censoring imputation sensitivity",
      estimate = imputed_result["logFC"],
      estimate_scale = "Adjusted log2 fold change after low-value imputation",
      p_value = imputed_result["p"],
      direction = ifelse(imputed_result["logFC"] > 0, "Higher in AAD", "Lower in AAD"),
      stringsAsFactors = FALSE
    )
  )
})
missingness_sensitivity <- do.call(rbind, missingness_rows)
write_csv(missingness_sensitivity, file.path(results_dir, "sdc1_missingness_mnar_sensitivity_v4.csv"))

# Rank-based GO enrichment keeps the complete ranked list and avoids an arbitrary DE cutoff.
run_gsea <- function(contrast_label) {
  differential <- read.csv(file.path(source_dir, "proteomics_adjusted_effects_all.csv"), stringsAsFactors = FALSE)
  differential <- differential[differential$contrast == contrast_label & differential$normalisation == "sample_median", ]
  differential <- differential[is.finite(differential$t) & !is.na(differential$gene_symbol) & nzchar(differential$gene_symbol), ]
  split_rows <- split(seq_len(nrow(differential)), differential$gene_symbol)
  best <- vapply(split_rows, function(index) index[which.max(abs(differential$t[index]))], integer(1))
  ranks <- differential$t[best]
  names(ranks) <- differential$gene_symbol[best]
  ranks <- sort(ranks, decreasing = TRUE)

  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(data.frame(contrast = contrast_label, status = "Required offline annotation package unavailable"))
  }
  enrichment <- suppressMessages(clusterProfiler::gseGO(
    geneList = ranks,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    eps = 0,
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  ))
  table <- as.data.frame(enrichment)
  if (!nrow(table)) return(data.frame(contrast = contrast_label, status = "No GO terms returned"))
  table$contrast <- contrast_label
  table$status <- "Rank-based GSEA; all eligible proteins retained"
  table
}

gsea_results <- rbind(
  run_gsea("AAD versus AMI"),
  run_gsea("AAD versus Healthy control")
)
write_csv(gsea_results, file.path(source_dir, "proteomics_rank_based_GO_GSEA_v4.csv"))

if ("p.adjust" %in% names(gsea_results)) {
  gsea_summary <- do.call(rbind, lapply(split(gsea_results, gsea_results$contrast), function(x) {
    x <- x[order(x$p.adjust, -abs(x$NES)), ]
    head(x, 15)
  }))
} else {
  gsea_summary <- gsea_results
}
write_csv(gsea_summary, file.path(results_dir, "proteomics_rank_based_GO_GSEA_top_v4.csv"))

writeLines(capture.output(sessionInfo()), file.path(logs_dir, "02b_sessionInfo.txt"))
cat("Advanced proteomic analyses complete: hurdle, processing grid, MNAR sensitivity and rank-based GO GSEA.\n")

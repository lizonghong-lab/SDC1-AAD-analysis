options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1L) as.integer(args[[1L]]) else 1000L
run_label <- if (length(args) >= 2L) args[[2L]] else "FINAL_STANDALONE_B1000_LOWESS_ABCLOC_LOCKED"
if (!is.finite(B) || B < 10L) stop("B must be at least 10")

root_dir <- normalizePath(getwd(), mustWork = TRUE)
project_dir <- normalizePath(
  file.path(root_dir, "SDC1_AAD_HARRELL_BOOTCAL_20260728"),
  mustWork = TRUE
)
.libPaths(c(file.path(project_dir, "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(rms)
  library(Hmisc)
})

input_file <- file.path(
  root_dir,
  "SDC1_AAD_REANALYSIS_FINAL_20260721",
  "02_analysis_data",
  "elisa_evaluation_primary.csv"
)
stat_dir <- file.path(project_dir, "04_results", "statistical_results", run_label)
source_dir <- file.path(project_dir, "04_results", "source_data", run_label)
log_dir <- file.path(project_dir, "04_results", "logs", run_label)
qc_dir <- file.path(project_dir, "07_quality_control", run_label)
for (d in c(stat_dir, source_dir, log_dir, qc_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "standalone_harrell_bootcal.log")
sink(log_file, split = TRUE)
on.exit({
  while (sink.number() > 0L) sink()
}, add = TRUE)

raw <- read.csv(input_file, check.names = FALSE)
required <- c(
  "group", "aad", "symptomatic", "log2_d_dimer",
  "log2_total_sdc1", "row_uid"
)
missing_columns <- setdiff(required, names(raw))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

eligible <- with(
  raw,
  symptomatic == 1 & !is.na(aad) &
    !is.na(log2_d_dimer) & !is.na(log2_total_sdc1)
)
d_sym <- raw[eligible, , drop = FALSE]
d_sym$aad <- as.integer(d_sym$aad)
d_pe <- d_sym[d_sym$group %in% c("AAD", "PE"), , drop = FALSE]

populations <- list(symptomatic = d_sym, aad_vs_pe = d_pe)
population_labels <- c(
  symptomatic = "AAD versus symptomatic non-AAD",
  aad_vs_pe = "AAD versus pulmonary embolism"
)
model_formulas <- list(
  D_dimer = aad ~ rcs(log2_d_dimer, 3),
  SDC1 = aad ~ rcs(log2_total_sdc1, 3)
)
model_labels <- c(D_dimer = "D-dimer", SDC1 = "SDC1")
population_seeds <- c(symptomatic = 202607285L, aad_vs_pe = 202607286L)

metric_measure <- function(
    x, y, fit, iter, evalfit = FALSE, pr = FALSE,
    penalty.matrix = NULL, kint = 1L, ...) {
  k <- fit$non.slopes
  null_model <- length(fit$coefficients) == k
  if (evalfit) {
    dxy <- unname(fit$stats["Dxy"])
    intercept <- if (null_model) NA_real_ else 0
    slope <- if (null_model) NA_real_ else 1
  } else {
    refit <- if (null_model) lrm.fit(y = y) else lrm.fit(x, y)
    kr <- refit$non.slopes
    dxy <- unname(refit$stats["Dxy"])
    intercept <- if (null_model) NA_real_ else unname(refit$coefficients[kint])
    slope <- if (null_model) NA_real_ else unname(refit$coefficients[kr + 1L])
  }
  probability <- pmin(pmax(plogis(x), 1e-15), 1 - 1e-15)
  outcome <- as.integer(y >= kint)
  c(
    Dxy = dxy,
    Intercept = intercept,
    Slope = slope,
    Brier = mean((outcome - probability)^2),
    LogLoss = -mean(
      outcome * log(probability) +
        (1 - outcome) * log(1 - probability)
    )
  )
}

lrm_refit <- function(
    x, y, penalty.matrix = NULL, xcol = NULL,
    strata, iter, ...) {
  if (length(xcol) && length(penalty.matrix) > 0L) {
    penalty.matrix <- penalty.matrix[xcol, xcol, drop = FALSE]
  }
  lrm.fit(x, y, penalty.matrix = penalty.matrix, ...)
}

validate_extended <- function(fit, B, group, seed) {
  fit_work <- fit
  y <- fit_work$y
  if (!is.factor(y)) y <- factor(y)
  fit_work$y <- unclass(y) - 1L
  set.seed(seed)
  z <- rms:::predab.resample(
    fit.orig = fit_work,
    fit = lrm_refit,
    measure = metric_measure,
    method = "boot",
    B = B,
    bw = FALSE,
    estimates = TRUE,
    non.slopes.in.x = FALSE,
    penalty.matrix = fit_work$penalty.matrix,
    kint = 1L,
    group = group,
    saveraw = TRUE
  )
  raw_draws <- get(".predab_raw.", envir = .GlobalEnv)
  colnames(z) <- c(
    "index.orig", "training", "test", "optimism",
    "index.corrected", "Lower", "Upper", "n"
  )[seq_len(ncol(z))]
  rownames(z) <- c("Dxy", "Intercept", "Slope", "Brier", "LogLoss")
  list(summary = z, raw = raw_draws)
}

extract_scalar <- function(validation, population, model) {
  z <- validation$summary
  pick <- function(row_name, metric, transform = identity) {
    data.frame(
      population = population,
      model = model,
      model_label = unname(model_labels[model]),
      metric = metric,
      apparent = transform(z[row_name, "index.orig"]),
      optimism = if (row_name == "Dxy") z[row_name, "optimism"] / 2 else z[row_name, "optimism"],
      corrected = transform(z[row_name, "index.corrected"]),
      lower = transform(z[row_name, "Lower"]),
      upper = transform(z[row_name, "Upper"]),
      successful_resamples = z[row_name, "n"],
      stringsAsFactors = FALSE
    )
  }
  rbind(
    pick("Dxy", "C-index", function(x) (x + 1) / 2),
    pick("Brier", "Brier score"),
    pick("LogLoss", "Log loss"),
    pick("Intercept", "Calibration intercept"),
    pick("Slope", "Calibration slope")
  )
}

raw_metric <- function(raw_object, metric, component) {
  raw_columns <- names(raw_object$orig)
  if (metric == "C-index") {
    return((raw_object[[component]][, match("Dxy", raw_columns)] + 1) / 2)
  }
  column <- switch(
    metric,
    "Brier score" = "Brier",
    "Log loss" = "LogLoss",
    stop("Unknown metric: ", metric)
  )
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
    apparent = orig,
    optimism = mean(optimism_draw),
    corrected = corrected,
    lower = corrected - z * unname(directional_sd["top"]),
    upper = corrected + z * unname(directional_sd["bottom"]),
    successful_resamples = length(btrain)
  )
}

fit_population <- function(data, population_key) {
  cat("Population", population_key, "n=", nrow(data), "events=", sum(data$aad), "\n")
  dd_local <- datadist(data[, c("aad", "log2_d_dimer", "log2_total_sdc1")])
  assign("dd_standalone", dd_local, envir = .GlobalEnv)
  options(datadist = "dd_standalone")

  fits <- lapply(model_formulas, function(f) lrm(f, data = data, x = TRUE, y = TRUE))
  validations <- list()
  calibrations <- list()
  scalar_rows <- list()
  calibration_rows <- list()
  base_seed <- unname(population_seeds[population_key])

  for (model in names(fits)) {
    cat("  validate", model, "\n")
    validations[[model]] <- validate_extended(
      fits[[model]], B = B, group = data$aad, seed = base_seed
    )
    scalar_rows[[model]] <- extract_scalar(validations[[model]], population_key, model)

    p <- plogis(predict(fits[[model]], type = "lp"))
    q <- seq(0.05, 0.95, length.out = 49L)
    grid <- unname(quantile(p, probs = q, type = 8, na.rm = TRUE))
    keep <- !duplicated(grid)
    grid <- grid[keep]
    q <- q[keep]
    cat("  calibrate", model, "\n")
    set.seed(base_seed + 10000L)
    cal <- calibrate(
      fits[[model]], method = "boot", B = B,
      predy = grid, smoother = "lowess", group = data$aad
    )
    calibrations[[model]] <- cal
    cal_df <- as.data.frame(unclass(cal))
    calibration_rows[[model]] <- data.frame(
      population = population_key,
      population_label = unname(population_labels[population_key]),
      model = model,
      model_label = unname(model_labels[model]),
      prediction_quantile = q,
      predicted_probability = cal_df$predy,
      apparent_observed = cal_df$calibrated.orig,
      corrected_observed = cal_df$calibrated.corrected,
      lower_error_scale = cal_df$Lower,
      upper_error_scale = cal_df$Upper,
      lower_observed = cal_df$predy + cal_df$Lower,
      upper_observed = cal_df$predy + cal_df$Upper,
      successful_resamples = cal_df$n,
      stringsAsFactors = FALSE
    )
  }

  metrics <- c("C-index", "Brier score", "Log loss")
  difference_rows <- lapply(metrics, function(metric) {
    s <- validations$SDC1$raw
    d <- validations$D_dimer$raw
    if (metric == "C-index") {
      ans <- abcloc_summary(
        orig_metric(s, metric) - orig_metric(d, metric),
        raw_metric(s, metric, "btrain") - raw_metric(d, metric, "btrain"),
        raw_metric(s, metric, "btest") - raw_metric(d, metric, "btest")
      )
    } else {
      ans <- abcloc_summary(
        orig_metric(d, metric) - orig_metric(s, metric),
        raw_metric(d, metric, "btrain") - raw_metric(s, metric, "btrain"),
        raw_metric(d, metric, "btest") - raw_metric(s, metric, "btest")
      )
    }
    data.frame(
      population = population_key,
      comparison = "SDC1_vs_Ddimer",
      metric = metric,
      positive_value_favours = "SDC1",
      corrected_improvement = ans["corrected"],
      lower = ans["lower"],
      upper = ans["upper"],
      successful_resamples = ans["successful_resamples"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  list(
    fits = fits,
    validations = validations,
    calibrations = calibrations,
    scalar = do.call(rbind, scalar_rows),
    calibration = do.call(rbind, calibration_rows),
    differences = do.call(rbind, difference_rows)
  )
}

results <- lapply(names(populations), function(pop) fit_population(populations[[pop]], pop))
names(results) <- names(populations)
scalar_all <- do.call(rbind, lapply(results, `[[`, "scalar"))
calibration_all <- do.call(rbind, lapply(results, `[[`, "calibration"))
differences_all <- do.call(rbind, lapply(results, `[[`, "differences"))
rownames(scalar_all) <- rownames(calibration_all) <- rownames(differences_all) <- NULL

write.csv(scalar_all, file.path(stat_dir, "standalone_harrell_abcloc_scalar_performance.csv"), row.names = FALSE, na = "")
write.csv(differences_all, file.path(stat_dir, "standalone_harrell_abcloc_paired_differences.csv"), row.names = FALSE, na = "")
write.csv(calibration_all, file.path(source_dir, "standalone_harrell_bootcal_calibration_source_data.csv"), row.names = FALSE, na = "")
saveRDS(results, file.path(source_dir, "standalone_harrell_bootcal_complete_results.rds"))

for (population in names(results)) {
  for (model in names(results[[population]]$validations)) {
    v <- results[[population]]$validations[[model]]$raw
    cols <- names(v$orig)
    take <- function(name, component) v[[component]][, match(name, cols)]
    out <- data.frame(
      bootstrap = seq_len(nrow(v$btrain)),
      population = population,
      model = model,
      train_C_index = (take("Dxy", "btrain") + 1) / 2,
      test_C_index = (take("Dxy", "btest") + 1) / 2,
      train_Brier = take("Brier", "btrain"),
      test_Brier = take("Brier", "btest"),
      train_Log_loss = take("LogLoss", "btrain"),
      test_Log_loss = take("LogLoss", "btest"),
      train_calibration_intercept = take("Intercept", "btrain"),
      test_calibration_intercept = take("Intercept", "btest"),
      train_calibration_slope = take("Slope", "btrain"),
      test_calibration_slope = take("Slope", "btest"),
      stringsAsFactors = FALSE
    )
    write.csv(
      out,
      file.path(source_dir, paste0("bootstrap_train_test_", population, "_", model, ".csv")),
      row.names = FALSE,
      na = ""
    )
  }
}

qc <- data.frame(
  check = c(
    "symptomatic_n_329", "aad_vs_pe_n_142",
    "all_scalar_resamples_1000", "all_difference_resamples_1000",
    "scalar_intervals_ordered", "calibration_corrected_present",
    "calibration_abcloc_limits_present"
  ),
  pass = c(
    nrow(d_sym) == 329L,
    nrow(d_pe) == 142L,
    all(scalar_all$successful_resamples == B),
    all(differences_all$successful_resamples == B),
    all(scalar_all$lower <= scalar_all$corrected & scalar_all$corrected <= scalar_all$upper),
    all(is.finite(calibration_all$corrected_observed)),
    any(is.finite(calibration_all$lower_observed) & is.finite(calibration_all$upper_observed))
  ),
  stringsAsFactors = FALSE
)
write.csv(qc, file.path(qc_dir, "standalone_harrell_bootcal_qa.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

cat("Standalone strict Harrell analysis complete. QA:", if (all(qc$pass)) "PASS" else "FAIL", "\n")
if (!all(qc$pass)) quit(status = 1)

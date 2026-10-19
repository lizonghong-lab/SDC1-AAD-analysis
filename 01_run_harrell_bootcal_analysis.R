options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1L) as.integer(args[[1L]]) else 1000L
run_label <- if (length(args) >= 2L) args[[2L]] else paste0("B", B)
calibration_smoother <- if (length(args) >= 3L) args[[3L]] else "lowess"
if (!is.finite(B) || B < 10L) stop("B must be an integer of at least 10")
if (!calibration_smoother %in% c("lowess", "hare")) {
  stop("Calibration smoother must be 'lowess' or 'hare'")
}

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
  library(scales)
})

required_versions <- c(rms = "8.1.1", Hmisc = "5.2.6")
for (pkg in names(required_versions)) {
  if (packageVersion(pkg) < package_version(required_versions[[pkg]])) {
    stop(pkg, " version ", required_versions[[pkg]], " or later is required")
  }
}

input_file <- file.path(
  root_dir,
  "SDC1_AAD_REANALYSIS_FINAL_20260721",
  "02_analysis_data",
  "elisa_evaluation_primary.csv"
)
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

stat_dir <- file.path(project_dir, "04_results", "statistical_results", run_label)
source_dir <- file.path(project_dir, "04_results", "source_data", run_label)
figure_dir <- file.path(project_dir, "04_results", "figures", run_label)
log_dir <- file.path(project_dir, "04_results", "logs", run_label)
qc_dir <- file.path(project_dir, "07_quality_control", run_label)
for (d in c(stat_dir, source_dir, figure_dir, log_dir, qc_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink_file <- file.path(log_dir, "harrell_bootcal_run.log")
sink(sink_file, split = TRUE)
on.exit({
  while (sink.number() > 0L) sink()
}, add = TRUE)

cat("Harrell bootcal analysis\n")
cat("Run label:", run_label, "\n")
cat("Bootstrap resamples:", B, "\n")
cat("Calibration smoother:", calibration_smoother, "\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Input MD5:", unname(tools::md5sum(input_file)), "\n")
cat("R:", R.version.string, "\n")
cat("rms:", as.character(packageVersion("rms")), "\n")
cat("Hmisc:", as.character(packageVersion("Hmisc")), "\n\n")

raw <- read.csv(input_file, check.names = FALSE)
required_columns <- c(
  "group", "aad", "symptomatic", "add_rs_corrected",
  "log2_d_dimer", "log2_total_sdc1", "row_uid"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

eligible <- with(
  raw,
  symptomatic == 1 &
    !is.na(aad) &
    !is.na(add_rs_corrected) &
    !is.na(log2_d_dimer) &
    !is.na(log2_total_sdc1)
)
d_sym <- raw[eligible, , drop = FALSE]
d_sym$aad <- as.integer(d_sym$aad)
d_pe <- d_sym[d_sym$group %in% c("AAD", "PE"), , drop = FALSE]

populations <- list(
  symptomatic = d_sym,
  aad_vs_pe = d_pe
)

model_formulas <- list(
  ADD_RS = aad ~ add_rs_corrected,
  ADD_RS_Ddimer = aad ~ add_rs_corrected + rcs(log2_d_dimer, 3),
  ADD_RS_SDC1 = aad ~ add_rs_corrected + rcs(log2_total_sdc1, 3),
  Full = aad ~ add_rs_corrected + rcs(log2_d_dimer, 3) +
    rcs(log2_total_sdc1, 3)
)

model_labels <- c(
  ADD_RS = "ADD-RS",
  ADD_RS_Ddimer = "ADD-RS + D-dimer",
  ADD_RS_SDC1 = "ADD-RS + SDC1",
  Full = "ADD-RS + D-dimer + SDC1"
)

model_colours <- c(
  ADD_RS = "#7F858B",
  ADD_RS_Ddimer = "#3F78A8",
  ADD_RS_SDC1 = "#D86613",
  Full = "#2A9D8F"
)

population_labels <- c(
  symptomatic = "AAD versus symptomatic non-AAD",
  aad_vs_pe = "AAD versus pulmonary embolism"
)

seed_by_population <- c(
  symptomatic = 202607281L,
  aad_vs_pe = 202607282L
)

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
    refit <- if (null_model) {
      lrm.fit(y = y)
    } else {
      lrm.fit(x, y)
    }
    kr <- refit$non.slopes
    dxy <- unname(refit$stats["Dxy"])
    intercept <- if (null_model) NA_real_ else unname(refit$coefficients[kint])
    slope <- if (null_model) NA_real_ else unname(refit$coefficients[kr + 1L])
  }
  probability <- pmin(pmax(plogis(x), 1e-15), 1 - 1e-15)
  outcome <- as.integer(y >= kint)
  brier <- mean((outcome - probability)^2)
  log_loss <- -mean(
    outcome * log(probability) +
      (1 - outcome) * log(1 - probability)
  )
  c(
    Dxy = dxy,
    Intercept = intercept,
    Slope = slope,
    Brier = brier,
    LogLoss = log_loss
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
  expected_columns <- c(
    "index.orig", "training", "test", "optimism",
    "index.corrected", "Lower", "Upper", "n"
  )
  colnames(z) <- expected_columns[seq_len(ncol(z))]
  rownames(z) <- c("Dxy", "Intercept", "Slope", "Brier", "LogLoss")
  list(summary = z, raw = raw_draws)
}

extract_scalar_summary <- function(v, population, model_key) {
  z <- v$summary
  pick <- function(row_name, label, transform = identity) {
    data.frame(
      population = population,
      model = model_key,
      model_label = unname(model_labels[model_key]),
      metric = label,
      apparent = transform(z[row_name, "index.orig"]),
      optimism = if (row_name == "Dxy") {
        z[row_name, "optimism"] / 2
      } else {
        z[row_name, "optimism"]
      },
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
    x <- raw_object[[component]][, match("Dxy", raw_columns), drop = TRUE]
    return((x + 1) / 2)
  }
  column <- switch(
    metric,
    "Brier score" = "Brier",
    "Log loss" = "LogLoss",
    "Calibration intercept" = "Intercept",
    "Calibration slope" = "Slope",
    stop("Unknown metric: ", metric)
  )
  raw_object[[component]][, match(column, raw_columns), drop = TRUE]
}

orig_metric <- function(raw_object, metric) {
  if (metric == "C-index") {
    return((unname(raw_object$orig["Dxy"]) + 1) / 2)
  }
  column <- switch(
    metric,
    "Brier score" = "Brier",
    "Log loss" = "LogLoss",
    "Calibration intercept" = "Intercept",
    "Calibration slope" = "Slope",
    stop("Unknown metric: ", metric)
  )
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

paired_difference <- function(
    raw_a, raw_b, population, comparison, model_a, model_b, metric,
    orientation = c("a_minus_b", "b_minus_a")) {
  orientation <- match.arg(orientation)
  oa <- orig_metric(raw_a, metric)
  ob <- orig_metric(raw_b, metric)
  ta <- raw_metric(raw_a, metric, "btrain")
  tb <- raw_metric(raw_b, metric, "btrain")
  wa <- raw_metric(raw_a, metric, "btest")
  wb <- raw_metric(raw_b, metric, "btest")
  if (orientation == "a_minus_b") {
    ans <- abcloc_summary(oa - ob, ta - tb, wa - wb)
  } else {
    ans <- abcloc_summary(ob - oa, tb - ta, wb - wa)
  }
  data.frame(
    population = population,
    comparison = comparison,
    model_a = model_a,
    model_b = model_b,
    metric = metric,
    positive_value_favours = model_a,
    apparent_improvement = ans["apparent"],
    optimism = ans["optimism"],
    corrected_improvement = ans["corrected"],
    lower = ans["lower"],
    upper = ans["upper"],
    successful_resamples = ans["successful_resamples"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

fit_population <- function(data, population_key, B) {
  cat("Population:", population_key, "n=", nrow(data),
      "events=", sum(data$aad), "non-events=", sum(1 - data$aad), "\n")
  dd_local <- datadist(data[, c(
    "aad", "add_rs_corrected", "log2_d_dimer", "log2_total_sdc1"
  )])
  assign("dd_bootcal", dd_local, envir = .GlobalEnv)
  options(datadist = "dd_bootcal")

  fits <- lapply(model_formulas, function(formula_object) {
    lrm(formula_object, data = data, x = TRUE, y = TRUE)
  })

  seed <- unname(seed_by_population[population_key])
  validations <- list()
  calibrations <- list()
  scalar_rows <- list()
  calibration_rows <- list()

  for (model_key in names(fits)) {
    cat("  validate:", model_key, "\n")
    validations[[model_key]] <- validate_extended(
      fits[[model_key]], B = B, group = data$aad, seed = seed
    )
    scalar_rows[[model_key]] <- extract_scalar_summary(
      validations[[model_key]], population_key, model_key
    )

    if (model_key != "ADD_RS") {
      cat("  calibrate:", model_key, "\n")
      apparent_probability <- plogis(predict(fits[[model_key]], type = "lp"))
      prediction_quantiles <- seq(0.05, 0.95, length.out = 49L)
      prediction_grid_raw <- unname(quantile(
        apparent_probability,
        probs = prediction_quantiles,
        na.rm = TRUE,
        type = 8
      ))
      keep_grid <- !duplicated(prediction_grid_raw)
      prediction_grid <- prediction_grid_raw[keep_grid]
      prediction_quantiles <- prediction_quantiles[keep_grid]
      set.seed(seed + 10000L)
      cal <- calibrate(
        fits[[model_key]],
        method = "boot",
        B = B,
        predy = prediction_grid,
        smoother = calibration_smoother,
        group = data$aad
      )
      calibrations[[model_key]] <- cal
      cal_df <- as.data.frame(unclass(cal))
      calibration_rows[[model_key]] <- data.frame(
        population = population_key,
        population_label = unname(population_labels[population_key]),
        model = model_key,
        model_label = unname(model_labels[model_key]),
        prediction_quantile = prediction_quantiles,
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
  }

  scalar <- do.call(rbind, scalar_rows)
  calibration <- do.call(rbind, calibration_rows)
  rownames(scalar) <- NULL
  rownames(calibration) <- NULL

  comparison_specs <- list(
    list("SDC1_vs_Ddimer", "ADD_RS_SDC1", "ADD_RS_Ddimer"),
    list("Add_SDC1_to_Ddimer", "Full", "ADD_RS_Ddimer"),
    list("Add_Ddimer_to_SDC1", "Full", "ADD_RS_SDC1"),
    list("Add_Ddimer_to_ADDRS", "ADD_RS_Ddimer", "ADD_RS"),
    list("Add_SDC1_to_ADDRS", "ADD_RS_SDC1", "ADD_RS")
  )
  difference_rows <- list()
  row_id <- 0L
  for (spec in comparison_specs) {
    comparison <- spec[[1L]]
    model_a <- spec[[2L]]
    model_b <- spec[[3L]]
    for (metric in c("C-index", "Brier score", "Log loss")) {
      row_id <- row_id + 1L
      orientation <- if (metric == "C-index") "a_minus_b" else "b_minus_a"
      difference_rows[[row_id]] <- paired_difference(
        validations[[model_a]]$raw,
        validations[[model_b]]$raw,
        population = population_key,
        comparison = comparison,
        model_a = model_a,
        model_b = model_b,
        metric = metric,
        orientation = orientation
      )
    }
  }
  differences <- do.call(rbind, difference_rows)
  rownames(differences) <- NULL

  list(
    fits = fits,
    validations = validations,
    calibrations = calibrations,
    scalar = scalar,
    calibration = calibration,
    differences = differences
  )
}

results <- list()
for (population_key in names(populations)) {
  results[[population_key]] <- fit_population(
    populations[[population_key]], population_key, B
  )
}

scalar_all <- do.call(rbind, lapply(results, `[[`, "scalar"))
calibration_all <- do.call(rbind, lapply(results, `[[`, "calibration"))
differences_all <- do.call(rbind, lapply(results, `[[`, "differences"))
rownames(scalar_all) <- NULL
rownames(calibration_all) <- NULL
rownames(differences_all) <- NULL

write.csv(
  scalar_all,
  file.path(stat_dir, "harrell_abcloc_scalar_performance.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  differences_all,
  file.path(stat_dir, "harrell_abcloc_paired_model_differences.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  calibration_all,
  file.path(source_dir, "harrell_bootcal_smooth_calibration_source_data.csv"),
  row.names = FALSE,
  na = ""
)

sample_flow <- do.call(rbind, lapply(names(populations), function(population_key) {
  d <- populations[[population_key]]
  data.frame(
    population = population_key,
    population_label = unname(population_labels[population_key]),
    n = nrow(d),
    aad = sum(d$aad),
    non_aad = sum(1 - d$aad),
    group_counts = paste(
      paste(names(table(d$group)), as.integer(table(d$group)), sep = "="),
      collapse = "; "
    ),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  sample_flow,
  file.path(qc_dir, "analysis_population_counts.csv"),
  row.names = FALSE
)

for (population_key in names(results)) {
  for (model_key in names(results[[population_key]]$validations)) {
    raw_v <- results[[population_key]]$validations[[model_key]]$raw
    raw_columns <- names(raw_v$orig)
    draw_column <- function(name, component) {
      raw_v[[component]][, match(name, raw_columns), drop = TRUE]
    }
    draw_df <- data.frame(
      bootstrap = seq_len(nrow(raw_v$btrain)),
      population = population_key,
      model = model_key,
      train_C_index = (draw_column("Dxy", "btrain") + 1) / 2,
      test_C_index = (draw_column("Dxy", "btest") + 1) / 2,
      train_Brier = draw_column("Brier", "btrain"),
      test_Brier = draw_column("Brier", "btest"),
      train_Log_loss = draw_column("LogLoss", "btrain"),
      test_Log_loss = draw_column("LogLoss", "btest"),
      train_calibration_intercept = draw_column("Intercept", "btrain"),
      test_calibration_intercept = draw_column("Intercept", "btest"),
      train_calibration_slope = draw_column("Slope", "btrain"),
      test_calibration_slope = draw_column("Slope", "btest"),
      stringsAsFactors = FALSE
    )
    write.csv(
      draw_df,
      file.path(
        source_dir,
        paste0("bootstrap_train_test_", population_key, "_", model_key, ".csv")
      ),
      row.names = FALSE,
      na = ""
    )
  }
}

cal_plot_data <- calibration_all[
  is.finite(calibration_all$corrected_observed) &
    is.finite(calibration_all$lower_observed) &
    is.finite(calibration_all$upper_observed),
]
cal_plot_data$lower_plot <- pmax(0, cal_plot_data$lower_observed)
cal_plot_data$upper_plot <- pmin(1, cal_plot_data$upper_observed)
cal_plot_data$model_label <- factor(
  cal_plot_data$model_label,
  levels = unname(model_labels)
)
cal_plot_data$population_label <- factor(
  cal_plot_data$population_label,
  levels = unname(population_labels)
)

theme_bootcal <- function() {
  theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      strip.background = element_rect(fill = "#F1F3F4", colour = NA),
      strip.text = element_text(size = 6.4, face = "bold"),
      panel.spacing = grid::unit(4, "mm"),
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      legend.key.width = grid::unit(9, "mm"),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6.3, colour = "#555B60"),
      plot.margin = margin(4, 5, 4, 4)
    )
}

make_calibration_figure <- function(plot_data, title, subtitle, ncol = 2L) {
  ggplot(
    plot_data,
    aes(x = predicted_probability, group = model)
  ) +
    geom_ribbon(
      aes(ymin = lower_plot, ymax = upper_plot, fill = model),
      alpha = 0.16,
      colour = NA,
      na.rm = TRUE
    ) +
    geom_abline(
      intercept = 0, slope = 1,
      linewidth = 0.45, linetype = "dashed", colour = "#8A9299"
    ) +
    geom_line(
      aes(y = apparent_observed),
      linewidth = 0.48, linetype = "dotted", colour = "#303438",
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = corrected_observed, colour = model),
      linewidth = 0.72,
      na.rm = TRUE
    ) +
    facet_wrap(~model_label, ncol = ncol) +
    scale_colour_manual(values = model_colours, guide = "none") +
    scale_fill_manual(values = model_colours, guide = "none") +
    scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      labels = label_percent(accuracy = 1), expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      labels = label_percent(accuracy = 1), expand = expansion(mult = 0.01)
    ) +
    labs(
      x = "Predicted probability of AAD",
      y = "Observed AAD proportion",
      title = title,
      subtitle = subtitle
    ) +
    coord_equal() +
    theme_bootcal()
}

save_figure <- function(plot, stem, width_mm, height_mm, dpi = 600L) {
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
    units = "in", res = dpi, compression = "lzw"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300
  )
  print(plot)
  dev.off()
}

for (population_key in names(populations)) {
  plot_data <- cal_plot_data[cal_plot_data$population == population_key, ]
  p <- make_calibration_figure(
    plot_data,
    title = paste0(
      "Harrell optimism-corrected calibration: ",
      unname(population_labels[population_key])
    ),
    subtitle = paste0(
      "n=", nrow(populations[[population_key]]),
      "; ", B,
      " outcome-stratified bootstrap resamples; ", calibration_smoother,
      " smooth ABCLOC 95% confidence bands"
    ),
    ncol = 2L
  )
  save_figure(
    p,
    file.path(figure_dir, paste0("Harrell_bootcal_", population_key, "_biomarker_models")),
    width_mm = 150,
    height_mm = 132
  )
}

biomarker_models <- c("ADD_RS_Ddimer", "ADD_RS_SDC1", "Full")
combined_data <- cal_plot_data[cal_plot_data$model %in% biomarker_models, ]
p_combined <- ggplot(
  combined_data,
  aes(x = predicted_probability, group = model)
) +
  geom_ribbon(
    aes(ymin = lower_plot, ymax = upper_plot, fill = model),
    alpha = 0.16, colour = NA, na.rm = TRUE
  ) +
  geom_abline(
    intercept = 0, slope = 1,
    linewidth = 0.45, linetype = "dashed", colour = "#8A9299"
  ) +
  geom_line(
    aes(y = apparent_observed),
    linewidth = 0.45, linetype = "dotted", colour = "#303438", na.rm = TRUE
  ) +
  geom_line(
    aes(y = corrected_observed, colour = model),
    linewidth = 0.72, na.rm = TRUE
  ) +
  facet_grid(population_label ~ model_label) +
  scale_colour_manual(values = model_colours, guide = "none") +
  scale_fill_manual(values = model_colours, guide = "none") +
  scale_x_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1), expand = expansion(mult = 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1), expand = expansion(mult = 0.01)
  ) +
  labs(
    x = "Predicted probability of AAD",
    y = "Observed AAD proportion",
    title = "Smooth optimism-corrected calibration with ABCLOC confidence bands",
    subtitle = paste0(
      B,
      " outcome-stratified bootstrap resamples; ", calibration_smoother,
      " smoother; solid, corrected; dotted, apparent; dashed, ideal"
    )
  ) +
  coord_equal() +
  theme_bootcal() +
  theme(legend.position = "none")

save_figure(
  p_combined,
  file.path(figure_dir, "Harrell_bootcal_biomarker_models_combined"),
  width_mm = 183,
  height_mm = 112
)

calibration_corrected_fraction <- mean(is.finite(calibration_all$corrected_observed))
calibration_limit_fraction <- mean(
  is.finite(calibration_all$lower_observed) &
    is.finite(calibration_all$upper_observed)
)
calibration_all$abcloc_interval_estimable <-
  is.finite(calibration_all$lower_observed) &
  is.finite(calibration_all$upper_observed)
calibration_model_coverage <- do.call(rbind, lapply(
  split(
    calibration_all,
    list(calibration_all$population, calibration_all$model),
    drop = TRUE
  ),
  function(x) {
    central <- x$prediction_quantile >= 0.25 & x$prediction_quantile <= 0.75
    data.frame(
      population = x$population[1L],
      model = x$model[1L],
      all_q05_q95_coverage = mean(x$abcloc_interval_estimable),
      central_q25_q75_coverage = mean(x$abcloc_interval_estimable[central]),
      total_grid_points = nrow(x),
      central_grid_points = sum(central),
      stringsAsFactors = FALSE
    )
  }
))
rownames(calibration_model_coverage) <- NULL
write.csv(
  calibration_model_coverage,
  file.path(qc_dir, "calibration_abcloc_interval_coverage.csv"),
  row.names = FALSE
)
qc <- data.frame(
  check = c(
    "symptomatic_n",
    "symptomatic_events",
    "aad_vs_pe_n",
    "aad_vs_pe_events",
    "all_scalar_resamples_complete",
    "all_scalar_intervals_ordered",
    "calibration_has_corrected_curve",
    "calibration_has_abcloc_limits"
  ),
  observed = c(
    nrow(d_sym),
    sum(d_sym$aad),
    nrow(d_pe),
    sum(d_pe$aad),
    all(scalar_all$successful_resamples == B),
    all(scalar_all$lower <= scalar_all$corrected &
          scalar_all$corrected <= scalar_all$upper),
    calibration_corrected_fraction >= 0.98,
    all(calibration_model_coverage$all_q05_q95_coverage >= 0.60) &&
      all(calibration_model_coverage$central_q25_q75_coverage >= 0.90)
  ),
  expected = c(329, 99, 142, 99, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
qc$pass <- as.character(qc$observed) == as.character(qc$expected)
write.csv(qc, file.path(qc_dir, "harrell_bootcal_qc_checks.csv"), row.names = FALSE)

if (!all(qc$pass)) {
  warning("One or more QC checks failed; see ", file.path(qc_dir, "harrell_bootcal_qc_checks.csv"))
}

capture.output(
  sessionInfo(),
  file = file.path(log_dir, "sessionInfo.txt")
)
saveRDS(
  results,
  file = file.path(source_dir, "harrell_bootcal_complete_results.rds")
)

cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Scalar results:", file.path(stat_dir, "harrell_abcloc_scalar_performance.csv"), "\n")
cat("Paired differences:", file.path(stat_dir, "harrell_abcloc_paired_model_differences.csv"), "\n")
cat("Calibration source:", file.path(source_dir, "harrell_bootcal_smooth_calibration_source_data.csv"), "\n")
cat("QC passed:", all(qc$pass), "\n")

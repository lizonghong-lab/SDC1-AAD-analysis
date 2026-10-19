options(stringsAsFactors = FALSE, scipen = 999)

root_dir <- normalizePath(getwd(), mustWork = TRUE)
project_dir <- file.path(root_dir, "SDC1_AAD_HARRELL_BOOTCAL_20260728")
.libPaths(c(file.path(project_dir, "R_library"), .libPaths()))

suppressPackageStartupMessages(library(Hmisc))

run_label <- "FINAL_B1000_LOWESS_ABCLOC_LOCKED"
stat_dir <- file.path(project_dir, "04_results", "statistical_results", run_label)
source_dir <- file.path(project_dir, "04_results", "source_data", run_label)
figure_dir <- file.path(project_dir, "04_results", "figures", run_label)
qc_dir <- file.path(project_dir, "07_quality_control", run_label)
doc_dir <- file.path(project_dir, "00_project_docs")

scalar_file <- file.path(stat_dir, "harrell_abcloc_scalar_performance.csv")
difference_file <- file.path(stat_dir, "harrell_abcloc_paired_model_differences.csv")
rds_file <- file.path(source_dir, "harrell_bootcal_complete_results.rds")

scalar <- read.csv(scalar_file, check.names = FALSE)
differences <- read.csv(difference_file, check.names = FALSE)
results <- readRDS(rds_file)

metric_column <- c(
  "C-index" = "Dxy",
  "Brier score" = "Brier",
  "Log loss" = "LogLoss",
  "Calibration intercept" = "Intercept",
  "Calibration slope" = "Slope"
)

recompute <- function(raw, metric) {
  column <- unname(metric_column[metric])
  column_index <- match(column, names(raw$orig))
  orig <- unname(raw$orig[column_index])
  train <- raw$btrain[, column_index]
  test <- raw$btest[, column_index]
  if (metric == "C-index") {
    orig <- (orig + 1) / 2
    train <- (train + 1) / 2
    test <- (test + 1) / 2
  }
  optimism <- mean(train - test)
  corrected <- orig - optimism
  s <- Hmisc::dualSD(train - 1.25 * test)
  z <- qnorm(0.975)
  c(
    apparent = orig,
    optimism = optimism,
    corrected = corrected,
    lower = corrected - z * unname(s["top"]),
    upper = corrected + z * unname(s["bottom"]),
    successful_resamples = length(train)
  )
}

audit_rows <- list()
i <- 0L
for (population in names(results)) {
  for (model in names(results[[population]]$validations)) {
    raw <- results[[population]]$validations[[model]]$raw
    for (metric in names(metric_column)) {
      i <- i + 1L
      recomputed <- recompute(raw, metric)
      reported <- scalar[
        scalar$population == population &
          scalar$model == model &
          scalar$metric == metric,
      ]
      audit_rows[[i]] <- data.frame(
        population = population,
        model = model,
        metric = metric,
        max_absolute_difference = max(abs(
          recomputed[c("apparent", "optimism", "corrected", "lower", "upper")] -
            unlist(reported[c("apparent", "optimism", "corrected", "lower", "upper")])
        )),
        resamples_match =
          recomputed["successful_resamples"] == reported$successful_resamples,
        stringsAsFactors = FALSE
      )
    }
  }
}
audit <- do.call(rbind, audit_rows)
rownames(audit) <- NULL
write.csv(
  audit,
  file.path(qc_dir, "independent_recalculation_audit.csv"),
  row.names = FALSE
)

knots <- list()
i <- 0L
for (population in names(results)) {
  for (model in names(results[[population]]$fits)) {
    parameters <- results[[population]]$fits[[model]]$Design$parms
    if (!length(parameters)) next
    for (variable in names(parameters)) {
      i <- i + 1L
      knots[[i]] <- data.frame(
        population = population,
        model = model,
        variable = variable,
        knot_1 = unname(parameters[[variable]][1L]),
        knot_2 = unname(parameters[[variable]][2L]),
        knot_3 = unname(parameters[[variable]][3L]),
        stringsAsFactors = FALSE
      )
    }
  }
}
knots <- do.call(rbind, knots)
rownames(knots) <- NULL
write.csv(
  knots,
  file.path(qc_dir, "locked_restricted_cubic_spline_knots.csv"),
  row.names = FALSE
)

figure_files <- list.files(
  figure_dir,
  pattern = "\\.(svg|pdf|tiff|png)$",
  full.names = TRUE
)
figure_audit <- data.frame(
  file = basename(figure_files),
  bytes = file.info(figure_files)$size,
  nonempty = file.info(figure_files)$size > 1000,
  stringsAsFactors = FALSE
)
svg_files <- figure_files[grepl("\\.svg$", figure_files, ignore.case = TRUE)]
figure_audit$editable_svg_text <- NA
for (svg in svg_files) {
  index <- match(basename(svg), figure_audit$file)
  svg_text <- paste(readLines(svg, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  figure_audit$editable_svg_text[index] <- grepl("<text", svg_text, fixed = TRUE)
}
write.csv(
  figure_audit,
  file.path(qc_dir, "figure_export_audit.csv"),
  row.names = FALSE,
  na = ""
)

summary_checks <- data.frame(
  check = c(
    "independent_recalculation_matches",
    "all_resample_counts_match",
    "all_figure_exports_nonempty",
    "all_svg_exports_have_editable_text"
  ),
  pass = c(
    max(audit$max_absolute_difference) < 1e-12,
    all(audit$resamples_match),
    all(figure_audit$nonempty),
    all(figure_audit$editable_svg_text[!is.na(figure_audit$editable_svg_text)])
  ),
  stringsAsFactors = FALSE
)
write.csv(
  summary_checks,
  file.path(qc_dir, "independent_audit_summary.csv"),
  row.names = FALSE
)

if (!all(summary_checks$pass)) stop("Independent audit failed")

format_ci <- function(point, lower, upper, digits = 3L) {
  paste0(
    formatC(point, format = "f", digits = digits),
    " (", formatC(lower, format = "f", digits = digits),
    " to ", formatC(upper, format = "f", digits = digits), ")"
  )
}

scalar$estimate_ci <- mapply(
  format_ci, scalar$corrected, scalar$lower, scalar$upper,
  MoreArgs = list(digits = 3L)
)
scalar_compact <- scalar[, c(
  "population", "model_label", "metric", "estimate_ci"
)]
write.csv(
  scalar_compact,
  file.path(stat_dir, "manuscript_ready_scalar_performance.csv"),
  row.names = FALSE
)

differences$estimate_ci <- mapply(
  format_ci,
  differences$corrected_improvement,
  differences$lower,
  differences$upper,
  MoreArgs = list(digits = 3L)
)
difference_compact <- differences[, c(
  "population", "comparison", "metric", "positive_value_favours", "estimate_ci"
)]
write.csv(
  difference_compact,
  file.path(stat_dir, "manuscript_ready_paired_differences.csv"),
  row.names = FALSE
)

cat("Independent audit PASS\n")
cat("Maximum numerical discrepancy:", max(audit$max_absolute_difference), "\n")
cat("Figure exports checked:", nrow(figure_audit), "\n")


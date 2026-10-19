options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(root, "02_analysis_data")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
log_dir <- file.path(root, "04_results", "logs")
for (directory in c(result_dir, source_dir, log_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")

auc_rank <- function(y, p) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_curve <- function(y, p) {
  order_index <- order(p, decreasing = TRUE)
  y <- y[order_index]
  p <- p[order_index]
  endpoint <- c(which(diff(p) != 0), length(p))
  data.frame(
    false_positive_rate = c(0, cumsum(y == 0)[endpoint] / sum(y == 0)),
    true_positive_rate = c(0, cumsum(y == 1)[endpoint] / sum(y == 1)),
    threshold = c(Inf, p[endpoint]),
    stringsAsFactors = FALSE
  )
}

cluster_outcome <- function(data) {
  split_y <- split(data$aad, data$cluster_id)
  inconsistent <- names(split_y)[vapply(split_y, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(inconsistent)) stop("Cluster IDs span outcomes")
  vapply(split_y, function(x) x[1], numeric(1))
}

bootstrap_cluster_indices <- function(data) {
  y_by_cluster <- cluster_outcome(data)
  selected <- unlist(lapply(sort(unique(y_by_cluster)), function(level) {
    ids <- names(y_by_cluster)[y_by_cluster == level]
    sample(ids, length(ids), replace = TRUE)
  }), use.names = FALSE)
  unlist(lapply(selected, function(id) which(data$cluster_id == id)), use.names = FALSE)
}

symptomatic <- read.csv(file.path(analysis_dir, "elisa_evaluation_primary.csv"), stringsAsFactors = FALSE, check.names = FALSE)
symptomatic$aad <- as.integer(symptomatic$aad)
symptomatic$cluster_id <- ifelse(is.na(symptomatic$cluster_id) | !nzchar(symptomatic$cluster_id), symptomatic$row_uid, symptomatic$cluster_id)
required <- c("aad", "sdc1_source_value", "d_dimer_mg_l", "group", "cluster_id", "row_uid")
data <- symptomatic[complete.cases(symptomatic[, required]) & symptomatic$sdc1_source_value > 0 & symptomatic$d_dimer_mg_l > 0, , drop = FALSE]
data$SDC1 <- log2(data$sdc1_source_value)
data$D_dimer <- log2(data$d_dimer_mg_l)

biomarkers <- c("D_dimer", "SDC1")
population_index <- list(
  "All symptomatic comparators" = rep(TRUE, nrow(data)),
  "AAD versus PE" = data$aad == 1 | data$group == "PE"
)

roc_source <- do.call(rbind, lapply(names(population_index), function(population_name) {
  index <- population_index[[population_name]]
  do.call(rbind, lapply(biomarkers, function(biomarker) {
    curve <- roc_curve(data$aad[index], data[[biomarker]][index])
    curve$population <- population_name
    curve$biomarker <- biomarker
    curve
  }))
}))

point_auc <- do.call(rbind, lapply(names(population_index), function(population_name) {
  index <- population_index[[population_name]]
  do.call(rbind, lapply(biomarkers, function(biomarker) data.frame(
    population = population_name,
    biomarker = biomarker,
    estimate = auc_rank(data$aad[index], data[[biomarker]][index]),
    n_total = sum(index),
    aad_cases = sum(data$aad[index] == 1),
    non_aad_controls = sum(data$aad[index] == 0),
    stringsAsFactors = FALSE
  )))
}))

B <- 2000L
bootstrap_draws <- vector("list", B)
set.seed(20284000)
for (b in seq_len(B)) {
  index <- bootstrap_cluster_indices(data)
  bootstrap_data <- data[index, , drop = FALSE]
  bootstrap_populations <- list(
    "All symptomatic comparators" = rep(TRUE, nrow(bootstrap_data)),
    "AAD versus PE" = bootstrap_data$aad == 1 | bootstrap_data$group == "PE"
  )
  bootstrap_draws[[b]] <- do.call(rbind, lapply(names(bootstrap_populations), function(population_name) {
    rows <- bootstrap_populations[[population_name]]
    do.call(rbind, lapply(biomarkers, function(biomarker) data.frame(
      bootstrap_replicate = b,
      population = population_name,
      biomarker = biomarker,
      auc = auc_rank(bootstrap_data$aad[rows], bootstrap_data[[biomarker]][rows]),
      stringsAsFactors = FALSE
    )))
  }))
  if (b %% 200L == 0L) cat("Standalone paired cluster bootstrap", b, "of", B, "\n")
}
bootstrap_draws <- do.call(rbind, bootstrap_draws)

auc_results <- point_auc
auc_results$lower_95 <- NA_real_
auc_results$upper_95 <- NA_real_
for (i in seq_len(nrow(auc_results))) {
  row <- auc_results[i, , drop = FALSE]
  draws <- bootstrap_draws$auc[bootstrap_draws$population == row$population & bootstrap_draws$biomarker == row$biomarker]
  auc_results$lower_95[i] <- quantile(draws, 0.025, na.rm = TRUE)
  auc_results$upper_95[i] <- quantile(draws, 0.975, na.rm = TRUE)
}
auc_results$validation <- "Prespecified continuous single-marker score; 2000 paired outcome-stratified author-resolved patient-key bootstrap resamples"

paired_auc_difference <- do.call(rbind, lapply(names(population_index), function(population_name) {
  sdc1_point <- auc_results$estimate[auc_results$population == population_name & auc_results$biomarker == "SDC1"]
  d_dimer_point <- auc_results$estimate[auc_results$population == population_name & auc_results$biomarker == "D_dimer"]
  sdc1_draws <- bootstrap_draws$auc[bootstrap_draws$population == population_name & bootstrap_draws$biomarker == "SDC1"]
  d_dimer_draws <- bootstrap_draws$auc[bootstrap_draws$population == population_name & bootstrap_draws$biomarker == "D_dimer"]
  difference_draws <- sdc1_draws - d_dimer_draws
  data.frame(
    population = population_name,
    comparison = "SDC1 minus D-dimer",
    difference = sdc1_point - d_dimer_point,
    lower_95 = quantile(difference_draws, 0.025, na.rm = TRUE),
    upper_95 = quantile(difference_draws, 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write_csv(roc_source, file.path(source_dir, "fig4_standalone_biomarker_roc_curves.csv"))
write_csv(bootstrap_draws, file.path(source_dir, "fig4_standalone_biomarker_roc_bootstrap_draws.csv"))
write_csv(data[, c("row_uid", "source_id", "cluster_id", "cluster_rule", "group", "aad", "sdc1_source_value", "d_dimer_mg_l")], file.path(source_dir, "fig4_standalone_biomarker_roc_source_data.csv"))
write_csv(auc_results, file.path(result_dir, "fig4_standalone_biomarker_roc_auc.csv"))
write_csv(paired_auc_difference, file.path(result_dir, "fig4_standalone_biomarker_roc_paired_difference.csv"))

log_lines <- c(
  paste("Fair same-patient standalone biomarker cohort:", nrow(data), "records;", sum(data$aad), "AAD"),
  paste("AAD versus PE subgroup:", sum(population_index[["AAD versus PE"]]), "records"),
  capture.output(print(auc_results[, c("population", "biomarker", "estimate", "lower_95", "upper_95")], row.names = FALSE)),
  capture.output(print(paired_auc_difference, row.names = FALSE)),
  "No clinical covariates, combined biomarker model, fitted diagnostic cut-point or Youden J statistic is used in these standalone ROC panels."
)
writeLines(log_lines, file.path(log_dir, "fig4_standalone_biomarker_roc_summary.txt"), useBytes = TRUE)
cat(paste(log_lines, collapse = "\n"), "\n")

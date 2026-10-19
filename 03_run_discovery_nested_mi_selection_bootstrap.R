options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1L) as.integer(args[[1L]]) else 300L
M <- if (length(args) >= 2L) as.integer(args[[2L]]) else 50L
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else 4L
run_label <- if (length(args) >= 4L) args[[4L]] else paste0("B", B, "_M", M)
if (!is.finite(B) || B < 1L) stop("B must be positive")
if (!is.finite(M) || M < 1L) stop("M must be positive")
if (!is.finite(workers) || workers < 1L) stop("workers must be positive")

root_dir <- normalizePath(getwd(), mustWork = TRUE)
project_dir <- file.path(root_dir, "SDC1_AAD_HARRELL_BOOTCAL_20260728")
old_project <- normalizePath(
  file.path(
    root_dir,
    "SDC1_AAD_REANALYSIS_V6_1_REVIEWER_REFINED_20260715"
  ),
  mustWork = TRUE
)
old_library <- file.path(old_project, "08_environment", "Rlib")
.libPaths(c(old_library, file.path(project_dir, "R_library"), .libPaths()))

suppressPackageStartupMessages(library(limma))

source_dir <- file.path(
  project_dir, "04_results", "source_data", "discovery_nested_selection", run_label
)
stat_dir <- file.path(
  project_dir, "04_results", "statistical_results", "discovery_nested_selection", run_label
)
log_dir <- file.path(
  project_dir, "04_results", "logs", "discovery_nested_selection", run_label
)
qc_dir <- file.path(
  project_dir, "07_quality_control", "discovery_nested_selection", run_label
)
for (d in c(source_dir, stat_dir, log_dir, qc_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "nested_mi_selection_bootstrap.log")
sink(log_file, split = TRUE)
on.exit({
  while (sink.number() > 0L) sink()
}, add = TRUE)

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

prepare_design <- function(metadata, positive_group, reference_group) {
  rows <- which(metadata$group %in% c(positive_group, reference_group))
  d <- metadata[rows, , drop = FALSE]
  complete <- complete.cases(d[, c(
    "group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l"
  )])
  rows <- rows[complete]
  d <- d[complete, , drop = FALSE]
  d$aad_indicator <- as.integer(d$group == positive_group)
  d$age_scaled <- safe_z(d$age_years)
  d$wbc_scaled <- safe_z(log1p(d$wbc_10e9_l))
  d$crp_scaled <- safe_z(log1p(d$crp_mg_l))
  design <- model.matrix(
    ~ aad_indicator + age_scaled + female + hypertension + wbc_scaled + crp_scaled,
    data = d
  )
  list(rows = rows, design = design)
}

outcome_blind_filter <- function(log2_matrix, minimum_detection = 0.50) {
  detection <- rowMeans(is.finite(log2_matrix))
  variability <- apply(log2_matrix, 1, IQR, na.rm = TRUE)
  is.finite(detection) & detection >= minimum_detection &
    is.finite(variability) & variability > 0
}

median_normalise <- function(log2_matrix) {
  sample_medians <- apply(log2_matrix, 2, median, na.rm = TRUE)
  grand_median <- median(sample_medians, na.rm = TRUE)
  sweep(log2_matrix, 2, sample_medians - grand_median, "-")
}

draw_left_censored <- function(matrix_value, shift = 1.8, scale = 0.3) {
  output <- matrix_value
  for (j in seq_len(ncol(matrix_value))) {
    observed <- matrix_value[, j]
    observed <- observed[is.finite(observed)]
    cutoff <- as.numeric(quantile(observed, 0.01, na.rm = TRUE, type = 8))
    robust_scale <- IQR(observed, na.rm = TRUE, type = 8) / 1.349
    if (!is.finite(robust_scale) || robust_scale <= 0.10) {
      robust_scale <- mad(observed, na.rm = TRUE)
    }
    if (!is.finite(robust_scale) || robust_scale <= 0.10) {
      robust_scale <- sd(observed, na.rm = TRUE)
    }
    if (!is.finite(robust_scale) || robust_scale <= 0.10) robust_scale <- 1
    mu <- cutoff - shift * robust_scale
    sigma <- max(scale * robust_scale, 0.05)
    missing <- which(!is.finite(matrix_value[, j]))
    if (length(missing)) {
      upper_probability <- pnorm(cutoff, mean = mu, sd = sigma)
      upper_probability <- min(max(upper_probability, 1e-10), 1 - 1e-10)
      u <- runif(length(missing), min = 1e-10, max = upper_probability)
      output[missing, j] <- qnorm(u, mean = mu, sd = sigma)
    }
  }
  output
}

fit_adjusted_limma <- function(matrix_value, metadata, positive_group, reference_group) {
  prepared <- prepare_design(metadata, positive_group, reference_group)
  fit <- limma::lmFit(matrix_value[, prepared$rows, drop = FALSE], prepared$design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  coefficient <- which(colnames(prepared$design) == "aad_indicator")
  beta <- fit$coefficients[, coefficient]
  se <- fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
  list(beta = as.numeric(beta), se2 = as.numeric(se^2))
}

pool_rubin <- function(beta_matrix, se2_matrix) {
  M_local <- ncol(beta_matrix)
  qbar <- rowMeans(beta_matrix)
  ubar <- rowMeans(se2_matrix)
  between <- if (M_local > 1L) {
    apply(beta_matrix, 1, var)
  } else {
    rep(0, nrow(beta_matrix))
  }
  total <- ubar + (1 + 1 / M_local) * between
  total[!is.finite(total) | total <= 0] <- NA_real_
  relative_increase <-
    (1 + 1 / M_local) * between / pmax(ubar, .Machine$double.eps)
  degrees_freedom <-
    (M_local - 1) * (1 + 1 / pmax(relative_increase, .Machine$double.eps))^2
  degrees_freedom[!is.finite(degrees_freedom) |
                    between <= .Machine$double.eps] <- 1e9
  statistic <- qbar / sqrt(total)
  list(
    beta = qbar,
    se = sqrt(total),
    p_one_sided_aad_higher = pt(
      statistic, df = degrees_freedom, lower.tail = FALSE
    )
  )
}

run_nested_mi_ranking <- function(
    log2_matrix, metadata, protein_ids, gene_symbols,
    M, seed, shift = 1.8, scale = 0.3) {
  set.seed(seed)
  matrix_norm <- median_normalise(log2_matrix)
  n_proteins <- nrow(matrix_norm)
  beta_ami <- se2_ami <- beta_hc <- se2_hc <-
    matrix(NA_real_, nrow = n_proteins, ncol = M)

  for (m in seq_len(M)) {
    imputed <- draw_left_censored(matrix_norm, shift = shift, scale = scale)
    ami <- fit_adjusted_limma(imputed, metadata, "AAD", "AMI")
    hc <- fit_adjusted_limma(imputed, metadata, "AAD", "Healthy control")
    beta_ami[, m] <- ami$beta
    se2_ami[, m] <- ami$se2
    beta_hc[, m] <- hc$beta
    se2_hc[, m] <- hc$se2
  }

  pooled_ami <- pool_rubin(beta_ami, se2_ami)
  pooled_hc <- pool_rubin(beta_hc, se2_hc)
  joint <- data.frame(
    protein_group_id = protein_ids,
    gene_symbol = gene_symbols,
    beta_ami = pooled_ami$beta,
    beta_hc = pooled_hc$beta,
    p_ami_greater = pooled_ami$p_one_sided_aad_higher,
    p_hc_greater = pooled_hc$p_one_sided_aad_higher,
    stringsAsFactors = FALSE
  )
  joint$joint_iut_p <- pmax(joint$p_ami_greater, joint$p_hc_greater)
  joint$minimum_effect <- pmin(joint$beta_ami, joint$beta_hc)
  joint <- joint[
    order(joint$joint_iut_p, -joint$minimum_effect, joint$protein_group_id),
    , drop = FALSE
  ]
  joint$rank <- seq_len(nrow(joint))
  joint
}

bootstrap_worker <- function(task) {
  b <- task$b
  sampled <- task$sampled
  boot_meta <- sample_metadata[sampled, , drop = FALSE]
  boot_matrix_full <- raw_log2[, sampled, drop = FALSE]
  colnames(boot_matrix_full) <- paste0(
    boot_meta$sample_id, "_boot", seq_along(sampled)
  )
  keep <- outcome_blind_filter(boot_matrix_full, minimum_detection = 0.50)
  boot_matrix <- boot_matrix_full[keep, , drop = FALSE]
  protein_ids <- all_protein_ids[keep]
  gene_symbols <- all_gene_symbols[keep]
  joint <- run_nested_mi_ranking(
    boot_matrix,
    boot_meta,
    protein_ids,
    gene_symbols,
    M = M_nested,
    seed = 202607280L + b,
    shift = 1.8,
    scale = 0.3
  )
  sdc1 <- joint[joint$protein_group_id == sdc1_id, , drop = FALSE]
  summary <- data.frame(
    bootstrap = b,
    proteins_retained = nrow(joint),
    imputations_repeated = M_nested,
    sdc1_retained = nrow(sdc1) == 1L,
    sdc1_rank = if (nrow(sdc1)) sdc1$rank else NA_integer_,
    sdc1_log2FC_aad_vs_ami = if (nrow(sdc1)) sdc1$beta_ami else NA_real_,
    sdc1_log2FC_aad_vs_healthy = if (nrow(sdc1)) sdc1$beta_hc else NA_real_,
    sdc1_aad_higher_both = if (nrow(sdc1)) {
      sdc1$beta_ami > 0 && sdc1$beta_hc > 0
    } else {
      FALSE
    },
    top_candidate_id = joint$protein_group_id[1L],
    top_candidate_gene = joint$gene_symbol[1L],
    stringsAsFactors = FALSE
  )
  top10 <- joint[seq_len(min(10L, nrow(joint))), c(
    "protein_group_id", "gene_symbol", "beta_ami", "beta_hc",
    "joint_iut_p", "minimum_effect", "rank"
  )]
  top10$bootstrap <- b
  list(summary = summary, top10 = top10)
}

input_file <- file.path(old_project, "02_analysis_data", "proteomics_raw_deidentified.rds")
raw <- readRDS(input_file)
raw_log2 <- log2(raw$intensity)
rownames(raw_log2) <- raw$protein_group_id
colnames(raw_log2) <- raw$sample_metadata$sample_id
sample_metadata <- raw$sample_metadata
all_protein_ids <- raw$protein_group_id
all_gene_symbols <- raw$gene_symbol
sdc1_candidates <- all_protein_ids[all_gene_symbols == "SDC1"]
if (length(sdc1_candidates) != 1L) stop("Expected one SDC1 protein group")
sdc1_id <- sdc1_candidates[[1L]]
M_nested <- M

set.seed(2026072803L)
group_indices <- split(seq_len(nrow(sample_metadata)), sample_metadata$group)
tasks <- lapply(seq_len(B), function(b) {
  sampled <- unlist(lapply(
    group_indices,
    function(index) sample(index, length(index), replace = TRUE)
  ), use.names = FALSE)
  list(b = b, sampled = sampled)
})
saveRDS(tasks, file.path(source_dir, "bootstrap_sample_indices.rds"))

cat("Nested discovery selection bootstrap\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Input:", input_file, "\n")
cat("Input MD5:", unname(tools::md5sum(input_file)), "\n")
cat("B:", B, "M per bootstrap:", M, "workers:", workers, "\n")
cat("Samples:", nrow(sample_metadata), "proteins before filtering:", nrow(raw_log2), "\n")
cat("limma:", as.character(packageVersion("limma")), "\n")

cluster <- NULL
if (workers > 1L) {
  old_library_worker <- old_library
  project_library_worker <- file.path(project_dir, "R_library")
  cluster <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterExport(
    cluster,
    varlist = c(
      "old_library_worker", "project_library_worker",
      "safe_z", "prepare_design", "outcome_blind_filter", "median_normalise",
      "draw_left_censored", "fit_adjusted_limma", "pool_rubin",
      "run_nested_mi_ranking", "bootstrap_worker", "sample_metadata",
      "raw_log2", "all_protein_ids", "all_gene_symbols", "sdc1_id", "M_nested"
    ),
    envir = environment()
  )
  parallel::clusterEvalQ(cluster, {
    .libPaths(c(old_library_worker, project_library_worker, .libPaths()))
    suppressPackageStartupMessages(library(limma))
    NULL
  })
}

chunk_size <- max(workers * 2L, 4L)
chunks <- split(seq_len(B), ceiling(seq_len(B) / chunk_size))
all_results <- vector("list", B)
start_time <- proc.time()[[3L]]
for (chunk_index in seq_along(chunks)) {
  ids <- chunks[[chunk_index]]
  if (workers > 1L) {
    chunk_results <- parallel::parLapplyLB(cluster, tasks[ids], bootstrap_worker)
  } else {
    chunk_results <- lapply(tasks[ids], bootstrap_worker)
  }
  all_results[ids] <- chunk_results
  elapsed <- proc.time()[[3L]] - start_time
  cat(
    "Completed", max(ids), "of", B,
    "elapsed_seconds", round(elapsed, 1), "\n"
  )
  flush.console()
}
if (!is.null(cluster)) {
  parallel::stopCluster(cluster)
  cluster <- NULL
}

summary_draws <- do.call(rbind, lapply(all_results, `[[`, "summary"))
top10_draws <- do.call(rbind, lapply(all_results, `[[`, "top10"))
rownames(summary_draws) <- NULL
rownames(top10_draws) <- NULL
write_csv(
  summary_draws,
  file.path(source_dir, "nested_mi_selection_bootstrap_draws.csv")
)
write_csv(
  top10_draws,
  file.path(source_dir, "nested_mi_selection_top10_draws.csv")
)

retained <- sum(summary_draws$sdc1_retained)
rank1_all <- sum(summary_draws$sdc1_rank == 1L, na.rm = TRUE)
rank1_retained <- sum(summary_draws$sdc1_rank == 1L, na.rm = TRUE)
retain_ci <- binom.test(retained, B)$conf.int
rank1_all_ci <- binom.test(rank1_all, B)$conf.int
conditional_ci <- if (retained > 0L) {
  binom.test(rank1_retained, retained)$conf.int
} else {
  c(NA_real_, NA_real_)
}

stability <- data.frame(
  bootstrap_repetitions = B,
  imputations_per_bootstrap = M,
  retained_count = retained,
  retained_frequency = retained / B,
  retained_ci_lower = retain_ci[1L],
  retained_ci_upper = retain_ci[2L],
  rank1_count_all = rank1_all,
  rank1_frequency_all = rank1_all / B,
  rank1_all_ci_lower = rank1_all_ci[1L],
  rank1_all_ci_upper = rank1_all_ci[2L],
  rank1_count_retained = rank1_retained,
  rank1_frequency_retained = if (retained) rank1_retained / retained else NA_real_,
  rank1_retained_ci_lower = conditional_ci[1L],
  rank1_retained_ci_upper = conditional_ci[2L],
  median_rank_retained = median(summary_draws$sdc1_rank, na.rm = TRUE),
  rank_q025_retained = unname(quantile(
    summary_draws$sdc1_rank, 0.025, na.rm = TRUE
  )),
  rank_q975_retained = unname(quantile(
    summary_draws$sdc1_rank, 0.975, na.rm = TRUE
  )),
  stringsAsFactors = FALSE
)
write_csv(
  stability,
  file.path(stat_dir, "nested_mi_sdc1_selection_stability.csv")
)

qa <- data.frame(
  check = c(
    "all_bootstraps_returned",
    "all_bootstraps_repeated_requested_imputations",
    "all_top_candidates_recorded",
    "sdc1_rank_only_when_retained",
    "all_strata_present_in_every_resample"
  ),
  pass = c(
    nrow(summary_draws) == B && length(unique(summary_draws$bootstrap)) == B,
    all(summary_draws$imputations_repeated == M),
    all(nzchar(summary_draws$top_candidate_id)),
    all(is.finite(summary_draws$sdc1_rank) == summary_draws$sdc1_retained),
    all(vapply(tasks, function(task) {
      all(table(sample_metadata$group[task$sampled]) > 0)
    }, logical(1)))
  ),
  stringsAsFactors = FALSE
)
write_csv(qa, file.path(qc_dir, "nested_mi_selection_bootstrap_qc.csv"))
if (!all(qa$pass)) stop("Nested selection bootstrap QA failed")

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
print(stability)
cat("QA passed: TRUE\n")

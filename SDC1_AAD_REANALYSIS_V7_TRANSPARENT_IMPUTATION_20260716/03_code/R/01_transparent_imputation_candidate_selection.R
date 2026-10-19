options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_V7_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source_root <- normalizePath(file.path(root, "..", "SDC1_AAD_REANALYSIS_V6_1_REVIEWER_REFINED_20260715"), winslash = "/", mustWork = TRUE)
source_lib <- file.path(source_root, "08_environment", "Rlib")
if (dir.exists(source_lib)) .libPaths(c(source_lib, .libPaths()))

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

set.seed(20260716)

M_primary <- as.integer(Sys.getenv("V7_MI_M", unset = "50"))
B_boot <- as.integer(Sys.getenv("V7_BOOT_B", unset = "300"))
M_sensitivity <- as.integer(Sys.getenv("V7_SENS_M", unset = "20"))
primary_shift <- 1.8
primary_scale <- 0.3

analysis_dir <- file.path(root, "02_analysis_data")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
main_figure_dir <- file.path(root, "04_results", "figures", "main")
supp_figure_dir <- file.path(root, "04_results", "figures", "supplementary")
log_dir <- file.path(root, "04_results", "logs")
for (directory in c(analysis_dir, result_dir, source_dir, main_figure_dir, supp_figure_dir, log_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

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
  complete <- complete.cases(d[, c("group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l")])
  rows <- rows[complete]
  d <- d[complete, , drop = FALSE]
  d$aad_indicator <- as.integer(d$group == positive_group)
  d$age_scaled <- safe_z(d$age_years)
  d$wbc_scaled <- safe_z(log1p(d$wbc_10e9_l))
  d$crp_scaled <- safe_z(log1p(d$crp_mg_l))
  design <- model.matrix(~ aad_indicator + age_scaled + female + hypertension + wbc_scaled + crp_scaled, data = d)
  list(rows = rows, metadata = d, design = design)
}

outcome_blind_filter <- function(log2_matrix, minimum_detection = 0.50) {
  detection <- rowMeans(is.finite(log2_matrix))
  variability <- apply(log2_matrix, 1, IQR, na.rm = TRUE)
  is.finite(detection) & detection >= minimum_detection & is.finite(variability) & variability > 0
}

median_normalise <- function(log2_matrix) {
  sample_medians <- apply(log2_matrix, 2, median, na.rm = TRUE)
  grand_median <- median(sample_medians, na.rm = TRUE)
  normalised <- sweep(log2_matrix, 2, sample_medians - grand_median, "-")
  list(matrix = normalised, sample_medians = sample_medians, grand_median = grand_median)
}

draw_left_censored <- function(matrix_value, shift = 1.8, scale = 0.3) {
  output <- matrix_value
  diagnostics <- vector("list", ncol(matrix_value))
  for (j in seq_len(ncol(matrix_value))) {
    observed <- matrix_value[, j]
    observed <- observed[is.finite(observed)]
    cutoff <- as.numeric(quantile(observed, 0.01, na.rm = TRUE, type = 8))
    robust_scale <- IQR(observed, na.rm = TRUE, type = 8) / 1.349
    if (!is.finite(robust_scale) || robust_scale <= 0.10) robust_scale <- mad(observed, na.rm = TRUE)
    if (!is.finite(robust_scale) || robust_scale <= 0.10) robust_scale <- sd(observed, na.rm = TRUE)
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
    diagnostics[[j]] <- data.frame(
      sample_id = colnames(matrix_value)[j],
      observed_n = length(observed),
      missing_n = length(missing),
      cutoff_1pct = cutoff,
      robust_scale = robust_scale,
      imputation_mean = mu,
      imputation_sd = sigma,
      stringsAsFactors = FALSE
    )
  }
  list(matrix = output, diagnostics = do.call(rbind, diagnostics))
}

fit_adjusted_limma <- function(matrix_value, metadata, positive_group, reference_group) {
  prepared <- prepare_design(metadata, positive_group, reference_group)
  fit <- lmFit(matrix_value[, prepared$rows, drop = FALSE], prepared$design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  coefficient <- which(colnames(prepared$design) == "aad_indicator")
  beta <- fit$coefficients[, coefficient]
  se <- fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
  data.frame(
    protein_group_id = rownames(matrix_value),
    beta = as.numeric(beta),
    se2 = as.numeric(se^2),
    stringsAsFactors = FALSE
  )
}

pool_rubin <- function(beta_matrix, se2_matrix) {
  M <- ncol(beta_matrix)
  qbar <- rowMeans(beta_matrix)
  ubar <- rowMeans(se2_matrix)
  between <- if (M > 1) apply(beta_matrix, 1, var) else rep(0, nrow(beta_matrix))
  total <- ubar + (1 + 1 / M) * between
  total[!is.finite(total) | total <= 0] <- NA_real_
  fraction_missing_information <- (1 + 1 / M) * between / total
  relative_increase <- (1 + 1 / M) * between / pmax(ubar, .Machine$double.eps)
  degrees_freedom <- (M - 1) * (1 + 1 / pmax(relative_increase, .Machine$double.eps))^2
  degrees_freedom[!is.finite(degrees_freedom) | between <= .Machine$double.eps] <- 1e9
  statistic <- qbar / sqrt(total)
  p_two_sided <- 2 * pt(-abs(statistic), df = degrees_freedom)
  p_greater <- pt(statistic, df = degrees_freedom, lower.tail = FALSE)
  data.frame(
    pooled_log2FC = qbar,
    pooled_se = sqrt(total),
    ci_lower_95 = qbar - qt(0.975, df = degrees_freedom) * sqrt(total),
    ci_upper_95 = qbar + qt(0.975, df = degrees_freedom) * sqrt(total),
    t_rubin = statistic,
    df_rubin = degrees_freedom,
    p_value = p_two_sided,
    p_one_sided_aad_higher = p_greater,
    fraction_missing_information = fraction_missing_information,
    stringsAsFactors = FALSE
  )
}

run_multiple_imputation <- function(log2_matrix, metadata, protein_lookup, shift, scale, M, seed, keep_example = FALSE) {
  set.seed(seed)
  normalisation <- median_normalise(log2_matrix)
  matrix_norm <- normalisation$matrix
  n_proteins <- nrow(matrix_norm)
  beta_ami <- se2_ami <- beta_hc <- se2_hc <- matrix(NA_real_, nrow = n_proteins, ncol = M)
  first_imputed <- NULL
  first_diagnostics <- NULL
  sdc1_draws <- vector("list", M)

  for (m in seq_len(M)) {
    imputed <- draw_left_censored(matrix_norm, shift = shift, scale = scale)
    if (keep_example && m == 1) {
      first_imputed <- imputed$matrix
      first_diagnostics <- imputed$diagnostics
    }
    ami <- fit_adjusted_limma(imputed$matrix, metadata, "AAD", "AMI")
    hc <- fit_adjusted_limma(imputed$matrix, metadata, "AAD", "Healthy control")
    beta_ami[, m] <- ami$beta
    se2_ami[, m] <- ami$se2
    beta_hc[, m] <- hc$beta
    se2_hc[, m] <- hc$se2
    sdc1_row <- which(protein_lookup$gene_symbol == "SDC1")
    sdc1_draws[[m]] <- data.frame(
      imputation = m,
      shift = shift,
      scale = scale,
      contrast = c("AAD versus AMI", "AAD versus Healthy control"),
      log2FC = c(ami$beta[sdc1_row], hc$beta[sdc1_row]),
      standard_error = sqrt(c(ami$se2[sdc1_row], hc$se2[sdc1_row])),
      stringsAsFactors = FALSE
    )
    if (m %% 10 == 0 || m == M) cat("Multiple imputation", m, "of", M, "\n")
  }

  pooled_ami <- cbind(protein_lookup, pool_rubin(beta_ami, se2_ami))
  pooled_hc <- cbind(protein_lookup, pool_rubin(beta_hc, se2_hc))
  pooled_ami$contrast <- "AAD versus AMI"
  pooled_hc$contrast <- "AAD versus Healthy control"
  pooled_ami$bh_fdr <- p.adjust(pooled_ami$p_value, method = "BH")
  pooled_hc$bh_fdr <- p.adjust(pooled_hc$p_value, method = "BH")

  joint <- merge(
    pooled_ami[, c("protein_group_id", "gene_symbol", "pooled_log2FC", "p_value", "p_one_sided_aad_higher", "bh_fdr")],
    pooled_hc[, c("protein_group_id", "pooled_log2FC", "p_value", "p_one_sided_aad_higher", "bh_fdr")],
    by = "protein_group_id", suffixes = c("_ami", "_healthy"), sort = FALSE
  )
  joint$joint_iut_p <- pmax(joint$p_one_sided_aad_higher_ami, joint$p_one_sided_aad_higher_healthy)
  joint$joint_iut_bh <- p.adjust(joint$joint_iut_p, method = "BH")
  joint$minimum_log2FC <- pmin(joint$pooled_log2FC_ami, joint$pooled_log2FC_healthy)
  joint$consistent_aad_higher <- joint$pooled_log2FC_ami > 0 & joint$pooled_log2FC_healthy > 0
  order_index <- order(joint$joint_iut_p, -joint$minimum_log2FC, joint$protein_group_id)
  joint <- joint[order_index, , drop = FALSE]
  joint$joint_rank <- seq_len(nrow(joint))

  list(
    pooled_ami = pooled_ami,
    pooled_hc = pooled_hc,
    joint = joint,
    sdc1_draws = do.call(rbind, sdc1_draws),
    first_imputed = first_imputed,
    first_diagnostics = first_diagnostics,
    normalised_matrix = matrix_norm,
    normalisation = normalisation
  )
}

raw <- readRDS(file.path(source_root, "02_analysis_data", "proteomics_raw_deidentified.rds"))
raw_log2 <- log2(raw$intensity)
rownames(raw_log2) <- raw$protein_group_id
colnames(raw_log2) <- raw$sample_metadata$sample_id
primary_keep <- outcome_blind_filter(raw_log2, minimum_detection = 0.50)
primary_matrix <- raw_log2[primary_keep, , drop = FALSE]
protein_lookup <- data.frame(
  protein_group_id = raw$protein_group_id[primary_keep],
  gene_symbol = raw$gene_symbol[primary_keep],
  stringsAsFactors = FALSE
)

if (sum(protein_lookup$gene_symbol == "SDC1") != 1L) stop("Expected exactly one eligible SDC1 protein group")

candidate_flow <- data.frame(
  stage = c("Locked protein groups", "Detection >=50% and non-zero IQR", "Primary imputation candidate universe"),
  n_protein_groups = c(nrow(raw_log2), sum(primary_keep), nrow(primary_matrix)),
  stringsAsFactors = FALSE
)
write_csv(candidate_flow, file.path(result_dir, "candidate_flow_v7.csv"))

cat("Primary candidate universe:", nrow(primary_matrix), "proteins\n")
primary <- run_multiple_imputation(
  log2_matrix = primary_matrix,
  metadata = raw$sample_metadata,
  protein_lookup = protein_lookup,
  shift = primary_shift,
  scale = primary_scale,
  M = M_primary,
  seed = 2026071601,
  keep_example = TRUE
)

write_csv(primary$pooled_ami, file.path(source_dir, "imputation_pooled_AAD_vs_AMI_v7.csv"))
write_csv(primary$pooled_hc, file.path(source_dir, "imputation_pooled_AAD_vs_healthy_v7.csv"))
write_csv(primary$joint, file.path(source_dir, "imputation_joint_candidate_ranking_v7.csv"))
write_csv(primary$sdc1_draws, file.path(source_dir, "sdc1_multiple_imputation_draws_v7.csv"))
write_csv(primary$first_diagnostics, file.path(source_dir, "imputation_sample_diagnostics_v7.csv"))

sdc1_id <- protein_lookup$protein_group_id[protein_lookup$gene_symbol == "SDC1"]
sdc1_raw_row <- match(sdc1_id, rownames(raw_log2))
sdc1_primary_row <- match(sdc1_id, rownames(primary$normalised_matrix))
sdc1_example <- data.frame(
  sample_id = raw$sample_metadata$sample_id,
  group = raw$sample_metadata$group,
  observed = is.finite(primary$normalised_matrix[sdc1_primary_row, ]),
  normalised_log2_value = primary$normalised_matrix[sdc1_primary_row, ],
  example_imputed_log2_value = primary$first_imputed[sdc1_primary_row, ],
  stringsAsFactors = FALSE
)
write_csv(sdc1_example, file.path(source_dir, "sdc1_observed_and_example_imputed_values_v7.csv"))

detect_table <- do.call(rbind, lapply(split(seq_len(nrow(raw$sample_metadata)), raw$sample_metadata$group), function(index) {
  detected <- is.finite(raw_log2[sdc1_raw_row, index])
  data.frame(
    group = raw$sample_metadata$group[index[1]],
    detected = sum(detected),
    total = length(index),
    detection_fraction = mean(detected),
    stringsAsFactors = FALSE
  )
}))
write_csv(detect_table, file.path(result_dir, "sdc1_detection_by_group_v7.csv"))

primary_sdc1_ami <- primary$pooled_ami[primary$pooled_ami$protein_group_id == sdc1_id, ]
primary_sdc1_hc <- primary$pooled_hc[primary$pooled_hc$protein_group_id == sdc1_id, ]
primary_sdc1_joint <- primary$joint[primary$joint$protein_group_id == sdc1_id, ]

# Pre-fixed sensitivity grid.
sensitivity_grid <- expand.grid(
  shift = c(1.5, 1.8, 2.1),
  scale = c(0.2, 0.3, 0.5),
  stringsAsFactors = FALSE
)
sensitivity_results <- vector("list", nrow(sensitivity_grid))
for (i in seq_len(nrow(sensitivity_grid))) {
  cat("Sensitivity", i, "of", nrow(sensitivity_grid), "shift", sensitivity_grid$shift[i], "scale", sensitivity_grid$scale[i], "\n")
  analysis <- run_multiple_imputation(
    log2_matrix = primary_matrix,
    metadata = raw$sample_metadata,
    protein_lookup = protein_lookup,
    shift = sensitivity_grid$shift[i],
    scale = sensitivity_grid$scale[i],
    M = M_sensitivity,
    seed = 2026071700 + i,
    keep_example = FALSE
  )
  row <- analysis$joint[analysis$joint$protein_group_id == sdc1_id, ]
  ami <- analysis$pooled_ami[analysis$pooled_ami$protein_group_id == sdc1_id, ]
  hc <- analysis$pooled_hc[analysis$pooled_hc$protein_group_id == sdc1_id, ]
  sensitivity_results[[i]] <- data.frame(
    shift = sensitivity_grid$shift[i],
    scale = sensitivity_grid$scale[i],
    imputations = M_sensitivity,
    sdc1_log2FC_aad_vs_ami = ami$pooled_log2FC,
    sdc1_p_aad_vs_ami = ami$p_value,
    sdc1_log2FC_aad_vs_healthy = hc$pooled_log2FC,
    sdc1_p_aad_vs_healthy = hc$p_value,
    sdc1_joint_iut_p = row$joint_iut_p,
    sdc1_joint_rank = row$joint_rank,
    top_candidate = analysis$joint$gene_symbol[1],
    stringsAsFactors = FALSE
  )
}
sensitivity_results <- do.call(rbind, sensitivity_results)
write_csv(sensitivity_results, file.path(source_dir, "sdc1_imputation_parameter_sensitivity_v7.csv"))

# Full-pipeline stratified bootstrap with a new stochastic imputation draw per resample.
set.seed(2026071801)
all_ids <- raw$protein_group_id
counter <- data.frame(
  protein_group_id = all_ids,
  retained = 0L,
  aad_higher_both = 0L,
  rank1 = 0L,
  top5 = 0L,
  top10 = 0L,
  top100 = 0L,
  stringsAsFactors = FALSE
)
rownames(counter) <- counter$protein_group_id
bootstrap_rows <- vector("list", B_boot)
group_indices <- split(seq_len(nrow(raw$sample_metadata)), raw$sample_metadata$group)

for (b in seq_len(B_boot)) {
  sampled <- unlist(lapply(group_indices, function(index) sample(index, length(index), replace = TRUE)), use.names = FALSE)
  boot_meta <- raw$sample_metadata[sampled, , drop = FALSE]
  boot_matrix_full <- raw_log2[, sampled, drop = FALSE]
  colnames(boot_matrix_full) <- paste0(boot_meta$sample_id, "_boot", seq_along(sampled))
  keep <- outcome_blind_filter(boot_matrix_full, minimum_detection = 0.50)
  boot_matrix <- boot_matrix_full[keep, , drop = FALSE]
  boot_lookup <- data.frame(
    protein_group_id = raw$protein_group_id[keep],
    gene_symbol = raw$gene_symbol[keep],
    stringsAsFactors = FALSE
  )
  normalised <- median_normalise(boot_matrix)$matrix
  imputed <- draw_left_censored(normalised, shift = primary_shift, scale = primary_scale)$matrix
  ami <- fit_adjusted_limma(imputed, boot_meta, "AAD", "AMI")
  hc <- fit_adjusted_limma(imputed, boot_meta, "AAD", "Healthy control")
  names(ami)[names(ami) == "beta"] <- "beta_ami"
  names(hc)[names(hc) == "beta"] <- "beta_hc"
  joint <- merge(ami[, c("protein_group_id", "beta_ami", "se2")], hc[, c("protein_group_id", "beta_hc", "se2")], by = "protein_group_id", suffixes = c("_ami", "_hc"), sort = FALSE)
  z_ami <- joint$beta_ami / sqrt(joint$se2_ami)
  z_hc <- joint$beta_hc / sqrt(joint$se2_hc)
  joint$p_ami_greater <- pnorm(z_ami, lower.tail = FALSE)
  joint$p_hc_greater <- pnorm(z_hc, lower.tail = FALSE)
  joint$joint_p <- pmax(joint$p_ami_greater, joint$p_hc_greater)
  joint$minimum_effect <- pmin(joint$beta_ami, joint$beta_hc)
  joint <- joint[order(joint$joint_p, -joint$minimum_effect, joint$protein_group_id), , drop = FALSE]
  joint$rank <- seq_len(nrow(joint))
  joint_ids <- joint$protein_group_id
  counter[joint_ids, "retained"] <- counter[joint_ids, "retained"] + 1L
  higher_ids <- joint$protein_group_id[joint$beta_ami > 0 & joint$beta_hc > 0]
  counter[higher_ids, "aad_higher_both"] <- counter[higher_ids, "aad_higher_both"] + 1L
  counter[joint_ids[joint$rank == 1], "rank1"] <- counter[joint_ids[joint$rank == 1], "rank1"] + 1L
  counter[joint_ids[joint$rank <= 5], "top5"] <- counter[joint_ids[joint$rank <= 5], "top5"] + 1L
  counter[joint_ids[joint$rank <= 10], "top10"] <- counter[joint_ids[joint$rank <= 10], "top10"] + 1L
  counter[joint_ids[joint$rank <= 100], "top100"] <- counter[joint_ids[joint$rank <= 100], "top100"] + 1L

  sdc1 <- joint[joint$protein_group_id == sdc1_id, , drop = FALSE]
  bootstrap_rows[[b]] <- data.frame(
    bootstrap = b,
    proteins_retained = nrow(joint),
    sdc1_retained = nrow(sdc1) == 1,
    sdc1_rank = if (nrow(sdc1)) sdc1$rank else NA_integer_,
    sdc1_log2FC_aad_vs_ami = if (nrow(sdc1)) sdc1$beta_ami else NA_real_,
    sdc1_log2FC_aad_vs_healthy = if (nrow(sdc1)) sdc1$beta_hc else NA_real_,
    sdc1_aad_higher_both = if (nrow(sdc1)) sdc1$beta_ami > 0 && sdc1$beta_hc > 0 else FALSE,
    top_candidate_id = joint$protein_group_id[1],
    top_candidate_gene = boot_lookup$gene_symbol[match(joint$protein_group_id[1], boot_lookup$protein_group_id)],
    stringsAsFactors = FALSE
  )
  if (b %% 25 == 0 || b == B_boot) cat("Bootstrap", b, "of", B_boot, "\n")
}

bootstrap_draws <- do.call(rbind, bootstrap_rows)
write_csv(bootstrap_draws, file.path(source_dir, "imputation_full_pipeline_bootstrap_draws_v7.csv"))

counter$gene_symbol <- raw$gene_symbol[match(counter$protein_group_id, raw$protein_group_id)]
counter$retained_frequency <- counter$retained / B_boot
counter$aad_higher_both_frequency <- counter$aad_higher_both / pmax(counter$retained, 1)
counter$rank1_frequency <- counter$rank1 / B_boot
counter$top5_frequency <- counter$top5 / B_boot
counter$top10_frequency <- counter$top10 / B_boot
counter$top100_frequency <- counter$top100 / B_boot
write_csv(counter, file.path(source_dir, "imputation_bootstrap_candidate_frequencies_v7.csv"))

sdc1_bootstrap_summary <- data.frame(
  protein_group_id = sdc1_id,
  gene_symbol = "SDC1",
  bootstrap_repetitions = B_boot,
  retained_count = sum(bootstrap_draws$sdc1_retained),
  retained_frequency = mean(bootstrap_draws$sdc1_retained),
  aad_higher_both_frequency_all = mean(bootstrap_draws$sdc1_aad_higher_both),
  aad_higher_both_frequency_retained = mean(bootstrap_draws$sdc1_aad_higher_both[bootstrap_draws$sdc1_retained]),
  rank1_frequency_all = sum(bootstrap_draws$sdc1_rank == 1, na.rm = TRUE) / B_boot,
  rank1_frequency_retained = mean(bootstrap_draws$sdc1_rank == 1, na.rm = TRUE),
  top5_frequency_all = sum(bootstrap_draws$sdc1_rank <= 5, na.rm = TRUE) / B_boot,
  top5_frequency_retained = mean(bootstrap_draws$sdc1_rank <= 5, na.rm = TRUE),
  top10_frequency_all = sum(bootstrap_draws$sdc1_rank <= 10, na.rm = TRUE) / B_boot,
  top10_frequency_retained = mean(bootstrap_draws$sdc1_rank <= 10, na.rm = TRUE),
  top100_frequency_all = sum(bootstrap_draws$sdc1_rank <= 100, na.rm = TRUE) / B_boot,
  top100_frequency_retained = mean(bootstrap_draws$sdc1_rank <= 100, na.rm = TRUE),
  median_rank = median(bootstrap_draws$sdc1_rank, na.rm = TRUE),
  rank_q025 = as.numeric(quantile(bootstrap_draws$sdc1_rank, 0.025, na.rm = TRUE)),
  rank_q975 = as.numeric(quantile(bootstrap_draws$sdc1_rank, 0.975, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write_csv(sdc1_bootstrap_summary, file.path(result_dir, "sdc1_imputation_bootstrap_stability_summary_v7.csv"))

sdc1_primary_summary <- rbind(
  data.frame(
    contrast = "AAD versus AMI",
    pooled_log2FC = primary_sdc1_ami$pooled_log2FC,
    ci_lower_95 = primary_sdc1_ami$ci_lower_95,
    ci_upper_95 = primary_sdc1_ami$ci_upper_95,
    p_value = primary_sdc1_ami$p_value,
    bh_fdr = primary_sdc1_ami$bh_fdr,
    fraction_missing_information = primary_sdc1_ami$fraction_missing_information,
    stringsAsFactors = FALSE
  ),
  data.frame(
    contrast = "AAD versus Healthy control",
    pooled_log2FC = primary_sdc1_hc$pooled_log2FC,
    ci_lower_95 = primary_sdc1_hc$ci_lower_95,
    ci_upper_95 = primary_sdc1_hc$ci_upper_95,
    p_value = primary_sdc1_hc$p_value,
    bh_fdr = primary_sdc1_hc$bh_fdr,
    fraction_missing_information = primary_sdc1_hc$fraction_missing_information,
    stringsAsFactors = FALSE
  )
)
sdc1_primary_summary$joint_iut_p <- primary_sdc1_joint$joint_iut_p
sdc1_primary_summary$joint_iut_bh <- primary_sdc1_joint$joint_iut_bh
sdc1_primary_summary$joint_rank <- primary_sdc1_joint$joint_rank
sdc1_primary_summary$imputations <- M_primary
sdc1_primary_summary$shift <- primary_shift
sdc1_primary_summary$scale <- primary_scale
write_csv(sdc1_primary_summary, file.path(result_dir, "sdc1_primary_imputation_summary_v7.csv"))

# Figure data.
plot_ami <- primary$pooled_ami
plot_hc <- primary$pooled_hc
joint_top <- head(primary$joint, 15)
highlight_ids <- joint_top$protein_group_id
plot_ami$display_group <- ifelse(plot_ami$protein_group_id == sdc1_id, "SDC1", ifelse(plot_ami$protein_group_id %in% highlight_ids, "Joint top 15", "Other proteins"))
plot_hc$display_group <- ifelse(plot_hc$protein_group_id == sdc1_id, "SDC1", ifelse(plot_hc$protein_group_id %in% highlight_ids, "Joint top 15", "Other proteins"))
plot_ami$negative_log10_p <- -log10(pmax(plot_ami$p_value, .Machine$double.xmin))
plot_hc$negative_log10_p <- -log10(pmax(plot_hc$p_value, .Machine$double.xmin))

palette <- c(
  AAD = "#BD3E4B", AMI = "#D8902E", Healthy = "#6F918B",
  SDC1 = "#BD3E4B", Joint = "#D8902E", Other = "#C9CDD1",
  Teal = "#2F8A7C", Blue = "#3677A8", Dark = "#202428"
)

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      plot.title = element_text(size = 8, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 6.4, colour = "#555555", margin = margin(b = 3)),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6.4, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 6.2),
      strip.text = element_text(size = 6.8, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

display_values <- c("Other proteins" = palette["Other"], "Joint top 15" = palette["Joint"], "SDC1" = palette["SDC1"])

p_a <- ggplot(detect_table, aes(x = factor(group, levels = c("Healthy control", "AMI", "AAD")), y = detection_fraction, fill = group)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.28) +
  geom_text(aes(label = paste0(detected, "/", total)), vjust = -0.35, size = 2.4) +
  scale_fill_manual(values = c("AAD" = palette["AAD"], "AMI" = palette["AMI"], "Healthy control" = palette["Healthy"])) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1.10), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "SDC1 detection", title = "Observed SDC1 detection") +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

make_volcano <- function(data, title) {
  label <- data[data$protein_group_id == sdc1_id, , drop = FALSE]
  ggplot(data, aes(pooled_log2FC, negative_log10_p, colour = display_group)) +
    geom_point(size = 0.65, alpha = 0.65) +
    geom_point(data = label, size = 2.1, colour = palette["SDC1"]) +
    geom_text_repel(data = label, aes(label = gene_symbol), size = 2.5, colour = palette["Dark"], min.segment.length = 0, seed = 17) +
    scale_colour_manual(values = display_values) +
    labs(x = "Pooled adjusted log2 fold change", y = "-log10 pooled P", title = title, subtitle = "50 stochastic imputations") +
    theme_pub() + theme(legend.position = "bottom")
}

p_b <- make_volcano(plot_ami, "AAD versus AMI")
p_c <- make_volcano(plot_hc, "AAD versus healthy controls")

joint_top$gene_label <- ifelse(is.na(joint_top$gene_symbol) | !nzchar(joint_top$gene_symbol), joint_top$protein_group_id, joint_top$gene_symbol)
joint_top$gene_label <- factor(joint_top$gene_label, levels = rev(joint_top$gene_label))
joint_top$joint_signal <- -log10(pmax(joint_top$joint_iut_p, .Machine$double.xmin))
p_d <- ggplot(joint_top, aes(joint_signal, gene_label, colour = protein_group_id == sdc1_id)) +
  geom_segment(aes(x = 0, xend = joint_signal, yend = gene_label), linewidth = 0.45, colour = "#C8CCD0") +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(`FALSE` = palette["Teal"], `TRUE` = palette["SDC1"]), guide = "none") +
  labs(x = "-log10 joint intersection-union P", y = NULL, title = "Cross-comparator candidate ranking", subtitle = "Continuous ranking; no stepwise selection") +
  theme_pub()

sdc1_forest <- sdc1_primary_summary
sdc1_forest$contrast <- factor(sdc1_forest$contrast, levels = rev(c("AAD versus AMI", "AAD versus Healthy control")))
p_e <- ggplot(sdc1_forest, aes(pooled_log2FC, contrast)) +
  geom_vline(xintercept = 0, colour = "#B8BDC2", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_lower_95, xmax = ci_upper_95), orientation = "y", width = 0.18, linewidth = 0.65, colour = palette["Dark"]) +
  geom_point(size = 2.8, colour = palette["SDC1"]) +
  labs(x = "Pooled adjusted log2 fold change (95% CI)", y = NULL, title = "SDC1 pooled effects", subtitle = paste0("Joint rank = ", primary_sdc1_joint$joint_rank)) +
  theme_pub()

stability_plot <- data.frame(
  measure = factor(c("Retained", "AAD higher", "Rank 1", "Top 5", "Top 10", "Top 100"), levels = c("Retained", "AAD higher", "Rank 1", "Top 5", "Top 10", "Top 100")),
  value = c(
    sdc1_bootstrap_summary$retained_frequency,
    sdc1_bootstrap_summary$aad_higher_both_frequency_all,
    sdc1_bootstrap_summary$rank1_frequency_all,
    sdc1_bootstrap_summary$top5_frequency_all,
    sdc1_bootstrap_summary$top10_frequency_all,
    sdc1_bootstrap_summary$top100_frequency_all
  )
)
p_f <- ggplot(stability_plot, aes(measure, value, fill = measure)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = paste0(round(100 * value), "%")), vjust = -0.35, size = 2.3) +
  scale_fill_manual(values = rep(c(palette["Teal"], palette["Blue"], palette["SDC1"]), length.out = 6), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1.10), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Frequency across all resamples", title = "Full-pipeline SDC1 stability", subtitle = paste0("Rank 1 in ", sum(bootstrap_draws$sdc1_rank == 1, na.rm = TRUE), "/", B_boot, " resamples")) +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))

main_figure <- (p_a | p_b | p_c) / (p_d | p_e | p_f) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 9, family = "Arial"))) +
  plot_layout(heights = c(1, 1.02))

save_figure <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(paste0(stem, ".pdf"), plot, width = width_in, height = height_in, device = cairo_pdf, family = "Arial")
  ggsave(paste0(stem, ".svg"), plot, width = width_in, height = height_in, device = svglite::svglite)
  ggsave(paste0(stem, ".png"), plot, width = width_in, height = height_in, dpi = 300, device = ragg::agg_png)
  ggsave(paste0(stem, ".tiff"), plot, width = width_in, height = height_in, dpi = dpi, compression = "lzw", device = ragg::agg_tiff)
}

save_figure(main_figure, file.path(main_figure_dir, "Figure_V7_transparent_imputation_candidate_selection"), 183, 145)

# Supplementary diagnostics.
example_long <- rbind(
  data.frame(group = sdc1_example$group, status = "Observed", value = sdc1_example$normalised_log2_value, stringsAsFactors = FALSE),
  data.frame(group = sdc1_example$group, status = "Imputed where missing", value = ifelse(sdc1_example$observed, NA, sdc1_example$example_imputed_log2_value), stringsAsFactors = FALSE)
)
example_long <- example_long[is.finite(example_long$value), ]
p_s1 <- ggplot(example_long, aes(value, fill = status, colour = status)) +
  geom_density(alpha = 0.28, linewidth = 0.55) +
  facet_wrap(~factor(group, levels = c("Healthy control", "AMI", "AAD")), scales = "free_y") +
  scale_fill_manual(values = c("Observed" = palette["Blue"], "Imputed where missing" = palette["SDC1"])) +
  scale_colour_manual(values = c("Observed" = palette["Blue"], "Imputed where missing" = palette["SDC1"])) +
  labs(x = "Normalised log2 SDC1", y = "Density", title = "Observed and example imputed SDC1 values", subtitle = "One stochastic draw shown for diagnostics only") +
  theme_pub() + theme(legend.position = "bottom")

p_s2 <- ggplot(primary$first_diagnostics, aes(cutoff_1pct, imputation_mean)) +
  geom_abline(slope = 1, intercept = 0, colour = "#AEB4B9", linewidth = 0.45) +
  geom_point(aes(size = missing_n, colour = robust_scale), alpha = 0.7) +
  scale_colour_gradient(low = palette["Teal"], high = palette["SDC1"]) +
  labs(x = "Sample 1% observed quantile", y = "Imputation distribution mean", size = "Missing proteins", colour = "Robust scale", title = "Sample-specific imputation parameters") +
  theme_pub() + theme(legend.position = "bottom")

p_s3 <- ggplot(primary$sdc1_draws, aes(contrast, log2FC, fill = contrast)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 0.75, alpha = 0.55) +
  scale_fill_manual(values = c("AAD versus AMI" = palette["AMI"], "AAD versus Healthy control" = palette["Healthy"])) +
  labs(x = NULL, y = "Adjusted log2 fold change", title = "Between-imputation SDC1 uncertainty") +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))

sensitivity_results$shift_label <- factor(paste0("shift ", sensitivity_results$shift), levels = paste0("shift ", c(1.5, 1.8, 2.1)))
sensitivity_results$scale_label <- factor(paste0("scale ", sensitivity_results$scale), levels = paste0("scale ", c(0.2, 0.3, 0.5)))
p_s4 <- ggplot(sensitivity_results, aes(scale_label, shift_label, fill = sdc1_joint_rank)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = paste0("rank ", sdc1_joint_rank)), size = 2.4) +
  scale_fill_gradient(low = "#F4F6F7", high = palette["SDC1"], trans = "reverse") +
  labs(x = NULL, y = NULL, fill = "SDC1 rank", title = "Pre-fixed imputation sensitivity grid") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_s5 <- ggplot(bootstrap_draws[is.finite(bootstrap_draws$sdc1_rank), ], aes(sdc1_rank)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = palette["Teal"], colour = "white", linewidth = 0.25) +
  scale_x_continuous(breaks = sort(unique(bootstrap_draws$sdc1_rank[is.finite(bootstrap_draws$sdc1_rank)]))) +
  labs(x = "SDC1 joint rank", y = "Bootstrap count", title = "Bootstrap SDC1 rank") +
  theme_pub()

supp_figure <- (p_s1 | p_s2) / (p_s3 | p_s4 | p_s5) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 9, family = "Arial"))) +
  plot_layout(heights = c(1, 1))
save_figure(supp_figure, file.path(supp_figure_dir, "Supplementary_Figure_V7_imputation_diagnostics_and_sensitivity"), 183, 150)

# Machine-readable run manifest and concise text log.
run_manifest <- data.frame(
  parameter = c(
    "analysis_date", "candidate_universe", "primary_imputations", "bootstrap_repetitions",
    "sensitivity_imputations", "primary_shift", "primary_scale", "candidate_rule"
  ),
  value = c(
    "2026-07-16", nrow(primary_matrix), M_primary, B_boot,
    M_sensitivity, primary_shift, primary_scale,
    "intersection-union P across AAD-vs-AMI and AAD-vs-healthy; continuous rank; no stepwise"
  ),
  stringsAsFactors = FALSE
)
write_csv(run_manifest, file.path(result_dir, "v7_imputation_run_manifest.csv"))

log_lines <- c(
  paste("Primary candidate universe:", nrow(primary_matrix)),
  paste("Primary top candidate:", primary$joint$gene_symbol[1], primary$joint$protein_group_id[1]),
  paste("SDC1 joint rank:", primary_sdc1_joint$joint_rank),
  paste("SDC1 joint IUT P:", format(primary_sdc1_joint$joint_iut_p, scientific = TRUE, digits = 5)),
  paste("SDC1 AAD-vs-AMI pooled log2FC:", round(primary_sdc1_ami$pooled_log2FC, 4)),
  paste("SDC1 AAD-vs-healthy pooled log2FC:", round(primary_sdc1_hc$pooled_log2FC, 4)),
  paste("SDC1 bootstrap retained frequency:", round(sdc1_bootstrap_summary$retained_frequency, 4)),
  paste("SDC1 bootstrap rank-1 frequency among all resamples:", round(sdc1_bootstrap_summary$rank1_frequency_all, 4)),
  paste("SDC1 bootstrap rank-1 frequency conditional on retention:", round(sdc1_bootstrap_summary$rank1_frequency_retained, 4)),
  paste("SDC1 bootstrap top-10 frequency among all resamples:", round(sdc1_bootstrap_summary$top10_frequency_all, 4)),
  paste("SDC1 bootstrap median rank:", sdc1_bootstrap_summary$median_rank),
  paste("Sensitivity SDC1 rank range:", paste(range(sensitivity_results$sdc1_joint_rank), collapse = " to "))
)
writeLines(log_lines, file.path(log_dir, "v7_imputation_summary.txt"), useBytes = TRUE)
cat(paste(log_lines, collapse = "\n"), "\n")

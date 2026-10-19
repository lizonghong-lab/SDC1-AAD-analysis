options(stringsAsFactors = FALSE, scipen = 999)

# Supplementary analysis contract
# Question: is SDC1 merely a proxy for a strongly collinear protein signal?
# Role: outcome-blind redundancy audit, not disease-guided candidate selection.
# Universe: the same 1,837 protein groups retained by the locked >=50%
# detection and non-zero-IQR rule.
# Missingness: 50 stochastic, sample-specific left-censoring imputations using
# the locked primary parameters (shift 1.8; scale 0.3).
# Estimand: Fisher-z pooled Spearman correlation between SDC1 and every other
# eligible protein; percentile ranges describe imputation sensitivity only.

set.seed(2026072204)

revision_root <- normalizePath(
  Sys.getenv(
    "SDC1_REVISION_PACKAGE",
    unset = file.path(getwd(), "SDC1_AAD_REVISED_SUBMISSION_FINAL_20260722")
  ),
  winslash = "/",
  mustWork = TRUE
)
source_root <- normalizePath(
  file.path(getwd(), "SDC1_AAD_REANALYSIS_V6_1_REVIEWER_REFINED_20260715"),
  winslash = "/",
  mustWork = TRUE
)
ranking_path <- normalizePath(
  file.path(getwd(), "SDC1_AAD_REANALYSIS_V7_TRANSPARENT_IMPUTATION_20260716", "04_results", "source_data", "imputation_joint_candidate_ranking_v7.csv"),
  winslash = "/",
  mustWork = TRUE
)

output_dir <- file.path(revision_root, "04_figures", "supplementary", "source_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
      output[missing, j] <- qnorm(
        runif(length(missing), min = 1e-10, max = upper_probability),
        mean = mu,
        sd = sigma
      )
    }
  }
  output
}

raw <- readRDS(file.path(source_root, "02_analysis_data", "proteomics_raw_deidentified.rds"))
log2_matrix <- log2(raw$intensity)
rownames(log2_matrix) <- raw$protein_group_id
colnames(log2_matrix) <- raw$sample_metadata$sample_id

keep <- outcome_blind_filter(log2_matrix)
analysis_matrix <- median_normalise(log2_matrix[keep, , drop = FALSE])
lookup <- data.frame(
  protein_group_id = raw$protein_group_id[keep],
  gene_symbol = raw$gene_symbol[keep],
  overall_detection_fraction = rowMeans(is.finite(log2_matrix[keep, , drop = FALSE])),
  stringsAsFactors = FALSE
)

sdc1_index <- which(toupper(lookup$gene_symbol) == "SDC1")
if (length(sdc1_index) != 1L) stop("Expected exactly one eligible SDC1 row")

M <- 50L
correlation_draws <- matrix(NA_real_, nrow = nrow(analysis_matrix), ncol = M)
for (m in seq_len(M)) {
  imputed <- draw_left_censored(analysis_matrix, shift = 1.8, scale = 0.3)
  correlation_draws[, m] <- apply(
    imputed,
    1,
    function(x) suppressWarnings(cor(x, imputed[sdc1_index, ], method = "spearman"))
  )
  if (m %% 10L == 0L) message("Correlation imputation ", m, "/", M)
}

clip_rho <- function(x) pmin(pmax(x, -0.999999), 0.999999)
pooled_rho <- tanh(rowMeans(atanh(clip_rho(correlation_draws))))
rho_lower <- apply(correlation_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
rho_upper <- apply(correlation_draws, 1, quantile, probs = 0.975, na.rm = TRUE)

ranking <- read.csv(ranking_path, stringsAsFactors = FALSE, check.names = FALSE)
audit <- merge(
  lookup,
  ranking[, c(
    "protein_group_id", "joint_rank", "pooled_log2FC_ami",
    "pooled_log2FC_healthy", "joint_iut_p"
  )],
  by = "protein_group_id",
  all.x = TRUE,
  sort = FALSE
)
audit$rho_with_sdc1 <- pooled_rho[match(audit$protein_group_id, lookup$protein_group_id)]
audit$rho_imputation_p025 <- rho_lower[match(audit$protein_group_id, lookup$protein_group_id)]
audit$rho_imputation_p975 <- rho_upper[match(audit$protein_group_id, lookup$protein_group_id)]
audit$absolute_rho_with_sdc1 <- abs(audit$rho_with_sdc1)
audit <- audit[order(-audit$absolute_rho_with_sdc1, audit$protein_group_id), , drop = FALSE]
audit$absolute_correlation_rank <- seq_len(nrow(audit))
audit$absolute_correlation_rank[audit$protein_group_id == lookup$protein_group_id[sdc1_index]] <- NA_integer_

write.csv(
  audit,
  file.path(output_dir, "sdc1_outcome_blind_correlation_audit_all_proteins.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

other <- audit[toupper(audit$gene_symbol) != "SDC1", , drop = FALSE]
top_neighbours <- head(other, 12)
write.csv(
  top_neighbours,
  file.path(output_dir, "sdc1_outcome_blind_correlation_neighbourhood.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

top_ranked <- ranking[ranking$joint_rank <= 15, c("protein_group_id", "gene_symbol", "joint_rank")]
top_ranked <- merge(
  top_ranked,
  audit[, c(
    "protein_group_id", "rho_with_sdc1", "rho_imputation_p025",
    "rho_imputation_p975", "overall_detection_fraction"
  )],
  by = "protein_group_id",
  all.x = TRUE,
  sort = FALSE
)
top_ranked <- top_ranked[order(top_ranked$joint_rank), ]
write.csv(
  top_ranked,
  file.path(output_dir, "sdc1_top15_candidate_correlation_audit.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

summary_table <- data.frame(
  item = c(
    "Eligible protein groups",
    "Imputation draws",
    "Largest absolute SDC1 correlation excluding SDC1",
    "Proteins with absolute rho >= 0.70",
    "Proteins with absolute rho >= 0.50",
    "Proteins with absolute rho >= 0.30",
    "Top correlated protein",
    "Top correlated protein joint candidate rank"
  ),
  value = c(
    nrow(audit),
    M,
    sprintf("%.3f", max(other$absolute_rho_with_sdc1, na.rm = TRUE)),
    sum(other$absolute_rho_with_sdc1 >= 0.70, na.rm = TRUE),
    sum(other$absolute_rho_with_sdc1 >= 0.50, na.rm = TRUE),
    sum(other$absolute_rho_with_sdc1 >= 0.30, na.rm = TRUE),
    top_neighbours$gene_symbol[1],
    top_neighbours$joint_rank[1]
  ),
  stringsAsFactors = FALSE
)
write.csv(
  summary_table,
  file.path(output_dir, "sdc1_outcome_blind_correlation_summary.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

message("Wrote outcome-blind SDC1 correlation audit to ", output_dir)

options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
revision_dir <- root
project_library <- file.path(revision_dir, "08_environment", "Rlib")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))

suppressPackageStartupMessages({
  library(limma)
})

set.seed(20260714)

analysis_dir <- file.path(revision_dir, "02_analysis_data")
results_dir <- file.path(revision_dir, "04_results", "statistical_results")
source_dir <- file.path(revision_dir, "04_results", "source_data")
logs_dir <- file.path(revision_dir, "04_results", "logs")
for (directory in c(results_dir, source_dir, logs_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

raw <- readRDS(file.path(analysis_dir, "proteomics_raw_deidentified.rds"))
intensity <- raw$intensity
sample_metadata <- raw$sample_metadata
protein_group_id <- raw$protein_group_id
gene_symbol <- raw$gene_symbol

stopifnot(identical(colnames(intensity), sample_metadata$sample_id))
stopifnot(length(protein_group_id) == nrow(intensity))

log2_intensity <- log2(intensity)
sample_median <- apply(log2_intensity, 2, median, na.rm = TRUE)
sample_iqr <- apply(log2_intensity, 2, IQR, na.rm = TRUE)
sample_detected <- colSums(!is.na(log2_intensity))

sample_qc <- cbind(
  sample_metadata,
  data.frame(
    detected_protein_groups = sample_detected,
    detected_fraction = sample_detected / nrow(log2_intensity),
    median_log2_intensity = sample_median,
    iqr_log2_intensity = sample_iqr,
    robust_median_z = as.numeric((sample_median - median(sample_median)) / mad(sample_median)),
    robust_detection_z = as.numeric((sample_detected - median(sample_detected)) / mad(sample_detected)),
    stringsAsFactors = FALSE
  )
)

group_levels <- c("AAD", "AMI", "Healthy control")
group_index <- lapply(group_levels, function(group_name) which(sample_metadata$group == group_name))
names(group_index) <- group_levels

overall_detection <- rowMeans(!is.na(log2_intensity))
detection_by_group <- sapply(group_levels, function(group_name) {
  rowMeans(!is.na(log2_intensity[, group_index[[group_name]], drop = FALSE]))
})
colnames(detection_by_group) <- c("AAD", "AMI", "Healthy_control")

protein_qc <- data.frame(
  protein_group_id = protein_group_id,
  gene_symbol = gene_symbol,
  overall_detection_fraction = overall_detection,
  aad_detection_fraction = detection_by_group[, "AAD"],
  ami_detection_fraction = detection_by_group[, "AMI"],
  healthy_detection_fraction = detection_by_group[, "Healthy_control"],
  median_log2_intensity = apply(log2_intensity, 1, median, na.rm = TRUE),
  iqr_log2_intensity = apply(log2_intensity, 1, IQR, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# The primary inferential filter is outcome-blinded: a protein group must be
# quantified in at least half of all 307 samples, regardless of disease group.
primary_keep <- overall_detection >= 0.50 & is.finite(protein_qc$iqr_log2_intensity) & protein_qc$iqr_log2_intensity > 0

log2_primary <- log2_intensity[primary_keep, , drop = FALSE]
primary_protein_id <- protein_group_id[primary_keep]
primary_gene <- gene_symbol[primary_keep]

# Median normalisation corrects sample-wide loading shifts without using the
# outcome. Missing values remain missing for protein-level inference.
primary_sample_median <- apply(log2_primary, 2, median, na.rm = TRUE)
grand_median <- median(primary_sample_median, na.rm = TRUE)
log2_median_normalised <- sweep(log2_primary, 2, primary_sample_median - grand_median, "-")

# Quantile-normalised values are retained only for a sensitivity analysis.
log2_quantile_normalised <- normalizeBetweenArrays(log2_primary, method = "quantile")
rownames(log2_quantile_normalised) <- rownames(log2_primary)
colnames(log2_quantile_normalised) <- colnames(log2_primary)

# PCA uses a stricter outcome-blinded completeness rule. Protein-median
# imputation is used only to display unsupervised structure and is never used
# for differential-abundance inference.
pca_keep_within_primary <- rowMeans(!is.na(log2_median_normalised)) >= 0.80
pca_matrix <- log2_median_normalised[pca_keep_within_primary, , drop = FALSE]
pca_protein_id <- primary_protein_id[pca_keep_within_primary]
pca_gene <- primary_gene[pca_keep_within_primary]

for (i in seq_len(nrow(pca_matrix))) {
  replacement <- median(pca_matrix[i, ], na.rm = TRUE)
  pca_matrix[i, is.na(pca_matrix[i, ])] <- replacement
}
pca_scaled <- t(scale(t(pca_matrix)))
pca_finite <- apply(pca_scaled, 1, function(x) all(is.finite(x)))
pca_scaled <- pca_scaled[pca_finite, , drop = FALSE]
pca_protein_id <- pca_protein_id[pca_finite]
pca_gene <- pca_gene[pca_finite]

pca_fit <- prcomp(t(pca_scaled), center = FALSE, scale. = FALSE)
pca_variance <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
n_pc_90 <- which(cumsum(pca_variance) >= 90)[1]

pca_scores <- cbind(
  sample_metadata,
  data.frame(
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    PC3 = pca_fit$x[, 3],
    stringsAsFactors = FALSE
  )
)
pca_loadings <- data.frame(
  protein_group_id = pca_protein_id,
  gene_symbol = pca_gene,
  PC1 = pca_fit$rotation[, 1],
  PC2 = pca_fit$rotation[, 2],
  PC3 = pca_fit$rotation[, 3],
  stringsAsFactors = FALSE
)
pca_variance_table <- data.frame(
  principal_component = paste0("PC", seq_along(pca_variance)),
  percent_variance = pca_variance,
  cumulative_percent_variance = cumsum(pca_variance),
  stringsAsFactors = FALSE
)

safe_cor <- function(x, y, method = "spearman") {
  if (sum(complete.cases(x, y)) < 5 || length(unique(na.omit(x))) < 2 || length(unique(na.omit(y))) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, use = "complete.obs", method = method))
}

pca_technical <- data.frame(
  principal_component = c("PC1", "PC2", "PC3"),
  correlation_with_export_order = c(
    safe_cor(pca_scores$PC1, pca_scores$export_column_order),
    safe_cor(pca_scores$PC2, pca_scores$export_column_order),
    safe_cor(pca_scores$PC3, pca_scores$export_column_order)
  ),
  correlation_with_detected_fraction = c(
    safe_cor(pca_scores$PC1, sample_qc$detected_fraction),
    safe_cor(pca_scores$PC2, sample_qc$detected_fraction),
    safe_cor(pca_scores$PC3, sample_qc$detected_fraction)
  ),
  correlation_with_sample_median = c(
    safe_cor(pca_scores$PC1, sample_qc$median_log2_intensity),
    safe_cor(pca_scores$PC2, sample_qc$median_log2_intensity),
    safe_cor(pca_scores$PC3, sample_qc$median_log2_intensity)
  ),
  stringsAsFactors = FALSE
)

# Correlation-based clustering is outcome-blinded. Restricting to the 300 most
# variable high-completeness proteins keeps the structure legible and avoids
# allowing thousands of nearly invariant proteins to dominate the dendrogram.
protein_variance <- apply(pca_scaled, 1, var)
cluster_order <- order(protein_variance, decreasing = TRUE)
cluster_order <- cluster_order[seq_len(min(300, length(cluster_order)))]
cluster_matrix <- pca_scaled[cluster_order, , drop = FALSE]
cluster_cor <- cor(t(cluster_matrix), method = "spearman")
cluster_distance <- as.dist(1 - abs(cluster_cor))
cluster_fit <- hclust(cluster_distance, method = "average")
cluster_pca <- prcomp(t(cluster_matrix), center = FALSE, scale. = FALSE)
cluster_pca_variance <- cluster_pca$sdev^2 / sum(cluster_pca$sdev^2)
cluster_count <- which(cumsum(cluster_pca_variance) >= 0.80)[1]
cluster_membership <- cutree(cluster_fit, k = cluster_count)

module_assignment <- data.frame(
  protein_group_id = pca_protein_id[cluster_order],
  gene_symbol = pca_gene[cluster_order],
  variance_rank = seq_along(cluster_order),
  module = paste0("M", cluster_membership),
  stringsAsFactors = FALSE
)

module_ids <- sort(unique(cluster_membership))
module_scores <- sapply(module_ids, function(module_id) {
  members <- which(cluster_membership == module_id)
  if (length(members) == 1) return(as.numeric(cluster_matrix[members, ]))
  first_pc <- prcomp(t(cluster_matrix[members, , drop = FALSE]), center = TRUE, scale. = FALSE)$x[, 1]
  mean_loading <- mean(cluster_matrix[members, ] %*% first_pc)
  if (is.finite(mean_loading) && mean_loading < 0) first_pc <- -first_pc
  as.numeric(first_pc)
})
if (is.vector(module_scores)) module_scores <- matrix(module_scores, ncol = 1)
colnames(module_scores) <- paste0("M", module_ids)
rownames(module_scores) <- colnames(cluster_matrix)
module_scores <- data.frame(sample_id = rownames(module_scores), module_scores, check.names = FALSE)

module_analysis <- merge(sample_metadata, module_scores, by = "sample_id", all.x = TRUE, sort = FALSE)
module_analysis <- module_analysis[match(sample_metadata$sample_id, module_analysis$sample_id), ]
module_names <- setdiff(names(module_scores), "sample_id")

fit_module_contrast <- function(positive_group, reference_group) {
  d <- module_analysis[module_analysis$group %in% c(positive_group, reference_group), , drop = FALSE]
  complete <- complete.cases(d[, c(
    "group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l"
  )])
  d <- d[complete, , drop = FALSE]
  d$aad_indicator <- as.integer(d$group == positive_group)
  d$age_scaled <- as.numeric(scale(d$age_years))
  d$wbc_scaled <- as.numeric(scale(log1p(d$wbc_10e9_l)))
  d$crp_scaled <- as.numeric(scale(log1p(d$crp_mg_l)))
  design <- model.matrix(
    ~ aad_indicator + age_scaled + female + hypertension + wbc_scaled + crp_scaled,
    data = d
  )
  module_matrix <- t(as.matrix(d[, module_names, drop = FALSE]))
  fit <- eBayes(lmFit(module_matrix, design), trend = TRUE, robust = TRUE)
  table <- topTable(fit, coef = "aad_indicator", number = Inf, sort.by = "none")
  standard_error <- abs(table$logFC / table$t)
  standard_error[!is.finite(standard_error)] <- NA_real_
  data.frame(
    contrast = paste(positive_group, "versus", reference_group),
    module = rownames(table),
    adjusted_module_score_difference = table$logFC,
    ci_lower_95 = table$logFC - qnorm(0.975) * standard_error,
    ci_upper_95 = table$logFC + qnorm(0.975) * standard_error,
    p_value = table$P.Value,
    bh_fdr = table$adj.P.Val,
    n_positive = sum(d$group == positive_group),
    n_reference = sum(d$group == reference_group),
    stringsAsFactors = FALSE
  )
}

module_associations <- rbind(
  fit_module_contrast("AAD", "AMI"),
  fit_module_contrast("AAD", "Healthy control")
)

sdc1_index <- which(toupper(gene_symbol) == "SDC1")
if (length(sdc1_index) != 1) stop("Expected exactly one SDC1 protein group")
sdc1_id <- protein_group_id[sdc1_index]
sdc1_eligibility <- data.frame(
  protein_group_id = sdc1_id,
  gene_symbol = gene_symbol[sdc1_index],
  overall_detection_fraction = protein_qc$overall_detection_fraction[sdc1_index],
  aad_detection_fraction = protein_qc$aad_detection_fraction[sdc1_index],
  ami_detection_fraction = protein_qc$ami_detection_fraction[sdc1_index],
  healthy_detection_fraction = protein_qc$healthy_detection_fraction[sdc1_index],
  primary_inference_eligible = sdc1_id %in% primary_protein_id,
  pca_80_percent_completeness_eligible = sdc1_id %in% pca_protein_id,
  top300_module_eligible = sdc1_id %in% module_assignment$protein_group_id,
  interpretation = paste(
    "SDC1 is retained for the >=50% primary inference universe but is not eligible",
    "for the >=80% PCA or high-completeness module analysis; detection is modelled separately"
  ),
  stringsAsFactors = FALSE
)

group_mean <- function(matrix_value, group_name) {
  rowMeans(matrix_value[, sample_metadata$group == group_name, drop = FALSE], na.rm = TRUE)
}

global_shift <- data.frame(
  protein_group_id = primary_protein_id,
  gene_symbol = primary_gene,
  raw_aad_minus_ami = group_mean(log2_primary, "AAD") - group_mean(log2_primary, "AMI"),
  raw_aad_minus_healthy = group_mean(log2_primary, "AAD") - group_mean(log2_primary, "Healthy control"),
  median_normalised_aad_minus_ami = group_mean(log2_median_normalised, "AAD") - group_mean(log2_median_normalised, "AMI"),
  median_normalised_aad_minus_healthy = group_mean(log2_median_normalised, "AAD") - group_mean(log2_median_normalised, "Healthy control"),
  stringsAsFactors = FALSE
)

normalisation_summary <- data.frame(
  contrast = c("AAD minus AMI", "AAD minus healthy control"),
  raw_median_log2_difference = c(
    median(global_shift$raw_aad_minus_ami, na.rm = TRUE),
    median(global_shift$raw_aad_minus_healthy, na.rm = TRUE)
  ),
  raw_fraction_positive = c(
    mean(global_shift$raw_aad_minus_ami > 0, na.rm = TRUE),
    mean(global_shift$raw_aad_minus_healthy > 0, na.rm = TRUE)
  ),
  normalised_median_log2_difference = c(
    median(global_shift$median_normalised_aad_minus_ami, na.rm = TRUE),
    median(global_shift$median_normalised_aad_minus_healthy, na.rm = TRUE)
  ),
  normalised_fraction_positive = c(
    mean(global_shift$median_normalised_aad_minus_ami > 0, na.rm = TRUE),
    mean(global_shift$median_normalised_aad_minus_healthy > 0, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

candidate_flow <- data.frame(
  stage = c(
    "Protein groups in raw Spectronaut export",
    "Protein groups in locked analysis matrix",
    "Unique non-empty gene annotations",
    "Protein groups retained for primary inference (overall detection >=50%)",
    "Protein groups retained for unsupervised PCA (overall detection >=80%)",
    "Most-variable high-completeness protein groups used for clustering",
    "SDC1 protein groups in the locked matrix"
  ),
  n = c(
    3820,
    nrow(intensity),
    length(unique(gene_symbol[gene_symbol != "Unannotated"])),
    sum(primary_keep),
    nrow(pca_scaled),
    nrow(cluster_matrix),
    sum(toupper(gene_symbol) == "SDC1")
  ),
  rule = c(
    "Source export row count",
    "Post-processing matrix provided for analysis",
    "Exact annotation count",
    "Outcome-blinded completeness and non-zero variability",
    "Outcome-blinded completeness; protein-median imputation for PCA only",
    "Outcome-blinded variance rank",
    "Exact case-insensitive gene-symbol match"
  ),
  stringsAsFactors = FALSE
)

qc_object <- list(
  primary_protein_id = primary_protein_id,
  primary_gene = primary_gene,
  primary_keep = primary_keep,
  log2_primary = log2_primary,
  log2_median_normalised = log2_median_normalised,
  log2_quantile_normalised = log2_quantile_normalised,
  sample_metadata = sample_metadata,
  sample_qc = sample_qc,
  protein_qc = protein_qc,
  pca_fit = pca_fit,
  pca_scores = pca_scores,
  pca_loadings = pca_loadings,
  pca_variance = pca_variance_table,
  pca_scaled = pca_scaled,
  cluster_fit = cluster_fit,
  cluster_cor = cluster_cor,
  module_assignment = module_assignment,
  module_scores = module_scores,
  module_associations = module_associations,
  sdc1_unsupervised_eligibility = sdc1_eligibility,
  cluster_count_rule = "Number of outcome-blinded PCs explaining 80% of variance among the 300 most-variable high-completeness proteins",
  cluster_count = cluster_count,
  candidate_flow = candidate_flow,
  global_shift = global_shift,
  normalisation_summary = normalisation_summary,
  n_pc_90 = n_pc_90
)
saveRDS(qc_object, file.path(analysis_dir, "proteomics_qc_reduction.rds"), compress = "xz")

write_csv(sample_qc, file.path(source_dir, "proteomics_sample_qc.csv"))
write_csv(protein_qc, file.path(source_dir, "proteomics_detection_missingness.csv"))
write_csv(pca_scores, file.path(source_dir, "proteomics_pca_scores.csv"))
write_csv(pca_loadings, file.path(source_dir, "proteomics_pca_loadings.csv"))
write_csv(pca_variance_table, file.path(source_dir, "proteomics_pca_variance.csv"))
write_csv(pca_technical, file.path(results_dir, "proteomics_pca_technical_associations.csv"))
write_csv(module_assignment, file.path(source_dir, "proteomics_unsupervised_modules.csv"))
write_csv(module_scores, file.path(source_dir, "proteomics_module_scores.csv"))
write_csv(module_associations, file.path(results_dir, "proteomics_module_associations_v4.csv"))
write_csv(sdc1_eligibility, file.path(results_dir, "sdc1_unsupervised_eligibility_v4.csv"))
write_csv(global_shift, file.path(source_dir, "proteomics_global_shift.csv"))
write_csv(normalisation_summary, file.path(results_dir, "proteomics_normalisation_summary.csv"))
write_csv(candidate_flow, file.path(results_dir, "proteomics_candidate_flow.csv"))

capture.output(sessionInfo(), file = file.path(logs_dir, "01_sessionInfo.txt"))

cat("Proteomic QC complete\n")
cat("Locked matrix:", nrow(intensity), "protein groups x", ncol(intensity), "samples\n")
cat("Primary outcome-blinded filter retained:", sum(primary_keep), "protein groups\n")
cat("PCA retained:", nrow(pca_scaled), "protein groups\n")
cat("PC1/PC2 variance:", sprintf("%.1f%% / %.1f%%", pca_variance[1], pca_variance[2]), "\n")
cat("PCs for 90% variance:", n_pc_90, "\n")
cat("Unsupervised modules among top-variable proteins:", length(module_ids), "\n")
cat("SDC1 eligible for PCA:", sdc1_eligibility$pca_80_percent_completeness_eligible, "\n")
print(normalisation_summary)

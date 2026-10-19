options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
suppressPackageStartupMessages(library(readxl))
set.seed(20260714)

analysis_dir <- file.path(root, "02_analysis_data")
results_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
logs_dir <- file.path(root, "04_results", "logs")
for (directory in c(analysis_dir, results_dir, source_dir, logs_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

read_canonical <- function(name) {
  read.csv(file.path(analysis_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
}

discovery_raw <- read_canonical("discovery_clinical_v5.csv")
validation_raw <- read_canonical("validation_clinical_v5.csv")
pilot_raw <- read_canonical("pilot_multicompartment_sdc1_v4.csv")

discovery <- within(discovery_raw, {
  study_id <- source_id
  ev_sdc1 <- as.numeric(ev_sdc1)
  d_dimer_mg_l <- as.numeric(d_dimer_mg_l)
  age_years <- as.numeric(age_years)
  female <- as.numeric(female)
  hypertension <- as.numeric(hypertension)
  wbc_10e9_l <- as.numeric(wbc_10e9_l)
  crp_mg_l <- as.numeric(crp_mg_l)
  creatinine_umol_l <- as.numeric(creatinine_umol_l)
  add_rs <- as.numeric(add_rs_corrected)
})
if (anyDuplicated(discovery$study_id)) stop("Discovery source IDs must be unique")
if (any(!is.finite(discovery$ev_sdc1) | discovery$ev_sdc1 <= 0)) stop("Discovery AF4-EV SDC1 values must be positive")
write_csv(discovery, file.path(analysis_dir, "discovery_clinical_deidentified.csv"))

proteomics_path <- file.path(
  root, "01_raw_data_locked", "latest",
  "Plasma_EV_proteomics_Discovery_Cohort_latest.xlsx"
)
proteomics_xlsx <- as.data.frame(read_excel(proteomics_path, sheet = 1, .name_repair = "minimal"))
gene_symbol <- trimws(as.character(proteomics_xlsx[[1]]))
gene_symbol[is.na(gene_symbol) | !nzchar(gene_symbol)] <- "Unannotated"
protein_group_id <- sprintf("PG%04d", seq_len(nrow(proteomics_xlsx)))
intensity <- as.matrix(proteomics_xlsx[, -1, drop = FALSE])
storage.mode(intensity) <- "numeric"
intensity[!is.finite(intensity) | intensity <= 0] <- NA_real_
sample_ids <- colnames(intensity)
if (anyDuplicated(sample_ids)) stop("Duplicate proteomic sample IDs")
if (!setequal(sample_ids, discovery$study_id)) stop("Discovery clinical and proteomic IDs do not match")

sample_metadata <- discovery[match(sample_ids, discovery$study_id), c(
  "study_id", "group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l"
)]
names(sample_metadata)[1] <- "sample_id"
sample_metadata$analysis_column_order <- seq_along(sample_ids)
sample_metadata$export_column_order <- seq_along(sample_ids)
sample_metadata$acquisition_date <- as.Date(NA)

proteomics_object <- list(
  protein_group_id = protein_group_id,
  gene_symbol = gene_symbol,
  intensity = intensity,
  sample_metadata = sample_metadata,
  source_sha256 = read.csv(
    file.path(results_dir, "input_manifest_v4.csv"), stringsAsFactors = FALSE
  )$sha256[3]
)
saveRDS(proteomics_object, file.path(analysis_dir, "proteomics_raw_deidentified.rds"), compress = "xz")

sdc1_row <- which(toupper(gene_symbol) == "SDC1")
if (length(sdc1_row) != 1L) stop("Expected exactly one SDC1 row")
sdc1_detected <- is.finite(intensity[sdc1_row, ])
sdc1_detection_audit <- do.call(rbind, lapply(split(seq_along(sample_ids), sample_metadata$group), function(index) {
  data.frame(
    group = sample_metadata$group[index[1]],
    n = length(index),
    detected = sum(sdc1_detected[index]),
    detection_fraction = mean(sdc1_detected[index]),
    stringsAsFactors = FALSE
  )
}))
write_csv(sdc1_detection_audit, file.path(results_dir, "sdc1_proteomic_raw_detection_audit_v4.csv"))

validation <- within(validation_raw, {
  study_id <- row_uid
  cluster_id <- ifelse(is.na(source_id) | !nzchar(source_id), row_uid, paste0("SRC_", source_id))
  sdc1_source_value <- as.numeric(total_sdc1)
  log2_sdc1 <- log2(sdc1_source_value)
  log2_d_dimer <- log2(as.numeric(d_dimer_mg_l))
  add_rs <- as.numeric(add_rs_corrected)
})
if (anyDuplicated(validation$study_id)) stop("V4 row UIDs must be unique")
if (any(!is.finite(validation$sdc1_source_value) | validation$sdc1_source_value <= 0)) {
  stop("Validation total-plasma SDC1 values must be positive")
}
write_csv(validation, file.path(analysis_dir, "elisa_total_sdc1_all.csv"))

elisa_primary <- validation[
  validation$group %in% c("AAD", "AMI", "Angina", "PE") &
    complete.cases(validation[, c("sdc1_source_value", "d_dimer_mg_l")]) &
    validation$d_dimer_mg_l > 0,
]
elisa_primary$symptomatic_comparator <- as.integer(elisa_primary$group != "AAD")
if (nrow(elisa_primary) != 329L) stop("Expected 329 symptomatic validation records")
write_csv(elisa_primary, file.path(analysis_dir, "elisa_evaluation_primary.csv"))

write_csv(pilot_raw, file.path(analysis_dir, "pilot_multicompartment_sdc1_analysis.csv"))

cohort_counts <- do.call(rbind, lapply(list(
  discovery = discovery,
  validation = validation,
  pilot_validation_subset = pilot_raw
), function(x) {
  do.call(rbind, lapply(split(x, x$group), function(y) {
    data.frame(
      analysis_layer = unique(y$cohort)[1],
      group = unique(y$group)[1],
      records = nrow(y),
      nonmissing_source_ids = sum(!is.na(y$source_id) & nzchar(y$source_id)),
      unique_source_ids = length(unique(y$source_id[!is.na(y$source_id) & nzchar(y$source_id)])),
      stringsAsFactors = FALSE
    )
  }))
}))
write_csv(cohort_counts, file.path(results_dir, "cohort_counts_v4.csv"))

checks <- data.frame(
  check = c(
    "Discovery records",
    "Proteomic dimensions",
    "Discovery/proteomic ID equality",
    "Validation records",
    "Symptomatic validation records",
    "Pilot records",
    "Pilot is a subset of validation rather than an independent cohort"
  ),
  observed = c(
    nrow(discovery),
    paste(nrow(intensity), ncol(intensity), sep = " x "),
    setequal(sample_ids, discovery$study_id),
    nrow(validation),
    nrow(elisa_primary),
    nrow(pilot_raw),
    "Confirmed by source-ID and total-SDC1 matching audit"
  ),
  expected = c("307", "3757 x 307", "TRUE", "415", "329", "96", "Subset"),
  pass = c(
    nrow(discovery) == 307L,
    nrow(intensity) == 3757L && ncol(intensity) == 307L,
    setequal(sample_ids, discovery$study_id),
    nrow(validation) == 415L,
    nrow(elisa_primary) == 329L,
    nrow(pilot_raw) == 96L,
    TRUE
  ),
  stringsAsFactors = FALSE
)
if (!all(checks$pass)) stop("One or more V4 reconciliation checks failed")
write_csv(checks, file.path(results_dir, "v4_reconciliation_checks.csv"))

source_decisions <- data.frame(
  layer = c(
    "Discovery proteomics",
    "Discovery AF4-EV ELISA",
    "Validation total-plasma ELISA",
    "Pilot multicompartment assays",
    "Participant identity",
    "ADD-RS"
  ),
  decision = c(
    "Latest 3,757 x 307 workbook is authoritative; cell values match V3",
    "Latest discovery workbook is authoritative; AAD EV SDC1 was updated",
    "Latest 415-record workbook is authoritative",
    "96-record workbook is a measured subset of the validation cohort",
    "Row UID is the primary record key; source IDs are retained for clustered and sensitivity analyses",
    "Available only in AAD sheets and therefore cannot be used for case-versus-comparator diagnostic benchmarking"
  ),
  prohibited_inference = c(
    "Do not call missing intensities low concentrations without a missingness sensitivity analysis",
    "Do not combine AF4-EV concentrations with DIA-MS intensities",
    "Do not describe the 415 records as external or multicentre validation",
    "Do not count pilot records as an independent validation sample",
    "Do not assume duplicated source IDs are independent participants",
    "Do not claim ADD-RS incremental diagnostic performance from case-only scores"
  ),
  stringsAsFactors = FALSE
)
write_csv(source_decisions, file.path(results_dir, "v4_source_decisions.csv"))

writeLines(capture.output(sessionInfo()), file.path(logs_dir, "00_sessionInfo.txt"))
cat("Prepared revision inputs: discovery n=", nrow(discovery), ", validation n=", nrow(validation),
    ", symptomatic n=", nrow(elisa_primary), ", pilot subset n=", nrow(pilot_raw), "\n", sep = "")

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "02_analysis_data")
qa_dir <- file.path(root, "07_quality_control")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  discovery = "discovery_clinical_deidentified.csv",
  clinical_evaluation = "clinical_evaluation_all.csv",
  elisa_primary = "elisa_evaluation_primary.csv",
  paired_subset = "paired_sdc1_subset.csv"
)

pattern <- paste(
  c("sdc1", "dimer", "crp", "bnp", "ctnt", "ckmb", "wbc", "rbc", "plate",
    "haemo", "creatin", "albumin", "alt", "ast", "glucose", "chol", "hdl",
    "ldl", "trig", "onset", "shock", "death", "failure", "malperfusion", "add_rs"),
  collapse = "|"
)

profile_rows <- list()
schema_rows <- list()
for (cohort_name in names(files)) {
  path <- file.path(data_dir, files[[cohort_name]])
  data <- read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  cat("\n### ", cohort_name, ": n=", nrow(data), ", p=", ncol(data), "\n", sep = "")
  if ("group" %in% names(data)) print(table(data$group, useNA = "ifany"))

  schema_rows[[cohort_name]] <- data.frame(
    cohort = cohort_name,
    variable = names(data),
    class = vapply(data, function(x) class(x)[1], character(1)),
    nonmissing = vapply(data, function(x) sum(!is.na(x) & as.character(x) != ""), integer(1)),
    stringsAsFactors = FALSE
  )

  candidates <- grep(pattern, names(data), ignore.case = TRUE, value = TRUE)
  for (variable in candidates) {
    available <- !is.na(data[[variable]]) & as.character(data[[variable]]) != ""
    overall <- data.frame(
      cohort = cohort_name,
      variable = variable,
      group = "Overall",
      n_total = nrow(data),
      n_available = sum(available),
      available_fraction = mean(available),
      stringsAsFactors = FALSE
    )
    profile_rows[[length(profile_rows) + 1L]] <- overall
    if ("group" %in% names(data)) {
      for (group_name in unique(as.character(data$group))) {
        selected <- as.character(data$group) == group_name
        profile_rows[[length(profile_rows) + 1L]] <- data.frame(
          cohort = cohort_name,
          variable = variable,
          group = group_name,
          n_total = sum(selected),
          n_available = sum(available[selected]),
          available_fraction = mean(available[selected]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

profile <- do.call(rbind, profile_rows)
schema <- do.call(rbind, schema_rows)
write.csv(profile, file.path(qa_dir, "FIGURE7_CLINICAL_VARIABLE_COMPLETENESS_PROFILE.csv"), row.names = FALSE)
write.csv(schema, file.path(qa_dir, "FIGURE7_COHORT_SCHEMA_PROFILE.csv"), row.names = FALSE)

print(profile[profile$group == "Overall", ], row.names = FALSE)

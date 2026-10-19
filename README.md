# SDC1 acute aortic dissection analysis code

This is a **code-only** release for the revision of “Plasma extracellular vesicles proteomics identifies circulating SDC1 as a diagnostic biomarker for acute aortic dissection.” It contains 35 R scripts and one Python preprocessing script. Patient-level data, raw MS files, processed proteomics tables, figures and local software libraries are not included.

The original analysis directory names have been retained for provenance. Packaging changed only machine-specific paths and moved the author-resolved patient-key mapping out of the Python code into a required private JSON configuration. No analysis was rerun or result recalculated while packaging.

## Code map

| Directory | Role | Status |
| --- | --- | --- |
| `SDC1_AAD_REANALYSIS_V7_TRANSPARENT_IMPUTATION_20260716/03_code/R/` | Proteomics filtering, left-censored multiple imputation, adjusted effects and candidate ranking | Current discovery method |
| `SDC1_AAD_HARRELL_BOOTCAL_20260728/03_code/R/03_*` and `04_*` | Nested full-pipeline candidate-rank bootstrap and audit | Current stability method |
| `SDC1_AAD_REANALYSIS_FINAL_20260721/03_code/R/` | Clinical group comparisons, paired compartments, clinical models, associations and figure scripts | Analysis source; some earlier display/performance versions were later superseded |
| `SDC1_AAD_HARRELL_BOOTCAL_20260728/03_code/R/01_*`, `02_*`, `05_*`, `06_*`, `08_*` | Harrell optimism correction, ABCLOC intervals, standalone models and updated panels | Current performance and calibration method |
| `SDC1_AAD_REVISED_SUBMISSION_FINAL_20260722/06_code/` | Supplementary PCA, protein-correlation audit and paired-assay display | Supporting code |
| `SDC1_AAD_PAIRED_AAD_ONLY_CORRELATIONS_20260727/03_code/R/` | AAD-only paired-assay correlations | Supporting code |
| `SDC1_AAD_REANALYSIS_V6_1_REVIEWER_REFINED_20260715/03_code/R/` | Input preparation and selected QC/power/theme analyses | Provenance; old primary protein effects were superseded by V7 |
| `preprocessing/build_final_analysis_data.py` | Builds final derived clinical inputs from the locked workbooks | Requires a private patient-key override file |

The paired-compartment source was updated in the **FINAL_20260727** lock; the other three source workbooks were unchanged from **FINAL_20260721**. Harrell/ABCLOC outputs must not be mixed with earlier repeated-CV/percentile-bootstrap outputs. `CODE_MANIFEST.json` lists SHA-256 hashes for all 36 packaged code files.

## Input and privacy boundary

The scripts expect their local project input folders, such as `02_analysis_data/`, to contain the source or derived data. In this release those folders contain only `.gitkeep`. The clinical preprocessing script also expects the four locked Excel workbooks under `data_lock/FINAL_20260727/` and requires the environment variable `SDC1_PATIENT_KEY_OVERRIDES` to point to a private JSON file with the following shape:

```json
{"overrides": [{"cohort": "clinical_evaluation", "group": "AAD", "source_row": 1, "analysis_key": "PRIVATE_KEY"}]}
```

The example is a format illustration, **not a study record**. The real mapping must not be committed. The data-lock workbooks and participant-level derived CSVs are likewise excluded. This code release documents the algorithms, but **does not by itself reproduce the numerical results**. A separate, approved and suitably de-identified data release with a data dictionary is required for independent execution.

This GitHub repository also does not replace the reviewer-requested public deposition of raw MS data, processed proteomics data and metadata in a specialist repository.

## Software

The original analyses used **R 4.6.0**. Packages used across the scripts include `limma`, `ggplot2`, `patchwork`, `ggrepel`, `glmnet`, `readxl`, `ragg`, `svglite`, `scales`, `clusterProfiler`, `org.Hs.eg.db`, `rms` and `Hmisc`. The Harrell script checks for at least `rms` 8.1.1 and `Hmisc` 5.2.6. The Python preprocessing script requires Python 3, `numpy`, `pandas` and `openpyxl`. The original full dependency lockfile is unavailable; exact software-version reproduction needs verification against run records.

The V7 script uses `SDC1_V7_ROOT`; clinical scripts use `SDC1_FINAL_ROOT`; several Harrell scripts use the repository root as their working directory. Random seeds and resampling defaults are preserved in the source. A full run may take substantial time and should only be attempted after the required data and patient-key configuration are available.

## Suggested reading order

1. Read the V7 script and Harrell `03_run_discovery_nested_mi_selection_bootstrap.R` for selection and stability.
2. Read the final clinical scripts for assay comparisons and clinical associations.
3. Read Harrell `01_run_harrell_bootcal_analysis.R` and `05_run_standalone_harrell_bootcal.R` for the current clinical performance estimation.

## Manuscript wording after publication

Only after the repository is public and its URL and archived version are verified may the manuscript say:

> The analysis code used for proteomic candidate prioritization, statistical analyses and figure generation is available at [GitHub repository URL]; the version used for this revision is archived at [DOI]. Access to participant-level data is described separately in the Data Availability statement.

The brackets are placeholders. This package has no assigned GitHub URL, DOI, license or MS-data accession.

from __future__ import annotations

import hashlib
import json
import math
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


WORKSPACE = Path(__file__).resolve().parents[1]
LOCK = WORKSPACE / "data_lock" / "FINAL_20260727"
PROJECT = WORKSPACE / "SDC1_AAD_REANALYSIS_FINAL_20260721"
ANALYSIS = PROJECT / "02_analysis_data"
AUDIT = PROJECT / "07_quality_control"
RESULTS = PROJECT / "04_results" / "statistical_results"
SOURCE_DATA = PROJECT / "04_results" / "source_data"
LOGS = PROJECT / "04_results" / "logs"

FILES = {
    "discovery": LOCK / "Discovery_Cohort_raw_data_FINAL_20260727.xlsx",
    "paired_subset": LOCK / "Paired_Compartment_Subset_raw_data_FINAL_20260727.xlsx",
    "clinical_evaluation": LOCK / "Clinical_Evaluation_Cohort_raw_data_FINAL_20260727.xlsx",
    "proteomics": LOCK / "Plasma_EV_proteomics_Discovery_Cohort_FINAL_20260727.xlsx",
}

GROUP_MAP = {
    "HC": "Healthy control",
    "AAD": "AAD",
    "AMI": "AMI",
    "Angina": "Angina",
    "PE": "PE",
}

def load_cluster_overrides() -> dict[tuple[str, str, int], str]:
    config = os.environ.get("SDC1_PATIENT_KEY_OVERRIDES")
    if not config:
        raise RuntimeError(
            "Set SDC1_PATIENT_KEY_OVERRIDES to a private JSON file containing "
            "the author-resolved patient-key mapping; do not commit that file."
        )
    records = json.loads(Path(config).read_text(encoding="utf-8"))["overrides"]
    return {
        (entry["cohort"], entry["group"], int(entry["source_row"])): entry["analysis_key"]
        for entry in records
    }


CLUSTER_OVERRIDES = load_cluster_overrides()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def norm_header(value: Any) -> str:
    text = str(value or "").strip().lower().replace("μ", "u").replace("^", "")
    text = re.sub(r"[（）()]", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def clean_id(value: Any) -> str | None:
    if pd.isna(value):
        return None
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(float(value)):
            return None
        return str(int(value)) if float(value).is_integer() else format(float(value), ".15g")
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "n/a", "none", "null", "/", "-", "--"}:
        return None
    if re.fullmatch(r"\d+\.0+", text):
        return text.split(".")[0]
    return re.sub(r"\s+", "", text)


def numeric(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip()
    text = text.replace({"": pd.NA, "/": pd.NA, "-": pd.NA, "--": pd.NA, "NA": pd.NA, "N/A": pd.NA})
    text = text.str.replace(",", "", regex=False)
    text = text.str.replace(r"^[<>]=?\s*", "", regex=True)
    return pd.to_numeric(text, errors="coerce")


def binary(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip().str.lower()
    mapped = text.map(
        {
            "yes": 1.0,
            "y": 1.0,
            "true": 1.0,
            "present": 1.0,
            "female": 1.0,
            "1": 1.0,
            "1.0": 1.0,
            "no": 0.0,
            "n": 0.0,
            "false": 0.0,
            "absent": 0.0,
            "male": 0.0,
            "0": 0.0,
            "0.0": 0.0,
        }
    )
    numeric_fallback = pd.to_numeric(text, errors="coerce")
    return mapped.fillna(numeric_fallback).astype(float)


def select_column(
    frame: pd.DataFrame,
    patterns: list[str],
    *,
    exclude: list[str] | None = None,
) -> pd.Series:
    exclude = exclude or []
    candidates: list[pd.Series] = []
    for position, column in enumerate(frame.columns):
        header = norm_header(column)
        if any(re.search(pattern, header) for pattern in patterns) and not any(
            re.search(pattern, header) for pattern in exclude
        ):
            candidate = frame.iloc[:, position]
            if isinstance(candidate, pd.DataFrame):
                for duplicate_position in range(candidate.shape[1]):
                    candidates.append(candidate.iloc[:, duplicate_position])
            else:
                candidates.append(candidate)
    if not candidates:
        return pd.Series(np.nan, index=frame.index)
    return max(candidates, key=lambda value: numeric(value).notna().sum())


def safe_log2(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float, na_value=np.nan)
    output = np.full(values.shape, np.nan, dtype=float)
    keep = np.isfinite(values) & (values > 0)
    output[keep] = np.log2(values[keep])
    return pd.Series(output, index=series.index)


def canonical_sheet(frame: pd.DataFrame, cohort: str, sheet: str) -> pd.DataFrame:
    group = GROUP_MAP[sheet]
    # Reset the row index before assembling columns so every derived value remains
    # bound to the same source row and patient ID even if blank rows are present.
    frame = frame.dropna(how="all").reset_index(drop=True).copy()
    source_rows = np.arange(2, len(frame) + 2)
    source_id = select_column(frame, [r"^(id|name)$"])
    output = pd.DataFrame(
        {
            "cohort": cohort,
            "source_sheet": sheet,
            "source_row": source_rows,
            "source_id": [clean_id(value) for value in source_id],
            "group": group,
            "aad": int(group == "AAD"),
            "symptomatic": int(group != "Healthy control"),
            "age_years": numeric(select_column(frame, [r"^age$"])),
            "female": binary(select_column(frame, [r"^gender"])),
            "bmi": numeric(select_column(frame, [r"^bmi$"])),
            "onset_hours": numeric(select_column(frame, [r"chest pain time"])),
            "total_sdc1": numeric(select_column(frame, [r"^total sdc1$"])),
            "ev_sdc1": numeric(select_column(frame, [r"^ev sdc1$"])),
            "shed_sdc1": numeric(select_column(frame, [r"^shed sdc1$"])),
            "hypertension": binary(select_column(frame, [r"^hypertension$"])),
            "diabetes": binary(select_column(frame, [r"^diabetes$"])),
            "stroke_history": binary(select_column(frame, [r"^stroke$", r"^non fatal stroke$"])),
            "smoking": binary(select_column(frame, [r"^smoking$"])),
            "hyperlipidaemia": binary(select_column(frame, [r"^hyperlipid"])),
            "wbc_10e9_l": numeric(select_column(frame, [r"^wbc "])),
            "rbc_10e9_l": numeric(select_column(frame, [r"^rbc "])),
            "platelets_10e9_l": numeric(select_column(frame, [r"^plt "])),
            "haemoglobin_g_l": numeric(select_column(frame, [r"^hb "])),
            "alt_u_l": numeric(select_column(frame, [r"^alt "])),
            "ast_u_l": numeric(select_column(frame, [r"^ast "])),
            "glucose_mmol_l": numeric(select_column(frame, [r"^glucose "])),
            "triglycerides_mmol_l": numeric(select_column(frame, [r"^tg "])),
            "total_cholesterol_mmol_l": numeric(select_column(frame, [r"^tc "])),
            "ldl_mmol_l": numeric(select_column(frame, [r"^ldl "])),
            "hdl_mmol_l": numeric(select_column(frame, [r"^hdl "])),
            "albumin_g_l": numeric(select_column(frame, [r"^alb "])),
            "creatinine_umol_l": numeric(select_column(frame, [r"^creatinine "])),
            "bnp_source_units": numeric(select_column(frame, [r"^bnp "])),
            "crp_mg_l": numeric(select_column(frame, [r"^crp "])),
            "ckmb_source_units": numeric(select_column(frame, [r"^ck mb "])),
            "ctnt_source_units": numeric(select_column(frame, [r"^ctnt "])),
            "d_dimer_mg_l": numeric(select_column(frame, [r"^d dimer "])),
            "add_rs_source": numeric(select_column(frame, [r"^add rs score$"])),
            "cardiogenic_shock": binary(select_column(frame, [r"^cardiac shock"])),
            "heart_failure": binary(select_column(frame, [r"^heart failure"])),
            "cardiac_death": binary(select_column(frame, [r"^cardiac death"])),
            "malperfusion": binary(select_column(frame, [r"^malperfusion"])),
        }
    )
    output["row_uid"] = (
        output["cohort"]
        + "_"
        + output["source_sheet"]
        + "_R"
        + output["source_row"].astype(str)
    )
    cluster_ids = []
    rules = []
    for row in output.itertuples(index=False):
        override = CLUSTER_OVERRIDES.get((cohort, sheet, int(row.source_row)))
        if override:
            cluster_ids.append(override)
            rules.append("author_resolved_row_override")
        elif row.source_id:
            cluster_ids.append(f"{sheet.upper()}_SRC_{row.source_id}")
            rules.append("group_qualified_nonmissing_source_id")
        else:
            cluster_ids.append(f"{sheet.upper()}_ROW{int(row.source_row)}_MISSING_ID")
            rules.append("stable_row_level_missing_id")
    output["cluster_id"] = cluster_ids
    output["cluster_rule"] = rules
    output["source_id_missing"] = output["source_id"].isna().astype(int)
    output["add_rs_corrected"] = output["add_rs_source"].replace({4.0: 3.0})
    output["add_rs_correction_flag"] = output["add_rs_source"].eq(4.0).fillna(False).astype(int)
    output["aad_stanford_type"] = np.where(output["group"] == "AAD", "A", "")
    output["ami_subtype"] = np.where(output["group"] == "AMI", "STEMI", "")
    output["ami_coronary_angiography_confirmed"] = np.where(
        output["group"] == "AMI", 1, np.nan
    )
    output["log2_total_sdc1"] = safe_log2(output["total_sdc1"])
    output["sdc1_source_value"] = output["total_sdc1"]
    output["log2_ev_sdc1"] = safe_log2(output["ev_sdc1"])
    output["log2_d_dimer"] = safe_log2(output["d_dimer_mg_l"])
    return output


def read_cohort(path: Path, cohort: str) -> pd.DataFrame:
    sheets = pd.read_excel(path, sheet_name=None, engine="openpyxl")
    frames = [canonical_sheet(sheets[name], cohort, name) for name in GROUP_MAP if name in sheets]
    output = pd.concat(frames, ignore_index=True)
    output["source_id_duplicate_raw"] = (
        output["source_id"].notna() & output.duplicated(["group", "source_id"], keep=False)
    ).astype(int)
    return output


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, encoding="utf-8-sig", na_rep="")


def inclusion_audit(discovery: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    rows = []
    specifications = {
        "fig3_ev_sdc1_group_adjusted": (
            discovery,
            ["ev_sdc1", "group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l"],
        ),
        "fig3_d_dimer_group_adjusted": (
            discovery,
            ["d_dimer_mg_l", "group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l"],
        ),
        "fig3_ctnt_group_adjusted": (
            discovery,
            ["ctnt_source_units", "group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l"],
        ),
        "fig4_total_sdc1_group_adjusted": (
            clinical,
            ["total_sdc1", "group", "age_years", "female", "hypertension", "diabetes", "smoking", "creatinine_umol_l", "wbc_10e9_l"],
        ),
        "fig4_d_dimer_group_adjusted_symptomatic": (
            clinical.loc[clinical["symptomatic"] == 1],
            ["d_dimer_mg_l", "group", "age_years", "female", "hypertension", "diabetes", "smoking", "creatinine_umol_l", "crp_mg_l"],
        ),
        "fig4_head_to_head_all_symptomatic": (
            clinical.loc[clinical["symptomatic"] == 1],
            ["aad", "total_sdc1", "d_dimer_mg_l", "cluster_id"],
        ),
        "fig4_head_to_head_aad_vs_pe": (
            clinical.loc[clinical["group"].isin(["AAD", "PE"])],
            ["aad", "total_sdc1", "d_dimer_mg_l", "cluster_id"],
        ),
    }
    for analysis, (frame, variables) in specifications.items():
        required = frame[variables].copy()
        complete = required.notna().all(axis=1)
        for group, subset in frame.assign(complete=complete).groupby("group", sort=False):
            rows.append(
                {
                    "analysis": analysis,
                    "group": group,
                    "available_records": len(subset),
                    "included_records": int(subset["complete"].sum()),
                    "excluded_for_missing_required_variable": int((~subset["complete"]).sum()),
                    "required_variables": "|".join(variables),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    for directory in [ANALYSIS, AUDIT, RESULTS, SOURCE_DATA, LOGS, PROJECT / "03_code" / "R"]:
        directory.mkdir(parents=True, exist_ok=True)

    discovery = read_cohort(FILES["discovery"], "discovery")
    paired = read_cohort(FILES["paired_subset"], "paired_subset")
    clinical = read_cohort(FILES["clinical_evaluation"], "clinical_evaluation")

    write_csv(discovery, ANALYSIS / "discovery_clinical_deidentified.csv")
    write_csv(paired, ANALYSIS / "paired_sdc1_subset.csv")
    write_csv(clinical, ANALYSIS / "clinical_evaluation_all.csv")

    symptomatic = clinical.loc[clinical["symptomatic"] == 1].copy()
    symptomatic["symptomatic_comparator"] = (symptomatic["group"] != "AAD").astype(int)
    write_csv(symptomatic, ANALYSIS / "elisa_evaluation_primary.csv")

    inclusion = inclusion_audit(discovery, clinical)
    write_csv(inclusion, AUDIT / "analysis_inclusion_audit.csv")
    write_csv(
        clinical.loc[
            clinical["cluster_rule"] != "group_qualified_nonmissing_source_id",
            ["row_uid", "group", "source_row", "source_id", "cluster_id", "cluster_rule"],
        ],
        AUDIT / "author_resolved_cluster_keys.csv",
    )
    cross_group_collisions = clinical.loc[
        clinical["source_id"].notna()
        & clinical.duplicated("source_id", keep=False)
        & ~clinical.duplicated(["group", "source_id"], keep=False),
        ["row_uid", "group", "source_row", "source_id", "cluster_id", "cluster_rule"],
    ].sort_values(["source_id", "group"])
    write_csv(cross_group_collisions, AUDIT / "cross_group_raw_id_collisions.csv")
    write_csv(
        pd.concat(
            [
                discovery.loc[discovery["add_rs_correction_flag"] == 1],
                clinical.loc[clinical["add_rs_correction_flag"] == 1],
            ]
        )[["cohort", "row_uid", "group", "source_id", "add_rs_source", "add_rs_corrected"]],
        AUDIT / "add_rs_correction_audit.csv",
    )

    proteomics_header = pd.read_excel(FILES["proteomics"], sheet_name=0, nrows=0, engine="openpyxl").columns.tolist()
    proteomics_ids = {clean_id(value) for value in proteomics_header[1:]}
    discovery_ids = set(discovery["source_id"].dropna())

    manifest = {
        "built_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "authoritative_lock": str(LOCK),
        "lock_version": "FINAL_20260727",
        "locked_input_hashes": {name: sha256(path) for name, path in FILES.items()},
        "record_counts": {
            "discovery": discovery.groupby("group").size().to_dict(),
            "paired_subset": paired.groupby("group").size().to_dict(),
            "clinical_evaluation": clinical.groupby("group").size().to_dict(),
            "symptomatic_modelling_set": symptomatic.groupby("group").size().to_dict(),
        },
        "clinical_analysis_cluster_count": int(clinical["cluster_id"].nunique()),
        "clinical_row_count": int(len(clinical)),
        "clinical_cluster_ids_unique": bool(clinical["cluster_id"].is_unique),
        "author_override_rows": int((clinical["cluster_rule"] == "author_resolved_row_override").sum()),
        "missing_id_row_clusters": int(clinical["source_id_missing"].sum()),
        "cross_group_raw_id_collision_rows": int(len(cross_group_collisions)),
        "add_rs_corrected_rows": int(
            discovery["add_rs_correction_flag"].sum() + clinical["add_rs_correction_flag"].sum()
        ),
        "proteomics_dimensions": [3757, len(proteomics_header) - 1],
        "discovery_proteomics_id_match": discovery_ids == proteomics_ids,
        "available_case_rule": "Each analysis includes all records with nonmissing values for that analysis's required variables; there is no global complete-case deletion.",
    }
    (AUDIT / "derived_data_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if len(discovery) != 307 or len(paired) != 96 or len(clinical) != 415 or len(symptomatic) != 329:
        raise RuntimeError("Unexpected cohort record count")
    if not clinical["cluster_id"].is_unique or clinical["cluster_id"].nunique() != 415:
        raise RuntimeError("Clinical analysis cluster IDs must be unique after author resolution")
    if discovery_ids != proteomics_ids:
        raise RuntimeError("Discovery clinical and proteomics sample IDs do not match")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

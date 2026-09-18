#!/usr/bin/env python

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CASE_ID = "example_cup_case"
OUTPUT = ROOT / "results" / "clinical_projection" / CASE_ID
REFERENCE = ROOT / "results" / "reference_atlas"


def fail(message: str) -> None:
    raise AssertionError(message)


paths = {
    "feature_vector": OUTPUT / "clinical_feature_vector.tsv",
    "scaled_vector": OUTPUT / "clinical_feature_vector_scaled.tsv",
    "missingness": OUTPUT / "clinical_feature_missingness.tsv",
    "neighbors": OUTPUT / "tcga_nearest_neighbors.tsv",
    "project_summary": OUTPUT / "tcga_nearest_project_summary.tsv",
    "group_summary": OUTPUT / "tcga_nearest_major_group_summary.tsv",
    "pca": OUTPUT / "clinical_tcga_pca_projection.tsv",
    "qc": OUTPUT / "clinical_projection_qc_summary.tsv",
    "pca_project_figure": OUTPUT / "clinical_tcga_pca_projection_by_project.pdf",
    "pca_group_figure": OUTPUT / "clinical_tcga_pca_projection_by_major_group.pdf",
    "neighbors_figure": OUTPUT / "clinical_nearest_neighbor_barplot.pdf",
    "missingness_figure": OUTPUT / "clinical_feature_missingness.pdf",
    "report": OUTPUT / "clinical_tcga_projection_report.md",
}

missing_outputs = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
if missing_outputs:
    fail(f"Missing or empty clinical projection outputs: {', '.join(missing_outputs)}")

feature_vector = pd.read_csv(paths["feature_vector"], sep="\t")
scaled_vector = pd.read_csv(paths["scaled_vector"], sep="\t")
missingness = pd.read_csv(paths["missingness"], sep="\t")
neighbors = pd.read_csv(paths["neighbors"], sep="\t")
project_summary = pd.read_csv(paths["project_summary"], sep="\t")
group_summary = pd.read_csv(paths["group_summary"], sep="\t")
pca = pd.read_csv(paths["pca"], sep="\t")
qc = pd.read_csv(paths["qc"], sep="\t")
reference_metadata = pd.read_csv(REFERENCE / "tcga_sample_reference_metadata.tsv", sep="\t")
with (REFERENCE / "tcga_reference_scaler_parameters.json").open() as handle:
    scaler = json.load(handle)
feature_order = scaler["feature_order"]

if len(feature_vector) != 1 or len(scaled_vector) != 1:
    fail("Example clinical feature vectors must contain exactly one query row")
if list(feature_vector.columns) != ["case_id", "sample_id"] + feature_order:
    fail("Harmonized clinical vector does not preserve locked feature order")
if list(scaled_vector.columns) != ["case_id", "sample_id"] + feature_order:
    fail("Scaled clinical vector does not preserve locked feature order")
if feature_vector["case_id"].iloc[0] != CASE_ID or feature_vector["sample_id"].iloc[0] != "SYNTHETIC-CUP-001":
    fail("Example clinical vector identifiers are incorrect")
if not np.isfinite(scaled_vector[feature_order].to_numpy(dtype=float)).all():
    fail("Scaled clinical vector contains non-finite values")

reconstructed = feature_vector[feature_order].copy()
for feature in feature_order:
    reconstructed[feature] = (
        reconstructed[feature] - scaler["scaling_mean"][feature]
    ) / scaler["scaling_std"][feature]
if not np.allclose(
    reconstructed.to_numpy(dtype=float),
    scaled_vector[feature_order].to_numpy(dtype=float),
    atol=1e-10,
    rtol=1e-10,
):
    fail("Example clinical scaling does not use saved TCGA reference parameters")

required_missingness_columns = {
    "case_id",
    "sample_id",
    "feature_name",
    "expected_in_reference",
    "feature_class",
    "clinical_modality",
    "source_file",
    "supplied",
    "imputed",
    "imputation_value",
    "value_used",
    "scaled_value",
    "status",
    "notes",
}
if not required_missingness_columns.issubset(missingness.columns):
    fail("Clinical feature missingness table lacks required columns")
expected_rows = missingness[missingness["expected_in_reference"]].copy()
extra_rows = missingness[~missingness["expected_in_reference"]].copy()
if expected_rows["feature_name"].tolist() != feature_order:
    fail("Clinical feature missingness table does not preserve locked feature order")
if len(expected_rows) != 97 or int(expected_rows["supplied"].sum()) != 77 or int(expected_rows["imputed"].sum()) != 20:
    fail("Example clinical feature coverage counts are incorrect")
missing_features = expected_rows.loc[expected_rows["imputed"], "feature_name"].tolist()
if missing_features != [f"expression_pca_PC{i}" for i in range(1, 21)]:
    fail("Example query should impute exactly the 20 expression PCA features")
for feature in missing_features:
    if not np.isclose(feature_vector[feature].iloc[0], scaler["imputation_median"][feature]):
        fail(f"Missing clinical feature was not imputed with the TCGA reference median: {feature}")
if len(extra_rows) != 4 or set(extra_rows["status"]) != {"extra_not_in_reference_ignored"}:
    fail("Optional non-reference clinical features were not preserved as ignored audit rows")

required_neighbor_columns = {
    "case_id",
    "query_sample_id",
    "rank",
    "tcga_sample_barcode",
    "tcga_project_id",
    "tcga_project_code",
    "tcga_disease_type",
    "tcga_primary_site",
    "tcga_major_cancer_group",
    "distance",
    "similarity_score",
    "distance_metric",
    "notes",
}
if set(neighbors.columns) != required_neighbor_columns:
    fail("Nearest-neighbor table columns do not match the required contract")
if len(neighbors) != 50 or neighbors["rank"].tolist() != list(range(1, 51)):
    fail("Nearest-neighbor table does not contain ranks 1 through 50")
if not neighbors["distance"].is_monotonic_increasing:
    fail("Nearest-neighbor distances are not sorted")
if not np.allclose(neighbors["similarity_score"], 1.0 / (1.0 + neighbors["distance"])):
    fail("Nearest-neighbor similarity scores do not match the documented transformation")
if set(neighbors["distance_metric"]) != {"euclidean"}:
    fail("Nearest-neighbor rows do not record the requested Euclidean metric")
if not set(neighbors["tcga_sample_barcode"]).issubset(set(reference_metadata["sample_barcode"])):
    fail("Nearest-neighbor output contains unknown TCGA sample barcodes")
valid_project_pairs = set(
    map(tuple, reference_metadata[["project_id", "project_code"]].drop_duplicates().to_numpy())
)
neighbor_project_pairs = set(
    map(tuple, neighbors[["tcga_project_id", "tcga_project_code"]].drop_duplicates().to_numpy())
)
if not neighbor_project_pairs.issubset(valid_project_pairs):
    fail("Nearest-neighbor output contains invalid TCGA project labels")

required_project_columns = {
    "project_code",
    "project_id",
    "n_neighbors",
    "median_distance",
    "min_distance",
    "weighted_similarity_score",
    "rank",
    "interpretation_note",
}
if set(project_summary.columns) != required_project_columns or project_summary.empty:
    fail("Project-level nearest-neighbor summary is missing or malformed")
if int(project_summary["n_neighbors"].sum()) != 50:
    fail("Project-level summary does not account for all nearest neighbors")
if not np.isclose(project_summary["weighted_similarity_score"].sum(), 1.0):
    fail("Project-level similarity weights do not sum to one")
if project_summary["rank"].tolist() != list(range(1, len(project_summary) + 1)):
    fail("Project-level summary ranks are malformed")

required_group_columns = {
    "major_cancer_group",
    "n_neighbors",
    "median_distance",
    "min_distance",
    "weighted_similarity_score",
    "rank",
    "interpretation_note",
}
if set(group_summary.columns) != required_group_columns or group_summary.empty:
    fail("Major-group nearest-neighbor summary is missing or malformed")
if int(group_summary["n_neighbors"].sum()) != 50:
    fail("Major-group summary does not account for all nearest neighbors")
if not np.isclose(group_summary["weighted_similarity_score"].sum(), 1.0):
    fail("Major-group similarity weights do not sum to one")

pc_columns = [column for column in pca.columns if re.fullmatch(r"PC\d+", column)]
if len(pc_columns) != 20 or not np.isfinite(pca[pc_columns].to_numpy(dtype=float)).all():
    fail("Clinical PCA projection does not contain 20 finite saved-model coordinates")
if "PCA was not refit" not in pca["projection_note"].iloc[0]:
    fail("Clinical PCA output does not document saved-model projection")

required_qc_metrics = {
    "number_of_expected_reference_features",
    "atlas_feature_count_used",
    "number_of_clinical_features_supplied",
    "number_of_missing_clinical_features",
    "percent_feature_coverage",
    "WES_features_present",
    "WTS_features_present",
    "projection_status",
    "warnings",
    "limitations",
}
if not required_qc_metrics.issubset(set(qc["metric"])):
    fail("Clinical projection QC summary lacks required metrics")
global_qc = qc[qc["qc_section"] == "global"]
qc_values = dict(zip(global_qc["metric"], global_qc["value"].astype(str), strict=False))
if qc_values["number_of_expected_reference_features"] != "97":
    fail("Clinical QC expected-feature count is incorrect")
if qc_values["atlas_feature_count_used"] != "97":
    fail("Clinical QC atlas-feature count is incorrect")
if qc_values["number_of_clinical_features_supplied"] != "77":
    fail("Clinical QC supplied-feature count is incorrect")
if qc_values["number_of_missing_clinical_features"] != "20":
    fail("Clinical QC missing-feature count is incorrect")
if qc_values["WES_features_present"] != "True" or qc_values["WTS_features_present"] != "True":
    fail("Clinical QC does not report both WES and WTS inputs")
if not qc_values["projection_status"].startswith("success"):
    fail("Clinical QC does not report successful projection")

report = paths["report"].read_text()
for phrase in [
    "research/prototype molecular projection",
    "not a validated clinical diagnostic classifier",
    "does not establish tissue of origin",
    "No supervised cancer-of-unknown-primary classifier",
    "does not process FASTQ, BAM, VCF, or raw RNA-seq files",
    "Locked TCGA atlas features used: 97",
    "## PCA Projection Summary",
]:
    if phrase not in report:
        fail(f"Clinical projection report lacks required caveat: {phrase}")
if "Artificial demonstration data only" not in report:
    fail("Clinical projection report does not preserve the synthetic-data label")

for key in ["pca_project_figure", "pca_group_figure", "neighbors_figure", "missingness_figure"]:
    with paths[key].open("rb") as handle:
        if handle.read(4) != b"%PDF":
            fail(f"Clinical projection figure is not a valid PDF: {paths[key]}")

print("Clinical projection output checks passed")

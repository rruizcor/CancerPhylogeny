#!/usr/bin/env python

from __future__ import annotations

import json
import os
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


paths = {
    "reference_matrix": ROOT / "results" / "reference_atlas" / "tcga_sample_reference_matrix.tsv.gz",
    "metadata": ROOT / "results" / "reference_atlas" / "tcga_sample_reference_metadata.tsv",
    "scaled_matrix": ROOT
    / "results"
    / "reference_atlas"
    / "tcga_sample_reference_scaled_matrix.tsv.gz",
    "feature_definitions": ROOT
    / "results"
    / "reference_atlas"
    / "tcga_reference_feature_definitions.yaml",
    "scaler": ROOT / "results" / "reference_atlas" / "tcga_reference_scaler_parameters.json",
    "feature_missingness": ROOT
    / "results"
    / "reference_atlas"
    / "tcga_reference_feature_missingness.tsv",
    "pca_coordinates": ROOT
    / "results"
    / "reference_atlas"
    / "tcga_reference_pca_coordinates.tsv",
    "pca_model": ROOT / "results" / "reference_atlas" / "tcga_reference_pca_model.pkl",
    "pca_variance": ROOT / "results" / "reference_atlas" / "tcga_reference_pca_variance.tsv",
    "nearest_neighbor_index": ROOT
    / "results"
    / "reference_atlas"
    / "tcga_reference_nearest_neighbor_index.pkl",
    "nearest_neighbor_metadata": ROOT
    / "results"
    / "reference_atlas"
    / "tcga_reference_nearest_neighbor_metadata.json",
    "qc": ROOT / "results" / "tables" / "tcga_sample_reference_atlas_qc_summary.tsv",
    "pca_project_figure": ROOT
    / "results"
    / "figures"
    / "tcga_sample_reference_pca_by_project.pdf",
    "pca_group_figure": ROOT
    / "results"
    / "figures"
    / "tcga_sample_reference_pca_by_major_group.pdf",
    "missingness_figure": ROOT
    / "results"
    / "figures"
    / "tcga_sample_reference_feature_missingness.pdf",
    "samples_by_project_figure": ROOT
    / "results"
    / "figures"
    / "tcga_sample_reference_samples_by_project.pdf",
}

missing_or_empty = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
if missing_or_empty:
    fail(f"Missing or empty TCGA reference atlas outputs: {', '.join(missing_or_empty)}")

reference = pd.read_csv(paths["reference_matrix"], sep="\t")
metadata = pd.read_csv(paths["metadata"], sep="\t")
scaled = pd.read_csv(paths["scaled_matrix"], sep="\t")
feature_missingness = pd.read_csv(paths["feature_missingness"], sep="\t")
pca_coordinates = pd.read_csv(paths["pca_coordinates"], sep="\t")
pca_variance = pd.read_csv(paths["pca_variance"], sep="\t")
qc = pd.read_csv(paths["qc"], sep="\t")
projects = pd.read_csv(ROOT / "data" / "interim" / "tcga_projects.tsv", sep="\t")
driver_genes = pd.read_csv(ROOT / "config" / "driver_genes.csv")

required_metadata_columns = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "disease_type",
    "primary_site",
    "major_cancer_group",
]
if not set(required_metadata_columns).issubset(metadata.columns):
    fail("Reference metadata lacks required columns")
if not set(required_metadata_columns).issubset(reference.columns):
    fail("Raw reference matrix lacks required metadata columns")
if metadata[required_metadata_columns].isna().any().any():
    fail("Required metadata columns contain missing values")
if metadata["sample_barcode"].duplicated().any():
    fail("Reference metadata contain duplicate sample barcodes")
if not metadata["sample_barcode"].equals(reference["sample_barcode"]):
    fail("Raw reference matrix and metadata use different sample ordering")
if not metadata["sample_barcode"].equals(scaled["sample_barcode"]):
    fail("Scaled reference matrix and metadata use different sample ordering")

if not (len(reference) == len(metadata) == len(scaled) == len(pca_coordinates)):
    fail("Reference matrix, metadata, scaled matrix, and PCA coordinates have inconsistent row counts")
if len(metadata) != 10201:
    fail(f"Expected 10,201 mutation-layer tumor samples in atlas, found {len(metadata)}")
if metadata["project_code"].nunique() != 33:
    fail("Reference atlas does not retain all 33 TCGA projects")

valid_project_pairs = set(map(tuple, projects[["project_id", "project_code"]].drop_duplicates().to_numpy()))
atlas_project_pairs = set(map(tuple, metadata[["project_id", "project_code"]].drop_duplicates().to_numpy()))
if not atlas_project_pairs.issubset(valid_project_pairs):
    fail("Reference atlas contains invalid TCGA project labels")
if len(atlas_project_pairs) != 33:
    fail("Reference atlas does not contain 33 valid project label pairs")

barcode_pattern = re.compile(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2})[A-Z0-9]$")
sample_type_codes = []
for barcode in metadata["sample_barcode"].astype(str):
    match = barcode_pattern.match(barcode)
    if not match:
        fail(f"Malformed TCGA sample barcode detected: {barcode}")
    sample_type_codes.append(int(match.group(1)))
if not all(1 <= code <= 9 for code in sample_type_codes):
    fail("Reference metadata contain non-tumor TCGA sample types")

required_feature_columns = {
    "nonsynonymous_count",
    "total_mutation_count",
    "purity",
    "ploidy",
    "aneuploidy_score",
    "arm_gain_count",
    "arm_loss_count",
    "total_arm_alteration_count",
    "proliferation_score",
    "immune_inflammatory_score",
}
if not required_feature_columns.issubset(reference.columns):
    fail("Raw reference matrix lacks expected molecular feature columns")

driver_binary_columns = {f"driver_mutation_{gene}" for gene in driver_genes["gene"].astype(str)}
driver_count_columns = {f"driver_mutation_count_{gene}" for gene in driver_genes["gene"].astype(str)}
if not driver_binary_columns.issubset(reference.columns):
    fail("Raw reference matrix lacks one or more configured driver mutation indicators")
if not driver_count_columns.issubset(reference.columns):
    fail("Raw reference matrix lacks one or more configured driver mutation counts")

with paths["scaler"].open() as handle:
    scaler = json.load(handle)
required_scaler_keys = {
    "atlas_schema_version",
    "raw_feature_order",
    "feature_order",
    "imputation_median",
    "imputation_values",
    "scaling_mean",
    "scaling_std",
    "dropped_features",
    "missing_fraction",
    "feature_contract_sha256",
    "sample_order_sha256",
    "serialization_runtime",
    "future_clinical_transform",
    "default_nearest_neighbor_metric",
    "available_nearest_neighbor_metrics",
}
missing_scaler_keys = required_scaler_keys - set(scaler)
if missing_scaler_keys:
    fail(f"Scaler JSON lacks required keys: {sorted(missing_scaler_keys)}")
if not re.fullmatch(r"[0-9a-f]{64}", scaler["feature_contract_sha256"]):
    fail("Scaler feature contract hash is malformed")
if not scaler["feature_order"]:
    fail("Scaler feature order is empty")
if list(reference.columns) != required_metadata_columns + scaler["raw_feature_order"]:
    fail("Raw reference matrix columns do not match saved raw feature order")
if list(scaled.columns) != ["sample_barcode"] + scaler["feature_order"]:
    fail("Scaled matrix columns do not match saved retained feature order")
if scaled[scaler["feature_order"]].isna().any().any():
    fail("Scaled matrix contains NA values")
if not np.isfinite(scaled[scaler["feature_order"]].to_numpy(dtype=float)).all():
    fail("Scaled matrix contains non-finite values")
if not all(pd.api.types.is_numeric_dtype(scaled[column]) for column in scaler["feature_order"]):
    fail("Scaled matrix contains non-numeric feature values")

reconstructed = reference[scaler["feature_order"]].copy()
for feature in scaler["feature_order"]:
    reconstructed[feature] = (
        reconstructed[feature].fillna(scaler["imputation_median"][feature])
        - scaler["scaling_mean"][feature]
    ) / scaler["scaling_std"][feature]
if not np.allclose(
    reconstructed.to_numpy(dtype=float),
    scaled[scaler["feature_order"]].to_numpy(dtype=float),
    atol=1e-9,
    rtol=1e-9,
):
    fail("Saved transform parameters do not reconstruct the scaled reference matrix")

required_missingness_columns = {
    "feature_name",
    "feature_type",
    "source_layer",
    "n_samples",
    "n_observed",
    "n_missing",
    "missing_fraction",
    "missingness_threshold",
    "retained",
    "drop_reason",
    "imputation_median",
}
if not required_missingness_columns.issubset(feature_missingness.columns):
    fail("Feature missingness table lacks required columns")
if feature_missingness["feature_name"].tolist() != scaler["raw_feature_order"]:
    fail("Feature missingness table does not preserve raw feature order")
retained_from_missingness = feature_missingness.loc[feature_missingness["retained"], "feature_name"].tolist()
if retained_from_missingness != scaler["feature_order"]:
    fail("Feature missingness retained flags do not match locked feature order")
if int(feature_missingness["n_missing"].sum()) < scaler["n_imputed_values"]:
    fail("Raw feature missingness is inconsistent with the imputation count")

with paths["feature_definitions"].open() as handle:
    definitions = yaml.safe_load(handle)
if definitions.get("feature_order") != scaler["feature_order"]:
    fail("Feature definitions do not preserve scaler feature order")
if definitions.get("raw_feature_order") != scaler["raw_feature_order"]:
    fail("Feature definitions do not preserve raw feature order")
if definitions.get("feature_contract_sha256") != scaler["feature_contract_sha256"]:
    fail("Feature definitions and scaler use different feature contracts")
feature_records = definitions.get("features", [])
if len(feature_records) != len(scaler["raw_feature_order"]):
    fail("Feature definitions do not contain one record per raw feature")
required_definition_fields = {
    "name",
    "feature_type",
    "source_layer",
    "description",
    "clinical_modality",
    "computable_from_clinical_wes",
    "computable_from_clinical_wts",
    "included_in_nearest_neighbor_index",
}
for record in feature_records:
    if not required_definition_fields.issubset(record):
        fail(f"Feature definition lacks required fields: {record.get('name', '<unknown>')}")
    if record["clinical_modality"] not in {"WES", "WTS", "both"}:
        fail(f"Invalid clinical modality in feature definition: {record['name']}")
if "future_case_rule" not in definitions.get("preprocessing", {}):
    fail("Feature definitions lack the future clinical transform rule")

if not {"PC1", "PC2", "sample_barcode", "project_code", "major_cancer_group"}.issubset(
    pca_coordinates.columns
):
    fail("PCA coordinate table lacks required columns")
if pca_coordinates[["PC1", "PC2"]].isna().any().any():
    fail("PCA coordinates contain NA values")
required_variance_columns = {
    "component",
    "explained_variance",
    "explained_variance_ratio",
    "cumulative_explained_variance_ratio",
    "singular_value",
}
if not required_variance_columns.issubset(pca_variance.columns):
    fail("PCA variance table lacks required columns")
if not pca_variance["cumulative_explained_variance_ratio"].is_monotonic_increasing:
    fail("PCA cumulative explained variance is not monotonic")
if pca_variance["cumulative_explained_variance_ratio"].iloc[-1] > 1.0000001:
    fail("PCA cumulative explained variance exceeds one")

with paths["pca_model"].open("rb") as handle:
    pca_model = pickle.load(handle)
if pca_model.get("feature_order") != scaler["feature_order"]:
    fail("PCA model artifact does not preserve feature order")
if pca_model.get("feature_contract_sha256") != scaler["feature_contract_sha256"]:
    fail("PCA model artifact does not preserve the feature contract")
if "model" not in pca_model:
    fail("PCA model artifact lacks the fitted model object")

with paths["nearest_neighbor_metadata"].open() as handle:
    neighbor_metadata = json.load(handle)
if neighbor_metadata.get("feature_order") != scaler["feature_order"]:
    fail("Nearest-neighbor metadata does not preserve feature order")
if neighbor_metadata.get("feature_contract_sha256") != scaler["feature_contract_sha256"]:
    fail("Nearest-neighbor metadata does not preserve the feature contract")
if set(neighbor_metadata.get("available_metrics", [])) != {"euclidean", "cosine"}:
    fail("Nearest-neighbor metadata does not include Euclidean and cosine options")
if not neighbor_metadata.get("query_contract"):
    fail("Nearest-neighbor metadata lacks future clinical transform/query instructions")
if not neighbor_metadata.get("serialization_runtime", {}).get("scikit_learn"):
    fail("Nearest-neighbor metadata lacks model-serialization runtime versions")

with paths["nearest_neighbor_index"].open("rb") as handle:
    neighbor_index = pickle.load(handle)
if neighbor_index.get("feature_order") != scaler["feature_order"]:
    fail("Nearest-neighbor artifact does not preserve feature order")
if set(neighbor_index.get("indices", {})) != {"euclidean", "cosine"}:
    fail("Nearest-neighbor artifact does not contain both requested metric models")
default_metric = neighbor_metadata["default_metric"]
if neighbor_index.get("default_metric") != default_metric:
    fail("Nearest-neighbor artifact and metadata disagree on the default metric")
distances, indices = neighbor_index["indices"][default_metric].kneighbors(
    scaled[scaler["feature_order"]].iloc[[0]].to_numpy(dtype=float), n_neighbors=5
)
if distances.shape != (1, 5) or indices.shape != (1, 5):
    fail("Nearest-neighbor index query returned unexpected shape")
if indices[0, 0] != 0 or not np.isclose(distances[0, 0], 0.0, atol=1e-6):
    fail("Nearest-neighbor self-query does not return the source row within text round-trip precision")

required_qc_metrics = {
    "number_of_samples_included",
    "number_of_patients_included",
    "number_of_projects_represented",
    "samples_per_project",
    "number_of_raw_features",
    "number_of_retained_features",
    "number_of_dropped_features",
    "dropped_features",
    "number_of_imputed_values",
    "expression_included",
    "pca_built",
    "nearest_neighbor_index_built",
    "limitations",
}
if not required_qc_metrics.issubset(set(qc["metric"])):
    fail("Atlas QC summary lacks required metrics")
metric_values = dict(zip(qc["metric"], qc["value"].astype(str), strict=False))
if metric_values["number_of_samples_included"] != "10201":
    fail("QC summary sample count does not match the reference matrix")
if metric_values["number_of_projects_represented"] != "33":
    fail("QC summary project count does not match the reference metadata")
if metric_values["expression_included"] != "True":
    fail("QC summary does not report expression inclusion")
if metric_values["nearest_neighbor_index_built"] != "True" or metric_values["pca_built"] != "True":
    fail("QC summary does not report PCA and nearest-neighbor construction")
project_qc_rows = qc[qc["qc_section"].astype(str).str.startswith("project:")]
if len(project_qc_rows) != 33:
    fail("QC summary does not include samples-per-project rows for all 33 projects")

for key in [
    "pca_project_figure",
    "pca_group_figure",
    "missingness_figure",
    "samples_by_project_figure",
]:
    with paths[key].open("rb") as handle:
        if handle.read(4) != b"%PDF":
            fail(f"Atlas figure is not a valid PDF: {paths[key]}")

print("TCGA sample reference atlas output checks passed")

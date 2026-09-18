#!/usr/bin/env python

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import pickle
import re
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcga_cancer_phylogeny_matplotlib"))
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.lines import Line2D

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


CONFIG = {
    "top_k": 50,
    "metric": None,
    "force_direct_distance": False,
    "max_abs_scaled_warning": 10.0,
}

CLINICAL_METADATA_REQUIRED_COLUMNS = [
    "case_id",
    "sample_id",
    "specimen_type",
    "submitted_diagnosis",
    "notes",
]

REFERENCE_METADATA_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "disease_type",
    "primary_site",
    "major_cancer_group",
]

NEIGHBOR_COLUMNS = [
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
]


@dataclass
class ProjectionPaths:
    reference_matrix: Path
    reference_metadata: Path
    reference_scaled_matrix: Path
    feature_definitions: Path
    scaler_parameters: Path
    pca_model: Path
    reference_pca_coordinates: Path
    nearest_neighbor_index: Path
    nearest_neighbor_metadata: Path
    clinical_metadata: Path
    wes_features: Path
    wts_features: Path
    output_dir: Path
    clinical_feature_vector: Path
    clinical_feature_vector_scaled: Path
    clinical_feature_missingness: Path
    nearest_neighbors: Path
    nearest_project_summary: Path
    nearest_major_group_summary: Path
    pca_projection: Path
    qc_summary: Path
    pca_by_project_figure: Path
    pca_by_major_group_figure: Path
    nearest_neighbor_figure: Path
    feature_missingness_figure: Path
    report: Path


def default_paths(
    root: Path,
    case_id: str,
    query_dir: Path | None = None,
    output_dir: Path | None = None,
) -> ProjectionPaths:
    query = query_dir or project_path("data", "clinical_queries", case_id, root=root)
    output = output_dir or project_path("results", "clinical_projection", case_id, root=root)
    return ProjectionPaths(
        reference_matrix=project_path(
            "results", "reference_atlas", "tcga_sample_reference_matrix.tsv.gz", root=root
        ),
        reference_metadata=project_path(
            "results", "reference_atlas", "tcga_sample_reference_metadata.tsv", root=root
        ),
        reference_scaled_matrix=project_path(
            "results", "reference_atlas", "tcga_sample_reference_scaled_matrix.tsv.gz", root=root
        ),
        feature_definitions=project_path(
            "results", "reference_atlas", "tcga_reference_feature_definitions.yaml", root=root
        ),
        scaler_parameters=project_path(
            "results", "reference_atlas", "tcga_reference_scaler_parameters.json", root=root
        ),
        pca_model=project_path("results", "reference_atlas", "tcga_reference_pca_model.pkl", root=root),
        reference_pca_coordinates=project_path(
            "results", "reference_atlas", "tcga_reference_pca_coordinates.tsv", root=root
        ),
        nearest_neighbor_index=project_path(
            "results", "reference_atlas", "tcga_reference_nearest_neighbor_index.pkl", root=root
        ),
        nearest_neighbor_metadata=project_path(
            "results", "reference_atlas", "tcga_reference_nearest_neighbor_metadata.json", root=root
        ),
        clinical_metadata=query / "clinical_metadata.tsv",
        wes_features=query / "wes_features.tsv",
        wts_features=query / "wts_features.tsv",
        output_dir=output,
        clinical_feature_vector=output / "clinical_feature_vector.tsv",
        clinical_feature_vector_scaled=output / "clinical_feature_vector_scaled.tsv",
        clinical_feature_missingness=output / "clinical_feature_missingness.tsv",
        nearest_neighbors=output / "tcga_nearest_neighbors.tsv",
        nearest_project_summary=output / "tcga_nearest_project_summary.tsv",
        nearest_major_group_summary=output / "tcga_nearest_major_group_summary.tsv",
        pca_projection=output / "clinical_tcga_pca_projection.tsv",
        qc_summary=output / "clinical_projection_qc_summary.tsv",
        pca_by_project_figure=output / "clinical_tcga_pca_projection_by_project.pdf",
        pca_by_major_group_figure=output / "clinical_tcga_pca_projection_by_major_group.pdf",
        nearest_neighbor_figure=output / "clinical_nearest_neighbor_barplot.pdf",
        feature_missingness_figure=output / "clinical_feature_missingness.pdf",
        report=output / "clinical_tcga_projection_report.md",
    )


def validate_case_id(case_id: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", case_id) or case_id in {".", ".."}:
        raise ValueError("case_id must contain only letters, numbers, periods, underscores, or hyphens")
    return case_id


def write_tsv(data: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_text(text: str, path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")
    logging.info("Wrote %s", path)


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def ordered_values_hash(values: list[str]) -> str:
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def require_columns(data: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise ValueError(f"{label} lacks required columns: {missing}")


def is_missing(value: Any) -> bool:
    return pd.isna(value) or (isinstance(value, str) and value.strip() == "")


def numeric_value(value: Any, feature_name: str, source_label: str) -> float | None:
    if is_missing(value):
        return None
    converted = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(converted) or not np.isfinite(float(converted)):
        raise ValueError(f"Non-numeric value for reference feature {feature_name} in {source_label}: {value!r}")
    return float(converted)


def load_clinical_metadata(path: Path, case_id: str) -> dict[str, Any]:
    require_file(path, "clinical metadata")
    metadata = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=True)
    require_columns(metadata, CLINICAL_METADATA_REQUIRED_COLUMNS, "clinical metadata")
    matching = metadata[metadata["case_id"].astype(str) == case_id]
    if len(matching) != 1:
        raise ValueError(f"clinical_metadata.tsv must contain exactly one row for case_id={case_id!r}")
    row = matching.iloc[0].to_dict()
    sample_id = "" if is_missing(row.get("sample_id")) else str(row["sample_id"]).strip()
    if not sample_id:
        raise ValueError("Clinical metadata sample_id is blank")
    row["case_id"] = case_id
    row["sample_id"] = sample_id
    for optional in ["differential_diagnosis", "tumor_purity"]:
        row.setdefault(optional, None)
    return row


def load_optional_feature_row(path: Path, sample_id: str, label: str) -> tuple[dict[str, Any] | None, bool]:
    if not path.exists():
        return None, False
    data = pd.read_csv(path, sep="\t", low_memory=False)
    require_columns(data, ["sample_id"], label)
    matching = data[data["sample_id"].astype(str) == sample_id]
    if len(matching) != 1:
        raise ValueError(f"{path.name} must contain exactly one row for sample_id={sample_id!r}")
    return matching.iloc[0].to_dict(), True


def validate_feature_value(value: float, feature: dict[str, Any]) -> None:
    feature_name = feature["name"]
    feature_type = feature.get("feature_type")
    if feature_type == "binary" and not (np.isclose(value, 0.0) or np.isclose(value, 1.0)):
        raise ValueError(f"Binary reference feature {feature_name} must be 0 or 1; found {value}")
    if feature_type == "count":
        if value < 0:
            raise ValueError(f"Count reference feature {feature_name} cannot be negative; found {value}")
        if not np.isclose(value, round(value)):
            raise ValueError(f"Count reference feature {feature_name} must be integer-like; found {value}")


def load_reference_contract(paths: ProjectionPaths) -> dict[str, Any]:
    required = [
        (paths.reference_matrix, "raw TCGA reference matrix"),
        (paths.reference_metadata, "TCGA reference metadata"),
        (paths.reference_scaled_matrix, "scaled TCGA reference matrix"),
        (paths.feature_definitions, "TCGA reference feature definitions"),
        (paths.scaler_parameters, "TCGA reference scaler parameters"),
        (paths.pca_model, "TCGA reference PCA model"),
        (paths.reference_pca_coordinates, "TCGA reference PCA coordinates"),
        (paths.nearest_neighbor_index, "TCGA nearest-neighbor index"),
        (paths.nearest_neighbor_metadata, "TCGA nearest-neighbor metadata"),
    ]
    for path, description in required:
        require_file(path, description)

    definitions = read_yaml(paths.feature_definitions)
    scaler = read_json(paths.scaler_parameters)
    neighbor_metadata = read_json(paths.nearest_neighbor_metadata)
    feature_order = list(scaler.get("feature_order", []))
    if not feature_order:
        raise ValueError("Reference scaler feature_order is empty")
    feature_records = definitions.get("features", [])
    feature_map = {record.get("name"): record for record in feature_records}
    if len(feature_map) != len(feature_records):
        raise ValueError("Reference feature definitions contain duplicate feature names")
    if definitions.get("feature_order") != feature_order:
        raise ValueError("Reference feature definitions and scaler disagree on feature order")
    if neighbor_metadata.get("feature_order") != feature_order:
        raise ValueError("Nearest-neighbor metadata and scaler disagree on feature order")
    if set(feature_order) - set(feature_map):
        raise ValueError("One or more locked reference features lack feature definitions")

    contract_hash = scaler.get("feature_contract_sha256")
    for label, artifact_hash in [
        ("feature definitions", definitions.get("feature_contract_sha256")),
        ("nearest-neighbor metadata", neighbor_metadata.get("feature_contract_sha256")),
    ]:
        if artifact_hash != contract_hash:
            raise ValueError(f"Reference {label} feature contract hash is incompatible with the scaler")

    raw_header = pd.read_csv(paths.reference_matrix, sep="\t", nrows=0)
    expected_raw_columns = list(definitions.get("metadata_columns", [])) + list(scaler.get("raw_feature_order", []))
    if list(raw_header.columns) != expected_raw_columns:
        raise ValueError("Raw TCGA reference matrix columns do not match the locked raw feature contract")

    metadata = pd.read_csv(paths.reference_metadata, sep="\t", low_memory=False)
    scaled = pd.read_csv(paths.reference_scaled_matrix, sep="\t", low_memory=False)
    require_columns(metadata, REFERENCE_METADATA_COLUMNS, "TCGA reference metadata")
    if metadata[REFERENCE_METADATA_COLUMNS].isna().any().any():
        raise ValueError("TCGA reference metadata contain missing required values")
    expected_scaled_columns = ["sample_barcode"] + feature_order
    if list(scaled.columns) != expected_scaled_columns:
        raise ValueError("Scaled TCGA reference matrix columns do not match locked feature order")
    if not scaled["sample_barcode"].equals(metadata["sample_barcode"]):
        raise ValueError("Scaled TCGA reference matrix and metadata use different sample ordering")
    scaled_values = scaled[feature_order].to_numpy(dtype=float)
    if not np.isfinite(scaled_values).all():
        raise ValueError("Scaled TCGA reference matrix contains non-finite values")
    sample_barcodes = metadata["sample_barcode"].astype(str).tolist()
    sample_hash = ordered_values_hash(sample_barcodes)
    if sample_hash != scaler.get("sample_order_sha256"):
        raise ValueError("TCGA reference metadata sample order hash does not match the scaler")
    if neighbor_metadata.get("sample_order_sha256") != sample_hash:
        raise ValueError("TCGA nearest-neighbor metadata sample order hash is incompatible")

    with paths.pca_model.open("rb") as handle:
        pca_package = pickle.load(handle)
    if pca_package.get("feature_order") != feature_order:
        raise ValueError("Saved PCA model feature order is incompatible with the scaler")
    if pca_package.get("feature_contract_sha256") != contract_hash:
        raise ValueError("Saved PCA model feature contract hash is incompatible")
    if pca_package.get("sample_order_sha256") != sample_hash:
        raise ValueError("Saved PCA model sample order hash is incompatible")
    pca_model = pca_package.get("model")
    if pca_model is None or int(getattr(pca_model, "n_features_in_", -1)) != len(feature_order):
        raise ValueError("Saved PCA model does not accept the locked reference feature count")

    reference_pca = pd.read_csv(paths.reference_pca_coordinates, sep="\t", low_memory=False)
    pca_columns = sorted(
        [column for column in reference_pca.columns if re.fullmatch(r"PC\d+", column)],
        key=lambda column: int(column.removeprefix("PC")),
    )
    expected_pca_components = int(getattr(pca_model, "n_components_", len(pca_columns)))
    require_columns(
        reference_pca,
        REFERENCE_METADATA_COLUMNS + [f"PC{i}" for i in range(1, expected_pca_components + 1)],
        "TCGA reference PCA coordinates",
    )
    for column in REFERENCE_METADATA_COLUMNS:
        if not reference_pca[column].astype(str).equals(metadata[column].astype(str)):
            raise ValueError(
                f"TCGA reference PCA coordinates and metadata disagree on {column} or sample ordering"
            )
    if pca_columns != [f"PC{i}" for i in range(1, expected_pca_components + 1)]:
        raise ValueError("TCGA reference PCA coordinate columns do not match the saved PCA model")
    if not np.isfinite(reference_pca[pca_columns].to_numpy(dtype=float)).all():
        raise ValueError("TCGA reference PCA coordinates contain non-finite values")

    neighbor_package = None
    neighbor_load_warning = ""
    try:
        with paths.nearest_neighbor_index.open("rb") as handle:
            neighbor_package = pickle.load(handle)
    except Exception as exc:  # pragma: no cover - depends on external pickle compatibility
        neighbor_load_warning = f"Saved nearest-neighbor index could not be loaded; direct distances will be used: {exc}"

    return {
        "definitions": definitions,
        "feature_map": feature_map,
        "scaler": scaler,
        "neighbor_metadata": neighbor_metadata,
        "metadata": metadata.reset_index(drop=True),
        "scaled": scaled.reset_index(drop=True),
        "scaled_values": scaled_values,
        "pca_package": pca_package,
        "pca_model": pca_model,
        "reference_pca": reference_pca[REFERENCE_METADATA_COLUMNS + pca_columns].reset_index(drop=True),
        "neighbor_package": neighbor_package,
        "neighbor_load_warning": neighbor_load_warning,
        "feature_order": feature_order,
        "feature_contract_sha256": contract_hash,
        "sample_order_sha256": sample_hash,
    }


def feature_sources_from_rows(
    feature_order: list[str],
    wes_row: dict[str, Any] | None,
    wts_row: dict[str, Any] | None,
) -> tuple[dict[str, tuple[Any, str]], list[dict[str, Any]]]:
    expected = set(feature_order)
    supplied: dict[str, tuple[Any, str]] = {}
    extras: list[dict[str, Any]] = []
    for row, source_file in [(wes_row, "wes_features.tsv"), (wts_row, "wts_features.tsv")]:
        if row is None:
            continue
        for column, raw_value in row.items():
            if column in {"case_id", "sample_id", "notes"}:
                continue
            if column not in expected:
                if not is_missing(raw_value):
                    extras.append(
                        {
                            "feature_name": column,
                            "source_file": source_file,
                            "raw_value": raw_value,
                        }
                    )
                continue
            if is_missing(raw_value):
                continue
            if column in supplied:
                previous_value = numeric_value(supplied[column][0], column, supplied[column][1])
                new_value = numeric_value(raw_value, column, source_file)
                if previous_value is None or new_value is None or not np.isclose(previous_value, new_value):
                    raise ValueError(f"Conflicting values supplied for reference feature {column}")
                supplied[column] = (raw_value, f"{supplied[column][1]};{source_file}")
            else:
                supplied[column] = (raw_value, source_file)
    return supplied, extras


def harmonize_clinical_features(
    case_id: str,
    clinical_metadata: dict[str, Any],
    wes_row: dict[str, Any] | None,
    wts_row: dict[str, Any] | None,
    reference: dict[str, Any],
    max_abs_scaled_warning: float,
) -> dict[str, Any]:
    feature_order = reference["feature_order"]
    feature_map = reference["feature_map"]
    scaler = reference["scaler"]
    sample_id = str(clinical_metadata["sample_id"])
    supplied, extras = feature_sources_from_rows(feature_order, wes_row, wts_row)

    if "purity" not in supplied and not is_missing(clinical_metadata.get("tumor_purity")):
        supplied["purity"] = (clinical_metadata["tumor_purity"], "clinical_metadata.tsv:tumor_purity")

    vector_values: dict[str, float] = {}
    scaled_values: dict[str, float] = {}
    audit_rows: list[dict[str, Any]] = []
    processing_warnings: list[str] = []
    for feature_name in feature_order:
        feature = feature_map[feature_name]
        imputation_value = float(scaler["imputation_median"][feature_name])
        mean = float(scaler["scaling_mean"][feature_name])
        std = float(scaler["scaling_std"][feature_name])
        if not np.isfinite(std) or std <= 0:
            raise ValueError(f"Reference scaling standard deviation is invalid for {feature_name}")

        if feature_name in supplied:
            raw_value, source_file = supplied[feature_name]
            value = numeric_value(raw_value, feature_name, source_file)
            if value is None:
                value = imputation_value
                supplied_flag = False
                imputed_flag = True
                status = "imputed_missing_value"
            else:
                validate_feature_value(value, feature)
                supplied_flag = True
                imputed_flag = False
                status = "supplied"
        else:
            raw_value = None
            source_file = "none"
            value = imputation_value
            supplied_flag = False
            imputed_flag = True
            status = "imputed_missing_feature"

        scaled_value = (value - mean) / std
        if not np.isfinite(scaled_value):
            raise ValueError(f"Clinical scaling produced a non-finite value for {feature_name}")
        vector_values[feature_name] = float(value)
        scaled_values[feature_name] = float(scaled_value)
        audit_rows.append(
            {
                "case_id": case_id,
                "sample_id": sample_id,
                "feature_name": feature_name,
                "expected_in_reference": True,
                "feature_class": feature.get("feature_class", "unknown"),
                "feature_type": feature.get("feature_type", "unknown"),
                "clinical_modality": feature.get("clinical_modality", "unknown"),
                "source_file": source_file,
                "supplied": supplied_flag,
                "raw_value": raw_value,
                "imputed": imputed_flag,
                "imputation_value": imputation_value,
                "value_used": float(value),
                "scaled_value": float(scaled_value),
                "status": status,
                "notes": "Exact reference feature name used" if supplied_flag else "TCGA reference median used",
            }
        )

    for extra in extras:
        audit_rows.append(
            {
                "case_id": case_id,
                "sample_id": sample_id,
                "feature_name": extra["feature_name"],
                "expected_in_reference": False,
                "feature_class": "extra_clinical_feature",
                "feature_type": "not_evaluated",
                "clinical_modality": "not_evaluated",
                "source_file": extra["source_file"],
                "supplied": True,
                "raw_value": extra["raw_value"],
                "imputed": False,
                "imputation_value": np.nan,
                "value_used": np.nan,
                "scaled_value": np.nan,
                "status": "extra_not_in_reference_ignored",
                "notes": "Preserved for audit but excluded from PCA and nearest-neighbor projection",
            }
        )

    feature_audit = pd.DataFrame(audit_rows)
    expected_audit = feature_audit[feature_audit["expected_in_reference"]].copy()
    supplied_count = int(expected_audit["supplied"].sum())
    missing_count = int(expected_audit["imputed"].sum())
    coverage = 100.0 * supplied_count / len(feature_order)
    if missing_count:
        processing_warnings.append(
            f"{missing_count} of {len(feature_order)} locked reference features were median-imputed"
        )
    if extras:
        processing_warnings.append(
            f"{len(extras)} supplied clinical features were not present in the locked atlas and were ignored"
        )
    high_z = expected_audit.loc[expected_audit["scaled_value"].abs() > max_abs_scaled_warning, "feature_name"].tolist()
    if high_z:
        processing_warnings.append(
            "Clinical values exceeded the configured absolute scaled-value warning threshold for: "
            + ", ".join(high_z)
        )
    if coverage < 50:
        processing_warnings.append(
            "Feature coverage is below 50%; similarity results are especially limited"
        )

    feature_vector = pd.DataFrame(
        [{"case_id": case_id, "sample_id": sample_id, **vector_values}],
        columns=["case_id", "sample_id"] + feature_order,
    )
    scaled_vector = pd.DataFrame(
        [{"case_id": case_id, "sample_id": sample_id, **scaled_values}],
        columns=["case_id", "sample_id"] + feature_order,
    )
    if supplied_count == 0:
        raise ValueError("No exact locked reference features were supplied by the clinical query")
    return {
        "feature_vector": feature_vector,
        "scaled_vector": scaled_vector,
        "feature_audit": feature_audit,
        "supplied_count": supplied_count,
        "missing_count": missing_count,
        "coverage_percent": coverage,
        "extra_count": len(extras),
        "warnings": processing_warnings,
    }


def neighbor_package_compatible(
    package: dict[str, Any] | None,
    reference: dict[str, Any],
    metric: str,
) -> tuple[bool, str]:
    if package is None:
        return False, "Saved nearest-neighbor package is unavailable"
    if package.get("feature_order") != reference["feature_order"]:
        return False, "Saved nearest-neighbor package feature order is incompatible"
    if package.get("feature_contract_sha256") != reference["feature_contract_sha256"]:
        return False, "Saved nearest-neighbor package feature contract is incompatible"
    if package.get("sample_order_sha256") != reference["sample_order_sha256"]:
        return False, "Saved nearest-neighbor package sample order is incompatible"
    if package.get("sample_barcodes") != reference["metadata"]["sample_barcode"].astype(str).tolist():
        return False, "Saved nearest-neighbor package barcode ordering is incompatible"
    model = package.get("indices", {}).get(metric)
    if model is None:
        return False, f"Saved nearest-neighbor package does not contain metric={metric}"
    fit_shape = tuple(getattr(model, "_fit_X", np.empty((0, 0))).shape)
    expected_shape = (len(reference["metadata"]), len(reference["feature_order"]))
    if fit_shape != expected_shape:
        return False, "Saved nearest-neighbor fitted matrix shape is incompatible"
    return True, ""


def direct_distances(reference_values: np.ndarray, query: np.ndarray, metric: str) -> np.ndarray:
    if metric == "euclidean":
        return np.linalg.norm(reference_values - query, axis=1)
    if metric == "cosine":
        query_norm = float(np.linalg.norm(query))
        if np.isclose(query_norm, 0.0):
            raise ValueError("Cosine distance is undefined for an all-zero scaled clinical vector")
        reference_norms = np.linalg.norm(reference_values, axis=1)
        similarities = np.zeros(len(reference_values), dtype=float)
        nonzero = reference_norms > 0
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            dot_products = np.sum(reference_values[nonzero] * query[np.newaxis, :], axis=1, dtype=np.float64)
            similarities[nonzero] = dot_products / (reference_norms[nonzero] * query_norm)
        if not np.isfinite(similarities[nonzero]).all():
            raise ValueError("Direct cosine calculation produced non-finite similarities")
        similarities = np.clip(similarities, -1.0, 1.0)
        return 1.0 - similarities
    raise ValueError(f"Unsupported nearest-neighbor metric: {metric}")


def find_nearest_neighbors(
    query_scaled: np.ndarray,
    reference: dict[str, Any],
    top_k: int,
    metric: str,
    force_direct_distance: bool,
) -> tuple[np.ndarray, np.ndarray, str, list[str]]:
    if top_k < 1:
        raise ValueError("top_k must be at least 1")
    top_k = min(int(top_k), len(reference["metadata"]))
    query_warnings: list[str] = []
    compatible, reason = neighbor_package_compatible(reference["neighbor_package"], reference, metric)
    if compatible and not force_direct_distance:
        try:
            model = reference["neighbor_package"]["indices"][metric]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                distances, indices = model.kneighbors(query_scaled.reshape(1, -1), n_neighbors=top_k)
            if not np.isfinite(distances).all():
                raise ValueError("saved nearest-neighbor index returned non-finite distances")
            return (
                distances[0].astype(float),
                indices[0].astype(int),
                "saved_nearest_neighbor_index",
                query_warnings,
            )
        except Exception as exc:  # pragma: no cover - environment-specific sklearn failures
            query_warnings.append(f"Saved nearest-neighbor query failed; direct distances used: {exc}")
    else:
        if force_direct_distance:
            query_warnings.append("Direct distance calculation was explicitly requested")
        elif reason:
            query_warnings.append(reason + "; direct distances used")

    distances_all = direct_distances(reference["scaled_values"], query_scaled, metric)
    order = np.argsort(distances_all, kind="stable")[:top_k]
    return distances_all[order].astype(float), order.astype(int), "direct_scaled_matrix_distance", query_warnings


def build_neighbor_table(
    case_id: str,
    sample_id: str,
    distances: np.ndarray,
    indices: np.ndarray,
    metadata: pd.DataFrame,
    metric: str,
    method: str,
) -> pd.DataFrame:
    selected = metadata.iloc[indices].reset_index(drop=True)
    similarity = 1.0 / (1.0 + distances)
    notes = (
        f"Research similarity only; metric={metric}; method={method}; "
        "similarity_score=1/(1+distance), not a diagnostic probability"
    )
    output = pd.DataFrame(
        {
            "case_id": case_id,
            "query_sample_id": sample_id,
            "rank": np.arange(1, len(selected) + 1),
            "tcga_sample_barcode": selected["sample_barcode"],
            "tcga_project_id": selected["project_id"],
            "tcga_project_code": selected["project_code"],
            "tcga_disease_type": selected["disease_type"],
            "tcga_primary_site": selected["primary_site"],
            "tcga_major_cancer_group": selected["major_cancer_group"],
            "distance": distances,
            "similarity_score": similarity,
            "distance_metric": metric,
            "notes": notes,
        }
    )
    return output[NEIGHBOR_COLUMNS]


def summarize_neighbors_by_project(neighbors: pd.DataFrame) -> pd.DataFrame:
    total_weight = float(neighbors["similarity_score"].sum())
    summary = (
        neighbors.groupby(["tcga_project_code", "tcga_project_id"], as_index=False)
        .agg(
            n_neighbors=("tcga_sample_barcode", "size"),
            median_distance=("distance", "median"),
            min_distance=("distance", "min"),
            similarity_weight=("similarity_score", "sum"),
        )
        .rename(columns={"tcga_project_code": "project_code", "tcga_project_id": "project_id"})
    )
    summary["weighted_similarity_score"] = summary["similarity_weight"] / total_weight if total_weight else 0.0
    summary = summary.sort_values(
        ["weighted_similarity_score", "min_distance", "project_code"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    summary["rank"] = np.arange(1, len(summary) + 1)
    summary["interpretation_note"] = (
        "Descriptive share of top-k molecular-similarity weight; not a tissue-of-origin probability or diagnosis"
    )
    return summary[
        [
            "project_code",
            "project_id",
            "n_neighbors",
            "median_distance",
            "min_distance",
            "weighted_similarity_score",
            "rank",
            "interpretation_note",
        ]
    ]


def summarize_neighbors_by_major_group(neighbors: pd.DataFrame) -> pd.DataFrame:
    total_weight = float(neighbors["similarity_score"].sum())
    summary = neighbors.groupby("tcga_major_cancer_group", as_index=False).agg(
        n_neighbors=("tcga_sample_barcode", "size"),
        median_distance=("distance", "median"),
        min_distance=("distance", "min"),
        similarity_weight=("similarity_score", "sum"),
    )
    summary = summary.rename(columns={"tcga_major_cancer_group": "major_cancer_group"})
    summary["weighted_similarity_score"] = summary["similarity_weight"] / total_weight if total_weight else 0.0
    summary = summary.sort_values(
        ["weighted_similarity_score", "min_distance", "major_cancer_group"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    summary["rank"] = np.arange(1, len(summary) + 1)
    summary["interpretation_note"] = (
        "Descriptive share of top-k molecular-similarity weight; not a lineage probability or diagnosis"
    )
    return summary[
        [
            "major_cancer_group",
            "n_neighbors",
            "median_distance",
            "min_distance",
            "weighted_similarity_score",
            "rank",
            "interpretation_note",
        ]
    ]


def project_pca(
    case_id: str,
    sample_id: str,
    query_scaled: np.ndarray,
    reference: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = reference["pca_model"]
    # NumPy 2.2 may emit spurious matmul RuntimeWarnings on macOS while returning
    # finite PCA coordinates. Suppress those warnings and enforce finiteness below.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        clinical_coordinates = model.transform(query_scaled.reshape(1, -1))
    if not np.isfinite(clinical_coordinates).all():
        raise ValueError("Saved PCA projection produced non-finite coordinates")
    pc_columns = [f"PC{i + 1}" for i in range(clinical_coordinates.shape[1])]
    clinical = pd.DataFrame(clinical_coordinates, columns=pc_columns)
    clinical.insert(0, "query_sample_id", sample_id)
    clinical.insert(0, "case_id", case_id)
    clinical["projection_note"] = "Projection with saved TCGA PCA model; PCA was not refit"
    reference_pca = reference["reference_pca"].copy()
    return clinical, reference_pca


def color_mapping(values: pd.Series) -> dict[str, Any]:
    labels = sorted(values.fillna("Unannotated").astype(str).unique())
    cmap = plt.get_cmap("nipy_spectral", max(len(labels), 1))
    return {label: cmap(i) for i, label in enumerate(labels)}


def plot_pca_projection(
    reference_pca: pd.DataFrame,
    clinical_pca: pd.DataFrame,
    color_column: str,
    path: Path,
    title: str,
    explained_variance_ratio: list[float],
) -> None:
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(13, 8))
    plot_data = reference_pca.copy()
    plot_data[color_column] = plot_data[color_column].fillna("Unannotated").astype(str)
    colors = color_mapping(plot_data[color_column])
    for label, group in plot_data.groupby(color_column, sort=True):
        ax.scatter(
            group["PC1"],
            group["PC2"],
            s=8,
            alpha=0.38,
            color=colors[label],
            label=label,
            linewidths=0,
            rasterized=True,
        )
    ax.scatter(
        clinical_pca["PC1"].iloc[0],
        clinical_pca["PC2"].iloc[0],
        marker="*",
        s=260,
        color="#D62728",
        edgecolor="#111111",
        linewidth=0.9,
        label="_nolegend_",
        zorder=10,
    )
    pc1_variance = 100 * explained_variance_ratio[0] if explained_variance_ratio else 0
    pc2_variance = 100 * explained_variance_ratio[1] if len(explained_variance_ratio) > 1 else 0
    ax.set_xlabel(f"TCGA reference PC1 ({pc1_variance:.1f}% variance)")
    ax.set_ylabel(f"TCGA reference PC2 ({pc2_variance:.1f}% variance)")
    ax.set_title(title)
    legend_columns = 2 if plot_data[color_column].nunique() > 18 else 1
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=colors[label],
            markeredgecolor="none",
            markersize=5,
            label=label,
        )
        for label in sorted(colors)
    ]
    legend_handles.append(
        Line2D(
            [0],
            [0],
            marker="*",
            linestyle="none",
            markerfacecolor="#D62728",
            markeredgecolor="#111111",
            markersize=11,
            label="Clinical query",
        )
    )
    ax.legend(
        handles=legend_handles,
        fontsize=7,
        ncol=legend_columns,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
    )
    ax.grid(color="#D9D9D9", linewidth=0.4, alpha=0.5)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_nearest_neighbors(neighbors: pd.DataFrame, path: Path, max_neighbors: int = 15) -> None:
    ensure_dir(path.parent)
    plotted = neighbors.head(max_neighbors).copy().sort_values("rank", ascending=False)
    labels = plotted.apply(
        lambda row: f"#{int(row['rank'])} {row['tcga_sample_barcode']} ({row['tcga_project_code']})", axis=1
    )
    fig, ax = plt.subplots(figsize=(11, max(6, 0.42 * len(plotted))))
    bars = ax.barh(labels, plotted["similarity_score"], color="#356A8A")
    ax.bar_label(bars, labels=[f"{value:.3f}" for value in plotted["similarity_score"]], padding=3, fontsize=8)
    ax.set_xlabel("Similarity score, 1 / (1 + distance)")
    ax.set_ylabel("")
    ax.set_title("Nearest TCGA reference samples for the clinical query")
    ax.set_xlim(0, min(1.0, max(0.1, float(plotted["similarity_score"].max()) * 1.18)))
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_feature_missingness(feature_audit: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    expected = feature_audit[feature_audit["expected_in_reference"]].copy()
    summary = (
        expected.groupby("feature_class", as_index=False)
        .agg(expected=("feature_name", "size"), supplied=("supplied", "sum"), imputed=("imputed", "sum"))
        .sort_values("expected", ascending=True)
    )
    fig, ax = plt.subplots(figsize=(11, max(6, 0.65 * len(summary))))
    ax.barh(summary["feature_class"], summary["supplied"], color="#3B7A57", label="Supplied")
    ax.barh(
        summary["feature_class"],
        summary["imputed"],
        left=summary["supplied"],
        color="#D9A441",
        label="TCGA median-imputed",
    )
    for y, row in enumerate(summary.itertuples(index=False)):
        ax.text(float(row.expected) + 0.3, y, f"{int(row.supplied)}/{int(row.expected)}", va="center", fontsize=8)
    coverage = 100 * float(expected["supplied"].sum()) / len(expected)
    ax.set_xlabel("Number of locked reference features")
    ax.set_ylabel("")
    ax.set_title(f"Clinical query feature coverage by class ({coverage:.1f}% overall)")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.set_xlim(0, max(summary["expected"].max() + 4, 10))
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def build_qc_summary(
    case_id: str,
    sample_id: str,
    harmonized: dict[str, Any],
    reference: dict[str, Any],
    wes_present: bool,
    wts_present: bool,
    top_k: int,
    metric: str,
    method: str,
    warnings: list[str],
) -> pd.DataFrame:
    expected = harmonized["feature_audit"][harmonized["feature_audit"]["expected_in_reference"]]
    wes_expected = int((expected["clinical_modality"] == "WES").sum())
    wts_expected = int((expected["clinical_modality"] == "WTS").sum())
    wes_supplied = int(((expected["clinical_modality"] == "WES") & expected["supplied"]).sum())
    wts_supplied = int(((expected["clinical_modality"] == "WTS") & expected["supplied"]).sum())
    coverage = float(harmonized["coverage_percent"])
    if harmonized["missing_count"] == 0:
        status = "success_complete_feature_coverage"
    elif coverage >= 50:
        status = "success_with_reference_median_imputation"
    else:
        status = "success_low_coverage_with_reference_median_imputation"
    rows = [
        ("global", "case_id", case_id),
        ("global", "query_sample_id", sample_id),
        ("global", "feature_contract_sha256", reference["feature_contract_sha256"]),
        ("global", "number_of_expected_reference_features", len(reference["feature_order"])),
        ("global", "atlas_feature_count_used", len(reference["feature_order"])),
        ("global", "number_of_clinical_features_supplied", harmonized["supplied_count"]),
        ("global", "number_of_missing_clinical_features", harmonized["missing_count"]),
        ("global", "percent_feature_coverage", round(coverage, 3)),
        ("global", "number_of_extra_clinical_features_ignored", harmonized["extra_count"]),
        ("global", "WES_features_present", wes_present),
        ("global", "WTS_features_present", wts_present),
        ("WES", "number_of_expected_reference_features", wes_expected),
        ("WES", "number_of_clinical_features_supplied", wes_supplied),
        ("WTS", "number_of_expected_reference_features", wts_expected),
        ("WTS", "number_of_clinical_features_supplied", wts_supplied),
        ("global", "nearest_neighbor_metric", metric),
        ("global", "nearest_neighbor_method", method),
        ("global", "top_k_neighbors_returned", top_k),
        ("global", "PCA_model_reused_without_refit", True),
        ("global", "projection_status", status),
        ("global", "warnings", "; ".join(warnings) if warnings else "none"),
        (
            "global",
            "limitations",
            "Research/prototype molecular projection only; not a validated clinical diagnostic classifier, not a final tissue-of-origin diagnosis, and not a supervised CUP classifier.",
        ),
    ]
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def markdown_escape(value: Any) -> str:
    if is_missing(value):
        return "Not supplied"
    return str(value).replace("|", "\\|").replace("\n", " ")


def markdown_table(data: pd.DataFrame, columns: list[tuple[str, str]]) -> str:
    headers = [label for _, label in columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in data.iterrows():
        values = []
        for column, _label in columns:
            value = row[column]
            if isinstance(value, float):
                values.append(f"{value:.4f}")
            else:
                values.append(markdown_escape(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_report(
    clinical_metadata: dict[str, Any],
    harmonized: dict[str, Any],
    project_summary: pd.DataFrame,
    group_summary: pd.DataFrame,
    pca_projection: pd.DataFrame,
    metric: str,
    method: str,
    warnings: list[str],
) -> str:
    case_id = str(clinical_metadata["case_id"])
    sample_id = str(clinical_metadata["sample_id"])
    expected = harmonized["feature_audit"][harmonized["feature_audit"]["expected_in_reference"]]
    modality_rows = []
    for modality in ["WES", "WTS"]:
        subset = expected[expected["clinical_modality"] == modality]
        supplied = int(subset["supplied"].sum())
        modality_rows.append(
            {
                "modality": modality,
                "expected": len(subset),
                "supplied": supplied,
                "imputed": int(subset["imputed"].sum()),
                "coverage": 100.0 * supplied / len(subset) if len(subset) else 0.0,
            }
        )
    modality_rows.append(
        {
            "modality": "Overall",
            "expected": len(expected),
            "supplied": harmonized["supplied_count"],
            "imputed": harmonized["missing_count"],
            "coverage": harmonized["coverage_percent"],
        }
    )
    coverage_table = markdown_table(
        pd.DataFrame(modality_rows),
        [
            ("modality", "Feature set"),
            ("expected", "Expected"),
            ("supplied", "Supplied"),
            ("imputed", "Imputed"),
            ("coverage", "Coverage (%)"),
        ],
    )
    top_projects = markdown_table(
        project_summary.head(10),
        [
            ("rank", "Rank"),
            ("project_code", "TCGA project"),
            ("n_neighbors", "Neighbors"),
            ("median_distance", "Median distance"),
            ("min_distance", "Minimum distance"),
            ("weighted_similarity_score", "Similarity-weight share"),
        ],
    )
    top_groups = markdown_table(
        group_summary.head(10),
        [
            ("rank", "Rank"),
            ("major_cancer_group", "Major cancer group"),
            ("n_neighbors", "Neighbors"),
            ("median_distance", "Median distance"),
            ("weighted_similarity_score", "Similarity-weight share"),
        ],
    )
    warning_lines = "\n".join(f"- {warning}" for warning in warnings) if warnings else "- No additional processing warnings."
    metadata_lines = "\n".join(
        [
            f"- Case ID: `{markdown_escape(case_id)}`",
            f"- Sample ID: `{markdown_escape(sample_id)}`",
            f"- Specimen type: {markdown_escape(clinical_metadata.get('specimen_type'))}",
            f"- Submitted diagnosis: {markdown_escape(clinical_metadata.get('submitted_diagnosis'))}",
            f"- Differential diagnosis: {markdown_escape(clinical_metadata.get('differential_diagnosis'))}",
            f"- Tumor purity metadata: {markdown_escape(clinical_metadata.get('tumor_purity'))}",
            f"- Notes: {markdown_escape(clinical_metadata.get('notes'))}",
        ]
    )
    pc1 = float(pca_projection["PC1"].iloc[0])
    pc2 = float(pca_projection["PC2"].iloc[0])
    return f"""# TCGA Molecular Projection Report: {case_id}

Generated by `scripts/21_project_clinical_case_to_tcga.py`.

## Research-Use Statement

This is a research/prototype molecular projection against the TCGA sample-level reference atlas. It is not a validated clinical diagnostic classifier, does not establish tissue of origin, and must not be used as a final cancer diagnosis. No supervised cancer-of-unknown-primary classifier was trained or applied.

## Clinical Metadata

{metadata_lines}

## Feature Completeness

Clinical inputs were matched by exact locked reference feature names. Missing locked features were filled with TCGA reference medians before applying the saved reference means and standard deviations. Imputed values support numerical projection but do not constitute observed molecular evidence.

Locked TCGA atlas features used: {len(expected)}.

{coverage_table}

Extra supplied clinical features not present in the locked atlas: {harmonized['extra_count']}. These values were preserved in `clinical_feature_missingness.tsv` and excluded from distance and PCA calculations.

## Projection Method

- Nearest-neighbor metric: `{metric}`
- Neighbor calculation: `{method}`
- PCA: saved TCGA reference PCA model reused without refitting
- Similarity score: `1 / (1 + distance)`; this is a display transformation, not a probability

## PCA Projection Summary

- PC1 coordinate: `{pc1:.4f}`
- PC2 coordinate: `{pc2:.4f}`

These coordinates were generated with the saved TCGA PCA model without refitting. They locate the query in the reference feature space but are not diagnostic scores or tissue-of-origin probabilities.

## Top Nearest TCGA Projects

{top_projects}

Project ranks summarize the composition and distance of the returned top-k neighbors. They are descriptive molecular-similarity results, not tissue-of-origin probabilities.

## Top Nearest Major Cancer Groups

{top_groups}

## Warnings

{warning_lines}

## Figures

- `clinical_tcga_pca_projection_by_project.pdf`
- `clinical_tcga_pca_projection_by_major_group.pdf`
- `clinical_nearest_neighbor_barplot.pdf`
- `clinical_feature_missingness.pdf`

## Caveats

- The script ingests already-computed harmonized WES/WTS features; it does not process FASTQ, BAM, VCF, or raw RNA-seq files.
- Mutation counts are TCGA-compatible mutation-count proxies rather than callable-territory-normalized TMB.
- TCGA is enriched for primary tumors and differs from metastatic, small-biopsy, cytology, and real-world clinical specimens.
- Missing features are median-imputed, which can reduce or distort similarity when an assay modality is incomplete.
- Expression pathway scores and expression PCs require preprocessing conventions compatible with the TCGA reference.
- The top-k summaries are sensitive to feature coverage, assay harmonization, distance metric, and the current feature balance.
- Clinical validation against known-primary cohorts and explicit uncertainty calibration are required before any diagnostic use.

## Conclusion

The query was projected into the locked TCGA molecular reference space for research comparison. The results identify nearby TCGA molecular contexts for review but do not provide a clinical diagnosis or a validated cancer-of-unknown-primary classification.
"""


def project_clinical_case(paths: ProjectionPaths, case_id: str, config: dict[str, Any] = CONFIG) -> dict[str, Any]:
    validate_case_id(case_id)
    reference = load_reference_contract(paths)
    clinical_metadata = load_clinical_metadata(paths.clinical_metadata, case_id)
    sample_id = str(clinical_metadata["sample_id"])
    wes_row, wes_present = load_optional_feature_row(paths.wes_features, sample_id, "WES feature table")
    wts_row, wts_present = load_optional_feature_row(paths.wts_features, sample_id, "WTS feature table")
    if not wes_present and not wts_present:
        raise FileNotFoundError("At least one of wes_features.tsv or wts_features.tsv is required")

    harmonized = harmonize_clinical_features(
        case_id,
        clinical_metadata,
        wes_row,
        wts_row,
        reference,
        float(config["max_abs_scaled_warning"]),
    )
    query_scaled = harmonized["scaled_vector"][reference["feature_order"]].iloc[0].to_numpy(dtype=float)
    metric = str(config.get("metric") or reference["neighbor_metadata"]["default_metric"])
    allowed_metrics = set(reference["neighbor_metadata"].get("available_metrics", [])) & {
        "euclidean",
        "cosine",
    }
    if metric not in allowed_metrics:
        raise ValueError(f"Requested metric is not supported by the reference contract: {metric}")

    distances, indices, method, neighbor_warnings = find_nearest_neighbors(
        query_scaled,
        reference,
        int(config["top_k"]),
        metric,
        bool(config["force_direct_distance"]),
    )
    warnings = list(harmonized["warnings"])
    if reference["neighbor_load_warning"]:
        warnings.append(reference["neighbor_load_warning"])
    warnings.extend(neighbor_warnings)
    neighbors = build_neighbor_table(
        case_id,
        sample_id,
        distances,
        indices,
        reference["metadata"],
        metric,
        method,
    )
    project_summary = summarize_neighbors_by_project(neighbors)
    group_summary = summarize_neighbors_by_major_group(neighbors)
    clinical_pca, reference_pca = project_pca(case_id, sample_id, query_scaled, reference)
    qc = build_qc_summary(
        case_id,
        sample_id,
        harmonized,
        reference,
        wes_present,
        wts_present,
        len(neighbors),
        metric,
        method,
        warnings,
    )
    report = build_report(
        clinical_metadata,
        harmonized,
        project_summary,
        group_summary,
        clinical_pca,
        metric,
        method,
        warnings,
    )

    ensure_dir(paths.output_dir)
    write_tsv(harmonized["feature_vector"], paths.clinical_feature_vector)
    write_tsv(harmonized["scaled_vector"], paths.clinical_feature_vector_scaled)
    write_tsv(harmonized["feature_audit"], paths.clinical_feature_missingness)
    write_tsv(neighbors, paths.nearest_neighbors)
    write_tsv(project_summary, paths.nearest_project_summary)
    write_tsv(group_summary, paths.nearest_major_group_summary)
    write_tsv(clinical_pca, paths.pca_projection)
    write_tsv(qc, paths.qc_summary)
    explained_variance_ratio = list(reference["pca_package"].get("explained_variance_ratio", []))
    plot_pca_projection(
        reference_pca,
        clinical_pca,
        "project_code",
        paths.pca_by_project_figure,
        "Clinical query projected into TCGA PCA space by project",
        explained_variance_ratio,
    )
    plot_pca_projection(
        reference_pca,
        clinical_pca,
        "major_cancer_group",
        paths.pca_by_major_group_figure,
        "Clinical query projected into TCGA PCA space by major cancer group",
        explained_variance_ratio,
    )
    plot_nearest_neighbors(neighbors, paths.nearest_neighbor_figure)
    plot_feature_missingness(harmonized["feature_audit"], paths.feature_missingness_figure)
    write_text(report, paths.report)
    return {
        "clinical_metadata": clinical_metadata,
        "reference": reference,
        "harmonized": harmonized,
        "neighbors": neighbors,
        "project_summary": project_summary,
        "major_group_summary": group_summary,
        "pca_projection": clinical_pca,
        "qc": qc,
        "metric": metric,
        "neighbor_method": method,
        "warnings": warnings,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Project an already-computed clinical WES/WTS feature case into the locked TCGA reference atlas."
    )
    parser.add_argument("--case-id", required=True, help="Clinical query case directory name and case_id value.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument("--query-dir", type=Path, default=None, help="Override data/clinical_queries/<case_id>.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Override results/clinical_projection/<case_id>.")
    parser.add_argument("--top-k", type=int, default=CONFIG["top_k"], help="Number of TCGA neighbors to return.")
    parser.add_argument("--metric", choices=["euclidean", "cosine"], default=None, help="Distance metric override.")
    parser.add_argument(
        "--force-direct-distance",
        action="store_true",
        help="Bypass the saved neighbor index and calculate distances directly from the scaled TCGA matrix.",
    )
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    case_id = validate_case_id(args.case_id)
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    paths = default_paths(root, case_id, args.query_dir, args.output_dir)
    config = dict(CONFIG)
    config.update(
        {
            "top_k": args.top_k,
            "metric": args.metric,
            "force_direct_distance": args.force_direct_distance,
        }
    )
    outputs = project_clinical_case(paths, case_id, config)
    logging.info(
        "Clinical projection complete: case=%s sample=%s coverage=%.3f%% top_project=%s metric=%s method=%s",
        case_id,
        outputs["clinical_metadata"]["sample_id"],
        outputs["harmonized"]["coverage_percent"],
        outputs["project_summary"].iloc[0]["project_code"],
        outputs["metric"],
        outputs["neighbor_method"],
    )


if __name__ == "__main__":
    main()

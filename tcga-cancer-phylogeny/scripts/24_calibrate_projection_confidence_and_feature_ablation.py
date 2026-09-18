#!/usr/bin/env python

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcga_cancer_phylogeny_matplotlib")
)
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.neighbors import NearestNeighbors

from lib.common import configure_logging, ensure_dir, find_project_root, require_file


FEATURE_SET_ORDER = [
    "all_features",
    "mutation_counts_only",
    "driver_indicators_only",
    "driver_counts_only",
    "mutation_features_combined",
    "purity_ploidy_only",
    "copy_number_aneuploidy_only",
    "pathway_scores_only",
    "expression_pcs_only",
    "expression_combined",
    "WES_like_features",
    "WTS_like_features",
    "WES_plus_WTS_combined",
]

FEATURE_CLASS_TO_PRIMARY_GROUP = {
    "mutation_count_proxy": "mutation_counts_only",
    "driver_mutation_binary": "driver_indicators_only",
    "driver_mutation_count": "driver_counts_only",
    "purity_ploidy": "purity_ploidy_only",
    "copy_number_aneuploidy": "copy_number_aneuploidy_only",
    "expression_pathway_score": "pathway_scores_only",
    "expression_pca": "expression_pcs_only",
}

AMBIGUITY_CONFIG = {
    "meaningful_project_score_min": 0.05,
    "uninterpretable_coverage_below": 0.50,
    "low_ambiguity_top1_min": 0.65,
    "low_ambiguity_margin_min": 0.25,
    "low_ambiguity_entropy_max": 0.60,
    "high_ambiguity_top1_below": 0.40,
    "high_ambiguity_margin_below": 0.10,
    "high_ambiguity_entropy_min": 0.80,
    "high_ambiguity_meaningful_projects_min": 8,
}

LEVEL_SPECS = {
    "project": {
        "metadata_column": "project_code",
        "correct_column": "correct_top_project",
        "predicted_column": "predicted_top_project",
    },
    "lineage": {
        "metadata_column": "cup_lineage_group",
        "correct_column": "correct_top_lineage",
        "predicted_column": "predicted_top_lineage",
    },
    "major_group": {
        "metadata_column": "major_cancer_group",
        "correct_column": "correct_top_major_group",
        "predicted_column": "predicted_top_major_group",
    },
}

CALIBRATION_METRICS = {
    "top1_top2_margin": ("margin", True),
    "entropy": ("entropy", False),
    "top1_score": ("top1_score", True),
    "neighbor_concentration": ("neighbor_concentration", True),
}

METADATA_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "disease_type",
    "primary_site",
    "major_cancer_group",
]


@dataclass(frozen=True)
class AnalysisPaths:
    predictions: Path
    overall_metrics: Path
    ambiguity_performance: Path
    scaled_matrix: Path
    metadata: Path
    feature_definitions: Path
    scaler_parameters: Path
    lineage_groups: Path
    cancer_group_map: Path
    output_dir: Path
    thresholds: Path
    calibrated_predictions: Path
    calibration_summary: Path
    feature_groups: Path
    ablation_summary: Path
    ablation_project: Path
    ablation_lineage: Path
    ablation_project_confusions: Path
    ablation_lineage_confusions: Path
    report: Path
    confidence_figure: Path
    margin_figure: Path
    entropy_figure: Path
    ablation_overall_figure: Path
    ablation_lineage_figure: Path
    ablation_project_heatmap: Path
    best_confusion_heatmap: Path


def default_paths(root: Path, output_dir: Path | None = None) -> AnalysisPaths:
    out = output_dir or root / "results" / "validation" / "calibration_feature_ablation"
    return AnalysisPaths(
        predictions=root
        / "results"
        / "validation"
        / "known_primary_projection"
        / "known_primary_projection_predictions.tsv",
        overall_metrics=root
        / "results"
        / "validation"
        / "known_primary_projection"
        / "known_primary_overall_metrics.tsv",
        ambiguity_performance=root
        / "results"
        / "validation"
        / "known_primary_projection"
        / "known_primary_ambiguity_performance.tsv",
        scaled_matrix=root
        / "results"
        / "reference_atlas"
        / "tcga_sample_reference_scaled_matrix.tsv.gz",
        metadata=root / "results" / "reference_atlas" / "tcga_sample_reference_metadata.tsv",
        feature_definitions=root
        / "results"
        / "reference_atlas"
        / "tcga_reference_feature_definitions.yaml",
        scaler_parameters=root
        / "results"
        / "reference_atlas"
        / "tcga_reference_scaler_parameters.json",
        lineage_groups=root / "config" / "cup_lineage_groups.yaml",
        cancer_group_map=root / "config" / "cancer_group_map.csv",
        output_dir=out,
        thresholds=out / "calibrated_confidence_thresholds.tsv",
        calibrated_predictions=out / "known_primary_predictions_calibrated.tsv",
        calibration_summary=out / "calibration_performance_summary.tsv",
        feature_groups=out / "feature_group_definitions.tsv",
        ablation_summary=out / "feature_ablation_performance_summary.tsv",
        ablation_project=out / "feature_ablation_project_performance.tsv",
        ablation_lineage=out / "feature_ablation_lineage_performance.tsv",
        ablation_project_confusions=out / "feature_ablation_project_confusion_matrices.tsv",
        ablation_lineage_confusions=out / "feature_ablation_lineage_confusion_matrices.tsv",
        report=out / "calibration_feature_ablation_report.md",
        confidence_figure=out / "calibration_accuracy_by_confidence_tier.pdf",
        margin_figure=out / "calibration_margin_vs_accuracy.pdf",
        entropy_figure=out / "calibration_entropy_vs_accuracy.pdf",
        ablation_overall_figure=out / "feature_ablation_overall_accuracy.pdf",
        ablation_lineage_figure=out / "feature_ablation_lineage_accuracy.pdf",
        ablation_project_heatmap=out / "feature_ablation_project_accuracy_heatmap.pdf",
        best_confusion_heatmap=out / "feature_ablation_confusion_heatmap_best_model.pdf",
    )


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected YAML mapping in {path}")
    return value


def read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def write_tsv(data: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_text(value: str, path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(value.rstrip() + "\n", encoding="utf-8")
    logging.info("Wrote %s", path)


def load_lineage_configuration(
    path: Path, project_codes: set[str]
) -> tuple[dict[str, str], list[str], dict[str, set[str]]]:
    config = read_yaml(path)
    records = config.get("lineages")
    if not isinstance(records, list) or not records:
        raise ValueError("CUP lineage configuration lacks a non-empty lineages list")
    project_to_lineage: dict[str, str] = {}
    lineage_order: list[str] = []
    lineage_major_groups: dict[str, set[str]] = {}
    for record in records:
        lineage = str(record.get("name", "")).strip()
        projects = [str(value).strip() for value in record.get("projects", [])]
        major_groups = {str(value).strip() for value in record.get("major_cancer_groups", [])}
        if not lineage or not projects or not major_groups:
            raise ValueError("Each CUP lineage must define name, projects, and major_cancer_groups")
        lineage_order.append(lineage)
        lineage_major_groups[lineage] = major_groups
        for project in projects:
            if project in project_to_lineage:
                raise ValueError(f"TCGA project {project} is assigned to multiple CUP lineages")
            project_to_lineage[project] = lineage
    missing = sorted(project_codes - set(project_to_lineage))
    extra = sorted(set(project_to_lineage) - project_codes)
    if missing or extra:
        raise ValueError(f"CUP lineage project mismatch; missing={missing}, extra={extra}")
    return project_to_lineage, lineage_order, lineage_major_groups


def load_inputs(paths: AnalysisPaths) -> dict[str, Any]:
    required = [
        (paths.predictions, "known-primary predictions"),
        (paths.overall_metrics, "known-primary overall metrics"),
        (paths.ambiguity_performance, "known-primary ambiguity summary"),
        (paths.scaled_matrix, "locked scaled reference matrix"),
        (paths.metadata, "reference metadata"),
        (paths.feature_definitions, "reference feature definitions"),
        (paths.scaler_parameters, "reference scaler parameters"),
        (paths.lineage_groups, "CUP lineage configuration"),
        (paths.cancer_group_map, "cancer group map"),
    ]
    for path, label in required:
        require_file(path, label)

    predictions = pd.read_csv(paths.predictions, sep="\t")
    overall_metrics = pd.read_csv(paths.overall_metrics, sep="\t")
    prior_ambiguity = pd.read_csv(paths.ambiguity_performance, sep="\t")
    scaled = pd.read_csv(paths.scaled_matrix, sep="\t")
    metadata = pd.read_csv(paths.metadata, sep="\t")
    definitions = read_yaml(paths.feature_definitions)
    scaler = read_json(paths.scaler_parameters)
    cancer_map = pd.read_csv(paths.cancer_group_map)

    missing_metadata = sorted(set(METADATA_COLUMNS) - set(metadata.columns))
    if missing_metadata or metadata[METADATA_COLUMNS].isna().any().any():
        raise ValueError(f"Reference metadata are incomplete; missing columns={missing_metadata}")
    if metadata["sample_barcode"].duplicated().any():
        raise ValueError("Reference metadata contain duplicate sample barcodes")
    if not metadata["sample_barcode"].astype(str).equals(scaled["sample_barcode"].astype(str)):
        raise ValueError("Scaled matrix row order does not match reference metadata")

    feature_order = [str(value) for value in definitions.get("feature_order", [])]
    records = definitions.get("features")
    if not feature_order or not isinstance(records, list) or len(records) != len(feature_order):
        raise ValueError("Feature definitions lack a complete locked feature order and records")
    if scaled.columns.tolist() != ["sample_barcode", *feature_order]:
        raise ValueError("Scaled matrix columns do not match the locked feature order")
    if list(scaler.get("feature_order", [])) != feature_order:
        raise ValueError("Scaler feature order does not match feature definitions")
    numeric = scaled[feature_order].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("Scaled matrix contains non-finite values")

    project_codes = set(metadata["project_code"].astype(str))
    project_to_lineage, lineage_order, lineage_major_groups = load_lineage_configuration(
        paths.lineage_groups, project_codes
    )
    metadata = metadata.copy()
    metadata["cup_lineage_group"] = metadata["project_code"].astype(str).map(project_to_lineage)

    required_map = {"project_code", "project_id", "broad_group"}
    if not required_map.issubset(cancer_map.columns):
        raise ValueError("Cancer group map lacks project_code, project_id, or broad_group")
    map_rows = cancer_map.set_index("project_code")[["project_id", "broad_group"]].to_dict("index")
    for row in metadata[["project_code", "project_id", "major_cancer_group"]].drop_duplicates().itertuples(
        index=False
    ):
        expected = map_rows.get(str(row.project_code))
        if expected is None or str(expected["project_id"]) != str(row.project_id):
            raise ValueError(f"Cancer group project mapping mismatch for {row.project_code}")
        if str(expected["broad_group"]) != str(row.major_cancer_group):
            raise ValueError(f"Cancer group label mismatch for {row.project_code}")

    required_prediction_columns = {
        "sample_barcode",
        "true_project_code",
        "true_project_id",
        "true_major_cancer_group",
        "true_cup_lineage_group",
        "predicted_top_project",
        "predicted_top_lineage",
        "predicted_top_major_group",
        "top_project_score",
        "top_lineage_score",
        "top_major_group_score",
        "top1_top2_project_margin",
        "top1_top2_lineage_margin",
        "ambiguity_class",
        "correct_top_project",
        "correct_top_lineage",
        "correct_top_major_group",
        "feature_coverage",
    }
    missing_predictions = sorted(required_prediction_columns - set(predictions.columns))
    if missing_predictions or predictions["sample_barcode"].duplicated().any():
        raise ValueError(f"Known-primary predictions are incomplete; missing={missing_predictions}")
    prediction_index = predictions.set_index("sample_barcode", drop=False)
    missing_samples = sorted(set(metadata["sample_barcode"].astype(str)) - set(prediction_index.index.astype(str)))
    if missing_samples:
        raise ValueError(f"Known-primary predictions lack {len(missing_samples)} atlas samples")
    aligned_predictions = prediction_index.loc[metadata["sample_barcode"].astype(str)].reset_index(drop=True)
    expected_labels = pd.DataFrame(
        {
            "true_project_code": metadata["project_code"].astype(str),
            "true_project_id": metadata["project_id"].astype(str),
            "true_major_cancer_group": metadata["major_cancer_group"].astype(str),
            "true_cup_lineage_group": metadata["cup_lineage_group"].astype(str),
        }
    )
    if not all(
        aligned_predictions[column].astype(str).equals(expected_labels[column])
        for column in expected_labels.columns
    ):
        raise ValueError("Known-primary true labels do not match locked reference metadata")

    return {
        "predictions": aligned_predictions,
        "overall_metrics": overall_metrics,
        "prior_ambiguity": prior_ambiguity,
        "matrix": numeric,
        "metadata": metadata,
        "feature_order": feature_order,
        "feature_records": records,
        "project_to_lineage": project_to_lineage,
        "lineage_order": lineage_order,
        "lineage_major_groups": lineage_major_groups,
    }


def build_feature_groups(
    feature_order: list[str], records: list[dict[str, Any]]
) -> tuple[dict[str, list[str]], pd.DataFrame, list[str]]:
    record_by_name = {str(record.get("name")): record for record in records}
    if set(record_by_name) != set(feature_order):
        raise ValueError("Feature records do not map one-to-one to the locked feature order")

    groups: dict[str, list[str]] = {name: [] for name in FEATURE_SET_ORDER}
    unclassified: list[str] = []
    for feature in feature_order:
        record = record_by_name[feature]
        feature_class = str(record.get("feature_class", ""))
        modality = str(record.get("clinical_modality", "")).upper()
        groups["all_features"].append(feature)
        primary = FEATURE_CLASS_TO_PRIMARY_GROUP.get(feature_class)
        if primary is None or modality not in {"WES", "WTS"}:
            unclassified.append(feature)
            continue
        groups[primary].append(feature)
        if primary in {"mutation_counts_only", "driver_indicators_only", "driver_counts_only"}:
            groups["mutation_features_combined"].append(feature)
        if primary in {"pathway_scores_only", "expression_pcs_only"}:
            groups["expression_combined"].append(feature)
        if modality == "WES":
            groups["WES_like_features"].append(feature)
        if modality == "WTS":
            groups["WTS_like_features"].append(feature)
        groups["WES_plus_WTS_combined"].append(feature)

    for group in FEATURE_SET_ORDER:
        groups[group] = [feature for feature in feature_order if feature in set(groups[group])]
        if not groups[group]:
            raise ValueError(f"Required feature group {group} is empty")

    rows: list[dict[str, Any]] = []
    for group in FEATURE_SET_ORDER:
        for feature in groups[group]:
            record = record_by_name[feature]
            rows.append(
                {
                    "feature_set": group,
                    "feature_name": feature,
                    "feature_class": record.get("feature_class"),
                    "clinical_modality": record.get("clinical_modality"),
                    "source_layer": record.get("source_layer"),
                    "n_features_in_set": len(groups[group]),
                    "classification_status": "classified",
                    "assignment_basis": "Locked feature_class and clinical_modality fields",
                    "notes": "Feature values remain in the locked TCGA standardized space.",
                }
            )
    for feature in unclassified:
        record = record_by_name[feature]
        rows.append(
            {
                "feature_set": "unclassified",
                "feature_name": feature,
                "feature_class": record.get("feature_class"),
                "clinical_modality": record.get("clinical_modality"),
                "source_layer": record.get("source_layer"),
                "n_features_in_set": len(unclassified),
                "classification_status": "unclassified",
                "assignment_basis": "No supported feature_class/modality assignment",
                "notes": "Excluded from named modality subsets; retained in no ablation group.",
            }
        )
    return groups, pd.DataFrame(rows), unclassified


def stratified_query_indices(
    metadata: pd.DataFrame, max_samples: int | None, seed: int
) -> np.ndarray:
    all_indices = np.arange(len(metadata), dtype=int)
    if max_samples is None or max_samples >= len(metadata):
        return all_indices
    if max_samples < metadata["project_code"].nunique():
        raise ValueError("max_samples must be at least the number of TCGA projects")
    counts = metadata["project_code"].value_counts().sort_index()
    allocation = pd.Series(1, index=counts.index, dtype=int)
    remaining = max_samples - int(allocation.sum())
    capacity = counts - allocation
    raw = capacity / capacity.sum() * remaining
    extra = np.floor(raw).astype(int).clip(upper=capacity)
    allocation += extra
    remaining -= int(extra.sum())
    fractional = (raw - np.floor(raw)).sort_values(ascending=False, kind="stable")
    while remaining > 0:
        progressed = False
        for project in fractional.index:
            if allocation[project] < counts[project]:
                allocation[project] += 1
                remaining -= 1
                progressed = True
                if remaining == 0:
                    break
        if not progressed:
            break
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    for project, group in metadata.groupby("project_code", sort=True):
        selected.extend(
            rng.choice(group.index.to_numpy(dtype=int), size=int(allocation[project]), replace=False).tolist()
        )
    return np.array(sorted(selected), dtype=int)


def compute_neighbor_indices(
    matrix: np.ndarray,
    query_indices: np.ndarray,
    top_k: int,
    metric: str,
    batch_size: int = 1024,
) -> tuple[np.ndarray, np.ndarray]:
    if matrix.ndim != 2 or matrix.shape[1] == 0 or not np.isfinite(matrix).all():
        raise ValueError("Ablation matrix must be finite, two-dimensional, and non-empty")
    if top_k < 1 or top_k >= len(matrix):
        raise ValueError("top_k must be between 1 and n_reference_samples - 1")
    if metric not in {"euclidean", "cosine"}:
        raise ValueError("metric must be euclidean or cosine")
    n_requested = min(len(matrix), top_k + 8)
    model = NearestNeighbors(
        metric=metric, algorithm="brute", n_jobs=-1, n_neighbors=n_requested
    ).fit(matrix)
    result_indices = np.empty((len(query_indices), top_k), dtype=int)
    result_distances = np.empty((len(query_indices), top_k), dtype=float)
    for start in range(0, len(query_indices), batch_size):
        stop = min(start + batch_size, len(query_indices))
        batch_queries = query_indices[start:stop]
        distances, indices = model.kneighbors(matrix[batch_queries], return_distance=True)
        for offset, query_index in enumerate(batch_queries):
            keep = indices[offset] != query_index
            filtered_indices = indices[offset][keep][:top_k]
            filtered_distances = distances[offset][keep][:top_k]
            if len(filtered_indices) != top_k:
                raise RuntimeError(f"Could not retain {top_k} non-self neighbors for row {query_index}")
            result_indices[start + offset] = filtered_indices
            result_distances[start + offset] = filtered_distances
        logging.info("Computed ablation neighbors for %d/%d queries", stop, len(query_indices))
    return result_indices, result_distances


def score_neighbor_labels(
    labels: np.ndarray,
    query_indices: np.ndarray,
    neighbor_indices: np.ndarray,
    neighbor_distances: np.ndarray,
    category_order: list[str] | None = None,
) -> dict[str, Any]:
    labels = labels.astype(str)
    categories = category_order or sorted(set(labels))
    category_to_code = {label: index for index, label in enumerate(categories)}
    codes = np.array([category_to_code[label] for label in labels], dtype=int)
    neighbor_codes = codes[neighbor_indices]
    weights = 1.0 / (1.0 + neighbor_distances)
    n_queries, top_k = neighbor_indices.shape
    n_categories = len(categories)
    scores = np.zeros((n_queries, n_categories), dtype=float)
    counts = np.zeros((n_queries, n_categories), dtype=int)
    row_indices = np.repeat(np.arange(n_queries), top_k)
    np.add.at(scores, (row_indices, neighbor_codes.ravel()), weights.ravel())
    np.add.at(counts, (row_indices, neighbor_codes.ravel()), 1)
    score_sums = scores.sum(axis=1, keepdims=True)
    scores = np.divide(scores, score_sums, out=np.zeros_like(scores), where=score_sums > 0)

    label_order = np.broadcast_to(np.arange(n_categories), scores.shape)
    order = np.lexsort((label_order, -counts, -scores), axis=1)
    top_codes = order[:, 0]
    second_codes = order[:, 1] if n_categories > 1 else np.full(n_queries, -1, dtype=int)
    top_scores = scores[np.arange(n_queries), top_codes]
    second_scores = (
        scores[np.arange(n_queries), second_codes] if n_categories > 1 else np.zeros(n_queries)
    )
    inverse_order = np.argsort(order, axis=1)
    true_codes = codes[query_indices]
    positive = counts > 0
    n_positive = positive.sum(axis=1)
    true_present = positive[np.arange(n_queries), true_codes]
    true_ranks = np.where(
        true_present,
        inverse_order[np.arange(n_queries), true_codes] + 1,
        n_positive + 1,
    )
    entropy_terms = np.zeros_like(scores)
    entropy_terms[positive] = -scores[positive] * np.log(scores[positive])
    entropy = np.zeros(n_queries, dtype=float)
    multi = n_positive > 1
    entropy[multi] = entropy_terms[multi].sum(axis=1) / np.log(n_positive[multi])
    concentration = counts[np.arange(n_queries), top_codes] / top_k

    return {
        "categories": categories,
        "scores": scores,
        "counts": counts,
        "order": order,
        "top_label": np.array(categories, dtype=object)[top_codes].astype(str),
        "top1_score": top_scores,
        "top2_score": second_scores,
        "margin": top_scores - second_scores,
        "entropy": entropy,
        "neighbor_concentration": concentration,
        "true_rank": true_ranks,
        "correct": top_codes == true_codes,
    }


def classify_ambiguity_vectorized(
    project_scores: dict[str, Any],
    lineage_scores: dict[str, Any],
    major_scores: dict[str, Any],
    feature_coverage: np.ndarray,
    lineage_major_groups: dict[str, set[str]],
) -> np.ndarray:
    cfg = AMBIGUITY_CONFIG
    meaningful_projects = (
        project_scores["scores"] >= float(cfg["meaningful_project_score_min"])
    ).sum(axis=1)
    agreement = np.array(
        [
            major in lineage_major_groups[lineage]
            for major, lineage in zip(
                major_scores["top_label"], lineage_scores["top_label"], strict=True
            )
        ],
        dtype=bool,
    )
    top1 = lineage_scores["top1_score"]
    margin = lineage_scores["margin"]
    entropy = lineage_scores["entropy"]
    result = np.full(len(top1), "moderate_ambiguity", dtype=object)
    uninterpretable = (feature_coverage < float(cfg["uninterpretable_coverage_below"])) | (
        top1 <= 0
    )
    low = (
        (top1 >= float(cfg["low_ambiguity_top1_min"]))
        & (margin >= float(cfg["low_ambiguity_margin_min"]))
        & (entropy <= float(cfg["low_ambiguity_entropy_max"]))
        & agreement
    )
    high = (
        (top1 < float(cfg["high_ambiguity_top1_below"]))
        | (margin < float(cfg["high_ambiguity_margin_below"]))
        | (entropy >= float(cfg["high_ambiguity_entropy_min"]))
        | (meaningful_projects >= int(cfg["high_ambiguity_meaningful_projects_min"]))
    )
    result[high] = "high_ambiguity"
    result[low] = "low_ambiguity"
    result[uninterpretable] = "uninterpretable_due_to_missing_features"
    return result.astype(str)


def evaluate_feature_matrix(
    matrix: np.ndarray,
    metadata: pd.DataFrame,
    query_indices: np.ndarray,
    feature_coverage: np.ndarray,
    lineage_order: list[str],
    lineage_major_groups: dict[str, set[str]],
    top_k: int,
    metric: str,
) -> pd.DataFrame:
    neighbor_indices, neighbor_distances = compute_neighbor_indices(
        matrix, query_indices, top_k=top_k, metric=metric
    )
    project = score_neighbor_labels(
        metadata["project_code"].astype(str).to_numpy(),
        query_indices,
        neighbor_indices,
        neighbor_distances,
        sorted(metadata["project_code"].astype(str).unique()),
    )
    lineage = score_neighbor_labels(
        metadata["cup_lineage_group"].astype(str).to_numpy(),
        query_indices,
        neighbor_indices,
        neighbor_distances,
        lineage_order,
    )
    major = score_neighbor_labels(
        metadata["major_cancer_group"].astype(str).to_numpy(),
        query_indices,
        neighbor_indices,
        neighbor_distances,
        sorted(metadata["major_cancer_group"].astype(str).unique()),
    )
    query_metadata = metadata.iloc[query_indices].reset_index(drop=True)
    ambiguity = classify_ambiguity_vectorized(
        project,
        lineage,
        major,
        feature_coverage[query_indices],
        lineage_major_groups,
    )
    result = pd.DataFrame(
        {
            "sample_barcode": query_metadata["sample_barcode"].astype(str),
            "true_project_code": query_metadata["project_code"].astype(str),
            "true_cup_lineage_group": query_metadata["cup_lineage_group"].astype(str),
            "true_major_cancer_group": query_metadata["major_cancer_group"].astype(str),
            "predicted_top_project": project["top_label"],
            "predicted_top_lineage": lineage["top_label"],
            "predicted_top_major_group": major["top_label"],
            "correct_top_project": project["correct"],
            "correct_top_lineage": lineage["correct"],
            "correct_top_major_group": major["correct"],
            "true_project_rank": project["true_rank"],
            "true_lineage_rank": lineage["true_rank"],
            "true_major_group_rank": major["true_rank"],
            "true_project_in_top3": project["true_rank"] <= 3,
            "true_project_in_top5": project["true_rank"] <= 5,
            "true_lineage_in_top3": lineage["true_rank"] <= 3,
            "heuristic_ambiguity_class": ambiguity,
            "feature_coverage": feature_coverage[query_indices],
        }
    )
    for level, scored in [("project", project), ("lineage", lineage), ("major_group", major)]:
        result[f"{level}_top1_score"] = scored["top1_score"]
        result[f"{level}_margin"] = scored["margin"]
        result[f"{level}_entropy"] = scored["entropy"]
        result[f"{level}_neighbor_concentration"] = scored["neighbor_concentration"]
    return result


def validate_all_feature_reproduction(
    evaluation: pd.DataFrame, source_predictions: pd.DataFrame
) -> None:
    source = source_predictions.set_index("sample_barcode").loc[evaluation["sample_barcode"]].reset_index()
    comparisons = [
        ("predicted_top_project", "predicted_top_project"),
        ("predicted_top_lineage", "predicted_top_lineage"),
        ("predicted_top_major_group", "predicted_top_major_group"),
    ]
    for calculated, expected in comparisons:
        if not evaluation[calculated].astype(str).equals(source[expected].astype(str)):
            mismatches = int(
                (evaluation[calculated].astype(str) != source[expected].astype(str)).sum()
            )
            raise ValueError(
                f"Vectorized all-feature predictions differ from locked validation for {mismatches} {calculated} rows"
            )
    numeric = [
        ("project_top1_score", "top_project_score"),
        ("lineage_top1_score", "top_lineage_score"),
        ("major_group_top1_score", "top_major_group_score"),
        ("project_margin", "top1_top2_project_margin"),
        ("lineage_margin", "top1_top2_lineage_margin"),
        ("lineage_entropy", "neighbor_entropy"),
    ]
    for calculated, expected in numeric:
        if not np.allclose(
            evaluation[calculated].to_numpy(dtype=float),
            source[expected].to_numpy(dtype=float),
            atol=1e-10,
            rtol=1e-10,
        ):
            raise ValueError(f"Vectorized all-feature metric {calculated} does not reproduce source output")


def deterministic_calibration_split(
    predictions: pd.DataFrame, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    calibration = np.zeros(len(predictions), dtype=bool)
    rng = np.random.default_rng(seed)
    for _, group in predictions.groupby("true_project_code", sort=True):
        indices = group.index.to_numpy(dtype=int).copy()
        rng.shuffle(indices)
        split = max(1, len(indices) // 2)
        calibration[indices[:split]] = True
    evaluation = ~calibration
    if calibration.sum() == 0 or evaluation.sum() == 0:
        raise ValueError("Calibration/evaluation split produced an empty partition")
    return calibration, evaluation


def _tail_accuracy(
    reliability: np.ndarray, correct: np.ndarray, cutoff: float, favorable: bool
) -> tuple[int, float]:
    selected = reliability >= cutoff if favorable else reliability <= cutoff
    n = int(selected.sum())
    return n, float(correct[selected].mean()) if n else np.nan


def derive_metric_thresholds(
    values: np.ndarray,
    correct: np.ndarray,
    higher_is_better: bool,
) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    correct = np.asarray(correct, dtype=bool)
    if not np.isfinite(values).all() or len(values) != len(correct):
        raise ValueError("Calibration metrics must be finite and aligned with correctness labels")
    reliability = values if higher_is_better else -values
    baseline = float(correct.mean())
    high_target = min(0.95, baseline + 0.20)
    low_target = max(0.25, baseline - 0.20)
    min_n = max(5, int(math.ceil(len(values) * 0.02)))

    high_candidates: list[tuple[float, float, int, float]] = []
    for quantile in np.arange(0.60, 0.951, 0.05):
        cutoff = float(np.quantile(reliability, quantile))
        n, accuracy = _tail_accuracy(reliability, correct, cutoff, favorable=True)
        if n >= min_n:
            high_candidates.append((float(quantile), cutoff, n, accuracy))
    qualifying_high = [row for row in high_candidates if row[3] >= high_target]
    if qualifying_high:
        high_quantile, high_cutoff, _, _ = qualifying_high[0]
    else:
        high_quantile, high_cutoff, _, _ = sorted(
            high_candidates, key=lambda row: (-row[3], row[0])
        )[0]

    low_candidates: list[tuple[float, float, int, float]] = []
    for quantile in np.arange(0.40, 0.049, -0.05):
        cutoff = float(np.quantile(reliability, quantile))
        n, accuracy = _tail_accuracy(reliability, correct, cutoff, favorable=False)
        if n >= min_n:
            low_candidates.append((float(quantile), cutoff, n, accuracy))
    qualifying_low = [row for row in low_candidates if row[3] <= low_target]
    if qualifying_low:
        low_quantile, low_cutoff, _, _ = qualifying_low[0]
    else:
        low_quantile, low_cutoff, _, _ = sorted(
            low_candidates, key=lambda row: (row[3], -row[0])
        )[0]

    if low_cutoff >= high_cutoff:
        low_quantile = 1 / 3
        high_quantile = 2 / 3
        low_cutoff = float(np.quantile(reliability, low_quantile))
        high_cutoff = float(np.quantile(reliability, high_quantile))
    median_cutoff = float(np.quantile(reliability, 0.50))
    return {
        "higher_is_better": bool(higher_is_better),
        "low_reliability_cutoff": low_cutoff,
        "moderate_reliability_anchor": median_cutoff,
        "high_reliability_cutoff": high_cutoff,
        "low_quantile": low_quantile,
        "high_quantile": high_quantile,
        "baseline_accuracy": baseline,
        "high_target_accuracy": high_target,
        "low_target_accuracy": low_target,
    }


def assign_metric_tier(values: np.ndarray, threshold: dict[str, float]) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    reliability = values if threshold["higher_is_better"] else -values
    result = np.full(len(values), "moderate_confidence", dtype=object)
    result[reliability <= threshold["low_reliability_cutoff"]] = "low_confidence"
    result[reliability >= threshold["high_reliability_cutoff"]] = "high_confidence"
    return result.astype(str)


def external_threshold_value(reliability_value: float, higher_is_better: bool) -> float:
    return float(reliability_value if higher_is_better else -reliability_value)


def calibrate_confidence(
    evaluation: pd.DataFrame,
    source_predictions: pd.DataFrame,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray]:
    source = (
        source_predictions.set_index("sample_barcode")
        .loc[evaluation["sample_barcode"]]
        .reset_index()
    )
    calibrated = source.copy()
    calibration_mask, holdout_mask = deterministic_calibration_split(source, seed)
    threshold_rows: list[dict[str, Any]] = []
    thresholds: dict[tuple[str, str], dict[str, float]] = {}

    for level, spec in LEVEL_SPECS.items():
        correct = evaluation[spec["correct_column"]].to_numpy(dtype=bool)
        for metric_name, (column_suffix, higher_is_better) in CALIBRATION_METRICS.items():
            column = f"{level}_{column_suffix}"
            threshold = derive_metric_thresholds(
                evaluation.loc[calibration_mask, column].to_numpy(dtype=float),
                correct[calibration_mask],
                higher_is_better,
            )
            thresholds[(level, metric_name)] = threshold
            holdout_tier = assign_metric_tier(
                evaluation.loc[holdout_mask, column].to_numpy(dtype=float), threshold
            )
            holdout_correct = correct[holdout_mask]
            accuracy_by_tier = {
                tier: (
                    float(holdout_correct[holdout_tier == tier].mean())
                    if (holdout_tier == tier).any()
                    else np.nan
                )
                for tier in ["low_confidence", "moderate_confidence", "high_confidence"]
            }
            threshold_rows.append(
                {
                    "prediction_level": level,
                    "metric": metric_name,
                    "threshold_low_confidence": external_threshold_value(
                        threshold["low_reliability_cutoff"], higher_is_better
                    ),
                    "threshold_moderate_confidence": external_threshold_value(
                        threshold["moderate_reliability_anchor"], higher_is_better
                    ),
                    "threshold_high_confidence": external_threshold_value(
                        threshold["high_reliability_cutoff"], higher_is_better
                    ),
                    "observed_accuracy_low_confidence": accuracy_by_tier["low_confidence"],
                    "observed_accuracy_moderate_confidence": accuracy_by_tier[
                        "moderate_confidence"
                    ],
                    "observed_accuracy_high_confidence": accuracy_by_tier["high_confidence"],
                    "recommended_use": (
                        "Equal-vote input to the combined validation-informed tier; high/reliable requires "
                        "a combined metric score of at least 6 of 8."
                    ),
                    "notes": (
                        f"Derived on deterministic internal calibration half; evaluated on held-out internal half; "
                        f"higher_is_better={higher_is_better}; low_tail_quantile={threshold['low_quantile']:.2f}; "
                        f"high_tail_quantile={threshold['high_quantile']:.2f}; empirical reliability only."
                    ),
                }
            )

    for level in LEVEL_SPECS:
        votes = np.zeros((len(evaluation), len(CALIBRATION_METRICS)), dtype=int)
        for metric_index, (metric_name, (column_suffix, _)) in enumerate(
            CALIBRATION_METRICS.items()
        ):
            metric_tier = assign_metric_tier(
                evaluation[f"{level}_{column_suffix}"].to_numpy(dtype=float),
                thresholds[(level, metric_name)],
            )
            votes[:, metric_index] = pd.Series(metric_tier).map(
                {"low_confidence": 0, "moderate_confidence": 1, "high_confidence": 2}
            )
        score = votes.sum(axis=1)
        combined = np.full(len(score), "moderate_confidence", dtype=object)
        combined[score <= 2] = "low_confidence"
        combined[score >= 6] = "high_confidence"
        calibrated[f"calibrated_{level}_confidence"] = combined.astype(str)
        calibrated[f"calibrated_{level}_reliable_binary"] = combined == "high_confidence"

    calibrated["calibration_notes"] = (
        "Validation-informed empirical reliability tier from four transparent metrics; thresholds were derived "
        "on a deterministic internal TCGA calibration half; reliable means high_confidence; similarity weights "
        "are not probabilities and external clinical validation remains required."
    )
    thresholds_frame = pd.DataFrame(threshold_rows)
    summary = calibration_performance_summary(calibrated, holdout_mask)
    return thresholds_frame, calibrated, summary, holdout_mask


def calibration_performance_summary(
    calibrated: pd.DataFrame, holdout_mask: np.ndarray
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    holdout = calibrated.loc[holdout_mask].copy()
    tier_order = ["low_confidence", "moderate_confidence", "high_confidence"]
    for level, spec in LEVEL_SPECS.items():
        tier_column = f"calibrated_{level}_confidence"
        reliable_column = f"calibrated_{level}_reliable_binary"
        correct = holdout[spec["correct_column"]].astype(bool)
        for tier in tier_order:
            subset = holdout[tier_column] == tier
            rows.append(
                {
                    "prediction_level": level,
                    "summary_type": "calibrated_confidence_tier",
                    "stratum": tier,
                    "n_samples": int(subset.sum()),
                    "fraction_samples": float(subset.mean()),
                    "observed_accuracy": float(correct[subset].mean()) if subset.any() else np.nan,
                    "sensitivity_for_correct_call": np.nan,
                    "specificity_for_incorrect_call": np.nan,
                    "positive_predictive_value": (
                        float(correct[subset].mean()) if tier == "high_confidence" and subset.any() else np.nan
                    ),
                    "negative_predictive_value": np.nan,
                    "comparison_reference": "held_out_internal_calibration_evaluation_half",
                    "notes": "Observed tier accuracy is empirical internal reliability, not a calibrated probability.",
                }
            )
        reliable = holdout[reliable_column].astype(bool)
        tp = int((reliable & correct).sum())
        fp = int((reliable & ~correct).sum())
        tn = int((~reliable & ~correct).sum())
        fn = int((~reliable & correct).sum())
        rows.append(
            {
                "prediction_level": level,
                "summary_type": "reliable_binary_performance",
                "stratum": "high_confidence_vs_not_high_confidence",
                "n_samples": len(holdout),
                "fraction_samples": float(reliable.mean()),
                "observed_accuracy": float(correct.mean()),
                "sensitivity_for_correct_call": tp / (tp + fn) if tp + fn else np.nan,
                "specificity_for_incorrect_call": tn / (tn + fp) if tn + fp else np.nan,
                "positive_predictive_value": tp / (tp + fp) if tp + fp else np.nan,
                "negative_predictive_value": tn / (tn + fn) if tn + fn else np.nan,
                "comparison_reference": "correct_top1_call_as_outcome",
                "notes": "Reliable is defined as high_confidence; sensitivity/specificity-style metrics are descriptive.",
            }
        )
        for ambiguity in ["low_ambiguity", "moderate_ambiguity", "high_ambiguity"]:
            subset = holdout["ambiguity_class"].astype(str) == ambiguity
            rows.append(
                {
                    "prediction_level": level,
                    "summary_type": "previous_heuristic_ambiguity",
                    "stratum": ambiguity,
                    "n_samples": int(subset.sum()),
                    "fraction_samples": float(subset.mean()),
                    "observed_accuracy": float(correct[subset].mean()) if subset.any() else np.nan,
                    "sensitivity_for_correct_call": np.nan,
                    "specificity_for_incorrect_call": np.nan,
                    "positive_predictive_value": np.nan,
                    "negative_predictive_value": np.nan,
                    "comparison_reference": "previous_known_primary_ambiguity_class",
                    "notes": "Previous heuristic ambiguity retained for comparison only.",
                }
            )
    return pd.DataFrame(rows)


def common_confusions(
    group: pd.DataFrame, true_column: str, predicted_column: str, limit: int = 3
) -> str:
    wrong = group[group[true_column].astype(str) != group[predicted_column].astype(str)]
    if wrong.empty:
        return "none"
    counts = wrong[predicted_column].astype(str).value_counts().head(limit)
    return "; ".join(f"{label}:{int(count)}" for label, count in counts.items())


def ablation_summary_row(
    feature_set: str,
    features: list[str],
    evaluation: pd.DataFrame,
    top_k: int,
    metric: str,
) -> dict[str, Any]:
    ambiguity_counts = evaluation["heuristic_ambiguity_class"].value_counts(normalize=True)
    return {
        "feature_set": feature_set,
        "n_features": len(features),
        "n_samples_validated": len(evaluation),
        "top1_project_accuracy": evaluation["correct_top_project"].mean(),
        "top3_project_accuracy": evaluation["true_project_in_top3"].mean(),
        "top5_project_accuracy": evaluation["true_project_in_top5"].mean(),
        "top1_lineage_accuracy": evaluation["correct_top_lineage"].mean(),
        "top3_lineage_accuracy": evaluation["true_lineage_in_top3"].mean(),
        "top1_major_group_accuracy": evaluation["correct_top_major_group"].mean(),
        "median_project_margin": evaluation["project_margin"].median(),
        "median_lineage_margin": evaluation["lineage_margin"].median(),
        "high_ambiguity_fraction": ambiguity_counts.get("high_ambiguity", 0.0),
        "low_ambiguity_fraction": ambiguity_counts.get("low_ambiguity", 0.0),
        "notes": (
            f"Internal TCGA pseudo-unknown ablation; self excluded; top_k={top_k}; metric={metric}; "
            "locked scaling retained; not external clinical validation."
        ),
    }


def project_ablation_rows(feature_set: str, evaluation: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for project, group in evaluation.groupby("true_project_code", sort=True):
        rows.append(
            {
                "feature_set": feature_set,
                "project_code": project,
                "n_samples": len(group),
                "top1_project_accuracy": group["correct_top_project"].mean(),
                "top3_project_accuracy": group["true_project_in_top3"].mean(),
                "top1_lineage_accuracy": group["correct_top_lineage"].mean(),
                "median_margin": group["project_margin"].median(),
                "common_confusions": common_confusions(
                    group, "true_project_code", "predicted_top_project"
                ),
                "notes": "Internal ablation comparison in the locked scaled TCGA atlas.",
            }
        )
    return pd.DataFrame(rows)


def lineage_ablation_rows(feature_set: str, evaluation: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for lineage, group in evaluation.groupby("true_cup_lineage_group", sort=True):
        rows.append(
            {
                "feature_set": feature_set,
                "cup_lineage_group": lineage,
                "n_samples": len(group),
                "top1_lineage_accuracy": group["correct_top_lineage"].mean(),
                "top3_lineage_accuracy": group["true_lineage_in_top3"].mean(),
                "top1_project_accuracy": group["correct_top_project"].mean(),
                "median_margin": group["lineage_margin"].median(),
                "common_confusions": common_confusions(
                    group, "true_cup_lineage_group", "predicted_top_lineage"
                ),
                "notes": "Configured CUP lineage grouping; internal ablation comparison only.",
            }
        )
    return pd.DataFrame(rows)


def confusion_long(
    feature_set: str,
    evaluation: pd.DataFrame,
    true_column: str,
    predicted_column: str,
) -> pd.DataFrame:
    return (
        evaluation.groupby([true_column, predicted_column], sort=True)
        .size()
        .rename("count")
        .reset_index()
        .rename(columns={true_column: "true_label", predicted_column: "predicted_label"})
        .assign(feature_set=feature_set)[["feature_set", "true_label", "predicted_label", "count"]]
    )


def run_feature_ablation(
    matrix: np.ndarray,
    feature_order: list[str],
    feature_groups: dict[str, list[str]],
    metadata: pd.DataFrame,
    query_indices: np.ndarray,
    feature_coverage: np.ndarray,
    lineage_order: list[str],
    lineage_major_groups: dict[str, set[str]],
    top_k: int,
    metric: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    feature_index = {feature: index for index, feature in enumerate(feature_order)}
    cache: dict[tuple[str, ...], pd.DataFrame] = {}
    evaluations: dict[str, pd.DataFrame] = {}
    summary_rows: list[dict[str, Any]] = []
    project_rows: list[pd.DataFrame] = []
    lineage_rows: list[pd.DataFrame] = []
    project_confusions: list[pd.DataFrame] = []
    lineage_confusions: list[pd.DataFrame] = []

    for group_number, feature_set in enumerate(FEATURE_SET_ORDER, start=1):
        features = feature_groups[feature_set]
        key = tuple(features)
        logging.info(
            "Evaluating feature set %s (%d/%d; %d features)",
            feature_set,
            group_number,
            len(FEATURE_SET_ORDER),
            len(features),
        )
        if key not in cache:
            columns = [feature_index[feature] for feature in features]
            cache[key] = evaluate_feature_matrix(
                matrix[:, columns],
                metadata,
                query_indices,
                feature_coverage,
                lineage_order,
                lineage_major_groups,
                top_k,
                metric,
            )
        evaluation = cache[key].copy()
        evaluations[feature_set] = evaluation
        summary_rows.append(ablation_summary_row(feature_set, features, evaluation, top_k, metric))
        project_rows.append(project_ablation_rows(feature_set, evaluation))
        lineage_rows.append(lineage_ablation_rows(feature_set, evaluation))
        project_confusions.append(
            confusion_long(
                feature_set, evaluation, "true_project_code", "predicted_top_project"
            )
        )
        lineage_confusions.append(
            confusion_long(
                feature_set,
                evaluation,
                "true_cup_lineage_group",
                "predicted_top_lineage",
            )
        )
    return (
        pd.DataFrame(summary_rows),
        pd.concat(project_rows, ignore_index=True),
        pd.concat(lineage_rows, ignore_index=True),
        pd.concat(project_confusions, ignore_index=True),
        pd.concat(lineage_confusions, ignore_index=True),
        evaluations,
    )


def figure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.grid": False,
        }
    )


def plot_calibration_accuracy(summary: pd.DataFrame, path: Path) -> None:
    data = summary[summary["summary_type"] == "calibrated_confidence_tier"].copy()
    tiers = ["low_confidence", "moderate_confidence", "high_confidence"]
    levels = ["project", "lineage", "major_group"]
    pivot = data.pivot(index="stratum", columns="prediction_level", values="observed_accuracy").reindex(
        tiers
    )
    x = np.arange(len(tiers))
    width = 0.25
    colors = {"project": "#B44D5E", "lineage": "#3E7898", "major_group": "#4F7F5C"}
    figure_style()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for offset, level in enumerate(levels):
        ax.bar(
            x + (offset - 1) * width,
            pivot[level],
            width,
            label=level.replace("_", " ").title(),
            color=colors[level],
        )
    ax.set_xticks(x, [tier.replace("_", " ") for tier in tiers])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Observed top-1 accuracy")
    ax.set_title("Held-out internal accuracy by validation-informed confidence tier")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def binned_accuracy(
    evaluation: pd.DataFrame,
    holdout_mask: np.ndarray,
    value_suffix: str,
    n_bins: int = 10,
) -> pd.DataFrame:
    rows = []
    for level, spec in LEVEL_SPECS.items():
        frame = evaluation.loc[holdout_mask, [f"{level}_{value_suffix}", spec["correct_column"]]].copy()
        frame.columns = ["value", "correct"]
        try:
            frame["bin"] = pd.qcut(frame["value"], q=n_bins, duplicates="drop")
        except ValueError:
            frame["bin"] = "all"
        for _, group in frame.groupby("bin", observed=True, sort=True):
            rows.append(
                {
                    "prediction_level": level,
                    "median_value": group["value"].median(),
                    "observed_accuracy": group["correct"].mean(),
                    "n_samples": len(group),
                }
            )
    return pd.DataFrame(rows)


def plot_metric_accuracy(
    data: pd.DataFrame, path: Path, title: str, xlabel: str
) -> None:
    colors = {"project": "#B44D5E", "lineage": "#3E7898", "major_group": "#4F7F5C"}
    figure_style()
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for level, group in data.groupby("prediction_level", sort=False):
        group = group.sort_values("median_value")
        ax.plot(
            group["median_value"],
            group["observed_accuracy"],
            marker="o",
            linewidth=1.8,
            markersize=4,
            label=level.replace("_", " ").title(),
            color=colors[level],
        )
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Observed top-1 accuracy")
    ax.set_title(title)
    ax.grid(color="#dddddd", linewidth=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_ablation_overall(summary: pd.DataFrame, path: Path) -> None:
    data = summary.copy()
    data["_order"] = data["feature_set"].map(
        {feature_set: index for index, feature_set in enumerate(FEATURE_SET_ORDER)}
    )
    data = data.sort_values("_order", ascending=False)
    y = np.arange(len(data))
    height = 0.24
    figure_style()
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.barh(y - height, data["top1_project_accuracy"], height, label="Project", color="#B44D5E")
    ax.barh(y, data["top1_lineage_accuracy"], height, label="Lineage", color="#3E7898")
    ax.barh(y + height, data["top1_major_group_accuracy"], height, label="Major group", color="#4F7F5C")
    ax.set_yticks(y, data["feature_set"])
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Top-1 accuracy")
    ax.set_title("Feature-ablation performance in the locked TCGA atlas")
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_ablation_lineage(summary: pd.DataFrame, path: Path) -> None:
    data = summary.copy()
    data["_order"] = data["feature_set"].map(
        {feature_set: index for index, feature_set in enumerate(FEATURE_SET_ORDER)}
    )
    data = data.sort_values("_order", ascending=False)
    y = np.arange(len(data))
    figure_style()
    fig, ax = plt.subplots(figsize=(11, 7.5))
    ax.barh(y - 0.17, data["top1_lineage_accuracy"], 0.34, label="Top 1", color="#3E7898")
    ax.barh(y + 0.17, data["top3_lineage_accuracy"], 0.34, label="Top 3", color="#4F7F5C")
    ax.set_yticks(y, data["feature_set"])
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Lineage accuracy")
    ax.set_title("Configured-lineage recovery by feature set")
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_project_accuracy_heatmap(project_performance: pd.DataFrame, path: Path) -> None:
    pivot = project_performance.pivot(
        index="feature_set", columns="project_code", values="top1_project_accuracy"
    ).reindex(FEATURE_SET_ORDER)
    figure_style()
    fig, ax = plt.subplots(figsize=(15, 7.5))
    image = ax.imshow(pivot.to_numpy(dtype=float), cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(pivot.columns)), labels=pivot.columns, rotation=90)
    ax.set_yticks(np.arange(len(pivot.index)), labels=pivot.index)
    ax.set_xlabel("TCGA project")
    ax.set_ylabel("Feature set")
    ax.set_title("Top-1 project accuracy across feature-ablation models")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Top-1 project accuracy")
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_best_confusion_heatmap(
    project_confusions: pd.DataFrame,
    best_feature_set: str,
    project_order: list[str],
    path: Path,
) -> None:
    subset = project_confusions[project_confusions["feature_set"] == best_feature_set]
    matrix = (
        subset.pivot(index="true_label", columns="predicted_label", values="count")
        .reindex(index=project_order, columns=project_order, fill_value=0)
        .fillna(0)
        .to_numpy(dtype=float)
    )
    row_totals = matrix.sum(axis=1, keepdims=True)
    normalized = np.divide(matrix, row_totals, out=np.zeros_like(matrix), where=row_totals > 0)
    figure_style()
    fig, ax = plt.subplots(figsize=(14, 11.5))
    image = ax.imshow(normalized, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(project_order)), labels=project_order, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(project_order)), labels=project_order, fontsize=7)
    ax.set_xlabel("Predicted project")
    ax.set_ylabel("Known project")
    ax.set_title(f"Best ablation model confusion matrix: {best_feature_set}")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    colorbar.set_label("Row-normalized fraction")
    fig.text(
        0.01,
        0.01,
        "Internal TCGA pseudo-unknown validation; not external clinical performance.",
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout(rect=[0, 0.025, 1, 1])
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def markdown_table(
    data: pd.DataFrame, columns: list[str], digits: int = 3, max_rows: int | None = None
) -> str:
    view = data[columns].head(max_rows).copy() if max_rows is not None else data[columns].copy()
    for column in view.select_dtypes(include=["float"]).columns:
        view[column] = view[column].map(
            lambda value: f"{value:.{digits}f}" if pd.notna(value) else "NA"
        )
    headers = [column.replace("_", " ").title() for column in columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in view.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def expression_dependency_tables(
    project_performance: pd.DataFrame, lineage_performance: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    project_pivot = project_performance.pivot(
        index="project_code", columns="feature_set", values="top1_project_accuracy"
    )
    project_dependency = pd.DataFrame(
        {
            "project_code": project_pivot.index,
            "combined_accuracy": project_pivot["WES_plus_WTS_combined"].to_numpy(),
            "WES_like_accuracy": project_pivot["WES_like_features"].to_numpy(),
            "WTS_like_accuracy": project_pivot["WTS_like_features"].to_numpy(),
        }
    )
    project_dependency["combined_minus_WES_like"] = (
        project_dependency["combined_accuracy"] - project_dependency["WES_like_accuracy"]
    )
    project_dependency["WTS_like_minus_WES_like"] = (
        project_dependency["WTS_like_accuracy"] - project_dependency["WES_like_accuracy"]
    )
    project_dependency = project_dependency.sort_values(
        ["combined_minus_WES_like", "project_code"], ascending=[False, True]
    ).reset_index(drop=True)

    lineage_pivot = lineage_performance.pivot(
        index="cup_lineage_group", columns="feature_set", values="top1_lineage_accuracy"
    )
    lineage_dependency = pd.DataFrame(
        {
            "cup_lineage_group": lineage_pivot.index,
            "combined_accuracy": lineage_pivot["WES_plus_WTS_combined"].to_numpy(),
            "WES_like_accuracy": lineage_pivot["WES_like_features"].to_numpy(),
            "WTS_like_accuracy": lineage_pivot["WTS_like_features"].to_numpy(),
        }
    )
    lineage_dependency["combined_minus_WES_like"] = (
        lineage_dependency["combined_accuracy"] - lineage_dependency["WES_like_accuracy"]
    )
    lineage_dependency["WTS_like_minus_WES_like"] = (
        lineage_dependency["WTS_like_accuracy"] - lineage_dependency["WES_like_accuracy"]
    )
    lineage_dependency = lineage_dependency.sort_values(
        ["combined_minus_WES_like", "cup_lineage_group"], ascending=[False, True]
    ).reset_index(drop=True)
    return project_dependency, lineage_dependency


def robust_and_unstable_groups(
    project_performance: pd.DataFrame, lineage_performance: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    project_robustness = (
        project_performance.groupby("project_code", as_index=False)
        .agg(
            mean_top1_accuracy=("top1_project_accuracy", "mean"),
            minimum_top1_accuracy=("top1_project_accuracy", "min"),
            maximum_top1_accuracy=("top1_project_accuracy", "max"),
        )
        .sort_values(["mean_top1_accuracy", "project_code"], ascending=[False, True])
    )
    lineage_robustness = (
        lineage_performance.groupby("cup_lineage_group", as_index=False)
        .agg(
            mean_top1_accuracy=("top1_lineage_accuracy", "mean"),
            minimum_top1_accuracy=("top1_lineage_accuracy", "min"),
            maximum_top1_accuracy=("top1_lineage_accuracy", "max"),
        )
        .sort_values(["mean_top1_accuracy", "cup_lineage_group"], ascending=[False, True])
    )
    return (
        project_robustness.head(8),
        project_robustness.sort_values(["mean_top1_accuracy", "project_code"]).head(8),
        lineage_robustness.head(8),
        lineage_robustness.sort_values(["mean_top1_accuracy", "cup_lineage_group"]).head(8),
    )


def top_confusions_for_model(
    project_confusions: pd.DataFrame, feature_set: str, limit: int = 12
) -> pd.DataFrame:
    subset = project_confusions[
        (project_confusions["feature_set"] == feature_set)
        & (project_confusions["true_label"] != project_confusions["predicted_label"])
    ].copy()
    return subset.sort_values(
        ["count", "true_label", "predicted_label"], ascending=[False, True, True]
    ).head(limit)


def build_report(
    thresholds: pd.DataFrame,
    calibration_summary: pd.DataFrame,
    ablation_summary: pd.DataFrame,
    project_performance: pd.DataFrame,
    lineage_performance: pd.DataFrame,
    project_confusions: pd.DataFrame,
    unclassified: list[str],
    calibration_n: int,
    holdout_n: int,
    top_k: int,
    metric: str,
) -> str:
    best = ablation_summary.sort_values(
        ["top1_project_accuracy", "top1_lineage_accuracy", "n_features"],
        ascending=[False, False, True],
    ).iloc[0]
    best_feature_set = str(best["feature_set"])
    tier_rows = calibration_summary[
        calibration_summary["summary_type"] == "calibrated_confidence_tier"
    ]
    modality_rows = ablation_summary[
        ablation_summary["feature_set"].isin(
            ["WES_like_features", "WTS_like_features", "WES_plus_WTS_combined"]
        )
    ]
    project_dependency, lineage_dependency = expression_dependency_tables(
        project_performance, lineage_performance
    )
    robust_projects, unstable_projects, robust_lineages, unstable_lineages = robust_and_unstable_groups(
        project_performance, lineage_performance
    )
    confusions = top_confusions_for_model(project_confusions, best_feature_set)
    threshold_view = thresholds.copy()
    threshold_view["metric"] = threshold_view["metric"].str.replace("_", " ")
    feature_rank = ablation_summary.sort_values(
        ["top1_project_accuracy", "top1_lineage_accuracy"], ascending=[False, False]
    )

    return f"""# Calibration and Feature-Ablation Validation

Generated by `scripts/24_calibrate_projection_confidence_and_feature_ablation.py`.

## Purpose

This analysis uses the locked 97-feature TCGA atlas and known-primary pseudo-unknown results to derive validation-informed confidence tiers and to measure projection performance after restricting the molecular feature space. Similarity weights are not probabilities, and the confidence tiers represent empirical internal reliability rather than clinical diagnostic certainty.

## Calibration Approach

- The complete query cohort was split deterministically within TCGA project into {calibration_n} calibration samples and {holdout_n} internal evaluation samples.
- Thresholds were derived only on the calibration partition and tier accuracy was measured on the held-out internal partition.
- Four transparent metrics were assessed separately for project, configured lineage, and broad major-group predictions: top-1/top-2 margin, entropy, top-1 similarity-weight share, and top-label neighbor concentration.
- Favorable and unfavorable metric tails were selected by observed calibration accuracy, subject to minimum tail size; the four metric tiers then contributed equal 0/1/2 votes.
- A combined score of at least 6 of 8 defines `high_confidence`; at most 2 defines `low_confidence`; intermediate scores define `moderate_confidence`. The reliable binary flag means `high_confidence` only.

This is a deterministic train/evaluation split inside TCGA, not independent external calibration. It reduces direct threshold-fit reporting bias but does not remove shared-cohort and preprocessing leakage.

## Recommended Confidence Thresholds

For margin, top-1 score, and neighbor concentration, larger values are more favorable. For entropy, smaller values are more favorable; its low-confidence threshold is therefore numerically higher than its high-confidence threshold.

{markdown_table(threshold_view, ["prediction_level", "metric", "threshold_low_confidence", "threshold_moderate_confidence", "threshold_high_confidence", "observed_accuracy_low_confidence", "observed_accuracy_moderate_confidence", "observed_accuracy_high_confidence"])}

These thresholds are recommended as validation-informed reporting gates for research projections. They should replace unsupported verbal certainty, but they must not be represented as calibrated probabilities or externally validated diagnostic cutoffs.

## Accuracy by Calibrated Confidence Tier

{markdown_table(tier_rows, ["prediction_level", "stratum", "n_samples", "fraction_samples", "observed_accuracy"])}

`calibration_performance_summary.tsv` also records sensitivity/specificity-style summaries for the high-confidence reliable flag and compares the new tiers with the prior heuristic ambiguity classes.

## Feature-Ablation Design

For each of 13 prespecified feature sets, the analysis retained the locked TCGA standardization, recomputed exact leave-one-out {metric} neighbors, excluded the query sample, summarized the top {top_k} neighbors, and measured the same project, lineage, and broad-group recovery endpoints. The WES-plus-WTS set and all-feature set are intentionally equivalent contract checks. No feature subset was rescaled or tuned separately.

Unclassified locked features: {len(unclassified)}{f" ({'; '.join(unclassified)})" if unclassified else ""}.

## WES-Like Versus WTS-Like Performance

{markdown_table(modality_rows, ["feature_set", "n_features", "top1_project_accuracy", "top3_project_accuracy", "top5_project_accuracy", "top1_lineage_accuracy", "top3_lineage_accuracy", "top1_major_group_accuracy"])}

The synthetic clinical example previously supplied all WES-compatible features but only 35.484% of WTS-compatible features. The ablation comparison estimates the information available from entire locked modality subsets; it does not directly simulate arbitrary partial WTS missingness or cross-platform clinical assay effects.

## Which Feature Sets Perform Best

The highest top-1 project accuracy was observed for `{best_feature_set}` ({float(best['top1_project_accuracy']):.3f}; {int(best['n_features'])} features). Results across all feature sets are:

{markdown_table(feature_rank, ["feature_set", "n_features", "top1_project_accuracy", "top3_project_accuracy", "top5_project_accuracy", "top1_lineage_accuracy", "top1_major_group_accuracy"])}

## Which Projects/Lineages Are Robust

Projects with the highest mean top-1 recovery across feature sets:

{markdown_table(robust_projects, ["project_code", "mean_top1_accuracy", "minimum_top1_accuracy", "maximum_top1_accuracy"])}

Configured lineages with the highest mean recovery across feature sets:

{markdown_table(robust_lineages, ["cup_lineage_group", "mean_top1_accuracy", "minimum_top1_accuracy", "maximum_top1_accuracy"])}

Projects with the lowest mean recovery across feature sets:

{markdown_table(unstable_projects, ["project_code", "mean_top1_accuracy", "minimum_top1_accuracy", "maximum_top1_accuracy"])}

Configured lineages with the lowest mean recovery across feature sets:

{markdown_table(unstable_lineages, ["cup_lineage_group", "mean_top1_accuracy", "minimum_top1_accuracy", "maximum_top1_accuracy"])}

## Which Projects/Lineages Are Unstable or Confused

Largest gains from adding expression features to the WES-like set:

{markdown_table(project_dependency, ["project_code", "combined_accuracy", "WES_like_accuracy", "WTS_like_accuracy", "combined_minus_WES_like"], max_rows=10)}

{markdown_table(lineage_dependency, ["cup_lineage_group", "combined_accuracy", "WES_like_accuracy", "WTS_like_accuracy", "combined_minus_WES_like"], max_rows=10)}

These deltas identify expression dependence within this locked TCGA analysis. A negative delta means the combined unweighted Euclidean space performed worse than WES-like features alone; it does not mean expression is biologically irrelevant.

## Implications for CUP Reporting

- Report `high_confidence`, `moderate_confidence`, or `low_confidence` as validation-informed empirical reliability, never as a posterior probability.
- Treat only high-confidence calls as `reliable` under the binary research flag; retain explicit abstention and competing-lineage language for all other calls.
- Report modality coverage alongside the tier. A WES-only or WTS-incomplete query should be compared with the corresponding ablation benchmark rather than the all-feature benchmark alone.
- Preserve project, lineage, and broad-group tiers separately because their empirical accuracies differ.
- Continue to integrate morphology, immunophenotype, imaging, specimen context, and clinical findings.

## Limitations

- This remains internal TCGA validation, not external clinical validation.
- Calibration and evaluation partitions share TCGA cohort selection, assays, feature engineering, global median imputation, and standard scaling.
- Feature ablation uses pre-imputed, pre-scaled TCGA values and therefore tests information content under the locked atlas contract, not clinical assay failure or missing-feature mechanisms.
- Thresholds may be optimistic and are not guaranteed to transfer to metastatic, treated, low-purity, small-biopsy, cytology, or independently processed specimens.
- Euclidean distance gives every retained standardized feature equal weight; duplicated biological information can influence modality comparisons.
- Project and lineage labels are imperfect surrogates for tissue-of-origin truth, and small projects yield less stable estimates.

## Recommended Changes to CUP Workflow

1. Load `calibrated_confidence_thresholds.tsv` as a versioned research configuration and calculate all four metrics for project, lineage, and major-group summaries.
2. Replace the current purely heuristic confidence labels with the validation-informed equal-vote tiers while retaining the previous ambiguity class as a secondary warning signal.
3. Record feature modality coverage and select the nearest applicable WES-like, WTS-like, or combined internal benchmark in every report.
4. Validate these thresholds without refitting on an independent known-primary clinical WES/WTS cohort before clinical use.
5. Reassess thresholds after assay harmonization, top-k/distance sensitivity analysis, patient-level exclusion, and explicit partial-modality simulations.

Most common errors for the best-performing feature set:

{markdown_table(confusions, ["true_label", "predicted_label", "count"], digits=0)}
"""


def generate_figures(
    paths: AnalysisPaths,
    calibration_summary: pd.DataFrame,
    all_feature_evaluation: pd.DataFrame,
    holdout_mask: np.ndarray,
    ablation_summary: pd.DataFrame,
    project_performance: pd.DataFrame,
    project_confusions: pd.DataFrame,
) -> None:
    plot_calibration_accuracy(calibration_summary, paths.confidence_figure)
    plot_metric_accuracy(
        binned_accuracy(all_feature_evaluation, holdout_mask, "margin"),
        paths.margin_figure,
        "Observed accuracy across margin bins",
        "Median top-1/top-2 margin within bin",
    )
    plot_metric_accuracy(
        binned_accuracy(all_feature_evaluation, holdout_mask, "entropy"),
        paths.entropy_figure,
        "Observed accuracy across entropy bins",
        "Median normalized entropy within bin",
    )
    plot_ablation_overall(ablation_summary, paths.ablation_overall_figure)
    plot_ablation_lineage(ablation_summary, paths.ablation_lineage_figure)
    plot_project_accuracy_heatmap(project_performance, paths.ablation_project_heatmap)
    best_feature_set = str(
        ablation_summary.sort_values(
            ["top1_project_accuracy", "top1_lineage_accuracy", "n_features"],
            ascending=[False, False, True],
        ).iloc[0]["feature_set"]
    )
    project_order = sorted(project_performance["project_code"].astype(str).unique())
    plot_best_confusion_heatmap(
        project_confusions, best_feature_set, project_order, paths.best_confusion_heatmap
    )


def run_analysis(
    root: Path,
    top_k: int = 50,
    metric: str = "euclidean",
    max_samples: int | None = None,
    seed: int = 123,
    output_dir: Path | None = None,
) -> dict[str, Path]:
    paths = default_paths(root, output_dir=output_dir)
    inputs = load_inputs(paths)
    feature_groups, feature_group_table, unclassified = build_feature_groups(
        inputs["feature_order"], inputs["feature_records"]
    )
    query_indices = stratified_query_indices(inputs["metadata"], max_samples, seed)
    feature_coverage = pd.to_numeric(
        inputs["predictions"]["feature_coverage"], errors="coerce"
    ).to_numpy(dtype=float)
    if not np.isfinite(feature_coverage).all():
        raise ValueError("Known-primary feature coverage contains non-finite values")
    logging.info(
        "Starting calibration/ablation: queries=%d, reference=%d, feature_sets=%d, top_k=%d, metric=%s",
        len(query_indices),
        len(inputs["metadata"]),
        len(FEATURE_SET_ORDER),
        top_k,
        metric,
    )
    (
        ablation_summary,
        project_performance,
        lineage_performance,
        project_confusions,
        lineage_confusions,
        evaluations,
    ) = run_feature_ablation(
        inputs["matrix"],
        inputs["feature_order"],
        feature_groups,
        inputs["metadata"],
        query_indices,
        feature_coverage,
        inputs["lineage_order"],
        inputs["lineage_major_groups"],
        top_k,
        metric,
    )
    if top_k == 50 and metric == "euclidean":
        calibration_evaluation = evaluations["all_features"]
    else:
        logging.info(
            "Recomputing locked 50-neighbor Euclidean all-feature results for calibration contract"
        )
        calibration_evaluation = evaluate_feature_matrix(
            inputs["matrix"],
            inputs["metadata"],
            query_indices,
            feature_coverage,
            inputs["lineage_order"],
            inputs["lineage_major_groups"],
            50,
            "euclidean",
        )
    validate_all_feature_reproduction(calibration_evaluation, inputs["predictions"])
    source_subset = (
        inputs["predictions"]
        .set_index("sample_barcode")
        .loc[calibration_evaluation["sample_barcode"]]
        .reset_index()
    )
    thresholds, calibrated_predictions, calibration_summary, holdout_mask = calibrate_confidence(
        calibration_evaluation, source_subset, seed
    )

    write_tsv(thresholds, paths.thresholds)
    write_tsv(calibrated_predictions, paths.calibrated_predictions)
    write_tsv(calibration_summary, paths.calibration_summary)
    write_tsv(feature_group_table, paths.feature_groups)
    write_tsv(ablation_summary, paths.ablation_summary)
    write_tsv(project_performance, paths.ablation_project)
    write_tsv(lineage_performance, paths.ablation_lineage)
    write_tsv(project_confusions, paths.ablation_project_confusions)
    write_tsv(lineage_confusions, paths.ablation_lineage_confusions)
    generate_figures(
        paths,
        calibration_summary,
        calibration_evaluation,
        holdout_mask,
        ablation_summary,
        project_performance,
        project_confusions,
    )
    write_text(
        build_report(
            thresholds,
            calibration_summary,
            ablation_summary,
            project_performance,
            lineage_performance,
            project_confusions,
            unclassified,
            int((~holdout_mask).sum()),
            int(holdout_mask.sum()),
            top_k,
            metric,
        ),
        paths.report,
    )
    logging.info("Calibration and feature-ablation validation complete: %s", paths.output_dir)
    return {
        field: getattr(paths, field)
        for field in paths.__dataclass_fields__
        if field not in {
            "predictions",
            "overall_metrics",
            "ambiguity_performance",
            "scaled_matrix",
            "metadata",
            "feature_definitions",
            "scaler_parameters",
            "lineage_groups",
            "cancer_group_map",
            "output_dir",
        }
    }


def parse_max_samples(value: str) -> int | None:
    if value.strip().lower() == "all":
        return None
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--max-samples must be 'all' or a positive integer") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("--max-samples must be 'all' or a positive integer")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Derive validation-informed projection confidence tiers and rerun TCGA feature-ablation validation."
        )
    )
    parser.add_argument("--top-k", type=int, default=50, help="Neighbors summarized per query.")
    parser.add_argument(
        "--metric",
        choices=["euclidean", "cosine"],
        default="euclidean",
        help="Nearest-neighbor distance metric.",
    )
    parser.add_argument(
        "--max-samples",
        type=parse_max_samples,
        default=None,
        metavar="all|N",
        help="Validate all samples (default) or a deterministic project-stratified query subset.",
    )
    parser.add_argument("--seed", type=int, default=123, help="Sampling and calibration-split seed.")
    parser.add_argument("--root", type=Path, default=None, help="Project root; auto-detected by default.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional alternate output directory.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root or Path(__file__).resolve())
    run_analysis(
        root,
        top_k=args.top_k,
        metric=args.metric,
        max_samples=args.max_samples,
        seed=args.seed,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()

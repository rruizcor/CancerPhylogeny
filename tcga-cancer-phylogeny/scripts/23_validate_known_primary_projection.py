#!/usr/bin/env python

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
import pickle
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcga_cancer_phylogeny_matplotlib"))
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.neighbors import NearestNeighbors

from lib.common import configure_logging, ensure_dir, find_project_root, require_file


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

METADATA_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "disease_type",
    "primary_site",
    "major_cancer_group",
]

PREDICTION_COLUMNS = [
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
    "top2_project",
    "top2_lineage",
    "top1_top2_project_margin",
    "top1_top2_lineage_margin",
    "ambiguity_class",
    "correct_top_project",
    "correct_top_lineage",
    "correct_top_major_group",
    "true_project_in_top3",
    "true_project_in_top5",
    "true_lineage_in_top3",
    "neighbor_entropy",
    "feature_coverage",
    "notes",
]


@dataclass(frozen=True)
class ValidationPaths:
    raw_matrix: Path
    metadata: Path
    scaled_matrix: Path
    feature_definitions: Path
    scaler_parameters: Path
    nearest_neighbor_index: Path
    nearest_neighbor_metadata: Path
    lineage_groups: Path
    cancer_group_map: Path
    output_dir: Path
    cohort: Path
    predictions: Path
    project_performance: Path
    lineage_performance: Path
    major_group_performance: Path
    overall_metrics: Path
    project_confusion: Path
    lineage_confusion: Path
    major_group_confusion: Path
    ambiguity_performance: Path
    project_confusion_figure: Path
    lineage_confusion_figure: Path
    topk_project_figure: Path
    lineage_accuracy_figure: Path
    ambiguity_accuracy_figure: Path
    margin_figure: Path
    entropy_figure: Path
    report: Path


def default_paths(root: Path, output_dir: Path | None = None) -> ValidationPaths:
    out = output_dir or root / "results" / "validation" / "known_primary_projection"
    return ValidationPaths(
        raw_matrix=root / "results" / "reference_atlas" / "tcga_sample_reference_matrix.tsv.gz",
        metadata=root / "results" / "reference_atlas" / "tcga_sample_reference_metadata.tsv",
        scaled_matrix=root / "results" / "reference_atlas" / "tcga_sample_reference_scaled_matrix.tsv.gz",
        feature_definitions=root / "results" / "reference_atlas" / "tcga_reference_feature_definitions.yaml",
        scaler_parameters=root / "results" / "reference_atlas" / "tcga_reference_scaler_parameters.json",
        nearest_neighbor_index=root
        / "results"
        / "reference_atlas"
        / "tcga_reference_nearest_neighbor_index.pkl",
        nearest_neighbor_metadata=root
        / "results"
        / "reference_atlas"
        / "tcga_reference_nearest_neighbor_metadata.json",
        lineage_groups=root / "config" / "cup_lineage_groups.yaml",
        cancer_group_map=root / "config" / "cancer_group_map.csv",
        output_dir=out,
        cohort=out / "known_primary_validation_cohort.tsv",
        predictions=out / "known_primary_projection_predictions.tsv",
        project_performance=out / "known_primary_project_performance.tsv",
        lineage_performance=out / "known_primary_lineage_performance.tsv",
        major_group_performance=out / "known_primary_major_group_performance.tsv",
        overall_metrics=out / "known_primary_overall_metrics.tsv",
        project_confusion=out / "known_primary_project_confusion_matrix.tsv",
        lineage_confusion=out / "known_primary_lineage_confusion_matrix.tsv",
        major_group_confusion=out / "known_primary_major_group_confusion_matrix.tsv",
        ambiguity_performance=out / "known_primary_ambiguity_performance.tsv",
        project_confusion_figure=out / "known_primary_project_confusion_matrix.pdf",
        lineage_confusion_figure=out / "known_primary_lineage_confusion_matrix.pdf",
        topk_project_figure=out / "known_primary_topk_accuracy_by_project.pdf",
        lineage_accuracy_figure=out / "known_primary_accuracy_by_lineage.pdf",
        ambiguity_accuracy_figure=out / "known_primary_accuracy_by_ambiguity_class.pdf",
        margin_figure=out / "known_primary_margin_vs_correctness.pdf",
        entropy_figure=out / "known_primary_entropy_vs_correctness.pdf",
        report=out / "known_primary_validation_report.md",
    )


def write_tsv(data: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_text(text: str, path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")
    logging.info("Wrote %s", path)


def ordered_values_hash(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(str(value) for value in values).encode("utf-8")).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a YAML mapping in {path}")
    return value


def load_lineage_mapping(path: Path, project_codes: set[str]) -> tuple[dict[str, str], list[str]]:
    config = read_yaml(path)
    lineages = config.get("lineages")
    if not isinstance(lineages, list) or not lineages:
        raise ValueError("CUP lineage configuration must contain a non-empty lineages list")
    project_to_lineage: dict[str, str] = {}
    lineage_order: list[str] = []
    for record in lineages:
        name = str(record.get("name", "")).strip()
        projects = [str(value).strip() for value in record.get("projects", [])]
        if not name or not projects:
            raise ValueError("Each CUP lineage must define a name and at least one project")
        lineage_order.append(name)
        for project in projects:
            if project in project_to_lineage:
                raise ValueError(f"TCGA project {project} is assigned to multiple CUP lineages")
            project_to_lineage[project] = name
    missing = sorted(project_codes - set(project_to_lineage))
    extra = sorted(set(project_to_lineage) - project_codes)
    if missing or extra:
        raise ValueError(f"CUP lineage mapping mismatch; missing={missing}, extra={extra}")
    return project_to_lineage, lineage_order


def load_reference_inputs(paths: ValidationPaths) -> dict[str, Any]:
    for path, description in [
        (paths.raw_matrix, "raw TCGA reference matrix"),
        (paths.metadata, "TCGA reference metadata"),
        (paths.scaled_matrix, "scaled TCGA reference matrix"),
        (paths.feature_definitions, "TCGA feature definitions"),
        (paths.scaler_parameters, "TCGA scaler parameters"),
        (paths.nearest_neighbor_index, "TCGA nearest-neighbor index"),
        (paths.nearest_neighbor_metadata, "TCGA nearest-neighbor metadata"),
        (paths.lineage_groups, "CUP lineage groups"),
        (paths.cancer_group_map, "cancer group map"),
    ]:
        require_file(path, description)

    metadata = pd.read_csv(paths.metadata, sep="\t")
    raw = pd.read_csv(paths.raw_matrix, sep="\t")
    scaled = pd.read_csv(paths.scaled_matrix, sep="\t")
    definitions = read_yaml(paths.feature_definitions)
    scaler = read_json(paths.scaler_parameters)
    neighbor_metadata = read_json(paths.nearest_neighbor_metadata)
    cancer_map = pd.read_csv(paths.cancer_group_map)
    with paths.nearest_neighbor_index.open("rb") as handle:
        neighbor_package = pickle.load(handle)

    missing_metadata = sorted(set(METADATA_COLUMNS) - set(metadata.columns))
    if missing_metadata:
        raise ValueError(f"Reference metadata lack required columns: {missing_metadata}")
    if metadata[METADATA_COLUMNS].isna().any().any():
        raise ValueError("Reference metadata contain missing required values")
    if metadata["sample_barcode"].duplicated().any():
        raise ValueError("Reference metadata contain duplicate sample barcodes")
    if not metadata["sample_barcode"].astype(str).equals(raw["sample_barcode"].astype(str)):
        raise ValueError("Raw reference matrix row order does not match metadata")
    if not metadata["sample_barcode"].astype(str).equals(scaled["sample_barcode"].astype(str)):
        raise ValueError("Scaled reference matrix row order does not match metadata")

    feature_order = [str(value) for value in definitions.get("feature_order", [])]
    if not feature_order:
        raise ValueError("Feature definitions do not provide a locked feature_order")
    for label, order in [
        ("scaler", scaler.get("feature_order")),
        ("neighbor metadata", neighbor_metadata.get("feature_order")),
        ("neighbor pickle", neighbor_package.get("feature_order")),
        ("scaled matrix", scaled.columns[1:].tolist()),
    ]:
        if list(order or []) != feature_order:
            raise ValueError(f"Locked feature order mismatch in {label}")
    if not set(feature_order).issubset(raw.columns):
        raise ValueError("Raw reference matrix lacks one or more locked features")

    contract_hash = definitions.get("feature_contract_sha256")
    for label, artifact in [
        ("scaler", scaler),
        ("neighbor metadata", neighbor_metadata),
        ("neighbor pickle", neighbor_package),
    ]:
        if artifact.get("feature_contract_sha256") != contract_hash:
            raise ValueError(f"Feature contract hash mismatch in {label}")
    sample_hash = ordered_values_hash(metadata["sample_barcode"].astype(str).tolist())
    for label, artifact in [
        ("scaler", scaler),
        ("neighbor metadata", neighbor_metadata),
        ("neighbor pickle", neighbor_package),
    ]:
        if artifact.get("sample_order_sha256") != sample_hash:
            raise ValueError(f"Sample-order hash mismatch in {label}")

    numeric = scaled[feature_order].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("Scaled reference matrix contains non-finite feature values")
    project_codes = set(metadata["project_code"].astype(str))
    project_to_lineage, lineage_order = load_lineage_mapping(paths.lineage_groups, project_codes)

    required_map_columns = {"project_code", "project_id", "broad_group"}
    if not required_map_columns.issubset(cancer_map.columns):
        raise ValueError("Cancer group map lacks required project and broad-group columns")
    map_pairs = cancer_map.set_index("project_code")[["project_id", "broad_group"]].to_dict("index")
    for row in metadata[["project_code", "project_id", "major_cancer_group"]].drop_duplicates().itertuples(index=False):
        expected = map_pairs.get(str(row.project_code))
        if expected is None:
            raise ValueError(f"Project {row.project_code} is absent from cancer_group_map.csv")
        if str(row.project_id) != str(expected["project_id"]):
            raise ValueError(f"Project ID mismatch for {row.project_code}")
        if str(row.major_cancer_group) != str(expected["broad_group"]):
            raise ValueError(f"Major cancer-group mismatch for {row.project_code}")

    feature_coverage = raw[feature_order].notna().mean(axis=1).to_numpy(dtype=float)
    metadata = metadata.copy()
    metadata["cup_lineage_group"] = metadata["project_code"].map(project_to_lineage)
    return {
        "metadata": metadata,
        "matrix": numeric,
        "feature_coverage": feature_coverage,
        "feature_order": feature_order,
        "project_to_lineage": project_to_lineage,
        "lineage_order": lineage_order,
        "contract_hash": contract_hash,
        "available_metrics": neighbor_metadata.get("available_metrics", []),
        "default_metric": neighbor_metadata.get("default_metric"),
    }


def stratified_query_indices(
    metadata: pd.DataFrame,
    eligible_indices: np.ndarray,
    max_samples: int | None,
    seed: int,
) -> np.ndarray:
    if max_samples is None or max_samples >= len(eligible_indices):
        return np.sort(eligible_indices.astype(int))
    if max_samples <= 0:
        raise ValueError("max_samples must be positive or 'all'")
    eligible = metadata.iloc[eligible_indices].copy()
    eligible["_row_index"] = eligible_indices
    groups = list(eligible.groupby("project_code", sort=True))
    if max_samples < len(groups):
        raise ValueError(f"max_samples={max_samples} cannot retain all {len(groups)} TCGA projects")

    counts = eligible["project_code"].value_counts().sort_index()
    allocation = pd.Series(1, index=counts.index, dtype=int)
    remaining = max_samples - int(allocation.sum())
    capacities = counts - allocation
    if remaining > 0:
        proportions = capacities / capacities.sum() if capacities.sum() else capacities.astype(float)
        raw_extra = proportions * remaining
        extra = np.floor(raw_extra).astype(int).clip(upper=capacities)
        allocation += extra
        remaining -= int(extra.sum())
        if remaining > 0:
            fractional = (raw_extra - np.floor(raw_extra)).sort_values(ascending=False, kind="stable")
            for project in fractional.index:
                if remaining == 0:
                    break
                if allocation[project] < counts[project]:
                    allocation[project] += 1
                    remaining -= 1
        while remaining > 0:
            progressed = False
            for project in allocation.index:
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
    for project, group in groups:
        candidates = group["_row_index"].to_numpy(dtype=int)
        n_select = int(allocation[project])
        selected.extend(rng.choice(candidates, size=n_select, replace=False).tolist())
    return np.array(sorted(selected), dtype=int)


def build_validation_cohort(
    metadata: pd.DataFrame,
    max_samples: int | None,
    seed: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    cohort = metadata[METADATA_COLUMNS + ["cup_lineage_group"]].copy()
    valid = cohort["cup_lineage_group"].notna().to_numpy()
    eligible_indices = np.flatnonzero(valid)
    query_indices = stratified_query_indices(metadata, eligible_indices, max_samples, seed)
    selected = np.zeros(len(cohort), dtype=bool)
    selected[query_indices] = True
    subset = len(query_indices) < len(eligible_indices)
    cohort["validation_split"] = (
        f"stratified_subset_{len(query_indices)}_queries_against_full_atlas"
        if subset
        else "full_leave_one_out_against_full_atlas"
    )
    cohort["included_in_validation"] = selected
    cohort["exclusion_reason"] = ""
    cohort.loc[~valid, "exclusion_reason"] = "missing_cup_lineage_mapping"
    cohort.loc[valid & ~selected, "exclusion_reason"] = "not_selected_for_stratified_max_samples_subset"
    cohort = cohort[
        [
            "sample_barcode",
            "patient_barcode",
            "project_id",
            "project_code",
            "disease_type",
            "primary_site",
            "major_cancer_group",
            "cup_lineage_group",
            "validation_split",
            "included_in_validation",
            "exclusion_reason",
        ]
    ]
    return cohort, query_indices


def compute_neighbor_indices(
    matrix: np.ndarray,
    query_indices: np.ndarray,
    top_k: int,
    metric: str,
    patient_ids: np.ndarray | None = None,
    exclude_same_patient: bool = False,
    batch_size: int = 512,
) -> tuple[np.ndarray, np.ndarray]:
    if matrix.ndim != 2 or not np.isfinite(matrix).all():
        raise ValueError("Neighbor matrix must be a finite two-dimensional array")
    if top_k < 1 or top_k >= len(matrix):
        raise ValueError("top_k must be between 1 and n_reference_samples - 1")
    if metric not in {"euclidean", "cosine"}:
        raise ValueError("metric must be euclidean or cosine")
    if exclude_same_patient and patient_ids is None:
        raise ValueError("patient_ids are required when excluding same-patient samples")

    max_patient_size = 1
    if exclude_same_patient and patient_ids is not None:
        max_patient_size = int(pd.Series(patient_ids).value_counts().max())
    n_requested = min(len(matrix), top_k + max(1, max_patient_size) + 5)
    model = NearestNeighbors(metric=metric, algorithm="brute", n_jobs=-1, n_neighbors=n_requested)
    model.fit(matrix)
    neighbor_indices = np.empty((len(query_indices), top_k), dtype=int)
    neighbor_distances = np.empty((len(query_indices), top_k), dtype=float)

    for start in range(0, len(query_indices), batch_size):
        stop = min(start + batch_size, len(query_indices))
        batch_queries = query_indices[start:stop]
        distances, indices = model.kneighbors(matrix[batch_queries], return_distance=True)
        for offset, query_index in enumerate(batch_queries):
            keep = indices[offset] != query_index
            if exclude_same_patient and patient_ids is not None:
                keep &= patient_ids[indices[offset]] != patient_ids[query_index]
            filtered_indices = indices[offset][keep][:top_k]
            filtered_distances = distances[offset][keep][:top_k]
            if len(filtered_indices) < top_k:
                raise RuntimeError(
                    f"Only {len(filtered_indices)} eligible neighbors remained for query row {query_index}"
                )
            neighbor_indices[start + offset] = filtered_indices
            neighbor_distances[start + offset] = filtered_distances
        logging.info("Computed neighbors for %d/%d validation queries", stop, len(query_indices))
    return neighbor_indices, neighbor_distances


def aggregate_scores(labels: np.ndarray, distances: np.ndarray) -> pd.DataFrame:
    weights = 1.0 / (1.0 + distances.astype(float))
    frame = pd.DataFrame({"label": labels.astype(str), "weight": weights})
    summary = (
        frame.groupby("label", sort=False)
        .agg(n_neighbors=("label", "size"), weight=("weight", "sum"))
        .reset_index()
    )
    total = float(summary["weight"].sum())
    summary["score"] = summary["weight"] / total if total > 0 else 0.0
    return summary.sort_values(
        ["score", "n_neighbors", "label"], ascending=[False, False, True], kind="stable"
    ).reset_index(drop=True)


def normalized_entropy(scores: np.ndarray) -> float:
    positive = scores[np.isfinite(scores) & (scores > 0)]
    if len(positive) <= 1:
        return 0.0
    positive = positive / positive.sum()
    return float(-np.sum(positive * np.log(positive)) / math.log(len(positive)))


def classify_ambiguity(
    top1_score: float,
    margin: float,
    entropy: float,
    meaningful_projects: int,
    coverage: float,
    major_group_agreement: bool,
    config: dict[str, float] | None = None,
) -> str:
    cfg = config or AMBIGUITY_CONFIG
    if coverage < float(cfg["uninterpretable_coverage_below"]) or top1_score <= 0:
        return "uninterpretable_due_to_missing_features"
    if (
        top1_score >= float(cfg["low_ambiguity_top1_min"])
        and margin >= float(cfg["low_ambiguity_margin_min"])
        and entropy <= float(cfg["low_ambiguity_entropy_max"])
        and major_group_agreement
    ):
        return "low_ambiguity"
    if (
        top1_score < float(cfg["high_ambiguity_top1_below"])
        or margin < float(cfg["high_ambiguity_margin_below"])
        or entropy >= float(cfg["high_ambiguity_entropy_min"])
        or meaningful_projects >= int(cfg["high_ambiguity_meaningful_projects_min"])
    ):
        return "high_ambiguity"
    return "moderate_ambiguity"


def label_rank(summary: pd.DataFrame, label: str) -> int:
    matches = np.flatnonzero(summary["label"].to_numpy(dtype=str) == str(label))
    return int(matches[0] + 1) if len(matches) else int(len(summary) + 1)


def evaluate_queries(
    metadata: pd.DataFrame,
    query_indices: np.ndarray,
    neighbor_indices: np.ndarray,
    neighbor_distances: np.ndarray,
    feature_coverage: np.ndarray,
    top_k: int,
    metric: str,
    exclude_same_patient: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    projects = metadata["project_code"].astype(str).to_numpy()
    project_ids = metadata["project_id"].astype(str).to_numpy()
    major_groups = metadata["major_cancer_group"].astype(str).to_numpy()
    lineages = metadata["cup_lineage_group"].astype(str).to_numpy()
    barcodes = metadata["sample_barcode"].astype(str).to_numpy()
    lineage_major_groups = (
        metadata.groupby("cup_lineage_group")["major_cancer_group"]
        .agg(lambda values: set(values.astype(str)))
        .to_dict()
    )
    rows: list[dict[str, Any]] = []
    internal: list[dict[str, Any]] = []

    for row_number, query_index in enumerate(query_indices):
        indices = neighbor_indices[row_number]
        distances = neighbor_distances[row_number]
        project_scores = aggregate_scores(projects[indices], distances)
        lineage_scores = aggregate_scores(lineages[indices], distances)
        major_scores = aggregate_scores(major_groups[indices], distances)

        top_project = project_scores.iloc[0]
        top2_project = project_scores.iloc[1] if len(project_scores) > 1 else None
        top_lineage = lineage_scores.iloc[0]
        top2_lineage = lineage_scores.iloc[1] if len(lineage_scores) > 1 else None
        top_major = major_scores.iloc[0]
        top2_major = major_scores.iloc[1] if len(major_scores) > 1 else None

        project_margin = float(top_project["score"]) - (
            float(top2_project["score"]) if top2_project is not None else 0.0
        )
        lineage_margin = float(top_lineage["score"]) - (
            float(top2_lineage["score"]) if top2_lineage is not None else 0.0
        )
        major_margin = float(top_major["score"]) - (
            float(top2_major["score"]) if top2_major is not None else 0.0
        )
        entropy = normalized_entropy(lineage_scores["score"].to_numpy(dtype=float))
        meaningful_projects = int(
            (project_scores["score"] >= AMBIGUITY_CONFIG["meaningful_project_score_min"]).sum()
        )
        lineage_major_agreement = (
            str(top_major["label"]) in lineage_major_groups[str(top_lineage["label"])]
        )
        correct_major_group = str(top_major["label"]) == major_groups[query_index]
        ambiguity = classify_ambiguity(
            float(top_lineage["score"]),
            lineage_margin,
            entropy,
            meaningful_projects,
            float(feature_coverage[query_index]),
            lineage_major_agreement,
        )
        true_project_rank = label_rank(project_scores, projects[query_index])
        true_lineage_rank = label_rank(lineage_scores, lineages[query_index])
        true_major_rank = label_rank(major_scores, major_groups[query_index])

        rows.append(
            {
                "sample_barcode": barcodes[query_index],
                "true_project_code": projects[query_index],
                "true_project_id": project_ids[query_index],
                "true_major_cancer_group": major_groups[query_index],
                "true_cup_lineage_group": lineages[query_index],
                "predicted_top_project": str(top_project["label"]),
                "predicted_top_lineage": str(top_lineage["label"]),
                "predicted_top_major_group": str(top_major["label"]),
                "top_project_score": float(top_project["score"]),
                "top_lineage_score": float(top_lineage["score"]),
                "top_major_group_score": float(top_major["score"]),
                "top2_project": str(top2_project["label"]) if top2_project is not None else "NA",
                "top2_lineage": str(top2_lineage["label"]) if top2_lineage is not None else "NA",
                "top1_top2_project_margin": project_margin,
                "top1_top2_lineage_margin": lineage_margin,
                "ambiguity_class": ambiguity,
                "correct_top_project": str(top_project["label"]) == projects[query_index],
                "correct_top_lineage": str(top_lineage["label"]) == lineages[query_index],
                "correct_top_major_group": correct_major_group,
                "true_project_in_top3": true_project_rank <= 3,
                "true_project_in_top5": true_project_rank <= 5,
                "true_lineage_in_top3": true_lineage_rank <= 3,
                "neighbor_entropy": entropy,
                "feature_coverage": float(feature_coverage[query_index]),
                "notes": (
                    f"Internal TCGA pseudo-unknown validation; self excluded; metric={metric}; top_k={top_k}; "
                    f"same_patient_excluded={exclude_same_patient}; scores are descriptive weights, not probabilities."
                ),
            }
        )
        internal.append(
            {
                "sample_barcode": barcodes[query_index],
                "true_project_rank": true_project_rank,
                "true_lineage_rank": true_lineage_rank,
                "true_major_group_rank": true_major_rank,
                "top1_top2_major_group_margin": major_margin,
                "meaningful_projects": meaningful_projects,
                "neighbor_indices": ";".join(str(value) for value in indices),
            }
        )
    return pd.DataFrame(rows, columns=PREDICTION_COLUMNS), pd.DataFrame(internal)


def common_confusions(group: pd.DataFrame, true_column: str, predicted_column: str, limit: int = 3) -> str:
    wrong = group[group[true_column].astype(str) != group[predicted_column].astype(str)]
    if wrong.empty:
        return "none"
    counts = wrong[predicted_column].astype(str).value_counts().head(limit)
    return "; ".join(f"{label}:{int(count)}" for label, count in counts.items())


def project_performance(predictions: pd.DataFrame, internal: pd.DataFrame) -> pd.DataFrame:
    joined = predictions.merge(internal, on="sample_barcode", validate="one_to_one")
    rows = []
    for project, group in joined.groupby("true_project_code", sort=True):
        rows.append(
            {
                "project_code": project,
                "n_samples": len(group),
                "top1_project_accuracy": group["correct_top_project"].mean(),
                "top3_project_accuracy": group["true_project_in_top3"].mean(),
                "top5_project_accuracy": group["true_project_in_top5"].mean(),
                "median_true_project_rank": group["true_project_rank"].median(),
                "top1_lineage_accuracy": group["correct_top_lineage"].mean(),
                "median_top1_margin": group["top1_top2_project_margin"].median(),
                "median_entropy": group["neighbor_entropy"].median(),
                "high_ambiguity_fraction": (group["ambiguity_class"] == "high_ambiguity").mean(),
                "common_confusions": common_confusions(
                    group, "true_project_code", "predicted_top_project"
                ),
                "notes": "Internal TCGA pseudo-unknown validation; project estimates are not external clinical performance.",
            }
        )
    return pd.DataFrame(rows)


def lineage_performance(predictions: pd.DataFrame, internal: pd.DataFrame) -> pd.DataFrame:
    joined = predictions.merge(internal, on="sample_barcode", validate="one_to_one")
    rows = []
    for lineage, group in joined.groupby("true_cup_lineage_group", sort=True):
        rows.append(
            {
                "cup_lineage_group": lineage,
                "n_samples": len(group),
                "top1_lineage_accuracy": group["correct_top_lineage"].mean(),
                "top3_lineage_accuracy": group["true_lineage_in_top3"].mean(),
                "median_margin": group["top1_top2_lineage_margin"].median(),
                "median_entropy": group["neighbor_entropy"].median(),
                "high_ambiguity_fraction": (group["ambiguity_class"] == "high_ambiguity").mean(),
                "common_confusions": common_confusions(
                    group, "true_cup_lineage_group", "predicted_top_lineage"
                ),
                "notes": "Broad configured lineage grouping; not a clinically validated tissue-of-origin category.",
            }
        )
    return pd.DataFrame(rows)


def major_group_performance(predictions: pd.DataFrame, internal: pd.DataFrame) -> pd.DataFrame:
    joined = predictions.merge(internal, on="sample_barcode", validate="one_to_one")
    rows = []
    for major_group, group in joined.groupby("true_major_cancer_group", sort=True):
        rows.append(
            {
                "major_cancer_group": major_group,
                "n_samples": len(group),
                "top1_major_group_accuracy": group["correct_top_major_group"].mean(),
                "median_top1_margin": group["top1_top2_major_group_margin"].median(),
                "median_entropy": group["neighbor_entropy"].median(),
                "high_ambiguity_fraction": (group["ambiguity_class"] == "high_ambiguity").mean(),
                "common_confusions": common_confusions(
                    group, "true_major_cancer_group", "predicted_top_major_group"
                ),
                "notes": "Broad major-group validation reuses the same neighbors and is not independent evidence.",
            }
        )
    return pd.DataFrame(rows)


def confusion_matrix_table(
    predictions: pd.DataFrame,
    true_column: str,
    predicted_column: str,
    row_label: str,
    category_order: list[str] | None = None,
) -> pd.DataFrame:
    categories = category_order or sorted(
        set(predictions[true_column].astype(str)) | set(predictions[predicted_column].astype(str))
    )
    matrix = pd.crosstab(
        predictions[true_column].astype(str),
        predictions[predicted_column].astype(str),
        dropna=False,
    ).reindex(index=categories, columns=categories, fill_value=0)
    matrix.index.name = row_label
    return matrix.reset_index()


def overall_metrics(
    predictions: pd.DataFrame,
    internal: pd.DataFrame,
    top_k: int,
    metric: str,
    exclude_same_patient: bool,
) -> pd.DataFrame:
    joined = predictions.merge(internal, on="sample_barcode", validate="one_to_one")
    ambiguity_counts = predictions["ambiguity_class"].value_counts(normalize=True)
    limitations = (
        "Internal TCGA leave-one-out validation only; not external clinical validation. Queries and references share "
        "cohort, assays, preprocessing, feature engineering, and global imputation/scaling parameters. Self-neighbors "
        "are excluded, but no independent cohort, prospective specimens, clinical assay harmonization, or diagnostic "
        "calibration is evaluated."
    )
    rows = [
        ("validation_scope", "internal_tcga_pseudo_unknown_leave_one_out", "Not external clinical validation."),
        ("n_samples_validated", len(predictions), "Number of pseudo-unknown query samples."),
        ("n_projects", predictions["true_project_code"].nunique(), "Known TCGA project labels."),
        ("n_lineage_groups", predictions["true_cup_lineage_group"].nunique(), "Configured CUP lineage groups."),
        ("top_k", top_k, "Neighbors summarized per query after exclusion."),
        ("distance_metric", metric, "Exact neighbors computed from the locked scaled matrix."),
        ("exclude_same_patient", exclude_same_patient, "Self-neighbor exclusion is always enabled."),
        ("top1_project_accuracy", predictions["correct_top_project"].mean(), "Exact project recovery."),
        ("top3_project_accuracy", predictions["true_project_in_top3"].mean(), "True project among top three project summaries."),
        ("top5_project_accuracy", predictions["true_project_in_top5"].mean(), "True project among top five project summaries."),
        ("top1_lineage_accuracy", predictions["correct_top_lineage"].mean(), "Configured lineage recovery."),
        ("top3_lineage_accuracy", predictions["true_lineage_in_top3"].mean(), "True lineage among top three lineage summaries."),
        ("top1_major_group_accuracy", predictions["correct_top_major_group"].mean(), "Broad major-group recovery."),
        ("median_project_margin", predictions["top1_top2_project_margin"].median(), "Median project top-1/top-2 score margin."),
        ("median_lineage_margin", predictions["top1_top2_lineage_margin"].median(), "Median lineage top-1/top-2 score margin."),
        ("median_true_project_rank", joined["true_project_rank"].median(), "Median rank of the known project."),
        ("high_ambiguity_fraction", ambiguity_counts.get("high_ambiguity", 0.0), "Fraction assigned high ambiguity."),
        ("moderate_ambiguity_fraction", ambiguity_counts.get("moderate_ambiguity", 0.0), "Fraction assigned moderate ambiguity."),
        ("low_ambiguity_fraction", ambiguity_counts.get("low_ambiguity", 0.0), "Fraction assigned low ambiguity."),
        (
            "uninterpretable_due_to_missing_features_fraction",
            ambiguity_counts.get("uninterpretable_due_to_missing_features", 0.0),
            "Fraction below the CUP feature-coverage threshold.",
        ),
        ("limitations", limitations, "Performance may be optimistic relative to external clinical use."),
    ]
    return pd.DataFrame(rows, columns=["metric", "value", "notes"])


def ambiguity_performance(predictions: pd.DataFrame) -> pd.DataFrame:
    frame = predictions.copy()
    frame["feature_coverage_bin"] = pd.cut(
        frame["feature_coverage"],
        bins=[-np.inf, 0.50, 0.80, 0.95, np.inf],
        labels=["below_50pct", "50_to_below_80pct", "80_to_below_95pct", "95_to_100pct"],
        right=False,
    )
    frame["project_margin_bin"] = pd.cut(
        frame["top1_top2_project_margin"],
        bins=[-np.inf, 0.05, 0.10, 0.20, np.inf],
        labels=["below_0.05", "0.05_to_below_0.10", "0.10_to_below_0.20", "0.20_or_higher"],
        right=False,
    )
    frame["entropy_bin"] = pd.cut(
        frame["neighbor_entropy"],
        bins=[-np.inf, 0.25, 0.50, 0.75, np.inf],
        labels=["below_0.25", "0.25_to_below_0.50", "0.50_to_below_0.75", "0.75_or_higher"],
        right=False,
    )
    dimensions = [
        ("ambiguity_class", "ambiguity_class"),
        ("feature_coverage_bin", "feature_coverage_bin"),
        ("top1_top2_project_margin_bin", "project_margin_bin"),
        ("neighbor_entropy_bin", "entropy_bin"),
    ]
    rows: list[dict[str, Any]] = []
    for dimension, column in dimensions:
        for stratum, group in frame.groupby(column, observed=True, sort=False):
            rows.append(
                {
                    "analysis_dimension": dimension,
                    "stratum": str(stratum),
                    "n_samples": len(group),
                    "top1_project_accuracy": group["correct_top_project"].mean(),
                    "top1_lineage_accuracy": group["correct_top_lineage"].mean(),
                    "top1_major_group_accuracy": group["correct_top_major_group"].mean(),
                    "median_project_margin": group["top1_top2_project_margin"].median(),
                    "median_lineage_margin": group["top1_top2_lineage_margin"].median(),
                    "median_entropy": group["neighbor_entropy"].median(),
                    "low_ambiguity_project_accuracy": pd.NA,
                    "high_ambiguity_project_accuracy": pd.NA,
                    "low_ambiguity_more_accurate_than_high_ambiguity": pd.NA,
                    "interpretation": "Internal descriptive calibration-style summary; not probability calibration.",
                }
            )
    low = frame[frame["ambiguity_class"] == "low_ambiguity"]
    high = frame[frame["ambiguity_class"] == "high_ambiguity"]
    comparison = pd.NA
    low_accuracy = np.nan
    high_accuracy = np.nan
    if len(low) and len(high):
        low_accuracy = float(low["correct_top_project"].mean())
        high_accuracy = float(high["correct_top_project"].mean())
        comparison = bool(low_accuracy > high_accuracy)
    rows.append(
        {
            "analysis_dimension": "ambiguity_comparison",
            "stratum": "low_vs_high_ambiguity",
            "n_samples": len(low) + len(high),
            "top1_project_accuracy": np.nan,
            "top1_lineage_accuracy": np.nan,
            "top1_major_group_accuracy": np.nan,
            "median_project_margin": np.nan,
            "median_lineage_margin": np.nan,
            "median_entropy": np.nan,
            "low_ambiguity_project_accuracy": low_accuracy,
            "high_ambiguity_project_accuracy": high_accuracy,
            "low_ambiguity_more_accurate_than_high_ambiguity": comparison,
            "interpretation": (
                "Direct comparison of top-1 project accuracy in low- and high-ambiguity strata."
            ),
        }
    )
    return pd.DataFrame(rows)


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


def plot_confusion_matrix(matrix: pd.DataFrame, path: Path, title: str) -> None:
    row_label = matrix.columns[0]
    values = matrix.iloc[:, 1:].to_numpy(dtype=float)
    row_totals = values.sum(axis=1, keepdims=True)
    normalized = np.divide(values, row_totals, out=np.zeros_like(values), where=row_totals > 0)
    labels = matrix[row_label].astype(str).tolist()
    size = max(8.0, min(15.0, len(labels) * 0.38 + 4.0))
    figure_style()
    fig, ax = plt.subplots(figsize=(size, size * 0.82))
    image = ax.imshow(normalized, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(labels)), labels=labels, rotation=90, fontsize=max(4, 8 - len(labels) // 10))
    ax.set_yticks(np.arange(len(labels)), labels=labels, fontsize=max(4, 8 - len(labels) // 10))
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("Known label")
    ax.set_title(title)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    colorbar.set_label("Row-normalized fraction")
    fig.text(
        0.01,
        0.01,
        "Internal TCGA pseudo-unknown validation; row normalization is descriptive, not external clinical performance.",
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout(rect=[0, 0.035, 1, 1])
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_project_topk(performance: pd.DataFrame, path: Path) -> None:
    data = performance.sort_values("top1_project_accuracy", ascending=True).reset_index(drop=True)
    y = np.arange(len(data))
    figure_style()
    fig, ax = plt.subplots(figsize=(10, max(8, len(data) * 0.28)))
    height = 0.24
    ax.barh(y - height, data["top1_project_accuracy"], height=height, label="Top 1", color="#B44D5E")
    ax.barh(y, data["top3_project_accuracy"], height=height, label="Top 3", color="#3E7898")
    ax.barh(y + height, data["top5_project_accuracy"], height=height, label="Top 5", color="#4F7F5C")
    ax.set_yticks(y, data["project_code"])
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Accuracy")
    ax.set_title("Known-project recovery by TCGA project")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_lineage_accuracy(performance: pd.DataFrame, path: Path) -> None:
    data = performance.sort_values("top1_lineage_accuracy", ascending=True).reset_index(drop=True)
    y = np.arange(len(data))
    figure_style()
    fig, ax = plt.subplots(figsize=(10, max(6, len(data) * 0.34)))
    ax.barh(y - 0.16, data["top1_lineage_accuracy"], height=0.32, label="Top 1", color="#3E7898")
    ax.barh(y + 0.16, data["top3_lineage_accuracy"], height=0.32, label="Top 3", color="#4F7F5C")
    ax.set_yticks(y, data["cup_lineage_group"])
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Accuracy")
    ax.set_title("Known-lineage recovery by configured CUP lineage")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_ambiguity_accuracy(summary: pd.DataFrame, path: Path) -> None:
    order = ["low_ambiguity", "moderate_ambiguity", "high_ambiguity", "uninterpretable_due_to_missing_features"]
    data = summary[summary["analysis_dimension"] == "ambiguity_class"].copy()
    data["_order"] = data["stratum"].map({value: index for index, value in enumerate(order)})
    data = data.sort_values("_order")
    x = np.arange(len(data))
    figure_style()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    width = 0.25
    ax.bar(x - width, data["top1_project_accuracy"], width, label="Project", color="#B44D5E")
    ax.bar(x, data["top1_lineage_accuracy"], width, label="Lineage", color="#3E7898")
    ax.bar(x + width, data["top1_major_group_accuracy"], width, label="Major group", color="#4F7F5C")
    ax.set_xticks(x, [value.replace("_", " ") for value in data["stratum"]], rotation=15, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Accuracy")
    ax.set_title("Accuracy by ambiguity class")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_metric_by_correctness(
    predictions: pd.DataFrame,
    value_column: str,
    path: Path,
    title: str,
    ylabel: str,
) -> None:
    incorrect = predictions.loc[~predictions["correct_top_project"], value_column].to_numpy(dtype=float)
    correct = predictions.loc[predictions["correct_top_project"], value_column].to_numpy(dtype=float)
    figure_style()
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    boxes = ax.boxplot(
        [incorrect, correct],
        tick_labels=[f"Incorrect\n(n={len(incorrect)})", f"Correct\n(n={len(correct)})"],
        patch_artist=True,
        showfliers=False,
    )
    for patch, color in zip(boxes["boxes"], ["#B44D5E", "#4F7F5C"], strict=True):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    fig.tight_layout()
    ensure_dir(path.parent)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def top_confusion_pairs(predictions: pd.DataFrame, limit: int = 12) -> pd.DataFrame:
    wrong = predictions[~predictions["correct_top_project"]]
    if wrong.empty:
        return pd.DataFrame(columns=["true_project", "predicted_project", "n_samples"])
    return (
        wrong.groupby(["true_project_code", "predicted_top_project"], sort=False)
        .size()
        .rename("n_samples")
        .reset_index()
        .rename(columns={"true_project_code": "true_project", "predicted_top_project": "predicted_project"})
        .sort_values(["n_samples", "true_project", "predicted_project"], ascending=[False, True, True])
        .head(limit)
        .reset_index(drop=True)
    )


def prioritized_confusion_analysis(predictions: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        ("COAD vs READ", {"COAD"}, {"READ"}),
        ("PAAD vs CHOL/STAD/COAD/READ", {"PAAD"}, {"CHOL", "STAD", "COAD", "READ"}),
        ("LUAD vs LUSC/HNSC", {"LUAD"}, {"LUSC", "HNSC"}),
        ("HNSC vs LUSC/CESC/ESCA", {"HNSC"}, {"LUSC", "CESC", "ESCA"}),
        ("KIRC vs KIRP/KICH", {"KIRC"}, {"KIRP", "KICH"}),
        ("BRCA vs OV/UCEC/UCS", {"BRCA"}, {"OV", "UCEC", "UCS"}),
        ("GBM vs LGG", {"GBM"}, {"LGG"}),
        ("SKCM vs UVM", {"SKCM"}, {"UVM"}),
        ("LAML vs DLBC", {"LAML"}, {"DLBC"}),
    ]
    rows = []
    true_labels = predictions["true_project_code"].astype(str)
    predicted_labels = predictions["predicted_top_project"].astype(str)
    for label, group_a, group_b in comparisons:
        a_as_b = int((true_labels.isin(group_a) & predicted_labels.isin(group_b)).sum())
        b_as_a = int((true_labels.isin(group_b) & predicted_labels.isin(group_a)).sum())
        n_queries = int(true_labels.isin(group_a | group_b).sum())
        rows.append(
            {
                "comparison": label,
                "group_a_to_group_b": a_as_b,
                "group_b_to_group_a": b_as_a,
                "n_queries_in_comparison": n_queries,
                "cross_confusion_fraction": (a_as_b + b_as_a) / n_queries if n_queries else np.nan,
            }
        )
    return pd.DataFrame(rows)


def markdown_table(frame: pd.DataFrame, columns: list[str], digits: int = 3) -> str:
    view = frame[columns].copy()
    for column in view.select_dtypes(include=["float"]).columns:
        view[column] = view[column].map(lambda value: f"{value:.{digits}f}" if pd.notna(value) else "NA")
    headers = [column.replace("_", " ").title() for column in columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in view.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def build_report(
    predictions: pd.DataFrame,
    project_summary: pd.DataFrame,
    lineage_summary: pd.DataFrame,
    major_summary: pd.DataFrame,
    metrics: pd.DataFrame,
    ambiguity_summary: pd.DataFrame,
    top_k: int,
    metric: str,
    exclude_same_patient: bool,
) -> str:
    metric_map = dict(zip(metrics["metric"].astype(str), metrics["value"].astype(str), strict=False))
    project_best = project_summary.sort_values(["top1_project_accuracy", "n_samples"], ascending=[False, False]).head(8)
    project_worst = project_summary.sort_values(["top1_project_accuracy", "n_samples"], ascending=[True, False]).head(8)
    lineage_ranked = lineage_summary.sort_values(["top1_lineage_accuracy", "n_samples"], ascending=[False, False])
    confusion_pairs = top_confusion_pairs(predictions)
    prioritized_confusions = prioritized_confusion_analysis(predictions)
    ambiguity_rows = ambiguity_summary[ambiguity_summary["analysis_dimension"] == "ambiguity_class"]
    comparison = ambiguity_summary[ambiguity_summary["analysis_dimension"] == "ambiguity_comparison"]
    comparison_text = "not assessable because one class was absent"
    if not comparison.empty and pd.notna(comparison.iloc[0]["low_ambiguity_more_accurate_than_high_ambiguity"]):
        comparison_text = str(bool(comparison.iloc[0]["low_ambiguity_more_accurate_than_high_ambiguity"]))

    return f"""# Known-Primary TCGA Projection Internal Validation

Generated by `scripts/23_validate_known_primary_projection.py`.

## Purpose

This analysis tests whether the locked 97-feature TCGA sample-level atlas and nearest-neighbor/CUP interpretation logic recover known TCGA project, configured CUP lineage, and broad major-group labels when atlas samples are treated as pseudo-unknown queries.

This is internal validation using TCGA-derived samples. It is not external clinical validation, does not establish clinical diagnostic accuracy, and does not validate performance on metastatic, treated, small-biopsy, cytology, or independently processed clinical specimens.

## Validation Design

- Reference and query source: the same locked TCGA atlas.
- Query design: leave-one-out pseudo-unknown projection against the full atlas.
- Self-neighbor exclusion: always enabled.
- Same-patient exclusion: `{exclude_same_patient}`.
- Distance metric: `{metric}`.
- Neighbors summarized per query: {top_k}.
- Feature coverage: observed pre-imputation fraction of the 97 locked features.
- Scoring: `1 / (1 + distance)` normalized across returned neighbors; scores are descriptive weights, not probabilities.

The global atlas medians and scaling parameters were not refit per query. This preserves the locked clinical-query contract but means the pseudo-unknown samples contributed to reference preprocessing estimates, a source of internal-validation leakage that may make performance optimistic.

## Cohort Composition

The validation included {metric_map.get("n_samples_validated", "NA")} samples from {metric_map.get("n_projects", "NA")} TCGA projects, {metric_map.get("n_lineage_groups", "NA")} configured CUP lineage groups, and {predictions["true_major_cancer_group"].nunique()} broad major cancer groups.

## Overall Performance

| Metric | Value |
| --- | ---: |
| Top-1 project accuracy | {float(metric_map.get("top1_project_accuracy", "nan")):.3f} |
| Top-3 project accuracy | {float(metric_map.get("top3_project_accuracy", "nan")):.3f} |
| Top-5 project accuracy | {float(metric_map.get("top5_project_accuracy", "nan")):.3f} |
| Top-1 lineage accuracy | {float(metric_map.get("top1_lineage_accuracy", "nan")):.3f} |
| Top-3 lineage accuracy | {float(metric_map.get("top3_lineage_accuracy", "nan")):.3f} |
| Top-1 major-group accuracy | {float(metric_map.get("top1_major_group_accuracy", "nan")):.3f} |
| Median project margin | {float(metric_map.get("median_project_margin", "nan")):.3f} |
| Median lineage margin | {float(metric_map.get("median_lineage_margin", "nan")):.3f} |

These values measure label recovery within the same TCGA molecular ecosystem. They must not be reported as external clinical accuracy.

## Project-Level Performance

Highest observed top-1 project recovery:

{markdown_table(project_best, ["project_code", "n_samples", "top1_project_accuracy", "top3_project_accuracy", "top5_project_accuracy", "common_confusions"])}

Lowest observed top-1 project recovery:

{markdown_table(project_worst, ["project_code", "n_samples", "top1_project_accuracy", "top3_project_accuracy", "top5_project_accuracy", "common_confusions"])}

## Lineage-Level Performance

{markdown_table(lineage_ranked, ["cup_lineage_group", "n_samples", "top1_lineage_accuracy", "top3_lineage_accuracy", "high_ambiguity_fraction", "common_confusions"])}

Lineage-level results are easier than exact project recovery when multiple related projects share one configured lineage. This grouping is heuristic and does not replace organ-specific pathology.

## Major-Group Performance

{markdown_table(major_summary.sort_values("top1_major_group_accuracy", ascending=False), ["major_cancer_group", "n_samples", "top1_major_group_accuracy", "high_ambiguity_fraction", "common_confusions"])}

Major-group performance reuses the same neighbors and should not be interpreted as independent validation.

## Common Confusions

{markdown_table(confusion_pairs, ["true_project", "predicted_project", "n_samples"], digits=0) if not confusion_pairs.empty else "No incorrect top-project predictions were observed."}

Prespecified directional confusion families (`group_a_to_group_b` and `group_b_to_group_a` follow the order shown in the comparison label):

{markdown_table(prioritized_confusions, ["comparison", "group_a_to_group_b", "group_b_to_group_a", "n_queries_in_comparison", "cross_confusion_fraction"])}

Project confusion may reflect genuine shared molecular programs, broad TCGA labels, sample composition, or feature limitations. Particular attention should be paid to colorectal, GI/pancreatobiliary, squamous, kidney, gynecologic, glioma, melanoma, and hematolymphoid project pairs.

## Ambiguity Analysis

{markdown_table(ambiguity_rows, ["stratum", "n_samples", "top1_project_accuracy", "top1_lineage_accuracy", "top1_major_group_accuracy", "median_project_margin", "median_entropy"])}

Low-ambiguity cases had higher top-1 project accuracy than high-ambiguity cases: `{comparison_text}`. This is a calibration-style descriptive check, not probability calibration.

## Failure Modes

- Exact project recovery is limited where related TCGA projects share molecular programs or configured lineages.
- Feature coverage varies because the raw atlas contains missing values subsequently filled by TCGA medians.
- Mutation-count features are proxies without callable-territory normalization.
- Broad pathway scores, arm-level copy-number burden, and driver counts omit variant pathogenicity, signatures, fusions, and allele-specific focal context.
- Same-cohort technical and biological similarity may inflate performance relative to independently processed clinical specimens.
- Rare projects and heterogeneous TCGA cohorts can produce unstable estimates despite full leave-one-out evaluation.

## Limitations

- This is internal TCGA validation, not external clinical validation.
- Query samples and references share cohort selection, assays, preprocessing, and global imputation/scaling parameters.
- Excluding self-neighbors prevents identity matching but does not remove cohort or preprocessing leakage.
- TCGA is primary-tumor enriched and incompletely represents metastatic and real-world diagnostic specimens.
- Project labels are not identical to histologic or tissue-of-origin truth in every case.
- CUP lineage definitions and ambiguity thresholds remain heuristic and uncalibrated.
- The workflow does not evaluate raw clinical data processing, clinical assay failure, actionability, or prospective outcomes.

## Recommended Improvements Before Clinical Use

1. Validate on independent known-primary WES/WTS cohorts processed with the intended clinical pipeline.
2. Include metastatic, small-biopsy, cytology, low-purity, and treatment-exposed specimens.
3. Repeat validation with training-only imputation/scaling or an independently frozen reference to remove preprocessing leakage.
4. Evaluate top-k, Euclidean/cosine, feature-class weighting, missing-modality, and patient-level exclusion sensitivity.
5. Calibrate abstention and ambiguity thresholds on held-out cohorts before any tissue-of-origin claim.
6. Add lineage-marker, fusion, signature, focal copy-number, and variant-level evidence under separate validation.
7. Integrate morphology, immunophenotype, imaging, and clinical findings; do not use this output as a final diagnosis.
"""


def generate_figures(
    paths: ValidationPaths,
    project_confusion: pd.DataFrame,
    lineage_confusion: pd.DataFrame,
    project_summary: pd.DataFrame,
    lineage_summary: pd.DataFrame,
    ambiguity_summary: pd.DataFrame,
    predictions: pd.DataFrame,
) -> None:
    plot_confusion_matrix(
        project_confusion,
        paths.project_confusion_figure,
        "Known-project confusion matrix (row-normalized)",
    )
    plot_confusion_matrix(
        lineage_confusion,
        paths.lineage_confusion_figure,
        "Known CUP-lineage confusion matrix (row-normalized)",
    )
    plot_project_topk(project_summary, paths.topk_project_figure)
    plot_lineage_accuracy(lineage_summary, paths.lineage_accuracy_figure)
    plot_ambiguity_accuracy(ambiguity_summary, paths.ambiguity_accuracy_figure)
    plot_metric_by_correctness(
        predictions,
        "top1_top2_project_margin",
        paths.margin_figure,
        "Project margin versus top-project correctness",
        "Top-1/top-2 project score margin",
    )
    plot_metric_by_correctness(
        predictions,
        "neighbor_entropy",
        paths.entropy_figure,
        "Lineage entropy versus top-project correctness",
        "Normalized lineage entropy",
    )


def run_validation(
    root: Path,
    top_k: int = 50,
    metric: str = "euclidean",
    max_samples: int | None = None,
    seed: int = 123,
    exclude_same_patient: bool = False,
    output_dir: Path | None = None,
) -> dict[str, Path]:
    paths = default_paths(root, output_dir=output_dir)
    reference = load_reference_inputs(paths)
    if metric not in set(reference["available_metrics"]):
        raise ValueError(
            f"Metric {metric!r} is not in locked atlas available metrics {reference['available_metrics']}"
        )
    cohort, query_indices = build_validation_cohort(reference["metadata"], max_samples, seed)
    logging.info(
        "Starting known-primary validation: queries=%d, reference=%d, features=%d, top_k=%d, metric=%s",
        len(query_indices),
        len(reference["metadata"]),
        len(reference["feature_order"]),
        top_k,
        metric,
    )
    neighbor_indices, neighbor_distances = compute_neighbor_indices(
        reference["matrix"],
        query_indices,
        top_k,
        metric,
        patient_ids=reference["metadata"]["patient_barcode"].astype(str).to_numpy(),
        exclude_same_patient=exclude_same_patient,
    )
    predictions, internal = evaluate_queries(
        reference["metadata"],
        query_indices,
        neighbor_indices,
        neighbor_distances,
        reference["feature_coverage"],
        top_k,
        metric,
        exclude_same_patient,
    )
    project_summary = project_performance(predictions, internal)
    lineage_summary = lineage_performance(predictions, internal)
    major_summary = major_group_performance(predictions, internal)
    metrics = overall_metrics(predictions, internal, top_k, metric, exclude_same_patient)
    project_order = sorted(reference["project_to_lineage"])
    lineage_order = reference["lineage_order"]
    major_order = sorted(reference["metadata"]["major_cancer_group"].astype(str).unique())
    project_confusion = confusion_matrix_table(
        predictions, "true_project_code", "predicted_top_project", "true_project_code", project_order
    )
    lineage_confusion = confusion_matrix_table(
        predictions,
        "true_cup_lineage_group",
        "predicted_top_lineage",
        "true_cup_lineage_group",
        lineage_order,
    )
    major_confusion = confusion_matrix_table(
        predictions,
        "true_major_cancer_group",
        "predicted_top_major_group",
        "true_major_cancer_group",
        major_order,
    )
    ambiguity_summary = ambiguity_performance(predictions)

    write_tsv(cohort, paths.cohort)
    write_tsv(predictions, paths.predictions)
    write_tsv(project_summary, paths.project_performance)
    write_tsv(lineage_summary, paths.lineage_performance)
    write_tsv(major_summary, paths.major_group_performance)
    write_tsv(metrics, paths.overall_metrics)
    write_tsv(project_confusion, paths.project_confusion)
    write_tsv(lineage_confusion, paths.lineage_confusion)
    write_tsv(major_confusion, paths.major_group_confusion)
    write_tsv(ambiguity_summary, paths.ambiguity_performance)
    generate_figures(
        paths,
        project_confusion,
        lineage_confusion,
        project_summary,
        lineage_summary,
        ambiguity_summary,
        predictions,
    )
    write_text(
        build_report(
            predictions,
            project_summary,
            lineage_summary,
            major_summary,
            metrics,
            ambiguity_summary,
            top_k,
            metric,
            exclude_same_patient,
        ),
        paths.report,
    )
    logging.info("Known-primary internal validation complete: %s", paths.output_dir)
    return {field: getattr(paths, field) for field in paths.__dataclass_fields__ if field != "output_dir"}


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
        description="Internally validate TCGA nearest-neighbor and CUP label recovery using pseudo-unknown samples."
    )
    parser.add_argument("--top-k", type=int, default=50, help="Neighbors summarized per query after exclusions.")
    parser.add_argument(
        "--metric", choices=["euclidean", "cosine"], default="euclidean", help="Exact neighbor distance metric."
    )
    parser.add_argument(
        "--max-samples",
        type=parse_max_samples,
        default=None,
        metavar="all|N",
        help="Validate all samples (default) or a deterministic project-stratified query subset.",
    )
    parser.add_argument("--seed", type=int, default=123, help="Seed for stratified query subsampling.")
    parser.add_argument(
        "--exclude-same-patient",
        action="store_true",
        help="Exclude every reference sample from the query patient's barcode, in addition to self exclusion.",
    )
    parser.add_argument("--root", type=Path, default=None, help="Project root; defaults to auto-detection.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional alternate output directory.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root or Path(__file__).resolve())
    run_validation(
        root,
        top_k=args.top_k,
        metric=args.metric,
        max_samples=args.max_samples,
        seed=args.seed,
        exclude_same_patient=args.exclude_same_patient,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()

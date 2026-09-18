#!/usr/bin/env python

from __future__ import annotations

import argparse
import logging
import math
import os
import re
import tempfile
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcga_cancer_phylogeny_matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


CONFIG = {
    "top_k": 50,
    "meaningful_project_score_min": 0.05,
    "uninterpretable_coverage_below": 0.50,
    "low_ambiguity_top1_min": 0.65,
    "low_ambiguity_margin_min": 0.25,
    "low_ambiguity_entropy_max": 0.60,
    "high_ambiguity_top1_below": 0.40,
    "high_ambiguity_margin_below": 0.10,
    "high_ambiguity_entropy_min": 0.80,
    "high_ambiguity_meaningful_projects_min": 8,
    "high_support_top1_min": 0.65,
    "high_support_margin_min": 0.25,
    "high_support_coverage_min": 0.85,
    "high_support_neighbor_concentration_min": 0.65,
    "moderate_support_top1_min": 0.45,
    "moderate_support_margin_min": 0.15,
    "moderate_support_coverage_min": 0.65,
    "moderate_support_neighbor_concentration_min": 0.40,
    "differential_relative_support_min": 0.50,
}

FEATURE_VECTOR_ID_COLUMNS = ["case_id", "sample_id"]
NEIGHBOR_REQUIRED_COLUMNS = [
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
]
PROJECT_SUMMARY_REQUIRED_COLUMNS = [
    "project_code",
    "project_id",
    "n_neighbors",
    "median_distance",
    "min_distance",
    "weighted_similarity_score",
    "rank",
]
MAJOR_GROUP_SUMMARY_REQUIRED_COLUMNS = [
    "major_cancer_group",
    "n_neighbors",
    "median_distance",
    "min_distance",
    "weighted_similarity_score",
    "rank",
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
RANKED_LINEAGE_COLUMNS = [
    "rank",
    "lineage_or_project_group",
    "supporting_tcga_projects",
    "similarity_score",
    "confidence_tier",
    "evidence_summary",
    "caveats",
]
EVIDENCE_TYPES = [
    "nearest-neighbor TCGA project evidence",
    "major cancer group evidence",
    "expression/pathway evidence",
    "driver mutation evidence",
    "copy-number/aneuploidy evidence",
    "immune/stromal/EMT evidence",
    "contradictory evidence",
    "missing evidence",
]
EVIDENCE_STRENGTH_ORDER = {"not assessed": 0, "weak": 1, "moderate": 2, "strong": 3}


@dataclass
class CUPPaths:
    clinical_feature_vector: Path
    clinical_feature_vector_scaled: Path
    clinical_feature_missingness: Path
    nearest_neighbors: Path
    nearest_project_summary: Path
    nearest_major_group_summary: Path
    projection_qc_summary: Path
    clinical_metadata: Path
    reference_matrix: Path
    reference_metadata: Path
    feature_definitions: Path
    lineage_groups: Path
    evidence_rules: Path
    output_dir: Path
    ranked_lineages: Path
    ambiguity_metrics: Path
    molecular_evidence: Path
    differential_summary: Path
    concordance_summary: Path
    report: Path
    top_project_figure: Path
    major_group_figure: Path
    ambiguity_figure: Path
    evidence_heatmap: Path


def default_paths(
    root: Path,
    case_id: str,
    projection_dir: Path | None = None,
    query_dir: Path | None = None,
    output_dir: Path | None = None,
    lineage_config: Path | None = None,
    evidence_config: Path | None = None,
) -> CUPPaths:
    projection = projection_dir or project_path("results", "clinical_projection", case_id, root=root)
    query = query_dir or project_path("data", "clinical_queries", case_id, root=root)
    output = output_dir or projection / "cup_interpretation"
    return CUPPaths(
        clinical_feature_vector=projection / "clinical_feature_vector.tsv",
        clinical_feature_vector_scaled=projection / "clinical_feature_vector_scaled.tsv",
        clinical_feature_missingness=projection / "clinical_feature_missingness.tsv",
        nearest_neighbors=projection / "tcga_nearest_neighbors.tsv",
        nearest_project_summary=projection / "tcga_nearest_project_summary.tsv",
        nearest_major_group_summary=projection / "tcga_nearest_major_group_summary.tsv",
        projection_qc_summary=projection / "clinical_projection_qc_summary.tsv",
        clinical_metadata=query / "clinical_metadata.tsv",
        reference_matrix=project_path(
            "results", "reference_atlas", "tcga_sample_reference_matrix.tsv.gz", root=root
        ),
        reference_metadata=project_path(
            "results", "reference_atlas", "tcga_sample_reference_metadata.tsv", root=root
        ),
        feature_definitions=project_path(
            "results", "reference_atlas", "tcga_reference_feature_definitions.yaml", root=root
        ),
        lineage_groups=lineage_config
        or project_path("config", "cup_lineage_groups.yaml", root=root),
        evidence_rules=evidence_config
        or project_path("config", "cup_feature_evidence_rules.yaml", root=root),
        output_dir=output,
        ranked_lineages=output / "cup_ranked_lineage_interpretation.tsv",
        ambiguity_metrics=output / "cup_ambiguity_metrics.tsv",
        molecular_evidence=output / "cup_molecular_evidence_table.tsv",
        differential_summary=output / "cup_differential_diagnosis_constrained_summary.tsv",
        concordance_summary=output / "cup_molecular_pathologic_concordance.tsv",
        report=output / "cup_molecular_interpretation_report.md",
        top_project_figure=output / "cup_top_project_similarity_barplot.pdf",
        major_group_figure=output / "cup_major_group_similarity_barplot.pdf",
        ambiguity_figure=output / "cup_ambiguity_summary.pdf",
        evidence_heatmap=output / "cup_evidence_heatmap.pdf",
    )


def validate_case_id(case_id: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", case_id) or case_id in {".", ".."}:
        raise ValueError("case_id must contain only letters, numbers, periods, underscores, or hyphens")
    return case_id


def read_yaml(path: Path, label: str) -> dict[str, Any]:
    require_file(path, label)
    with path.open(encoding="utf-8") as handle:
        parsed = yaml.safe_load(handle)
    if not isinstance(parsed, dict):
        raise ValueError(f"{label} must contain a YAML mapping")
    return parsed


def read_tsv(path: Path, label: str, dtype: Any = None) -> pd.DataFrame:
    require_file(path, label)
    data = pd.read_csv(path, sep="\t", dtype=dtype, keep_default_na=True)
    if data.empty:
        raise ValueError(f"{label} is empty: {path}")
    return data


def write_tsv(data: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_text(text: str, path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")
    logging.info("Wrote %s", path)


def require_columns(data: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise ValueError(f"{label} lacks required columns: {missing}")


def is_missing(value: Any) -> bool:
    return pd.isna(value) or (isinstance(value, str) and value.strip() == "")


def bool_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if is_missing(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def finite_numeric(data: pd.DataFrame, columns: list[str], label: str) -> None:
    for column in columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    if data[columns].isna().any().any() or not np.isfinite(data[columns].to_numpy(dtype=float)).all():
        raise ValueError(f"{label} contains missing or non-finite numeric values in {columns}")


def one_case_row(data: pd.DataFrame, case_id: str, sample_id: str | None, label: str) -> pd.Series:
    require_columns(data, ["case_id", "sample_id"], label)
    matching = data[data["case_id"].astype(str) == case_id]
    if sample_id is not None:
        matching = matching[matching["sample_id"].astype(str) == sample_id]
    if len(matching) != 1:
        raise ValueError(f"{label} must contain exactly one matching case/sample row")
    return matching.iloc[0]


def feature_name(record: dict[str, Any]) -> str:
    name = record.get("name", record.get("feature_name"))
    if not name:
        raise ValueError("Feature definition lacks name/feature_name")
    return str(name)


def load_inputs(paths: CUPPaths, case_id: str, top_k: int) -> dict[str, Any]:
    metadata = read_tsv(paths.clinical_metadata, "clinical metadata", dtype=str)
    require_columns(
        metadata,
        ["case_id", "sample_id", "specimen_type", "submitted_diagnosis", "notes"],
        "clinical metadata",
    )
    matching_metadata = metadata[metadata["case_id"].astype(str) == case_id]
    if len(matching_metadata) != 1:
        raise ValueError(f"clinical_metadata.tsv must contain exactly one row for case_id={case_id!r}")
    clinical_metadata = matching_metadata.iloc[0].to_dict()
    sample_id = str(clinical_metadata["sample_id"]).strip()
    if not sample_id:
        raise ValueError("Clinical metadata sample_id is blank")
    clinical_metadata.setdefault("differential_diagnosis", None)
    clinical_metadata.setdefault("tumor_purity", None)

    definitions = read_yaml(paths.feature_definitions, "TCGA reference feature definitions")
    records = definitions.get("features")
    if not isinstance(records, list) or not records:
        raise ValueError("TCGA reference feature definitions lack a non-empty features list")
    feature_records = [dict(record) for record in records]
    feature_order = [feature_name(record) for record in feature_records]
    if len(feature_order) != len(set(feature_order)):
        raise ValueError("TCGA reference feature definitions contain duplicate feature names")

    raw_vector = read_tsv(paths.clinical_feature_vector, "clinical feature vector")
    scaled_vector = read_tsv(paths.clinical_feature_vector_scaled, "scaled clinical feature vector")
    raw_row = one_case_row(raw_vector, case_id, sample_id, "clinical feature vector")
    scaled_row = one_case_row(scaled_vector, case_id, sample_id, "scaled clinical feature vector")
    missing_raw = [feature for feature in feature_order if feature not in raw_vector.columns]
    missing_scaled = [feature for feature in feature_order if feature not in scaled_vector.columns]
    if missing_raw or missing_scaled:
        raise ValueError(
            f"Clinical feature vectors do not match the locked feature contract; raw missing={missing_raw}, "
            f"scaled missing={missing_scaled}"
        )
    raw_values = pd.to_numeric(raw_row[feature_order], errors="coerce")
    scaled_values = pd.to_numeric(scaled_row[feature_order], errors="coerce")
    if raw_values.isna().any() or scaled_values.isna().any():
        raise ValueError("Clinical feature vectors contain missing or non-numeric locked feature values")
    if not np.isfinite(raw_values.to_numpy(dtype=float)).all() or not np.isfinite(
        scaled_values.to_numpy(dtype=float)
    ).all():
        raise ValueError("Clinical feature vectors contain non-finite values")

    neighbors_all = read_tsv(paths.nearest_neighbors, "TCGA nearest-neighbor table")
    require_columns(neighbors_all, NEIGHBOR_REQUIRED_COLUMNS, "TCGA nearest-neighbor table")
    neighbors_all = neighbors_all[
        (neighbors_all["case_id"].astype(str) == case_id)
        & (neighbors_all["query_sample_id"].astype(str) == sample_id)
    ].copy()
    if neighbors_all.empty:
        raise ValueError("TCGA nearest-neighbor table has no rows for the requested case/sample")
    finite_numeric(neighbors_all, ["rank", "distance", "similarity_score"], "TCGA nearest-neighbor table")
    neighbors_all["rank"] = neighbors_all["rank"].astype(int)
    neighbors_all = neighbors_all.sort_values("rank", kind="stable").reset_index(drop=True)
    if neighbors_all["rank"].duplicated().any() or neighbors_all["tcga_sample_barcode"].duplicated().any():
        raise ValueError("TCGA nearest-neighbor table contains duplicate ranks or sample barcodes")
    expected_similarity = 1.0 / (1.0 + neighbors_all["distance"].to_numpy(dtype=float))
    if not np.allclose(
        expected_similarity,
        neighbors_all["similarity_score"].to_numpy(dtype=float),
        rtol=1e-8,
        atol=1e-10,
    ):
        raise ValueError("TCGA nearest-neighbor similarity scores do not equal 1/(1+distance)")
    if top_k <= 0:
        raise ValueError("top_k must be positive")
    neighbors = neighbors_all.head(min(top_k, len(neighbors_all))).copy()

    project_summary = read_tsv(paths.nearest_project_summary, "TCGA nearest-project summary")
    require_columns(project_summary, PROJECT_SUMMARY_REQUIRED_COLUMNS, "TCGA nearest-project summary")
    finite_numeric(
        project_summary,
        ["n_neighbors", "median_distance", "min_distance", "weighted_similarity_score", "rank"],
        "TCGA nearest-project summary",
    )
    project_summary["n_neighbors"] = project_summary["n_neighbors"].astype(int)
    project_summary["rank"] = project_summary["rank"].astype(int)

    major_summary = read_tsv(paths.nearest_major_group_summary, "TCGA nearest major-group summary")
    require_columns(major_summary, MAJOR_GROUP_SUMMARY_REQUIRED_COLUMNS, "TCGA nearest major-group summary")
    finite_numeric(
        major_summary,
        ["n_neighbors", "median_distance", "min_distance", "weighted_similarity_score", "rank"],
        "TCGA nearest major-group summary",
    )
    major_summary["n_neighbors"] = major_summary["n_neighbors"].astype(int)
    major_summary["rank"] = major_summary["rank"].astype(int)

    qc = read_tsv(paths.projection_qc_summary, "clinical projection QC summary")
    require_columns(qc, ["qc_section", "metric", "value"], "clinical projection QC summary")
    global_qc = (
        qc[qc["qc_section"].astype(str) == "global"]
        .drop_duplicates("metric", keep="last")
        .set_index("metric")["value"]
        .to_dict()
    )
    qc_case = str(global_qc.get("case_id", ""))
    qc_sample = str(global_qc.get("query_sample_id", ""))
    if qc_case != case_id or qc_sample != sample_id:
        raise ValueError("Clinical projection QC case/sample identifiers do not match clinical metadata")
    distance_metrics = set(neighbors_all["distance_metric"].dropna().astype(str))
    if len(distance_metrics) != 1:
        raise ValueError("TCGA nearest-neighbor table must contain exactly one distance metric")
    qc_metric = str(global_qc.get("nearest_neighbor_metric", "")).strip()
    if qc_metric and distance_metrics != {qc_metric}:
        raise ValueError("Nearest-neighbor table distance metric disagrees with projection QC")

    reference_metadata = read_tsv(paths.reference_metadata, "TCGA reference metadata")
    require_columns(reference_metadata, REFERENCE_METADATA_COLUMNS, "TCGA reference metadata")
    if reference_metadata["sample_barcode"].duplicated().any():
        raise ValueError("TCGA reference metadata contains duplicate sample barcodes")
    reference_matrix = read_tsv(paths.reference_matrix, "TCGA raw reference matrix")
    require_columns(
        reference_matrix,
        REFERENCE_METADATA_COLUMNS + feature_order,
        "TCGA raw reference matrix",
    )
    if reference_matrix["sample_barcode"].astype(str).tolist() != reference_metadata[
        "sample_barcode"
    ].astype(str).tolist():
        raise ValueError("TCGA reference matrix and metadata sample orders do not match")
    unknown_neighbors = set(neighbors_all["tcga_sample_barcode"].astype(str)) - set(
        reference_metadata["sample_barcode"].astype(str)
    )
    if unknown_neighbors:
        raise ValueError(f"Nearest-neighbor table contains samples absent from reference metadata: {sorted(unknown_neighbors)}")

    feature_audit: pd.DataFrame | None = None
    supplied_features: set[str] = set()
    missing_features: list[str] = []
    if paths.clinical_feature_missingness.exists():
        feature_audit = pd.read_csv(paths.clinical_feature_missingness, sep="\t", keep_default_na=True)
        require_columns(
            feature_audit,
            ["case_id", "sample_id", "feature_name", "expected_in_reference", "supplied", "imputed"],
            "clinical feature missingness audit",
        )
        feature_audit = feature_audit[
            (feature_audit["case_id"].astype(str) == case_id)
            & (feature_audit["sample_id"].astype(str) == sample_id)
            & feature_audit["expected_in_reference"].map(bool_value)
        ].copy()
        supplied_features = set(
            feature_audit.loc[feature_audit["supplied"].map(bool_value), "feature_name"].astype(str)
        )
        missing_features = feature_audit.loc[
            feature_audit["imputed"].map(bool_value), "feature_name"
        ].astype(str).tolist()

    return {
        "clinical_metadata": clinical_metadata,
        "sample_id": sample_id,
        "feature_records": feature_records,
        "feature_order": feature_order,
        "raw_values": raw_values.astype(float).to_dict(),
        "scaled_values": scaled_values.astype(float).to_dict(),
        "neighbors_all": neighbors_all,
        "neighbors": neighbors,
        "provided_project_summary": project_summary,
        "provided_major_summary": major_summary,
        "projection_qc": qc,
        "global_qc": global_qc,
        "reference_metadata": reference_metadata,
        "reference_matrix": reference_matrix,
        "feature_audit": feature_audit,
        "supplied_features": supplied_features,
        "missing_features": missing_features,
        "feature_availability_audited": feature_audit is not None,
    }


def load_lineage_groups(path: Path, reference_metadata: pd.DataFrame) -> dict[str, Any]:
    config = read_yaml(path, "CUP lineage-group configuration")
    lineages = config.get("lineages")
    if not isinstance(lineages, list) or not lineages:
        raise ValueError("CUP lineage-group configuration lacks a non-empty lineages list")
    names: list[str] = []
    assigned_projects: list[str] = []
    normalized: list[dict[str, Any]] = []
    for record in lineages:
        if not isinstance(record, dict) or not record.get("name"):
            raise ValueError("Each CUP lineage-group record requires a name")
        name = str(record["name"])
        projects = [str(value) for value in record.get("projects", [])]
        if not projects:
            raise ValueError(f"CUP lineage group {name} has no projects")
        names.append(name)
        assigned_projects.extend(projects)
        normalized.append(
            {
                **record,
                "name": name,
                "display_name": str(record.get("display_name", name)),
                "projects": projects,
                "aliases": [str(value) for value in record.get("aliases", [])],
                "major_cancer_groups": [str(value) for value in record.get("major_cancer_groups", [])],
                "description": str(record.get("description", "")),
                "caveat": str(record.get("caveat", "")),
            }
        )
    duplicate_names = sorted({name for name in names if names.count(name) > 1})
    duplicate_projects = sorted({project for project in assigned_projects if assigned_projects.count(project) > 1})
    if duplicate_names or duplicate_projects:
        raise ValueError(
            f"CUP lineage configuration has duplicate names/projects: names={duplicate_names}, "
            f"projects={duplicate_projects}"
        )
    reference_projects = set(reference_metadata["project_code"].astype(str))
    configured_projects = set(assigned_projects)
    if reference_projects != configured_projects:
        raise ValueError(
            "CUP lineage configuration must assign every reference project exactly once; "
            f"missing={sorted(reference_projects - configured_projects)}, "
            f"extra={sorted(configured_projects - reference_projects)}"
        )
    expected_count = int(config.get("expected_tcga_project_count", len(reference_projects)))
    if expected_count != len(reference_projects):
        raise ValueError(
            f"CUP lineage configuration expects {expected_count} projects but the reference has "
            f"{len(reference_projects)}"
        )
    project_to_lineage = {
        project: record["name"] for record in normalized for project in record["projects"]
    }
    by_name = {record["name"]: record for record in normalized}
    return {
        "raw": config,
        "lineages": normalized,
        "by_name": by_name,
        "project_to_lineage": project_to_lineage,
    }


def condition_features(condition: dict[str, Any]) -> list[str]:
    clauses = condition.get("all", condition.get("any", []))
    if not isinstance(clauses, list):
        raise ValueError("Evidence-rule condition all/any value must be a list")
    return [str(clause.get("feature", "")) for clause in clauses]


def load_evidence_rules(
    path: Path,
    available_features: set[str],
    lineage_names: set[str],
) -> dict[str, Any]:
    config = read_yaml(path, "CUP feature-evidence rule configuration")
    rules = config.get("rules")
    if not isinstance(rules, list):
        raise ValueError("CUP feature-evidence rule configuration lacks a rules list")
    rule_ids: list[str] = []
    for rule in rules:
        if not isinstance(rule, dict) or not rule.get("rule_id"):
            raise ValueError("Each CUP feature-evidence rule requires rule_id")
        rule_id = str(rule["rule_id"])
        rule_ids.append(rule_id)
        condition = rule.get("condition", {})
        if ("all" in condition) == ("any" in condition):
            raise ValueError(f"Evidence rule {rule_id} must define exactly one of condition.all or condition.any")
        unknown_features = set(condition_features(condition)) - available_features
        if unknown_features:
            raise ValueError(f"Evidence rule {rule_id} references unavailable features: {sorted(unknown_features)}")
        for target in list(rule.get("supports", [])) + list(rule.get("argues_against", [])):
            target = str(target)
            if target in lineage_names or target.startswith("major_group:") or target.endswith("_phenotype"):
                continue
            raise ValueError(f"Evidence rule {rule_id} has unknown target {target!r}")
        if str(rule.get("strength", "")) not in {"weak", "moderate", "strong"}:
            raise ValueError(f"Evidence rule {rule_id} has unsupported strength")
    if len(rule_ids) != len(set(rule_ids)):
        raise ValueError("CUP feature-evidence rules contain duplicate rule_id values")
    return config


def recompute_project_summary(neighbors: pd.DataFrame) -> pd.DataFrame:
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
    summary["weighted_similarity_score"] = (
        summary["similarity_weight"] / total_weight if total_weight > 0 else 0.0
    )
    summary = summary.sort_values(
        ["weighted_similarity_score", "min_distance", "project_code"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    summary["rank"] = np.arange(1, len(summary) + 1)
    summary["interpretation_note"] = (
        "Normalized share of top-k molecular-similarity weight; not a tissue-of-origin probability"
    )
    return summary[PROJECT_SUMMARY_REQUIRED_COLUMNS + ["interpretation_note"]]


def recompute_major_group_summary(neighbors: pd.DataFrame) -> pd.DataFrame:
    total_weight = float(neighbors["similarity_score"].sum())
    summary = neighbors.groupby("tcga_major_cancer_group", as_index=False).agg(
        n_neighbors=("tcga_sample_barcode", "size"),
        median_distance=("distance", "median"),
        min_distance=("distance", "min"),
        similarity_weight=("similarity_score", "sum"),
    )
    summary = summary.rename(columns={"tcga_major_cancer_group": "major_cancer_group"})
    summary["weighted_similarity_score"] = (
        summary["similarity_weight"] / total_weight if total_weight > 0 else 0.0
    )
    summary = summary.sort_values(
        ["weighted_similarity_score", "min_distance", "major_cancer_group"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    summary["rank"] = np.arange(1, len(summary) + 1)
    summary["interpretation_note"] = (
        "Normalized share of top-k molecular-similarity weight; not a lineage probability"
    )
    return summary[MAJOR_GROUP_SUMMARY_REQUIRED_COLUMNS + ["interpretation_note"]]


def resolve_similarity_summaries(inputs: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    neighbors = inputs["neighbors"]
    project_recomputed = recompute_project_summary(neighbors)
    major_recomputed = recompute_major_group_summary(neighbors)
    using_all_neighbors = len(neighbors) == len(inputs["neighbors_all"])
    if not using_all_neighbors:
        return project_recomputed, major_recomputed, "recomputed_from_requested_top_k"

    project_provided = inputs["provided_project_summary"].copy()
    project_comparison = project_recomputed.merge(
        project_provided[["project_code", "n_neighbors", "weighted_similarity_score"]],
        on="project_code",
        how="outer",
        suffixes=("_calculated", "_provided"),
        indicator=True,
    )
    if not (project_comparison["_merge"] == "both").all():
        raise ValueError("Provided project summary does not contain the same projects as the neighbor table")
    if not np.array_equal(
        project_comparison["n_neighbors_calculated"].to_numpy(dtype=int),
        project_comparison["n_neighbors_provided"].to_numpy(dtype=int),
    ) or not np.allclose(
        project_comparison["weighted_similarity_score_calculated"].to_numpy(dtype=float),
        project_comparison["weighted_similarity_score_provided"].to_numpy(dtype=float),
        atol=1e-8,
        rtol=1e-8,
    ):
        raise ValueError("Provided project summary is inconsistent with the nearest-neighbor table")

    major_provided = inputs["provided_major_summary"].copy()
    major_comparison = major_recomputed.merge(
        major_provided[["major_cancer_group", "n_neighbors", "weighted_similarity_score"]],
        on="major_cancer_group",
        how="outer",
        suffixes=("_calculated", "_provided"),
        indicator=True,
    )
    if not (major_comparison["_merge"] == "both").all():
        raise ValueError("Provided major-group summary is inconsistent with the neighbor table")
    if not np.array_equal(
        major_comparison["n_neighbors_calculated"].to_numpy(dtype=int),
        major_comparison["n_neighbors_provided"].to_numpy(dtype=int),
    ) or not np.allclose(
        major_comparison["weighted_similarity_score_calculated"].to_numpy(dtype=float),
        major_comparison["weighted_similarity_score_provided"].to_numpy(dtype=float),
        atol=1e-8,
        rtol=1e-8,
    ):
        raise ValueError("Provided major-group summary is inconsistent with the nearest-neighbor table")

    project_provided = project_provided.sort_values("rank", kind="stable").reset_index(drop=True)
    major_provided = major_provided.sort_values("rank", kind="stable").reset_index(drop=True)
    return project_provided, major_provided, "validated_projection_summary"


def feature_coverage(global_qc: dict[str, Any], feature_count: int) -> float:
    if "percent_feature_coverage" in global_qc:
        coverage = float(global_qc["percent_feature_coverage"]) / 100.0
    else:
        supplied = float(global_qc.get("number_of_clinical_features_supplied", 0))
        coverage = supplied / feature_count if feature_count else 0.0
    return float(np.clip(coverage, 0.0, 1.0))


def modality_feature_coverage(inputs: dict[str, Any], modality: str) -> float:
    qc = inputs["projection_qc"]
    modality_qc = (
        qc[qc["qc_section"].astype(str).str.upper() == modality.upper()]
        .drop_duplicates("metric", keep="last")
        .set_index("metric")["value"]
        .to_dict()
    )
    if {
        "number_of_expected_reference_features",
        "number_of_clinical_features_supplied",
    }.issubset(modality_qc):
        expected = float(modality_qc["number_of_expected_reference_features"])
        supplied = float(modality_qc["number_of_clinical_features_supplied"])
        if expected > 0:
            return float(np.clip(supplied / expected, 0.0, 1.0))

    modality_features = {
        feature_name(record)
        for record in inputs["feature_records"]
        if str(record.get("clinical_modality", "")).upper() in {modality.upper(), "BOTH"}
    }
    if not modality_features:
        return 0.0
    if inputs["feature_availability_audited"]:
        return float(
            np.clip(
                len(modality_features & inputs["supplied_features"]) / len(modality_features),
                0.0,
                1.0,
            )
        )
    return feature_coverage(inputs["global_qc"], len(inputs["feature_order"]))


def initial_lineage_scores(
    lineages: dict[str, Any],
    project_summary: pd.DataFrame,
    neighbors: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    project_to_lineage = lineages["project_to_lineage"]
    neighbor_lineages = neighbors["tcga_project_code"].astype(str).map(project_to_lineage)
    for order, record in enumerate(lineages["lineages"]):
        subset = project_summary[project_summary["project_code"].isin(record["projects"])].copy()
        subset = subset.sort_values(
            ["weighted_similarity_score", "project_code"], ascending=[False, True], kind="stable"
        )
        supporting_projects = "; ".join(
            f"{row.project_code} ({float(row.weighted_similarity_score):.4f})"
            for row in subset.itertuples(index=False)
            if float(row.weighted_similarity_score) > 0
        )
        rows.append(
            {
                "config_order": order,
                "lineage_or_project_group": record["name"],
                "display_name": record["display_name"],
                "supporting_tcga_projects": supporting_projects,
                "similarity_score": float(subset["weighted_similarity_score"].sum()),
                "n_neighbors": int((neighbor_lineages == record["name"]).sum()),
                "description": record["description"],
                "caveats": record["caveat"],
                "major_cancer_groups": record["major_cancer_groups"],
            }
        )
    scores = pd.DataFrame(rows).sort_values(
        ["similarity_score", "config_order"], ascending=[False, True], kind="stable"
    ).reset_index(drop=True)
    scores["rank"] = np.arange(1, len(scores) + 1)
    score_sum = float(scores["similarity_score"].sum())
    if not np.isclose(score_sum, 1.0, atol=1e-8):
        raise ValueError(f"Lineage scores must sum to one after complete project mapping; observed {score_sum}")
    return scores


def normalized_entropy(scores: np.ndarray) -> float:
    positive = scores[np.isfinite(scores) & (scores > 0)]
    if len(positive) <= 1:
        return 0.0
    positive = positive / positive.sum()
    return float(-np.sum(positive * np.log(positive)) / math.log(len(positive)))


def build_ambiguity_metrics(
    case_id: str,
    sample_id: str,
    lineage_scores: pd.DataFrame,
    project_summary: pd.DataFrame,
    major_summary: pd.DataFrame,
    coverage: float,
    wes_coverage: float,
    wts_coverage: float,
    top_k: int,
    config: dict[str, Any],
) -> pd.DataFrame:
    top1 = lineage_scores.iloc[0]
    top2 = lineage_scores.iloc[1] if len(lineage_scores) > 1 else None
    top1_score = float(top1["similarity_score"])
    top2_score = float(top2["similarity_score"]) if top2 is not None else 0.0
    margin = top1_score - top2_score
    entropy = normalized_entropy(lineage_scores["similarity_score"].to_numpy(dtype=float))
    meaningful_projects = int(
        (project_summary["weighted_similarity_score"] >= float(config["meaningful_project_score_min"])).sum()
    )
    top_major_group = str(major_summary.iloc[0]["major_cancer_group"])
    top_major_group_score = float(major_summary.iloc[0]["weighted_similarity_score"])
    agreement = top_major_group in set(top1["major_cancer_groups"])
    neighbor_concentration = float(top1["n_neighbors"]) / top_k if top_k else 0.0

    if coverage < float(config["uninterpretable_coverage_below"]) or top1_score <= 0:
        ambiguity_class = "uninterpretable_due_to_missing_features"
    elif (
        top1_score >= float(config["low_ambiguity_top1_min"])
        and margin >= float(config["low_ambiguity_margin_min"])
        and entropy <= float(config["low_ambiguity_entropy_max"])
        and agreement
    ):
        ambiguity_class = "low_ambiguity"
    elif (
        top1_score < float(config["high_ambiguity_top1_below"])
        or margin < float(config["high_ambiguity_margin_below"])
        or entropy >= float(config["high_ambiguity_entropy_min"])
        or meaningful_projects >= int(config["high_ambiguity_meaningful_projects_min"])
    ):
        ambiguity_class = "high_ambiguity"
    else:
        ambiguity_class = "moderate_ambiguity"

    note = (
        "Scores are normalized shares of top-k neighbor similarity weight, not diagnostic probabilities. "
        "Entropy is normalized over lineages with nonzero neighbor weight."
    )
    return pd.DataFrame(
        [
            {
                "case_id": case_id,
                "query_sample_id": sample_id,
                "top1_lineage": str(top1["lineage_or_project_group"]),
                "top2_lineage": str(top2["lineage_or_project_group"]) if top2 is not None else "NA",
                "top1_score": top1_score,
                "top2_score": top2_score,
                "top1_top2_margin": margin,
                "entropy_like_score": entropy,
                "number_of_projects_with_meaningful_similarity": meaningful_projects,
                "feature_coverage": coverage,
                "wes_feature_coverage": wes_coverage,
                "wts_feature_coverage": wts_coverage,
                "top_major_cancer_group": top_major_group,
                "top_major_group_score": top_major_group_score,
                "major_group_agreement": agreement,
                "neighbor_concentration_top_lineage": neighbor_concentration,
                "ambiguity_class": ambiguity_class,
                "interpretation_note": note,
            }
        ]
    )


def confidence_tier(
    rank: int,
    score: float,
    ambiguity: pd.Series,
    config: dict[str, Any],
) -> str:
    if str(ambiguity["ambiguity_class"]) == "uninterpretable_due_to_missing_features":
        return "indeterminate"
    if rank == 1:
        if (
            score >= float(config["high_support_top1_min"])
            and float(ambiguity["top1_top2_margin"]) >= float(config["high_support_margin_min"])
            and float(ambiguity["feature_coverage"]) >= float(config["high_support_coverage_min"])
            and float(ambiguity["neighbor_concentration_top_lineage"])
            >= float(config["high_support_neighbor_concentration_min"])
            and bool(ambiguity["major_group_agreement"])
        ):
            return "high_support"
        if (
            score >= float(config["moderate_support_top1_min"])
            and float(ambiguity["top1_top2_margin"]) >= float(config["moderate_support_margin_min"])
            and float(ambiguity["feature_coverage"]) >= float(config["moderate_support_coverage_min"])
            and float(ambiguity["neighbor_concentration_top_lineage"])
            >= float(config["moderate_support_neighbor_concentration_min"])
            and bool(ambiguity["major_group_agreement"])
        ):
            return "moderate_support"
        if score >= 0.20 and float(ambiguity["feature_coverage"]) >= 0.50:
            return "weak_support"
        return "indeterminate"
    if score >= 0.10:
        return "weak_support"
    return "indeterminate"


def compare_value(value: float, operator: str, threshold: float) -> bool:
    operators = {
        "ge": lambda: value >= threshold,
        "gt": lambda: value > threshold,
        "le": lambda: value <= threshold,
        "lt": lambda: value < threshold,
        "eq": lambda: np.isclose(value, threshold),
    }
    if operator not in operators:
        raise ValueError(f"Unsupported evidence-rule operator: {operator}")
    return bool(operators[operator]())


def evaluate_rule(
    rule: dict[str, Any],
    raw_values: dict[str, float],
    scaled_values: dict[str, float],
    supplied_features: set[str],
    availability_audited: bool,
) -> dict[str, Any]:
    condition = rule["condition"]
    mode = "all" if "all" in condition else "any"
    results: list[bool | None] = []
    details: list[str] = []
    for clause in condition[mode]:
        feature = str(clause["feature"])
        source = str(clause.get("value_source", "raw"))
        if availability_audited and feature not in supplied_features:
            results.append(None)
            details.append(f"{feature}=not supplied")
            continue
        if not availability_audited:
            results.append(None)
            details.append(f"{feature}=supply status unavailable")
            continue
        values = raw_values if source == "raw" else scaled_values
        value = float(values[feature])
        threshold = float(clause["threshold"])
        operator = str(clause["operator"])
        results.append(compare_value(value, operator, threshold))
        details.append(f"{feature}={value:.4g} ({source}, {operator} {threshold:.4g})")
    if mode == "all":
        evaluable = all(result is not None for result in results)
        triggered = evaluable and all(bool(result) for result in results)
    else:
        if any(result is True for result in results):
            evaluable = True
            triggered = True
        elif any(result is None for result in results):
            evaluable = False
            triggered = False
        else:
            evaluable = True
            triggered = False
    return {
        **rule,
        "evaluable": evaluable,
        "triggered": triggered,
        "condition_details": "; ".join(details),
    }


def evaluate_rules(evidence_config: dict[str, Any], inputs: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        evaluate_rule(
            rule,
            inputs["raw_values"],
            inputs["scaled_values"],
            inputs["supplied_features"],
            inputs["feature_availability_audited"],
        )
        for rule in evidence_config.get("rules", [])
    ]


def finish_ranked_lineages(
    lineage_scores: pd.DataFrame,
    ambiguity_metrics: pd.DataFrame,
    evaluated_rules: list[dict[str, Any]],
    config: dict[str, Any],
) -> pd.DataFrame:
    ambiguity = ambiguity_metrics.iloc[0]
    triggered = [rule for rule in evaluated_rules if rule["triggered"]]
    evidence_by_lineage: dict[str, list[str]] = {}
    for rule in triggered:
        for target in rule.get("supports", []):
            target = str(target)
            evidence_by_lineage.setdefault(target, []).append(str(rule["finding"]))
    output = lineage_scores.copy()
    output["confidence_tier"] = [
        confidence_tier(int(row.rank), float(row.similarity_score), ambiguity, config)
        for row in output.itertuples(index=False)
    ]
    summaries: list[str] = []
    caveats: list[str] = []
    for row in output.itertuples(index=False):
        parts: list[str] = []
        if float(row.similarity_score) > 0:
            parts.append(
                f"Top-k neighbor weight {float(row.similarity_score):.4f} from "
                f"{row.supporting_tcga_projects}."
            )
        else:
            parts.append("No returned top-k neighbors mapped to this lineage group.")
        rule_findings = evidence_by_lineage.get(str(row.lineage_or_project_group), [])
        if rule_findings:
            parts.append("Heuristic feature observations: " + " ".join(rule_findings))
        else:
            parts.append("No lineage-specific heuristic feature rule triggered.")
        summaries.append(" ".join(parts))
        caveats.append(
            f"{row.caveats} Similarity scores are normalized descriptive weights, not probabilities or diagnoses."
        )
    output["evidence_summary"] = summaries
    output["caveats"] = caveats
    return output[RANKED_LINEAGE_COLUMNS]


def tcga_percentile(reference_matrix: pd.DataFrame, feature: str, value: float) -> float | None:
    reference = pd.to_numeric(reference_matrix[feature], errors="coerce").dropna().to_numpy(dtype=float)
    reference = reference[np.isfinite(reference)]
    if len(reference) == 0:
        return None
    return float(100.0 * np.mean(reference <= value))


def feature_summary(
    features: list[str],
    inputs: dict[str, Any],
    max_features: int = 5,
) -> str:
    candidates = [
        feature
        for feature in features
        if feature in inputs["supplied_features"] and feature in inputs["scaled_values"]
    ]
    candidates = sorted(candidates, key=lambda feature: abs(inputs["scaled_values"][feature]), reverse=True)
    entries: list[str] = []
    for feature in candidates[:max_features]:
        raw = float(inputs["raw_values"][feature])
        scaled = float(inputs["scaled_values"][feature])
        percentile = tcga_percentile(inputs["reference_matrix"], feature, raw)
        percentile_text = f", TCGA percentile={percentile:.1f}" if percentile is not None else ""
        entries.append(f"{feature}={raw:.4g} (z={scaled:.2f}{percentile_text})")
    return "; ".join(entries) if entries else "No supplied features were available for this evidence category."


def maximum_strength(rules: list[dict[str, Any]], default: str = "not assessed") -> str:
    strengths = [str(rule["strength"]) for rule in rules if rule.get("triggered")]
    if not strengths:
        return default
    return max(strengths, key=lambda value: EVIDENCE_STRENGTH_ORDER[value])


def joined_targets(rules: list[dict[str, Any]], field: str) -> str:
    targets = sorted({str(target) for rule in rules if rule.get("triggered") for target in rule.get(field, [])})
    return "; ".join(targets) if targets else "None"


def build_molecular_evidence(
    inputs: dict[str, Any],
    project_summary: pd.DataFrame,
    major_summary: pd.DataFrame,
    ranked_lineages: pd.DataFrame,
    ambiguity_metrics: pd.DataFrame,
    evaluated_rules: list[dict[str, Any]],
    evidence_config: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    top_projects = "; ".join(
        f"{row.project_code}={float(row.weighted_similarity_score):.4f} ({int(row.n_neighbors)} neighbors)"
        for row in project_summary.head(6).itertuples(index=False)
    )
    top_lineage = ranked_lineages.iloc[0]
    project_strength = {
        "high_support": "strong",
        "moderate_support": "moderate",
        "weak_support": "weak",
        "indeterminate": "not assessed",
    }[str(top_lineage["confidence_tier"])]
    rows.append(
        {
            "evidence_type": EVIDENCE_TYPES[0],
            "finding": f"Top-k project weights: {top_projects}.",
            "supports": str(top_lineage["lineage_or_project_group"]),
            "argues_against": "None",
            "strength": project_strength,
            "notes": "Project weights are the primary molecular-similarity evidence and are not tissue-of-origin probabilities.",
        }
    )

    top_major = major_summary.iloc[0]
    rows.append(
        {
            "evidence_type": EVIDENCE_TYPES[1],
            "finding": (
                f"{top_major['major_cancer_group']} accounts for "
                f"{float(top_major['weighted_similarity_score']):.4f} of top-k similarity weight "
                f"({int(top_major['n_neighbors'])} neighbors)."
            ),
            "supports": str(top_major["major_cancer_group"]),
            "argues_against": "None",
            "strength": "moderate" if float(top_major["weighted_similarity_score"]) >= 0.50 else "weak",
            "notes": "Major-group evidence is derived from the same neighbors as project evidence and is not independent validation.",
        }
    )

    category_features = {
        EVIDENCE_TYPES[2]: [
            "proliferation_score",
            "epithelial_score",
            "hypoxia_score",
            "cell_cycle_score",
            "DNA_repair_score",
            "angiogenesis_score",
        ],
        EVIDENCE_TYPES[4]: [
            "ploidy",
            "aneuploidy_score",
            "arm_gain_count",
            "arm_loss_count",
            "total_arm_alteration_count",
        ],
        EVIDENCE_TYPES[5]: [
            "immune_inflammatory_score",
            "interferon_gamma_score",
            "cytotoxic_t_cell_score",
            "stromal_score",
            "EMT_score",
        ],
    }
    triggered_by_type = {
        evidence_type: [
            rule for rule in evaluated_rules if rule.get("triggered") and rule["evidence_type"] == evidence_type
        ]
        for evidence_type in EVIDENCE_TYPES
    }
    for evidence_type in [EVIDENCE_TYPES[2], EVIDENCE_TYPES[4], EVIDENCE_TYPES[5]]:
        hits = triggered_by_type[evidence_type]
        hit_text = " ".join(str(rule["finding"]) for rule in hits)
        summary = feature_summary(category_features[evidence_type], inputs)
        finding = f"{hit_text} Observed context: {summary}" if hit_text else f"Observed context: {summary}"
        rows.append(
            {
                "evidence_type": evidence_type,
                "finding": finding,
                "supports": joined_targets(hits, "supports"),
                "argues_against": joined_targets(hits, "argues_against"),
                "strength": maximum_strength(hits, default="weak" if "No supplied" not in summary else "not assessed"),
                "notes": (
                    "Pathway scores are broad phenotype summaries, not validated lineage classifiers."
                    if evidence_type == EVIDENCE_TYPES[2]
                    else "Broad copy-number burden is not tissue-specific."
                    if evidence_type == EVIDENCE_TYPES[4]
                    else "Immune, stromal, and EMT signals may reflect microenvironment as well as tumor cells."
                ),
            }
        )

    driver_hits = triggered_by_type[EVIDENCE_TYPES[3]]
    driver_features = [
        feature
        for feature in inputs["feature_order"]
        if feature.startswith("driver_mutation_")
        and not feature.startswith("driver_mutation_count_")
        and feature in inputs["supplied_features"]
        and float(inputs["raw_values"][feature]) >= 1
    ]
    driver_names = [feature.removeprefix("driver_mutation_") for feature in driver_features]
    driver_finding = (
        "Observed configured driver indicators: " + ", ".join(driver_names) + "."
        if driver_names
        else "No positive supplied configured driver indicators were observed."
    )
    if driver_hits:
        driver_finding += " Triggered heuristic observations: " + " ".join(
            str(rule["finding"]) for rule in driver_hits
        )
    rows.insert(
        3,
        {
            "evidence_type": EVIDENCE_TYPES[3],
            "finding": driver_finding,
            "supports": joined_targets(driver_hits, "supports"),
            "argues_against": joined_targets(driver_hits, "argues_against"),
            "strength": maximum_strength(driver_hits, default="not assessed"),
            "notes": "Binary/count driver features omit variant pathogenicity, clonality, signatures, and many lineage-relevant genes.",
        },
    )

    moderate_lineage_targets = sorted(
        {
            str(target)
            for rule in evaluated_rules
            if rule.get("triggered") and str(rule.get("strength")) in {"moderate", "strong"}
            for target in rule.get("supports", [])
            if str(target) in set(ranked_lineages["lineage_or_project_group"])
        }
    )
    competing_lineages = ranked_lineages[ranked_lineages["similarity_score"] >= 0.15][
        "lineage_or_project_group"
    ].astype(str).tolist()
    contradiction_parts: list[str] = []
    if len(competing_lineages) > 1:
        contradiction_parts.append(
            "Substantial top-k neighbor weight is split across " + ", ".join(competing_lineages) + "."
        )
    if len(moderate_lineage_targets) > 1:
        contradiction_parts.append(
            "Moderate heuristic feature patterns support more than one lineage: "
            + ", ".join(moderate_lineage_targets)
            + "."
        )
    contradiction = " ".join(contradiction_parts) or "No explicit cross-lineage contradiction was identified."
    rows.append(
        {
            "evidence_type": EVIDENCE_TYPES[6],
            "finding": contradiction,
            "supports": "Maintaining an appropriately broad differential",
            "argues_against": "Overconfident single-lineage interpretation",
            "strength": "moderate" if contradiction_parts else "not assessed",
            "notes": "Competing molecular patterns may reflect shared biology, mixed evidence, or incomplete assay harmonization.",
        }
    )

    missing_count = int(float(inputs["global_qc"].get("number_of_missing_clinical_features", 0)))
    if inputs["missing_features"]:
        feature_classes = {
            feature_name(record): str(record.get("feature_class", record.get("source_layer", "unclassified")))
            for record in inputs["feature_records"]
        }
        grouped: dict[str, int] = {}
        for feature in inputs["missing_features"]:
            group = feature_classes.get(feature, "unclassified")
            grouped[group] = grouped.get(group, 0) + 1
        missing_detail = ", ".join(f"{group}={count}" for group, count in sorted(grouped.items()))
    elif missing_count:
        missing_detail = "specific feature identities were not available from the optional feature audit"
    else:
        missing_detail = "none"
    unavailable = evidence_config.get("unavailable_evidence", [])
    unavailable_labels = "; ".join(str(item.get("label")) for item in unavailable)
    rows.append(
        {
            "evidence_type": EVIDENCE_TYPES[7],
            "finding": (
                f"{missing_count} locked query features were median-imputed ({missing_detail}). "
                f"Evidence not represented in the atlas includes: {unavailable_labels}."
            ),
            "supports": "Cautious interpretation and orthogonal review",
            "argues_against": "High-confidence tissue-of-origin assignment",
            "strength": "strong" if missing_count > 0 or unavailable else "not assessed",
            "notes": "Imputed features support numerical projection but are not observed molecular evidence.",
        }
    )

    evidence = pd.DataFrame(rows)
    evidence["evidence_type"] = pd.Categorical(evidence["evidence_type"], EVIDENCE_TYPES, ordered=True)
    evidence = evidence.sort_values("evidence_type").reset_index(drop=True)
    evidence["evidence_type"] = evidence["evidence_type"].astype(str)
    return evidence[["evidence_type", "finding", "supports", "argues_against", "strength", "notes"]]


def normalize_text(value: Any) -> str:
    if is_missing(value):
        return ""
    normalized = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    normalized = normalized.lower().replace("_", " ").replace("-", " ")
    normalized = re.sub(r"[^a-z0-9]+", " ", normalized)
    return re.sub(r"\s+", " ", normalized).strip()


def split_differential(value: Any) -> list[str]:
    if is_missing(value):
        return []
    raw = str(value).strip()
    pieces = re.split(r"\s+(?:vs\.?|versus|or)\s+|[,;/|]", raw, flags=re.IGNORECASE)
    stop_words = {"synthetic", "differential", "diagnosis", "possible", "consider", "consideration"}
    cleaned: list[str] = []
    for piece in pieces:
        words = [word for word in normalize_text(piece).split() if word not in stop_words]
        candidate = " ".join(words).strip()
        if candidate and candidate not in cleaned:
            cleaned.append(candidate)
    return cleaned


def match_candidate_lineages(
    candidate: str,
    lineages: dict[str, Any],
    reference_metadata: pd.DataFrame,
) -> list[str]:
    normalized_candidate = normalize_text(candidate)
    matches: set[str] = set()
    for record in lineages["lineages"]:
        aliases = [record["name"], record["display_name"], *record["aliases"]]
        if any(
            normalize_text(alias)
            and re.search(rf"\b{re.escape(normalize_text(alias))}\b", normalized_candidate)
            for alias in aliases
        ):
            matches.add(record["name"])
        for project in record["projects"]:
            if re.search(rf"\b{re.escape(normalize_text(project))}\b", normalized_candidate):
                matches.add(record["name"])
    major_groups = sorted(reference_metadata["major_cancer_group"].dropna().astype(str).unique())
    for major_group in major_groups:
        normalized_group = normalize_text(major_group)
        if normalized_group and re.search(rf"\b{re.escape(normalized_group)}\b", normalized_candidate):
            for record in lineages["lineages"]:
                if major_group in record["major_cancer_groups"]:
                    matches.add(record["name"])
    return sorted(matches, key=lambda name: list(lineages["by_name"]).index(name))


def build_differential_summary(
    clinical_metadata: dict[str, Any],
    ranked_lineages: pd.DataFrame,
    lineages: dict[str, Any],
    reference_metadata: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    differential = clinical_metadata.get("differential_diagnosis")
    candidates = split_differential(differential)
    columns = [
        "candidate_rank",
        "candidate_text",
        "matched_lineages_or_projects",
        "constrained_similarity_score",
        "best_supported_lineage",
        "overall_lineage_rank",
        "interpretation",
        "caveats",
    ]
    if not candidates:
        return pd.DataFrame(
            [
                {
                    "candidate_rank": pd.NA,
                    "candidate_text": "No differential diagnosis supplied",
                    "matched_lineages_or_projects": "None",
                    "constrained_similarity_score": pd.NA,
                    "best_supported_lineage": "Not assessed",
                    "overall_lineage_rank": pd.NA,
                    "interpretation": "No constrained differential analysis was performed.",
                    "caveats": "The unconstrained TCGA similarity interpretation remains research-only.",
                }
            ],
            columns=columns,
        )
    score_map = ranked_lineages.set_index("lineage_or_project_group")["similarity_score"].to_dict()
    rank_map = ranked_lineages.set_index("lineage_or_project_group")["rank"].to_dict()
    rows: list[dict[str, Any]] = []
    for original_order, candidate in enumerate(candidates):
        matched = match_candidate_lineages(candidate, lineages, reference_metadata)
        constrained_score = float(sum(float(score_map.get(name, 0.0)) for name in matched))
        best = max(matched, key=lambda name: float(score_map.get(name, 0.0))) if matched else None
        if not matched:
            interpretation = "Candidate text did not map to a configured TCGA project or lineage group."
        elif constrained_score >= float(config["differential_relative_support_min"]):
            interpretation = "Relatively supported within the supplied differential; this is not a diagnosis."
        elif constrained_score >= 0.20:
            interpretation = "Some relative molecular support within the supplied differential."
        elif constrained_score > 0:
            interpretation = "Limited relative molecular support within the supplied differential."
        else:
            interpretation = "No returned top-k neighbor weight supports the mapped candidate group."
        rows.append(
            {
                "original_order": original_order,
                "candidate_text": candidate,
                "matched_lineages_or_projects": "; ".join(matched) if matched else "None",
                "constrained_similarity_score": constrained_score,
                "best_supported_lineage": best or "Not assessed",
                "overall_lineage_rank": int(rank_map[best]) if best is not None else pd.NA,
                "interpretation": interpretation,
                "caveats": (
                    "Differential constraints highlight supplied candidates but do not change the unconstrained "
                    "molecular ranking or constitute clinical validation."
                ),
            }
        )
    output = pd.DataFrame(rows).sort_values(
        ["constrained_similarity_score", "original_order"], ascending=[False, True], kind="stable"
    ).reset_index(drop=True)
    output["candidate_rank"] = np.arange(1, len(output) + 1)
    return output[columns]


def molecular_pathologic_concordance(
    case_id: str,
    sample_id: str,
    submitted_diagnosis: Any,
    top_lineage: str,
    top_major_group: str,
    lineages: dict[str, Any],
    reference_metadata: pd.DataFrame,
) -> pd.DataFrame:
    normalized = normalize_text(submitted_diagnosis)
    if not normalized:
        concordance_class = "not_assessable"
        rationale = "No submitted diagnosis was supplied."
    elif "unknown primary" in normalized or re.search(r"\bcup\b", normalized):
        concordance_class = "not_assessable"
        rationale = "The submitted diagnosis does not specify a tissue lineage for comparison."
    else:
        matched = match_candidate_lineages(normalized, lineages, reference_metadata)
        if top_lineage in matched:
            concordance_class = "concordant"
            rationale = "Submitted diagnosis text heuristically maps to the top molecular lineage group."
        elif matched:
            concordance_class = "discordant"
            rationale = "Submitted diagnosis text maps to a different configured lineage than the top result."
        elif normalize_text(top_major_group) and normalize_text(top_major_group) in normalized:
            concordance_class = "partially_concordant"
            rationale = "Submitted diagnosis and molecular result share only a broad major cancer group."
        else:
            concordance_class = "not_assessable"
            rationale = "Submitted diagnosis text could not be mapped to the configured lineage vocabulary."
    return pd.DataFrame(
        [
            {
                "case_id": case_id,
                "query_sample_id": sample_id,
                "submitted_diagnosis": submitted_diagnosis,
                "top_lineage_or_project_group": top_lineage,
                "top_major_cancer_group": top_major_group,
                "concordance_class": concordance_class,
                "rationale": rationale,
                "caveats": (
                    "Heuristic text comparison only; definitive concordance requires morphology, "
                    "immunophenotype, imaging, specimen context, and expert pathologic review."
                ),
            }
        ]
    )


def markdown_escape(value: Any) -> str:
    if is_missing(value):
        return "Not supplied"
    return str(value).replace("|", "\\|").replace("\n", " ")


def markdown_table(data: pd.DataFrame, columns: list[tuple[str, str]], max_rows: int | None = None) -> str:
    shown = data.head(max_rows) if max_rows is not None else data
    headers = [label for _, label in columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in shown.iterrows():
        values: list[str] = []
        for column, _label in columns:
            value = row[column]
            if isinstance(value, (float, np.floating)) and np.isfinite(value):
                values.append(f"{float(value):.4f}")
            else:
                values.append(markdown_escape(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def plot_project_similarity(project_summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    plotted = project_summary.head(12).sort_values("weighted_similarity_score", ascending=True)
    colors = ["#B04A5A" if int(rank) == 1 else "#3D6F8E" for rank in plotted["rank"]]
    fig, ax = plt.subplots(figsize=(10, max(5.5, 0.55 * len(plotted))))
    bars = ax.barh(plotted["project_code"], plotted["weighted_similarity_score"], color=colors)
    ax.bar_label(
        bars,
        labels=[f"{value:.3f}" for value in plotted["weighted_similarity_score"]],
        padding=3,
        fontsize=8,
    )
    ax.set_xlabel("Normalized share of top-k similarity weight")
    ax.set_ylabel("")
    ax.set_title("TCGA project composition of nearest neighbors")
    ax.set_xlim(0, max(0.10, float(plotted["weighted_similarity_score"].max()) * 1.22))
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.text(
        0,
        -0.14,
        "Descriptive molecular similarity only; bars are not tissue-of-origin probabilities.",
        transform=ax.transAxes,
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_major_group_similarity(major_summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    plotted = major_summary.head(12).sort_values("weighted_similarity_score", ascending=True)
    fig, ax = plt.subplots(figsize=(10, max(4.8, 0.58 * len(plotted))))
    bars = ax.barh(plotted["major_cancer_group"], plotted["weighted_similarity_score"], color="#4D7A55")
    ax.bar_label(
        bars,
        labels=[f"{value:.3f}" for value in plotted["weighted_similarity_score"]],
        padding=3,
        fontsize=8,
    )
    ax.set_xlabel("Normalized share of top-k similarity weight")
    ax.set_ylabel("")
    ax.set_title("Major cancer-group composition of nearest neighbors")
    ax.set_xlim(0, max(0.10, float(plotted["weighted_similarity_score"].max()) * 1.15))
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.text(
        0,
        -0.16,
        "Major-group weights reuse the same neighbors and are not independent validation.",
        transform=ax.transAxes,
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_ambiguity(ambiguity_metrics: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    row = ambiguity_metrics.iloc[0]
    labels = ["Top lineage", "Top1-top2 margin", "Normalized entropy", "Feature coverage"]
    values = [
        float(row["top1_score"]),
        float(row["top1_top2_margin"]),
        float(row["entropy_like_score"]),
        float(row["feature_coverage"]),
    ]
    colors = ["#3D6F8E", "#4D7A55", "#B06A3C", "#7A5A8C"]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    bars = ax.bar(labels, values, color=colors, width=0.66)
    ax.bar_label(bars, labels=[f"{value:.3f}" for value in values], padding=4, fontsize=9)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Score on 0-1 scale")
    ambiguity_label = str(row["ambiguity_class"]).replace("_", " ")
    ax.set_title(f"CUP ambiguity assessment: {ambiguity_label}")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.text(
        0,
        -0.17,
        "Entropy is normalized over represented lineages; these metrics are heuristic and not clinically calibrated.",
        transform=ax.transAxes,
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def evidence_matrix(
    ranked_lineages: pd.DataFrame,
    major_summary: pd.DataFrame,
    lineages: dict[str, Any],
    evaluated_rules: list[dict[str, Any]],
    strength_weights: dict[str, float],
) -> tuple[list[str], list[str], np.ndarray]:
    represented = ranked_lineages[ranked_lineages["similarity_score"] > 0].head(8)
    if represented.empty:
        represented = ranked_lineages.head(5)
    columns = represented["lineage_or_project_group"].astype(str).tolist()
    rows = EVIDENCE_TYPES[:6]
    matrix = np.zeros((len(rows), len(columns)), dtype=float)
    top_score = max(float(represented["similarity_score"].max()), 1e-12)
    score_map = represented.set_index("lineage_or_project_group")["similarity_score"].to_dict()
    for column_index, lineage in enumerate(columns):
        matrix[0, column_index] = float(score_map[lineage]) / top_score
        expected_groups = set(lineages["by_name"][lineage]["major_cancer_groups"])
        matrix[1, column_index] = float(
            major_summary.loc[
                major_summary["major_cancer_group"].isin(expected_groups), "weighted_similarity_score"
            ].sum()
        )
    row_lookup = {name: index for index, name in enumerate(rows)}
    for rule in evaluated_rules:
        if not rule.get("triggered") or rule.get("evidence_type") not in row_lookup:
            continue
        weight = float(strength_weights.get(str(rule["strength"]), 0.0))
        row_index = row_lookup[str(rule["evidence_type"])]
        for target in rule.get("supports", []):
            if target in columns:
                matrix[row_index, columns.index(target)] += weight
        for target in rule.get("argues_against", []):
            if target in columns:
                matrix[row_index, columns.index(target)] -= weight
    matrix = np.clip(matrix, -1.0, 1.0)
    return rows, columns, matrix


def plot_evidence_heatmap(
    ranked_lineages: pd.DataFrame,
    major_summary: pd.DataFrame,
    lineages: dict[str, Any],
    evaluated_rules: list[dict[str, Any]],
    evidence_config: dict[str, Any],
    path: Path,
) -> None:
    ensure_dir(path.parent)
    rows, columns, matrix = evidence_matrix(
        ranked_lineages,
        major_summary,
        lineages,
        evaluated_rules,
        {key: float(value) for key, value in evidence_config.get("strength_weights", {}).items()},
    )
    fig, ax = plt.subplots(figsize=(max(9, 1.5 * len(columns)), 7))
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    x_labels = [value.replace("_", "\n") for value in columns]
    ax.set_xticks(np.arange(len(columns)), x_labels, rotation=0, ha="center")
    ax.set_yticks(np.arange(len(rows)), rows)
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix[row_index, column_index]
            ax.text(
                column_index,
                row_index,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if abs(value) >= 0.60 else "#222222",
            )
    ax.set_title("Heuristic evidence map across represented lineage groups")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03)
    colorbar.set_label("Relative heuristic support (-1 to +1)")
    fig.text(
        0.5,
        0.02,
        "Nearest-neighbor and major-group rows reuse the same top-k samples. Rule rows are conservative heuristics, not probabilities.",
        ha="center",
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def build_report(
    case_id: str,
    inputs: dict[str, Any],
    project_summary: pd.DataFrame,
    major_summary: pd.DataFrame,
    ranked_lineages: pd.DataFrame,
    ambiguity_metrics: pd.DataFrame,
    evidence: pd.DataFrame,
    differential: pd.DataFrame,
    summary_source: str,
    concordance: pd.DataFrame,
) -> str:
    metadata = inputs["clinical_metadata"]
    ambiguity = ambiguity_metrics.iloc[0]
    top_lineage = ranked_lineages.iloc[0]
    coverage_percent = 100.0 * float(ambiguity["feature_coverage"])
    expected_features = int(float(inputs["global_qc"].get("number_of_expected_reference_features", len(inputs["feature_order"]))))
    supplied_features = int(float(inputs["global_qc"].get("number_of_clinical_features_supplied", 0)))
    missing_features = int(float(inputs["global_qc"].get("number_of_missing_clinical_features", 0)))
    project_table = markdown_table(
        project_summary,
        [
            ("rank", "Rank"),
            ("project_code", "TCGA project"),
            ("n_neighbors", "Neighbors"),
            ("median_distance", "Median distance"),
            ("weighted_similarity_score", "Similarity-weight share"),
        ],
        max_rows=10,
    )
    lineage_table = markdown_table(
        ranked_lineages[ranked_lineages["similarity_score"] > 0],
        [
            ("rank", "Rank"),
            ("lineage_or_project_group", "Lineage group"),
            ("supporting_tcga_projects", "Supporting TCGA projects"),
            ("similarity_score", "Similarity-weight share"),
            ("confidence_tier", "Confidence tier"),
        ],
        max_rows=10,
    )
    major_table = markdown_table(
        major_summary,
        [
            ("rank", "Rank"),
            ("major_cancer_group", "Major group"),
            ("n_neighbors", "Neighbors"),
            ("weighted_similarity_score", "Similarity-weight share"),
        ],
        max_rows=10,
    )
    ambiguity_table = markdown_table(
        ambiguity_metrics,
        [
            ("top1_score", "Top-1 score"),
            ("top2_score", "Top-2 score"),
            ("top1_top2_margin", "Margin"),
            ("entropy_like_score", "Normalized entropy"),
            ("number_of_projects_with_meaningful_similarity", "Meaningful projects"),
            ("feature_coverage", "Feature coverage"),
            ("wes_feature_coverage", "WES coverage"),
            ("wts_feature_coverage", "WTS coverage"),
            ("ambiguity_class", "Ambiguity class"),
        ],
    )
    evidence_table = markdown_table(
        evidence,
        [
            ("evidence_type", "Evidence type"),
            ("finding", "Finding"),
            ("supports", "Supports"),
            ("argues_against", "Argues against"),
            ("strength", "Strength"),
            ("notes", "Notes"),
        ],
    )
    differential_table = markdown_table(
        differential,
        [
            ("candidate_rank", "Rank"),
            ("candidate_text", "Candidate"),
            ("matched_lineages_or_projects", "Matched groups"),
            ("constrained_similarity_score", "Constrained score"),
            ("best_supported_lineage", "Best matched lineage"),
            ("interpretation", "Interpretation"),
        ],
    )
    concordance_row = concordance.iloc[0]
    return f"""# CUP Molecular Interpretation Report: {case_id}

Generated by `scripts/22_cup_nearest_neighbor_classifier.py` from the locked TCGA projection outputs.

## Case Metadata

- Case ID: `{markdown_escape(case_id)}`
- Sample ID: `{markdown_escape(inputs['sample_id'])}`
- Specimen type: {markdown_escape(metadata.get('specimen_type'))}
- Submitted diagnosis: {markdown_escape(metadata.get('submitted_diagnosis'))}
- Differential diagnosis: {markdown_escape(metadata.get('differential_diagnosis'))}
- Tumor purity metadata: {markdown_escape(metadata.get('tumor_purity'))}
- Notes: {markdown_escape(metadata.get('notes'))}

## Feature Completeness

- Expected locked reference features: {expected_features}
- Clinical features supplied: {supplied_features}
- Median-imputed features: {missing_features}
- Feature coverage: {coverage_percent:.3f}%
- WES feature coverage: {100.0 * float(ambiguity['wes_feature_coverage']):.3f}%
- WTS feature coverage: {100.0 * float(ambiguity['wts_feature_coverage']):.3f}%

Median-imputed values permit numerical projection but are not treated as observed molecular evidence by the heuristic rules. Feature completeness, assay harmonization, and the TCGA reference composition directly affect this interpretation.

## Top TCGA-Like Projects

{project_table}

The project summary source was `{summary_source}`. Similarity-weight shares are normalized descriptions of the returned top-k neighbors, not tissue-of-origin probabilities.

## Top Lineage/Major Group Interpretation

Top lineage group: `{top_lineage['lineage_or_project_group']}` with a similarity-weight share of {float(top_lineage['similarity_score']):.4f} and confidence tier `{top_lineage['confidence_tier']}`.

{lineage_table}

{major_table}

Lineage grouping makes project-level neighbor composition easier to review but does not convert the unsupervised search into a validated classifier.

## Ambiguity/Confidence Assessment

{ambiguity_table}

Confidence tiers use transparent thresholds based on top-lineage weight, top-1/top-2 margin, neighbor concentration, feature coverage, and agreement with the broad major-group summary. They have not been clinically calibrated.

## Molecular Evidence Table

{evidence_table}

Feature rules are conservative heuristics. Untriggered rules and absent mutations are not interpreted as evidence against a lineage.

## Differential Diagnosis Constrained Interpretation

{differential_table}

The constrained analysis highlights molecular support among supplied candidates. It does not change the unconstrained ranking and cannot exclude diagnoses not represented by TCGA or the locked feature set.

## Actionability Placeholder

Actionability integration is planned but not implemented in this prototype. No OncoKB, CIViC, therapeutic-response, or trial-matching interpretation was performed.

## Molecular-Pathologic Concordance

- Submitted diagnosis: {markdown_escape(concordance_row['submitted_diagnosis'])}
- Top molecular lineage group: `{markdown_escape(concordance_row['top_lineage_or_project_group'])}`
- Heuristic concordance class: `{markdown_escape(concordance_row['concordance_class'])}`
- Interpretation: {markdown_escape(concordance_row['rationale'])}

{markdown_escape(concordance_row['caveats'])}

## Figures

- `cup_top_project_similarity_barplot.pdf`
- `cup_major_group_similarity_barplot.pdf`
- `cup_ambiguity_summary.pdf`
- `cup_evidence_heatmap.pdf`

## Limitations

- The workflow reinterprets an unsupervised nearest-neighbor projection; it does not train or apply a supervised CUP classifier.
- TCGA is enriched for primary tumors and may not represent metastatic, treated, poorly differentiated, small-biopsy, or rare clinical tumors.
- Similarity depends on top-k choice, feature balance, exact assay harmonization, and the Euclidean or cosine projection metric used upstream.
- Mutation counts are proxies without callable-territory normalization, and driver indicators omit variant-level pathogenicity and clonality.
- Expression pathway scores are broad first-pass signatures. Missing expression features are median-imputed for projection and excluded from heuristic evidence.
- Arm-level copy-number summaries are not allele-specific integer copy-number calls and are not tissue-specific.
- Differential and molecular-pathologic comparisons use keyword heuristics and can miss synonyms or nuanced diagnostic categories.
- No independent known-primary validation, uncertainty calibration, prospective validation, or clinical outcome validation has been performed.

## Research-Use Caveat

This research/prototype analysis compares the query tumor to TCGA molecular reference profiles. It is not a clinically validated tissue-of-origin assay and should be interpreted only in conjunction with morphology, immunophenotype, imaging, and clinical findings.

The output must not be used as a final clinical diagnosis. The top lineage and project summaries are hypotheses for multidisciplinary review and future validation only.
"""


def interpret_cup_case(paths: CUPPaths, case_id: str, config: dict[str, Any] = CONFIG) -> dict[str, Any]:
    validate_case_id(case_id)
    top_k = int(config["top_k"])
    inputs = load_inputs(paths, case_id, top_k)
    lineages = load_lineage_groups(paths.lineage_groups, inputs["reference_metadata"])
    evidence_config = load_evidence_rules(
        paths.evidence_rules,
        set(inputs["feature_order"]),
        set(lineages["by_name"]),
    )
    project_summary, major_summary, summary_source = resolve_similarity_summaries(inputs)
    lineage_scores = initial_lineage_scores(lineages, project_summary, inputs["neighbors"])
    coverage = feature_coverage(inputs["global_qc"], len(inputs["feature_order"]))
    wes_coverage = modality_feature_coverage(inputs, "WES")
    wts_coverage = modality_feature_coverage(inputs, "WTS")
    ambiguity_metrics = build_ambiguity_metrics(
        case_id,
        inputs["sample_id"],
        lineage_scores,
        project_summary,
        major_summary,
        coverage,
        wes_coverage,
        wts_coverage,
        len(inputs["neighbors"]),
        config,
    )
    evaluated_rules = evaluate_rules(evidence_config, inputs)
    ranked_lineages = finish_ranked_lineages(lineage_scores, ambiguity_metrics, evaluated_rules, config)
    evidence = build_molecular_evidence(
        inputs,
        project_summary,
        major_summary,
        ranked_lineages,
        ambiguity_metrics,
        evaluated_rules,
        evidence_config,
    )
    differential = build_differential_summary(
        inputs["clinical_metadata"],
        ranked_lineages,
        lineages,
        inputs["reference_metadata"],
        config,
    )
    concordance = molecular_pathologic_concordance(
        case_id,
        inputs["sample_id"],
        inputs["clinical_metadata"].get("submitted_diagnosis"),
        str(ranked_lineages.iloc[0]["lineage_or_project_group"]),
        str(major_summary.iloc[0]["major_cancer_group"]),
        lineages,
        inputs["reference_metadata"],
    )
    report = build_report(
        case_id,
        inputs,
        project_summary,
        major_summary,
        ranked_lineages,
        ambiguity_metrics,
        evidence,
        differential,
        summary_source,
        concordance,
    )

    ensure_dir(paths.output_dir)
    write_tsv(ranked_lineages, paths.ranked_lineages)
    write_tsv(ambiguity_metrics, paths.ambiguity_metrics)
    write_tsv(evidence, paths.molecular_evidence)
    write_tsv(differential, paths.differential_summary)
    write_tsv(concordance, paths.concordance_summary)
    plot_project_similarity(project_summary, paths.top_project_figure)
    plot_major_group_similarity(major_summary, paths.major_group_figure)
    plot_ambiguity(ambiguity_metrics, paths.ambiguity_figure)
    plot_evidence_heatmap(
        ranked_lineages,
        major_summary,
        lineages,
        evaluated_rules,
        evidence_config,
        paths.evidence_heatmap,
    )
    write_text(report, paths.report)
    return {
        "inputs": inputs,
        "lineages": lineages,
        "evidence_config": evidence_config,
        "evaluated_rules": evaluated_rules,
        "project_summary": project_summary,
        "major_summary": major_summary,
        "ranked_lineages": ranked_lineages,
        "ambiguity_metrics": ambiguity_metrics,
        "molecular_evidence": evidence,
        "differential_summary": differential,
        "summary_source": summary_source,
        "concordance": concordance,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a research/prototype CUP interpretation from an existing locked TCGA clinical projection."
        )
    )
    parser.add_argument("--case-id", required=True, help="Clinical query case directory name and case_id value.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument(
        "--projection-dir",
        type=Path,
        default=None,
        help="Override results/clinical_projection/<case_id>.",
    )
    parser.add_argument(
        "--query-dir", type=Path, default=None, help="Override data/clinical_queries/<case_id>."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Override results/clinical_projection/<case_id>/cup_interpretation.",
    )
    parser.add_argument(
        "--lineage-config", type=Path, default=None, help="Override config/cup_lineage_groups.yaml."
    )
    parser.add_argument(
        "--evidence-config", type=Path, default=None, help="Override config/cup_feature_evidence_rules.yaml."
    )
    parser.add_argument("--top-k", type=int, default=CONFIG["top_k"], help="Number of ranked neighbors to interpret.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    case_id = validate_case_id(args.case_id)
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    paths = default_paths(
        root,
        case_id,
        projection_dir=args.projection_dir,
        query_dir=args.query_dir,
        output_dir=args.output_dir,
        lineage_config=args.lineage_config,
        evidence_config=args.evidence_config,
    )
    config = {**CONFIG, "top_k": args.top_k}
    outputs = interpret_cup_case(paths, case_id, config)
    top_lineage = outputs["ranked_lineages"].iloc[0]
    ambiguity = outputs["ambiguity_metrics"].iloc[0]
    logging.info(
        "CUP interpretation complete: case=%s sample=%s top_lineage=%s score=%.4f "
        "confidence=%s ambiguity=%s",
        case_id,
        outputs["inputs"]["sample_id"],
        top_lineage["lineage_or_project_group"],
        float(top_lineage["similarity_score"]),
        top_lineage["confidence_tier"],
        ambiguity["ambiguity_class"],
    )


if __name__ == "__main__":
    main()

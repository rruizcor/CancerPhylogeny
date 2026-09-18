#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "23_validate_known_primary_projection.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("known_primary_validation", SCRIPT)
validation = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["known_primary_validation"] = validation
spec.loader.exec_module(validation)


def fail(message: str) -> None:
    raise AssertionError(message)


# Verify exact self exclusion and the optional patient-level exclusion path.
matrix = np.array(
    [
        [0.0, 0.0],
        [0.0, 0.1],
        [0.0, 0.2],
        [2.0, 2.0],
        [3.0, 3.0],
        [4.0, 4.0],
    ]
)
patients = np.array(["PAT-A", "PAT-A", "PAT-B", "PAT-C", "PAT-D", "PAT-E"])
indices, distances = validation.compute_neighbor_indices(
    matrix,
    np.array([0]),
    top_k=2,
    metric="euclidean",
    patient_ids=patients,
)
if 0 in indices[0] or indices[0].tolist() != [1, 2]:
    fail("Nearest-neighbor validation did not exclude the query sample itself")
if not np.allclose(distances[0], [0.1, 0.2]):
    fail("Synthetic Euclidean neighbor distances are incorrect")

patient_indices, _ = validation.compute_neighbor_indices(
    matrix,
    np.array([0]),
    top_k=2,
    metric="euclidean",
    patient_ids=patients,
    exclude_same_patient=True,
)
if 0 in patient_indices[0] or 1 in patient_indices[0] or patient_indices[0, 0] != 2:
    fail("Same-patient exclusion did not remove every sample from the query patient")


# Verify that the lineage configuration is complete, unique, and parsed as expected.
tmp = Path(tempfile.mkdtemp(prefix="known_primary_validation_synthetic_"))
lineage_path = tmp / "cup_lineage_groups.yaml"
lineage_path.write_text(
    yaml.safe_dump(
        {
            "lineages": [
                {"name": "lineage_1", "projects": ["P1", "P2"]},
                {"name": "lineage_2", "projects": ["P3", "P4"]},
                {"name": "lineage_3", "projects": ["P5", "P6"]},
                {"name": "lineage_4", "projects": ["P7", "P8", "P9"]},
            ]
        },
        sort_keys=False,
    )
)
project_to_lineage, lineage_order = validation.load_lineage_mapping(
    lineage_path, {f"P{index}" for index in range(1, 10)}
)
if project_to_lineage["P3"] != "lineage_2" or lineage_order != [
    "lineage_1",
    "lineage_2",
    "lineage_3",
    "lineage_4",
]:
    fail("CUP lineage mapping was not parsed deterministically")


# Construct explicit neighbor rankings so project top-1/top-3/top-5 behavior is testable.
projects = ["P1", "P1", "P2", "P3", "P3", "P4", "P5", "P5", "P6", "P7", "P8", "P9"]
lineages = [project_to_lineage[project] for project in projects]
major_by_lineage = {
    "lineage_1": "major_1",
    "lineage_2": "major_2",
    "lineage_3": "major_3",
    "lineage_4": "major_4",
}
metadata = pd.DataFrame(
    {
        "sample_barcode": [f"S{index:02d}" for index in range(len(projects))],
        "patient_barcode": [f"PAT-{index:02d}" for index in range(len(projects))],
        "project_id": [f"TCGA-{project}" for project in projects],
        "project_code": projects,
        "disease_type": [f"Synthetic {project}" for project in projects],
        "primary_site": [f"Site {project}" for project in projects],
        "major_cancer_group": [major_by_lineage[lineage] for lineage in lineages],
        "cup_lineage_group": lineages,
    }
)
query_indices = np.array([0, 3, 6])
neighbor_indices = np.array(
    [
        [1, 2, 4, 5, 7],       # True P1 ranks first.
        [2, 5, 4, 7, 8],       # True P3 ranks third.
        [2, 5, 8, 9, 7],       # True P5 ranks fifth.
    ]
)
neighbor_distances = np.tile(np.array([0.1, 0.2, 0.3, 0.4, 0.5]), (3, 1))
predictions, internal = validation.evaluate_queries(
    metadata,
    query_indices,
    neighbor_indices,
    neighbor_distances,
    np.ones(len(metadata)),
    top_k=5,
    metric="euclidean",
    exclude_same_patient=False,
)

if predictions["correct_top_project"].tolist() != [True, False, False]:
    fail("Synthetic top-1 project accuracy flags are incorrect")
if predictions["true_project_in_top3"].tolist() != [True, True, False]:
    fail("Synthetic top-3 project accuracy flags are incorrect")
if predictions["true_project_in_top5"].tolist() != [True, True, True]:
    fail("Synthetic top-5 project accuracy flags are incorrect")
if not predictions["correct_top_lineage"].all() or not predictions["true_lineage_in_top3"].all():
    fail("Synthetic lineage recovery did not use the configured project grouping")
if not predictions["correct_top_major_group"].all():
    fail("Synthetic broad major-group accuracy flags are incorrect")

allowed_ambiguity = {
    "low_ambiguity",
    "moderate_ambiguity",
    "high_ambiguity",
    "uninterpretable_due_to_missing_features",
}
if predictions["ambiguity_class"].isna().any() or not set(predictions["ambiguity_class"]).issubset(
    allowed_ambiguity
):
    fail("Ambiguity classes are missing or invalid")
if not predictions["neighbor_entropy"].between(0, 1).all():
    fail("Normalized neighbor entropy is outside [0, 1]")

project_summary = validation.project_performance(predictions, internal)
p3 = project_summary.set_index("project_code").loc["P3"]
if p3["common_confusions"] != "P2:1":
    fail("Common project confusion was not summarized correctly")

project_confusion = validation.confusion_matrix_table(
    predictions,
    "true_project_code",
    "predicted_top_project",
    "true_project_code",
)
if project_confusion.empty or project_confusion.iloc[:, 1:].to_numpy().sum() != len(predictions):
    fail("Project confusion matrix does not contain all synthetic queries")

prioritized = validation.prioritized_confusion_analysis(
    pd.DataFrame(
        {
            "true_project_code": ["COAD", "READ", "COAD"],
            "predicted_top_project": ["READ", "COAD", "COAD"],
        }
    )
).set_index("comparison")
coad_read = prioritized.loc["COAD vs READ"]
if (
    coad_read["group_a_to_group_b"] != 1
    or coad_read["group_b_to_group_a"] != 1
    or not np.isclose(coad_read["cross_confusion_fraction"], 2 / 3)
):
    fail("Prespecified directional confusion analysis is incorrect")

ambiguity_summary = validation.ambiguity_performance(predictions)
required_dimensions = {
    "ambiguity_class",
    "feature_coverage_bin",
    "top1_top2_project_margin_bin",
    "neighbor_entropy_bin",
    "ambiguity_comparison",
}
if not required_dimensions.issubset(set(ambiguity_summary["analysis_dimension"])):
    fail("Ambiguity performance lacks one or more required calibration-style summaries")

print("Known-primary validation synthetic checks passed")

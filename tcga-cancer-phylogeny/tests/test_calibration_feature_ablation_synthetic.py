#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "24_calibrate_projection_confidence_and_feature_ablation.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("calibration_ablation", SCRIPT)
analysis = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["calibration_ablation"] = analysis
spec.loader.exec_module(analysis)


def fail(message: str) -> None:
    raise AssertionError(message)


# Feature grouping must be driven by locked feature_class and modality metadata.
records = [
    {"name": "mutation_count", "feature_class": "mutation_count_proxy", "clinical_modality": "WES", "source_layer": "mutation"},
    {"name": "driver_binary", "feature_class": "driver_mutation_binary", "clinical_modality": "WES", "source_layer": "mutation"},
    {"name": "driver_count", "feature_class": "driver_mutation_count", "clinical_modality": "WES", "source_layer": "mutation"},
    {"name": "purity", "feature_class": "purity_ploidy", "clinical_modality": "WES", "source_layer": "purity"},
    {"name": "aneuploidy", "feature_class": "copy_number_aneuploidy", "clinical_modality": "WES", "source_layer": "copy_number"},
    {"name": "pathway", "feature_class": "expression_pathway_score", "clinical_modality": "WTS", "source_layer": "expression"},
    {"name": "expression_pc", "feature_class": "expression_pca", "clinical_modality": "WTS", "source_layer": "expression"},
]
feature_order = [record["name"] for record in records]
groups, definitions, unclassified = analysis.build_feature_groups(feature_order, records)
if unclassified:
    fail("Supported synthetic features were incorrectly marked unclassified")
if set(groups) != set(analysis.FEATURE_SET_ORDER):
    fail("Feature grouping does not contain all required ablation sets")
expected_counts = {
    "all_features": 7,
    "mutation_counts_only": 1,
    "driver_indicators_only": 1,
    "driver_counts_only": 1,
    "mutation_features_combined": 3,
    "purity_ploidy_only": 1,
    "copy_number_aneuploidy_only": 1,
    "pathway_scores_only": 1,
    "expression_pcs_only": 1,
    "expression_combined": 2,
    "WES_like_features": 5,
    "WTS_like_features": 2,
    "WES_plus_WTS_combined": 7,
}
if {name: len(values) for name, values in groups.items()} != expected_counts:
    fail("Feature group membership counts are incorrect")
if set(definitions["feature_set"]) != set(analysis.FEATURE_SET_ORDER):
    fail("Feature-group definition output lacks a required feature set")


# Exact nearest neighbors must exclude the query sample itself.
matrix = np.array(
    [
        [0.0, 0.0],
        [0.0, 0.1],
        [0.1, 0.0],
        [5.0, 5.0],
        [5.0, 5.1],
        [5.1, 5.0],
        [10.0, 10.0],
        [10.0, 10.1],
        [10.1, 10.0],
    ]
)
query_indices = np.arange(len(matrix))
neighbor_indices, _ = analysis.compute_neighbor_indices(
    matrix, query_indices, top_k=2, metric="euclidean"
)
if any(query in neighbors for query, neighbors in zip(query_indices, neighbor_indices, strict=True)):
    fail("Feature-ablation nearest-neighbor search retained a self neighbor")

metadata = pd.DataFrame(
    {
        "sample_barcode": [f"S{index}" for index in range(9)],
        "patient_barcode": [f"P{index}" for index in range(9)],
        "project_id": ["TCGA-A"] * 3 + ["TCGA-B"] * 3 + ["TCGA-C"] * 3,
        "project_code": ["A"] * 3 + ["B"] * 3 + ["C"] * 3,
        "disease_type": ["Synthetic"] * 9,
        "primary_site": ["Synthetic"] * 9,
        "major_cancer_group": ["M1"] * 3 + ["M2"] * 6,
        "cup_lineage_group": ["L1"] * 3 + ["L2"] * 6,
    }
)
evaluation = analysis.evaluate_feature_matrix(
    matrix,
    metadata,
    query_indices,
    np.ones(9),
    ["L1", "L2"],
    {"L1": {"M1"}, "L2": {"M2"}},
    top_k=2,
    metric="euclidean",
)
if not evaluation["correct_top_project"].all():
    fail("Synthetic clustered matrix did not recover its project labels")
if not evaluation["correct_top_lineage"].all() or not evaluation["correct_top_major_group"].all():
    fail("Synthetic clustered matrix did not recover lineage/major-group labels")
if not evaluation["project_neighbor_concentration"].eq(1.0).all():
    fail("Synthetic top-label neighbor concentration is incorrect")

# A true class absent from the neighbors ranks after represented classes, not
# according to the lexical order of zero-score categories.
absent_rank = analysis.score_neighbor_labels(
    np.array(["A", "B", "C", "D"]),
    np.array([3]),
    np.array([[0, 1]]),
    np.array([[0.1, 0.2]]),
)
if absent_rank["true_rank"][0] != 3:
    fail("An absent true class did not rank immediately after represented neighbor classes")


# Threshold calibration should orient favorable metrics correctly and assign populated tiers.
n = 180
x = np.linspace(0.0, 1.0, n)
correct = x >= 0.38
calibration_evaluation = pd.DataFrame(
    {
        "sample_barcode": [f"Q{index:03d}" for index in range(n)],
        "correct_top_project": correct,
        "correct_top_lineage": correct,
        "correct_top_major_group": correct,
    }
)
for level in analysis.LEVEL_SPECS:
    calibration_evaluation[f"{level}_margin"] = x
    calibration_evaluation[f"{level}_entropy"] = 1.0 - x
    calibration_evaluation[f"{level}_top1_score"] = x
    calibration_evaluation[f"{level}_neighbor_concentration"] = x

source = pd.DataFrame(
    {
        "sample_barcode": calibration_evaluation["sample_barcode"],
        "true_project_code": ["A"] * n,
        "true_project_id": ["TCGA-A"] * n,
        "true_major_cancer_group": ["M1"] * n,
        "true_cup_lineage_group": ["L1"] * n,
        "correct_top_project": correct,
        "correct_top_lineage": correct,
        "correct_top_major_group": correct,
        "ambiguity_class": np.where(x >= 0.75, "low_ambiguity", np.where(x < 0.35, "high_ambiguity", "moderate_ambiguity")),
    }
)
thresholds, calibrated, performance, holdout = analysis.calibrate_confidence(
    calibration_evaluation, source, seed=123
)
if len(thresholds) != 12:
    fail("Calibration did not produce four metric thresholds for all three prediction levels")
for level in analysis.LEVEL_SPECS:
    margin = thresholds[(thresholds["prediction_level"] == level) & (thresholds["metric"] == "top1_top2_margin")].iloc[0]
    entropy = thresholds[(thresholds["prediction_level"] == level) & (thresholds["metric"] == "entropy")].iloc[0]
    if margin["threshold_high_confidence"] <= margin["threshold_low_confidence"]:
        fail("Higher-is-better margin thresholds are reversed")
    if entropy["threshold_high_confidence"] >= entropy["threshold_low_confidence"]:
        fail("Lower-is-better entropy thresholds are reversed")
    tiers = set(calibrated[f"calibrated_{level}_confidence"])
    if not {"low_confidence", "high_confidence"}.issubset(tiers):
        fail("Calibrated confidence tiers were not assigned across favorable/unfavorable cases")
    tier_rows = performance[
        (performance["prediction_level"] == level)
        & (performance["summary_type"] == "calibrated_confidence_tier")
    ].set_index("stratum")
    if tier_rows.loc["high_confidence", "observed_accuracy"] <= tier_rows.loc["low_confidence", "observed_accuracy"]:
        fail("High-confidence synthetic calls are not more accurate than low-confidence calls")
if holdout.sum() == 0 or "sample_barcode" not in calibrated.columns:
    fail("Calibration output lost its held-out partition or sample identifier")

print("Calibration and feature-ablation synthetic checks passed")

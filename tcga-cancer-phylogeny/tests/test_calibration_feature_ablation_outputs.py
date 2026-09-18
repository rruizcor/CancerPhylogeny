#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "validation" / "calibration_feature_ablation"


def fail(message: str) -> None:
    raise AssertionError(message)


required_tables = [
    "calibrated_confidence_thresholds.tsv",
    "known_primary_predictions_calibrated.tsv",
    "calibration_performance_summary.tsv",
    "feature_group_definitions.tsv",
    "feature_ablation_performance_summary.tsv",
    "feature_ablation_project_performance.tsv",
    "feature_ablation_lineage_performance.tsv",
    "feature_ablation_project_confusion_matrices.tsv",
    "feature_ablation_lineage_confusion_matrices.tsv",
]
required_figures = [
    "calibration_accuracy_by_confidence_tier.pdf",
    "calibration_margin_vs_accuracy.pdf",
    "calibration_entropy_vs_accuracy.pdf",
    "feature_ablation_overall_accuracy.pdf",
    "feature_ablation_lineage_accuracy.pdf",
    "feature_ablation_project_accuracy_heatmap.pdf",
    "feature_ablation_confusion_heatmap_best_model.pdf",
]
report_path = OUTPUT / "calibration_feature_ablation_report.md"
for name in required_tables + required_figures:
    path = OUTPUT / name
    if not path.exists() or path.stat().st_size == 0:
        fail(f"Required calibration/ablation output missing or empty: {path}")
for name in required_figures:
    if (OUTPUT / name).read_bytes()[:5] != b"%PDF-":
        fail(f"Calibration/ablation figure is not a valid PDF: {name}")
if not report_path.exists() or report_path.stat().st_size == 0:
    fail("Calibration/ablation report is missing or empty")

thresholds = pd.read_csv(OUTPUT / "calibrated_confidence_thresholds.tsv", sep="\t")
required_threshold_columns = {
    "prediction_level",
    "metric",
    "threshold_low_confidence",
    "threshold_moderate_confidence",
    "threshold_high_confidence",
    "observed_accuracy_low_confidence",
    "observed_accuracy_moderate_confidence",
    "observed_accuracy_high_confidence",
    "recommended_use",
    "notes",
}
if set(thresholds.columns) != required_threshold_columns or len(thresholds) != 12:
    fail("Calibrated threshold table does not match the required 3-level by 4-metric contract")
if set(thresholds["prediction_level"]) != {"project", "lineage", "major_group"}:
    fail("Calibrated thresholds lack a prediction level")
if set(thresholds["metric"]) != {
    "top1_top2_margin",
    "entropy",
    "top1_score",
    "neighbor_concentration",
}:
    fail("Calibrated thresholds lack a required metric")

calibrated = pd.read_csv(OUTPUT / "known_primary_predictions_calibrated.tsv", sep="\t")
required_calibrated = {
    "sample_barcode",
    "calibrated_project_confidence",
    "calibrated_lineage_confidence",
    "calibrated_major_group_confidence",
    "calibrated_project_reliable_binary",
    "calibrated_lineage_reliable_binary",
    "calibrated_major_group_reliable_binary",
    "calibration_notes",
}
if not required_calibrated.issubset(calibrated.columns):
    fail("Calibrated prediction table lacks required identifier/tier fields")
if len(calibrated) != 10201 or calibrated["sample_barcode"].duplicated().any():
    fail("Calibrated prediction table does not contain all 10,201 unique validation samples")
allowed_tiers = {"low_confidence", "moderate_confidence", "high_confidence"}
for level in ["project", "lineage", "major_group"]:
    if not set(calibrated[f"calibrated_{level}_confidence"]).issubset(allowed_tiers):
        fail(f"Invalid calibrated {level} confidence tier")

summary = pd.read_csv(OUTPUT / "calibration_performance_summary.tsv", sep="\t")
if not {"calibrated_confidence_tier", "reliable_binary_performance", "previous_heuristic_ambiguity"}.issubset(
    set(summary["summary_type"])
):
    fail("Calibration summary lacks tier, reliable-binary, or heuristic comparison rows")

feature_groups = pd.read_csv(OUTPUT / "feature_group_definitions.tsv", sep="\t")
required_feature_sets = {
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
}
if not required_feature_sets.issubset(set(feature_groups["feature_set"])):
    fail("Feature-group definitions lack one or more required groups")

ablation = pd.read_csv(OUTPUT / "feature_ablation_performance_summary.tsv", sep="\t")
if set(ablation["feature_set"]) != required_feature_sets or len(ablation) != 13:
    fail("Feature-ablation summary does not contain exactly the 13 required feature sets")
if not ablation["n_samples_validated"].eq(10201).all():
    fail("Full ablation did not validate all 10,201 samples for every feature set")
expected_counts = {
    "all_features": 97,
    "mutation_counts_only": 2,
    "driver_indicators_only": 29,
    "driver_counts_only": 29,
    "mutation_features_combined": 60,
    "purity_ploidy_only": 2,
    "copy_number_aneuploidy_only": 4,
    "pathway_scores_only": 11,
    "expression_pcs_only": 20,
    "expression_combined": 31,
    "WES_like_features": 66,
    "WTS_like_features": 31,
    "WES_plus_WTS_combined": 97,
}
observed_counts = ablation.set_index("feature_set")["n_features"].to_dict()
if observed_counts != expected_counts:
    fail("Feature-ablation feature counts do not match the locked 97-feature taxonomy")

for name in [
    "feature_ablation_project_confusion_matrices.tsv",
    "feature_ablation_lineage_confusion_matrices.tsv",
]:
    confusion = pd.read_csv(OUTPUT / name, sep="\t")
    if confusion.empty or not {"feature_set", "true_label", "predicted_label", "count"}.issubset(
        confusion.columns
    ):
        fail(f"Ablation confusion table is missing or malformed: {name}")

report = report_path.read_text()
required_sections = [
    "## Purpose",
    "## Calibration Approach",
    "## Recommended Confidence Thresholds",
    "## Accuracy by Calibrated Confidence Tier",
    "## Feature-Ablation Design",
    "## WES-Like Versus WTS-Like Performance",
    "## Which Feature Sets Perform Best",
    "## Which Projects/Lineages Are Robust",
    "## Which Projects/Lineages Are Unstable or Confused",
    "## Implications for CUP Reporting",
    "## Limitations",
    "## Recommended Changes to CUP Workflow",
]
for section in required_sections:
    if section not in report:
        fail(f"Calibration/ablation report lacks required section: {section}")
if "This remains internal TCGA validation, not external clinical validation." not in report:
    fail("Calibration/ablation report lacks the internal-validation caveat")
if "Similarity weights are not probabilities" not in report:
    fail("Calibration/ablation report lacks the required score-interpretation boundary")

text_outputs = "\n".join((OUTPUT / name).read_text(errors="ignore") for name in required_tables) + report
unsupported_claims = [
    "similarity scores are probabilities",
    "similarity weights are probabilities",
    "clinically calibrated confidence",
    "external clinical validation demonstrates",
]
if any(claim in text_outputs.lower() for claim in unsupported_claims):
    fail("Calibration/ablation outputs make an unsupported probability or clinical-validation claim")

print("Calibration and feature-ablation output checks passed")

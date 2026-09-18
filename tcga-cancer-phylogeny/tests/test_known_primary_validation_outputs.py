#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "validation" / "known_primary_projection"


def fail(message: str) -> None:
    raise AssertionError(message)


required_tables = [
    "known_primary_validation_cohort.tsv",
    "known_primary_projection_predictions.tsv",
    "known_primary_project_performance.tsv",
    "known_primary_lineage_performance.tsv",
    "known_primary_major_group_performance.tsv",
    "known_primary_overall_metrics.tsv",
    "known_primary_project_confusion_matrix.tsv",
    "known_primary_lineage_confusion_matrix.tsv",
    "known_primary_major_group_confusion_matrix.tsv",
    "known_primary_ambiguity_performance.tsv",
]
required_figures = [
    "known_primary_project_confusion_matrix.pdf",
    "known_primary_lineage_confusion_matrix.pdf",
    "known_primary_topk_accuracy_by_project.pdf",
    "known_primary_accuracy_by_lineage.pdf",
    "known_primary_accuracy_by_ambiguity_class.pdf",
    "known_primary_margin_vs_correctness.pdf",
    "known_primary_entropy_vs_correctness.pdf",
]
report_path = OUTPUT / "known_primary_validation_report.md"

for name in required_tables + required_figures:
    path = OUTPUT / name
    if not path.exists() or path.stat().st_size == 0:
        fail(f"Required known-primary validation output missing or empty: {path}")
for name in required_figures:
    if (OUTPUT / name).read_bytes()[:5] != b"%PDF-":
        fail(f"Known-primary validation figure is not a valid PDF: {name}")
if not report_path.exists() or report_path.stat().st_size == 0:
    fail("Known-primary validation report is missing or empty")

cohort = pd.read_csv(OUTPUT / "known_primary_validation_cohort.tsv", sep="\t")
required_cohort_columns = {
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
}
if set(cohort.columns) != required_cohort_columns:
    fail("Validation cohort columns do not match the required contract")
if len(cohort) != 10201 or cohort["sample_barcode"].duplicated().any():
    fail("Validation cohort does not preserve all 10,201 unique atlas samples")
if not cohort["included_in_validation"].all() or set(cohort["validation_split"]) != {
    "full_leave_one_out_against_full_atlas"
}:
    fail("Default output is not the requested full leave-one-out validation cohort")

predictions = pd.read_csv(OUTPUT / "known_primary_projection_predictions.tsv", sep="\t")
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
}
if set(predictions.columns) != required_prediction_columns:
    fail("Per-sample prediction columns do not match the required contract")
if len(predictions) != 10201 or predictions["sample_barcode"].duplicated().any():
    fail("Full validation prediction table does not contain 10,201 unique queries")
if predictions["true_project_code"].nunique() != 33:
    fail("Full validation predictions do not represent all 33 TCGA projects")
if not predictions["feature_coverage"].between(0, 1).all():
    fail("Feature coverage values are outside [0, 1]")

metrics = pd.read_csv(OUTPUT / "known_primary_overall_metrics.tsv", sep="\t")
metric_names = set(metrics["metric"].astype(str))
required_metrics = {
    "n_samples_validated",
    "n_projects",
    "n_lineage_groups",
    "top1_project_accuracy",
    "top3_project_accuracy",
    "top5_project_accuracy",
    "top1_lineage_accuracy",
    "top3_lineage_accuracy",
    "top1_major_group_accuracy",
    "median_project_margin",
    "median_lineage_margin",
    "high_ambiguity_fraction",
    "moderate_ambiguity_fraction",
    "low_ambiguity_fraction",
    "limitations",
}
if not required_metrics.issubset(metric_names):
    fail("Overall validation metrics lack one or more required fields")

for name, row_label in [
    ("known_primary_project_confusion_matrix.tsv", "true_project_code"),
    ("known_primary_lineage_confusion_matrix.tsv", "true_cup_lineage_group"),
    ("known_primary_major_group_confusion_matrix.tsv", "true_major_cancer_group"),
]:
    matrix = pd.read_csv(OUTPUT / name, sep="\t")
    if matrix.empty or matrix.columns[0] != row_label:
        fail(f"Confusion matrix is missing or malformed: {name}")
    if matrix.iloc[:, 1:].to_numpy().sum() != len(predictions):
        fail(f"Confusion matrix does not account for all validation queries: {name}")

ambiguity = pd.read_csv(OUTPUT / "known_primary_ambiguity_performance.tsv", sep="\t")
required_dimensions = {
    "ambiguity_class",
    "feature_coverage_bin",
    "top1_top2_project_margin_bin",
    "neighbor_entropy_bin",
    "ambiguity_comparison",
}
if not required_dimensions.issubset(set(ambiguity["analysis_dimension"])):
    fail("Ambiguity output lacks required accuracy/calibration-style dimensions")

report = report_path.read_text()
required_sections = [
    "## Purpose",
    "## Validation Design",
    "## Cohort Composition",
    "## Overall Performance",
    "## Project-Level Performance",
    "## Lineage-Level Performance",
    "## Major-Group Performance",
    "## Common Confusions",
    "## Ambiguity Analysis",
    "## Failure Modes",
    "## Limitations",
    "## Recommended Improvements Before Clinical Use",
]
for section in required_sections:
    if section not in report:
        fail(f"Known-primary validation report lacks required section: {section}")
if "This is internal validation using TCGA-derived samples. It is not external clinical validation" not in report:
    fail("Validation report lacks the explicit internal-validation caveat")
unsupported_claims = [
    "This is external clinical validation",
    "externally validated diagnostic accuracy",
    "clinically validated tissue-of-origin assay",
]
if any(claim.lower() in report.lower() for claim in unsupported_claims):
    fail("Validation report makes an unsupported external clinical-validation claim")

print("Known-primary validation output checks passed")

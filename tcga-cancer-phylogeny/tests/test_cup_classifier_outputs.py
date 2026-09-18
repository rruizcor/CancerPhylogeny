#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CASE_ID = "example_cup_case"
OUTPUT = ROOT / "results" / "clinical_projection" / CASE_ID / "cup_interpretation"


def fail(message: str) -> None:
    raise AssertionError(message)


required_files = [
    OUTPUT / "cup_ranked_lineage_interpretation.tsv",
    OUTPUT / "cup_ambiguity_metrics.tsv",
    OUTPUT / "cup_molecular_evidence_table.tsv",
    OUTPUT / "cup_differential_diagnosis_constrained_summary.tsv",
    OUTPUT / "cup_molecular_pathologic_concordance.tsv",
    OUTPUT / "cup_molecular_interpretation_report.md",
    OUTPUT / "cup_top_project_similarity_barplot.pdf",
    OUTPUT / "cup_major_group_similarity_barplot.pdf",
    OUTPUT / "cup_ambiguity_summary.pdf",
    OUTPUT / "cup_evidence_heatmap.pdf",
]
for path in required_files:
    if not path.exists() or path.stat().st_size == 0:
        fail(f"Required CUP interpretation output missing or empty: {path}")

ranked = pd.read_csv(OUTPUT / "cup_ranked_lineage_interpretation.tsv", sep="\t")
required_ranked_columns = {
    "rank",
    "lineage_or_project_group",
    "supporting_tcga_projects",
    "similarity_score",
    "confidence_tier",
    "evidence_summary",
    "caveats",
}
if set(ranked.columns) != required_ranked_columns:
    fail("Ranked CUP lineage table columns do not match the required contract")
if ranked.iloc[0]["lineage_or_project_group"] != "lower_GI":
    fail("Example CUP top lineage is not lower_GI")
if not np.isclose(ranked["similarity_score"].sum(), 1.0, atol=1e-8):
    fail("Ranked CUP lineage scores do not sum to one")
allowed_tiers = {"high_support", "moderate_support", "weak_support", "indeterminate"}
if ranked["confidence_tier"].isna().any() or not set(ranked["confidence_tier"]).issubset(allowed_tiers):
    fail("Ranked CUP lineage confidence tiers are missing or invalid")

ambiguity = pd.read_csv(OUTPUT / "cup_ambiguity_metrics.tsv", sep="\t")
required_ambiguity_columns = {
    "top1_score",
    "top2_score",
    "top1_top2_margin",
    "entropy_like_score",
    "number_of_projects_with_meaningful_similarity",
    "feature_coverage",
    "wes_feature_coverage",
    "wts_feature_coverage",
    "ambiguity_class",
}
if not required_ambiguity_columns.issubset(ambiguity.columns):
    fail("CUP ambiguity metrics lack required fields")
allowed_ambiguity = {
    "low_ambiguity",
    "moderate_ambiguity",
    "high_ambiguity",
    "uninterpretable_due_to_missing_features",
}
if ambiguity.iloc[0]["ambiguity_class"] not in allowed_ambiguity:
    fail("CUP ambiguity class is invalid")
if ambiguity.iloc[0]["ambiguity_class"] != "moderate_ambiguity":
    fail("Example CUP ambiguity class changed unexpectedly")
if not np.isclose(float(ambiguity.iloc[0]["feature_coverage"]), 0.79381, atol=1e-6):
    fail("Example CUP feature coverage does not match the projection QC")
if not np.isclose(float(ambiguity.iloc[0]["wes_feature_coverage"]), 1.0, atol=1e-8):
    fail("Example CUP WES coverage does not match the projection QC")
if not np.isclose(float(ambiguity.iloc[0]["wts_feature_coverage"]), 11 / 31, atol=1e-8):
    fail("Example CUP WTS coverage does not match the projection QC")

evidence = pd.read_csv(OUTPUT / "cup_molecular_evidence_table.tsv", sep="\t")
expected_evidence_types = {
    "nearest-neighbor TCGA project evidence",
    "major cancer group evidence",
    "expression/pathway evidence",
    "driver mutation evidence",
    "copy-number/aneuploidy evidence",
    "immune/stromal/EMT evidence",
    "contradictory evidence",
    "missing evidence",
}
if set(evidence["evidence_type"]) != expected_evidence_types or len(evidence) != 8:
    fail("CUP molecular evidence table does not contain the eight required evidence types")
if evidence[["finding", "strength", "notes"]].isna().any().any():
    fail("CUP molecular evidence table contains missing required interpretation fields")

differential = pd.read_csv(
    OUTPUT / "cup_differential_diagnosis_constrained_summary.tsv", sep="\t"
)
if differential.empty or differential.iloc[0]["candidate_text"] != "gastrointestinal":
    fail("Example differential-constrained summary did not identify gastrointestinal as the top candidate")
if "lower_GI" not in str(differential.iloc[0]["matched_lineages_or_projects"]):
    fail("Example differential-constrained summary did not map gastrointestinal to lower_GI")

concordance = pd.read_csv(OUTPUT / "cup_molecular_pathologic_concordance.tsv", sep="\t")
required_concordance_columns = {
    "case_id",
    "query_sample_id",
    "submitted_diagnosis",
    "top_lineage_or_project_group",
    "top_major_cancer_group",
    "concordance_class",
    "rationale",
    "caveats",
}
if set(concordance.columns) != required_concordance_columns or len(concordance) != 1:
    fail("Molecular-pathologic concordance table is missing or malformed")
if concordance.iloc[0]["concordance_class"] != "not_assessable":
    fail("Synthetic CUP diagnosis should be not_assessable for concordance")
if "heuristic" not in concordance.iloc[0]["caveats"].lower():
    fail("Molecular-pathologic concordance table lacks its heuristic caveat")

report = (OUTPUT / "cup_molecular_interpretation_report.md").read_text()
required_sections = [
    "## Case Metadata",
    "## Feature Completeness",
    "## Top TCGA-Like Projects",
    "## Top Lineage/Major Group Interpretation",
    "## Ambiguity/Confidence Assessment",
    "## Molecular Evidence Table",
    "## Differential Diagnosis Constrained Interpretation",
    "## Actionability Placeholder",
    "## Molecular-Pathologic Concordance",
    "## Limitations",
    "## Research-Use Caveat",
]
for section in required_sections:
    if section not in report:
        fail(f"CUP report lacks required section: {section}")
required_caveat = (
    "This research/prototype analysis compares the query tumor to TCGA molecular reference profiles. "
    "It is not a clinically validated tissue-of-origin assay"
)
if required_caveat not in report:
    fail("CUP report lacks the required clinical-use caveat")
if "Actionability integration is planned but not implemented in this prototype." not in report:
    fail("CUP report lacks the actionability placeholder")
if "This is a clinically validated tissue-of-origin assay" in report or "The diagnosis is" in report:
    fail("CUP report makes an unsupported validated-diagnosis claim")

for pdf in required_files[-4:]:
    if pdf.read_bytes()[:5] != b"%PDF-":
        fail(f"CUP figure is not a valid PDF file: {pdf}")

lineage_config = yaml.safe_load((ROOT / "config" / "cup_lineage_groups.yaml").read_text())
assigned_projects = [project for group in lineage_config["lineages"] for project in group["projects"]]
reference_projects = set(
    pd.read_csv(
        ROOT / "results" / "reference_atlas" / "tcga_sample_reference_metadata.tsv",
        sep="\t",
    )["project_code"]
)
if len(assigned_projects) != len(set(assigned_projects)):
    fail("CUP lineage configuration assigns a TCGA project more than once")
if set(assigned_projects) != reference_projects:
    fail("CUP lineage configuration does not cover all reference TCGA projects")
if "other_or_mixed" not in {group["name"] for group in lineage_config["lineages"]}:
    fail("CUP lineage configuration lacks the required other_or_mixed group")

print("CUP interpretation output checks passed")

#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


paths = {
    "main_report": ROOT / "results" / "reports" / "integrated_tcga_cancer_phylogeny_report.md",
    "executive_summary": ROOT / "results" / "reports" / "integrated_tcga_cancer_phylogeny_executive_summary.md",
    "methods": ROOT / "results" / "reports" / "integrated_tcga_cancer_phylogeny_methods.md",
    "limitations": ROOT / "results" / "reports" / "integrated_tcga_cancer_phylogeny_limitations_next_steps.md",
    "html_report": ROOT / "results" / "reports" / "integrated_tcga_cancer_phylogeny_report.html",
    "summary": ROOT / "results" / "tables" / "integrated_project_summary.tsv",
    "outputs_index": ROOT / "results" / "tables" / "integrated_project_outputs_index.tsv",
}

missing_or_empty = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
if missing_or_empty:
    fail(f"Missing or empty integrated report outputs: {', '.join(missing_or_empty)}")

report = paths["main_report"].read_text(encoding="utf-8")
required_phrases = [
    "TCGA Multilevel Molecular Similarity, Clonal-Structure, and Clinical Projection Framework",
    "## Level 1 Results",
    "## Level 2 Results",
    "## Level 3 Results",
    "## TCGA Sample-Level Reference Atlas",
    "## Clinical WES/WTS Projection Prototype",
    "## CUP Molecular Interpretation Prototype",
    "## Current PyClone-VI/PhyloWGS Status",
    "## Limitations",
    "97 retained features",
    "descriptive similarity weights, not tissue-of-origin probabilities",
    "not a full copy-number-aware clonal phylogeny",
    "does not infer definitive phylogenetic branching",
]
for phrase in required_phrases:
    if phrase not in report:
        fail(f"Main report is missing required phrase: {phrase}")

false_claims = [
    "is a full copy-number-aware clonal phylogeny",
    "full copy-number-aware clonal phylogeny was inferred",
    "definitive branching order was inferred",
    "definitive phylogenetic branching was inferred",
]
lower_report = report.lower()
for phrase in false_claims:
    if phrase in lower_report:
        fail(f"Main report appears to make a false clonal-phylogeny claim: {phrase}")

if "70 retained features" in report or "70-feature" in report:
    fail("Main report contains an outdated 70-feature atlas claim")

summary = pd.read_csv(paths["summary"], sep="\t")
required_summary_columns = {"level", "metric", "value", "source_file", "notes"}
if not required_summary_columns.issubset(summary.columns):
    fail("Integrated project summary lacks required columns")
required_metrics = {
    "tcga_projects_represented",
    "biologic_features_used",
    "groups_with_project_level_trees",
    "total_selected_pilot_samples",
    "total_clustered_mutations",
    "multicluster_subclonal_like",
    "number_of_retained_features",
    "percent_feature_coverage",
    "top_lineage",
    "ambiguity_class",
    "samples_ready_for_pyclone_vi",
    "samples_ready_for_phylowgs",
}
if not required_metrics.issubset(set(summary["metric"].astype(str))):
    fail("Integrated project summary lacks expected metrics")

outputs_index = pd.read_csv(paths["outputs_index"], sep="\t")
required_index_columns = {"level", "output_type", "file_path", "description", "exists", "notes"}
if not required_index_columns.issubset(outputs_index.columns):
    fail("Integrated output index lacks required columns")

expected_files = {
    "results/reports/integrated_tcga_cancer_phylogeny_report.md",
    "results/reports/integrated_tcga_cancer_phylogeny_executive_summary.md",
    "results/reports/integrated_tcga_cancer_phylogeny_methods.md",
    "results/reports/integrated_tcga_cancer_phylogeny_limitations_next_steps.md",
    "results/tables/integrated_project_summary.tsv",
    "results/tables/integrated_project_outputs_index.tsv",
    "results/tables/level1_feature_matrix.tsv",
    "results/tables/level2_group_definitions.tsv",
    "results/tables/level3_limited_vaf_clustering_qc_summary.tsv",
    "results/tables/level3_limited_vaf_clonal_complexity_by_sample.tsv",
    "results/tables/level3_allele_specific_cn_readiness_qc_summary.tsv",
    "results/reference_atlas/tcga_sample_reference_matrix.tsv.gz",
    "results/clinical_projection/example_cup_case/clinical_tcga_projection_report.md",
    "results/clinical_projection/example_cup_case/tcga_nearest_project_summary.tsv",
    "results/clinical_projection/example_cup_case/cup_interpretation/cup_molecular_interpretation_report.md",
    "results/clinical_projection/example_cup_case/cup_interpretation/cup_ranked_lineage_interpretation.tsv",
}
indexed_files = set(outputs_index["file_path"].astype(str))
missing_index_entries = expected_files - indexed_files
if missing_index_entries:
    fail(f"Output index lacks expected files: {', '.join(sorted(missing_index_entries))}")

exists_map = dict(zip(outputs_index["file_path"].astype(str), outputs_index["exists"].astype(str), strict=False))
for file_path in expected_files:
    if exists_map[file_path].lower() != "true":
        fail(f"Expected output index file is not marked as existing: {file_path}")

expected_levels = {
    "level1",
    "level2",
    "level3",
    "level3_readiness",
    "translational",
    "clinical_projection",
    "cup_interpretation",
    "integrated",
}
if not expected_levels.issubset(set(outputs_index["level"].astype(str))):
    fail("Output index lacks one or more required analysis/output levels")

atlas_rows = summary[
    (summary["level"].astype(str) == "translational_reference_atlas")
    & (summary["metric"].astype(str) == "number_of_retained_features")
]
if atlas_rows.empty or str(atlas_rows.iloc[0]["value"]) != "97":
    fail("Integrated project summary does not record the current 97-feature atlas")

print("Integrated project report output checks passed")

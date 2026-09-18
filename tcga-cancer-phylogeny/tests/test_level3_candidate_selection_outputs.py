#!/usr/bin/env python

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


paths = {
    "candidates": ROOT / "results" / "tables" / "level3_candidate_samples.tsv",
    "summary": ROOT / "results" / "tables" / "level3_candidate_summary_by_project.tsv",
    "pilot": ROOT / "results" / "tables" / "level3_pilot_cohort.tsv",
    "excluded": ROOT / "results" / "tables" / "level3_excluded_samples.tsv",
    "qc": ROOT / "results" / "tables" / "level3_candidate_selection_qc_summary.tsv",
    "mutation_purity_figure": ROOT / "results" / "figures" / "level3_candidate_mutation_count_vs_purity.pdf",
    "counts_figure": ROOT / "results" / "figures" / "level3_candidate_counts_by_project.pdf",
    "score_figure": ROOT / "results" / "figures" / "level3_candidate_score_by_project.pdf",
    "pilot_figure": ROOT / "results" / "figures" / "level3_pilot_cohort_overview.pdf",
}

missing_or_empty = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
if missing_or_empty:
    fail(f"Missing or empty Level 3 candidate-selection outputs: {', '.join(missing_or_empty)}")

candidates = pd.read_csv(paths["candidates"], sep="\t")
summary = pd.read_csv(paths["summary"], sep="\t")
pilot = pd.read_csv(paths["pilot"], sep="\t")
excluded = pd.read_csv(paths["excluded"], sep="\t")
qc = pd.read_csv(paths["qc"], sep="\t")

if candidates.empty:
    fail("Candidate table is empty")
if summary.empty:
    fail("Project summary is empty")
if excluded.empty:
    fail("Exclusion table is empty")
if qc.empty:
    fail("QC summary is empty")

required_candidate_cols = {
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "nonsynonymous_count",
    "purity",
    "ploidy",
    "ref_alt_counts_available",
    "candidate_class",
    "candidate_score",
    "inclusion_status",
    "exclusion_reason",
}
if not required_candidate_cols.issubset(candidates.columns):
    fail("Candidate table lacks required columns")

if not candidates["sample_barcode"].astype(str).str.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]").all():
    fail("Candidate sample barcodes are not valid TCGA sample barcodes")
if pd.to_numeric(candidates["candidate_score"], errors="coerce").isna().any():
    fail("candidate_score is not numeric for all candidates")
if candidates["inclusion_status"].isna().any() or (candidates["inclusion_status"].astype(str) == "").any():
    fail("inclusion_status is not populated")

excluded_reasons = excluded["exclusion_reason"].fillna("").astype(str)
if (excluded_reasons == "").any() or (excluded_reasons == "none").any():
    fail("Excluded samples lack populated exclusion reasons")

eligible_count = candidates["candidate_class"].isin(
    ["copy_number_aware_clonal_phylogeny_candidate", "clonal_clustering_candidate_limited_cn"]
).sum()
if eligible_count > 0 and pilot.empty:
    fail("Pilot cohort is empty despite eligible samples")
if not pilot.empty and not set(pilot["sample_barcode"]).issubset(set(candidates["sample_barcode"])):
    fail("Pilot cohort contains samples absent from candidate table")

required_qc_metrics = {
    "total_mutation_layer_samples",
    "samples_with_ref_alt_counts",
    "samples_with_purity",
    "samples_with_ploidy",
    "samples_with_copy_number_or_aneuploidy",
    "samples_passing_minimum_filters",
    "samples_selected_for_pilot_cohort",
}
if not required_qc_metrics.issubset(set(qc["metric"])):
    fail("QC summary lacks required metrics")

print("Level 3 candidate-selection output checks passed")

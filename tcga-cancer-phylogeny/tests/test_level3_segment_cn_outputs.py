#!/usr/bin/env python

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


paths = {
    "segments_gz": ROOT / "data" / "processed" / "copy_number" / "segments_by_sample.tsv.gz",
    "segment_summary": ROOT / "data" / "processed" / "copy_number" / "segment_level_cn_summary_by_sample.tsv",
    "candidates": ROOT / "results" / "tables" / "level3_candidate_samples_with_segments.tsv",
    "pilot": ROOT / "results" / "tables" / "level3_pilot_cohort_with_segments.tsv",
    "qc": ROOT / "results" / "tables" / "level3_segment_cn_qc_summary.tsv",
    "availability": ROOT / "results" / "tables" / "level3_segment_cn_availability_by_project.tsv",
    "availability_figure": ROOT / "results" / "figures" / "level3_segment_cn_availability_by_project.pdf",
    "pilot_coverage_figure": ROOT / "results" / "figures" / "level3_segment_cn_pilot_coverage.pdf",
    "upgrade_by_project_figure": ROOT / "results" / "figures" / "level3_segment_cn_upgrade_by_project.pdf",
    "count_figure": ROOT / "results" / "figures" / "level3_segment_count_distribution.pdf",
    "upgrade_figure": ROOT / "results" / "figures" / "level3_candidate_upgrade_status.pdf",
}

missing_or_empty = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
if missing_or_empty:
    fail(f"Missing or empty Level 3 segment-CN outputs: {', '.join(missing_or_empty)}")

candidates = pd.read_csv(paths["candidates"], sep="\t")
pilot = pd.read_csv(paths["pilot"], sep="\t")
qc = pd.read_csv(paths["qc"], sep="\t")
availability = pd.read_csv(paths["availability"], sep="\t")
segments = pd.read_csv(paths["segments_gz"], sep="\t")
summary = pd.read_csv(paths["segment_summary"], sep="\t")

if candidates.empty:
    fail("Candidate-with-segments table is empty")
if pilot.empty:
    fail("Pilot-with-segments table is empty")
if qc.empty:
    fail("Segment CN QC summary is empty")
if availability.empty:
    fail("Segment CN availability table is empty")

required_candidate_cols = {"sample_barcode", "candidate_class", "has_segment_level_cn", "inclusion_status", "notes"}
if not required_candidate_cols.issubset(candidates.columns):
    fail("Candidate-with-segments table lacks required columns")
if not candidates["sample_barcode"].astype(str).str.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]").all():
    fail("Malformed sample barcode detected")

segment_metrics = {
    "total_segment_rows_processed",
    "samples_with_segment_level_cn",
    "candidates_upgraded_to_copy_number_aware",
    "pilot_samples_total",
    "pilot_samples_with_sample_level_match",
    "pilot_samples_with_selected_best_file",
    "pilot_samples_still_lacking_segment_cn",
    "pilot_missing_reason_no_sample_level_match",
    "pilot_missing_reason_download_failed",
    "pilot_missing_reason_parse_failed",
    "pilot_missing_reason_file_already_absent",
    "pilot_missing_reason_patient_level_only",
    "pilot_missing_reason_other",
}
if not segment_metrics.issubset(set(qc["metric"])):
    fail("Segment CN QC summary lacks required metrics")

metric_values = dict(zip(qc["metric"], qc["value"].astype(str), strict=False))
segment_rows = int(metric_values["total_segment_rows_processed"])
upgraded = int(metric_values["candidates_upgraded_to_copy_number_aware"])
if segment_rows > 0:
    if summary.empty:
        fail("Segment rows exist but segment summary is empty")
    numeric_cols = [
        "n_segments",
        "n_autosomal_segments",
        "median_segment_mean",
        "sd_segment_mean",
        "max_abs_segment_mean",
        "fraction_segments_gain_like",
        "fraction_segments_loss_like",
    ]
    for col in numeric_cols:
        if pd.to_numeric(summary[col], errors="coerce").isna().all():
            fail(f"Segment summary numeric field is not numeric: {col}")
    if upgraded <= 0:
        fail("Segment data exist but no candidates were upgraded")
else:
    limitations = metric_values.get("limitations", "")
    if "No parseable segment-level CN rows" not in limitations:
        fail("No-source fallback did not record explicit manual-download limitation")
    if upgraded != 0:
        fail("No segment data available but candidates were upgraded")

valid_classes = {
    "copy_number_aware_clonal_phylogeny_candidate",
    "clonal_clustering_candidate_limited_cn",
    "descriptive_only",
    "excluded",
}
if not set(candidates["candidate_class"]).issubset(valid_classes):
    fail("Invalid candidate_class after segment CN join")

print("Level 3 segment-CN output checks passed")

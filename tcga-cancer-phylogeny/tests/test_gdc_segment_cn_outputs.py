#!/usr/bin/env python

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


paths = {
    "manifest": ROOT / "results" / "tables" / "gdc_segment_cn_manifest_pilot.tsv",
    "best_files": ROOT / "results" / "tables" / "gdc_segment_cn_best_file_per_pilot_sample.tsv",
    "download_log": ROOT / "results" / "tables" / "gdc_segment_cn_download_log.tsv",
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
    fail(f"Missing or empty GDC segment-CN outputs: {', '.join(missing_or_empty)}")

manifest = pd.read_csv(paths["manifest"], sep="\t")
best_files = pd.read_csv(paths["best_files"], sep="\t")
download_log = pd.read_csv(paths["download_log"], sep="\t")
candidates = pd.read_csv(paths["candidates"], sep="\t")
pilot = pd.read_csv(paths["pilot"], sep="\t")
qc = pd.read_csv(paths["qc"], sep="\t")
availability = pd.read_csv(paths["availability"], sep="\t")
summary = pd.read_csv(paths["segment_summary"], sep="\t")
segments = pd.read_csv(paths["segments_gz"], sep="\t")

required_manifest_cols = {
    "gdc_file_id",
    "file_name",
    "data_category",
    "data_type",
    "access",
    "project_id",
    "sample_submitter_ids",
    "aliquot_barcodes",
    "matched_sample_barcode",
    "matched_patient_barcode",
    "match_level",
    "download_status",
    "local_path",
    "parse_status",
    "notes",
}
if not required_manifest_cols.issubset(manifest.columns):
    fail("GDC manifest lacks required columns")

required_best_cols = {
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "selected_gdc_file_id",
    "selected_file_name",
    "selected_data_type",
    "match_level",
    "file_size",
    "n_candidate_files",
    "selection_reason",
    "selection_rank",
    "notes",
}
if not required_best_cols.issubset(best_files.columns):
    fail("Best-file-per-pilot-sample table lacks required columns")
if best_files["sample_barcode"].duplicated().any():
    fail("Best-file table should contain at most one selected file per pilot sample")
if not set(best_files["match_level"].dropna()).issubset({"sample", "patient"}):
    fail("Best-file table contains invalid match levels")

required_download_cols = {
    "sample_barcode",
    "gdc_file_id",
    "file_name",
    "download_attempted",
    "download_status",
    "local_path",
    "error_message",
    "file_size",
}
if not required_download_cols.issubset(download_log.columns):
    fail("GDC download log lacks required columns")

if candidates.empty or pilot.empty or qc.empty or availability.empty:
    fail("One or more GDC segment-CN output tables are unexpectedly empty")

if len(pilot) != 84:
    fail("Pilot-with-segments table should preserve 84 Level 3 pilot samples")
if not pilot["sample_barcode"].astype(str).str.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]").all():
    fail("Malformed TCGA sample barcode detected in pilot-with-segments table")

candidate_required = {"sample_barcode", "candidate_class", "has_segment_level_cn", "segment_cn_match_level", "notes"}
if not candidate_required.issubset(candidates.columns):
    fail("Candidate-with-segments table lacks GDC segment annotations")

availability_required = {
    "project_id",
    "project_code",
    "n_level3_candidates",
    "n_pilot_samples",
    "n_pilot_samples_with_segments",
    "n_gdc_files_found_for_project",
    "n_gdc_files_matched_to_pilot",
}
if not availability_required.issubset(availability.columns):
    fail("GDC availability-by-project table lacks required fields")

required_metrics = {
    "gdc_query_filters",
    "gdc_segment_files_found",
    "gdc_segment_files_matched_to_pilot",
    "gdc_unique_files_downloaded",
    "gdc_unique_files_parsed_successfully",
    "total_segment_rows_processed",
    "pilot_samples_with_segment_level_cn",
    "candidates_upgraded_to_copy_number_aware",
    "sample_level_manifest_matches",
    "patient_level_manifest_matches_lower_confidence",
    "pilot_samples_total",
    "pilot_samples_with_any_gdc_match",
    "pilot_samples_with_sample_level_match",
    "pilot_samples_with_selected_best_file",
    "pilot_samples_downloaded",
    "pilot_samples_parsed",
    "pilot_samples_upgraded",
    "pilot_samples_still_lacking_segment_cn",
    "pilot_missing_reason_no_sample_level_match",
    "pilot_missing_reason_download_failed",
    "pilot_missing_reason_parse_failed",
    "pilot_missing_reason_file_already_absent",
    "pilot_missing_reason_patient_level_only",
    "pilot_missing_reason_other",
    "pilot_projects_with_no_matching_segment_files",
    "limitations",
}
if not required_metrics.issubset(set(qc["metric"])):
    fail("GDC segment-CN QC summary lacks required metrics")

metric_values = dict(zip(qc["metric"], qc["value"].astype(str), strict=False))
try:
    found = int(metric_values["gdc_segment_files_found"])
    matched = int(metric_values["gdc_segment_files_matched_to_pilot"])
    downloaded = int(metric_values["gdc_unique_files_downloaded"])
    parsed_files = int(metric_values["gdc_unique_files_parsed_successfully"])
    segment_rows = int(metric_values["total_segment_rows_processed"])
    pilot_with_segments = int(metric_values["pilot_samples_with_segment_level_cn"])
    upgraded = int(metric_values["candidates_upgraded_to_copy_number_aware"])
    selected_best = int(metric_values["pilot_samples_with_selected_best_file"])
except ValueError as exc:
    fail(f"QC metric expected to be integer was not parseable: {exc}")

if matched > found:
    fail("GDC matched-file count cannot exceed found-file count")
if parsed_files > downloaded + int(metric_values.get("gdc_unique_files_already_present", "0")):
    fail("Parsed-file count cannot exceed downloaded/existing matched files")
if selected_best != len(best_files):
    fail("QC selected-best-file count does not match best-file table rows")

with gzip.open(paths["segments_gz"], "rt") as handle:
    header = handle.readline().rstrip("\n").split("\t")
required_segment_cols = {
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "chromosome",
    "start",
    "end",
    "num_probes",
    "segment_mean",
    "source_file",
    "gdc_file_id",
    "match_level",
    "notes",
}
if not required_segment_cols.issubset(header):
    fail("Processed segment table lacks required GDC segment columns")

if segment_rows > 0:
    if segments.empty or summary.empty:
        fail("QC reports segment rows but processed segment outputs are empty")
    if pilot_with_segments <= 0:
        fail("Segment rows exist but no pilot samples were marked with segment-level CN")
    numeric_cols = ["start", "end", "segment_mean"]
    for col in numeric_cols:
        if pd.to_numeric(segments[col], errors="coerce").isna().any():
            fail(f"Processed segment field is not fully numeric: {col}")
    valid_match_levels = {"sample", "patient"}
    if not set(segments["match_level"]).issubset(valid_match_levels):
        fail("Processed segment rows contain invalid match levels")
    if upgraded <= 0:
        fail("Segment-level CN rows exist but no eligible candidates were upgraded")
else:
    if not summary.empty:
        fail("No segment rows were processed but summary table is not empty")
    limitations = metric_values.get("limitations", "")
    if "No usable GDC pilot segment-level CN rows were parsed" not in limitations:
        fail("No-source GDC fallback did not record explicit limitations/manual instructions")
    if upgraded != 0:
        fail("No segment data were parsed but candidates were upgraded")

if not set(manifest["match_level"].dropna()).issubset({"sample", "patient", "unmatched"}):
    fail("Manifest contains unexpected match levels")
if not re.search(r"Copy Number", metric_values["gdc_query_filters"]):
    fail("QC filters do not record copy-number segment query terms")

print("GDC segment-CN output checks passed")

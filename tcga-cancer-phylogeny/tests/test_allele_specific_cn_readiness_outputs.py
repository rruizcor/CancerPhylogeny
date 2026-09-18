#!/usr/bin/env python

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


availability_path = ROOT / "results" / "tables" / "level3_allele_specific_cn_availability_by_sample.tsv"
standardized_path = (
    ROOT
    / "data"
    / "processed"
    / "allele_specific_cn"
    / "allele_specific_segments_by_sample.tsv.gz"
)
mutation_readiness_path = (
    ROOT / "results" / "tables" / "level3_mutation_to_allele_specific_cn_readiness.tsv"
)
manifest_path = (
    ROOT / "results" / "tables" / "level3_allele_specific_cn_external_inference_manifest.tsv"
)
qc_path = ROOT / "results" / "tables" / "level3_allele_specific_cn_readiness_qc_summary.tsv"
report_path = ROOT / "results" / "reports" / "level3_allele_specific_cn_and_pyclone_plan.md"
figure_paths = [
    ROOT / "results" / "figures" / "level3_allele_specific_cn_readiness_by_project.pdf",
    ROOT / "results" / "figures" / "level3_allele_specific_cn_missing_inputs.pdf",
]

for path in [
    availability_path,
    standardized_path,
    mutation_readiness_path,
    manifest_path,
    qc_path,
    report_path,
    *figure_paths,
]:
    if not path.exists() or path.stat().st_size == 0:
        fail(f"Required ASCN readiness output missing or empty: {path}")

availability = pd.read_csv(availability_path, sep="\t")
required_availability_columns = {
    "sample_barcode",
    "has_segment_mean_cn",
    "has_allele_specific_cn",
    "has_total_cn",
    "has_major_cn",
    "has_minor_cn",
    "can_prepare_cn_aware_pyclone_input",
    "recommended_next_action",
}
if not required_availability_columns.issubset(availability.columns):
    fail("ASCN availability table lacks required columns")
if len(availability) != 84 or availability["sample_barcode"].nunique() != 84:
    fail("ASCN availability table does not contain all 84 unique pilot samples")
if int(availability["in_limited_vaf_pilot"].sum()) != 30:
    fail("ASCN availability table does not identify all 30 limited-VAF pilot samples")
if not availability["has_segment_mean_cn"].all():
    fail("Current ASCN availability table lost known segment-mean CN coverage")
if availability["has_allele_specific_cn"].any():
    fail("Current project falsely reports local allele-specific CN")
if availability[["has_major_cn", "has_minor_cn"]].any().any():
    fail("Current project falsely reports major/minor CN")
if availability["can_prepare_cn_aware_pyclone_input"].any():
    fail("Current project falsely reports PyClone-VI-ready samples")

with gzip.open(standardized_path, "rt") as handle:
    standardized = pd.read_csv(handle, sep="\t")
required_segment_columns = {
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "chromosome",
    "start",
    "end",
    "total_cn",
    "major_cn",
    "minor_cn",
    "loh_status",
    "source_file",
    "source_tool",
    "notes",
}
if not required_segment_columns.issubset(standardized.columns):
    fail("Empty standardized ASCN table lacks the required header")
if not standardized.empty:
    fail("Current standardized ASCN table should be empty because no local files were found")

mutation_readiness = pd.read_csv(mutation_readiness_path, sep="\t")
if len(mutation_readiness) != 84:
    fail("Mutation-to-ASCN readiness table does not contain 84 pilot samples")
if mutation_readiness["can_build_pyclone_vi_input"].any():
    fail("Current mutation readiness falsely reports PyClone-VI-ready samples")
if mutation_readiness["can_build_phylowgs_input"].any():
    fail("Current mutation readiness falsely reports PhyloWGS-ready samples")
if (mutation_readiness["n_mutations_with_allele_specific_cn"] != 0).any():
    fail("Current mutation readiness falsely reports mutation-to-ASCN overlaps")

manifest = pd.read_csv(manifest_path, sep="\t")
if len(manifest) != 84 or not (manifest["recommended_tool"] == "FACETS").all():
    fail("External ASCN inference manifest does not cover all pilot samples with the planned tool")
if not manifest["missing_required_inputs"].str.contains("tumor_bam", regex=False).all():
    fail("External ASCN inference manifest does not record missing tumor inputs")
if not manifest["missing_required_inputs"].str.contains("matched_normal", regex=False).all():
    fail("External ASCN inference manifest does not record missing matched-normal inputs")

qc = pd.read_csv(qc_path, sep="\t")
global_qc = qc[qc["qc_section"] == "global"].set_index("metric")["value"].to_dict()
expected_zero_metrics = [
    "candidate_local_allele_specific_cn_files_found",
    "samples_with_existing_allele_specific_cn",
    "samples_with_major_cn",
    "samples_with_minor_cn",
    "samples_ready_for_pyclone_vi",
    "samples_ready_for_phylowgs",
]
for metric in expected_zero_metrics:
    if int(global_qc[metric]) != 0:
        fail(f"ASCN QC metric should be zero for the current local data: {metric}")
if int(global_qc["pilot_samples_evaluated"]) != 84:
    fail("ASCN QC pilot sample count is incorrect")
if int(global_qc["samples_requiring_external_allele_specific_cn_inference"]) != 84:
    fail("ASCN QC external-inference count is incorrect")

report = report_path.read_text()
if "segment_mean is not allele-specific integer copy number" not in report:
    fail("ASCN planning report lacks the required segment-mean caveat")
if "Recommended next action" not in report:
    fail("ASCN planning report lacks a recommended next action")
if "does not run PyClone-VI or PhyloWGS" not in report:
    fail("ASCN planning report does not state that inference tools were not run")

for path in figure_paths:
    if path.read_bytes()[:5] != b"%PDF-":
        fail(f"ASCN readiness figure is not a valid PDF: {path}")

print("Allele-specific CN readiness output checks passed")

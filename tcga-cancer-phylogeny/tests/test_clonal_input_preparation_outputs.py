#!/usr/bin/env python

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


paths = {
    "copy_neutral": ROOT / "data" / "processed" / "clonal" / "level3_copy_neutral_vaf_clustering_input.tsv.gz",
    "segment_annotated": ROOT / "data" / "processed" / "clonal" / "level3_segment_annotated_clonal_input.tsv.gz",
    "readiness": ROOT / "results" / "tables" / "level3_clonal_input_readiness_by_sample.tsv",
    "pilot": ROOT / "results" / "tables" / "level3_limited_vaf_clustering_pilot_samples.tsv",
    "qc": ROOT / "results" / "tables" / "level3_clonal_input_preparation_qc_summary.tsv",
    "sample_fig": ROOT / "results" / "figures" / "level3_copy_neutral_candidate_counts_by_sample.pdf",
    "project_fig": ROOT / "results" / "figures" / "level3_copy_neutral_candidate_counts_by_project.pdf",
    "vaf_fig": ROOT / "results" / "figures" / "level3_vaf_distribution_copy_neutral_candidates.pdf",
    "readiness_fig": ROOT / "results" / "figures" / "level3_readiness_class_by_project.pdf",
}

missing_or_empty = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
if missing_or_empty:
    fail(f"Missing or empty clonal input preparation outputs: {', '.join(missing_or_empty)}")

with gzip.open(paths["copy_neutral"], "rt") as handle:
    copy_neutral_header = handle.readline().rstrip("\n").split("\t")
with gzip.open(paths["segment_annotated"], "rt") as handle:
    segment_header = handle.readline().rstrip("\n").split("\t")

required_copy_neutral = {
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "mutation_id",
    "chromosome",
    "position",
    "Hugo_Symbol",
    "Variant_Classification",
    "t_ref_count",
    "t_alt_count",
    "total_depth",
    "observed_vaf",
    "purity",
    "ploidy",
    "local_segment_mean",
    "local_cn_status",
    "depth_filter_pass",
    "vaf_filter_pass",
    "copy_neutral_filter_pass",
    "notes",
}
if not required_copy_neutral.issubset(copy_neutral_header):
    fail("Copy-neutral input lacks required columns")

required_segment = {
    "sample_barcode",
    "mutation_id",
    "observed_vaf",
    "local_segment_mean",
    "local_cn_status",
    "segment_annotation_label",
    "copy_number_aware_status",
    "notes",
}
if not required_segment.issubset(segment_header):
    fail("Segment-annotated input lacks required columns")

copy_neutral = pd.read_csv(paths["copy_neutral"], sep="\t")
segment_annotated = pd.read_csv(paths["segment_annotated"], sep="\t")
readiness = pd.read_csv(paths["readiness"], sep="\t")
pilot = pd.read_csv(paths["pilot"], sep="\t")
qc = pd.read_csv(paths["qc"], sep="\t")

if copy_neutral.empty:
    fail("Copy-neutral input is unexpectedly empty")
if segment_annotated.empty:
    fail("Segment-annotated input is unexpectedly empty")
if readiness.empty or qc.empty:
    fail("Readiness or QC table is unexpectedly empty")

if not copy_neutral["sample_barcode"].astype(str).str.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]").all():
    fail("Malformed sample barcode detected in copy-neutral input")

if not copy_neutral["local_segment_mean"].between(-0.15, 0.15, inclusive="both").all():
    fail("Copy-neutral input contains mutations outside neutral segment thresholds")
if not copy_neutral["notes"].astype(str).str.contains("limited_vaf_based_clonal_clustering").all():
    fail("Copy-neutral input is not labeled for limited VAF clustering")
if copy_neutral["notes"].astype(str).str.contains("copy_number_aware_phylogeny").any():
    fail("Copy-neutral input should not claim copy-number-aware phylogeny readiness")

vaf = pd.to_numeric(copy_neutral["observed_vaf"], errors="coerce")
if ((vaf <= 0) | (vaf > 1)).any():
    fail("Copy-neutral VAF values must be >0 and <=1")
if (pd.to_numeric(copy_neutral["total_depth"], errors="coerce") < 20).any():
    fail("Copy-neutral input contains rows below default depth threshold")

if not segment_annotated["copy_number_aware_status"].astype(str).str.contains("requires_allele_specific_cn|copy_number_aware_candidate").all():
    fail("Segment-annotated input lacks copy-number-aware status labels")

required_readiness = {
    "sample_barcode",
    "n_total_mutations",
    "n_candidate_mutations_with_local_cn",
    "n_copy_neutral_candidate_mutations",
    "eligible_for_limited_vaf_clustering",
    "eligible_for_copy_number_aware_clustering",
    "reason_not_copy_number_aware",
    "readiness_class",
}
if not required_readiness.issubset(readiness.columns):
    fail("Readiness table lacks required columns")
if readiness["eligible_for_copy_number_aware_clustering"].astype(str).str.lower().isin({"true", "1"}).any():
    fail("Copy-number-aware eligibility should currently remain false without major/minor CN")
if not readiness["reason_not_copy_number_aware"].astype(str).str.contains("requires_allele_specific_cn|major_minor_cn_available").all():
    fail("Readiness table does not record allele-specific CN requirement")

valid_classes = {
    "ready_for_limited_vaf_clustering",
    "insufficient_neutral_mutations",
    "requires_allele_specific_cn",
    "insufficient_depth",
    "excluded",
}
if not set(readiness["readiness_class"].dropna()).issubset(valid_classes):
    fail("Unexpected readiness class detected")

if not pilot.empty:
    if pilot.groupby("project_code")["sample_barcode"].nunique().max() > 3:
        fail("Selected limited VAF pilot should respect default per-project cap")
    if len(pilot) > 30:
        fail("Selected limited VAF pilot should respect default maximum sample count")

required_qc = {
    "total_pilot_samples_evaluated",
    "total_annotated_mutations",
    "total_copy_neutral_candidate_mutations",
    "samples_eligible_for_limited_vaf_clustering",
    "samples_eligible_for_copy_number_aware_clustering",
    "samples_requiring_allele_specific_cn",
    "median_copy_neutral_candidate_mutations_per_sample",
    "projects_represented_among_limited_vaf_candidates",
    "projects_represented_in_selected_limited_vaf_pilot",
    "limitations",
}
if not required_qc.issubset(set(qc["metric"])):
    fail("Clonal input preparation QC summary lacks required metrics")

metric_values = dict(zip(qc["metric"], qc["value"].astype(str), strict=False))
if int(float(metric_values["samples_eligible_for_copy_number_aware_clustering"])) != 0:
    fail("QC should report zero copy-number-aware eligible samples at this stage")
if "requires_allele_specific_cn" not in metric_values["limitations"]:
    fail("QC limitations should record allele-specific CN requirement")

print("Clonal input preparation output checks passed")

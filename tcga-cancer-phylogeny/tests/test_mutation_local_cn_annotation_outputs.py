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
    "annotated": ROOT / "data" / "processed" / "clonal" / "level3_mutations_with_local_cn.tsv.gz",
    "summary": ROOT / "results" / "tables" / "level3_mutation_local_cn_summary_by_sample.tsv",
    "preview": ROOT / "data" / "processed" / "clonal" / "level3_pyclone_input_preview.tsv",
    "qc": ROOT / "results" / "tables" / "level3_mutation_local_cn_annotation_qc_summary.tsv",
    "coverage_figure": ROOT / "results" / "figures" / "level3_mutation_cn_annotation_coverage_by_sample.pdf",
    "vaf_figure": ROOT / "results" / "figures" / "level3_mutation_vaf_vs_local_segment_mean.pdf",
    "depth_figure": ROOT / "results" / "figures" / "level3_mutation_depth_distribution.pdf",
}

missing_or_empty = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
if missing_or_empty:
    fail(f"Missing or empty mutation local-CN outputs: {', '.join(missing_or_empty)}")

with gzip.open(paths["annotated"], "rt") as handle:
    header = handle.readline().rstrip("\n").split("\t")

required_annotated_cols = {
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "chromosome",
    "position",
    "start_position",
    "end_position",
    "Hugo_Symbol",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "t_ref_count",
    "t_alt_count",
    "total_depth",
    "observed_vaf",
    "purity",
    "ploidy",
    "local_cn_segment_chromosome",
    "local_cn_segment_start",
    "local_cn_segment_end",
    "local_cn_num_probes",
    "local_segment_mean",
    "local_cn_match_status",
    "local_cn_match_notes",
}
if not required_annotated_cols.issubset(header):
    fail("Annotated mutation table lacks required columns")

annotated = pd.read_csv(paths["annotated"], sep="\t")
summary = pd.read_csv(paths["summary"], sep="\t")
preview = pd.read_csv(paths["preview"], sep="\t")
qc = pd.read_csv(paths["qc"], sep="\t")

if annotated.empty:
    fail("Annotated mutation table is unexpectedly empty")
if summary.empty or qc.empty:
    fail("Summary or QC output is unexpectedly empty")

if not annotated["sample_barcode"].astype(str).str.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]").all():
    fail("Malformed sample barcode detected in annotated mutation table")

vaf = pd.to_numeric(annotated["observed_vaf"], errors="coerce").dropna()
if ((vaf < 0) | (vaf > 1)).any():
    fail("Observed VAF values must be between 0 and 1")

valid_status = {"matched", "unmatched"}
if not set(annotated["local_cn_match_status"].dropna().astype(str)).issubset(valid_status):
    fail("Unexpected local CN match status detected")

matched = annotated[annotated["local_cn_match_status"] == "matched"]
if matched.empty:
    fail("Segment data exist but no mutations were annotated with local CN")
for column in ["local_cn_segment_start", "local_cn_segment_end", "local_segment_mean"]:
    if pd.to_numeric(matched[column], errors="coerce").isna().any():
        fail(f"Matched mutations have missing segment field: {column}")

unmatched = annotated[annotated["local_cn_match_status"] == "unmatched"]
if not unmatched.empty and unmatched["local_cn_match_notes"].astype(str).eq("").any():
    fail("Unmatched mutations should include an explanatory match note")

required_summary_cols = {
    "sample_barcode",
    "n_total_mutations",
    "n_mutations_with_local_cn",
    "pct_mutations_with_local_cn",
    "eligible_for_copy_number_aware_clonal_input",
    "exclusion_reason_or_warning",
}
if not required_summary_cols.issubset(summary.columns):
    fail("Sample-level mutation/CN summary lacks required columns")
if summary["eligible_for_copy_number_aware_clonal_input"].isna().any():
    fail("Eligibility flag is not populated for all pilot samples")

required_qc_metrics = {
    "pilot_samples_evaluated",
    "pilot_samples_with_segment_cn",
    "pilot_samples_with_mutations",
    "pilot_samples_with_mutation_to_segment_annotation",
    "total_mutations_evaluated",
    "total_mutations_annotated_with_local_cn",
    "percentage_mutations_annotated",
    "samples_eligible_for_downstream_clonal_input",
    "chromosome_naming_issues",
    "coordinate_mismatch_issues",
    "non_allele_specific_cn_warning",
}
if not required_qc_metrics.issubset(set(qc["metric"])):
    fail("Mutation local-CN QC summary lacks required metrics")

if not preview.empty:
    required_preview_cols = {
        "sample_id",
        "mutation_id",
        "ref_counts",
        "var_counts",
        "normal_cn",
        "major_cn",
        "minor_cn",
        "local_segment_mean",
        "purity",
        "notes",
    }
    if not required_preview_cols.issubset(preview.columns):
        fail("PyClone-style preview lacks required columns")
    if preview["major_cn"].notna().any() or preview["minor_cn"].notna().any():
        fail("PyClone preview should not invent major/minor allele-specific CN")
    if not preview["notes"].astype(str).str.contains("preview_only|segment_mean_not_allele_specific", regex=True).all():
        fail("PyClone preview notes should label the table as non-final")

metric_values = dict(zip(qc["metric"], qc["value"].astype(str), strict=False))
if int(float(metric_values["total_mutations_annotated_with_local_cn"])) > len(annotated):
    fail("QC annotated mutation count exceeds annotated table rows")
if not re.search("allele-specific", metric_values["non_allele_specific_cn_warning"]):
    fail("QC summary should warn that segment means are not allele-specific CN")

print("Mutation local-CN annotation output checks passed")

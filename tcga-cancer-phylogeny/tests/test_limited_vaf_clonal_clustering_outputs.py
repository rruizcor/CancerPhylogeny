#!/usr/bin/env python

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


paths = {
    "assignments": ROOT / "results" / "clonal" / "limited_vaf" / "level3_limited_vaf_cluster_assignments.tsv.gz",
    "sample_summary": ROOT / "results" / "tables" / "level3_limited_vaf_clonal_complexity_by_sample.tsv",
    "project_summary": ROOT / "results" / "tables" / "level3_limited_vaf_clonal_complexity_by_project.tsv",
    "qc": ROOT / "results" / "tables" / "level3_limited_vaf_clustering_qc_summary.tsv",
    "density_fig": ROOT / "results" / "figures" / "level3_limited_vaf_cluster_vaf_density_by_sample.pdf",
    "counts_fig": ROOT / "results" / "figures" / "level3_limited_vaf_cluster_counts_by_project.pdf",
    "heatmap_fig": ROOT / "results" / "figures" / "level3_limited_vaf_complexity_heatmap.pdf",
    "examples_fig": ROOT / "results" / "figures" / "level3_limited_vaf_example_samples.pdf",
}

missing_or_empty = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
if missing_or_empty:
    fail(f"Missing or empty limited VAF clustering outputs: {', '.join(missing_or_empty)}")

with gzip.open(paths["assignments"], "rt") as handle:
    assignment_header = handle.readline().rstrip("\n").split("\t")

required_assignment = {
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
    "cluster_id",
    "cluster_mean_vaf",
    "cluster_median_vaf",
    "cluster_size",
    "cluster_interpretation",
    "notes",
}
if not required_assignment.issubset(assignment_header):
    fail("Cluster assignment table lacks required columns")

assignments = pd.read_csv(paths["assignments"], sep="\t")
sample_summary = pd.read_csv(paths["sample_summary"], sep="\t")
project_summary = pd.read_csv(paths["project_summary"], sep="\t")
qc = pd.read_csv(paths["qc"], sep="\t")

if assignments.empty:
    fail("Cluster assignment table is unexpectedly empty")
if sample_summary.empty or project_summary.empty or qc.empty:
    fail("One or more limited VAF summary tables are unexpectedly empty")

vaf = pd.to_numeric(assignments["observed_vaf"], errors="coerce")
if vaf.isna().any() or ((vaf <= 0) | (vaf > 1)).any():
    fail("Observed VAF values must be populated and between 0 and 1")
if assignments["cluster_id"].isna().any() or assignments["cluster_id"].astype(str).eq("").any():
    fail("Cluster IDs must be populated")
valid_cluster_labels = {"clonal_like", "subclonal_like", "uncertain"}
if not set(assignments["cluster_interpretation"].dropna()).issubset(valid_cluster_labels):
    fail("Unexpected cluster interpretation label")
if assignments["notes"].astype(str).str.contains("definitive|full_copy_number_aware", case=False, regex=True).any():
    fail("Assignment notes should not claim definitive/full CN-aware phylogeny")
if not assignments["notes"].astype(str).str.contains("limited_vaf_based_clonal_clustering").all():
    fail("Assignments should be labeled as limited VAF-based clustering")

required_sample = {
    "sample_barcode",
    "n_input_mutations",
    "n_clusters",
    "dominant_cluster_id",
    "dominant_cluster_mean_vaf",
    "dominant_cluster_size",
    "dominant_cluster_fraction",
    "subclonal_cluster_count",
    "subclonal_mutation_fraction",
    "vaf_entropy",
    "vaf_dispersion",
    "interpretation_class",
    "warnings",
}
if not required_sample.issubset(sample_summary.columns):
    fail("Sample complexity summary lacks required columns")
valid_classes = {"predominantly_clonal_like", "oligoclonal_like", "multicluster_subclonal_like", "low_confidence"}
if not set(sample_summary["interpretation_class"].dropna()).issubset(valid_classes):
    fail("Unexpected sample interpretation class")

required_project = {
    "project_code",
    "n_samples",
    "median_n_input_mutations",
    "median_n_clusters",
    "median_dominant_cluster_fraction",
    "median_subclonal_mutation_fraction",
    "median_vaf_entropy",
    "interpretation_summary",
}
if not required_project.issubset(project_summary.columns):
    fail("Project summary lacks required columns")

required_qc = {
    "total_selected_pilot_samples",
    "samples_successfully_clustered",
    "samples_skipped",
    "total_input_mutations",
    "total_clustered_mutations",
    "clustering_method_used",
    "parameter_settings",
    "limitations",
}
if not required_qc.issubset(set(qc["metric"])):
    fail("Limited VAF QC summary lacks required metrics")
metric_values = dict(zip(qc["metric"], qc["value"].astype(str), strict=False))
if int(float(metric_values["total_clustered_mutations"])) != len(assignments):
    fail("QC clustered mutation count does not match assignments")
if "not_full_copy_number_aware_phylogeny" not in metric_values["limitations"]:
    fail("QC limitations must state this is not full CN-aware phylogeny")
if "definitive" in " ".join(qc["value"].astype(str)).lower().replace("does_not_infer_definitive", ""):
    fail("QC should not claim definitive branching")

print("Limited VAF clonal clustering output checks passed")

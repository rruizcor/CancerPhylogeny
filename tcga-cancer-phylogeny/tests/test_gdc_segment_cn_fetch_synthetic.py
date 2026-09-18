#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "09c_fetch_gdc_segment_cn.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("gdc_segment_cn", SCRIPT)
gdc_segment_cn = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["gdc_segment_cn"] = gdc_segment_cn
spec.loader.exec_module(gdc_segment_cn)


def fail(message: str) -> None:
    raise AssertionError(message)


pilot = pd.DataFrame(
    {
        "sample_barcode": ["TCGA-AA-0001-01A", "TCGA-BB-0001-01A"],
        "patient_barcode": ["TCGA-AA-0001", "TCGA-BB-0001"],
        "project_id": ["TCGA-SKCM", "TCGA-UCEC"],
        "project_code": ["SKCM", "UCEC"],
        "candidate_class": [
            "clonal_clustering_candidate_limited_cn",
            "clonal_clustering_candidate_limited_cn",
        ],
        "inclusion_status": ["selected_for_pilot", "selected_for_pilot"],
    }
)
candidates = pilot.assign(
    ref_alt_counts_available=True,
    purity_available=True,
    ploidy_available=True,
    copy_number_available=True,
    has_segment_level_cn=False,
    exclusion_reason="none",
    notes="single_bulk_tcga_sample",
)

hits = [
    {
        "file_id": "file-sample",
        "file_name": "TCGA-AA-0001-01A-11D.segment.tsv",
        "data_category": "Copy Number Variation",
        "data_type": "Masked Copy Number Segment",
        "experimental_strategy": "Genotyping Array",
        "access": "open",
        "file_size": 123,
        "cases": [
            {
                "submitter_id": "TCGA-AA-0001",
                "project": {"project_id": "TCGA-SKCM"},
                "samples": [{"submitter_id": "TCGA-AA-0001-01A-11D"}],
            }
        ],
    },
    {
        "file_id": "file-sample-filtered",
        "file_name": "TCGA-AA-0001-01A.filtered.segment.tsv",
        "data_category": "Copy Number Variation",
        "data_type": "Filtered Copy Number Segment",
        "experimental_strategy": "Genotyping Array",
        "access": "open",
        "file_size": 100,
        "cases": [
            {
                "submitter_id": "TCGA-AA-0001",
                "project": {"project_id": "TCGA-SKCM"},
                "samples": [{"submitter_id": "TCGA-AA-0001-01A"}],
            }
        ],
    },
    {
        "file_id": "file-sample-patient-lower-confidence",
        "file_name": "TCGA-AA-0001.patient.segment.tsv",
        "data_category": "Copy Number Variation",
        "data_type": "Masked Copy Number Segment",
        "experimental_strategy": "Genotyping Array",
        "access": "open",
        "file_size": 50,
        "cases": [
            {
                "submitter_id": "TCGA-AA-0001",
                "project": {"project_id": "TCGA-SKCM"},
                "samples": [{"submitter_id": "TCGA-AA-0001-10A"}],
            }
        ],
    },
    {
        "file_id": "file-patient",
        "file_name": "patient_level.segment.tsv",
        "data_category": "Copy Number Variation",
        "data_type": "Copy Number Segment",
        "experimental_strategy": "Genotyping Array",
        "access": "open",
        "file_size": 456,
        "cases": [
            {
                "submitter_id": "TCGA-BB-0001",
                "project": {"project_id": "TCGA-UCEC"},
                "samples": [{"submitter_id": "TCGA-BB-0001-10A"}],
            }
        ],
    },
    {
        "file_id": "file-unmatched",
        "file_name": "TCGA-CC-0001-01A.segment.tsv",
        "data_category": "Copy Number Variation",
        "data_type": "Filtered Copy Number Segment",
        "experimental_strategy": "Genotyping Array",
        "access": "open",
        "file_size": 789,
        "cases": [
            {
                "submitter_id": "TCGA-CC-0001",
                "project": {"project_id": "TCGA-LUAD"},
                "samples": [{"submitter_id": "TCGA-CC-0001-01A"}],
            }
        ],
    },
]

manifest = gdc_segment_cn.manifest_from_hits(hits, pilot)
if len(manifest) != 5:
    fail("Synthetic GDC manifest should have one row per file/match")
levels = dict(zip(manifest["gdc_file_id"], manifest["match_level"], strict=False))
if levels["file-sample"] != "sample":
    fail("Sample-level GDC file was not matched at sample level")
if levels["file-sample-filtered"] != "sample":
    fail("Second sample-level GDC file was not matched at sample level")
if levels["file-sample-patient-lower-confidence"] != "patient":
    fail("Patient-level file for a sample with sample-level options was not retained as patient match")
if levels["file-patient"] != "patient":
    fail("Patient-level GDC file was not matched as lower-confidence patient match")
if levels["file-unmatched"] != "unmatched":
    fail("Unmatched GDC file should remain unmatched")
if manifest.loc[manifest["gdc_file_id"] == "file-sample", "aliquot_barcodes"].iloc[0] == "NA":
    fail("Aliquot barcode was not preserved in the manifest")

tmp = Path(tempfile.mkdtemp(prefix="gdc_segment_cn_"))
best_default = gdc_segment_cn.select_best_files_per_pilot_sample(manifest, pilot, candidates, tmp)
if len(best_default) != 1:
    fail("Default best-file selection should select only sample-level matched pilot samples")
selected_aa = best_default[best_default["sample_barcode"] == "TCGA-AA-0001-01A"].iloc[0]
if selected_aa["selected_gdc_file_id"] != "file-sample":
    fail("Best-file selection should prefer masked sample-level segment files over filtered or patient-level files")
best_with_patient = gdc_segment_cn.select_best_files_per_pilot_sample(
    manifest,
    pilot,
    candidates,
    tmp,
    {**gdc_segment_cn.CONFIG, "allow_patient_level_matches": True},
)
if len(best_with_patient) != 2:
    fail("Allowing patient-level matches should select a best file for the patient-only pilot sample")
capped = gdc_segment_cn.apply_download_mode(
    best_with_patient,
    {**gdc_segment_cn.CONFIG, "download_mode": "pilot_capped", "max_downloads": 1},
)
if len(capped) != 1:
    fail("pilot_capped mode should cap selected best files")

sample_file = tmp / "file-sample.tsv"
sample_file.write_text(
    "GDC_Aliquot\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    "TCGA-AA-0001-01A-11D\t1\t100\t200\t20\t0.30\n"
    "TCGA-AA-0001-01A-11D\t2\t250\t400\t10\t-0.40\n",
    encoding="utf-8",
)
expected_present = gdc_segment_cn.expected_local_path(tmp, "file-sample", "TCGA-AA-0001-01A-11D.segment.tsv")
expected_present.write_text(sample_file.read_text(encoding="utf-8"), encoding="utf-8")
patient_file = tmp / "file-patient.tsv"
patient_file.write_text(
    "GDC_Aliquot\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    "TCGA-BB-0001-10A\t3\t500\t700\t15\t0.10\n",
    encoding="utf-8",
)

download_manifest, download_log, download_warnings = gdc_segment_cn.download_selected_files(
    manifest,
    best_default,
    tmp,
    {**gdc_segment_cn.CONFIG, "download_mode": "pilot_full", "overwrite_existing": False},
)
if download_warnings:
    fail(f"Existing-file skip generated unexpected warnings: {download_warnings}")
if download_log.loc[download_log["gdc_file_id"] == "file-sample", "download_status"].iloc[0] != "already_present":
    fail("Existing selected file should be skipped when overwrite is false")

manifest.loc[manifest["gdc_file_id"] == "file-sample", "download_status"] = "already_present"
manifest.loc[manifest["gdc_file_id"] == "file-sample", "local_path"] = str(sample_file)
manifest.loc[manifest["gdc_file_id"] == "file-patient", "download_status"] = "already_present"
manifest.loc[manifest["gdc_file_id"] == "file-patient", "local_path"] = str(patient_file)

segments, parsed_manifest, warnings = gdc_segment_cn.parse_downloaded_segments(manifest, candidates)
if warnings:
    fail(f"Unexpected synthetic parse warnings: {warnings}")
if len(segments) != 3:
    fail("Downloaded synthetic segment files were not parsed correctly")
if not set(gdc_segment_cn.SEGMENT_COLUMNS).issubset(segments.columns):
    fail("Parsed segments lack required GDC segment columns")
if set(parsed_manifest[parsed_manifest["download_status"] == "already_present"]["parse_status"]) != {"parsed"}:
    fail("Synthetic manifest parse status was not updated")

summary = gdc_segment_cn.segment_summary(segments, candidates)
if len(summary) != 2:
    fail("Segment summary should contain two pilot samples")
aa = summary[summary["sample_barcode"] == "TCGA-AA-0001-01A"].iloc[0]
if aa["n_segments"] != 2 or aa["segment_cn_match_level"] != "sample":
    fail("Sample-level segment summary is incorrect")
bb = summary[summary["sample_barcode"] == "TCGA-BB-0001-01A"].iloc[0]
if bb["segment_cn_match_level"] != "patient":
    fail("Patient-level lower-confidence match was not retained in summary")

updated = gdc_segment_cn.update_candidates_with_segments(candidates, summary)
if int(updated["has_segment_level_cn"].sum()) != 2:
    fail("Candidates were not marked as having segment-level CN")
sample_level_class = updated.loc[updated["sample_barcode"] == "TCGA-AA-0001-01A", "candidate_class"].iloc[0]
patient_level_class = updated.loc[updated["sample_barcode"] == "TCGA-BB-0001-01A", "candidate_class"].iloc[0]
if sample_level_class != "copy_number_aware_clonal_phylogeny_candidate":
    fail("Eligible sample-level segment CN candidate was not upgraded")
if patient_level_class != "clonal_clustering_candidate_limited_cn":
    fail("Patient-level lower-confidence segment CN should not upgrade to full copy-number-aware class")
if "patient_level_segment_cn_lower_confidence" not in updated.loc[
    updated["sample_barcode"] == "TCGA-BB-0001-01A", "notes"
].iloc[0]:
    fail("Patient-level caveat was not added to candidate notes")

filters = gdc_segment_cn.build_gdc_filters(["TCGA-SKCM"], use_file_prefix=True, include_data_category=True)
qc_best = best_with_patient
qc_log = pd.DataFrame(
    [
        {
            "sample_barcode": "TCGA-AA-0001-01A",
            "gdc_file_id": "file-sample",
            "file_name": "TCGA-AA-0001-01A-11D.segment.tsv",
            "download_attempted": False,
            "download_status": "already_present",
            "local_path": str(sample_file),
            "error_message": "NA",
            "file_size": 123,
        },
        {
            "sample_barcode": "TCGA-BB-0001-01A",
            "gdc_file_id": "file-patient",
            "file_name": "patient_level.segment.tsv",
            "download_attempted": False,
            "download_status": "already_present",
            "local_path": str(patient_file),
            "error_message": "NA",
            "file_size": 456,
        },
    ]
)
qc = gdc_segment_cn.qc_summary(
    parsed_manifest,
    qc_best,
    qc_log,
    segments,
    summary,
    candidates,
    updated,
    pilot,
    gdc_segment_cn.update_pilot_with_segments(pilot, updated),
    filters,
    [],
)
required_metrics = {
    "gdc_query_filters",
    "gdc_segment_files_found",
    "gdc_segment_files_matched_to_pilot",
    "pilot_samples_with_selected_best_file",
    "pilot_missing_reason_no_sample_level_match",
}
if not required_metrics.issubset(set(qc["metric"])):
    fail("Synthetic QC summary lacks required GDC metrics")

print("Synthetic GDC segment-CN fetch checks passed")

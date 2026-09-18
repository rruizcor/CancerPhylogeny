#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "09b_download_segment_level_cn.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("level3_segment_cn", SCRIPT)
segment_cn = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["level3_segment_cn"] = segment_cn
spec.loader.exec_module(segment_cn)


def fail(message: str) -> None:
    raise AssertionError(message)


assert segment_cn.normalize_sample_barcode("TCGA-AA-0001-01A-01D-1234-01") == "TCGA-AA-0001-01A"
assert segment_cn.normalize_patient_barcode("TCGA-AA-0001-01A-01D-1234-01") == "TCGA-AA-0001"

tmp = Path(tempfile.mkdtemp(prefix="level3_segment_cn_"))
source_a = tmp / "synthetic_a.seg"
source_b = tmp / "synthetic_b.tsv"
source_a.write_text(
    "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    "TCGA-AA-0001-01A-01D-1111-01\t1\t100\t200\t20\t0.35\n"
    "TCGA-AA-0001-01A-01D-1111-01\t2\t300\t450\t15\t-0.30\n",
    encoding="utf-8",
)
source_b.write_text(
    "ID\tchrom\tloc.start\tloc.end\tseg.mean\n"
    "TCGA-BB-0001-01A\t3\t500\t900\t0.10\n"
    "TCGA-BB-0001-01A\tX\t100\t200\t-0.50\n",
    encoding="utf-8",
)

candidates = pd.DataFrame(
    {
        "sample_barcode": ["TCGA-AA-0001-01A", "TCGA-BB-0001-01A", "TCGA-CC-0001-01A"],
        "patient_barcode": ["TCGA-AA-0001", "TCGA-BB-0001", "TCGA-CC-0001"],
        "project_id": ["TCGA-SKCM", "TCGA-UCEC", "TCGA-LUAD"],
        "project_code": ["SKCM", "UCEC", "LUAD"],
        "nonsynonymous_count": [500, 300, 20],
        "ref_alt_counts_available": [True, True, True],
        "purity_available": [True, True, True],
        "ploidy_available": [True, True, True],
        "copy_number_available": [True, True, True],
        "has_segment_level_cn": [False, False, False],
        "candidate_class": [
            "clonal_clustering_candidate_limited_cn",
            "clonal_clustering_candidate_limited_cn",
            "excluded",
        ],
        "candidate_score": [80, 70, 20],
        "inclusion_status": ["eligible_not_selected", "selected_for_pilot", "excluded"],
        "exclusion_reason": ["none", "none", "too_few_mutations_lt_100"],
        "notes": ["single_bulk_tcga_sample", "single_bulk_tcga_sample", "single_bulk_tcga_sample"],
    }
)
lookup = segment_cn.project_lookup_from_candidates(candidates)
raw_a = segment_cn.read_segment_file(source_a)
parsed_a, warnings_a = segment_cn.standardize_segment_table(raw_a, source_a, lookup)
raw_b = segment_cn.read_segment_file(source_b)
parsed_b, warnings_b = segment_cn.standardize_segment_table(raw_b, source_b, lookup)
if warnings_a or warnings_b:
    fail(f"Unexpected synthetic parsing warnings: {warnings_a + warnings_b}")
if len(parsed_a) != 2 or len(parsed_b) != 2:
    fail("Synthetic segment files were not parsed correctly")
if parsed_a["sample_barcode"].iloc[0] != "TCGA-AA-0001-01A":
    fail("Aliquot barcode was not normalized to sample barcode")

segments = pd.concat([parsed_a, parsed_b], ignore_index=True)
summary = segment_cn.segment_summary(segments, candidates)
if len(summary) != 2:
    fail("Sample-level segment summary has wrong row count")
aa = summary[summary["sample_barcode"] == "TCGA-AA-0001-01A"].iloc[0]
if aa["n_segments"] != 2 or aa["n_autosomal_segments"] != 2:
    fail("Segment counts were not summarized correctly")
if abs(aa["fraction_segments_gain_like"] - 0.5) > 1e-8:
    fail("Gain-like segment fraction was not calculated correctly")

updated = segment_cn.update_candidates_with_segments(candidates, summary)
upgraded = updated[updated["sample_barcode"] == "TCGA-AA-0001-01A"].iloc[0]
if upgraded["candidate_class"] != "copy_number_aware_clonal_phylogeny_candidate":
    fail("Candidate with segment CN was not upgraded")
not_upgraded = updated[updated["sample_barcode"] == "TCGA-CC-0001-01A"].iloc[0]
if not_upgraded["candidate_class"] != "excluded":
    fail("Excluded sample should remain excluded")

empty_summary = segment_cn.segment_summary(segment_cn.empty_segments(), candidates)
fallback = segment_cn.update_candidates_with_segments(candidates, empty_summary)
if fallback["has_segment_level_cn"].sum() != 0:
    fail("Empty segment fallback incorrectly marked segment CN availability")

print("Synthetic Level 3 segment-CN checks passed")

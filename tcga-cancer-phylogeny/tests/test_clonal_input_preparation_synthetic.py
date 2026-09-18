#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "11_prepare_clonal_input_sets.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("clonal_input_sets", SCRIPT)
clonal_input_sets = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["clonal_input_sets"] = clonal_input_sets
spec.loader.exec_module(clonal_input_sets)


def fail(message: str) -> None:
    raise AssertionError(message)


def mutation_row(sample: str, project: str, pos: int, seg_mean: float, depth: int, alt: int, note: str = "sample_level_segment_overlap") -> dict:
    return {
        "sample_barcode": sample,
        "patient_barcode": sample[:12],
        "project_id": f"TCGA-{project}",
        "project_code": project,
        "chromosome": "1",
        "position": pos,
        "Hugo_Symbol": f"GENE{pos}",
        "Variant_Classification": "Missense_Mutation",
        "Variant_Type": "SNP",
        "Reference_Allele": "A",
        "Tumor_Seq_Allele2": "T",
        "HGVSp_Short": "p.X1Y",
        "t_ref_count": depth - alt,
        "t_alt_count": alt,
        "total_depth": depth,
        "observed_vaf": alt / depth if depth else 0,
        "purity": 0.70,
        "ploidy": 2.0,
        "local_segment_mean": seg_mean,
        "local_cn_status": "neutral_like",
        "local_cn_match_status": "matched",
        "local_cn_match_notes": note,
        "is_nonsynonymous": True,
    }


rows: list[dict] = []
for i in range(60):
    rows.append(mutation_row("TCGA-AA-0001-01A", "SKCM", 1000 + i, 0.02, 40, 10))
rows.append(mutation_row("TCGA-AA-0001-01A", "SKCM", 2000, 0.40, 40, 10))
rows.append(mutation_row("TCGA-AA-0001-01A", "SKCM", 2001, -0.40, 40, 10))
rows.append(mutation_row("TCGA-AA-0001-01A", "SKCM", 2002, 0.00, 10, 3))
rows.append(mutation_row("TCGA-AA-0001-01A", "SKCM", 2003, 0.00, 40, 0))
rows.append(mutation_row("TCGA-AA-0001-01A", "SKCM", 2004, 0.00, 40, 10, "patient_level_segment_overlap"))
for i in range(10):
    rows.append(mutation_row("TCGA-BB-0001-01A", "LUAD", 3000 + i, 0.01, 35, 7))
for i in range(55):
    rows.append(mutation_row("TCGA-CC-0001-01A", "SKCM", 4000 + i, 0.01, 35, 7))
for i in range(55):
    rows.append(mutation_row("TCGA-DD-0001-01A", "SKCM", 5000 + i, 0.01, 35, 7))
for i in range(55):
    rows.append(mutation_row("TCGA-EE-0001-01A", "SKCM", 6000 + i, 0.01, 35, 7))
annotated = pd.DataFrame(rows)

pilot = annotated[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates()
config = {
    **clonal_input_sets.CONFIG,
    "min_total_depth": 20,
    "neutral_segment_mean_min": -0.15,
    "neutral_segment_mean_max": 0.15,
    "min_copy_neutral_candidate_mutations": 50,
    "min_median_depth": 20,
    "pilot_per_project_cap": 3,
    "pilot_max_samples": 4,
}

flagged = clonal_input_sets.prepare_mutation_flags(annotated, config)
if clonal_input_sets.normalize_local_cn_status(0.0, config) != "neutral_like":
    fail("Neutral-like CN status threshold failed")
if clonal_input_sets.normalize_local_cn_status(0.20, config) != "gain_like":
    fail("Gain-like CN status threshold failed")
if clonal_input_sets.normalize_local_cn_status(-0.20, config) != "loss_like":
    fail("Loss-like CN status threshold failed")

copy_neutral = clonal_input_sets.build_copy_neutral_input(flagged)
aa = copy_neutral[copy_neutral["sample_barcode"] == "TCGA-AA-0001-01A"]
if len(aa) != 60:
    fail("Copy-neutral set should exclude gain/loss, low-depth, zero-VAF, and patient-level-like rows")
if not aa["notes"].astype(str).str.contains("limited_vaf_based_clonal_clustering").all():
    fail("Copy-neutral set must be honestly labeled as limited VAF clustering")
if not ((copy_neutral["observed_vaf"] > 0) & (copy_neutral["total_depth"] >= 20)).all():
    fail("Depth and VAF filters were not applied")

segment_annotated = clonal_input_sets.build_segment_annotated_input(flagged)
if len(segment_annotated) <= len(copy_neutral):
    fail("Segment-annotated set should retain non-neutral gain/loss candidates")
if not segment_annotated["copy_number_aware_status"].astype(str).str.contains("requires_allele_specific_cn").all():
    fail("Segment-annotated set should not claim copy-number-aware status without major/minor CN")

readiness = clonal_input_sets.build_readiness_table(flagged, pilot, config)
aa_ready = readiness[readiness["sample_barcode"] == "TCGA-AA-0001-01A"].iloc[0]
bb_ready = readiness[readiness["sample_barcode"] == "TCGA-BB-0001-01A"].iloc[0]
if not bool(aa_ready["eligible_for_limited_vaf_clustering"]):
    fail("Sample with enough neutral mutations should be eligible for limited VAF clustering")
if bool(aa_ready["eligible_for_copy_number_aware_clustering"]):
    fail("Copy-number-aware eligibility should remain false without major/minor CN")
if aa_ready["reason_not_copy_number_aware"] != "requires_allele_specific_cn":
    fail("Missing major/minor CN should be recorded as requiring allele-specific CN")
if bb_ready["readiness_class"] != "insufficient_neutral_mutations":
    fail("Sample with too few neutral candidates should be classified correctly")

selected = clonal_input_sets.select_limited_vaf_pilot(readiness, config)
skcm_selected = selected[selected["project_code"] == "SKCM"]
if len(skcm_selected) > 3:
    fail("Pilot selection should respect per-project cap")
if selected.empty or not selected["selection_reason"].astype(str).str.contains("project_diversity_cap_3").all():
    fail("Pilot selection reason should record diversity cap")

qc = clonal_input_sets.build_qc_summary(flagged, copy_neutral, readiness, selected)
required_metrics = {
    "total_pilot_samples_evaluated",
    "total_copy_neutral_candidate_mutations",
    "samples_eligible_for_limited_vaf_clustering",
    "samples_eligible_for_copy_number_aware_clustering",
    "samples_requiring_allele_specific_cn",
}
if not required_metrics.issubset(set(qc["metric"])):
    fail("Synthetic QC summary lacks required metrics")

print("Synthetic clonal input preparation checks passed")

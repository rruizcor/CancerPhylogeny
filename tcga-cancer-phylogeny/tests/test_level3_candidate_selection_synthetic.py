#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "09_select_level3_clonal_candidates.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("level3_candidates", SCRIPT)
level3 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["level3_candidates"] = level3
spec.loader.exec_module(level3)


def fail(message: str) -> None:
    raise AssertionError(message)


samples = pd.DataFrame(
    {
        "sample_barcode": [
            "TCGA-AA-0001-01A",
            "TCGA-AA-0002-01A",
            "TCGA-BB-0001-01A",
            "TCGA-BB-0002-01A",
            "TCGA-CC-0001-01A",
            "TCGA-DD-0001-01A",
        ],
        "patient_barcode": [
            "TCGA-AA-0001",
            "TCGA-AA-0002",
            "TCGA-BB-0001",
            "TCGA-BB-0002",
            "TCGA-CC-0001",
            "TCGA-DD-0001",
        ],
        "project_id": ["TCGA-SKCM", "TCGA-SKCM", "TCGA-UCEC", "TCGA-UCEC", "TCGA-LUAD", "TCGA-LAML"],
        "project_code": ["SKCM", "SKCM", "UCEC", "UCEC", "LUAD", "LAML"],
        "nonsynonymous_count": [500, 20, 800, 200, 300, 80],
        "total_mutation_count": [700, 30, 1000, 250, 400, 100],
    }
)
purity = pd.DataFrame(
    {
        "sample_barcode": [
            "TCGA-AA-0001-01A",
            "TCGA-AA-0002-01A",
            "TCGA-BB-0001-01A",
            "TCGA-BB-0002-01A",
            "TCGA-CC-0001-01A",
            "TCGA-DD-0001-01A",
        ],
        "purity": [0.80, 0.90, None, 0.20, 0.70, 0.60],
        "ploidy": [2.0, 2.1, 2.0, 3.0, 2.5, 2.0],
    }
)
aneuploidy = pd.DataFrame(
    {
        "sample_barcode": samples["sample_barcode"],
        "aneuploidy_score": [10, 2, 5, 8, 12, 1],
        "arm_gain_count": [5, 1, 2, 4, 6, 0],
        "arm_loss_count": [5, 1, 3, 4, 6, 1],
        "total_arm_alteration_count": [10, 2, 5, 8, 12, 1],
    }
)
projects = pd.DataFrame(
    {
        "project_id": ["TCGA-SKCM", "TCGA-UCEC", "TCGA-LUAD", "TCGA-LAML"],
        "project_code": ["SKCM", "UCEC", "LUAD", "LAML"],
        "disease_type": ["Melanoma", "Endometrial", "Lung", "Leukemia"],
        "primary_site": ["Skin", "Uterus", "Lung", "Blood"],
    }
)
group_map = pd.DataFrame(
    {
        "project_id": ["TCGA-SKCM", "TCGA-UCEC", "TCGA-LUAD", "TCGA-LAML"],
        "project_code": ["SKCM", "UCEC", "LUAD", "LAML"],
        "broad_group": ["Melanocytic", "Carcinoma", "Carcinoma", "Hematolymphoid"],
    }
)
ref_alt = pd.DataFrame(
    {
        "sample_barcode": samples["sample_barcode"],
        "n_mutation_records": [500, 20, 800, 200, 300, 80],
        "n_mutations_with_ref_alt_counts": [500, 20, 800, 200, 300, 0],
        "n_mutations_with_positive_depth": [500, 20, 800, 200, 300, 0],
        "ref_alt_counts_available": [True, True, True, True, True, False],
    }
)

config = dict(level3.CONFIG)
config["min_nonsynonymous_solid"] = 100
config["min_nonsynonymous_hematolymphoid"] = 50
config["min_purity"] = 0.30
config["pilot_project_priority"] = ["SKCM", "UCEC", "LUAD", "LAML"]
config["pilot_per_project_cap"] = 1
config["pilot_total_min"] = 2
config["pilot_total_max"] = 3

candidates = level3.build_candidate_table(
    samples,
    purity,
    aneuploidy,
    projects,
    group_map,
    ref_alt,
    segment_samples={"TCGA-CC-0001-01A"},
    arm_call_samples=set(samples["sample_barcode"]),
    config=config,
)

high = candidates[candidates["sample_barcode"] == "TCGA-AA-0001-01A"].iloc[0]
if high["candidate_class"] != "clonal_clustering_candidate_limited_cn":
    fail("High mutation/high purity sample should be a limited-CN clonal clustering candidate")
if high["candidate_score"] <= 50:
    fail("High mutation/high purity sample did not score highly")

missing_purity = candidates[candidates["sample_barcode"] == "TCGA-BB-0001-01A"].iloc[0]
if missing_purity["inclusion_status"] != "excluded" or "missing_purity" not in missing_purity["exclusion_reason"]:
    fail("Sample with missing purity was not excluded correctly")

low_mut = candidates[candidates["sample_barcode"] == "TCGA-AA-0002-01A"].iloc[0]
if "too_few_mutations" not in low_mut["exclusion_reason"]:
    fail("Low mutation-count sample was not excluded correctly")

low_purity = candidates[candidates["sample_barcode"] == "TCGA-BB-0002-01A"].iloc[0]
if "low_purity" not in low_purity["exclusion_reason"]:
    fail("Low purity sample was not excluded correctly")

full_cn = candidates[candidates["sample_barcode"] == "TCGA-CC-0001-01A"].iloc[0]
if full_cn["candidate_class"] != "copy_number_aware_clonal_phylogeny_candidate":
    fail("Segment-level CN sample was not classified as copy-number-aware")

pilot = level3.select_pilot_cohort(candidates, config)
if pilot.empty:
    fail("Synthetic pilot cohort was empty")
if pilot["project_code"].value_counts().max() > config["pilot_per_project_cap"]:
    fail("Pilot cohort ignored per-project cap")

candidates = level3.apply_pilot_status(candidates, pilot)
summary = level3.summarize_by_project(candidates)
if summary.empty or "n_samples_passing_minimum_filters" not in summary.columns:
    fail("Project-level summary was not generated")

print("Synthetic Level 3 candidate-selection checks passed")

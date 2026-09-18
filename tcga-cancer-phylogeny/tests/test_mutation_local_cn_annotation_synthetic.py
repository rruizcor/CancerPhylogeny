#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "10_annotate_mutations_with_local_cn.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("mutation_local_cn", SCRIPT)
mutation_local_cn = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["mutation_local_cn"] = mutation_local_cn
spec.loader.exec_module(mutation_local_cn)


def fail(message: str) -> None:
    raise AssertionError(message)


tmp = Path(tempfile.mkdtemp(prefix="mutation_local_cn_"))

if mutation_local_cn.normalize_chromosome("chr1") != "1":
    fail("chr-prefixed chromosome was not normalized")
if mutation_local_cn.normalize_chromosome("23") != "X":
    fail("Numeric chromosome 23 should normalize to X")
if mutation_local_cn.normalize_chromosome("chrM") != "MT":
    fail("Mitochondrial chromosome should normalize to MT")

pilot = pd.DataFrame(
    {
        "sample_barcode": ["TCGA-AA-0001-01A", "TCGA-BB-0001-01A"],
        "patient_barcode": ["TCGA-AA-0001", "TCGA-BB-0001"],
        "project_id": ["TCGA-SKCM", "TCGA-UCEC"],
        "project_code": ["SKCM", "UCEC"],
        "has_segment_level_cn": [True, True],
    }
)
purity = pd.DataFrame(
    {
        "sample_barcode": ["TCGA-AA-0001-01A", "TCGA-BB-0001-01A"],
        "purity": [0.75, 0.60],
        "ploidy": [2.1, 2.5],
    }
)
segment_summary = pd.DataFrame(
    {
        "sample_barcode": ["TCGA-AA-0001-01A", "TCGA-BB-0001-01A"],
        "has_segment_level_cn": [True, True],
    }
)
mutations = pd.DataFrame(
    {
        "sample_barcode": [
            "TCGA-AA-0001-01A",
            "TCGA-AA-0001-01A",
            "TCGA-AA-0001-01A",
            "TCGA-AA-0001-01A",
            "TCGA-BB-0001-01A",
        ],
        "patient_barcode": [
            "TCGA-AA-0001",
            "TCGA-AA-0001",
            "TCGA-AA-0001",
            "TCGA-AA-0001",
            "TCGA-BB-0001",
        ],
        "project_id": ["TCGA-SKCM", "TCGA-SKCM", "TCGA-SKCM", "TCGA-SKCM", "TCGA-UCEC"],
        "project_code": ["SKCM", "SKCM", "SKCM", "SKCM", "UCEC"],
        "Chromosome": ["chr1", "1", "chrX", "chrM", "1"],
        "Start_Position": [150, 350, 60, 10, 150],
        "End_Position": [150, 350, 60, 10, 150],
        "Hugo_Symbol": ["A", "B", "C", "MT", "D"],
        "Variant_Classification": ["Missense_Mutation", "Missense_Mutation", "Silent", "Missense_Mutation", "Missense_Mutation"],
        "Variant_Type": ["SNP", "SNP", "SNP", "SNP", "SNP"],
        "Reference_Allele": ["A", "G", "C", "T", "A"],
        "Tumor_Seq_Allele2": ["T", "A", "T", "C", "G"],
        "HGVSp_Short": ["p.A1T", "p.B1A", "p.C1C", "p.M1T", "p.D1G"],
        "t_ref_count": [30, 10, 20, None, 30],
        "t_alt_count": [10, 0, 20, 5, 15],
        "is_nonsynonymous": [True, True, False, True, True],
    }
)
mutations["chromosome"] = mutations["Chromosome"].map(mutation_local_cn.normalize_chromosome)
mutations["position"] = pd.to_numeric(mutations["Start_Position"])
mutations["start_position"] = pd.to_numeric(mutations["Start_Position"])
mutations["end_position"] = pd.to_numeric(mutations["End_Position"])
mutations["total_depth"] = pd.to_numeric(mutations["t_ref_count"]) + pd.to_numeric(mutations["t_alt_count"])
mutations["observed_vaf"] = mutations["t_alt_count"] / mutations["total_depth"]

segments_file = tmp / "segments.tsv"
pd.DataFrame(
    {
        "sample_barcode": ["TCGA-AA-0001-01A", "TCGA-AA-0001-01A", "TCGA-BB-0001-01A"],
        "patient_barcode": ["TCGA-AA-0001", "TCGA-AA-0001", "TCGA-BB-0001"],
        "project_id": ["TCGA-SKCM", "TCGA-SKCM", "TCGA-UCEC"],
        "project_code": ["SKCM", "SKCM", "UCEC"],
        "chromosome": ["1", "chrX", "1"],
        "start": [100, 50, 100],
        "end": [200, 80, 200],
        "num_probes": [10, 8, 12],
        "segment_mean": [0.30, -0.30, 0.10],
        "match_level": ["sample", "sample", "patient"],
    }
).to_csv(segments_file, sep="\t", index=False)

segments = mutation_local_cn.load_segments(segments_file, set(pilot["sample_barcode"]))
if set(segments["sample_barcode"]) != {"TCGA-AA-0001-01A"}:
    fail("Patient-level segments should not be loaded by default")

annotated = mutation_local_cn.annotate_mutations_with_segments(mutations, segments)
annotated = mutation_local_cn.attach_purity_ploidy_and_metadata(annotated, purity, pilot)
annotated = mutation_local_cn.flag_clonal_input_candidates(
    annotated,
    {**mutation_local_cn.CONFIG, "min_total_depth": 20, "candidate_nonsynonymous_only": True},
)

first = annotated[(annotated["sample_barcode"] == "TCGA-AA-0001-01A") & (annotated["position"] == 150)].iloc[0]
if first["local_cn_match_status"] != "matched":
    fail("Coordinate overlap did not annotate matching mutation")
if abs(float(first["observed_vaf"]) - 0.25) > 1e-9:
    fail("Observed VAF calculation is incorrect")
if first["local_cn_status"] != "gain_like":
    fail("Gain-like segment status was not assigned")
if not bool(first["is_clonal_input_candidate"]):
    fail("Matched nonsynonymous mutation with adequate depth should be a preview candidate")

outside = annotated[(annotated["sample_barcode"] == "TCGA-AA-0001-01A") & (annotated["position"] == 350)].iloc[0]
if outside["local_cn_match_status"] != "unmatched":
    fail("Mutation outside segment should be unmatched")
if bool(outside["is_clonal_input_candidate"]):
    fail("Unmatched mutation should not be a preview candidate")

silent = annotated[(annotated["sample_barcode"] == "TCGA-AA-0001-01A") & (annotated["chromosome"] == "X")].iloc[0]
if silent["local_cn_match_status"] != "matched" or bool(silent["is_clonal_input_candidate"]):
    fail("Silent matched mutation should be annotated but excluded from default candidate subset")

bb = annotated[annotated["sample_barcode"] == "TCGA-BB-0001-01A"].iloc[0]
if bb["local_cn_match_status"] != "unmatched":
    fail("Patient-level segment should not be used to annotate mutation by default")

summary = mutation_local_cn.build_sample_summary(
    annotated,
    pilot,
    segment_summary,
    {
        **mutation_local_cn.CONFIG,
        "min_candidate_mutations_with_local_cn": 1,
        "min_annotation_rate": 0.40,
        "min_median_depth": 20,
    },
)
aa_summary = summary[summary["sample_barcode"] == "TCGA-AA-0001-01A"].iloc[0]
if not bool(aa_summary["eligible_for_copy_number_aware_clonal_input"]):
    fail("Sample eligibility logic should pass with relaxed synthetic thresholds")
bb_summary = summary[summary["sample_barcode"] == "TCGA-BB-0001-01A"].iloc[0]
if bool(bb_summary["eligible_for_copy_number_aware_clonal_input"]):
    fail("Sample without local CN mutation annotations should not be eligible")

preview = mutation_local_cn.build_pyclone_preview(annotated)
if len(preview) != 1:
    fail("PyClone preview should contain exactly one candidate mutation")
if preview["major_cn"].notna().any() or preview["minor_cn"].notna().any():
    fail("Major/minor CN should remain NA when allele-specific CN is unavailable")
if "placeholder" not in preview["notes"].iloc[0]:
    fail("PyClone preview should label normal_cn as a placeholder")

qc = mutation_local_cn.qc_summary(annotated, summary, pilot, segments)
required_metrics = {
    "pilot_samples_evaluated",
    "total_mutations_evaluated",
    "total_mutations_annotated_with_local_cn",
    "non_allele_specific_cn_warning",
}
if not required_metrics.issubset(set(qc["metric"])):
    fail("QC summary lacks required synthetic metrics")

print("Synthetic mutation local-CN annotation checks passed")

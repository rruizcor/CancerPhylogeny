#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "12b_allele_specific_cn_readiness.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("allele_specific_cn_readiness", SCRIPT)
ascn = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["allele_specific_cn_readiness"] = ascn
spec.loader.exec_module(ascn)


def fail(message: str) -> None:
    raise AssertionError(message)


tmp = Path(tempfile.mkdtemp(prefix="allele_specific_cn_readiness_synthetic_"))
paths = ascn.default_paths(tmp)
samples = [f"TCGA-AA-000{i}-01A" for i in range(1, 5)]
patients = [sample[:12] for sample in samples]
projects = ["P1", "P1", "P2", "P2"]

pilot = pd.DataFrame(
    {
        "sample_barcode": samples,
        "patient_barcode": patients,
        "project_id": [f"TCGA-{project}" for project in projects],
        "project_code": projects,
        "has_segment_level_cn": [True] * 4,
    }
)
paths.pilot_with_segments.parent.mkdir(parents=True)
pilot.to_csv(paths.pilot_with_segments, sep="\t", index=False)

pd.DataFrame(
    {
        "sample_barcode": samples,
        "project_code": projects,
        "n_copy_neutral_candidate_mutations": [100] * 4,
        "median_depth": [40.0] * 4,
        "median_observed_vaf": [0.3] * 4,
        "purity": [0.8] * 4,
        "ploidy": [2.0] * 4,
        "selection_rank": [1, 2, 3, 4],
        "selection_reason": ["synthetic"] * 4,
    }
).to_csv(paths.limited_vaf_pilot, sep="\t", index=False)

mutation_rows = []
for sample in samples:
    for index in range(1, 101):
        position = index * 10
        mutation_rows.append(
            {
                "sample_barcode": sample,
                "chromosome": "1",
                "position": position,
                "t_ref_count": 30,
                "t_alt_count": 10,
                "total_depth": 40,
                "observed_vaf": 0.25,
                "local_segment_mean": 0.0,
                "local_cn_match_status": "matched",
                "is_clonal_input_candidate": True,
            }
        )
mutations = pd.DataFrame(mutation_rows)
paths.mutations_with_local_cn.parent.mkdir(parents=True)
mutations.to_csv(paths.mutations_with_local_cn, sep="\t", index=False, compression="gzip")

segment_mean = pd.DataFrame(
    {
        "sample_barcode": samples,
        "patient_barcode": patients,
        "project_id": [f"TCGA-{project}" for project in projects],
        "project_code": projects,
        "chromosome": ["1"] * 4,
        "start": [1] * 4,
        "end": [2000] * 4,
        "segment_mean": [0.0] * 4,
    }
)
paths.segment_mean_cn.parent.mkdir(parents=True)
segment_mean.to_csv(paths.segment_mean_cn, sep="\t", index=False, compression="gzip")
pd.DataFrame(
    {"sample_barcode": samples, "has_segment_level_cn": [True] * 4}
).to_csv(paths.segment_mean_summary, sep="\t", index=False)

paths.purity_ploidy.parent.mkdir(parents=True)
pd.DataFrame(
    {"sample_barcode": samples, "purity": [0.8] * 4, "ploidy": [2.0] * 4}
).to_csv(paths.purity_ploidy, sep="\t", index=False)
pd.DataFrame(
    {
        "sample_barcode": samples,
        "eligible_for_copy_number_aware_clustering": [False] * 4,
        "reason_not_copy_number_aware": ["requires_allele_specific_cn"] * 4,
    }
).to_csv(paths.clonal_input_readiness, sep="\t", index=False)
pd.DataFrame(
    {
        "sample_barcode": samples,
        "n_input_mutations": [100] * 4,
        "n_clusters": [3] * 4,
        "interpretation_class": ["synthetic_multicluster"] * 4,
    }
).to_csv(paths.limited_vaf_complexity, sep="\t", index=False)

facets_dir = tmp / "data" / "raw" / "facets"
sequenza_dir = tmp / "data" / "raw" / "sequenza"
ascat_dir = tmp / "data" / "raw" / "ascat"
generic_dir = tmp / "data" / "raw" / "allele_specific_cn"
for directory in [facets_dir, sequenza_dir, ascat_dir, generic_dir]:
    directory.mkdir(parents=True)

pd.DataFrame(
    {
        "sample": [samples[0], samples[0]],
        "chrom": ["chr1", "chr1"],
        "loc.start": [1, 501],
        "loc.end": [500, 2000],
        "tcn": [3, 2],
        "lcn": [1, 0],
        "cf": [0.9, 0.9],
    }
).to_csv(facets_dir / "synthetic_facets_cncf.tsv", sep="\t", index=False)

pd.DataFrame(
    {
        "sample_id": [samples[1]],
        "chromosome": ["1"],
        "start": [1],
        "end": [2000],
        "CNt": [3],
        "A": [1],
        "B": [2],
    }
).to_csv(sequenza_dir / "synthetic_sequenza_segments.tsv", sep="\t", index=False)

pd.DataFrame(
    {
        "Tumor_Sample_Barcode": [samples[2]],
        "Chromosome": ["1"],
        "Start": [1],
        "End": [2000],
        "total_copy_number": [2],
        "nMajor": [2],
        "nMinor": [0],
        "LOH": [1],
    }
).to_csv(ascat_dir / "synthetic_ascat.tsv", sep="\t", index=False)

pd.DataFrame(
    {
        "Sample": [samples[3]],
        "Chromosome": ["1"],
        "Start": [1],
        "End": [2000],
        "segment_mean": [0.0],
    }
).to_csv(generic_dir / "synthetic_segment_mean_only_segments.txt", sep="\t", index=False)

outputs = ascn.run_allele_specific_cn_readiness(paths)
segments = outputs["standardized_segments"]
availability = outputs["availability"].set_index("sample_barcode")
readiness = outputs["mutation_readiness"].set_index("sample_barcode")

if len(outputs["candidate_files"]) != 4:
    fail("Synthetic local ASCN discovery did not find all candidate files")
if len(outputs["parsed_files"]) != 3:
    fail("Synthetic segment-mean-only file was incorrectly parsed as allele-specific CN")
if set(segments["source_tool"]) != {"FACETS", "Sequenza", "ASCAT"}:
    fail("FACETS-, Sequenza-, or ASCAT-like source detection failed")

facets = segments[segments["sample_barcode"] == samples[0]]
if facets.empty or not ((facets["major_cn"] == facets["total_cn"] - facets["minor_cn"]).all()):
    fail("FACETS tcn/lcn fields were not standardized to major/minor CN")
sequenza = segments[segments["sample_barcode"] == samples[1]].iloc[0]
if float(sequenza["major_cn"]) != 2 or float(sequenza["minor_cn"]) != 1:
    fail("Sequenza A/B fields were not ordered into major/minor CN")
ascat_row = segments[segments["sample_barcode"] == samples[2]].iloc[0]
if float(ascat_row["major_cn"]) != 2 or float(ascat_row["minor_cn"]) != 0:
    fail("ASCAT major/minor fields were not parsed correctly")
if str(ascat_row["loh_status"]) != "LOH":
    fail("ASCAT LOH field was not standardized")

for sample in samples[:3]:
    if not bool(availability.loc[sample, "has_allele_specific_cn"]):
        fail(f"Synthetic ASCN availability was not detected for {sample}")
    if not bool(readiness.loc[sample, "can_build_pyclone_vi_input"]):
        fail(f"Synthetic PyClone-VI readiness was not detected for {sample}")
    if not bool(readiness.loc[sample, "can_build_phylowgs_input"]):
        fail(f"Synthetic cautious PhyloWGS readiness was not detected for {sample}")

if bool(availability.loc[samples[3], "has_allele_specific_cn"]):
    fail("Segment-mean-only synthetic sample falsely has allele-specific CN")
if bool(readiness.loc[samples[3], "can_build_pyclone_vi_input"]):
    fail("Sample without major/minor CN was falsely marked PyClone-VI ready")
if bool(readiness.loc[samples[3], "can_build_phylowgs_input"]):
    fail("Sample without major/minor CN was falsely marked PhyloWGS ready")

expected_outputs = [
    paths.standardized_segments,
    paths.availability_by_sample,
    paths.mutation_readiness,
    paths.external_manifest,
    paths.qc_summary,
    paths.planning_report,
    paths.readiness_figure,
    paths.missing_inputs_figure,
]
for path in expected_outputs:
    if not path.exists() or path.stat().st_size == 0:
        fail(f"Synthetic ASCN readiness output missing or empty: {path}")

print("Synthetic allele-specific CN readiness checks passed")

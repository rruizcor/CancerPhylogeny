#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "12_limited_vaf_clonal_clustering.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("limited_vaf", SCRIPT)
limited_vaf = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["limited_vaf"] = limited_vaf
spec.loader.exec_module(limited_vaf)


def fail(message: str) -> None:
    raise AssertionError(message)


def make_rows(sample: str, project: str, vafs: list[float], purity: float = 0.8) -> pd.DataFrame:
    rows = []
    for i, vaf in enumerate(vafs):
        depth = 100
        alt = max(1, int(round(vaf * depth)))
        rows.append(
            {
                "sample_barcode": sample,
                "patient_barcode": sample[:12],
                "project_id": f"TCGA-{project}",
                "project_code": project,
                "mutation_id": f"{sample}|1|{1000+i}|A|T|G{i}",
                "chromosome": "1",
                "position": 1000 + i,
                "Hugo_Symbol": f"G{i}",
                "Variant_Classification": "Missense_Mutation",
                "HGVSp_Short": "p.X1Y",
                "t_ref_count": depth - alt,
                "t_alt_count": alt,
                "total_depth": depth,
                "observed_vaf": alt / depth,
                "purity": purity,
                "ploidy": 2.0,
                "local_segment_mean": 0.0,
            }
        )
    return pd.DataFrame(rows)


rng = np.random.default_rng(1)
one_cluster = np.clip(rng.normal(0.38, 0.015, 120), 0.01, 0.99).tolist()
two_cluster = np.concatenate([rng.normal(0.38, 0.015, 80), rng.normal(0.15, 0.015, 55)])
two_cluster = np.clip(two_cluster, 0.01, 0.99).tolist()
small_sample = np.clip(rng.normal(0.35, 0.02, 20), 0.01, 0.99).tolist()

copy_neutral = pd.concat(
    [
        make_rows("TCGA-AA-0001-01A", "SKCM", one_cluster),
        make_rows("TCGA-BB-0001-01A", "LUAD", two_cluster),
        make_rows("TCGA-CC-0001-01A", "UCEC", small_sample),
    ],
    ignore_index=True,
)
selected = pd.DataFrame(
    {
        "sample_barcode": ["TCGA-AA-0001-01A", "TCGA-BB-0001-01A", "TCGA-CC-0001-01A"],
        "project_code": ["SKCM", "LUAD", "UCEC"],
        "purity": [0.8, 0.8, 0.8],
        "ploidy": [2.0, 2.0, 2.0],
    }
)
config = {**limited_vaf.CONFIG, "min_input_mutations": 50, "min_cluster_size": 10, "max_clusters": 4, "random_state": 1}

assignments, sample_summary = limited_vaf.run_clustering(copy_neutral, selected, config)
if assignments.empty:
    fail("Synthetic clustering produced no assignments")

aa_summary = sample_summary[sample_summary["sample_barcode"] == "TCGA-AA-0001-01A"].iloc[0]
bb_summary = sample_summary[sample_summary["sample_barcode"] == "TCGA-BB-0001-01A"].iloc[0]
cc_summary = sample_summary[sample_summary["sample_barcode"] == "TCGA-CC-0001-01A"].iloc[0]

if int(aa_summary["n_clusters"]) != 1:
    fail("GMM/BIC should recover a simple one-cluster VAF distribution")
if int(bb_summary["n_clusters"]) < 2:
    fail("GMM/BIC should recover a two-cluster VAF distribution")
if cc_summary["interpretation_class"] != "low_confidence":
    fail("Low mutation count sample should be low confidence")

bb_assign = assignments[assignments["sample_barcode"] == "TCGA-BB-0001-01A"]
if not {"clonal_like", "subclonal_like"}.issubset(set(bb_assign["cluster_interpretation"])):
    fail("Two-cluster synthetic sample should include clonal-like and subclonal-like labels")
if bb_assign.groupby("cluster_id").size().min() < config["min_cluster_size"]:
    fail("Minimum cluster size filtering failed")

project_summary = limited_vaf.build_project_summary(sample_summary)
if project_summary.empty or not {"SKCM", "LUAD", "UCEC"}.issubset(set(project_summary["project_code"])):
    fail("Project summary is missing synthetic projects")

qc = limited_vaf.build_qc(selected, assignments, sample_summary, config)
required_metrics = {
    "total_selected_pilot_samples",
    "samples_successfully_clustered",
    "samples_skipped",
    "total_clustered_mutations",
    "clustering_method_used",
    "limitations",
}
if not required_metrics.issubset(set(qc["metric"])):
    fail("Synthetic QC lacks required metrics")
limitations = qc.loc[qc["metric"] == "limitations", "value"].iloc[0]
if "not_full_copy_number_aware_phylogeny" not in limitations:
    fail("QC limitations must avoid claiming full CN-aware phylogeny")

print("Synthetic limited VAF clonal clustering checks passed")

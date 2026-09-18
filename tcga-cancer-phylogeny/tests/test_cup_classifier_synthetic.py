#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "22_cup_nearest_neighbor_classifier.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("cup_interpretation", SCRIPT)
cup = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["cup_interpretation"] = cup
spec.loader.exec_module(cup)


def fail(message: str) -> None:
    raise AssertionError(message)


tmp = Path(tempfile.mkdtemp(prefix="cup_classifier_synthetic_"))
case_id = "synthetic_cup_interpretation"
projection_dir = tmp / "results" / "clinical_projection" / case_id
query_dir = tmp / "data" / "clinical_queries" / case_id
output_dir = projection_dir / "cup_interpretation"
lineage_config = tmp / "config" / "cup_lineage_groups.yaml"
evidence_config = tmp / "config" / "cup_feature_evidence_rules.yaml"
paths = cup.default_paths(
    tmp,
    case_id,
    projection_dir=projection_dir,
    query_dir=query_dir,
    output_dir=output_dir,
    lineage_config=lineage_config,
    evidence_config=evidence_config,
)
projection_dir.mkdir(parents=True)
query_dir.mkdir(parents=True)
paths.reference_matrix.parent.mkdir(parents=True)
lineage_config.parent.mkdir(parents=True)

feature_order = [
    "driver_mutation_APC",
    "driver_mutation_KRAS",
    "epithelial_score",
    "imputed_marker",
]
reference_metadata = pd.DataFrame(
    {
        "sample_barcode": ["TCGA-AA-0001-01A", "TCGA-AA-0002-01A", "TCGA-AA-0003-01A", "TCGA-BB-0001-01A"],
        "patient_barcode": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003", "TCGA-BB-0001"],
        "project_id": ["TCGA-P1", "TCGA-P1", "TCGA-P1", "TCGA-P2"],
        "project_code": ["P1", "P1", "P1", "P2"],
        "disease_type": ["Synthetic carcinoma alpha"] * 3 + ["Synthetic carcinoma beta"],
        "primary_site": ["Synthetic site alpha"] * 3 + ["Synthetic site beta"],
        "major_cancer_group": ["Carcinoma"] * 4,
    }
)
reference_features = pd.DataFrame(
    {
        "driver_mutation_APC": [1, 1, 1, 0],
        "driver_mutation_KRAS": [1, 1, 0, 0],
        "epithelial_score": [0.8, 0.7, 0.6, -0.5],
        "imputed_marker": [0.0, 0.0, 0.0, 1.0],
    }
)
reference_matrix = pd.concat([reference_metadata, reference_features], axis=1)
reference_metadata.to_csv(paths.reference_metadata, sep="\t", index=False)
reference_matrix.to_csv(paths.reference_matrix, sep="\t", index=False)

definitions = {
    "feature_order": feature_order,
    "features": [
        {
            "name": "driver_mutation_APC",
            "feature_type": "binary",
            "feature_class": "driver_mutation_binary",
            "clinical_modality": "WES",
        },
        {
            "name": "driver_mutation_KRAS",
            "feature_type": "binary",
            "feature_class": "driver_mutation_binary",
            "clinical_modality": "WES",
        },
        {
            "name": "epithelial_score",
            "feature_type": "continuous",
            "feature_class": "expression_pathway_score",
            "clinical_modality": "WTS",
        },
        {
            "name": "imputed_marker",
            "feature_type": "continuous",
            "feature_class": "synthetic_imputed_feature",
            "clinical_modality": "WTS",
        },
    ],
}
with paths.feature_definitions.open("w") as handle:
    yaml.safe_dump(definitions, handle, sort_keys=False)

lineage_definition = {
    "schema_version": "1.0.0",
    "expected_tcga_project_count": 2,
    "lineages": [
        {
            "name": "alpha_lineage",
            "display_name": "Alpha lineage",
            "projects": ["P1"],
            "major_cancer_groups": ["Carcinoma"],
            "aliases": ["alpha"],
            "description": "Synthetic alpha group",
            "caveat": "Synthetic grouping only.",
        },
        {
            "name": "beta_lineage",
            "display_name": "Beta lineage",
            "projects": ["P2"],
            "major_cancer_groups": ["Carcinoma"],
            "aliases": ["beta"],
            "description": "Synthetic beta group",
            "caveat": "Synthetic grouping only.",
        },
    ],
}
with paths.lineage_groups.open("w") as handle:
    yaml.safe_dump(lineage_definition, handle, sort_keys=False)

evidence_definition = {
    "schema_version": "1.0.0",
    "strength_weights": {"weak": 0.33, "moderate": 0.67, "strong": 1.0},
    "rules": [
        {
            "rule_id": "synthetic_alpha_driver_pattern",
            "evidence_type": "driver mutation evidence",
            "condition": {
                "all": [
                    {
                        "feature": "driver_mutation_APC",
                        "value_source": "raw",
                        "operator": "ge",
                        "threshold": 1,
                    },
                    {
                        "feature": "driver_mutation_KRAS",
                        "value_source": "raw",
                        "operator": "ge",
                        "threshold": 1,
                    },
                ]
            },
            "finding": "Synthetic APC/KRAS pattern is present.",
            "supports": ["alpha_lineage"],
            "argues_against": [],
            "strength": "moderate",
            "notes": "Synthetic rule only.",
        },
        {
            "rule_id": "imputed_beta_marker_must_not_trigger",
            "evidence_type": "expression/pathway evidence",
            "condition": {
                "all": [
                    {
                        "feature": "imputed_marker",
                        "value_source": "raw",
                        "operator": "ge",
                        "threshold": 1,
                    }
                ]
            },
            "finding": "Synthetic imputed beta marker is present.",
            "supports": ["beta_lineage"],
            "argues_against": [],
            "strength": "strong",
            "notes": "This rule must not trigger because the feature was imputed.",
        },
    ],
    "unavailable_evidence": [
        {"label": "Synthetic unavailable evidence", "notes": "Synthetic test limitation"}
    ],
}
with paths.evidence_rules.open("w") as handle:
    yaml.safe_dump(evidence_definition, handle, sort_keys=False)

pd.DataFrame(
    {
        "case_id": [case_id],
        "sample_id": ["QUERY-1"],
        "specimen_type": ["synthetic_specimen"],
        "submitted_diagnosis": ["Synthetic carcinoma of unknown primary"],
        "differential_diagnosis": ["Synthetic alpha vs beta differential"],
        "tumor_purity": [0.7],
        "notes": ["Artificial test data only"],
    }
).to_csv(paths.clinical_metadata, sep="\t", index=False)

raw_query = pd.DataFrame(
    {
        "case_id": [case_id],
        "sample_id": ["QUERY-1"],
        "driver_mutation_APC": [1.0],
        "driver_mutation_KRAS": [1.0],
        "epithelial_score": [0.75],
        "imputed_marker": [1.0],
    }
)
scaled_query = raw_query.copy()
scaled_query["epithelial_score"] = 1.0
scaled_query["imputed_marker"] = 2.0
raw_query.to_csv(paths.clinical_feature_vector, sep="\t", index=False)
scaled_query.to_csv(paths.clinical_feature_vector_scaled, sep="\t", index=False)

pd.DataFrame(
    {
        "case_id": [case_id] * 4,
        "sample_id": ["QUERY-1"] * 4,
        "feature_name": feature_order,
        "expected_in_reference": [True] * 4,
        "supplied": [True, True, True, False],
        "imputed": [False, False, False, True],
    }
).to_csv(paths.clinical_feature_missingness, sep="\t", index=False)

distances = np.array([0.8, 0.9, 1.0, 2.0])
similarity = 1.0 / (1.0 + distances)
neighbors = pd.DataFrame(
    {
        "case_id": [case_id] * 4,
        "query_sample_id": ["QUERY-1"] * 4,
        "rank": [1, 2, 3, 4],
        "tcga_sample_barcode": reference_metadata["sample_barcode"],
        "tcga_project_id": reference_metadata["project_id"],
        "tcga_project_code": reference_metadata["project_code"],
        "tcga_disease_type": reference_metadata["disease_type"],
        "tcga_primary_site": reference_metadata["primary_site"],
        "tcga_major_cancer_group": reference_metadata["major_cancer_group"],
        "distance": distances,
        "similarity_score": similarity,
        "distance_metric": ["euclidean"] * 4,
        "notes": ["Synthetic research-only neighbor"] * 4,
    }
)
neighbors.to_csv(paths.nearest_neighbors, sep="\t", index=False)
cup.recompute_project_summary(neighbors).to_csv(paths.nearest_project_summary, sep="\t", index=False)
cup.recompute_major_group_summary(neighbors).to_csv(
    paths.nearest_major_group_summary, sep="\t", index=False
)

pd.DataFrame(
    [
        ("global", "case_id", case_id),
        ("global", "query_sample_id", "QUERY-1"),
        ("global", "number_of_expected_reference_features", 4),
        ("global", "number_of_clinical_features_supplied", 3),
        ("global", "number_of_missing_clinical_features", 1),
        ("global", "percent_feature_coverage", 75.0),
        ("global", "nearest_neighbor_metric", "euclidean"),
        ("global", "top_k_neighbors_returned", 4),
    ],
    columns=["qc_section", "metric", "value"],
).to_csv(paths.projection_qc_summary, sep="\t", index=False)

outputs = cup.interpret_cup_case(paths, case_id, {**cup.CONFIG, "top_k": 4})

ranked = outputs["ranked_lineages"]
if ranked.iloc[0]["lineage_or_project_group"] != "alpha_lineage":
    fail("Synthetic CUP interpretation did not rank the expected lineage first")
if ranked.iloc[0]["confidence_tier"] not in {
    "high_support",
    "moderate_support",
    "weak_support",
    "indeterminate",
}:
    fail("Synthetic CUP confidence tier is not populated")

ambiguity = outputs["ambiguity_metrics"].iloc[0]
if not np.isclose(float(ambiguity["wes_feature_coverage"]), 1.0):
    fail("Synthetic WES coverage is incorrect")
if not np.isclose(float(ambiguity["wts_feature_coverage"]), 0.5):
    fail("Synthetic WTS coverage is incorrect")

rule_status = {rule["rule_id"]: rule for rule in outputs["evaluated_rules"]}
if not rule_status["synthetic_alpha_driver_pattern"]["triggered"]:
    fail("Supplied synthetic driver evidence did not trigger")
if rule_status["imputed_beta_marker_must_not_trigger"]["triggered"]:
    fail("Median-imputed synthetic feature was incorrectly used as molecular evidence")
if rule_status["imputed_beta_marker_must_not_trigger"]["evaluable"]:
    fail("Median-imputed synthetic feature was incorrectly marked evaluable")

differential = outputs["differential_summary"]
if set(differential["candidate_text"]) != {"alpha", "beta"}:
    fail("Synthetic differential diagnosis candidates were not parsed")
if differential.iloc[0]["candidate_text"] != "alpha":
    fail("Synthetic differential-constrained summary did not rank alpha first")

concordance = outputs["concordance"]
if concordance.iloc[0]["concordance_class"] != "not_assessable":
    fail("Synthetic unknown-primary diagnosis should not receive an assessable concordance class")

expected_outputs = [
    paths.ranked_lineages,
    paths.ambiguity_metrics,
    paths.molecular_evidence,
    paths.differential_summary,
    paths.concordance_summary,
    paths.report,
    paths.top_project_figure,
    paths.major_group_figure,
    paths.ambiguity_figure,
    paths.evidence_heatmap,
]
for path in expected_outputs:
    if not path.exists() or path.stat().st_size == 0:
        fail(f"Synthetic CUP output missing or empty: {path}")

report = paths.report.read_text()
required_caveat = "not a clinically validated tissue-of-origin assay"
if required_caveat not in report:
    fail("Synthetic CUP report lacks the required clinical-use caveat")
if "Actionability integration is planned but not implemented" not in report:
    fail("Synthetic CUP report lacks the actionability placeholder")
if "## Molecular-Pathologic Concordance" not in report:
    fail("Synthetic CUP report lacks the molecular-pathologic concordance section")
if "This is a clinically validated" in report or "The diagnosis is" in report:
    fail("Synthetic CUP report makes an unsupported validated-diagnosis claim")

print("Synthetic CUP interpretation checks passed")

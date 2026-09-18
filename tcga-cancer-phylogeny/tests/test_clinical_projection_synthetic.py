#!/usr/bin/env python

from __future__ import annotations

import importlib.util
import json
import pickle
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "21_project_clinical_case_to_tcga.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("clinical_projection", SCRIPT)
projection = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["clinical_projection"] = projection
spec.loader.exec_module(projection)


def fail(message: str) -> None:
    raise AssertionError(message)


tmp = Path(tempfile.mkdtemp(prefix="clinical_projection_synthetic_"))
case_id = "synthetic_projection_case"
query_dir = tmp / "data" / "clinical_queries" / case_id
output_dir = tmp / "results" / "clinical_projection" / case_id
paths = projection.default_paths(tmp, case_id, query_dir=query_dir, output_dir=output_dir)
query_dir.mkdir(parents=True)
paths.reference_matrix.parent.mkdir(parents=True)

feature_order = ["mutation_count", "driver_flag", "pathway_score", "expression_pc1"]
metadata = pd.DataFrame(
    {
        "sample_barcode": [f"TCGA-AA-000{i}-01A" for i in range(1, 7)],
        "patient_barcode": [f"TCGA-AA-000{i}" for i in range(1, 7)],
        "project_id": ["TCGA-P1", "TCGA-P1", "TCGA-P1", "TCGA-P2", "TCGA-P2", "TCGA-P2"],
        "project_code": ["P1", "P1", "P1", "P2", "P2", "P2"],
        "disease_type": ["Synthetic disease 1"] * 3 + ["Synthetic disease 2"] * 3,
        "primary_site": ["Synthetic site 1"] * 3 + ["Synthetic site 2"] * 3,
        "major_cancer_group": ["Synthetic group 1"] * 3 + ["Synthetic group 2"] * 3,
    }
)
raw_features = pd.DataFrame(
    {
        "mutation_count": [0, 1, 2, 3, 4, 5],
        "driver_flag": [0, 0, 0, 1, 1, 1],
        "pathway_score": [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5],
        "expression_pc1": [2.0, 2.5, 3.0, 3.5, 4.0, 4.5],
    }
)
medians = raw_features.median()
means = raw_features.mean()
stds = raw_features.std(ddof=0)
scaled_features = (raw_features - means) / stds
raw_matrix = pd.concat([metadata, raw_features], axis=1)
scaled_matrix = pd.concat([metadata[["sample_barcode"]], scaled_features], axis=1)

contract_hash = "a" * 64
sample_hash = projection.ordered_values_hash(metadata["sample_barcode"].tolist())
feature_records = [
    {
        "name": "mutation_count",
        "feature_type": "count",
        "feature_class": "mutation_count_proxy",
        "clinical_modality": "WES",
        "description": "Synthetic mutation count",
        "retained": True,
        "included_in_nearest_neighbor_index": True,
    },
    {
        "name": "driver_flag",
        "feature_type": "binary",
        "feature_class": "driver_mutation_binary",
        "clinical_modality": "WES",
        "description": "Synthetic driver indicator",
        "retained": True,
        "included_in_nearest_neighbor_index": True,
    },
    {
        "name": "pathway_score",
        "feature_type": "continuous",
        "feature_class": "expression_pathway_score",
        "clinical_modality": "WTS",
        "description": "Synthetic pathway score",
        "retained": True,
        "included_in_nearest_neighbor_index": True,
    },
    {
        "name": "expression_pc1",
        "feature_type": "continuous",
        "feature_class": "expression_pca",
        "clinical_modality": "WTS",
        "description": "Synthetic expression PC",
        "retained": True,
        "included_in_nearest_neighbor_index": True,
    },
]
definitions = {
    "artifact": "synthetic_feature_definitions",
    "feature_contract_sha256": contract_hash,
    "metadata_columns": projection.REFERENCE_METADATA_COLUMNS,
    "raw_feature_order": feature_order,
    "feature_order": feature_order,
    "features": feature_records,
}
scaler = {
    "artifact": "synthetic_scaler",
    "feature_contract_sha256": contract_hash,
    "sample_order_sha256": sample_hash,
    "raw_feature_order": feature_order,
    "feature_order": feature_order,
    "imputation_median": {feature: float(medians[feature]) for feature in feature_order},
    "scaling_mean": {feature: float(means[feature]) for feature in feature_order},
    "scaling_std": {feature: float(stds[feature]) for feature in feature_order},
}
pca_model = PCA(n_components=2, svd_solver="full").fit(scaled_features.to_numpy())
models = {}
for metric in ["euclidean", "cosine"]:
    model = NearestNeighbors(metric=metric, algorithm="brute" if metric == "cosine" else "auto")
    model.fit(scaled_features.to_numpy())
    models[metric] = model
pca_package = {
    "model": pca_model,
    "feature_order": feature_order,
    "feature_contract_sha256": contract_hash,
    "sample_order_sha256": sample_hash,
    "sample_barcodes": metadata["sample_barcode"].tolist(),
    "explained_variance_ratio": pca_model.explained_variance_ratio_.tolist(),
}
reference_pca = pd.concat(
    [
        metadata,
        pd.DataFrame(
            pca_model.transform(scaled_features.to_numpy()),
            columns=["PC1", "PC2"],
        ),
    ],
    axis=1,
)
neighbor_package = {
    "indices": models,
    "index": models["euclidean"],
    "default_metric": "euclidean",
    "available_metrics": ["euclidean", "cosine"],
    "feature_order": feature_order,
    "feature_contract_sha256": contract_hash,
    "sample_order_sha256": sample_hash,
    "sample_barcodes": metadata["sample_barcode"].tolist(),
}
neighbor_metadata = {
    "feature_contract_sha256": contract_hash,
    "sample_order_sha256": sample_hash,
    "feature_order": feature_order,
    "default_metric": "euclidean",
    "available_metrics": ["euclidean", "cosine"],
}

raw_matrix.to_csv(paths.reference_matrix, sep="\t", index=False)
metadata.to_csv(paths.reference_metadata, sep="\t", index=False)
scaled_matrix.to_csv(paths.reference_scaled_matrix, sep="\t", index=False)
with paths.feature_definitions.open("w") as handle:
    yaml.safe_dump(definitions, handle, sort_keys=False)
with paths.scaler_parameters.open("w") as handle:
    json.dump(scaler, handle)
with paths.pca_model.open("wb") as handle:
    pickle.dump(pca_package, handle)
reference_pca.to_csv(paths.reference_pca_coordinates, sep="\t", index=False)
with paths.nearest_neighbor_index.open("wb") as handle:
    pickle.dump(neighbor_package, handle)
with paths.nearest_neighbor_metadata.open("w") as handle:
    json.dump(neighbor_metadata, handle)

pd.DataFrame(
    {
        "case_id": [case_id],
        "sample_id": ["QUERY-1"],
        "specimen_type": ["synthetic_specimen"],
        "submitted_diagnosis": ["Synthetic unknown primary"],
        "differential_diagnosis": ["Synthetic P1 versus P2"],
        "tumor_purity": [0.7],
        "notes": ["Artificial test data only"],
    }
).to_csv(paths.clinical_metadata, sep="\t", index=False)
pd.DataFrame(
    {
        "sample_id": ["QUERY-1"],
        "mutation_count": [3],
        "driver_flag": [1],
        "extra_msi_like_feature": [0.5],
    }
).to_csv(paths.wes_features, sep="\t", index=False)
pd.DataFrame({"sample_id": ["QUERY-1"], "pathway_score": [0.6]}).to_csv(
    paths.wts_features, sep="\t", index=False
)

outputs = projection.project_clinical_case(
    paths,
    case_id,
    {**projection.CONFIG, "top_k": 3, "metric": "euclidean"},
)

if outputs["neighbor_method"] != "saved_nearest_neighbor_index":
    fail("Synthetic projection did not use the compatible saved neighbor index")
if outputs["harmonized"]["supplied_count"] != 3 or outputs["harmonized"]["missing_count"] != 1:
    fail("Synthetic feature coverage counts are incorrect")
if not np.isclose(outputs["harmonized"]["feature_vector"]["expression_pc1"].iloc[0], medians["expression_pc1"]):
    fail("Missing synthetic feature was not imputed with the reference median")

scaled_query = outputs["harmonized"]["scaled_vector"][feature_order].iloc[0]
expected_scaled = (
    outputs["harmonized"]["feature_vector"][feature_order].iloc[0] - means
) / stds
if not np.allclose(scaled_query.to_numpy(dtype=float), expected_scaled.to_numpy(dtype=float)):
    fail("Synthetic scaling did not use saved reference parameters")
if len(outputs["neighbors"]) != 3 or outputs["project_summary"].empty:
    fail("Synthetic nearest neighbors or project summary are missing")
if set(outputs["neighbors"]["distance_metric"]) != {"euclidean"}:
    fail("Synthetic nearest-neighbor rows do not record the selected metric")
if not {"PC1", "PC2"}.issubset(outputs["pca_projection"].columns):
    fail("Synthetic PCA projection is missing coordinates")

direct_distances, direct_indices, direct_method, _warnings = projection.find_nearest_neighbors(
    scaled_query.to_numpy(dtype=float),
    outputs["reference"],
    top_k=3,
    metric="euclidean",
    force_direct_distance=True,
)
if direct_method != "direct_scaled_matrix_distance":
    fail("Direct-distance fallback did not activate")
if not np.allclose(direct_distances, outputs["neighbors"]["distance"].to_numpy(dtype=float)):
    fail("Direct-distance fallback and saved index return different distances")
if direct_indices.tolist() != [
    outputs["reference"]["metadata"].index[
        outputs["reference"]["metadata"]["sample_barcode"] == barcode
    ][0]
    for barcode in outputs["neighbors"]["tcga_sample_barcode"]
]:
    fail("Direct-distance fallback and saved index return different neighbor rows")

cosine_distances, cosine_indices, cosine_method, _cosine_warnings = projection.find_nearest_neighbors(
    scaled_query.to_numpy(dtype=float),
    outputs["reference"],
    top_k=3,
    metric="cosine",
    force_direct_distance=False,
)
direct_cosine_distances, direct_cosine_indices, direct_cosine_method, _direct_cosine_warnings = (
    projection.find_nearest_neighbors(
        scaled_query.to_numpy(dtype=float),
        outputs["reference"],
        top_k=3,
        metric="cosine",
        force_direct_distance=True,
    )
)
if cosine_method != "saved_nearest_neighbor_index" or direct_cosine_method != "direct_scaled_matrix_distance":
    fail("Synthetic cosine projection did not exercise saved-index and direct-distance paths")
if not np.allclose(cosine_distances, direct_cosine_distances) or not np.array_equal(
    cosine_indices, direct_cosine_indices
):
    fail("Saved-index and direct cosine results disagree")

qc_values = dict(zip(outputs["qc"]["metric"], outputs["qc"]["value"].astype(str), strict=False))
if qc_values.get("atlas_feature_count_used") != "4":
    fail("Synthetic QC does not record the locked atlas feature count")

expected_outputs = [
    paths.clinical_feature_vector,
    paths.clinical_feature_vector_scaled,
    paths.clinical_feature_missingness,
    paths.nearest_neighbors,
    paths.nearest_project_summary,
    paths.nearest_major_group_summary,
    paths.pca_projection,
    paths.qc_summary,
    paths.pca_by_project_figure,
    paths.pca_by_major_group_figure,
    paths.nearest_neighbor_figure,
    paths.feature_missingness_figure,
    paths.report,
]
for path in expected_outputs:
    if not path.exists() or path.stat().st_size == 0:
        fail(f"Synthetic clinical projection output missing or empty: {path}")

report = paths.report.read_text()
if "not a validated clinical diagnostic classifier" not in report:
    fail("Synthetic clinical report lacks the required clinical-use caveat")
if "No supervised cancer-of-unknown-primary classifier" not in report:
    fail("Synthetic clinical report does not state that no CUP classifier was run")
if "## PCA Projection Summary" not in report:
    fail("Synthetic clinical report lacks the PCA projection summary")

print("Synthetic clinical projection checks passed")

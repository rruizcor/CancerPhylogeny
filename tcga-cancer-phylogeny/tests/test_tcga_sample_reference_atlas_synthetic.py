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
import pyarrow as pa
import pyarrow.parquet as pq
import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "20_build_tcga_sample_reference_atlas.py"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("tcga_sample_reference_atlas", SCRIPT)
atlas = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["tcga_sample_reference_atlas"] = atlas
spec.loader.exec_module(atlas)


def fail(message: str) -> None:
    raise AssertionError(message)


tmp = Path(tempfile.mkdtemp(prefix="tcga_reference_atlas_"))
paths = atlas.AtlasPaths(
    mutation_parquet=tmp / "mutations.parquet",
    tmb_by_sample=tmp / "tmb_by_sample.tsv",
    purity_ploidy_by_sample=tmp / "purity_ploidy.tsv",
    aneuploidy_by_sample=tmp / "aneuploidy.tsv",
    expression_pathway_scores_by_sample=tmp / "expression.tsv",
    expression_pca_by_sample=tmp / "expression_pca.tsv",
    cancer_group_map=tmp / "cancer_group_map.csv",
    projects=tmp / "projects.tsv",
    driver_genes=tmp / "driver_genes.csv",
    reference_matrix=tmp / "results" / "reference_atlas" / "tcga_sample_reference_matrix.tsv.gz",
    metadata=tmp / "results" / "reference_atlas" / "tcga_sample_reference_metadata.tsv",
    scaled_matrix=tmp / "results" / "reference_atlas" / "tcga_sample_reference_scaled_matrix.tsv.gz",
    feature_definitions=tmp / "results" / "reference_atlas" / "tcga_reference_feature_definitions.yaml",
    scaler_parameters=tmp / "results" / "reference_atlas" / "tcga_reference_scaler_parameters.json",
    feature_missingness=tmp / "results" / "reference_atlas" / "tcga_reference_feature_missingness.tsv",
    pca_coordinates=tmp / "results" / "reference_atlas" / "tcga_reference_pca_coordinates.tsv",
    pca_model=tmp / "results" / "reference_atlas" / "tcga_reference_pca_model.pkl",
    pca_variance=tmp / "results" / "reference_atlas" / "tcga_reference_pca_variance.tsv",
    nearest_neighbor_index=tmp / "results" / "reference_atlas" / "tcga_reference_nearest_neighbor_index.pkl",
    nearest_neighbor_metadata=tmp
    / "results"
    / "reference_atlas"
    / "tcga_reference_nearest_neighbor_metadata.json",
    qc_summary=tmp / "results" / "tables" / "tcga_sample_reference_atlas_qc_summary.tsv",
    pca_by_project_figure=tmp / "results" / "figures" / "tcga_sample_reference_pca_by_project.pdf",
    pca_by_major_group_figure=tmp
    / "results"
    / "figures"
    / "tcga_sample_reference_pca_by_major_group.pdf",
    feature_missingness_figure=tmp
    / "results"
    / "figures"
    / "tcga_sample_reference_feature_missingness.pdf",
    samples_by_project_figure=tmp
    / "results"
    / "figures"
    / "tcga_sample_reference_samples_by_project.pdf",
)

samples = [
    "TCGA-AA-0001-01A",
    "TCGA-AA-0002-01A",
    "TCGA-BB-0001-01A",
    "TCGA-CC-0001-06A",
    "TCGA-ZZ-0001-11A",
]
tmb = pd.DataFrame(
    {
        "sample_barcode": samples,
        "patient_barcode": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-BB-0001", "TCGA-CC-0001", "TCGA-ZZ-0001"],
        "project_id": ["TCGA-SKCM", "TCGA-SKCM", "TCGA-UCEC", "TCGA-LUAD", "TCGA-SKCM"],
        "project_code": ["SKCM", "SKCM", "UCEC", "LUAD", "SKCM"],
        "nonsynonymous_count": [500, 25, 300, 200, 5],
        "total_mutation_count": [700, 40, 450, 260, 8],
        "metric_label": ["mutation_count_proxy_no_callable_territory"] * 5,
    }
)
purity = pd.DataFrame(
    {
        "sample_barcode": samples[:4],
        "purity": [0.80, None, 0.55, 0.60],
        "ploidy": [2.1, 2.0, None, 3.4],
    }
)
aneuploidy = pd.DataFrame(
    {
        "sample_barcode": samples[:4],
        "aneuploidy_score": [10, 2, None, 8],
        "arm_gain_count": [5, 1, 3, 4],
        "arm_loss_count": [5, 1, None, 4],
        "total_arm_alteration_count": [10, 2, 6, 8],
    }
)
expression = pd.DataFrame(
    {
        "sample_barcode": samples[:4],
        "proliferation_score": [1.0, 0.2, -0.4, 0.7],
        "immune_inflammatory_score": [0.3, None, -0.1, 0.4],
        "mostly_missing_expression_score": [1.0, None, None, None],
        "lineage_or_tissue_PC1": [2.0, 1.0, -1.0, -2.0],
        "lineage_or_tissue_PC2": [0.5, -0.5, 0.4, -0.4],
        "source": ["synthetic"] * 4,
        "notes": ["synthetic"] * 4,
    }
)
expression_pca = pd.DataFrame(
    {
        "sample_barcode": samples[:4],
        "PC1": [2.0, 1.0, -1.0, -2.0],
        "PC2": [0.5, -0.5, 0.4, -0.4],
    }
)
projects = pd.DataFrame(
    {
        "project_id": ["TCGA-SKCM", "TCGA-UCEC", "TCGA-LUAD"],
        "project_code": ["SKCM", "UCEC", "LUAD"],
        "disease_type": ["Melanoma", "Endometrial", "Lung"],
        "primary_site": ["Skin", "Uterus", "Lung"],
    }
)
group_map = pd.DataFrame(
    {
        "project_id": ["TCGA-SKCM", "TCGA-UCEC", "TCGA-LUAD"],
        "project_code": ["SKCM", "UCEC", "LUAD"],
        "broad_group": ["Melanocytic", "Carcinoma", "Carcinoma"],
    }
)
drivers = pd.DataFrame(
    {"gene": ["TP53", "KRAS", "BRAF"], "feature_group": ["tsg", "oncogene", "oncogene"]}
)
mutations = pd.DataFrame(
    {
        "sample_barcode": [
            "TCGA-AA-0001-01A",
            "TCGA-AA-0001-01A",
            "TCGA-AA-0001-01A",
            "TCGA-BB-0001-01A",
            "TCGA-CC-0001-06A",
            "TCGA-ZZ-0001-11A",
        ],
        "Hugo_Symbol": ["TP53", "TP53", "KRAS", "BRAF", "TP53", "TP53"],
        "is_nonsynonymous": [True, True, False, True, True, True],
    }
)

tmb.to_csv(paths.tmb_by_sample, sep="\t", index=False)
purity.to_csv(paths.purity_ploidy_by_sample, sep="\t", index=False)
aneuploidy.to_csv(paths.aneuploidy_by_sample, sep="\t", index=False)
expression.to_csv(paths.expression_pathway_scores_by_sample, sep="\t", index=False)
expression_pca.to_csv(paths.expression_pca_by_sample, sep="\t", index=False)
projects.to_csv(paths.projects, sep="\t", index=False)
group_map.to_csv(paths.cancer_group_map, index=False)
drivers.to_csv(paths.driver_genes, index=False)
pq.write_table(pa.Table.from_pandas(mutations), paths.mutation_parquet)

outputs = atlas.build_reference_atlas(
    paths,
    {
        **atlas.CONFIG,
        "max_feature_missing_fraction": 0.30,
        "pca_max_components": 3,
    },
)

reference = outputs["reference_matrix"]
scaled = outputs["scaled_matrix"]
metadata = outputs["metadata"]
missingness = outputs["feature_missingness"]
if len(reference) != 4:
    fail("Reference atlas did not preserve exactly the four tumor samples")
if "TCGA-ZZ-0001-11A" in set(reference["sample_barcode"]):
    fail("Normal sample barcode was not filtered out")
if "TCGA-CC-0001-06A" not in set(reference["sample_barcode"]):
    fail("Metastatic TCGA tumor sample type 06 was incorrectly removed")
if not set(atlas.METADATA_COLUMNS).issubset(metadata.columns):
    fail("Metadata table lacks required metadata columns")
if metadata[atlas.METADATA_COLUMNS].isna().any().any():
    fail("Required synthetic metadata contain missing values")

sample_one = reference["sample_barcode"] == "TCGA-AA-0001-01A"
if reference.loc[sample_one, "driver_mutation_TP53"].iloc[0] != 1:
    fail("Driver mutation binary indicator was not set")
if reference.loc[sample_one, "driver_mutation_count_TP53"].iloc[0] != 2:
    fail("Driver mutation count did not retain both nonsynonymous records")
if reference.loc[sample_one, "driver_mutation_KRAS"].iloc[0] != 0:
    fail("Synonymous driver mutation should not set the binary indicator")
if reference.loc[sample_one, "driver_mutation_count_KRAS"].iloc[0] != 0:
    fail("Synonymous driver mutation should not contribute to the count feature")

if "mostly_missing_expression_score" not in reference.columns:
    fail("Raw matrix should retain candidate features before the missingness filter")
if "mostly_missing_expression_score" in scaled.columns:
    fail("Feature with excessive missingness was not removed from the scaled matrix")
if any(column.startswith("lineage_or_tissue_PC") for column in reference.columns):
    fail("Duplicated pathway-table expression PC columns were not excluded")
if not {"expression_pca_PC1", "expression_pca_PC2"}.issubset(reference.columns):
    fail("Expression PCA features were not included from the dedicated PCA table")
if not reference["purity"].isna().any():
    fail("Raw matrix should preserve observed missing values before imputation")
if scaled.drop(columns=["sample_barcode"]).isna().any().any():
    fail("Scaled matrix contains NA values after median imputation")
if not all(
    pd.api.types.is_numeric_dtype(scaled[column])
    for column in scaled.columns
    if column != "sample_barcode"
):
    fail("Scaled matrix contains non-numeric feature columns")

dropped_row = missingness[missingness["feature_name"] == "mostly_missing_expression_score"]
if len(dropped_row) != 1 or bool(dropped_row["retained"].iloc[0]):
    fail("Feature missingness table does not record the dropped synthetic feature")
if float(dropped_row["missing_fraction"].iloc[0]) != 0.75:
    fail("Feature missingness fraction is incorrect")

for path in [
    paths.reference_matrix,
    paths.metadata,
    paths.scaled_matrix,
    paths.feature_definitions,
    paths.scaler_parameters,
    paths.feature_missingness,
    paths.pca_coordinates,
    paths.pca_model,
    paths.pca_variance,
    paths.nearest_neighbor_index,
    paths.nearest_neighbor_metadata,
    paths.qc_summary,
    paths.pca_by_project_figure,
    paths.pca_by_major_group_figure,
    paths.feature_missingness_figure,
    paths.samples_by_project_figure,
]:
    if not path.exists() or path.stat().st_size == 0:
        fail(f"Expected atlas artifact missing or empty: {path}")

with paths.scaler_parameters.open() as handle:
    scaler = json.load(handle)
for key in [
    "feature_order",
    "raw_feature_order",
    "imputation_median",
    "scaling_mean",
    "scaling_std",
    "feature_contract_sha256",
    "future_clinical_transform",
    "available_nearest_neighbor_metrics",
]:
    if key not in scaler:
        fail(f"Scaler parameter JSON lacks required field: {key}")
if scaler["available_nearest_neighbor_metrics"] != ["euclidean", "cosine"]:
    fail("Scaler metadata does not preserve both nearest-neighbor metric options")

reconstructed = reference[scaler["feature_order"]].copy()
for feature in scaler["feature_order"]:
    reconstructed[feature] = (
        reconstructed[feature].fillna(scaler["imputation_median"][feature])
        - scaler["scaling_mean"][feature]
    ) / scaler["scaling_std"][feature]
if not np.allclose(
    reconstructed.to_numpy(dtype=float),
    scaled[scaler["feature_order"]].to_numpy(dtype=float),
    atol=1e-12,
):
    fail("Saved imputation and scaling parameters do not reconstruct the scaled matrix")

with paths.feature_definitions.open() as handle:
    definitions = yaml.safe_load(handle)
if definitions.get("feature_order") != scaler["feature_order"]:
    fail("Feature definition YAML does not preserve the locked feature order")
required_definition_fields = {
    "name",
    "feature_type",
    "source_layer",
    "description",
    "clinical_modality",
    "computable_from_clinical_wes",
    "computable_from_clinical_wts",
    "included_in_nearest_neighbor_index",
}
if not definitions.get("features") or not required_definition_fields.issubset(definitions["features"][0]):
    fail("Feature definition YAML lacks required per-feature transform metadata")

pca_coordinates = pd.read_csv(paths.pca_coordinates, sep="\t")
pca_variance = pd.read_csv(paths.pca_variance, sep="\t")
if not {"PC1", "PC2"}.issubset(pca_coordinates.columns):
    fail("Synthetic PCA coordinates lack PC1/PC2")
if not {
    "component",
    "explained_variance",
    "explained_variance_ratio",
    "cumulative_explained_variance_ratio",
}.issubset(pca_variance.columns):
    fail("PCA variance table lacks required fields")

with paths.pca_model.open("rb") as handle:
    pca_model = pickle.load(handle)
if pca_model.get("feature_contract_sha256") != scaler["feature_contract_sha256"]:
    fail("PCA model does not preserve the feature contract hash")

with paths.nearest_neighbor_metadata.open() as handle:
    neighbor_metadata = json.load(handle)
if neighbor_metadata.get("feature_order") != scaler["feature_order"]:
    fail("Nearest-neighbor metadata does not preserve feature order")
if neighbor_metadata.get("default_metric") != "euclidean":
    fail("Nearest-neighbor metadata does not record the default metric")

with paths.nearest_neighbor_index.open("rb") as handle:
    neighbor_index = pickle.load(handle)
if set(neighbor_index.get("indices", {})) != {"euclidean", "cosine"}:
    fail("Nearest-neighbor artifact does not include Euclidean and cosine models")
distances, indices = neighbor_index["index"].kneighbors(
    scaled.drop(columns=["sample_barcode"]).iloc[[0]].to_numpy(dtype=float), n_neighbors=2
)
if distances.shape != (1, 2) or indices.shape != (1, 2):
    fail("Nearest-neighbor index did not return expected query shape")

qc = pd.read_csv(paths.qc_summary, sep="\t")
qc_values = dict(zip(qc["metric"], qc["value"].astype(str), strict=False))
if qc_values.get("number_of_patients_included") != "4":
    fail("QC summary does not include the synthetic patient count")

print("Synthetic TCGA sample reference atlas checks passed")

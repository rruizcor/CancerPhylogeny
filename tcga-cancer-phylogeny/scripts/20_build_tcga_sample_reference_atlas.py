#!/usr/bin/env python

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
import pickle
import platform
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcga_cancer_phylogeny_matplotlib"))
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import sklearn
import yaml
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


ATLAS_SCHEMA_VERSION = "1.0.0"

CONFIG = {
    "max_feature_missing_fraction": 0.30,
    "pca_max_components": 20,
    "nearest_neighbor_metric": "euclidean",
    "nearest_neighbor_metrics": ("euclidean", "cosine"),
}

METADATA_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "disease_type",
    "primary_site",
    "major_cancer_group",
]

MUTATION_COUNT_COLUMNS = ["nonsynonymous_count", "total_mutation_count"]
PURITY_PLOIDY_COLUMNS = ["purity", "ploidy"]
COPY_NUMBER_COLUMNS = ["aneuploidy_score", "arm_gain_count", "arm_loss_count", "total_arm_alteration_count"]
EXCLUDE_FEATURE_COLUMNS = set(METADATA_COLUMNS + ["source", "notes", "metric_label"])


@dataclass
class AtlasPaths:
    mutation_parquet: Path
    tmb_by_sample: Path
    purity_ploidy_by_sample: Path
    aneuploidy_by_sample: Path
    expression_pathway_scores_by_sample: Path
    expression_pca_by_sample: Path
    cancer_group_map: Path
    projects: Path
    driver_genes: Path
    reference_matrix: Path
    metadata: Path
    scaled_matrix: Path
    feature_definitions: Path
    scaler_parameters: Path
    feature_missingness: Path
    pca_coordinates: Path
    pca_model: Path
    pca_variance: Path
    nearest_neighbor_index: Path
    nearest_neighbor_metadata: Path
    qc_summary: Path
    pca_by_project_figure: Path
    pca_by_major_group_figure: Path
    feature_missingness_figure: Path
    samples_by_project_figure: Path


def default_paths(root: Path) -> AtlasPaths:
    return AtlasPaths(
        mutation_parquet=project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root=root),
        tmb_by_sample=project_path("data", "processed", "features", "tmb_by_sample.tsv", root=root),
        purity_ploidy_by_sample=project_path(
            "data", "processed", "features", "purity_ploidy_by_sample.tsv", root=root
        ),
        aneuploidy_by_sample=project_path("data", "processed", "features", "aneuploidy_by_sample.tsv", root=root),
        expression_pathway_scores_by_sample=project_path(
            "data", "processed", "features", "expression_pathway_scores_by_sample.tsv", root=root
        ),
        expression_pca_by_sample=project_path(
            "data", "processed", "features", "expression_pca_by_sample.tsv", root=root
        ),
        cancer_group_map=project_path("config", "cancer_group_map.csv", root=root),
        projects=project_path("data", "interim", "tcga_projects.tsv", root=root),
        driver_genes=project_path("config", "driver_genes.csv", root=root),
        reference_matrix=project_path(
            "results", "reference_atlas", "tcga_sample_reference_matrix.tsv.gz", root=root
        ),
        metadata=project_path("results", "reference_atlas", "tcga_sample_reference_metadata.tsv", root=root),
        scaled_matrix=project_path(
            "results", "reference_atlas", "tcga_sample_reference_scaled_matrix.tsv.gz", root=root
        ),
        feature_definitions=project_path(
            "results", "reference_atlas", "tcga_reference_feature_definitions.yaml", root=root
        ),
        scaler_parameters=project_path(
            "results", "reference_atlas", "tcga_reference_scaler_parameters.json", root=root
        ),
        feature_missingness=project_path(
            "results", "reference_atlas", "tcga_reference_feature_missingness.tsv", root=root
        ),
        pca_coordinates=project_path(
            "results", "reference_atlas", "tcga_reference_pca_coordinates.tsv", root=root
        ),
        pca_model=project_path("results", "reference_atlas", "tcga_reference_pca_model.pkl", root=root),
        pca_variance=project_path(
            "results", "reference_atlas", "tcga_reference_pca_variance.tsv", root=root
        ),
        nearest_neighbor_index=project_path(
            "results", "reference_atlas", "tcga_reference_nearest_neighbor_index.pkl", root=root
        ),
        nearest_neighbor_metadata=project_path(
            "results", "reference_atlas", "tcga_reference_nearest_neighbor_metadata.json", root=root
        ),
        qc_summary=project_path("results", "tables", "tcga_sample_reference_atlas_qc_summary.tsv", root=root),
        pca_by_project_figure=project_path(
            "results", "figures", "tcga_sample_reference_pca_by_project.pdf", root=root
        ),
        pca_by_major_group_figure=project_path(
            "results", "figures", "tcga_sample_reference_pca_by_major_group.pdf", root=root
        ),
        feature_missingness_figure=project_path(
            "results", "figures", "tcga_sample_reference_feature_missingness.pdf", root=root
        ),
        samples_by_project_figure=project_path(
            "results", "figures", "tcga_sample_reference_samples_by_project.pdf", root=root
        ),
    )


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_json(data: dict[str, Any], path: Path) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")
    logging.info("Wrote %s", path)


def write_yaml(data: dict[str, Any], path: Path) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=False)
    logging.info("Wrote %s", path)


def write_pickle(data: Any, path: Path) -> None:
    ensure_dir(path.parent)
    with path.open("wb") as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info("Wrote %s", path)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def normalize_sample_barcode(value: Any) -> str:
    text = "" if pd.isna(value) else str(value).strip().upper()
    parts = text.split("-")
    if len(parts) >= 4 and parts[0] == "TCGA":
        return "-".join(parts[:4])
    return text


def normalize_patient_barcode(value: Any) -> str:
    text = "" if pd.isna(value) else str(value).strip().upper()
    parts = text.split("-")
    if len(parts) >= 3:
        return "-".join(parts[:3])
    return text[:12]


def is_tcga_tumor_sample(sample_barcode: str) -> bool:
    parts = normalize_sample_barcode(sample_barcode).split("-")
    if len(parts) < 4:
        return False
    try:
        sample_type_code = int(parts[3][:2])
    except ValueError:
        return False
    return 1 <= sample_type_code <= 9


def require_columns(data: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise ValueError(f"{label} lacks required columns: {missing}")


def harmonize_sample_table(data: pd.DataFrame, label: str) -> pd.DataFrame:
    require_columns(data, ["sample_barcode"], label)
    out = data.copy()
    out["sample_barcode"] = out["sample_barcode"].map(normalize_sample_barcode)
    if (out["sample_barcode"] == "").any():
        raise ValueError(f"{label} contains blank sample barcodes")
    if "patient_barcode" in out.columns:
        out["patient_barcode"] = out["patient_barcode"].map(normalize_patient_barcode)
    for column in ["project_id", "project_code"]:
        if column in out.columns:
            out[column] = out[column].astype(str).str.strip()
    out = out.drop_duplicates()
    duplicates = out.loc[out["sample_barcode"].duplicated(keep=False), "sample_barcode"].unique().tolist()
    if duplicates:
        preview = ", ".join(duplicates[:5])
        raise ValueError(f"{label} has conflicting duplicate sample rows after barcode harmonization: {preview}")
    return out


def parse_boolean(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    normalized = series.astype(str).str.strip().str.lower()
    return normalized.isin({"true", "t", "1", "yes", "y"})


def load_driver_genes(path: Path) -> list[str]:
    drivers = pd.read_csv(path)
    if "gene" not in drivers.columns:
        raise ValueError(f"Driver gene file lacks required 'gene' column: {path}")
    genes = drivers["gene"].dropna().astype(str).str.upper().drop_duplicates().tolist()
    if not genes:
        raise ValueError("Driver gene list is empty")
    return genes


def driver_binary_feature_name(gene: str) -> str:
    return f"driver_mutation_{gene}"


def driver_count_feature_name(gene: str) -> str:
    return f"driver_mutation_count_{gene}"


def driver_mutation_features_from_mutations(mutations: pd.DataFrame, driver_genes: list[str]) -> pd.DataFrame:
    require_columns(mutations, ["sample_barcode", "Hugo_Symbol"], "mutation table")
    columns = ["sample_barcode", "Hugo_Symbol"]
    if "is_nonsynonymous" in mutations.columns:
        columns.append("is_nonsynonymous")
    data = mutations[columns].copy()
    data["sample_barcode"] = data["sample_barcode"].map(normalize_sample_barcode)
    data["Hugo_Symbol"] = data["Hugo_Symbol"].astype(str).str.upper()
    data = data[data["sample_barcode"].map(is_tcga_tumor_sample)]
    data = data[data["Hugo_Symbol"].isin(set(driver_genes))]
    if "is_nonsynonymous" in data.columns:
        data = data[parse_boolean(data["is_nonsynonymous"])]

    binary_columns = [driver_binary_feature_name(gene) for gene in driver_genes]
    count_columns = [driver_count_feature_name(gene) for gene in driver_genes]
    if data.empty:
        return pd.DataFrame(columns=["sample_barcode"] + binary_columns + count_columns)

    counts = data.groupby(["sample_barcode", "Hugo_Symbol"], sort=False).size().rename("mutation_count").reset_index()
    count_wide = (
        counts.pivot(index="sample_barcode", columns="Hugo_Symbol", values="mutation_count")
        .fillna(0)
        .astype(int)
        .reset_index()
        .rename_axis(None, axis=1)
    )
    out = count_wide[["sample_barcode"]].copy()
    for gene in driver_genes:
        gene_counts = (
            count_wide[gene]
            if gene in count_wide.columns
            else pd.Series(0, index=count_wide.index, dtype=int)
        )
        out[driver_binary_feature_name(gene)] = (gene_counts > 0).astype(int)
        out[driver_count_feature_name(gene)] = gene_counts.astype(int)
    return out[["sample_barcode"] + binary_columns + count_columns]


def load_driver_mutation_features(mutation_parquet: Path, driver_genes: list[str]) -> pd.DataFrame:
    require_file(mutation_parquet, "standardized mutation parquet")
    available_columns = set(pq.ParquetFile(mutation_parquet).schema.names)
    columns = ["sample_barcode", "Hugo_Symbol"]
    if "is_nonsynonymous" in available_columns:
        columns.append("is_nonsynonymous")
    table = pq.read_table(mutation_parquet, columns=columns)
    return driver_mutation_features_from_mutations(table.to_pandas(), driver_genes)


def validate_project_labels(tmb: pd.DataFrame, projects: pd.DataFrame, group_map: pd.DataFrame) -> None:
    require_columns(tmb, ["project_id", "project_code"], "mutation-count table")
    require_columns(projects, ["project_id", "project_code", "disease_type", "primary_site"], "TCGA project metadata")
    require_columns(group_map, ["project_id", "project_code"], "cancer group map")

    valid_pairs = set(map(tuple, projects[["project_id", "project_code"]].drop_duplicates().to_numpy()))
    observed_pairs = set(map(tuple, tmb[["project_id", "project_code"]].drop_duplicates().to_numpy()))
    invalid_pairs = sorted(observed_pairs - valid_pairs)
    if invalid_pairs:
        raise ValueError(f"Mutation-count table contains project labels absent from TCGA project metadata: {invalid_pairs}")

    group_pairs = set(map(tuple, group_map[["project_id", "project_code"]].drop_duplicates().to_numpy()))
    missing_group_pairs = sorted(observed_pairs - group_pairs)
    if missing_group_pairs:
        raise ValueError(f"Cancer group map lacks observed project labels: {missing_group_pairs}")


def metadata_from_inputs(tmb: pd.DataFrame, projects: pd.DataFrame, group_map: pd.DataFrame) -> pd.DataFrame:
    require_columns(
        tmb,
        ["sample_barcode", "patient_barcode", "project_id", "project_code"] + MUTATION_COUNT_COLUMNS,
        "mutation-count table",
    )
    validate_project_labels(tmb, projects, group_map)
    out = tmb[["sample_barcode", "patient_barcode", "project_id", "project_code"] + MUTATION_COUNT_COLUMNS].copy()
    out = out[out["sample_barcode"].map(is_tcga_tumor_sample)].copy()
    if out.empty:
        raise ValueError("No TCGA tumor samples remain after sample-type filtering")

    project_meta = projects[["project_id", "project_code", "disease_type", "primary_site"]].drop_duplicates()
    if project_meta.duplicated(["project_id", "project_code"]).any():
        raise ValueError("TCGA project metadata has duplicate project label pairs")
    out = out.merge(project_meta, on=["project_id", "project_code"], how="left", validate="many_to_one")

    group_column = "major_cancer_group" if "major_cancer_group" in group_map.columns else "broad_group"
    if group_column not in group_map.columns:
        raise ValueError("Cancer group map lacks 'major_cancer_group' or 'broad_group'")
    group_meta = group_map[["project_id", "project_code", group_column]].drop_duplicates()
    if group_meta.duplicated(["project_id", "project_code"]).any():
        raise ValueError("Cancer group map has duplicate project label pairs")
    group_meta = group_meta.rename(columns={group_column: "major_cancer_group"})
    out = out.merge(group_meta, on=["project_id", "project_code"], how="left", validate="many_to_one")

    missing_patient = out["patient_barcode"].isna() | (out["patient_barcode"].astype(str).str.strip() == "")
    out.loc[missing_patient, "patient_barcode"] = out.loc[missing_patient, "sample_barcode"].map(
        normalize_patient_barcode
    )
    out["patient_barcode"] = out["patient_barcode"].map(normalize_patient_barcode)

    required_metadata = out[METADATA_COLUMNS]
    if required_metadata.isna().any().any() or (required_metadata.astype(str).apply(lambda x: x.str.strip()) == "").any().any():
        raise ValueError("Required atlas metadata contain missing or blank values after project annotation")
    if out["sample_barcode"].duplicated().any():
        raise ValueError("Tumor sample universe contains duplicate sample barcodes")
    return out.sort_values("sample_barcode", kind="stable").reset_index(drop=True)


def selected_numeric_columns(data: pd.DataFrame, exclude: set[str] = EXCLUDE_FEATURE_COLUMNS) -> list[str]:
    columns: list[str] = []
    for column in data.columns:
        if column in exclude or column.endswith("_source") or column.endswith("_notes"):
            continue
        if column == "sample_barcode":
            continue
        if numeric(data[column]).notna().any():
            columns.append(column)
    return columns


def prepare_expression_pathway_features(expression: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    score_columns = [
        column for column in selected_numeric_columns(expression) if column.endswith("_score")
    ]
    out = expression[["sample_barcode"] + score_columns].copy()
    for column in score_columns:
        out[column] = numeric(out[column])
    return out, score_columns


def pc_sort_key(column: str) -> tuple[int, str]:
    match = re.search(r"(\d+)$", column)
    return (int(match.group(1)) if match else 10**9, column)


def prepare_expression_pca_features(
    expression_pca: pd.DataFrame | None,
    expression_pathways: pd.DataFrame,
) -> tuple[pd.DataFrame | None, list[str]]:
    if expression_pca is not None and not expression_pca.empty:
        pca_columns = sorted(
            [column for column in expression_pca.columns if re.fullmatch(r"PC\d+", column) and numeric(expression_pca[column]).notna().any()],
            key=pc_sort_key,
        )
        out = expression_pca[["sample_barcode"] + pca_columns].copy()
        rename = {column: f"expression_pca_{column}" for column in pca_columns}
        out = out.rename(columns=rename)
        output_columns = [rename[column] for column in pca_columns]
        for column in output_columns:
            out[column] = numeric(out[column])
        return out, output_columns

    fallback_columns = sorted(
        [
            column
            for column in expression_pathways.columns
            if re.fullmatch(r"lineage_or_tissue_PC\d+", column)
            and numeric(expression_pathways[column]).notna().any()
        ],
        key=pc_sort_key,
    )
    if not fallback_columns:
        return None, []
    out = expression_pathways[["sample_barcode"] + fallback_columns].copy()
    rename = {
        column: f"expression_pca_PC{re.search(r'(\d+)$', column).group(1)}" for column in fallback_columns
    }
    out = out.rename(columns=rename)
    output_columns = [rename[column] for column in fallback_columns]
    for column in output_columns:
        out[column] = numeric(out[column])
    return out, output_columns


def merge_reference_features(
    tmb: pd.DataFrame,
    projects: pd.DataFrame,
    group_map: pd.DataFrame,
    drivers: pd.DataFrame,
    purity_ploidy: pd.DataFrame,
    aneuploidy: pd.DataFrame,
    expression_pathways: pd.DataFrame,
    expression_pca: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    base = metadata_from_inputs(tmb, projects, group_map)
    base[MUTATION_COUNT_COLUMNS] = base[MUTATION_COUNT_COLUMNS].apply(numeric)

    driver_binary_columns = [column for column in drivers.columns if column.startswith("driver_mutation_") and not column.startswith("driver_mutation_count_")]
    driver_count_columns = [column for column in drivers.columns if column.startswith("driver_mutation_count_")]
    base = base.merge(drivers, on="sample_barcode", how="left", validate="one_to_one")
    for column in driver_binary_columns + driver_count_columns:
        base[column] = numeric(base[column]).fillna(0).astype(int)

    purity_columns = [column for column in PURITY_PLOIDY_COLUMNS if column in purity_ploidy.columns]
    base = base.merge(purity_ploidy[["sample_barcode"] + purity_columns], on="sample_barcode", how="left", validate="one_to_one")
    for column in purity_columns:
        base[column] = numeric(base[column])

    copy_number_columns = [column for column in COPY_NUMBER_COLUMNS if column in aneuploidy.columns]
    base = base.merge(aneuploidy[["sample_barcode"] + copy_number_columns], on="sample_barcode", how="left", validate="one_to_one")
    for column in copy_number_columns:
        base[column] = numeric(base[column])

    pathway_features, pathway_columns = prepare_expression_pathway_features(expression_pathways)
    base = base.merge(pathway_features, on="sample_barcode", how="left", validate="one_to_one")

    expression_pca_features, expression_pca_columns = prepare_expression_pca_features(
        expression_pca, expression_pathways
    )
    if expression_pca_features is not None:
        base = base.merge(expression_pca_features, on="sample_barcode", how="left", validate="one_to_one")

    feature_groups = {
        "mutation_count_proxy": MUTATION_COUNT_COLUMNS,
        "driver_mutation_binary": driver_binary_columns,
        "driver_mutation_count": driver_count_columns,
        "purity_ploidy": purity_columns,
        "copy_number_aneuploidy": copy_number_columns,
        "expression_pathway_score": pathway_columns,
        "expression_pca": expression_pca_columns,
    }
    return base, feature_groups


def all_feature_columns(feature_groups: dict[str, list[str]]) -> list[str]:
    ordered: list[str] = []
    for group_columns in feature_groups.values():
        for column in group_columns:
            if column in ordered:
                raise ValueError(f"Feature appears in more than one feature group: {column}")
            ordered.append(column)
    return ordered


def preprocess_numeric_features(
    matrix: pd.DataFrame,
    feature_columns: list[str],
    max_missing_fraction: float = CONFIG["max_feature_missing_fraction"],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    if not 0 <= max_missing_fraction <= 1:
        raise ValueError("max_feature_missing_fraction must be between 0 and 1")
    if not feature_columns:
        raise ValueError("No molecular features were assembled for the reference atlas")

    raw = matrix[feature_columns].copy()
    for column in feature_columns:
        raw[column] = numeric(raw[column])
    missing_fraction = raw.isna().mean()
    missing_count = raw.isna().sum()
    dropped = missing_fraction[missing_fraction > max_missing_fraction].index.tolist()
    retained = [column for column in feature_columns if column not in dropped]
    if not retained:
        raise ValueError("All molecular features were removed by the missingness filter")

    retained_raw = raw[retained].copy()
    medians = retained_raw.median(skipna=True)
    if medians.isna().any():
        invalid = medians[medians.isna()].index.tolist()
        raise ValueError(f"Retained features lack a finite imputation median: {invalid}")
    imputed_values = int(missing_count[retained].sum())
    imputed = retained_raw.fillna(medians)
    means = imputed.mean()
    observed_stds = imputed.std(ddof=0)
    zero_variance = observed_stds[(observed_stds == 0) | observed_stds.isna()].index.tolist()
    scaling_stds = observed_stds.replace(0, 1).fillna(1)
    scaled = (imputed - means) / scaling_stds
    if not np.isfinite(scaled.to_numpy(dtype=float)).all():
        raise ValueError("Scaled atlas matrix contains non-finite values")

    parameters: dict[str, Any] = {
        "max_feature_missing_fraction": float(max_missing_fraction),
        "raw_feature_order": feature_columns,
        "feature_order": retained,
        "dropped_features": dropped,
        "missing_count": {column: int(missing_count[column]) for column in feature_columns},
        "missing_fraction": {column: float(missing_fraction[column]) for column in feature_columns},
        "imputation_median": {column: float(medians[column]) for column in retained},
        "imputation_values": {column: float(medians[column]) for column in retained},
        "scaling_mean": {column: float(means[column]) for column in retained},
        "scaling_std": {column: float(scaling_stds[column]) for column in retained},
        "zero_variance_features": zero_variance,
        "n_imputed_values": imputed_values,
        "standardization": "z_score_after_reference_median_imputation",
        "scaling_ddof": 0,
    }
    return raw, imputed, scaled, parameters


def feature_class_for_column(column: str, feature_groups: dict[str, list[str]]) -> str:
    for feature_class, columns in feature_groups.items():
        if column in columns:
            return feature_class
    raise ValueError(f"Feature is not assigned to a feature class: {column}")


def source_for_feature_class(feature_class: str) -> str:
    return {
        "mutation_count_proxy": "data/processed/features/tmb_by_sample.tsv",
        "driver_mutation_binary": "data/processed/mutations/mc3_somatic_mutations.parquet plus config/driver_genes.csv",
        "driver_mutation_count": "data/processed/mutations/mc3_somatic_mutations.parquet plus config/driver_genes.csv",
        "purity_ploidy": "data/processed/features/purity_ploidy_by_sample.tsv",
        "copy_number_aneuploidy": "data/processed/features/aneuploidy_by_sample.tsv",
        "expression_pathway_score": "data/processed/features/expression_pathway_scores_by_sample.tsv",
        "expression_pca": "data/processed/features/expression_pca_by_sample.tsv; fallback lineage PCs from data/processed/features/expression_pathway_scores_by_sample.tsv when the dedicated PCA table is unavailable",
    }[feature_class]


def source_layer_for_feature_class(feature_class: str) -> str:
    return {
        "mutation_count_proxy": "somatic_mutation",
        "driver_mutation_binary": "somatic_mutation",
        "driver_mutation_count": "somatic_mutation",
        "purity_ploidy": "purity_ploidy",
        "copy_number_aneuploidy": "aneuploidy_arm_level_copy_number",
        "expression_pathway_score": "expression_pathway",
        "expression_pca": "expression_pca",
    }[feature_class]


def feature_type_for_column(column: str, feature_class: str) -> str:
    if feature_class == "driver_mutation_binary":
        return "binary"
    if feature_class in {"driver_mutation_count", "mutation_count_proxy"}:
        return "count"
    if column in {"arm_gain_count", "arm_loss_count", "total_arm_alteration_count"}:
        return "count"
    return "continuous"


def clinical_modality_for_feature_class(feature_class: str) -> str:
    if feature_class in {"expression_pathway_score", "expression_pca"}:
        return "WTS"
    return "WES"


def description_for_feature(column: str, feature_class: str) -> str:
    descriptions = {
        "nonsynonymous_count": "Number of nonsynonymous MC3 mutation records for the tumor sample; a mutation-count proxy without callable-territory normalization.",
        "total_mutation_count": "Total number of MC3 mutation records for the tumor sample; a mutation-count proxy without callable-territory normalization.",
        "purity": "Tumor purity estimate harmonized from the project purity/ploidy layer.",
        "ploidy": "Tumor ploidy estimate harmonized from the project purity/ploidy layer.",
        "aneuploidy_score": "Sample-level aneuploidy score from the harmonized arm-level copy-number layer.",
        "arm_gain_count": "Number of chromosome arms called gained in the sample.",
        "arm_loss_count": "Number of chromosome arms called lost in the sample.",
        "total_arm_alteration_count": "Total number of chromosome arms called gained or lost in the sample.",
    }
    if column in descriptions:
        return descriptions[column]
    if feature_class == "driver_mutation_binary":
        gene = column.removeprefix("driver_mutation_")
        return f"Binary indicator that at least one nonsynonymous MC3 mutation record affects configured driver gene {gene}."
    if feature_class == "driver_mutation_count":
        gene = column.removeprefix("driver_mutation_count_")
        return f"Count of nonsynonymous MC3 mutation records affecting configured driver gene {gene}."
    if feature_class == "expression_pathway_score":
        label = column.removesuffix("_score").replace("_", " ")
        return f"Sample-level {label} expression score computed with the project pathway-scoring workflow."
    if feature_class == "expression_pca":
        component = column.removeprefix("expression_pca_")
        return f"Sample coordinate on {component} from the project expression principal-component model."
    return f"Harmonized sample-level molecular feature: {column}."


def clinical_requirement(feature_class: str) -> str:
    if feature_class in {"mutation_count_proxy", "driver_mutation_binary", "driver_mutation_count"}:
        return "Clinical WES or an equivalent somatic mutation callset processed with matching mutation filters."
    if feature_class in {"purity_ploidy", "copy_number_aneuploidy"}:
        return "Clinical WES-derived or orthogonal purity and copy-number estimates computed with compatible conventions."
    if feature_class in {"expression_pathway_score", "expression_pca"}:
        return "Clinical WTS/RNA-seq processed with the same gene, normalization, pathway-score, and PCA conventions."
    return "Compute with the feature-specific harmonization documented for the TCGA reference atlas."


def contract_hash(parameters: dict[str, Any]) -> str:
    payload = {
        "atlas_schema_version": ATLAS_SCHEMA_VERSION,
        "feature_order": parameters["feature_order"],
        "imputation_median": parameters["imputation_median"],
        "scaling_mean": parameters["scaling_mean"],
        "scaling_std": parameters["scaling_std"],
        "standardization": parameters["standardization"],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def ordered_values_hash(values: list[str]) -> str:
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def serialization_runtime() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
    }


def finalize_scaler_parameters(
    parameters: dict[str, Any],
    sample_barcodes: list[str],
    config: dict[str, Any],
) -> dict[str, Any]:
    output = dict(parameters)
    default_metric = str(config["nearest_neighbor_metric"])
    available_metrics = list(dict.fromkeys([str(metric) for metric in config["nearest_neighbor_metrics"]]))
    if default_metric not in available_metrics:
        raise ValueError("Default nearest-neighbor metric must be listed in nearest_neighbor_metrics")
    output.update(
        {
            "artifact": "tcga_sample_reference_atlas_scaler_parameters",
            "atlas_schema_version": ATLAS_SCHEMA_VERSION,
            "generated_by": "scripts/20_build_tcga_sample_reference_atlas.py",
            "feature_order_locked": True,
            "raw_matrix_semantics": "observed_pre_imputation_values_for_all_candidate_features",
            "scaled_matrix_semantics": "retained_features_after_reference_median_imputation_and_z_score_scaling",
            "sample_order_sha256": ordered_values_hash(sample_barcodes),
            "default_nearest_neighbor_metric": default_metric,
            "available_nearest_neighbor_metrics": available_metrics,
            "serialization_runtime": serialization_runtime(),
            "future_clinical_transform": {
                "status": "reference_transform_only_no_clinical_ingestion_or_classifier",
                "required_feature_order_key": "feature_order",
                "missing_value_policy": "replace missing retained features with imputation_median",
                "scaling_formula": "(value_after_imputation - scaling_mean) / scaling_std",
                "unexpected_feature_policy": "ignore features not present in feature_order",
                "required_model_input": "finite numeric matrix in exact feature_order",
                "pca_model_file": "results/reference_atlas/tcga_reference_pca_model.pkl",
                "nearest_neighbor_index_file": "results/reference_atlas/tcga_reference_nearest_neighbor_index.pkl",
                "sample_metadata_file": "results/reference_atlas/tcga_sample_reference_metadata.tsv",
            },
        }
    )
    output["feature_contract_sha256"] = contract_hash(output)
    return output


def build_feature_definitions(
    feature_columns: list[str],
    feature_groups: dict[str, list[str]],
    scaler_parameters: dict[str, Any],
) -> dict[str, Any]:
    retained = set(scaler_parameters["feature_order"])
    definitions: list[dict[str, Any]] = []
    for column in feature_columns:
        feature_class = feature_class_for_column(column, feature_groups)
        clinical_modality = clinical_modality_for_feature_class(feature_class)
        definitions.append(
            {
                "name": column,
                "feature_type": feature_type_for_column(column, feature_class),
                "feature_class": feature_class,
                "source_layer": source_layer_for_feature_class(feature_class),
                "source_file": source_for_feature_class(feature_class),
                "description": description_for_feature(column, feature_class),
                "clinical_modality": clinical_modality,
                "computable_from_clinical_wes": clinical_modality in {"WES", "both"},
                "computable_from_clinical_wts": clinical_modality in {"WTS", "both"},
                "clinical_case_requirement": clinical_requirement(feature_class),
                "retained": column in retained,
                "included_in_nearest_neighbor_index": column in retained,
                "missing_fraction": scaler_parameters["missing_fraction"].get(column),
                "imputation_median": scaler_parameters["imputation_median"].get(column),
                "scaling_mean": scaler_parameters["scaling_mean"].get(column),
                "scaling_std": scaler_parameters["scaling_std"].get(column),
            }
        )
    return {
        "artifact": "tcga_sample_reference_atlas_feature_definitions",
        "atlas_schema_version": ATLAS_SCHEMA_VERSION,
        "generated_by": "scripts/20_build_tcga_sample_reference_atlas.py",
        "feature_contract_sha256": scaler_parameters["feature_contract_sha256"],
        "metadata_columns": METADATA_COLUMNS,
        "raw_feature_order": feature_columns,
        "feature_order": scaler_parameters["feature_order"],
        "dropped_features": scaler_parameters["dropped_features"],
        "nearest_neighbor_feature_order": scaler_parameters["feature_order"],
        "preprocessing": {
            "raw_matrix": "Observed, pre-imputation numeric values for all candidate features.",
            "missingness_filter": (
                "Drop features with missing fraction greater than "
                f"{scaler_parameters['max_feature_missing_fraction']}."
            ),
            "imputation": "Reference-median imputation using values saved in tcga_reference_scaler_parameters.json.",
            "scaling": "Population-standard-deviation z-score scaling using saved TCGA means and standard deviations.",
            "future_case_rule": "Compute retained features in exact feature_order, impute with saved medians, scale with saved means and standard deviations, then apply the saved PCA or nearest-neighbor model.",
        },
        "scope_limitations": [
            "Reference construction only; no supervised CUP classifier is trained.",
            "Mutation burden features are mutation-count proxies without callable-territory normalization.",
            "Clinical WES/WTS projection requires assay harmonization and external validation.",
        ],
        "features": definitions,
    }


def build_feature_missingness_table(
    raw: pd.DataFrame,
    feature_groups: dict[str, list[str]],
    scaler_parameters: dict[str, Any],
) -> pd.DataFrame:
    retained = set(scaler_parameters["feature_order"])
    rows = []
    for column in scaler_parameters["raw_feature_order"]:
        feature_class = feature_class_for_column(column, feature_groups)
        missing_count = int(raw[column].isna().sum())
        missing_fraction = float(raw[column].isna().mean())
        rows.append(
            {
                "feature_name": column,
                "feature_type": feature_type_for_column(column, feature_class),
                "feature_class": feature_class,
                "source_layer": source_layer_for_feature_class(feature_class),
                "n_samples": int(len(raw)),
                "n_observed": int(len(raw) - missing_count),
                "n_missing": missing_count,
                "missing_fraction": missing_fraction,
                "missingness_threshold": float(scaler_parameters["max_feature_missing_fraction"]),
                "retained": column in retained,
                "drop_reason": "" if column in retained else "missing_fraction_exceeds_threshold",
                "imputation_median": scaler_parameters["imputation_median"].get(column),
            }
        )
    return pd.DataFrame(rows)


def build_pca(
    scaled_features: pd.DataFrame,
    metadata: pd.DataFrame,
    max_components: int = CONFIG["pca_max_components"],
) -> tuple[PCA, pd.DataFrame, pd.DataFrame]:
    if len(scaled_features) < 2:
        raise ValueError("At least two tumor samples are required to build PCA")
    n_components = min(int(max_components), scaled_features.shape[0] - 1, scaled_features.shape[1])
    if n_components < 1:
        raise ValueError("At least one retained feature is required to build PCA")
    values = scaled_features.to_numpy(dtype=np.float64)
    pca = PCA(n_components=n_components, svd_solver="full")
    coordinates = pca.fit_transform(values)
    coord_columns = [f"PC{i + 1}" for i in range(n_components)]
    coords = pd.concat(
        [metadata[METADATA_COLUMNS].reset_index(drop=True), pd.DataFrame(coordinates, columns=coord_columns)],
        axis=1,
    )
    variance = pd.DataFrame(
        {
            "component": coord_columns,
            "explained_variance": pca.explained_variance_,
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative_explained_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
            "singular_value": pca.singular_values_,
        }
    )
    return pca, coords, variance


def build_nearest_neighbor_indices(
    scaled_features: pd.DataFrame,
    metrics: list[str],
) -> dict[str, NearestNeighbors]:
    values = np.ascontiguousarray(scaled_features.to_numpy(dtype=np.float64))
    indices: dict[str, NearestNeighbors] = {}
    for metric in metrics:
        algorithm = "brute" if metric == "cosine" else "auto"
        model = NearestNeighbors(metric=metric, algorithm=algorithm, n_jobs=1)
        model.fit(values)
        indices[metric] = model
    return indices


def color_mapping(values: pd.Series) -> dict[str, Any]:
    labels = sorted(values.fillna("Unannotated").astype(str).unique())
    cmap = plt.get_cmap("nipy_spectral", max(len(labels), 1))
    return {label: cmap(i) for i, label in enumerate(labels)}


def plot_pca(
    coords: pd.DataFrame,
    pca_variance: pd.DataFrame,
    color_column: str,
    path: Path,
    title: str,
) -> None:
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(13, 8))
    if {"PC1", "PC2"}.issubset(coords.columns):
        plot_data = coords.copy()
        plot_data[color_column] = plot_data[color_column].fillna("Unannotated").astype(str)
        colors = color_mapping(plot_data[color_column])
        for label, group in plot_data.groupby(color_column, sort=True):
            ax.scatter(
                group["PC1"],
                group["PC2"],
                s=9,
                alpha=0.58,
                color=colors[label],
                label=label,
                linewidths=0,
                rasterized=True,
            )
        variance = dict(zip(pca_variance["component"], pca_variance["explained_variance_ratio"], strict=False))
        ax.set_xlabel(f"PC1 ({100 * variance.get('PC1', 0):.1f}% variance)")
        ax.set_ylabel(f"PC2 ({100 * variance.get('PC2', 0):.1f}% variance)")
        ax.set_title(title)
        legend_columns = 2 if plot_data[color_column].nunique() > 18 else 1
        ax.legend(
            markerscale=2,
            fontsize=7,
            ncol=legend_columns,
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
        )
        ax.grid(color="#D9D9D9", linewidth=0.4, alpha=0.5)
    else:
        ax.text(0.5, 0.5, "PCA has fewer than two components", ha="center", va="center")
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_feature_missingness(missingness: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    plotted = missingness[missingness["missing_fraction"] > 0].copy()
    if plotted.empty:
        plotted = missingness.nlargest(min(20, len(missingness)), "missing_fraction").copy()
    plotted = plotted.sort_values(["missing_fraction", "feature_name"], ascending=[True, True])
    height = max(6.0, min(13.0, 0.28 * len(plotted)))
    fig, ax = plt.subplots(figsize=(11, height))
    colors = np.where(plotted["retained"], "#3B7A57", "#B84A4A")
    ax.barh(plotted["feature_name"], plotted["missing_fraction"], color=colors)
    threshold = float(missingness["missingness_threshold"].iloc[0])
    ax.axvline(threshold, color="#222222", linestyle="--", linewidth=1.0, label=f"Drop threshold ({threshold:.0%})")
    ax.set_xlabel("Missing fraction")
    ax.set_ylabel("")
    ax.set_xlim(0, max(threshold * 1.12, plotted["missing_fraction"].max() * 1.12, 0.05))
    zero_count = int((missingness["missing_fraction"] == 0).sum())
    ax.set_title(f"TCGA reference feature missingness ({zero_count} zero-missingness features not shown)")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_samples_by_project(metadata: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    counts = metadata["project_code"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(13, 6.5))
    bars = ax.bar(counts.index, counts.values, color="#356A8A")
    ax.bar_label(bars, labels=[str(value) for value in counts.values], fontsize=6, padding=2, rotation=90)
    ax.set_xlabel("TCGA project")
    ax.set_ylabel("Tumor samples")
    ax.set_title("TCGA reference atlas samples by project")
    ax.tick_params(axis="x", labelrotation=90, labelsize=8)
    ax.margins(x=0.01)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def build_nearest_neighbor_metadata(
    metadata: pd.DataFrame,
    scaler_parameters: dict[str, Any],
    metrics: list[str],
    default_metric: str,
) -> dict[str, Any]:
    return {
        "artifact": "tcga_sample_reference_atlas_nearest_neighbor_metadata",
        "atlas_schema_version": ATLAS_SCHEMA_VERSION,
        "generated_by": "scripts/20_build_tcga_sample_reference_atlas.py",
        "feature_contract_sha256": scaler_parameters["feature_contract_sha256"],
        "sample_order_sha256": scaler_parameters["sample_order_sha256"],
        "default_metric": default_metric,
        "available_metrics": metrics,
        "serialization_runtime": serialization_runtime(),
        "n_samples": int(len(metadata)),
        "n_features": int(len(scaler_parameters["feature_order"])),
        "feature_order": scaler_parameters["feature_order"],
        "model_file": "results/reference_atlas/tcga_reference_nearest_neighbor_index.pkl",
        "sample_metadata_file": "results/reference_atlas/tcga_sample_reference_metadata.tsv",
        "scaler_parameters_file": "results/reference_atlas/tcga_reference_scaler_parameters.json",
        "query_contract": {
            "input": "One or more finite numeric rows after saved-median imputation and saved-parameter scaling.",
            "column_order": "Use feature_order exactly.",
            "neighbor_row_mapping": "Returned integer indices map to row order in tcga_sample_reference_metadata.tsv and sample_barcodes in the pickle artifact.",
            "distance_interpretation": "Smaller values indicate greater molecular similarity within the selected metric; distances are descriptive and not diagnostic probabilities.",
        },
        "scope": "Unsupervised reference search only; no clinical case ingestion and no supervised CUP classifier.",
    }


def qc_summary(
    matrix: pd.DataFrame,
    feature_columns: list[str],
    scaler_parameters: dict[str, Any],
    expression_included: bool,
    pca_built: bool,
    nn_built: bool,
) -> pd.DataFrame:
    project_counts = matrix["project_code"].value_counts().sort_index()
    rows: list[tuple[str, str, Any]] = [
        ("global", "atlas_schema_version", ATLAS_SCHEMA_VERSION),
        ("global", "number_of_samples_included", len(matrix)),
        ("global", "number_of_patients_included", matrix["patient_barcode"].nunique()),
        ("global", "number_of_projects_represented", matrix["project_code"].nunique()),
        (
            "global",
            "samples_per_project",
            ";".join(f"{project_code}={int(count)}" for project_code, count in project_counts.items()),
        ),
        ("global", "number_of_raw_features", len(feature_columns)),
        ("global", "number_of_retained_features", len(scaler_parameters["feature_order"])),
        ("global", "number_of_dropped_features", len(scaler_parameters["dropped_features"])),
        (
            "global",
            "dropped_features",
            ";".join(scaler_parameters["dropped_features"]) if scaler_parameters["dropped_features"] else "none",
        ),
        ("global", "number_of_imputed_values", scaler_parameters["n_imputed_values"]),
        ("global", "expression_included", expression_included),
        ("global", "pca_built", pca_built),
        ("global", "nearest_neighbor_index_built", nn_built),
        (
            "global",
            "nearest_neighbor_metrics",
            ";".join(scaler_parameters["available_nearest_neighbor_metrics"]),
        ),
        (
            "global",
            "default_nearest_neighbor_metric",
            scaler_parameters["default_nearest_neighbor_metric"],
        ),
        ("global", "feature_contract_sha256", scaler_parameters["feature_contract_sha256"]),
        (
            "global",
            "limitations",
            "Reference construction only; no supervised CUP classifier or clinical validation. TCGA is primary-tumor enriched; mutation burden is a mutation-count proxy; median imputation and cross-assay WES/WTS harmonization require external validation before clinical projection.",
        ),
    ]
    for project_code, count in project_counts.items():
        rows.append((f"project:{project_code}", "samples_included", int(count)))
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def validate_atlas_outputs_in_memory(
    raw_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    scaled_matrix: pd.DataFrame,
    scaler_parameters: dict[str, Any],
) -> None:
    if not (len(raw_matrix) == len(metadata) == len(scaled_matrix)):
        raise ValueError("Reference, metadata, and scaled matrices have inconsistent row counts")
    if metadata["sample_barcode"].duplicated().any():
        raise ValueError("Reference metadata contain duplicate sample barcodes")
    if not metadata["sample_barcode"].map(is_tcga_tumor_sample).all():
        raise ValueError("Reference metadata contain a non-tumor TCGA sample barcode")
    if metadata[METADATA_COLUMNS].isna().any().any():
        raise ValueError("Reference metadata contain missing required values")
    expected_scaled_columns = ["sample_barcode"] + scaler_parameters["feature_order"]
    if list(scaled_matrix.columns) != expected_scaled_columns:
        raise ValueError("Scaled matrix columns do not match the locked feature order")
    if not raw_matrix["sample_barcode"].equals(metadata["sample_barcode"]):
        raise ValueError("Raw matrix and metadata sample ordering differ")
    if not scaled_matrix["sample_barcode"].equals(metadata["sample_barcode"]):
        raise ValueError("Scaled matrix and metadata sample ordering differ")


def build_reference_atlas(paths: AtlasPaths, config: dict[str, Any] = CONFIG) -> dict[str, pd.DataFrame]:
    required = [
        (paths.mutation_parquet, "standardized mutation parquet"),
        (paths.tmb_by_sample, "sample mutation-count proxy table"),
        (paths.purity_ploidy_by_sample, "sample purity/ploidy table"),
        (paths.aneuploidy_by_sample, "sample aneuploidy table"),
        (paths.expression_pathway_scores_by_sample, "sample expression pathway-score table"),
        (paths.cancer_group_map, "cancer group map"),
        (paths.projects, "TCGA project metadata"),
        (paths.driver_genes, "driver gene config"),
    ]
    for path, description in required:
        require_file(path, description)

    tmb = harmonize_sample_table(read_tsv(paths.tmb_by_sample), "mutation-count table")
    purity_ploidy = harmonize_sample_table(read_tsv(paths.purity_ploidy_by_sample), "purity/ploidy table")
    aneuploidy = harmonize_sample_table(read_tsv(paths.aneuploidy_by_sample), "aneuploidy table")
    expression = harmonize_sample_table(
        read_tsv(paths.expression_pathway_scores_by_sample), "expression pathway-score table"
    )
    expression_pca = (
        harmonize_sample_table(read_tsv(paths.expression_pca_by_sample), "expression PCA table")
        if paths.expression_pca_by_sample.exists()
        else None
    )
    projects = read_tsv(paths.projects)
    group_map = pd.read_csv(paths.cancer_group_map)
    for table in [projects, group_map]:
        for column in ["project_id", "project_code"]:
            if column in table.columns:
                table[column] = table[column].astype(str).str.strip()
    driver_genes = load_driver_genes(paths.driver_genes)
    drivers = load_driver_mutation_features(paths.mutation_parquet, driver_genes)

    matrix, feature_groups = merge_reference_features(
        tmb,
        projects,
        group_map,
        drivers,
        purity_ploidy,
        aneuploidy,
        expression,
        expression_pca,
    )
    feature_columns = all_feature_columns(feature_groups)
    raw, _imputed, scaled, scaler_parameters = preprocess_numeric_features(
        matrix, feature_columns, float(config["max_feature_missing_fraction"])
    )
    metadata = matrix[METADATA_COLUMNS].copy().reset_index(drop=True)
    sample_barcodes = metadata["sample_barcode"].astype(str).tolist()
    scaler_parameters = finalize_scaler_parameters(scaler_parameters, sample_barcodes, config)
    retained_features = scaler_parameters["feature_order"]

    raw_matrix = pd.concat(
        [metadata.reset_index(drop=True), raw[feature_columns].reset_index(drop=True)],
        axis=1,
    )
    scaled_matrix = pd.concat(
        [metadata[["sample_barcode"]].reset_index(drop=True), scaled[retained_features].reset_index(drop=True)],
        axis=1,
    )
    pca, pca_coordinates, pca_variance = build_pca(
        scaled[retained_features], metadata, int(config["pca_max_components"])
    )
    metrics = scaler_parameters["available_nearest_neighbor_metrics"]
    default_metric = scaler_parameters["default_nearest_neighbor_metric"]
    nn_indices = build_nearest_neighbor_indices(scaled[retained_features], metrics)

    feature_definitions = build_feature_definitions(feature_columns, feature_groups, scaler_parameters)
    feature_missingness = build_feature_missingness_table(raw, feature_groups, scaler_parameters)
    nearest_neighbor_metadata = build_nearest_neighbor_metadata(
        metadata, scaler_parameters, metrics, default_metric
    )
    expression_included = bool(feature_groups["expression_pathway_score"] or feature_groups["expression_pca"])
    qc = qc_summary(matrix, feature_columns, scaler_parameters, expression_included, True, True)
    validate_atlas_outputs_in_memory(raw_matrix, metadata, scaled_matrix, scaler_parameters)

    write_tsv_gz(raw_matrix, paths.reference_matrix)
    write_tsv(metadata, paths.metadata)
    write_tsv_gz(scaled_matrix, paths.scaled_matrix)
    write_yaml(feature_definitions, paths.feature_definitions)
    write_json(scaler_parameters, paths.scaler_parameters)
    write_tsv(feature_missingness, paths.feature_missingness)
    write_tsv(pca_coordinates, paths.pca_coordinates)
    write_pickle(
        {
            "artifact": "tcga_sample_reference_atlas_pca_model",
            "atlas_schema_version": ATLAS_SCHEMA_VERSION,
            "feature_contract_sha256": scaler_parameters["feature_contract_sha256"],
            "sample_order_sha256": scaler_parameters["sample_order_sha256"],
            "serialization_runtime": serialization_runtime(),
            "model": pca,
            "feature_order": retained_features,
            "sample_barcodes": sample_barcodes,
            "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        },
        paths.pca_model,
    )
    write_tsv(pca_variance, paths.pca_variance)
    write_pickle(
        {
            "artifact": "tcga_sample_reference_atlas_nearest_neighbor_index",
            "atlas_schema_version": ATLAS_SCHEMA_VERSION,
            "feature_contract_sha256": scaler_parameters["feature_contract_sha256"],
            "sample_order_sha256": scaler_parameters["sample_order_sha256"],
            "serialization_runtime": serialization_runtime(),
            "index": nn_indices[default_metric],
            "indices": nn_indices,
            "metric": default_metric,
            "default_metric": default_metric,
            "available_metrics": metrics,
            "feature_order": retained_features,
            "sample_barcodes": sample_barcodes,
        },
        paths.nearest_neighbor_index,
    )
    write_json(nearest_neighbor_metadata, paths.nearest_neighbor_metadata)
    write_tsv(qc, paths.qc_summary)
    plot_pca(
        pca_coordinates,
        pca_variance,
        "project_code",
        paths.pca_by_project_figure,
        "TCGA sample reference atlas PCA by project",
    )
    plot_pca(
        pca_coordinates,
        pca_variance,
        "major_cancer_group",
        paths.pca_by_major_group_figure,
        "TCGA sample reference atlas PCA by major cancer group",
    )
    plot_feature_missingness(feature_missingness, paths.feature_missingness_figure)
    plot_samples_by_project(metadata, paths.samples_by_project_figure)
    return {
        "reference_matrix": raw_matrix,
        "metadata": metadata,
        "scaled_matrix": scaled_matrix,
        "feature_missingness": feature_missingness,
        "pca_coordinates": pca_coordinates,
        "pca_variance": pca_variance,
        "qc": qc,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a locked sample-level TCGA molecular reference atlas.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument(
        "--max-feature-missing-fraction",
        type=float,
        default=CONFIG["max_feature_missing_fraction"],
        help="Drop numeric features with missing fraction greater than this threshold.",
    )
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    config = dict(CONFIG)
    config["max_feature_missing_fraction"] = args.max_feature_missing_fraction
    paths = default_paths(root)
    outputs = build_reference_atlas(paths, config)
    logging.info(
        "TCGA sample reference atlas complete: %d samples, %d projects, %d retained features",
        len(outputs["metadata"]),
        outputs["metadata"]["project_code"].nunique(),
        len(outputs["scaled_matrix"].columns) - 1,
    )


if __name__ == "__main__":
    main()

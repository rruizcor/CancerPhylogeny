#!/usr/bin/env python

from __future__ import annotations

import argparse
import html
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from lib.common import configure_logging, ensure_dir, find_project_root, project_path


TITLE = "TCGA Multilevel Molecular Similarity, Clonal-Structure, and Clinical Projection Framework"


ANNOTATION_COLUMNS = {
    "project_id",
    "project_code",
    "cancer_name",
    "broad_group",
    "level2_group",
    "lineage_notes",
}


@dataclass(frozen=True)
class ReportPaths:
    main_report: Path
    executive_summary: Path
    methods: Path
    limitations_next_steps: Path
    html_report: Path
    project_summary: Path
    outputs_index: Path


def default_report_paths(root: Path) -> ReportPaths:
    return ReportPaths(
        main_report=project_path("results", "reports", "integrated_tcga_cancer_phylogeny_report.md", root=root),
        executive_summary=project_path(
            "results", "reports", "integrated_tcga_cancer_phylogeny_executive_summary.md", root=root
        ),
        methods=project_path("results", "reports", "integrated_tcga_cancer_phylogeny_methods.md", root=root),
        limitations_next_steps=project_path(
            "results", "reports", "integrated_tcga_cancer_phylogeny_limitations_next_steps.md", root=root
        ),
        html_report=project_path("results", "reports", "integrated_tcga_cancer_phylogeny_report.html", root=root),
        project_summary=project_path("results", "tables", "integrated_project_summary.tsv", root=root),
        outputs_index=project_path("results", "tables", "integrated_project_outputs_index.tsv", root=root),
    )


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        logging.warning("Missing optional report input: %s", path)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def clean_value(value: object, default: str = "not_available") -> str:
    if value is None:
        return default
    try:
        if pd.isna(value):
            return default
    except TypeError:
        pass
    text = str(value)
    if text.lower() == "nan":
        return default
    return text


def metric_value(df: pd.DataFrame, metric: str, default: str = "not_available") -> str:
    if df.empty or "metric" not in df.columns or "value" not in df.columns:
        return default
    rows = df[df["metric"].astype(str) == metric]
    if rows.empty:
        return default
    return clean_value(rows.iloc[0]["value"], default=default)


def section_metric_value(
    df: pd.DataFrame,
    section: str,
    metric: str,
    default: str = "not_available",
) -> str:
    if df.empty or not {"qc_section", "metric", "value"}.issubset(df.columns):
        return default
    rows = df[
        (df["qc_section"].astype(str) == section)
        & (df["metric"].astype(str) == metric)
    ]
    if rows.empty:
        return default
    return clean_value(rows.iloc[0]["value"], default=default)


def first_value(df: pd.DataFrame, column: str, default: str = "not_available") -> str:
    if df.empty or column not in df.columns:
        return default
    return clean_value(df.iloc[0][column], default=default)


def format_float(value: object, digits: int = 4, default: str = "not_available") -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return default


def percentage(numerator: object, denominator: object, default: str = "not_available") -> str:
    try:
        denominator_value = float(denominator)
        if denominator_value == 0:
            return default
        return f"{100.0 * float(numerator) / denominator_value:.3f}"
    except (TypeError, ValueError):
        return default


def boolish(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def rel(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")
    logging.info("Wrote %s", path)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def output_record(root: Path, level: str, output_type: str, file_path: str, description: str, notes: str = "") -> dict:
    path = root / file_path
    return {
        "level": level,
        "output_type": output_type,
        "file_path": file_path,
        "description": description,
        "exists": bool(path.exists()),
        "notes": notes,
    }


def figure_records(root: Path, level: str, pattern: str, description_prefix: str) -> list[dict]:
    records = []
    for path in sorted((root / "results" / "figures").glob(pattern)):
        records.append(
            output_record(
                root,
                level,
                "figure",
                rel(path, root),
                f"{description_prefix}: {path.name}",
                "auto-indexed existing figure",
            )
        )
    return records


def directory_records(
    root: Path,
    level: str,
    directory: str,
    pattern: str,
    output_type: str,
    description_prefix: str,
) -> list[dict]:
    records = []
    for path in sorted((root / directory).glob(pattern)):
        if not path.is_file():
            continue
        records.append(
            output_record(
                root,
                level,
                output_type,
                rel(path, root),
                f"{description_prefix}: {path.name}",
                "auto-indexed existing output",
            )
        )
    return records


def build_outputs_index(root: Path, paths: ReportPaths) -> pd.DataFrame:
    records = [
        output_record(root, "integrated", "report", rel(paths.main_report, root), "Main integrated project report"),
        output_record(root, "integrated", "report", rel(paths.executive_summary, root), "Executive summary for collaborators"),
        output_record(root, "integrated", "report", rel(paths.methods, root), "Concise methods summary"),
        output_record(
            root,
            "integrated",
            "report",
            rel(paths.limitations_next_steps, root),
            "Limitations and next recommended development steps",
        ),
        output_record(root, "integrated", "report", rel(paths.html_report, root), "Optional HTML rendering of main report"),
        output_record(
            root,
            "integrated",
            "table",
            rel(paths.project_summary, root),
            "Machine-readable integrated project summary",
        ),
        output_record(root, "integrated", "table", rel(paths.outputs_index, root), "Machine-readable project output index"),
        output_record(root, "level1", "table", "results/tables/level1_feature_matrix.tsv", "Level 1 project feature matrix"),
        output_record(root, "level1", "table", "results/tables/level1_tree_qc_summary.tsv", "Level 1 tree QC summary"),
        output_record(
            root,
            "level1",
            "table",
            "results/tables/level1_bootstrap_cluster_stability.tsv",
            "Level 1 feature-bootstrap co-clustering stability",
        ),
        output_record(
            root,
            "level1",
            "tree",
            "results/trees/level1_pan_cancer_gower_hclust_tree.nwk",
            "Level 1 Gower average-linkage Newick tree",
        ),
        output_record(
            root,
            "level1",
            "tree",
            "results/trees/level1_pan_cancer_correlation_hclust_tree.nwk",
            "Level 1 correlation average-linkage Newick tree",
        ),
        output_record(root, "level2", "table", "results/tables/level2_group_definitions.tsv", "Level 2 group definitions"),
        output_record(root, "level2", "table", "results/tables/level2_tree_qc_summary.tsv", "Level 2 tree QC summary"),
        output_record(root, "level2", "report", "results/reports/level2_interpretation_summary.md", "Level 2 interpretation summary"),
        output_record(
            root,
            "level3",
            "table",
            "results/tables/level3_candidate_selection_qc_summary.tsv",
            "Level 3 candidate-selection QC summary",
        ),
        output_record(
            root,
            "level3",
            "table",
            "results/tables/level3_candidate_summary_by_project.tsv",
            "Level 3 candidate summary by project",
        ),
        output_record(root, "level3", "table", "results/tables/level3_pilot_cohort.tsv", "Level 3 pilot cohort"),
        output_record(
            root,
            "level3",
            "table",
            "results/tables/level3_segment_cn_qc_summary.tsv",
            "Level 3 segment-level copy-number QC summary",
        ),
        output_record(
            root,
            "level3",
            "table",
            "results/tables/level3_mutation_local_cn_annotation_qc_summary.tsv",
            "Level 3 mutation local-CN annotation QC summary",
        ),
        output_record(
            root,
            "level3",
            "table",
            "results/tables/level3_clonal_input_preparation_qc_summary.tsv",
            "Level 3 clonal-input preparation QC summary",
        ),
        output_record(
            root,
            "level3",
            "table",
            "results/tables/level3_limited_vaf_clustering_qc_summary.tsv",
            "Level 3 limited VAF clustering QC summary",
        ),
        output_record(
            root,
            "level3",
            "table",
            "results/tables/level3_limited_vaf_clonal_complexity_by_sample.tsv",
            "Level 3 limited VAF complexity by sample",
        ),
        output_record(
            root,
            "level3",
            "table",
            "results/tables/level3_limited_vaf_clonal_complexity_by_project.tsv",
            "Level 3 limited VAF complexity by project",
        ),
        output_record(
            root,
            "translational",
            "table",
            "results/tables/tcga_sample_reference_atlas_qc_summary.tsv",
            "TCGA sample-level reference atlas QC summary",
        ),
        output_record(
            root,
            "translational",
            "artifact",
            "results/reference_atlas/tcga_reference_feature_definitions.yaml",
            "Locked reference-atlas feature definitions used by clinical projection",
        ),
        output_record(
            root,
            "translational",
            "artifact",
            "results/reference_atlas/tcga_reference_scaler_parameters.json",
            "Locked reference-atlas imputation and scaling parameters",
        ),
        output_record(
            root,
            "level3_readiness",
            "table",
            "results/tables/level3_allele_specific_cn_readiness_qc_summary.tsv",
            "Allele-specific CN readiness QC summary",
        ),
        output_record(
            root,
            "level3_readiness",
            "table",
            "results/tables/level3_allele_specific_cn_availability_by_sample.tsv",
            "Allele-specific CN availability by pilot sample",
        ),
        output_record(
            root,
            "level3_readiness",
            "table",
            "results/tables/level3_mutation_to_allele_specific_cn_readiness.tsv",
            "Mutation overlap readiness for allele-specific CN",
        ),
        output_record(
            root,
            "level3_readiness",
            "table",
            "results/tables/level3_allele_specific_cn_external_inference_manifest.tsv",
            "External allele-specific CN inference manifest",
        ),
        output_record(
            root,
            "level3_readiness",
            "report",
            "results/reports/level3_allele_specific_cn_and_pyclone_plan.md",
            "Allele-specific CN acquisition and PyClone-VI planning report",
        ),
        output_record(
            root,
            "translational",
            "matrix",
            "results/reference_atlas/tcga_sample_reference_matrix.tsv.gz",
            "Observed pre-imputation TCGA sample reference matrix",
        ),
        output_record(
            root,
            "translational",
            "table",
            "results/reference_atlas/tcga_sample_reference_metadata.tsv",
            "TCGA sample reference metadata",
        ),
        output_record(
            root,
            "translational",
            "matrix",
            "results/reference_atlas/tcga_sample_reference_scaled_matrix.tsv.gz",
            "Median-imputed and scaled TCGA sample reference matrix",
        ),
        output_record(
            root,
            "translational",
            "table",
            "results/reference_atlas/tcga_reference_feature_missingness.tsv",
            "Reference-atlas feature missingness",
        ),
        output_record(
            root,
            "translational",
            "table",
            "results/reference_atlas/tcga_reference_pca_coordinates.tsv",
            "Reference-atlas PCA coordinates",
        ),
        output_record(
            root,
            "translational",
            "model",
            "results/reference_atlas/tcga_reference_pca_model.pkl",
            "Locked reference PCA model",
        ),
        output_record(
            root,
            "translational",
            "table",
            "results/reference_atlas/tcga_reference_pca_variance.tsv",
            "Reference PCA explained variance",
        ),
        output_record(
            root,
            "translational",
            "model",
            "results/reference_atlas/tcga_reference_nearest_neighbor_index.pkl",
            "Default nearest-neighbor index",
        ),
        output_record(
            root,
            "translational",
            "artifact",
            "results/reference_atlas/tcga_reference_nearest_neighbor_metadata.json",
            "Nearest-neighbor index metadata and feature contract",
        ),
    ]
    records += figure_records(root, "level1", "level1*.pdf", "Level 1 figure")
    records += figure_records(root, "level2", "level2*.pdf", "Level 2 figure")
    records += figure_records(root, "level3", "level3*.pdf", "Level 3 figure")
    records += figure_records(root, "translational", "tcga_sample_reference*.pdf", "Reference atlas figure")
    records += directory_records(
        root,
        "clinical_projection",
        "results/clinical_projection/example_cup_case",
        "*",
        "clinical_projection_output",
        "Synthetic clinical projection output",
    )
    records += directory_records(
        root,
        "cup_interpretation",
        "results/clinical_projection/example_cup_case/cup_interpretation",
        "*",
        "cup_interpretation_output",
        "Synthetic CUP interpretation output",
    )
    return pd.DataFrame(records, columns=["level", "output_type", "file_path", "description", "exists", "notes"])


def list_join(values: Iterable[str], empty: str = "none") -> str:
    filtered = [str(value) for value in values if str(value) and str(value).lower() != "nan"]
    return ", ".join(filtered) if filtered else empty


def build_context(root: Path) -> dict:
    tables = {
        "level1_matrix": read_tsv(root / "results" / "tables" / "level1_feature_matrix.tsv"),
        "level1_qc": read_tsv(root / "results" / "tables" / "level1_tree_qc_summary.tsv"),
        "level1_bootstrap": read_tsv(root / "results" / "tables" / "level1_bootstrap_cluster_stability.tsv"),
        "level2_definitions": read_tsv(root / "results" / "tables" / "level2_group_definitions.tsv"),
        "level2_qc": read_tsv(root / "results" / "tables" / "level2_tree_qc_summary.tsv"),
        "level3_candidate_qc": read_tsv(root / "results" / "tables" / "level3_candidate_selection_qc_summary.tsv"),
        "level3_candidate_project": read_tsv(root / "results" / "tables" / "level3_candidate_summary_by_project.tsv"),
        "level3_pilot": read_tsv(root / "results" / "tables" / "level3_pilot_cohort.tsv"),
        "level3_segment_qc": read_tsv(root / "results" / "tables" / "level3_segment_cn_qc_summary.tsv"),
        "level3_annotation_qc": read_tsv(
            root / "results" / "tables" / "level3_mutation_local_cn_annotation_qc_summary.tsv"
        ),
        "level3_input_qc": read_tsv(root / "results" / "tables" / "level3_clonal_input_preparation_qc_summary.tsv"),
        "level3_vaf_qc": read_tsv(root / "results" / "tables" / "level3_limited_vaf_clustering_qc_summary.tsv"),
        "level3_vaf_sample": read_tsv(
            root / "results" / "tables" / "level3_limited_vaf_clonal_complexity_by_sample.tsv"
        ),
        "level3_vaf_project": read_tsv(
            root / "results" / "tables" / "level3_limited_vaf_clonal_complexity_by_project.tsv"
        ),
        "level3_ascn_qc": read_tsv(
            root / "results" / "tables" / "level3_allele_specific_cn_readiness_qc_summary.tsv"
        ),
        "reference_atlas_qc": read_tsv(root / "results" / "tables" / "tcga_sample_reference_atlas_qc_summary.tsv"),
        "clinical_qc": read_tsv(
            root
            / "results"
            / "clinical_projection"
            / "example_cup_case"
            / "clinical_projection_qc_summary.tsv"
        ),
        "clinical_project": read_tsv(
            root
            / "results"
            / "clinical_projection"
            / "example_cup_case"
            / "tcga_nearest_project_summary.tsv"
        ),
        "clinical_major_group": read_tsv(
            root
            / "results"
            / "clinical_projection"
            / "example_cup_case"
            / "tcga_nearest_major_group_summary.tsv"
        ),
        "cup_ranked": read_tsv(
            root
            / "results"
            / "clinical_projection"
            / "example_cup_case"
            / "cup_interpretation"
            / "cup_ranked_lineage_interpretation.tsv"
        ),
        "cup_ambiguity": read_tsv(
            root
            / "results"
            / "clinical_projection"
            / "example_cup_case"
            / "cup_interpretation"
            / "cup_ambiguity_metrics.tsv"
        ),
        "cup_evidence": read_tsv(
            root
            / "results"
            / "clinical_projection"
            / "example_cup_case"
            / "cup_interpretation"
            / "cup_molecular_evidence_table.tsv"
        ),
        "cup_concordance": read_tsv(
            root
            / "results"
            / "clinical_projection"
            / "example_cup_case"
            / "cup_interpretation"
            / "cup_molecular_pathologic_concordance.tsv"
        ),
    }
    level1_matrix = tables["level1_matrix"]
    level2_defs = tables["level2_definitions"]
    vaf_project = tables["level3_vaf_project"]
    vaf_sample = tables["level3_vaf_sample"]

    context = {
        "tables": tables,
        "level1_projects": metric_value(tables["level1_qc"], "included_tcga_projects", str(len(level1_matrix))),
        "level1_features": metric_value(
            tables["level1_qc"],
            "features_retained_for_tree",
            str(len([col for col in level1_matrix.columns if col not in ANNOTATION_COLUMNS])),
        ),
        "level1_distance_methods": metric_value(tables["level1_qc"], "distance_methods_generated", "Gower;correlation"),
        "level1_bootstrap_iterations": metric_value(tables["level1_qc"], "bootstrap_iterations_completed", "not_recorded"),
        "level2_tree_groups": [],
        "level2_underpowered_groups": [],
        "level3_top_candidate_projects": [],
        "level3_vaf_project_lines": [],
        "clinical_project_lines": [],
        "cup_evidence_types": [],
    }

    if not level2_defs.empty:
        attempted = level2_defs[level2_defs["tree_attempted"].map(boolish)]
        underpowered = level2_defs[~level2_defs["tree_attempted"].map(boolish)]
        context["level2_tree_groups"] = attempted["group_name"].astype(str).tolist()
        context["level2_underpowered_groups"] = [
            f"{row.group_name} ({clean_value(row.reason_if_not_attempted, default='not_attempted')})"
            for row in underpowered.itertuples(index=False)
        ]

    candidate_project = tables["level3_candidate_project"]
    if not candidate_project.empty and "n_samples_passing_minimum_filters" in candidate_project.columns:
        top = candidate_project.sort_values("n_samples_passing_minimum_filters", ascending=False).head(8)
        context["level3_top_candidate_projects"] = [
            f"{row.project_code}: {int(row.n_samples_passing_minimum_filters)}" for row in top.itertuples(index=False)
        ]

    if not vaf_project.empty:
        context["level3_vaf_project_lines"] = [
            f"{row.project_code}: n={int(row.n_samples)}, median clusters={float(row.median_n_clusters):.1f}"
            for row in vaf_project.itertuples(index=False)
        ]

    if not vaf_sample.empty and "interpretation_class" in vaf_sample.columns:
        context["vaf_class_counts"] = vaf_sample["interpretation_class"].value_counts().to_dict()
    else:
        context["vaf_class_counts"] = {}

    clinical_project = tables["clinical_project"]
    if not clinical_project.empty:
        clinical_project = clinical_project.sort_values("rank")
        context["clinical_project_lines"] = [
            (
                f"{row.project_code}: {int(row.n_neighbors)} neighbors, "
                f"{float(row.weighted_similarity_score):.4f} similarity-weight share"
            )
            for row in clinical_project.itertuples(index=False)
        ]

    cup_evidence = tables["cup_evidence"]
    if not cup_evidence.empty and "evidence_type" in cup_evidence.columns:
        context["cup_evidence_types"] = cup_evidence["evidence_type"].astype(str).tolist()
    return context


def add_summary_row(rows: list[dict], level: str, metric: str, value: object, source_file: str, notes: str = "") -> None:
    rows.append(
        {
            "level": level,
            "metric": metric,
            "value": clean_value(value),
            "source_file": source_file,
            "notes": notes,
        }
    )


def build_project_summary(context: dict) -> pd.DataFrame:
    tables = context["tables"]
    rows: list[dict] = []
    add_summary_row(rows, "overall", "tcga_projects_represented", context["level1_projects"], "results/tables/level1_tree_qc_summary.tsv")
    add_summary_row(rows, "level1", "biologic_features_used", context["level1_features"], "results/tables/level1_tree_qc_summary.tsv")
    add_summary_row(
        rows,
        "level1",
        "distance_methods_generated",
        context["level1_distance_methods"],
        "results/tables/level1_tree_qc_summary.tsv",
    )
    add_summary_row(
        rows,
        "level1",
        "bootstrap_iterations_completed",
        context["level1_bootstrap_iterations"],
        "results/tables/level1_tree_qc_summary.tsv",
    )
    add_summary_row(
        rows,
        "level2",
        "groups_with_project_level_trees",
        list_join(context["level2_tree_groups"]),
        "results/tables/level2_group_definitions.tsv",
    )
    add_summary_row(
        rows,
        "level2",
        "underpowered_groups",
        list_join(context["level2_underpowered_groups"]),
        "results/tables/level2_group_definitions.tsv",
    )
    candidate_qc = tables["level3_candidate_qc"]
    for metric in [
        "total_mutation_layer_samples",
        "samples_with_ref_alt_counts",
        "samples_with_purity",
        "samples_with_ploidy",
        "samples_with_copy_number_or_aneuploidy",
        "samples_passing_minimum_filters",
        "samples_selected_for_pilot_cohort",
    ]:
        add_summary_row(rows, "level3_candidate_selection", metric, metric_value(candidate_qc, metric), "results/tables/level3_candidate_selection_qc_summary.tsv")
    segment_qc = tables["level3_segment_qc"]
    for metric in [
        "gdc_segment_files_found",
        "gdc_segment_files_matched_to_pilot",
        "selected_best_files",
        "total_segment_rows_processed",
        "pilot_samples_with_segment_level_cn",
    ]:
        add_summary_row(rows, "level3_segment_cn", metric, metric_value(segment_qc, metric), "results/tables/level3_segment_cn_qc_summary.tsv")
    annotation_qc = tables["level3_annotation_qc"]
    for metric in [
        "total_mutations_evaluated",
        "total_mutations_annotated_with_local_cn",
        "percentage_mutations_annotated",
        "samples_eligible_for_downstream_clonal_input",
    ]:
        add_summary_row(
            rows,
            "level3_mutation_local_cn",
            metric,
            metric_value(annotation_qc, metric),
            "results/tables/level3_mutation_local_cn_annotation_qc_summary.tsv",
        )
    input_qc = tables["level3_input_qc"]
    for metric in [
        "total_copy_neutral_candidate_mutations",
        "samples_eligible_for_limited_vaf_clustering",
        "samples_eligible_for_copy_number_aware_clustering",
        "samples_requiring_allele_specific_cn",
    ]:
        add_summary_row(rows, "level3_clonal_input", metric, metric_value(input_qc, metric), "results/tables/level3_clonal_input_preparation_qc_summary.tsv")
    vaf_qc = tables["level3_vaf_qc"]
    for metric in [
        "total_selected_pilot_samples",
        "samples_successfully_clustered",
        "samples_skipped",
        "total_clustered_mutations",
        "predominantly_clonal_like",
        "oligoclonal_like",
        "multicluster_subclonal_like",
    ]:
        add_summary_row(rows, "level3_limited_vaf", metric, metric_value(vaf_qc, metric), "results/tables/level3_limited_vaf_clustering_qc_summary.tsv")
    atlas_qc = tables["reference_atlas_qc"]
    for metric in [
        "number_of_samples_included",
        "number_of_patients_included",
        "number_of_projects_represented",
        "number_of_raw_features",
        "number_of_retained_features",
        "number_of_dropped_features",
        "dropped_features",
        "number_of_imputed_values",
        "nearest_neighbor_metrics",
        "default_nearest_neighbor_metric",
        "nearest_neighbor_index_built",
        "pca_built",
    ]:
        add_summary_row(rows, "translational_reference_atlas", metric, metric_value(atlas_qc, metric), "results/tables/tcga_sample_reference_atlas_qc_summary.tsv")

    clinical_qc = tables["clinical_qc"]
    for metric in [
        "number_of_expected_reference_features",
        "number_of_clinical_features_supplied",
        "number_of_missing_clinical_features",
        "percent_feature_coverage",
        "number_of_extra_clinical_features_ignored",
        "nearest_neighbor_metric",
        "nearest_neighbor_method",
        "top_k_neighbors_returned",
        "PCA_model_reused_without_refit",
        "projection_status",
    ]:
        add_summary_row(
            rows,
            "clinical_projection_example",
            metric,
            section_metric_value(clinical_qc, "global", metric),
            "results/clinical_projection/example_cup_case/clinical_projection_qc_summary.tsv",
        )
    for modality in ["WES", "WTS"]:
        for metric in ["number_of_expected_reference_features", "number_of_clinical_features_supplied"]:
            add_summary_row(
                rows,
                "clinical_projection_example",
                f"{modality.lower()}_{metric}",
                section_metric_value(clinical_qc, modality, metric),
                "results/clinical_projection/example_cup_case/clinical_projection_qc_summary.tsv",
            )

    cup_ranked = tables["cup_ranked"]
    cup_ambiguity = tables["cup_ambiguity"]
    cup_concordance = tables["cup_concordance"]
    for metric, value, source_file in [
        (
            "top_lineage",
            first_value(cup_ranked, "lineage_or_project_group"),
            "results/clinical_projection/example_cup_case/cup_interpretation/cup_ranked_lineage_interpretation.tsv",
        ),
        (
            "top_lineage_score",
            first_value(cup_ranked, "similarity_score"),
            "results/clinical_projection/example_cup_case/cup_interpretation/cup_ranked_lineage_interpretation.tsv",
        ),
        (
            "top_lineage_confidence_tier",
            first_value(cup_ranked, "confidence_tier"),
            "results/clinical_projection/example_cup_case/cup_interpretation/cup_ranked_lineage_interpretation.tsv",
        ),
        (
            "ambiguity_class",
            first_value(cup_ambiguity, "ambiguity_class"),
            "results/clinical_projection/example_cup_case/cup_interpretation/cup_ambiguity_metrics.tsv",
        ),
        (
            "top1_top2_margin",
            first_value(cup_ambiguity, "top1_top2_margin"),
            "results/clinical_projection/example_cup_case/cup_interpretation/cup_ambiguity_metrics.tsv",
        ),
        (
            "entropy_like_score",
            first_value(cup_ambiguity, "entropy_like_score"),
            "results/clinical_projection/example_cup_case/cup_interpretation/cup_ambiguity_metrics.tsv",
        ),
        (
            "molecular_pathologic_concordance",
            first_value(cup_concordance, "concordance_class"),
            "results/clinical_projection/example_cup_case/cup_interpretation/cup_molecular_pathologic_concordance.tsv",
        ),
        (
            "evidence_categories_generated",
            len(context["cup_evidence_types"]),
            "results/clinical_projection/example_cup_case/cup_interpretation/cup_molecular_evidence_table.tsv",
        ),
    ]:
        add_summary_row(
            rows,
            "cup_interpretation_example",
            metric,
            value,
            source_file,
        )

    ascn_qc = tables["level3_ascn_qc"]
    for metric in [
        "candidate_local_allele_specific_cn_files_found",
        "samples_with_major_cn",
        "samples_with_minor_cn",
        "samples_ready_for_pyclone_vi",
        "samples_ready_for_phylowgs",
        "samples_requiring_external_allele_specific_cn_inference",
    ]:
        add_summary_row(
            rows,
            "level3_allele_specific_cn_readiness",
            metric,
            metric_value(ascn_qc, metric),
            "results/tables/level3_allele_specific_cn_readiness_qc_summary.tsv",
        )
    return pd.DataFrame(rows, columns=["level", "metric", "value", "source_file", "notes"])


def bullet_list(items: Iterable[str]) -> str:
    return "\n".join(f"- {item}" for item in items)


def build_main_report(root: Path, context: dict) -> str:
    tables = context["tables"]
    candidate_qc = tables["level3_candidate_qc"]
    segment_qc = tables["level3_segment_qc"]
    annotation_qc = tables["level3_annotation_qc"]
    input_qc = tables["level3_input_qc"]
    vaf_qc = tables["level3_vaf_qc"]
    atlas_qc = tables["reference_atlas_qc"]
    clinical_qc = tables["clinical_qc"]
    clinical_major_group = tables["clinical_major_group"]
    cup_ranked = tables["cup_ranked"]
    cup_ambiguity = tables["cup_ambiguity"]
    cup_concordance = tables["cup_concordance"]
    ascn_qc = tables["level3_ascn_qc"]

    level2_groups = list_join(context["level2_tree_groups"])
    underpowered = list_join(context["level2_underpowered_groups"])
    top_projects = list_join(context["level3_top_candidate_projects"])
    vaf_project_lines = bullet_list(context["level3_vaf_project_lines"])
    clinical_project_lines = bullet_list(context["clinical_project_lines"])

    wes_expected = section_metric_value(clinical_qc, "WES", "number_of_expected_reference_features")
    wes_supplied = section_metric_value(clinical_qc, "WES", "number_of_clinical_features_supplied")
    wts_expected = section_metric_value(clinical_qc, "WTS", "number_of_expected_reference_features")
    wts_supplied = section_metric_value(clinical_qc, "WTS", "number_of_clinical_features_supplied")
    top_lineage = first_value(cup_ranked, "lineage_or_project_group")
    top_lineage_score = format_float(first_value(cup_ranked, "similarity_score"))
    top_confidence = first_value(cup_ranked, "confidence_tier")
    if len(cup_ranked) > 1:
        second_lineage = clean_value(cup_ranked.iloc[1]["lineage_or_project_group"])
        second_lineage_score = format_float(cup_ranked.iloc[1]["similarity_score"])
    else:
        second_lineage = "not_available"
        second_lineage_score = "not_available"

    return f"""# {TITLE}

Generated by `scripts/13_generate_integrated_project_report.py`.

## Project Overview

This project implements a multilevel TCGA molecular analysis and translational projection framework.

- Level 1 builds pan-cancer molecular similarity trees across TCGA cancer types using project-level molecular traits.
- Level 2 builds within-group molecular similarity trees for biologically motivated cancer groups.
- Level 3 builds a limited within-tumor clonal-structure prototype using copy-neutral mutation VAF clustering in selected high-information samples.
- A locked 97-feature TCGA sample-level atlas preserves tumor-level heterogeneity for research projection.
- Clinical WES/WTS and CUP prototypes harmonize an already-computed query feature vector to that atlas and produce uncertainty-aware, non-diagnostic molecular comparisons.

Levels 1 and 2 are molecular similarity, taxonomic, and trait-mapping trees. They are not literal species-like evolutionary phylogenies and should not be used to infer directional ancestry between cancer types. Level 3 is closer to tumor evolutionary inference because it works within individual tumor samples, but the current implementation is limited VAF-based clustering. It does not run PyClone-VI or PhyloWGS, does not infer allele-specific integer copy number, and does not infer definitive phylogenetic branching.

The clinical projection and CUP outputs are research prototypes. Their scores are descriptive similarity weights, not tissue-of-origin probabilities or a clinical diagnosis.

## Data Resources and Current Scope

The project integrates MC3 somatic mutation calls, tumor ref/alt read counts, purity and ploidy estimates, aneuploidy and arm-level copy-number calls, GDC segment-level copy-number files for the Level 3 pilot cohort, and RNA expression/pathway features.

Key current data-resource counts:

- TCGA projects represented in Level 1: {context["level1_projects"]} / 33.
- Mutation-layer tumor samples evaluated for Level 3: {metric_value(candidate_qc, "total_mutation_layer_samples")}.
- Samples with ref/alt counts: {metric_value(candidate_qc, "samples_with_ref_alt_counts")}.
- Samples with purity and ploidy: {metric_value(candidate_qc, "samples_with_purity")} with purity and {metric_value(candidate_qc, "samples_with_ploidy")} with ploidy.
- Samples with copy-number or aneuploidy data: {metric_value(candidate_qc, "samples_with_copy_number_or_aneuploidy")}.
- Segment-level CN files found for pilot projects: {metric_value(segment_qc, "gdc_segment_files_found")}; pilot-matched rows: {metric_value(segment_qc, "gdc_segment_files_matched_to_pilot")}; selected best files: {metric_value(segment_qc, "selected_best_files")}.
- Reference-atlas samples: {metric_value(atlas_qc, "number_of_samples_included")} tumors from {metric_value(atlas_qc, "number_of_patients_included")} patients across {metric_value(atlas_qc, "number_of_projects_represented")} projects.
- Reference-atlas features: {metric_value(atlas_qc, "number_of_raw_features")} raw and {metric_value(atlas_qc, "number_of_retained_features")} retained; {metric_value(atlas_qc, "number_of_dropped_features")} dropped; {metric_value(atlas_qc, "number_of_imputed_values")} median-imputed modeling values.

## Level 1 Results: Pan-Cancer Molecular Similarity

Level 1 includes {context["level1_projects"]} TCGA projects and {context["level1_features"]} retained biologic features spanning mutation-count proxies, driver mutation prevalence, purity, ploidy, aneuploidy/copy-number burden, and expression pathway scores. Generated distance and tree views include {context["level1_distance_methods"]}. Bootstrap iterations completed: {context["level1_bootstrap_iterations"]}.

The Gower tree emphasizes scaled absolute molecular trait differences. In the current interpretation, it separates a broad lower-mutation or lineage-distinct cluster containing kidney, CNS, endocrine/adrenal, liver/biliary, mesothelioma, prostate, thyroid, thymoma, uveal melanoma, and sarcoma projects; a broad carcinoma-enriched cluster including BLCA, BRCA, CESC, COAD/READ, ESCA/HNSC, LUAD/LUSC, OV, and STAD; and singleton or outlier-like DLBC, SKCM, UCEC, and UCS groups.

The correlation tree emphasizes relative feature-profile shape rather than absolute scaled trait levels. It therefore differs from Gower by grouping projects with similar molecular profiles even when their overall mutation burden, copy-number burden, or expression-program magnitude differs.

Key Level 1 outputs:

- `results/tables/level1_feature_matrix.tsv`
- `results/tables/level1_tree_qc_summary.tsv`
- `results/tables/level1_bootstrap_cluster_stability.tsv`
- `results/trees/level1_pan_cancer_gower_hclust_tree.nwk`
- `results/trees/level1_pan_cancer_correlation_hclust_tree.nwk`
- `results/figures/level1_pan_cancer_gower_tree.pdf`
- `results/figures/level1_pan_cancer_correlation_tree.pdf`
- `results/figures/level1_pan_cancer_tree_with_trait_heatmap.pdf`
- `results/figures/level1_pan_cancer_feature_pca.pdf`

## Level 2 Results: Within-Group Molecular Similarity

Level 2 reuses the project-level molecular trait space but filters and scales features within biologically motivated groups. Groups with generated trees are: {level2_groups}. Underpowered groups are: {underpowered}.

The carcinoma-all Gower tree groups a broad carcinoma cluster containing BLCA, BRCA, CESC, ESCA, HNSC, LUAD, LUSC, OV, and STAD; a CHOL/KICH/KIRP/LIHC/PRAD/THCA cluster; a COAD/READ pair; and separate KIRC, PAAD, and UCEC branches in the current k-level summary. Expected subgroup signals appear cautiously: HNSC/ESCA cluster within the squamous-enriched set; COAD/READ/STAD and CHOL/LIHC cluster within the GI/hepatobiliary set; KIRC/KIRP cluster separately from KICH in kidney; and BRCA/OV cluster within the gynecologic/breast set, with UCEC, CESC, and UCS separate.

Level 2 remains descriptive. Small group sizes, especially kidney with three projects and underpowered two-project groups, limit robust topology interpretation.

Key Level 2 outputs:

- `results/tables/level2_group_definitions.tsv`
- `results/tables/level2_tree_qc_summary.tsv`
- `results/reports/level2_interpretation_summary.md`
- `results/figures/level2_carcinoma_all_gower_tree.pdf`
- `results/figures/level2_pan_squamous_gower_tree.pdf`
- `results/figures/level2_gi_pancancreatobiliary_gower_tree.pdf`
- `results/figures/level2_kidney_gower_tree.pdf`
- `results/figures/level2_gynecologic_breast_gower_tree.pdf`

## Level 3 Results: Limited Within-Tumor Clonal Structure Prototype

Candidate selection evaluated {metric_value(candidate_qc, "total_mutation_layer_samples")} mutation-layer tumor samples. {metric_value(candidate_qc, "samples_passing_minimum_filters")} samples passed minimum filters and {metric_value(candidate_qc, "samples_selected_for_pilot_cohort")} samples were selected for the pilot cohort. Top candidate projects by passing samples include: {top_projects}.

The segment-CN acquisition step found {metric_value(segment_qc, "gdc_segment_files_found")} GDC segment-CN files across pilot projects, matched {metric_value(segment_qc, "gdc_segment_files_matched_to_pilot")} rows to pilot samples or patients, selected {metric_value(segment_qc, "selected_best_files")} best files, parsed {metric_value(segment_qc, "gdc_unique_files_parsed_successfully")} files, and processed {metric_value(segment_qc, "total_segment_rows_processed")} segment rows. {metric_value(segment_qc, "pilot_samples_with_segment_level_cn")} pilot samples have segment-level CN in the current processed outputs.

Mutation-to-segment annotation evaluated {metric_value(annotation_qc, "total_mutations_evaluated")} pilot mutations and annotated {metric_value(annotation_qc, "total_mutations_annotated_with_local_cn")} with local segment CN ({metric_value(annotation_qc, "percentage_mutations_annotated")}%). {metric_value(annotation_qc, "samples_eligible_for_downstream_clonal_input")} pilot samples were eligible for downstream clonal-input preparation.

Clonal-input stratification produced {metric_value(input_qc, "total_copy_neutral_candidate_mutations")} copy-neutral candidate mutations and {metric_value(input_qc, "total_segment_annotated_candidate_mutations")} broader segment-annotated candidate mutations. {metric_value(input_qc, "samples_eligible_for_limited_vaf_clustering")} samples were eligible for limited VAF clustering, while {metric_value(input_qc, "samples_eligible_for_copy_number_aware_clustering")} were eligible for fully copy-number-aware clustering because allele-specific major/minor integer CN is absent.

The limited VAF clustering prototype clustered {metric_value(vaf_qc, "total_selected_pilot_samples")} selected pilot samples and {metric_value(vaf_qc, "total_clustered_mutations")} copy-neutral mutations, with {metric_value(vaf_qc, "samples_skipped")} skipped samples. Current interpretation classes are {metric_value(vaf_qc, "predominantly_clonal_like")} predominantly clonal-like, {metric_value(vaf_qc, "oligoclonal_like")} oligoclonal-like, and {metric_value(vaf_qc, "multicluster_subclonal_like")} multicluster subclonal-like.

Project-level limited VAF summaries:

{vaf_project_lines}

This Level 3 result is not a full copy-number-aware clonal phylogeny. It is a descriptive VAF-cluster prototype restricted to copy-neutral candidate regions and does not infer branching order.

Key Level 3 outputs:

- `results/tables/level3_candidate_selection_qc_summary.tsv`
- `results/tables/level3_candidate_summary_by_project.tsv`
- `results/tables/level3_pilot_cohort.tsv`
- `results/tables/level3_segment_cn_qc_summary.tsv`
- `results/tables/level3_mutation_local_cn_annotation_qc_summary.tsv`
- `results/tables/level3_clonal_input_preparation_qc_summary.tsv`
- `results/tables/level3_limited_vaf_clustering_qc_summary.tsv`
- `results/tables/level3_limited_vaf_clonal_complexity_by_sample.tsv`
- `results/tables/level3_limited_vaf_clonal_complexity_by_project.tsv`
- `results/figures/level3_limited_vaf_cluster_vaf_density_by_sample.pdf`
- `results/figures/level3_limited_vaf_complexity_heatmap.pdf`

## TCGA Sample-Level Reference Atlas

The hardened atlas contains {metric_value(atlas_qc, "number_of_samples_included")} TCGA tumor samples from {metric_value(atlas_qc, "number_of_patients_included")} patients and preserves all {metric_value(atlas_qc, "number_of_projects_represented")} TCGA projects. All {metric_value(atlas_qc, "number_of_raw_features")} candidate features passed the 30% missingness threshold, leaving {metric_value(atlas_qc, "number_of_retained_features")} retained features and no dropped features. The modeling matrix contains {metric_value(atlas_qc, "number_of_imputed_values")} median-imputed values.

The locked feature composition is 2 mutation-count proxies, 29 driver-gene mutation indicators, 29 driver-gene mutation counts, 2 purity/ploidy features, 4 copy-number/aneuploidy features, 11 expression pathway scores, and 20 expression principal components. PCA and nearest-neighbor models were built successfully. Euclidean distance is the recorded default, with cosine also available. Saved feature definitions, medians, means, standard deviations, PCA loadings, model metadata, and feature-contract hashes allow future query samples to be transformed without refitting the reference.

This atlas is an unsupervised molecular comparator, not a supervised tissue-of-origin classifier and not a clinically validated diagnostic resource.

## Clinical WES/WTS Projection Prototype

The synthetic `example_cup_case` demonstrates projection of already-computed, harmonized WES/WTS features into the locked TCGA atlas. The query supplied {section_metric_value(clinical_qc, "global", "number_of_clinical_features_supplied")} of {section_metric_value(clinical_qc, "global", "number_of_expected_reference_features")} expected features ({section_metric_value(clinical_qc, "global", "percent_feature_coverage")}% coverage). WES coverage was {wes_supplied}/{wes_expected} ({percentage(wes_supplied, wes_expected)}%), while WTS coverage was {wts_supplied}/{wts_expected} ({percentage(wts_supplied, wts_expected)}%). The 20 missing expression-PC features were imputed with locked TCGA medians, and 4 non-reference supplied features were ignored.

The saved PCA model was reused without refitting, and the saved Euclidean nearest-neighbor index returned {section_metric_value(clinical_qc, "global", "top_k_neighbors_returned")} neighbors. The query coordinates were PC1 = 0.9452 and PC2 = -2.0524. The returned project composition was:

{clinical_project_lines}

All {section_metric_value(clinical_qc, "global", "top_k_neighbors_returned")} neighbors mapped to the broad {first_value(clinical_major_group, "major_cancer_group")} group. These results are descriptive similarity weights, not tissue-of-origin probabilities. The case is artificial and demonstrates workflow behavior only.

Key projection outputs include `clinical_tcga_projection_report.md`, the harmonized and scaled query vectors, feature-missingness audit, TCGA neighbor tables, PCA coordinates, and four PDF figures under `results/clinical_projection/example_cup_case/`.

## CUP Molecular Interpretation Prototype

The CUP interpretation layer aggregates project-level neighbor weights into cautious lineage groups, computes ambiguity metrics, summarizes eight evidence categories, applies optional differential-diagnosis constraints, and reports a heuristic molecular-pathologic concordance class. It does not train or apply a supervised classifier.

For the synthetic example, the top lineage was `{top_lineage}` with score {top_lineage_score} and `{top_confidence}`. The second lineage was `{second_lineage}` with score {second_lineage_score}. The ambiguity class was `{first_value(cup_ambiguity, "ambiguity_class")}`, with a top-1/top-2 margin of {format_float(first_value(cup_ambiguity, "top1_top2_margin"))} and normalized entropy of {format_float(first_value(cup_ambiguity, "entropy_like_score"))}. Overall, WES, and WTS feature coverage were {percentage(first_value(cup_ambiguity, "feature_coverage"), 1)}%, {percentage(first_value(cup_ambiguity, "wes_feature_coverage"), 1)}%, and {percentage(first_value(cup_ambiguity, "wts_feature_coverage"), 1)}%, respectively.

The evidence table contains: {list_join(context["cup_evidence_types"])}. It records both supporting and contradictory observations and treats median-imputed features as missing evidence rather than observed evidence. Molecular-pathologic concordance was `{first_value(cup_concordance, "concordance_class")}` because the synthetic submitted diagnosis did not specify a tissue lineage.

The CUP scores are normalized descriptive similarity weights, not diagnostic probabilities. Confidence rules are transparent but heuristic and uncalibrated. The report is a research/prototype interpretation that must be integrated with morphology, immunophenotype, imaging, and clinical findings.

## Current PyClone-VI/PhyloWGS Status

The allele-specific CN readiness assessment evaluated {metric_value(ascn_qc, "pilot_samples_evaluated")} pilot samples, including {metric_value(ascn_qc, "limited_vaf_pilot_samples_evaluated")} limited-VAF pilot samples. It found {metric_value(ascn_qc, "candidate_local_allele_specific_cn_files_found")} local allele-specific CN files, {metric_value(ascn_qc, "samples_with_major_cn")} samples with major CN, and {metric_value(ascn_qc, "samples_with_minor_cn")} samples with minor CN. Consequently, {metric_value(ascn_qc, "samples_ready_for_pyclone_vi")} samples are ready for PyClone-VI and {metric_value(ascn_qc, "samples_ready_for_phylowgs")} are ready for PhyloWGS; {metric_value(ascn_qc, "samples_requiring_external_allele_specific_cn_inference")} require external allele-specific CN inference.

`segment_mean` is not allele-specific integer copy number. It does not provide validated `total_cn`, `major_cn`, `minor_cn`, mutation multiplicity, or LOH state. PyClone-VI and PhyloWGS therefore remain paused. The recommended path is to obtain paired tumor/normal BAMs or validated allele-count inputs for the 30 limited-VAF pilot samples, run one allele-specific CN method consistently, standardize accepted outputs, reannotate mutations with major/minor CN, and only then prepare PyClone-VI inputs. PhyloWGS should be considered later on a high-confidence subset, with the single-bulk-sample limitation stated explicitly.

No PyClone-VI or PhyloWGS analysis has been run in this project.

## Limitations

- TCGA is enriched for primary tumors and may not represent clinical metastatic, relapse, small-biopsy, or cytology specimens.
- TCGA sample types and processing differ from real-world clinical WES/WTS samples.
- Levels 1 and 2 are project-level molecular similarity analyses, not within-patient tumor phylogenies.
- Level 1 and Level 2 trees should not be interpreted as literal species-like evolutionary trees or directional ancestry between cancer types.
- Current mutation burden is a mutation-count proxy until callable territory is added.
- Expression pathway scores are simple first-pass signatures, not model-based pathway activity estimates.
- Level 3 lacks allele-specific major/minor integer copy number.
- Segment mean is not allele-specific copy number and cannot substitute for integer major/minor CN in CN-aware clonal tools.
- Limited VAF clustering cannot infer definitive branching order.
- Single bulk TCGA samples limit clonal evolution inference.
- The synthetic projection is not evidence of performance on real clinical specimens.
- Twenty expression PCs were imputed in the example, leaving WTS coverage at 35.484%.
- Nearest-neighbor and CUP scores are descriptive similarity weights, not diagnostic probabilities.
- CUP lineage and confidence rules are heuristic and uncalibrated.
- Clinical actionability integration is not implemented.
- Clinical classification, CUP support, or diagnostic use requires validation on known-primary cohorts and prospective safeguards.

## Recommended Next Steps

1. Validate the projection and CUP interpretation workflows on independent known-primary clinical WES/WTS cases, including assay-stratified feature coverage, top-k and metric sensitivity, and blinded molecular-pathologic concordance review.
2. Add separately validated actionability integration after the tissue-context workflow is benchmarked; no OncoKB, CIViC, or AMP interpretation is currently implemented.
3. Obtain paired tumor/normal or validated allele-count inputs for the 30 limited-VAF pilot samples and generate reviewed allele-specific integer CN before any PyClone-VI or PhyloWGS execution.
4. Calibrate ambiguity and abstention thresholds and evaluate TCGA primary-tumor bias using metastatic, small-biopsy, and low-purity validation cohorts.
5. Replace the pending conceptual figure and complete formal reference formatting before manuscript submission.

## Integrated Report Artifacts

- `results/reports/integrated_tcga_cancer_phylogeny_report.md`
- `results/reports/integrated_tcga_cancer_phylogeny_executive_summary.md`
- `results/reports/integrated_tcga_cancer_phylogeny_methods.md`
- `results/reports/integrated_tcga_cancer_phylogeny_limitations_next_steps.md`
- `results/reports/integrated_tcga_cancer_phylogeny_report.html`
- `results/tables/integrated_project_summary.tsv`
- `results/tables/integrated_project_outputs_index.tsv`
"""


def build_executive_summary(context: dict) -> str:
    tables = context["tables"]
    candidate_qc = tables["level3_candidate_qc"]
    vaf_qc = tables["level3_vaf_qc"]
    atlas_qc = tables["reference_atlas_qc"]
    clinical_qc = tables["clinical_qc"]
    cup_ranked = tables["cup_ranked"]
    cup_ambiguity = tables["cup_ambiguity"]
    ascn_qc = tables["level3_ascn_qc"]
    return f"""# {TITLE}: Executive Summary

This project provides an end-to-end computational proof of concept spanning TCGA pan-cancer molecular similarity, group-specific molecular trait dendrograms, limited within-tumor VAF clustering, a locked sample-level reference atlas, and research-only clinical WES/WTS and CUP projection.

Level 1 builds pan-cancer molecular similarity trees across {context["level1_projects"]} TCGA cancer types using {context["level1_features"]} biologic features. The Gower and correlation trees provide complementary views: Gower emphasizes absolute scaled trait differences, while correlation emphasizes relative feature-profile shape.

Level 2 builds group-specific molecular similarity trees. Trees were generated for {list_join(context["level2_tree_groups"])}. CNS/glial, melanocytic, hematolymphoid, and sarcoma project-level analyses remain underpowered in the current project-level design.

Level 3 evaluates within-tumor clonal structure in a limited prototype. Candidate selection screened {metric_value(candidate_qc, "total_mutation_layer_samples")} tumor samples, selected {metric_value(candidate_qc, "samples_selected_for_pilot_cohort")} pilot samples, prepared copy-neutral mutation sets, and clustered {metric_value(vaf_qc, "total_clustered_mutations")} mutations from {metric_value(vaf_qc, "total_selected_pilot_samples")} selected samples. Interpretation classes are {metric_value(vaf_qc, "predominantly_clonal_like")} predominantly clonal-like, {metric_value(vaf_qc, "oligoclonal_like")} oligoclonal-like, and {metric_value(vaf_qc, "multicluster_subclonal_like")} multicluster subclonal-like.

The main scientific caveat is that Levels 1 and 2 are molecular similarity and trait-mapping trees, not literal organismal phylogenies. Level 3 is closer to tumor evolution, but the current implementation is limited VAF-based clustering in copy-neutral candidate regions. It does not run PyClone-VI or PhyloWGS, does not infer allele-specific integer copy number, and does not infer definitive branching order.

The hardened sample-level atlas contains {metric_value(atlas_qc, "number_of_samples_included")} tumors from {metric_value(atlas_qc, "number_of_patients_included")} patients across {metric_value(atlas_qc, "number_of_projects_represented")} projects. It retains all {metric_value(atlas_qc, "number_of_retained_features")} features, records {metric_value(atlas_qc, "number_of_imputed_values")} median-imputed modeling values, and provides locked PCA and Euclidean/cosine nearest-neighbor artifacts.

The artificial `example_cup_case` supplied {section_metric_value(clinical_qc, "global", "number_of_clinical_features_supplied")} of {section_metric_value(clinical_qc, "global", "number_of_expected_reference_features")} features ({section_metric_value(clinical_qc, "global", "percent_feature_coverage")}% coverage). Its top CUP lineage was `{first_value(cup_ranked, "lineage_or_project_group")}` with score {format_float(first_value(cup_ranked, "similarity_score"))} and `{first_value(cup_ambiguity, "ambiguity_class")}`. These are descriptive similarity weights from a synthetic workflow demonstration, not diagnostic probabilities.

PyClone-VI and PhyloWGS remain paused: {metric_value(ascn_qc, "samples_ready_for_pyclone_vi")} and {metric_value(ascn_qc, "samples_ready_for_phylowgs")} samples, respectively, are ready because allele-specific integer major/minor CN is unavailable. The next priorities are known-primary clinical validation, later actionability integration, and external allele-specific CN acquisition for the Level 3 pilot.
"""


def build_methods_summary(context: dict) -> str:
    return f"""# Integrated Methods Summary

## Data Acquisition

The workflow retrieves or reads TCGA project metadata, MC3 somatic mutation calls, PanCanAtlas/GDC purity and ploidy resources, aneuploidy and arm-level copy-number calls, GDC segment-level copy-number files for Level 3 pilot samples, and PanCanAtlas expression data. Local source-file fallbacks are used where possible, and missing-source states are recorded explicitly in QC tables.

## Feature Engineering and Reference Atlas

Mutation features include mutation-count proxies and configured driver-gene mutation indicators and counts. Copy-number features include purity, ploidy, aneuploidy score, and arm gain/loss burden. Expression features include transparent first-pass pathway scores and expression principal components. Level 1 uses project-level aggregated traits. The locked sample-level atlas retains 97 features, uses TCGA-reference median imputation and standard scaling, and saves feature definitions, transformation parameters, PCA, and Euclidean/cosine nearest-neighbor artifacts.

## Level 1 Tree Construction

Level 1 reads the project-level feature matrix, removes annotation fields from tree calculations, median-imputes remaining numeric missing values, scales biologic features, and generates Gower and correlation average-linkage molecular similarity dendrograms. Bootstrap feature resampling summarizes project-pair co-clustering stability.

## Level 2 Group Definitions

Level 2 defines biologically motivated project groups, including carcinoma_all, pan_squamous, gi_pancancreatobiliary, kidney, and gynecologic_breast. It filters features for within-group variability, scales features within group, and builds Gower and correlation trees where project counts are sufficient. Underpowered groups are retained in QC tables rather than forced into trees.

## Level 3 Candidate Selection

Level 3 candidate selection evaluates mutation counts, ref/alt count availability, purity, ploidy, copy-number availability, and candidate project priorities. It selects a high-information pilot cohort for clonal-structure prototyping.

## Segment CN Processing

The GDC-backed segment-CN step queries open-access TCGA copy-number segment files for pilot projects, matches files to pilot samples where possible, selects one best sample-level file per pilot sample, parses segment intervals, and updates candidate/pilot copy-number availability tables.

## Mutation-to-Segment Annotation

The mutation local-CN annotation step overlaps somatic mutations with sample-level segment intervals, computes depth and observed VAF, records local segment coordinates and segment mean, and writes a non-final PyClone-style preview. Segment mean is treated conservatively and is not allele-specific integer copy number.

## Limited VAF Clustering

The limited VAF clustering prototype uses copy-neutral candidate mutations from selected pilot samples. It fits one-dimensional Gaussian mixture models to observed VAFs per sample, selects cluster counts by BIC with a conservative simpler-model tolerance, and falls back to quantile binning when needed. Cluster labels are descriptive VAF-structure summaries.

## Clinical WES/WTS Projection

The clinical projection module ingests already-computed harmonized feature tables rather than raw FASTQ, BAM, VCF, or RNA-seq data. Exact locked feature names are matched, absent features are filled with saved TCGA medians, and saved TCGA means and standard deviations are applied without refitting. The saved PCA model is reused, and nearest neighbors are queried from the compatible saved index with direct scaled-matrix distances as a fallback.

## CUP Molecular Interpretation

The CUP layer validates the top-k project and broad-group summaries, aggregates project weights into predefined cautious lineage groups, and computes top-lineage share, top-1/top-2 margin, normalized entropy, neighbor concentration, feature coverage, and broad-group agreement. Configurable heuristic rules summarize project, major-group, expression, driver, copy-number, microenvironment, contradictory, and missing evidence. Differential diagnosis and molecular-pathologic comparisons are keyword-based aids, not diagnostic inference.

## Allele-Specific CN Readiness

The readiness workflow searches supported local allele-specific CN outputs, validates integer total/major/minor copy-number fields, and measures mutation overlap. It does not infer major/minor CN from segment mean. No PyClone-VI or PhyloWGS execution is performed; those tools remain gated on reviewed allele-specific integer CN, purity, mutation read counts, and sufficient mutation overlap.

## Statistical Caveats

Mutation burden remains a count proxy rather than true mutations/Mb. Expression pathway scores are simple gene-set summaries. Levels 1 and 2 are project-level molecular similarity trees rather than literal evolutionary phylogenies. Level 3 does not infer allele-specific integer CN or branching order, and single bulk samples limit clonal-evolution inference. Clinical and CUP scores are descriptive similarity summaries, not probabilities, and have not been clinically calibrated.
"""


def build_limitations_next_steps(context: dict) -> str:
    return f"""# Limitations and Next Steps

## Key Limitations

- TCGA is mostly primary tumors and does not fully represent metastatic, relapse, small-biopsy, or cytology clinical specimens.
- Levels 1 and 2 summarize project-level molecular traits, not within-patient phylogeny.
- Level 1 and Level 2 trees are molecular similarity trees and should not be interpreted as literal species-like evolutionary phylogenies.
- Mutation burden is currently a mutation-count proxy without callable-territory normalization.
- Segment mean is not allele-specific major/minor integer copy number.
- Level 3 lacks allele-specific CN and therefore does not produce CN-aware PyClone-VI or PhyloWGS inputs.
- Limited VAF clustering cannot infer definitive branching order.
- Single bulk TCGA samples limit robust clonal-evolution inference.
- The example clinical query is synthetic, supplies only 11 of 31 WTS features, and median-imputes 20 expression PCs.
- Clinical projection and CUP scores are descriptive similarity weights, not tissue-of-origin probabilities.
- CUP confidence and evidence rules are heuristic and uncalibrated, and no supervised CUP classifier was trained.
- Clinical actionability is not implemented.
- Clinical use requires independent known-primary validation, assay harmonization, uncertainty calibration, and prospective reporting safeguards.

## Next Recommended Steps

1. Validate projection and CUP interpretation on independent known-primary WES/WTS cohorts, including metastatic and limited-specimen subsets.
2. Calibrate top-k, distance-metric, ambiguity, and abstention behavior before any supervised classifier development.
3. Add a separately validated actionability layer after the molecular-context workflow is benchmarked.
4. Obtain paired tumor/normal or validated allele-count inputs for the 30 limited-VAF pilot samples and run one allele-specific CN method consistently.
5. Prepare PyClone-VI inputs only after reviewed integer major/minor CN and mutation-overlap QC are available; consider PhyloWGS later on a high-confidence subset.
6. Create the pending conceptual manuscript figure, format references, and obtain external scientific and clinical review.
"""


def simple_markdown_to_html(markdown_text: str, title: str) -> str:
    lines = markdown_text.splitlines()
    out = [
        "<!doctype html>",
        "<html>",
        "<head>",
        '<meta charset="utf-8">',
        f"<title>{html.escape(title)}</title>",
        "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;line-height:1.5;max-width:980px;margin:40px auto;padding:0 24px;color:#202124}code{background:#f1f3f4;padding:1px 4px;border-radius:4px}h1,h2,h3{line-height:1.2}li{margin:4px 0}</style>",
        "</head>",
        "<body>",
    ]
    in_list = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            if in_list:
                out.append("</ul>")
                in_list = False
            continue
        if stripped.startswith("- "):
            if not in_list:
                out.append("<ul>")
                in_list = True
            out.append(f"<li>{html.escape(stripped[2:])}</li>")
            continue
        if in_list:
            out.append("</ul>")
            in_list = False
        if stripped.startswith("### "):
            out.append(f"<h3>{html.escape(stripped[4:])}</h3>")
        elif stripped.startswith("## "):
            out.append(f"<h2>{html.escape(stripped[3:])}</h2>")
        elif stripped.startswith("# "):
            out.append(f"<h1>{html.escape(stripped[2:])}</h1>")
        else:
            out.append(f"<p>{html.escape(stripped)}</p>")
    if in_list:
        out.append("</ul>")
    out.extend(["</body>", "</html>"])
    return "\n".join(out)


def run_report(root: Path) -> dict[str, Path]:
    paths = default_report_paths(root)
    context = build_context(root)
    main_report = build_main_report(root, context)
    executive_summary = build_executive_summary(context)
    methods = build_methods_summary(context)
    limitations = build_limitations_next_steps(context)

    write_text(paths.main_report, main_report)
    write_text(paths.executive_summary, executive_summary)
    write_text(paths.methods, methods)
    write_text(paths.limitations_next_steps, limitations)
    write_text(paths.html_report, simple_markdown_to_html(main_report, TITLE))

    project_summary = build_project_summary(context)
    write_tsv(project_summary, paths.project_summary)
    outputs_index = build_outputs_index(root, paths)
    outputs_index.loc[outputs_index["file_path"] == rel(paths.outputs_index, root), "exists"] = True
    write_tsv(outputs_index, paths.outputs_index)
    return {
        "main_report": paths.main_report,
        "executive_summary": paths.executive_summary,
        "methods": paths.methods,
        "limitations_next_steps": paths.limitations_next_steps,
        "html_report": paths.html_report,
        "project_summary": paths.project_summary,
        "outputs_index": paths.outputs_index,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate integrated TCGA cancer phylogeny project reports.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    outputs = run_report(root)
    logging.info("Integrated report generation complete: %s", ", ".join(rel(path, root) for path in outputs.values()))


if __name__ == "__main__":
    main()

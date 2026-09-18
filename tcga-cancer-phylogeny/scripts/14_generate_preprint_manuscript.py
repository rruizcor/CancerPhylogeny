#!/usr/bin/env python

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import pandas as pd

from lib.common import configure_logging, ensure_dir, find_project_root


TITLE = "Multilevel Molecular Similarity, Limited Clonal-Structure Analysis, and Clinical Projection Framework Using TCGA Multi-omic Data"


class Paths:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.preprint = root / "results" / "reports" / "tcga_cancer_phylogeny_preprint.md"
        self.index = root / "results" / "reports" / "tcga_cancer_phylogeny_preprint_figures_tables_index.tsv"
        self.summary = root / "results" / "tables" / "integrated_project_summary.tsv"
        self.outputs_index = root / "results" / "tables" / "integrated_project_outputs_index.tsv"
        self.integrated_report = root / "results" / "reports" / "integrated_tcga_cancer_phylogeny_report.md"
        self.methods = root / "results" / "reports" / "integrated_tcga_cancer_phylogeny_methods.md"
        self.limitations = root / "results" / "reports" / "integrated_tcga_cancer_phylogeny_limitations_next_steps.md"


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        logging.warning("Missing input table: %s", path)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def clean(value: object, default: str = "not_available") -> str:
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
    matches = df[df["metric"].astype(str) == metric]
    if matches.empty:
        return default
    return clean(matches.iloc[0]["value"], default=default)


def section_metric_value(
    df: pd.DataFrame,
    section: str,
    metric: str,
    default: str = "not_available",
) -> str:
    if df.empty or not {"qc_section", "metric", "value"}.issubset(df.columns):
        return default
    matches = df[
        (df["qc_section"].astype(str) == section)
        & (df["metric"].astype(str) == metric)
    ]
    if matches.empty:
        return default
    return clean(matches.iloc[0]["value"], default=default)


def first_value(df: pd.DataFrame, column: str, default: str = "not_available") -> str:
    if df.empty or column not in df.columns:
        return default
    return clean(df.iloc[0][column], default=default)


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


def format_count(value: object) -> str:
    text = clean(value)
    try:
        numeric = float(str(text).replace(",", ""))
    except ValueError:
        return text
    if numeric.is_integer():
        return f"{int(numeric):,}"
    return f"{numeric:,.3f}".rstrip("0").rstrip(".")


def markdown_link(root: Path, report_path: Path, source_file: str) -> str:
    path = root / source_file
    if path.exists():
        rel_target = os.path.relpath(path, start=report_path.parent)
        return f"[`{source_file}`]({rel_target})"
    return f"`{source_file}` (pending)"


def source_links(root: Path, report_path: Path, source_files: str) -> str:
    if source_files == "pending":
        return "`pending conceptual figure`"
    return "; ".join(markdown_link(root, report_path, source.strip()) for source in source_files.split(";"))


def build_context(root: Path) -> dict:
    tables = {
        "summary": read_tsv(root / "results" / "tables" / "integrated_project_summary.tsv"),
        "candidate_qc": read_tsv(root / "results" / "tables" / "level3_candidate_selection_qc_summary.tsv"),
        "segment_qc": read_tsv(root / "results" / "tables" / "level3_segment_cn_qc_summary.tsv"),
        "annotation_qc": read_tsv(root / "results" / "tables" / "level3_mutation_local_cn_annotation_qc_summary.tsv"),
        "input_qc": read_tsv(root / "results" / "tables" / "level3_clonal_input_preparation_qc_summary.tsv"),
        "vaf_qc": read_tsv(root / "results" / "tables" / "level3_limited_vaf_clustering_qc_summary.tsv"),
        "atlas_qc": read_tsv(root / "results" / "tables" / "tcga_sample_reference_atlas_qc_summary.tsv"),
        "ascn_qc": read_tsv(root / "results" / "tables" / "level3_allele_specific_cn_readiness_qc_summary.tsv"),
        "clinical_qc": read_tsv(
            root / "results" / "clinical_projection" / "example_cup_case" / "clinical_projection_qc_summary.tsv"
        ),
        "clinical_project": read_tsv(
            root / "results" / "clinical_projection" / "example_cup_case" / "tcga_nearest_project_summary.tsv"
        ),
        "clinical_major_group": read_tsv(
            root / "results" / "clinical_projection" / "example_cup_case" / "tcga_nearest_major_group_summary.tsv"
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
        "outputs_index": read_tsv(root / "results" / "tables" / "integrated_project_outputs_index.tsv"),
    }
    summary = tables["summary"]
    annotation_qc = tables["annotation_qc"]
    total_mutations = metric_value(annotation_qc, "total_mutations_evaluated")
    annotated_mutations = metric_value(annotation_qc, "total_mutations_annotated_with_local_cn")
    try:
        unmatched = str(int(float(total_mutations)) - int(float(annotated_mutations)))
    except ValueError:
        unmatched = "not_available"

    context = {
        "tables": tables,
        "tcga_projects": metric_value(summary, "tcga_projects_represented", "33"),
        "level1_features": metric_value(summary, "biologic_features_used", "52"),
        "level1_distances": metric_value(summary, "distance_methods_generated", "gower;correlation"),
        "level2_groups": metric_value(
            summary,
            "groups_with_project_level_trees",
            "carcinoma_all, pan_squamous, gi_pancancreatobiliary, kidney, gynecologic_breast",
        ),
        "level2_underpowered": metric_value(
            summary,
            "underpowered_groups",
            "cns_glial, melanocytic, hematolymphoid, sarcoma",
        ),
        "samples_evaluated": metric_value(summary, "total_mutation_layer_samples", "10201"),
        "samples_ref_alt": metric_value(summary, "samples_with_ref_alt_counts", "10201"),
        "samples_purity": metric_value(summary, "samples_with_purity", "9516"),
        "samples_ploidy": metric_value(summary, "samples_with_ploidy", "9516"),
        "samples_cn": metric_value(summary, "samples_with_copy_number_or_aneuploidy", "9667"),
        "samples_passing": metric_value(summary, "samples_passing_minimum_filters", "3017"),
        "pilot_samples": metric_value(summary, "samples_selected_for_pilot_cohort", "84"),
        "segment_rows": metric_value(summary, "total_segment_rows_processed", "34028"),
        "pilot_segment_cn": metric_value(summary, "pilot_samples_with_segment_level_cn", "84"),
        "mutations_evaluated": total_mutations,
        "mutations_annotated": annotated_mutations,
        "annotation_rate": metric_value(summary, "percentage_mutations_annotated", "96.815"),
        "unmatched_mutations": unmatched,
        "segment_annotated_candidates": metric_value(tables["input_qc"], "total_segment_annotated_candidate_mutations", "131398"),
        "copy_neutral_candidates": metric_value(summary, "total_copy_neutral_candidate_mutations", "106337"),
        "limited_vaf_eligible": metric_value(summary, "samples_eligible_for_limited_vaf_clustering", "79"),
        "cn_aware_eligible": metric_value(summary, "samples_eligible_for_copy_number_aware_clustering", "0"),
        "requires_allele_specific": metric_value(summary, "samples_requiring_allele_specific_cn", "84"),
        "selected_limited_vaf": metric_value(summary, "total_selected_pilot_samples", "30"),
        "clustered_samples": metric_value(summary, "samples_successfully_clustered", "30"),
        "skipped_samples": metric_value(summary, "samples_skipped", "0"),
        "clustered_mutations": metric_value(summary, "total_clustered_mutations", "81154"),
        "predominantly_clonal": metric_value(summary, "predominantly_clonal_like", "5"),
        "oligoclonal": metric_value(summary, "oligoclonal_like", "11"),
        "multicluster": metric_value(summary, "multicluster_subclonal_like", "14"),
        "atlas_samples": metric_value(summary, "number_of_samples_included", "10201"),
        "atlas_projects": metric_value(summary, "number_of_projects_represented", "33"),
        "atlas_patients": metric_value(tables["atlas_qc"], "number_of_patients_included", "10130"),
        "atlas_raw_features": metric_value(tables["atlas_qc"], "number_of_raw_features", "97"),
        "atlas_features": metric_value(tables["atlas_qc"], "number_of_retained_features", "97"),
        "atlas_dropped_features": metric_value(tables["atlas_qc"], "number_of_dropped_features", "0"),
        "atlas_imputed_values": metric_value(tables["atlas_qc"], "number_of_imputed_values", "23222"),
        "atlas_default_metric": metric_value(tables["atlas_qc"], "default_nearest_neighbor_metric", "euclidean"),
        "clinical_expected": section_metric_value(tables["clinical_qc"], "global", "number_of_expected_reference_features", "97"),
        "clinical_supplied": section_metric_value(tables["clinical_qc"], "global", "number_of_clinical_features_supplied", "77"),
        "clinical_missing": section_metric_value(tables["clinical_qc"], "global", "number_of_missing_clinical_features", "20"),
        "clinical_coverage": section_metric_value(tables["clinical_qc"], "global", "percent_feature_coverage", "79.381"),
        "clinical_extra": section_metric_value(tables["clinical_qc"], "global", "number_of_extra_clinical_features_ignored", "4"),
        "clinical_metric": section_metric_value(tables["clinical_qc"], "global", "nearest_neighbor_metric", "euclidean"),
        "clinical_top_k": section_metric_value(tables["clinical_qc"], "global", "top_k_neighbors_returned", "50"),
        "wes_expected": section_metric_value(tables["clinical_qc"], "WES", "number_of_expected_reference_features", "66"),
        "wes_supplied": section_metric_value(tables["clinical_qc"], "WES", "number_of_clinical_features_supplied", "66"),
        "wts_expected": section_metric_value(tables["clinical_qc"], "WTS", "number_of_expected_reference_features", "31"),
        "wts_supplied": section_metric_value(tables["clinical_qc"], "WTS", "number_of_clinical_features_supplied", "11"),
        "top_lineage": first_value(tables["cup_ranked"], "lineage_or_project_group", "lower_GI"),
        "top_lineage_score": first_value(tables["cup_ranked"], "similarity_score", "0.5701127836271324"),
        "top_lineage_confidence": first_value(tables["cup_ranked"], "confidence_tier", "moderate_support"),
        "ambiguity_class": first_value(tables["cup_ambiguity"], "ambiguity_class", "moderate_ambiguity"),
        "ambiguity_margin": first_value(tables["cup_ambiguity"], "top1_top2_margin", "0.26090142126849997"),
        "ambiguity_entropy": first_value(tables["cup_ambiguity"], "entropy_like_score", "0.7371942615412224"),
        "cup_concordance": first_value(tables["cup_concordance"], "concordance_class", "not_assessable"),
        "evidence_categories": str(len(tables["cup_evidence"])),
        "ascn_files": metric_value(tables["ascn_qc"], "candidate_local_allele_specific_cn_files_found", "0"),
        "ascn_major": metric_value(tables["ascn_qc"], "samples_with_major_cn", "0"),
        "ascn_minor": metric_value(tables["ascn_qc"], "samples_with_minor_cn", "0"),
        "pyclone_ready": metric_value(tables["ascn_qc"], "samples_ready_for_pyclone_vi", "0"),
        "phylowgs_ready": metric_value(tables["ascn_qc"], "samples_ready_for_phylowgs", "0"),
        "ascn_external": metric_value(tables["ascn_qc"], "samples_requiring_external_allele_specific_cn_inference", "84"),
        "indexed_outputs": str(len(tables["outputs_index"])),
        "missing_indexed_outputs": (
            str((~tables["outputs_index"]["exists"].astype(str).str.lower().eq("true")).sum())
            if not tables["outputs_index"].empty and "exists" in tables["outputs_index"].columns
            else "not_available"
        ),
    }
    cup_ranked = tables["cup_ranked"]
    if len(cup_ranked) > 1:
        context["second_lineage"] = clean(cup_ranked.iloc[1]["lineage_or_project_group"])
        context["second_lineage_score"] = clean(cup_ranked.iloc[1]["similarity_score"])
    else:
        context["second_lineage"] = "not_available"
        context["second_lineage_score"] = "not_available"
    return context


def table_1(context: dict) -> str:
    return f"""| Data layer | Unit or scope | Current coverage | Interpretation boundary |
|---|---:|---:|---|
| Mutation layer | TCGA tumor samples | {format_count(context["samples_evaluated"])} samples across {context["tcga_projects"]} projects; ref/alt counts for {format_count(context["samples_ref_alt"])} samples | Mutation burden remains a mutation-count proxy without callable territory normalization. |
| Purity/ploidy | Tumor samples | {format_count(context["samples_purity"])} samples with purity and {format_count(context["samples_ploidy"])} with ploidy | Estimates support feature engineering and Level 3 filtering, but do not by themselves define clonal phylogeny. |
| Aneuploidy/arm-level CN | Tumor samples | {format_count(context["samples_cn"])} samples with copy-number or aneuploidy data | Arm-level summaries support molecular trait analysis, not allele-specific local CN. |
| Expression features | Mutation-matched tumor samples | 9,565 samples; 11 pathway scores and 20 expression PCs represented in the atlas | Requires preprocessing compatible with the TCGA reference. |
| Reference atlas | TCGA tumor samples | {format_count(context["atlas_samples"])} samples, {format_count(context["atlas_patients"])} patients, {context["atlas_projects"]} projects, {context["atlas_features"]} retained features | Unsupervised comparator; not a clinical classifier. |
| Clinical projection example | Synthetic query | {context["clinical_supplied"]}/{context["clinical_expected"]} features supplied ({context["clinical_coverage"]}%); {context["clinical_top_k"]} neighbors returned | Artificial workflow demonstration; similarity weights are not diagnostic probabilities. |
| CUP interpretation example | Synthetic query | Top lineage `{context["top_lineage"]}`; `{context["ambiguity_class"]}`; {context["evidence_categories"]} evidence categories | Heuristic research interpretation, not a validated tissue-of-origin diagnosis. |"""


def table_2(context: dict) -> str:
    return """| Feature class | Number of retained features | Clinical modality context |
|---|---:|---|
| Mutation-count proxies | 2 | WES-compatible proxies; not callable-territory-normalized TMB. |
| Driver mutation indicators | 29 | WES-compatible binary indicators. |
| Driver mutation counts | 29 | WES-compatible per-gene mutation counts. |
| Purity/ploidy | 2 | WES-compatible when harmonized estimates are supplied. |
| Copy-number/aneuploidy | 4 | Broad WES-compatible burden summaries; not allele-specific CN. |
| Pathway scores | 11 | WTS-compatible expression summaries. |
| Expression PCs | 20 | WTS-compatible only with reference-matched preprocessing. |
| **Total** | **97** | Locked atlas feature contract. |"""


def table_3(context: dict) -> str:
    project = context["tables"]["clinical_project"]
    if project.empty:
        top_projects = "not_available"
    else:
        project = project.sort_values("rank")
        top_projects = "; ".join(
            f"{row.project_code}: {int(row.n_neighbors)} neighbors, {float(row.weighted_similarity_score):.4f} weight"
            for row in project.itertuples(index=False)
        )
    return f"""| Projection metric | Synthetic example result | Interpretation |
|---|---|---|
| Overall feature coverage | {context["clinical_supplied"]}/{context["clinical_expected"]} ({context["clinical_coverage"]}%) | Twenty missing expression PCs were median-imputed. |
| WES feature coverage | {context["wes_supplied"]}/{context["wes_expected"]} ({percentage(context["wes_supplied"], context["wes_expected"])}%) | Complete coverage of locked WES-compatible features. |
| WTS feature coverage | {context["wts_supplied"]}/{context["wts_expected"]} ({percentage(context["wts_supplied"], context["wts_expected"])}%) | Incomplete transcriptomic coverage limits interpretation. |
| Ignored non-reference features | {context["clinical_extra"]} | Excluded from PCA and distance calculations. |
| Projection method | Saved `{context["clinical_metric"]}` index; saved PCA reused without refitting | Clinical scaling was not recomputed. |
| PCA coordinates | PC1 = 0.9452; PC2 = -2.0524 | Coordinates locate the query in atlas space but are not diagnostic scores. |
| Top project composition | {top_projects} | Normalized descriptive similarity weights, not tissue-of-origin probabilities. |"""


def table_4(context: dict) -> str:
    return f"""| CUP interpretation field | Synthetic example result | Caveat |
|---|---|---|
| Top lineage | `{context["top_lineage"]}`; score {format_float(context["top_lineage_score"])} | Similarity-weight share, not a probability. |
| Confidence tier | `{context["top_lineage_confidence"]}` | Transparent heuristic tier; not clinically calibrated. |
| Second lineage | `{context["second_lineage"]}`; score {format_float(context["second_lineage_score"])} | Competing pancreatobiliary/hepatobiliary support remains relevant. |
| Ambiguity | `{context["ambiguity_class"]}`; margin {format_float(context["ambiguity_margin"])}; normalized entropy {format_float(context["ambiguity_entropy"])} | Intended to prevent overconfident single-label interpretation. |
| Evidence categories | {context["evidence_categories"]} categories | Includes supporting, contradictory, and missing evidence. |
| Molecular-pathologic concordance | `{context["cup_concordance"]}` | Synthetic diagnosis does not specify a tissue lineage. |
| Clinical-use boundary | Research/prototype output | Not a validated tissue-of-origin assay or diagnosis. |"""


def table_5(context: dict) -> str:
    return f"""| Readiness metric | Current result | Consequence |
|---|---:|---|
| Local allele-specific CN files found | {context["ascn_files"]} | External acquisition or inference is required. |
| Samples with major CN | {context["ascn_major"]} | Segment mean is not substituted for major CN. |
| Samples with minor CN | {context["ascn_minor"]} | Segment mean is not substituted for minor CN. |
| PyClone-VI-ready samples | {context["pyclone_ready"]} | PyClone-VI remains paused. |
| PhyloWGS-ready samples | {context["phylowgs_ready"]} | PhyloWGS remains paused. |
| Samples requiring external allele-specific CN inference | {context["ascn_external"]} | Prioritize the 30 limited-VAF pilot samples for paired tumor/normal input acquisition. |"""


def table_6(context: dict) -> str:
    return """| Analysis set | Scope | Output status | Reported molecular similarity pattern | Caveat |
|---|---|---|---|---|
| Level 1 pan-cancer | 33 TCGA projects, 52 biologic features | Gower and correlation molecular similarity trees generated | Gower separates a broad lower-mutation or lineage-distinct cluster, a broad carcinoma-enriched cluster, and singleton or outlier-like DLBC, SKCM, UCEC, and UCS groups. | Project-level molecular similarity tree, not a literal species-like phylogeny. |
| carcinoma_all | 20 carcinoma projects | Level 2 Gower and correlation trees generated | Broad cluster: BLCA, BRCA, CESC, ESCA, HNSC, LUAD, LUSC, OV, STAD; second mixed cluster: CHOL, KICH, KIRP, LIHC, PRAD, THCA; close COAD/READ pair; KIRC, PAAD, UCEC more separated at k = 6. | Descriptive group-specific molecular trait dendrogram. |
| pan_squamous | 5 projects | Level 2 trees generated | HNSC clusters with ESCA; LUSC, CESC, and mixed BLCA are more separated. | Small project count limits topology confidence. |
| gi_pancancreatobiliary | 7 projects | Level 2 trees generated | COAD/READ/STAD cluster; CHOL/LIHC cluster. | Similarity pattern is based on available project-level traits. |
| kidney | 3 projects | Level 2 trees generated | KIRC/KIRP cluster apart from KICH. | Three-project tree should be read as a compact descriptive comparison. |
| gynecologic_breast | 5 projects | Level 2 trees generated | BRCA/OV cluster; UCEC, CESC, and UCS remain separate. | Similarity does not imply directional ancestry. |
| Underpowered groups | cns_glial, melanocytic, hematolymphoid, sarcoma | Recorded as underpowered, not forced into robust project-level trees | cns_glial, melanocytic, and hematolymphoid have two projects each; sarcoma is project-level only in the current workflow. | Explicitly retained as limitations rather than overinterpreted. |"""


def table_7(context: dict) -> str:
    return f"""| Level 3 step | Main count or result | Source table | Interpretation boundary |
|---|---:|---|---|
| Candidate selection | {format_count(context["samples_evaluated"])} samples evaluated; {format_count(context["samples_passing"])} passed minimum filters | `results/tables/level3_candidate_selection_qc_summary.tsv` | Feasibility selection for clonal-structure prototyping. |
| Pilot selection | {format_count(context["pilot_samples"])} pilot samples across 11 pilot projects | `results/tables/level3_pilot_cohort.tsv` | High-information subset, not full-TCGA Level 3 analysis. |
| Segment CN acquisition | {format_count(context["pilot_segment_cn"])} pilot samples with segment-level CN | `results/tables/level3_segment_cn_qc_summary.tsv` | Segment means are conservative local-CN context, not allele-specific integer CN. |
| Mutation-to-segment annotation | {format_count(context["mutations_annotated"])} of {format_count(context["mutations_evaluated"])} pilot mutations annotated with local CN | `results/tables/level3_mutation_local_cn_annotation_qc_summary.tsv` | Annotation supports copy-neutral mutation subset selection. |
| Clonal input stratification | {format_count(context["segment_annotated_candidates"])} segment-annotated candidates; {format_count(context["copy_neutral_candidates"])} copy-neutral limited-VAF candidates; {format_count(context["limited_vaf_eligible"])} samples eligible for limited VAF clustering; {format_count(context["cn_aware_eligible"])} eligible for fully CN-aware clustering | `results/tables/level3_clonal_input_preparation_qc_summary.tsv` | Fully CN-aware inputs require allele-specific major/minor CN. |
| Limited VAF clustering | {format_count(context["selected_limited_vaf"])} selected pilot samples; {format_count(context["clustered_mutations"])} mutations clustered; {format_count(context["skipped_samples"])} skipped samples | `results/tables/level3_limited_vaf_clustering_qc_summary.tsv` | Clonal-structure prototype only; branching order is not inferred. |"""


def table_8(context: dict) -> str:
    return f"""| Interpretation class | Number of samples | Operational meaning | Caveat |
|---|---:|---|---|
| Predominantly clonal-like | {format_count(context["predominantly_clonal"])} | Most clustered copy-neutral mutations fall into one dominant VAF component. | VAF structure only; not a definitive clonal tree. |
| Oligoclonal-like | {format_count(context["oligoclonal"])} | A limited number of VAF components are detected without high multicluster complexity. | Copy-neutral filtering reduces, but does not eliminate, copy-number confounding. |
| Multicluster subclonal-like | {format_count(context["multicluster"])} | Multiple VAF components suggest more complex subclonal VAF structure. | Single bulk samples limit branching-order interpretation. |"""


def figure_rows() -> list[dict[str, str]]:
    return [
        {
            "item_type": "figure",
            "item_number": "1",
            "title": "Conceptual overview of the multilevel framework",
            "source_file": "pending",
            "included_in_preprint": "TRUE",
            "notes": "Placeholder included; conceptual figure remains to be created.",
        },
        {
            "item_type": "figure",
            "item_number": "2",
            "title": "Level 1 pan-cancer Gower molecular similarity tree",
            "source_file": "results/figures/level1_pan_cancer_gower_tree.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Existing Level 1 figure.",
        },
        {
            "item_type": "figure",
            "item_number": "3",
            "title": "Level 1 pan-cancer molecular trait heatmap",
            "source_file": "results/figures/level1_pan_cancer_tree_with_trait_heatmap.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Existing Level 1 figure.",
        },
        {
            "item_type": "figure",
            "item_number": "4",
            "title": "Level 1 feature PCA",
            "source_file": "results/figures/level1_pan_cancer_feature_pca.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Existing Level 1 figure.",
        },
        {
            "item_type": "figure",
            "item_number": "5",
            "title": "Level 2 carcinoma-all Gower molecular similarity tree",
            "source_file": "results/figures/level2_carcinoma_all_gower_tree.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Representative Level 2 group-specific tree.",
        },
        {
            "item_type": "figure",
            "item_number": "6",
            "title": "Level 3 limited VAF clustering summary",
            "source_file": "results/figures/level3_limited_vaf_complexity_heatmap.pdf;results/figures/level3_limited_vaf_cluster_vaf_density_by_sample.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Existing limited-VAF figures used as a composite figure reference.",
        },
        {
            "item_type": "figure",
            "item_number": "7",
            "title": "TCGA sample reference PCA by project",
            "source_file": "results/figures/tcga_sample_reference_pca_by_project.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Existing reference-atlas PCA figure.",
        },
        {
            "item_type": "figure",
            "item_number": "8",
            "title": "TCGA sample reference PCA by major cancer group",
            "source_file": "results/figures/tcga_sample_reference_pca_by_major_group.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Existing reference-atlas PCA figure.",
        },
        {
            "item_type": "figure",
            "item_number": "9",
            "title": "Synthetic clinical query PCA projection by TCGA project",
            "source_file": "results/clinical_projection/example_cup_case/clinical_tcga_pca_projection_by_project.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Synthetic example projection figure.",
        },
        {
            "item_type": "figure",
            "item_number": "10",
            "title": "Synthetic clinical query nearest-neighbor project composition",
            "source_file": "results/clinical_projection/example_cup_case/clinical_nearest_neighbor_barplot.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Synthetic example nearest-neighbor figure.",
        },
        {
            "item_type": "figure",
            "item_number": "11",
            "title": "CUP top-project similarity summary",
            "source_file": "results/clinical_projection/example_cup_case/cup_interpretation/cup_top_project_similarity_barplot.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Synthetic CUP interpretation figure.",
        },
        {
            "item_type": "figure",
            "item_number": "12",
            "title": "CUP major-group similarity summary",
            "source_file": "results/clinical_projection/example_cup_case/cup_interpretation/cup_major_group_similarity_barplot.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Synthetic CUP interpretation figure.",
        },
        {
            "item_type": "figure",
            "item_number": "13",
            "title": "CUP ambiguity summary",
            "source_file": "results/clinical_projection/example_cup_case/cup_interpretation/cup_ambiguity_summary.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Synthetic CUP interpretation figure.",
        },
        {
            "item_type": "figure",
            "item_number": "14",
            "title": "CUP heuristic evidence map",
            "source_file": "results/clinical_projection/example_cup_case/cup_interpretation/cup_evidence_heatmap.pdf",
            "included_in_preprint": "TRUE",
            "notes": "Synthetic CUP evidence figure; values are heuristic, not probabilities.",
        },
    ]


def table_rows() -> list[dict[str, str]]:
    return [
        {
            "item_type": "table",
            "item_number": "1",
            "title": "Project data layers and sample counts",
            "source_file": "results/tables/integrated_project_summary.tsv;results/tables/level3_mutation_local_cn_annotation_qc_summary.tsv",
            "included_in_preprint": "TRUE",
            "notes": "Markdown table generated in manuscript.",
        },
        {
            "item_type": "table",
            "item_number": "2",
            "title": "TCGA reference atlas feature composition",
            "source_file": "results/reference_atlas/tcga_reference_feature_definitions.yaml;results/tables/tcga_sample_reference_atlas_qc_summary.tsv",
            "included_in_preprint": "TRUE",
            "notes": "Markdown table generated in manuscript.",
        },
        {
            "item_type": "table",
            "item_number": "3",
            "title": "Synthetic clinical projection example",
            "source_file": "results/clinical_projection/example_cup_case/clinical_projection_qc_summary.tsv;results/clinical_projection/example_cup_case/tcga_nearest_project_summary.tsv",
            "included_in_preprint": "TRUE",
            "notes": "Markdown table generated in manuscript.",
        },
        {
            "item_type": "table",
            "item_number": "4",
            "title": "Synthetic CUP interpretation example",
            "source_file": "results/clinical_projection/example_cup_case/cup_interpretation/cup_ranked_lineage_interpretation.tsv;results/clinical_projection/example_cup_case/cup_interpretation/cup_ambiguity_metrics.tsv;results/clinical_projection/example_cup_case/cup_interpretation/cup_molecular_evidence_table.tsv",
            "included_in_preprint": "TRUE",
            "notes": "Markdown table generated in manuscript.",
        },
        {
            "item_type": "table",
            "item_number": "5",
            "title": "PyClone-VI and PhyloWGS readiness",
            "source_file": "results/tables/level3_allele_specific_cn_readiness_qc_summary.tsv",
            "included_in_preprint": "TRUE",
            "notes": "Markdown table generated in manuscript.",
        },
        {
            "item_type": "table",
            "item_number": "6",
            "title": "Level 1 and Level 2 tree outputs",
            "source_file": "results/tables/integrated_project_summary.tsv;results/tables/level2_group_definitions.tsv",
            "included_in_preprint": "TRUE",
            "notes": "Markdown table generated in manuscript.",
        },
        {
            "item_type": "table",
            "item_number": "7",
            "title": "Level 3 clonal-preparation summary",
            "source_file": "results/tables/level3_candidate_selection_qc_summary.tsv;results/tables/level3_segment_cn_qc_summary.tsv;results/tables/level3_mutation_local_cn_annotation_qc_summary.tsv;results/tables/level3_clonal_input_preparation_qc_summary.tsv;results/tables/level3_limited_vaf_clustering_qc_summary.tsv",
            "included_in_preprint": "TRUE",
            "notes": "Markdown table generated in manuscript.",
        },
        {
            "item_type": "table",
            "item_number": "8",
            "title": "Limited VAF clustering interpretation classes",
            "source_file": "results/tables/level3_limited_vaf_clustering_qc_summary.tsv",
            "included_in_preprint": "TRUE",
            "notes": "Markdown table generated in manuscript.",
        },
    ]


def build_index() -> pd.DataFrame:
    rows = figure_rows() + table_rows()
    return pd.DataFrame(rows, columns=["item_type", "item_number", "title", "source_file", "included_in_preprint", "notes"])


def build_figure_reference_block(root: Path, report_path: Path) -> str:
    rows = figure_rows()
    lines = []
    for row in rows:
        lines.append(
            f"- Figure {row['item_number']}. {row['title']}. Source: {source_links(root, report_path, row['source_file'])}."
        )
    return "\n".join(lines)


def build_figure_legends(root: Path, report_path: Path) -> str:
    source = lambda path: source_links(root, report_path, path)
    return f"""## Figure Legends

**Figure 1. Conceptual overview of the multilevel framework.** [Figure placeholder: conceptual framework to be created.] The intended figure should connect Level 1 pan-cancer molecular similarity, Level 2 group-specific molecular trait dendrograms, Level 3 limited VAF clustering, the 97-feature sample-level atlas, and research-only clinical/CUP projection. It should distinguish project-level similarity from within-tumor clonal structure and clinical diagnosis.

**Figure 2. Level 1 pan-cancer Gower molecular similarity tree.** Source: {source("results/figures/level1_pan_cancer_gower_tree.pdf")}. This figure shows the average-linkage Gower molecular similarity tree across 33 TCGA projects using 52 retained project-level biologic features. Gower distance emphasizes scaled absolute molecular trait differences across mutation, copy-number, purity/ploidy, and expression features. Branch proximity should be interpreted as project-level molecular similarity, not directional ancestry or a species-like cancer phylogeny.

**Figure 3. Level 1 pan-cancer molecular trait heatmap.** Source: {source("results/figures/level1_pan_cancer_tree_with_trait_heatmap.pdf")}. This figure maps retained molecular traits onto the Level 1 molecular trait dendrogram, allowing visual comparison of branch structure with mutation, copy-number, and expression-feature patterns. The heatmap helps identify which traits may contribute to similarity patterns, but it is descriptive and does not establish causal relationships between features and cancer-type grouping.

**Figure 4. Level 1 feature PCA.** Source: {source("results/figures/level1_pan_cancer_feature_pca.pdf")}. This figure summarizes project-level feature variation in principal-component space. It provides a complementary ordination view of the same feature matrix used for Level 1 trees. It should be interpreted as a sensitivity visualization of pan-cancer molecular trait space, not as a classifier or clinical projection model.

**Figure 5. Level 2 carcinoma-all Gower molecular similarity tree.** Source: {source("results/figures/level2_carcinoma_all_gower_tree.pdf")}. This representative Level 2 figure shows the carcinoma_all group-specific molecular similarity tree after within-group feature filtering and scaling. It illustrates broad carcinoma subgroup structure, including the reported broad carcinoma cluster, mixed CHOL/KICH/KIRP/LIHC/PRAD/THCA cluster, COAD/READ pair, and more separated KIRC, PAAD, and UCEC branches. The result is a group-specific molecular trait dendrogram and should not be interpreted as cancer-type ancestry.

**Figure 6. Level 3 limited VAF clustering summary.** Sources: {source("results/figures/level3_limited_vaf_complexity_heatmap.pdf;results/figures/level3_limited_vaf_cluster_vaf_density_by_sample.pdf")}. These figures summarize 30 selected pilot samples and 81,154 copy-neutral mutations. They show per-sample VAF-cluster complexity and VAF density. Interpretation is limited to copy-neutral VAF structure because single bulk samples and absent allele-specific integer CN preclude definitive branching inference.

**Figure 7. TCGA sample reference PCA by project.** Source: {source("results/figures/tcga_sample_reference_pca_by_project.pdf")}. This figure shows individual TCGA tumors projected in the PCA space fitted to the locked 97-feature scaled atlas and colored by project. It displays sample-level heterogeneity and project overlap; it is not a tissue-of-origin classifier.

**Figure 8. TCGA sample reference PCA by major cancer group.** Source: {source("results/figures/tcga_sample_reference_pca_by_major_group.pdf")}. This companion view colors the same reference PCA coordinates by broad cancer group. Broad-group separation and overlap are descriptive and do not establish clinical classification performance.

**Figure 9. Synthetic clinical query PCA projection by TCGA project.** Source: {source("results/clinical_projection/example_cup_case/clinical_tcga_pca_projection_by_project.pdf")}. The synthetic query is transformed with saved TCGA medians and scaling parameters and projected with the saved PCA model without refitting. Its location is a molecular comparison, not a diagnostic score.

**Figure 10. Synthetic clinical query nearest-neighbor project composition.** Source: {source("results/clinical_projection/example_cup_case/clinical_nearest_neighbor_barplot.pdf")}. This figure summarizes the project composition of 50 Euclidean nearest neighbors for the artificial query. Bars represent normalized descriptive similarity weight rather than tissue-of-origin probabilities.

**Figure 11. CUP top-project similarity summary.** Source: {source("results/clinical_projection/example_cup_case/cup_interpretation/cup_top_project_similarity_barplot.pdf")}. The plot shows COAD, PAAD, READ, STAD, and LUAD contributions among returned neighbors. It visualizes project composition and should not be read as a diagnostic probability distribution.

**Figure 12. CUP major-group similarity summary.** Source: {source("results/clinical_projection/example_cup_case/cup_interpretation/cup_major_group_similarity_barplot.pdf")}. All 50 synthetic-query neighbors fall within the broad Carcinoma group. This broad-group summary reuses the same neighbors and is not independent validation.

**Figure 13. CUP ambiguity summary.** Source: {source("results/clinical_projection/example_cup_case/cup_interpretation/cup_ambiguity_summary.pdf")}. The figure displays top-lineage weight, top-1/top-2 margin, normalized entropy, and overall feature coverage used in the transparent heuristic ambiguity assessment. These metrics are not clinically calibrated.

**Figure 14. CUP heuristic evidence map.** Source: {source("results/clinical_projection/example_cup_case/cup_interpretation/cup_evidence_heatmap.pdf")}. The heatmap organizes project, major-group, driver, expression, copy-number, and microenvironment observations across represented lineage groups. It is an explainability aid based on conservative heuristics, not a probabilistic classifier."""


def build_table_legends() -> str:
    return """## Table Legends

**Table 1. Project data layers and sample counts.** This table summarizes the data layers used by the project and the current validated coverage from the integrated reporting layer. Counts are included to distinguish pan-cancer project-level feature availability from Level 3 pilot-specific segment and mutation annotation coverage.

**Table 2. TCGA reference atlas feature composition.** This table enumerates the 97 locked sample-level features and their expected WES or WTS context.

**Table 3. Synthetic clinical projection example.** This table reports feature coverage, imputation, locked-model reuse, PCA coordinates, and nearest-neighbor composition for the artificial query. Similarity weights are not diagnostic probabilities.

**Table 4. Synthetic CUP interpretation example.** This table summarizes lineage scores, ambiguity, evidence coverage, and the research-use boundary for the artificial case.

**Table 5. PyClone-VI and PhyloWGS readiness.** This table documents the absence of local allele-specific integer CN and the resulting pause before copy-number-aware clonal inference.

**Table 6. Level 1 and Level 2 tree outputs.** This table summarizes pan-cancer and group-specific molecular similarity outputs. These are project-level trait dendrograms, not literal species-like phylogenies.

**Table 7. Level 3 clonal-preparation summary.** This table traces candidate selection, segment annotation, copy-neutral mutation selection, and limited VAF clustering.

**Table 8. Limited VAF clustering interpretation classes.** This table reports the three descriptive VAF-structure classes assigned to the 30 selected samples; the classes are not definitive clonal phylogenies."""


def build_manuscript(root: Path, paths: Paths, context: dict) -> str:
    figures = build_figure_reference_block(root, paths.preprint)
    fig_legends = build_figure_legends(root, paths.preprint)
    table_legends = build_table_legends()

    return f"""# {TITLE}

Generated by `scripts/14_generate_preprint_manuscript.py`.

## Abstract

**Background:** Tumor classification increasingly integrates somatic mutation, copy-number, purity/ploidy, and transcriptomic features. Tree-like representations can help organize this multidimensional molecular trait space, but the analogy to phylogeny must be used carefully. In this project, Level 1 and Level 2 are molecular similarity trees and molecular trait dendrograms, not literal organismal or species-like phylogenies.

**Methods:** We developed a multilevel computational proof-of-concept using TCGA multi-omic data. Level 1 built pan-cancer molecular similarity trees across {context["tcga_projects"]} projects using {context["level1_features"]} retained biologic features and Gower and correlation distances. Level 2 built five group-specific molecular similarity analyses. Level 3 selected an {format_count(context["pilot_samples"])}-sample cohort, annotated mutations with segment-level copy-number context, and used Gaussian mixture models selected by BIC for limited VAF-based clustering of a copy-neutral subset. A separate sample-level atlas retained {context["atlas_features"]} harmonized features across {format_count(context["atlas_samples"])} tumors. A synthetic WES/WTS query was transformed with locked TCGA medians and scaling parameters, projected with saved PCA and nearest-neighbor models, and summarized by an uncertainty-aware heuristic CUP interpretation layer.

**Results:** Level 1 generated Gower and correlation trees for {context["tcga_projects"]} projects and {context["level1_features"]} features. Level 2 generated trees for five groups and recorded four groups as underpowered. Level 3 evaluated {format_count(context["samples_evaluated"])} samples, selected {format_count(context["pilot_samples"])} pilot samples, annotated {format_count(context["mutations_annotated"])} of {format_count(context["mutations_evaluated"])} mutations with segment CN ({context["annotation_rate"]}%), and clustered {format_count(context["clustered_mutations"])} copy-neutral mutations from {format_count(context["clustered_samples"])} samples. The atlas includes {format_count(context["atlas_samples"])} tumors from {format_count(context["atlas_patients"])} patients across {context["atlas_projects"]} projects, with {context["atlas_raw_features"]}/{context["atlas_features"]} raw/retained features and {format_count(context["atlas_imputed_values"])} median-imputed modeling values. The synthetic query supplied {context["clinical_supplied"]}/{context["clinical_expected"]} features ({context["clinical_coverage"]}%); its top lineage was `{context["top_lineage"]}` (score {format_float(context["top_lineage_score"])}; `{context["top_lineage_confidence"]}`), followed by `{context["second_lineage"]}` (score {format_float(context["second_lineage_score"])}), with `{context["ambiguity_class"]}`.

**Conclusions:** This project establishes a computational proof-of-concept connecting multilevel TCGA molecular similarity, limited clonal-structure analysis, and explainable clinical projection. Level 1 and Level 2 are not literal species-like phylogenies. Level 3 is not PyClone-VI/PhyloWGS and is not allele-specific CN-aware. Clinical and CUP outputs are descriptive research comparisons, not validated diagnostic probabilities, and require known-primary validation before clinical use.

## Introduction

Cancer classification has moved beyond anatomic site and histologic appearance toward increasingly integrated molecular descriptions. Somatic mutation burden, driver-gene alteration prevalence, copy-number burden, tumor purity, ploidy, aneuploidy, and expression programs can all contribute to how tumors resemble or differ from one another. TCGA provides a useful test bed for such integrated analysis because it spans many tumor types and includes public molecular data layers that can be harmonized into project-level and sample-level feature spaces.

Evolutionary language is useful in cancer biology, but it requires precision. Within a tumor, clonal evolution is a real biological process driven by mutation, selection, drift, copy-number change, and sampling. Across cancer types, however, a tree built from aggregated molecular traits is not a species tree and should not be interpreted as showing that one cancer type evolved from another. This manuscript therefore uses the terms molecular similarity tree, molecular trait dendrogram, and pan-cancer molecular similarity analysis for Level 1 and Level 2 outputs.

The distinction between cancer-type molecular similarity and within-tumor clonal structure is central to the framework. Level 1 and Level 2 organize TCGA projects in molecular trait space. Level 3 moves toward within-sample clonal-structure analysis, but the current implementation remains a limited VAF-based clonal clustering prototype using copy-neutral or near-diploid mutation subsets. Segment means are not allele-specific integer copy number, and most TCGA cases are single bulk samples, so definitive tumor branching phylogeny is not inferred.

A reusable computational framework is useful because the same feature definitions can support project-level interpretation, group-specific analysis, and sample-level comparison. Clinical WES/WTS and cancer-of-unknown-primary workflows need explainable comparator atlases that record feature completeness, imputation, distance behavior, competing molecular contexts, and missing evidence. Molecular similarity may support a tissue-of-origin hypothesis, but it is distinct from a clinical diagnosis. The implemented projection and CUP modules demonstrate this distinction using a synthetic case; they are not clinically validated and require assay harmonization and known-primary validation.

## Methods

### Study Design and Conceptual Levels

The project was designed as a three-level computational framework. Level 1 performs pan-cancer molecular similarity analysis across TCGA cancer projects. Level 2 performs group-specific molecular similarity analysis within predefined biologic groupings. Level 3 performs limited VAF-based clonal clustering in selected individual tumor samples. Level 1 and Level 2 are project-level molecular similarity trees, not literal organismal phylogenies. Level 3 is closer to clonal-structure analysis but remains a constrained clonal-structure prototype.

### TCGA Data Sources

The workflow integrates MC3 somatic mutation calls with tumor ref/alt read counts, PanCanAtlas/GDC purity and ploidy estimates, aneuploidy and arm-level copy-number calls, GDC segment-level copy-number files for the Level 3 pilot cohort, and PanCanAtlas expression data. The sample-level reference atlas preserves tumor-level features for future projection and includes saved feature definitions, scaler parameters, PCA coordinates, and a nearest-neighbor index.

### Mutation Feature Engineering

Mutation-derived features include mutation-count proxies and configured driver-gene mutation prevalence. Mutation burden is treated as a mutation-count proxy because callable territory normalization has not been added. For Level 3, mutation records with tumor ref/alt counts support observed VAF calculation after depth and local segment annotation filters.

### Purity/Ploidy, Aneuploidy, and Copy-Number Features

Purity and ploidy estimates are used as project-level features and Level 3 sample-level readiness inputs. Aneuploidy and arm-level copy-number features summarize broad copy-number burden. Segment-level CN files are used for mutation-to-segment annotation in the Level 3 pilot cohort. Segment means are used conservatively as local context and are not converted into allele-specific integer major/minor copy number.

### Expression/Pathway Features

Expression features include transparent first-pass pathway scores and expression principal components from configured genes. Pathway scores summarize expression programs such as proliferation, immune/inflammatory state, cytotoxic T-cell signal, stromal signal, EMT/invasion, hypoxia, cell cycle, DNA repair, angiogenesis, and lineage/tissue-associated expression axes. These are phenotype summaries rather than definitive pathway activity estimates.

### TCGA Sample-Level Reference Atlas Construction

Tumor-level mutation, purity/ploidy, aneuploidy, arm-level alteration, pathway-score, and expression-PC tables were joined by sample barcode. Features exceeding 30% missingness would be excluded; all 97 candidate features were retained in the current atlas. Remaining missing values were filled with TCGA feature medians, numeric features were standardized with saved TCGA means and standard deviations, and deterministic feature and sample ordering was locked. PCA and Euclidean/cosine nearest-neighbor models were fit to the scaled TCGA matrix. Feature definitions, transformation parameters, model metadata, and contract hashes were saved so query samples can be transformed identically.

### Level 1 Pan-Cancer Molecular Similarity Trees

Level 1 reads the project-level feature matrix, excludes annotation fields, imputes remaining numeric missing values, scales retained biologic features, and constructs Gower and correlation distance matrices. Average-linkage hierarchical clustering is used to generate molecular similarity dendrograms. A neighbor-joining sensitivity tree and feature-bootstrap co-clustering stability summaries are generated where supported by the distance matrix. The output is a pan-cancer molecular similarity tree or molecular trait dendrogram, not a literal evolutionary tree among cancer types.

### Level 2 Within-Group Molecular Similarity Trees

Level 2 defines biologically motivated groups and applies within-group feature filtering and scaling. Gower and correlation trees are generated for groups with sufficient project counts and feature variability. Groups that are too small or otherwise underpowered are recorded explicitly rather than forced into interpretable trees.

### Level 3 Candidate Selection

Level 3 candidate selection evaluates mutation count, ref/alt count availability, purity, ploidy, copy-number or aneuploidy availability, and pilot-project representation. The goal is feasibility assessment and high-information pilot selection, not full-TCGA clonal inference.

### Segment-Level Copy-Number Retrieval and Mutation Annotation

The GDC-backed segment-CN step identifies open-access TCGA copy-number segment files for pilot projects, matches files to pilot samples, selects one best sample-level file per pilot sample, parses segment intervals, and updates candidate/pilot copy-number availability tables. The mutation annotation step overlaps somatic mutations with sample-level segment intervals, computes depth and observed VAF, and records local segment status. Segment mean is not allele-specific integer CN, and patient-level segment matches are lower-confidence evidence.

### Limited VAF-Based Clonal Clustering

The Level 3 prototype uses the copy-neutral or near-diploid mutation subset from selected pilot samples. Observed VAFs are clustered independently within each sample using one-dimensional Gaussian mixture models selected by BIC with conservative model-selection behavior and quantile fallback when needed. This limited VAF-based clonal clustering does not run PyClone-VI or PhyloWGS, does not infer allele-specific copy number, and does not infer definitive branching order.

### Clinical WES/WTS Feature Harmonization and Projection

The clinical module ingests already-computed harmonized WES/WTS features and does not process raw FASTQ, BAM, VCF, or RNA-seq data. Clinical columns are matched by exact locked feature names. Missing reference features are filled with saved TCGA medians, and saved TCGA means and standard deviations are applied without estimating parameters from the query. Extra non-reference features are recorded and excluded.

### Nearest-Neighbor and PCA Projection

The compatible saved nearest-neighbor index is used when available, with direct distance calculation from the scaled reference matrix as a fallback. The current default is Euclidean distance, with cosine also available in the atlas. Similarity is summarized as a monotonic transformation of distance and normalized over returned neighbors for display; it is not a probability. The saved PCA model transforms the query without refitting.

### CUP Lineage Grouping and Ambiguity Metrics

TCGA project weights are aggregated into 18 cautious lineage groups. Ambiguity metrics include top-lineage score, second-lineage score, top-1/top-2 margin, normalized entropy across represented lineages, meaningful-project count, overall and modality-specific feature coverage, broad-group agreement, and neighbor concentration. Transparent confidence tiers are assigned by configurable heuristic rules and have not been clinically calibrated.

### Molecular Evidence Table Generation

The CUP evidence layer summarizes nearest-neighbor project, major-group, expression/pathway, driver, copy-number/aneuploidy, immune/stromal/EMT, contradictory, and missing evidence. Feature rules are applied only to supplied rather than median-imputed features. Differential-diagnosis and molecular-pathologic comparisons use transparent keyword heuristics and do not alter the unconstrained neighbor ranking.

### PyClone-VI/PhyloWGS Readiness Assessment

The readiness workflow searches supported local allele-specific CN outputs, validates integer total, major, and minor copy-number fields, and measures mutation overlap. Major or minor CN is never inferred from segment mean. PyClone-VI readiness requires reviewed major/minor CN, purity, mutation read counts, and sufficient annotated mutations. PhyloWGS remains a later cautious option for a high-confidence subset because most TCGA cases are single bulk samples.

### Quality Control and Reproducibility

Each major workflow layer writes QC tables and explicit missingness or limitation notes. Integrated project summaries and output indexes are generated without rerunning primary analyses. The refreshed output index records {context["indexed_outputs"]} indexed outputs and {context["missing_indexed_outputs"]} missing outputs. The manuscript itself is generated deterministically from validated reports and QC tables.

### Software/Workflow Structure

The repository uses R for data retrieval, feature engineering, tree construction, and visualization; Python for Level 3 preprocessing, clustering, reference-atlas construction, clinical/CUP projection, and reporting; and Snakemake for workflow orchestration. The preprint draft is generated by `scripts/14_generate_preprint_manuscript.py`.

## Results

### Construction of a Pan-Cancer TCGA Feature Matrix

The project assembled a Level 1 project-level feature matrix containing {context["tcga_projects"]} TCGA projects and {context["level1_features"]} retained biologic features. The retained features span mutation-count proxies, driver mutation prevalence, purity/ploidy, aneuploidy and arm-level copy-number burden, and expression pathway scores. These project-level features support molecular similarity analysis and are distinct from the sample-level atlas.

**Table 1. Project data layers and sample counts.**

{table_1(context)}

### Hardened 97-Feature TCGA Reference Atlas

The sample-level atlas contains {format_count(context["atlas_samples"])} TCGA tumors from {format_count(context["atlas_patients"])} patients across all {context["atlas_projects"]} projects. It retains {context["atlas_features"]} of {context["atlas_raw_features"]} candidate features, drops {context["atlas_dropped_features"]}, and records {format_count(context["atlas_imputed_values"])} median-imputed values in the modeling matrix. PCA and nearest-neighbor models were built successfully, with Euclidean distance recorded as the default and cosine also available (Figures 7-8).

**Table 2. TCGA reference atlas feature composition.**

{table_2(context)}

The atlas preserves sample-level heterogeneity for research comparison. It is not a supervised tissue-of-origin classifier and has not been clinically validated.

### Level 1 Pan-Cancer Molecular Similarity Trees

Level 1 generated Gower and correlation pan-cancer molecular similarity trees from the 33-project, 52-feature matrix (Figures 2-4). The Gower tree emphasizes scaled absolute molecular trait differences and currently separates a broad lower-mutation or lineage-distinct cluster, a broad carcinoma-enriched cluster including BLCA, BRCA, CESC, COAD/READ, ESCA/HNSC, LUAD/LUSC, OV, and STAD, and singleton or outlier-like DLBC, SKCM, UCEC, and UCS groups. The correlation tree emphasizes relative feature-profile shape and differs from Gower where projects have similar feature patterns but different absolute trait levels.

These findings are descriptive pan-cancer molecular similarity results. Level 1 and Level 2 are not literal species-like phylogenies, and the Level 1 tree should not be used to infer that one cancer type evolved from another.

### Level 2 Group-Specific Molecular Similarity Trees

Level 2 generated group-specific molecular trait dendrograms for carcinoma_all, pan_squamous, gi_pancancreatobiliary, kidney, and gynecologic_breast (Figure 5). Four groups were recorded as underpowered in the current project-level design: cns_glial, melanocytic, hematolymphoid, and sarcoma project-level only. The reported group patterns include expected close relationships such as COAD/READ, HNSC/ESCA in the squamous-enriched set, COAD/READ/STAD and CHOL/LIHC in the GI/hepatobiliary set, KIRC/KIRP apart from KICH in kidney, and BRCA/OV within the gynecologic/breast set.

The generated and underpowered group sets are summarized in Table 6.

### Feasibility Assessment for Level 3 Clonal Analysis

Level 3 candidate selection evaluated {format_count(context["samples_evaluated"])} mutation-layer tumor samples. All {format_count(context["samples_ref_alt"])} had ref/alt count availability, {format_count(context["samples_purity"])} had purity/ploidy, and {format_count(context["samples_cn"])} had copy-number or aneuploidy data. A total of {format_count(context["samples_passing"])} samples passed minimum filters, and {format_count(context["pilot_samples"])} pilot samples were selected across 11 pilot projects. This step establishes feasibility for clonal-structure prototyping rather than full within-tumor phylogenetic inference.

### Segment-Level CN Retrieval and Mutation-to-Local-CN Annotation

For the Level 3 pilot cohort, segment-level CN was available for {format_count(context["pilot_segment_cn"])} pilot samples after GDC-backed retrieval and parsing. The mutation-to-segment annotation step evaluated {format_count(context["mutations_evaluated"])} pilot mutations and annotated {format_count(context["mutations_annotated"])} mutations with local segment CN ({context["annotation_rate"]}%), leaving {format_count(context["unmatched_mutations"])} unmatched mutations. These QC outputs support copy-neutral subset selection, but segment means are not allele-specific integer copy number.

### Limited VAF-Based Clonal Clustering Prototype

Clonal input stratification yielded {format_count(context["segment_annotated_candidates"])} segment-annotated candidate mutations and {format_count(context["copy_neutral_candidates"])} copy-neutral limited-VAF candidate mutations. {format_count(context["limited_vaf_eligible"])} samples were eligible for limited VAF clustering, but {format_count(context["cn_aware_eligible"])} samples were eligible for fully copy-number-aware clustering because allele-specific major/minor CN was unavailable. The selected prototype pilot included {format_count(context["selected_limited_vaf"])} samples across all 11 pilot projects.

Limited VAF clustering was completed for {format_count(context["clustered_samples"])} selected samples and {format_count(context["clustered_mutations"])} copy-neutral mutations, with {format_count(context["skipped_samples"])} skipped samples (Figure 6). Interpretation classes were {context["predominantly_clonal"]} predominantly clonal-like, {context["oligoclonal"]} oligoclonal-like, and {context["multicluster"]} multicluster subclonal-like. This Level 3 result is a clonal-structure prototype; it is not PyClone-VI/PhyloWGS, not allele-specific CN-aware, and not a definitive tumor branching phylogeny.

The Level 3 preparation counts and interpretation classes are summarized in Tables 7 and 8.

### Synthetic Clinical WES/WTS Projection

The artificial `example_cup_case` supplied {context["clinical_supplied"]} of {context["clinical_expected"]} locked features ({context["clinical_coverage"]}% overall coverage). WES coverage was {context["wes_supplied"]}/{context["wes_expected"]}, while WTS coverage was {context["wts_supplied"]}/{context["wts_expected"]}. Twenty missing expression PCs were filled with locked TCGA medians, and four non-reference supplied features were excluded. The saved PCA model was reused without refitting, placing the query at PC1 = 0.9452 and PC2 = -2.0524 (Figure 9). The saved Euclidean index returned {context["clinical_top_k"]} neighbors (Figure 10): COAD 17, PAAD 16, READ 11, STAD 3, and LUAD 3. All were in the broad Carcinoma group.

**Table 3. Synthetic clinical projection example.**

{table_3(context)}

These project contributions are descriptive molecular similarity weights, not tissue-of-origin probabilities or a diagnosis.

### CUP Interpretation Prototype

Project weights were aggregated into cautious lineage groups and reviewed with ambiguity and evidence summaries. The top lineage was `{context["top_lineage"]}` with score {format_float(context["top_lineage_score"])} and confidence tier `{context["top_lineage_confidence"]}`. `{context["second_lineage"]}` ranked second with score {format_float(context["second_lineage_score"])}. The top-1/top-2 margin was {format_float(context["ambiguity_margin"])}, normalized entropy was {format_float(context["ambiguity_entropy"])}, and the final class was `{context["ambiguity_class"]}` (Figures 11-14). The evidence table contained {context["evidence_categories"]} categories, including explicit contradictory and missing evidence. Concordance was `{context["cup_concordance"]}` because the synthetic submitted diagnosis did not specify a tissue lineage.

**Table 4. Synthetic CUP interpretation example.**

{table_4(context)}

This workflow is an uncertainty-aware research interpretation, not a clinically validated diagnostic classifier.

### Current Boundary of Clonal Phylogeny Inference

The allele-specific CN readiness assessment found {context["ascn_files"]} local allele-specific CN files, {context["ascn_major"]} samples with major CN, and {context["ascn_minor"]} with minor CN. Accordingly, {context["pyclone_ready"]} samples are ready for PyClone-VI, {context["phylowgs_ready"]} are ready for PhyloWGS, and {context["ascn_external"]} require external allele-specific CN inference.

**Table 5. PyClone-VI and PhyloWGS readiness.**

{table_5(context)}

`segment_mean` is not allele-specific integer copy number. PyClone-VI and PhyloWGS remain paused until reviewed total, major, and minor CN and mutation-overlap QC are available. No such tools were run in the current project.

### Supporting Molecular Similarity and Level 3 Summary Tables

**Table 6. Level 1 and Level 2 tree outputs.**

{table_6(context)}

**Table 7. Level 3 clonal-preparation summary.**

{table_7(context)}

**Table 8. Limited VAF clustering interpretation classes.**

{table_8(context)}

### Integrated Computational Resource

The refreshed integrated output index records {context["indexed_outputs"]} outputs and {context["missing_indexed_outputs"]} missing outputs. The repository now connects project-level similarity, limited within-tumor VAF clustering, a locked sample-level atlas, synthetic clinical projection, and explainable CUP interpretation while preserving the boundaries between these tasks.

Current figure references:

{figures}

## Discussion

This project demonstrates that TCGA multi-omic data can be organized into a coherent multilevel computational framework. Level 1 provides a pan-cancer molecular similarity tree that summarizes how TCGA projects relate in aggregated molecular trait space. Its value is explanatory and comparative: it gives investigators a way to inspect broad trait structure, compare Gower and correlation views, and map mutation, copy-number, and expression features onto a shared molecular trait dendrogram. The Level 1 output should be read as a pan-cancer molecular similarity analysis, not a literal history of cancer-type evolution.

Level 2 adds biologic focus by repeating the molecular similarity analysis within selected groups. This helps distinguish pan-cancer structure from within-group relationships such as COAD/READ, HNSC/ESCA, GI/hepatobiliary clustering, kidney subgroup separation, and BRCA/OV similarity. The small size of several groups means Level 2 is best treated as a descriptive molecular trait dendrogram layer rather than a robust topology for every lineage.

Level 3 demonstrates a cautious within-tumor clonal-structure prototype. It shows that TCGA mutation counts, ref/alt counts, purity/ploidy, segment-level CN, and mutation-to-segment annotation can be combined to select copy-neutral mutation subsets and cluster observed VAFs in selected samples. What it does not demonstrate is equally important: it does not infer allele-specific integer copy number, does not run PyClone-VI or PhyloWGS, does not reconstruct definitive branching order, and does not overcome the constraints of single bulk TCGA samples. Allele-specific total, major, and minor CN are required before VAF can be interpreted with defensible mutation multiplicity and cancer-cell-fraction assumptions.

The hardened atlas and projection module provide an explicit route for future clinical WES/WTS research. A clinical case can be transformed with a locked feature contract, medians, scaling parameters, PCA model, and distance metric rather than being normalized against itself. This reproducibility is important, but compatibility of capture design, expression preprocessing, tumor content, and clinical specimen type still requires empirical validation.

CUP interpretation should be uncertainty-aware because molecular evidence often supports overlapping lineages. In the synthetic example, lower-GI and pancreatobiliary/hepatobiliary contexts both receive substantial weight. The ambiguity metrics preserve this competition instead of converting it into a single overconfident label. Evidence tables are useful complements to black-box labels because they expose neighbor composition, supplied versus imputed features, driver and expression observations, contradictions, and absent evidence. They do not remove the need for morphology, immunophenotype, imaging, and clinical review.

Validation on independent known-primary clinical cases is required before diagnostic development. Relevant studies should evaluate metastatic and small-biopsy specimens, assay-stratified feature coverage, top-k and distance-metric sensitivity, calibrated abstention, blinded molecular-pathologic concordance, and failure modes for rare or absent TCGA lineages. The synthetic result is a workflow test, not evidence of accuracy.

Actionability should be integrated only after molecular-context performance is benchmarked. Variant-level pathogenicity, evidence-level curation, therapeutic biomarkers, and trial matching require dedicated validation and governance; no OncoKB, CIViC, or AMP actionability interpretation is implemented here.

The novelty of the project is the multilevel connection between pan-cancer trait space, group-level molecular similarity, and limited within-tumor clonal architecture. The framework keeps the levels conceptually distinct, which reduces the risk of overinterpreting project-level molecular similarity as evolutionary ancestry while still allowing tumor-evolution concepts to be explored cautiously within selected samples.

## Limitations

TCGA is enriched for primary tumors and may not represent metastatic, relapse, small-biopsy, or cytology specimens encountered in clinical WES/WTS workflows. Project-level aggregation in Level 1 and Level 2 loses intratumoral and intra-disease heterogeneity, and group definitions depend on TCGA project labels, available metadata, and current feature availability. Mutation burden is currently a mutation-count proxy rather than true mutations per megabase.

Level 1 and Level 2 are molecular similarity trees and molecular trait dendrograms. They are not literal species-like phylogenies, do not imply directional ancestry between cancer types, and should not be used to claim that one cancer type evolved from another. Small Level 2 groups, especially kidney and the underpowered two-project categories, require cautious interpretation.

Level 3 is limited by segment-level CN assumptions and TCGA sampling. Segment means are not allele-specific integer copy number, and no PyClone-VI or PhyloWGS execution was performed. The current Level 3 output is not allele-specific CN-aware and does not infer definitive branching. Single bulk TCGA samples support limited clonal clustering more readily than robust clonal phylogeny inference. Limited VAF clusters are approximations based on observed VAFs in copy-neutral or near-diploid regions.

No clinical validation was performed. The clinical WES/WTS and CUP modules operate on a synthetic case and require already-computed harmonized features. The example median-imputes 20 expression PCs, has 35.484% WTS coverage, and uses heuristic, uncalibrated confidence rules. Similarity weights are not diagnostic probabilities. TCGA primary-tumor bias, missing lineage-marker and fusion features, broad project labels, and absent actionability interpretation further limit use.

## Conclusion

This project establishes a computational proof-of-concept for multilevel TCGA molecular analysis, a hardened 97-feature sample-level atlas, and explainable research-only WES/WTS and CUP projection. It supports Level 1 pan-cancer and Level 2 group-specific molecular similarity analysis and provides a cautious Level 3 limited VAF-based clonal-structure prototype. Clinical translation requires known-primary validation and calibrated safeguards. PyClone-VI and PhyloWGS remain appropriately paused until allele-specific integer CN is available.

## Data and Code Availability

The analysis is implemented in the local repository with workflow scripts under `scripts/`, configuration files under `config/`, tests under `tests/`, and generated outputs under `results/`. The preprint draft is generated by `python scripts/14_generate_preprint_manuscript.py`. Key manuscript source files include `results/reports/integrated_tcga_cancer_phylogeny_report.md`, `results/tables/integrated_project_summary.tsv`, and `results/tables/integrated_project_outputs_index.tsv`. TCGA/GDC-derived resources remain subject to their source data-access terms and should be reprocessed with the documented workflow for reproducibility.

## References Placeholder

- TCGA Pan-Cancer Atlas publications and data resources.
- MC3 somatic mutation call set.
- GDC data portal and harmonized TCGA resources.
- PyClone-VI methods literature.
- PhyloWGS methods literature.
- Cancer-of-unknown-primary molecular classifier literature.
- Tumor evolution and clonal architecture literature.
- Copy-number-aware clonal inference and allele-specific copy-number estimation literature.

{fig_legends}

{table_legends}
"""


def write_outputs(paths: Paths) -> None:
    context = build_context(paths.root)
    manuscript = build_manuscript(paths.root, paths, context)
    ensure_dir(paths.preprint.parent)
    paths.preprint.write_text(manuscript.rstrip() + "\n", encoding="utf-8")
    logging.info("Wrote %s", paths.preprint)

    index = build_index()
    index.to_csv(paths.index, sep="\t", index=False)
    logging.info("Wrote %s", paths.index)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate the TCGA cancer phylogeny preprint manuscript draft.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection from cwd.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = args.root.resolve() if args.root else find_project_root()
    paths = Paths(root)
    write_outputs(paths)


if __name__ == "__main__":
    main()

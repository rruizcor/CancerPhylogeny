#!/usr/bin/env python

from __future__ import annotations

import argparse
import gzip
import logging
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcga_cancer_phylogeny_matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


CONFIG = {
    "allow_patient_level_segment_matching": False,
    "min_total_depth": 20,
    "min_observed_vaf": 0.0,
    "neutral_segment_mean_min": -0.15,
    "neutral_segment_mean_max": 0.15,
    "candidate_nonsynonymous_only": True,
    "min_copy_neutral_candidate_mutations": 50,
    "min_median_depth": 20,
    "min_cn_aware_mutations": 50,
    "pilot_min_samples": 20,
    "pilot_max_samples": 30,
    "pilot_per_project_cap": 3,
    "priority_projects": ["SKCM", "LUAD", "LUSC", "UCEC", "COAD", "READ", "BLCA", "HNSC", "SARC", "UCS", "DLBC"],
}

COPY_NEUTRAL_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "mutation_id",
    "chromosome",
    "position",
    "Hugo_Symbol",
    "Variant_Classification",
    "HGVSp_Short",
    "t_ref_count",
    "t_alt_count",
    "total_depth",
    "observed_vaf",
    "purity",
    "ploidy",
    "local_segment_mean",
    "local_cn_status",
    "depth_filter_pass",
    "vaf_filter_pass",
    "copy_neutral_filter_pass",
    "notes",
]

SEGMENT_ANNOTATED_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "mutation_id",
    "chromosome",
    "position",
    "Hugo_Symbol",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "HGVSp_Short",
    "t_ref_count",
    "t_alt_count",
    "total_depth",
    "observed_vaf",
    "purity",
    "ploidy",
    "local_segment_mean",
    "local_cn_status",
    "segment_annotation_label",
    "copy_number_aware_status",
    "notes",
]

READINESS_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "n_total_mutations",
    "n_candidate_mutations_with_local_cn",
    "n_copy_neutral_candidate_mutations",
    "median_depth",
    "median_observed_vaf",
    "purity",
    "ploidy",
    "median_local_segment_mean",
    "pct_mutations_neutral_like",
    "pct_mutations_gain_like",
    "pct_mutations_loss_like",
    "eligible_for_limited_vaf_clustering",
    "eligible_for_copy_number_aware_clustering",
    "reason_not_copy_number_aware",
    "readiness_class",
    "notes",
]

PILOT_COLUMNS = [
    "sample_barcode",
    "project_code",
    "n_copy_neutral_candidate_mutations",
    "median_depth",
    "median_observed_vaf",
    "purity",
    "ploidy",
    "selection_rank",
    "selection_reason",
]


@dataclass
class ClonalInputPaths:
    annotated_mutations: Path
    mutation_local_cn_summary: Path
    pilot_with_segments: Path
    candidates_with_segments: Path
    purity_ploidy: Path
    copy_neutral_input: Path
    segment_annotated_input: Path
    readiness_by_sample: Path
    limited_vaf_pilot_samples: Path
    qc_summary: Path
    neutral_counts_by_sample_figure: Path
    neutral_counts_by_project_figure: Path
    vaf_distribution_figure: Path
    readiness_by_project_figure: Path


def default_paths(root: Path) -> ClonalInputPaths:
    return ClonalInputPaths(
        annotated_mutations=project_path("data", "processed", "clonal", "level3_mutations_with_local_cn.tsv.gz", root=root),
        mutation_local_cn_summary=project_path(
            "results", "tables", "level3_mutation_local_cn_summary_by_sample.tsv", root=root
        ),
        pilot_with_segments=project_path("results", "tables", "level3_pilot_cohort_with_segments.tsv", root=root),
        candidates_with_segments=project_path("results", "tables", "level3_candidate_samples_with_segments.tsv", root=root),
        purity_ploidy=project_path("data", "processed", "features", "purity_ploidy_by_sample.tsv", root=root),
        copy_neutral_input=project_path(
            "data", "processed", "clonal", "level3_copy_neutral_vaf_clustering_input.tsv.gz", root=root
        ),
        segment_annotated_input=project_path(
            "data", "processed", "clonal", "level3_segment_annotated_clonal_input.tsv.gz", root=root
        ),
        readiness_by_sample=project_path("results", "tables", "level3_clonal_input_readiness_by_sample.tsv", root=root),
        limited_vaf_pilot_samples=project_path(
            "results", "tables", "level3_limited_vaf_clustering_pilot_samples.tsv", root=root
        ),
        qc_summary=project_path("results", "tables", "level3_clonal_input_preparation_qc_summary.tsv", root=root),
        neutral_counts_by_sample_figure=project_path(
            "results", "figures", "level3_copy_neutral_candidate_counts_by_sample.pdf", root=root
        ),
        neutral_counts_by_project_figure=project_path(
            "results", "figures", "level3_copy_neutral_candidate_counts_by_project.pdf", root=root
        ),
        vaf_distribution_figure=project_path(
            "results", "figures", "level3_vaf_distribution_copy_neutral_candidates.pdf", root=root
        ),
        readiness_by_project_figure=project_path("results", "figures", "level3_readiness_class_by_project.pdf", root=root),
    )


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["NA", ""])


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def is_autosome_or_sex(chromosome: str) -> bool:
    text = str(chromosome)
    return text in {str(i) for i in range(1, 23)} | {"X", "Y"}


def normalize_local_cn_status(segment_mean: float | int | str | None, config: dict = CONFIG) -> str:
    value = pd.to_numeric(pd.Series([segment_mean]), errors="coerce").iloc[0]
    if pd.isna(value):
        return "unknown"
    if float(config["neutral_segment_mean_min"]) <= value <= float(config["neutral_segment_mean_max"]):
        return "neutral_like"
    if value > float(config["neutral_segment_mean_max"]):
        return "gain_like"
    return "loss_like"


def stable_mutation_id(data: pd.DataFrame) -> pd.Series:
    return (
        data["sample_barcode"].astype(str)
        + "|"
        + data["chromosome"].astype(str)
        + "|"
        + data["position"].round().astype("Int64").astype(str)
        + "|"
        + data["Reference_Allele"].astype(str)
        + "|"
        + data["Tumor_Seq_Allele2"].astype(str)
        + "|"
        + data["Hugo_Symbol"].astype(str)
    )


def has_major_minor_cn(data: pd.DataFrame) -> pd.Series:
    if "major_cn" not in data.columns or "minor_cn" not in data.columns:
        return pd.Series(False, index=data.index)
    return numeric(data["major_cn"]).notna() & numeric(data["minor_cn"]).notna()


def prepare_mutation_flags(annotated: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    data = annotated.copy()
    for col in ["t_ref_count", "t_alt_count", "total_depth", "observed_vaf", "purity", "ploidy", "local_segment_mean"]:
        if col in data.columns:
            data[col] = numeric(data[col])
    if "is_nonsynonymous" in data.columns:
        coding_mask = as_bool(data["is_nonsynonymous"])
    else:
        coding_mask = data["Variant_Classification"].astype(str).ne("Silent")
    data["mutation_id"] = stable_mutation_id(data)
    data["local_cn_status"] = data["local_segment_mean"].map(lambda x: normalize_local_cn_status(x, config))
    data["has_ref_alt_counts"] = data["t_ref_count"].notna() & data["t_alt_count"].notna()
    data["depth_filter_pass"] = data["has_ref_alt_counts"] & (data["total_depth"] >= int(config["min_total_depth"]))
    data["vaf_filter_pass"] = data["observed_vaf"].notna() & (data["observed_vaf"] > float(config["min_observed_vaf"]))
    data["local_cn_matched"] = data["local_cn_match_status"].astype(str).eq("matched")
    data["sample_level_cn_annotation"] = data["local_cn_match_notes"].astype(str).str.contains("sample_level", na=False)
    if bool(config.get("allow_patient_level_segment_matching", False)):
        data["cn_annotation_allowed"] = data["local_cn_matched"]
    else:
        data["cn_annotation_allowed"] = data["local_cn_matched"] & data["sample_level_cn_annotation"]
    data["copy_neutral_filter_pass"] = data["local_segment_mean"].between(
        float(config["neutral_segment_mean_min"]),
        float(config["neutral_segment_mean_max"]),
        inclusive="both",
    )
    data["coding_filter_pass"] = coding_mask if bool(config.get("candidate_nonsynonymous_only", True)) else True
    data["chromosome_filter_pass"] = data["chromosome"].map(is_autosome_or_sex)
    data["has_major_minor_cn"] = has_major_minor_cn(data)
    data["basic_clonal_candidate"] = (
        data["cn_annotation_allowed"]
        & data["depth_filter_pass"]
        & data["vaf_filter_pass"]
        & data["coding_filter_pass"]
        & data["chromosome_filter_pass"]
    )
    data["copy_neutral_candidate"] = data["basic_clonal_candidate"] & data["copy_neutral_filter_pass"]
    return data


def build_copy_neutral_input(data: pd.DataFrame) -> pd.DataFrame:
    out = data[data["copy_neutral_candidate"]].copy()
    out["notes"] = (
        "limited_vaf_based_clonal_clustering;"
        "copy_neutral_near_diploid_by_segment_mean;"
        "not_copy_number_aware;"
        "requires_allele_specific_cn_for_cn_aware_phylogeny"
    )
    for column in COPY_NEUTRAL_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    return out[COPY_NEUTRAL_COLUMNS]


def build_segment_annotated_input(data: pd.DataFrame) -> pd.DataFrame:
    out = data[data["basic_clonal_candidate"]].copy()
    out["segment_annotation_label"] = "segment_annotated_not_allele_specific"
    out["copy_number_aware_status"] = np.where(
        out["has_major_minor_cn"],
        "copy_number_aware_candidate_major_minor_cn_available",
        "not_copy_number_aware_requires_allele_specific_cn",
    )
    out["notes"] = (
        "exploratory_visualization_or_future_cn_aware_preparation;"
        "local_segment_mean_retained;"
        "major_minor_cn_not_inferred"
    )
    for column in SEGMENT_ANNOTATED_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    return out[SEGMENT_ANNOTATED_COLUMNS]


def build_readiness_table(data: pd.DataFrame, pilot: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    data = data.copy()
    data["neutral_like"] = data["local_cn_status"].eq("neutral_like") & data["cn_annotation_allowed"]
    data["gain_like"] = data["local_cn_status"].eq("gain_like") & data["cn_annotation_allowed"]
    data["loss_like"] = data["local_cn_status"].eq("loss_like") & data["cn_annotation_allowed"]
    grouped = data.groupby("sample_barcode", dropna=False)
    summary = grouped.agg(
        patient_barcode=("patient_barcode", "first"),
        project_id=("project_id", "first"),
        project_code=("project_code", "first"),
        n_total_mutations=("sample_barcode", "size"),
        n_candidate_mutations_with_local_cn=("basic_clonal_candidate", "sum"),
        n_copy_neutral_candidate_mutations=("copy_neutral_candidate", "sum"),
        median_depth=("total_depth", "median"),
        median_observed_vaf=("observed_vaf", "median"),
        purity=("purity", "first"),
        ploidy=("ploidy", "first"),
        median_local_segment_mean=("local_segment_mean", lambda x: x.dropna().median() if x.notna().any() else np.nan),
        n_neutral_like=("neutral_like", "sum"),
        n_gain_like=("gain_like", "sum"),
        n_loss_like=("loss_like", "sum"),
        n_major_minor_cn=("has_major_minor_cn", "sum"),
    ).reset_index()
    denom = summary["n_neutral_like"] + summary["n_gain_like"] + summary["n_loss_like"]
    summary["pct_mutations_neutral_like"] = np.where(denom > 0, 100 * summary["n_neutral_like"] / denom, 0)
    summary["pct_mutations_gain_like"] = np.where(denom > 0, 100 * summary["n_gain_like"] / denom, 0)
    summary["pct_mutations_loss_like"] = np.where(denom > 0, 100 * summary["n_loss_like"] / denom, 0)
    summary["eligible_for_limited_vaf_clustering"] = (
        (summary["n_copy_neutral_candidate_mutations"] >= int(config["min_copy_neutral_candidate_mutations"]))
        & (summary["median_depth"] >= float(config["min_median_depth"]))
        & summary["purity"].notna()
    )
    summary["eligible_for_copy_number_aware_clustering"] = (
        (summary["n_major_minor_cn"] >= int(config["min_cn_aware_mutations"]))
        & (summary["median_depth"] >= float(config["min_median_depth"]))
        & summary["purity"].notna()
        & summary["ploidy"].notna()
    )
    summary["reason_not_copy_number_aware"] = np.where(
        summary["eligible_for_copy_number_aware_clustering"],
        "major_minor_cn_available",
        "requires_allele_specific_cn",
    )

    readiness_class: list[str] = []
    notes: list[str] = []
    for row in summary.itertuples(index=False):
        if pd.isna(row.purity) or pd.isna(row.ploidy):
            readiness_class.append("excluded")
            notes.append("missing_purity_or_ploidy")
        elif pd.isna(row.median_depth) or float(row.median_depth) < float(config["min_median_depth"]):
            readiness_class.append("insufficient_depth")
            notes.append(f"median_depth_lt_{config['min_median_depth']}")
        elif int(row.n_copy_neutral_candidate_mutations) < int(config["min_copy_neutral_candidate_mutations"]):
            readiness_class.append("insufficient_neutral_mutations")
            notes.append(f"copy_neutral_candidates_lt_{config['min_copy_neutral_candidate_mutations']}")
        elif not bool(row.eligible_for_copy_number_aware_clustering):
            readiness_class.append("ready_for_limited_vaf_clustering")
            notes.append("limited_vaf_based_clonal_clustering;requires_allele_specific_cn_for_cn_aware")
        else:
            readiness_class.append("requires_allele_specific_cn")
            notes.append("copy_number_aware_inputs_need_review")
    summary["readiness_class"] = readiness_class
    summary["notes"] = notes

    pilot_meta = pilot[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates()
    summary = pilot_meta.merge(summary, on=["sample_barcode", "patient_barcode", "project_id", "project_code"], how="left")
    fill_zero_cols = [
        "n_total_mutations",
        "n_candidate_mutations_with_local_cn",
        "n_copy_neutral_candidate_mutations",
        "pct_mutations_neutral_like",
        "pct_mutations_gain_like",
        "pct_mutations_loss_like",
    ]
    for column in fill_zero_cols:
        summary[column] = summary[column].fillna(0)
    summary["eligible_for_limited_vaf_clustering"] = summary["eligible_for_limited_vaf_clustering"].fillna(False).astype(bool)
    summary["eligible_for_copy_number_aware_clustering"] = (
        summary["eligible_for_copy_number_aware_clustering"].fillna(False).astype(bool)
    )
    summary["reason_not_copy_number_aware"] = summary["reason_not_copy_number_aware"].fillna("requires_allele_specific_cn")
    summary["readiness_class"] = summary["readiness_class"].fillna("excluded")
    summary["notes"] = summary["notes"].fillna("no_annotated_mutations")
    return summary[READINESS_COLUMNS]


def select_limited_vaf_pilot(readiness: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    eligible = readiness[readiness["eligible_for_limited_vaf_clustering"]].copy()
    if eligible.empty:
        return pd.DataFrame(columns=PILOT_COLUMNS)
    eligible = eligible.sort_values(
        ["n_copy_neutral_candidate_mutations", "median_depth", "purity"],
        ascending=[False, False, False],
    )
    selected_rows: list[pd.Series] = []
    selected_samples: set[str] = set()
    project_counts: dict[str, int] = {}
    cap = int(config["pilot_per_project_cap"])

    for project in config["priority_projects"]:
        project_rows = eligible[(eligible["project_code"] == project) & (~eligible["sample_barcode"].isin(selected_samples))]
        if project_rows.empty:
            continue
        row = project_rows.iloc[0]
        selected_rows.append(row)
        selected_samples.add(str(row["sample_barcode"]))
        project_counts[project] = project_counts.get(project, 0) + 1

    for _, row in eligible.iterrows():
        if len(selected_rows) >= int(config["pilot_max_samples"]):
            break
        sample = str(row["sample_barcode"])
        project = str(row["project_code"])
        if sample in selected_samples or project_counts.get(project, 0) >= cap:
            continue
        selected_rows.append(row)
        selected_samples.add(sample)
        project_counts[project] = project_counts.get(project, 0) + 1

    selected = pd.DataFrame(selected_rows)
    if selected.empty:
        return pd.DataFrame(columns=PILOT_COLUMNS)
    selected = selected.sort_values(
        ["n_copy_neutral_candidate_mutations", "median_depth"],
        ascending=[False, False],
    ).reset_index(drop=True)
    selected["selection_rank"] = selected.index + 1
    selected["selection_reason"] = selected.apply(
        lambda row: (
            "limited_vaf_based_clonal_clustering_prototype;"
            f"project_diversity_cap_{cap};"
            f"copy_neutral_candidates={int(row['n_copy_neutral_candidate_mutations'])}"
        ),
        axis=1,
    )
    return selected[PILOT_COLUMNS]


def build_qc_summary(
    data: pd.DataFrame,
    copy_neutral: pd.DataFrame,
    readiness: pd.DataFrame,
    selected_pilot: pd.DataFrame,
) -> pd.DataFrame:
    eligible_limited = readiness[readiness["eligible_for_limited_vaf_clustering"]]
    eligible_cn = readiness[readiness["eligible_for_copy_number_aware_clustering"]]
    requires_cn = readiness[readiness["reason_not_copy_number_aware"] == "requires_allele_specific_cn"]
    rows = [
        ("global", "total_pilot_samples_evaluated", len(readiness)),
        ("global", "total_annotated_mutations", len(data)),
        ("global", "total_segment_annotated_candidate_mutations", int(data["basic_clonal_candidate"].sum())),
        ("global", "total_copy_neutral_candidate_mutations", len(copy_neutral)),
        ("global", "samples_eligible_for_limited_vaf_clustering", len(eligible_limited)),
        ("global", "samples_eligible_for_copy_number_aware_clustering", len(eligible_cn)),
        ("global", "samples_requiring_allele_specific_cn", len(requires_cn)),
        (
            "global",
            "median_copy_neutral_candidate_mutations_per_sample",
            float(readiness["n_copy_neutral_candidate_mutations"].median()) if not readiness.empty else 0,
        ),
        (
            "global",
            "projects_represented_among_limited_vaf_candidates",
            ";".join(sorted(eligible_limited["project_code"].astype(str).unique())) if not eligible_limited.empty else "none",
        ),
        (
            "global",
            "projects_represented_in_selected_limited_vaf_pilot",
            ";".join(sorted(selected_pilot["project_code"].astype(str).unique())) if not selected_pilot.empty else "none",
        ),
        (
            "global",
            "limitations",
            "limited_vaf_based_clonal_clustering_only;segment_annotated_not_allele_specific;not_copy_number_aware;requires_allele_specific_cn_for_pyClone_vi_or_phylowgs_cn_aware_inputs",
        ),
    ]
    for readiness_class, count in readiness["readiness_class"].value_counts().items():
        rows.append(("readiness_class", str(readiness_class), int(count)))
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def plot_neutral_counts_by_sample(readiness: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = readiness.sort_values("n_copy_neutral_candidate_mutations", ascending=False)
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.bar(data["sample_barcode"], data["n_copy_neutral_candidate_mutations"], color="#4c78a8")
    ax.axhline(y=CONFIG["min_copy_neutral_candidate_mutations"], color="#e45756", linestyle="--", linewidth=1)
    ax.set_ylabel("Copy-neutral candidate mutations")
    ax.set_title("Level 3 copy-neutral VAF clustering candidates by sample")
    ax.tick_params(axis="x", rotation=90, labelsize=6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_neutral_counts_by_project(readiness: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = (
        readiness.groupby("project_code", dropna=False)
        .agg(n_copy_neutral_candidate_mutations=("n_copy_neutral_candidate_mutations", "sum"), n_samples=("sample_barcode", "size"))
        .reset_index()
        .sort_values("n_copy_neutral_candidate_mutations", ascending=False)
    )
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.bar(data["project_code"], data["n_copy_neutral_candidate_mutations"], color="#54a24b")
    ax.set_ylabel("Copy-neutral candidate mutations")
    ax.set_title("Copy-neutral VAF clustering candidates by project")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_vaf_distribution(copy_neutral: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(9, 6))
    values = copy_neutral["observed_vaf"].dropna()
    if values.empty:
        ax.text(0.5, 0.5, "No copy-neutral candidate mutations", ha="center", va="center")
        ax.axis("off")
    else:
        ax.hist(values, bins=60, color="#4c78a8")
        ax.set_xlabel("Observed VAF")
        ax.set_ylabel("Mutations")
        ax.set_title("VAF distribution for copy-neutral candidate mutations")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_readiness_by_project(readiness: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    counts = readiness.groupby(["project_code", "readiness_class"], dropna=False).size().unstack(fill_value=0)
    counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(12, 7))
    bottom = np.zeros(len(counts))
    colors = {
        "ready_for_limited_vaf_clustering": "#54a24b",
        "insufficient_neutral_mutations": "#f58518",
        "requires_allele_specific_cn": "#b279a2",
        "insufficient_depth": "#eeca3b",
        "excluded": "#bab0ac",
    }
    for column in counts.columns:
        ax.bar(counts.index, counts[column], bottom=bottom, label=column, color=colors.get(column, None))
        bottom += counts[column].to_numpy()
    ax.set_ylabel("Pilot samples")
    ax.set_title("Level 3 clonal input readiness by project")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def run_clonal_input_preparation(paths: ClonalInputPaths, config: dict = CONFIG) -> dict[str, pd.DataFrame]:
    for path, description in [
        (paths.annotated_mutations, "mutation local-CN annotation table"),
        (paths.mutation_local_cn_summary, "mutation local-CN sample summary"),
        (paths.pilot_with_segments, "pilot cohort with segment CN"),
        (paths.candidates_with_segments, "candidate table with segment CN"),
        (paths.purity_ploidy, "purity/ploidy by sample"),
    ]:
        require_file(path, description)

    annotated = pd.read_csv(paths.annotated_mutations, sep="\t")
    pilot = read_tsv(paths.pilot_with_segments)
    data = prepare_mutation_flags(annotated, config)
    copy_neutral = build_copy_neutral_input(data)
    segment_annotated = build_segment_annotated_input(data)
    readiness = build_readiness_table(data, pilot, config)
    selected_pilot = select_limited_vaf_pilot(readiness, config)
    qc = build_qc_summary(data, copy_neutral, readiness, selected_pilot)

    write_tsv_gz(copy_neutral, paths.copy_neutral_input)
    write_tsv_gz(segment_annotated, paths.segment_annotated_input)
    write_tsv(readiness, paths.readiness_by_sample)
    write_tsv(selected_pilot, paths.limited_vaf_pilot_samples)
    write_tsv(qc, paths.qc_summary)
    plot_neutral_counts_by_sample(readiness, paths.neutral_counts_by_sample_figure)
    plot_neutral_counts_by_project(readiness, paths.neutral_counts_by_project_figure)
    plot_vaf_distribution(copy_neutral, paths.vaf_distribution_figure)
    plot_readiness_by_project(readiness, paths.readiness_by_project_figure)

    return {
        "copy_neutral": copy_neutral,
        "segment_annotated": segment_annotated,
        "readiness": readiness,
        "selected_pilot": selected_pilot,
        "qc": qc,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare honest Level 3 clonal input sets and readiness tables.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument("--min-total-depth", type=int, default=CONFIG["min_total_depth"])
    parser.add_argument("--neutral-segment-mean-min", type=float, default=CONFIG["neutral_segment_mean_min"])
    parser.add_argument("--neutral-segment-mean-max", type=float, default=CONFIG["neutral_segment_mean_max"])
    parser.add_argument("--min-copy-neutral-candidates", type=int, default=CONFIG["min_copy_neutral_candidate_mutations"])
    parser.add_argument("--pilot-max-samples", type=int, default=CONFIG["pilot_max_samples"])
    parser.add_argument("--pilot-per-project-cap", type=int, default=CONFIG["pilot_per_project_cap"])
    parser.add_argument("--include-silent", action="store_true", help="Include synonymous mutations in candidate sets.")
    parser.add_argument(
        "--allow-patient-level-segment-matching",
        action="store_true",
        help="Allow lower-confidence patient-level segment annotations. Off by default.",
    )
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    config = dict(CONFIG)
    config["min_total_depth"] = args.min_total_depth
    config["neutral_segment_mean_min"] = args.neutral_segment_mean_min
    config["neutral_segment_mean_max"] = args.neutral_segment_mean_max
    config["min_copy_neutral_candidate_mutations"] = args.min_copy_neutral_candidates
    config["pilot_max_samples"] = args.pilot_max_samples
    config["pilot_per_project_cap"] = args.pilot_per_project_cap
    config["candidate_nonsynonymous_only"] = not args.include_silent
    config["allow_patient_level_segment_matching"] = bool(args.allow_patient_level_segment_matching)
    outputs = run_clonal_input_preparation(default_paths(root), config)
    logging.info(
        "Clonal input preparation complete: %d copy-neutral candidates, %d segment-annotated candidates, %d selected pilot samples",
        len(outputs["copy_neutral"]),
        len(outputs["segment_annotated"]),
        len(outputs["selected_pilot"]),
    )


if __name__ == "__main__":
    main()

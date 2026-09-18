#!/usr/bin/env python

from __future__ import annotations

import argparse
import logging
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcga_cancer_phylogeny_matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pyarrow.parquet as pq

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


CONFIG = {
    "min_nonsynonymous_solid": 100,
    "min_nonsynonymous_hematolymphoid": 50,
    "min_purity": 0.30,
    "pilot_project_priority": [
        "SKCM",
        "UCEC",
        "COAD",
        "READ",
        "LUAD",
        "LUSC",
        "BLCA",
        "HNSC",
        "SARC",
        "UCS",
        "LAML",
        "DLBC",
    ],
    "pilot_per_project_cap": 8,
    "pilot_total_min": 50,
    "pilot_total_max": 100,
    "mutation_score_cap": 1000,
}


@dataclass
class Level3Paths:
    mutation_parquet: Path
    tmb_by_sample: Path
    purity_ploidy_by_sample: Path
    aneuploidy_by_sample: Path
    projects: Path
    cancer_group_map: Path
    segments_tsv: Path
    segments_parquet: Path
    arm_calls: Path
    candidate_samples: Path
    summary_by_project: Path
    pilot_cohort: Path
    excluded_samples: Path
    qc_summary: Path
    mutation_purity_figure: Path
    counts_by_project_figure: Path
    score_by_project_figure: Path
    pilot_overview_figure: Path


def default_paths(root: Path) -> Level3Paths:
    return Level3Paths(
        mutation_parquet=project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root=root),
        tmb_by_sample=project_path("data", "processed", "features", "tmb_by_sample.tsv", root=root),
        purity_ploidy_by_sample=project_path("data", "processed", "features", "purity_ploidy_by_sample.tsv", root=root),
        aneuploidy_by_sample=project_path("data", "processed", "features", "aneuploidy_by_sample.tsv", root=root),
        projects=project_path("data", "interim", "tcga_projects.tsv", root=root),
        cancer_group_map=project_path("config", "cancer_group_map.csv", root=root),
        segments_tsv=project_path("data", "processed", "copy_number", "segments_by_sample.tsv", root=root),
        segments_parquet=project_path("data", "processed", "copy_number", "segments_by_sample.parquet", root=root),
        arm_calls=project_path("data", "processed", "copy_number", "arm_level_calls_by_sample.tsv", root=root),
        candidate_samples=project_path("results", "tables", "level3_candidate_samples.tsv", root=root),
        summary_by_project=project_path("results", "tables", "level3_candidate_summary_by_project.tsv", root=root),
        pilot_cohort=project_path("results", "tables", "level3_pilot_cohort.tsv", root=root),
        excluded_samples=project_path("results", "tables", "level3_excluded_samples.tsv", root=root),
        qc_summary=project_path("results", "tables", "level3_candidate_selection_qc_summary.tsv", root=root),
        mutation_purity_figure=project_path("results", "figures", "level3_candidate_mutation_count_vs_purity.pdf", root=root),
        counts_by_project_figure=project_path("results", "figures", "level3_candidate_counts_by_project.pdf", root=root),
        score_by_project_figure=project_path("results", "figures", "level3_candidate_score_by_project.pdf", root=root),
        pilot_overview_figure=project_path("results", "figures", "level3_pilot_cohort_overview.pdf", root=root),
    )


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["NA", ""])


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def patient_from_sample(sample_barcode: str) -> str:
    parts = str(sample_barcode).split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else str(sample_barcode)[:12]


def mutation_ref_alt_summary(mutation_parquet: Path) -> pd.DataFrame:
    require_file(mutation_parquet, "standardized mutation parquet")
    columns = ["sample_barcode", "t_ref_count", "t_alt_count"]
    table = pq.read_table(mutation_parquet, columns=columns)
    mutations = table.to_pandas()
    mutations["t_ref_count"] = numeric(mutations["t_ref_count"])
    mutations["t_alt_count"] = numeric(mutations["t_alt_count"])
    mutations["has_ref_alt"] = mutations["t_ref_count"].notna() & mutations["t_alt_count"].notna()
    mutations["has_positive_depth"] = mutations["has_ref_alt"] & ((mutations["t_ref_count"] + mutations["t_alt_count"]) > 0)
    summary = (
        mutations.groupby("sample_barcode", as_index=False)
        .agg(
            n_mutation_records=("sample_barcode", "size"),
            n_mutations_with_ref_alt_counts=("has_ref_alt", "sum"),
            n_mutations_with_positive_depth=("has_positive_depth", "sum"),
        )
        .reset_index(drop=True)
    )
    summary["ref_alt_counts_available"] = summary["n_mutations_with_positive_depth"] > 0
    return summary


def sample_set_from_tsv(path: Path, sample_column: str = "sample_barcode") -> set[str]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    try:
        data = pd.read_csv(path, sep="\t", usecols=[sample_column], dtype=str)
    except (ValueError, pd.errors.EmptyDataError):
        return set()
    return set(data[sample_column].dropna().astype(str))


def sample_set_from_parquet(path: Path, sample_column: str = "sample_barcode") -> set[str]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    try:
        table = pq.read_table(path, columns=[sample_column])
    except Exception:
        return set()
    data = table.to_pandas()
    return set(data[sample_column].dropna().astype(str))


def load_optional_copy_number_sets(paths: Level3Paths) -> tuple[set[str], set[str]]:
    segment_samples = sample_set_from_tsv(paths.segments_tsv) | sample_set_from_parquet(paths.segments_parquet)
    arm_call_samples = sample_set_from_tsv(paths.arm_calls)
    return segment_samples, arm_call_samples


def attach_metadata(tmb: pd.DataFrame, projects: pd.DataFrame, group_map: pd.DataFrame) -> pd.DataFrame:
    project_cols = [c for c in ["project_id", "project_code", "disease_type", "primary_site"] if c in projects.columns]
    group_cols = [c for c in ["project_id", "project_code", "broad_group", "major_cancer_group"] if c in group_map.columns]
    out = tmb.merge(projects[project_cols].drop_duplicates(), on=["project_id", "project_code"], how="left")
    group_meta = group_map[group_cols].drop_duplicates()
    if "major_cancer_group" not in group_meta.columns and "broad_group" in group_meta.columns:
        group_meta = group_meta.assign(major_cancer_group=group_meta["broad_group"])
    keep = [c for c in ["project_id", "project_code", "major_cancer_group"] if c in group_meta.columns]
    out = out.merge(group_meta[keep].drop_duplicates(), on=["project_id", "project_code"], how="left")
    for col in ["disease_type", "primary_site", "major_cancer_group"]:
        if col not in out.columns:
            out[col] = pd.NA
    out["major_cancer_group"] = out["major_cancer_group"].fillna("Unannotated")
    return out


def min_mutation_threshold(row: pd.Series, config: dict = CONFIG) -> int:
    if str(row.get("major_cancer_group", "")) == "Hematolymphoid":
        return int(config["min_nonsynonymous_hematolymphoid"])
    return int(config["min_nonsynonymous_solid"])


def exclusion_reasons(row: pd.Series, config: dict = CONFIG) -> list[str]:
    reasons: list[str] = []
    threshold = min_mutation_threshold(row, config)
    if not bool(row.get("ref_alt_counts_available", False)):
        reasons.append("missing_ref_alt_counts")
    if not bool(row.get("purity_available", False)):
        reasons.append("missing_purity")
    if not bool(row.get("ploidy_available", False)):
        reasons.append("missing_ploidy")
    if pd.isna(row.get("nonsynonymous_count")) or float(row["nonsynonymous_count"]) < threshold:
        reasons.append(f"too_few_mutations_lt_{threshold}")
    if bool(row.get("purity_available", False)) and float(row["purity"]) < float(config["min_purity"]):
        reasons.append(f"low_purity_lt_{config['min_purity']}")
    if not bool(row.get("copy_number_available", False)):
        reasons.append("missing_copy_number_data")
    return reasons


def candidate_class(row: pd.Series, reasons: list[str]) -> str:
    hard_failures = [
        reason
        for reason in reasons
        if reason
        in {
            "missing_ref_alt_counts",
            "missing_purity",
            "missing_ploidy",
            "missing_copy_number_data",
        }
        or reason.startswith("too_few_mutations")
        or reason.startswith("low_purity")
    ]
    if hard_failures:
        if reasons == ["missing_copy_number_data"]:
            return "descriptive_only"
        return "excluded"
    if bool(row.get("has_segment_level_cn", False)):
        return "copy_number_aware_clonal_phylogeny_candidate"
    return "clonal_clustering_candidate_limited_cn"


def score_candidate(row: pd.Series, reasons: list[str], config: dict = CONFIG) -> float:
    # Transparent 0-100 score:
    # 40 points mutation burden, 25 purity, 25 availability of ref/alt + purity/ploidy/copy-number,
    # 10 segment-level CN bonus. Penalties reduce low-information samples but are not hidden filters.
    mutation_cap = float(config["mutation_score_cap"])
    nonsyn = 0.0 if pd.isna(row.get("nonsynonymous_count")) else max(0.0, float(row["nonsynonymous_count"]))
    mutation_component = 40.0 * min(1.0, math.log1p(nonsyn) / math.log1p(mutation_cap))
    purity = 0.0 if pd.isna(row.get("purity")) else max(0.0, min(1.0, float(row["purity"])))
    purity_component = 25.0 * purity
    availability_component = 0.0
    availability_component += 8.0 if bool(row.get("ref_alt_counts_available", False)) else 0.0
    availability_component += 5.0 if bool(row.get("purity_available", False)) else 0.0
    availability_component += 5.0 if bool(row.get("ploidy_available", False)) else 0.0
    availability_component += 7.0 if bool(row.get("copy_number_available", False)) else 0.0
    segment_component = 10.0 if bool(row.get("has_segment_level_cn", False)) else 0.0
    penalty = 0.0
    if not bool(row.get("has_segment_level_cn", False)):
        penalty += 8.0
    if any(reason.startswith("too_few_mutations") for reason in reasons):
        penalty += 25.0
    if any(reason.startswith("low_purity") for reason in reasons):
        penalty += 25.0
    if any(reason in {"missing_ref_alt_counts", "missing_purity", "missing_ploidy"} for reason in reasons):
        penalty += 35.0
    if "missing_copy_number_data" in reasons:
        penalty += 12.0
    return round(max(0.0, min(100.0, mutation_component + purity_component + availability_component + segment_component - penalty)), 3)


def build_candidate_table(
    tmb_by_sample: pd.DataFrame,
    purity_ploidy: pd.DataFrame,
    aneuploidy: pd.DataFrame,
    projects: pd.DataFrame,
    cancer_group_map: pd.DataFrame,
    ref_alt_summary: pd.DataFrame,
    segment_samples: set[str],
    arm_call_samples: set[str],
    config: dict = CONFIG,
) -> pd.DataFrame:
    tmb = tmb_by_sample.copy()
    tmb["nonsynonymous_count"] = numeric(tmb["nonsynonymous_count"])
    tmb["total_mutation_count"] = numeric(tmb["total_mutation_count"])
    tmb["patient_barcode"] = tmb.get("patient_barcode", tmb["sample_barcode"].map(patient_from_sample))
    candidates = attach_metadata(tmb, projects, cancer_group_map)

    purity_cols = [c for c in ["sample_barcode", "purity", "ploidy"] if c in purity_ploidy.columns]
    purity = purity_ploidy[purity_cols].copy()
    if "purity" in purity:
        purity["purity"] = numeric(purity["purity"])
    if "ploidy" in purity:
        purity["ploidy"] = numeric(purity["ploidy"])
    candidates = candidates.merge(purity.drop_duplicates("sample_barcode"), on="sample_barcode", how="left")

    aneuploidy_cols = [
        c
        for c in [
            "sample_barcode",
            "aneuploidy_score",
            "arm_gain_count",
            "arm_loss_count",
            "total_arm_alteration_count",
        ]
        if c in aneuploidy.columns
    ]
    cn = aneuploidy[aneuploidy_cols].copy()
    for col in set(aneuploidy_cols) - {"sample_barcode"}:
        cn[col] = numeric(cn[col])
    candidates = candidates.merge(cn.drop_duplicates("sample_barcode"), on="sample_barcode", how="left")

    candidates = candidates.merge(ref_alt_summary, on="sample_barcode", how="left")
    candidates["ref_alt_counts_available"] = candidates["ref_alt_counts_available"].fillna(False).astype(bool)
    candidates["purity_available"] = candidates["purity"].notna()
    candidates["ploidy_available"] = candidates["ploidy"].notna()
    candidates["has_segment_level_cn"] = candidates["sample_barcode"].isin(segment_samples)
    candidates["copy_number_available"] = (
        candidates["aneuploidy_score"].notna()
        | candidates["arm_gain_count"].notna()
        | candidates["arm_loss_count"].notna()
        | candidates["sample_barcode"].isin(arm_call_samples)
        | candidates["has_segment_level_cn"]
    )

    candidates["mutation_count_proxy_rank_within_project"] = (
        candidates.groupby("project_code")["nonsynonymous_count"].rank(method="first", ascending=False).astype("Int64")
    )
    candidates["mutation_count_proxy_rank_pan_cancer"] = candidates["nonsynonymous_count"].rank(method="first", ascending=False).astype("Int64")

    reasons = candidates.apply(lambda row: exclusion_reasons(row, config), axis=1)
    candidates["exclusion_reason"] = reasons.map(lambda values: ";".join(values) if values else "none")
    candidates["candidate_class"] = [
        candidate_class(row, row_reasons) for (_, row), row_reasons in zip(candidates.iterrows(), reasons, strict=False)
    ]
    candidates["candidate_score"] = [
        score_candidate(row, row_reasons, config) for (_, row), row_reasons in zip(candidates.iterrows(), reasons, strict=False)
    ]
    candidates["inclusion_status"] = candidates["candidate_class"].map(
        {
            "copy_number_aware_clonal_phylogeny_candidate": "eligible_not_selected",
            "clonal_clustering_candidate_limited_cn": "eligible_not_selected",
            "descriptive_only": "excluded",
            "excluded": "excluded",
        }
    )
    candidates["notes"] = candidates.apply(candidate_notes, axis=1)

    ordered_cols = [
        "sample_barcode",
        "patient_barcode",
        "project_id",
        "project_code",
        "disease_type",
        "primary_site",
        "major_cancer_group",
        "nonsynonymous_count",
        "total_mutation_count",
        "mutation_count_proxy_rank_within_project",
        "mutation_count_proxy_rank_pan_cancer",
        "purity",
        "ploidy",
        "aneuploidy_score",
        "arm_gain_count",
        "arm_loss_count",
        "total_arm_alteration_count",
        "ref_alt_counts_available",
        "purity_available",
        "ploidy_available",
        "copy_number_available",
        "has_segment_level_cn",
        "candidate_class",
        "candidate_score",
        "inclusion_status",
        "exclusion_reason",
        "notes",
    ]
    return candidates[ordered_cols].sort_values(["candidate_score", "nonsynonymous_count"], ascending=[False, False]).reset_index(drop=True)


def candidate_notes(row: pd.Series) -> str:
    notes = ["single_bulk_tcga_sample"]
    if row.get("candidate_class") == "clonal_clustering_candidate_limited_cn":
        notes.append("segment_level_cn_missing_suitable_for_clonal_clustering_prototype_not_full_cn_aware_phylogeny")
    if row.get("candidate_class") == "copy_number_aware_clonal_phylogeny_candidate":
        notes.append("segment_level_cn_available")
    if row.get("project_code") in {"UCEC", "COAD", "READ"} and row.get("nonsynonymous_count", 0) >= 500:
        notes.append("very_high_mutation_count_possible_msi_pole_like_proxy")
    if row.get("project_code") in {"SKCM", "LUAD", "LUSC", "BLCA", "HNSC"} and row.get("nonsynonymous_count", 0) >= 250:
        notes.append("high_mutation_count_priority_project")
    return ";".join(notes)


def select_pilot_cohort(candidates: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    eligible = candidates[candidates["inclusion_status"] == "eligible_not_selected"].copy()
    if eligible.empty:
        return eligible
    eligible["_priority"] = eligible["project_code"].map(
        {project: i for i, project in enumerate(config["pilot_project_priority"])}
    )
    eligible["_priority"] = eligible["_priority"].fillna(len(config["pilot_project_priority"]) + 1)
    eligible = eligible.sort_values(["_priority", "candidate_score", "nonsynonymous_count"], ascending=[True, False, False])

    selected_parts: list[pd.DataFrame] = []
    selected_indices: set[int] = set()
    per_project_cap = int(config["pilot_per_project_cap"])
    total_max = int(config["pilot_total_max"])
    total_min = int(config["pilot_total_min"])

    for project in config["pilot_project_priority"]:
        project_candidates = eligible[eligible["project_code"] == project].head(per_project_cap)
        if project_candidates.empty:
            continue
        remaining = total_max - sum(len(part) for part in selected_parts)
        if remaining <= 0:
            break
        project_candidates = project_candidates.head(remaining)
        selected_parts.append(project_candidates)
        selected_indices.update(project_candidates.index)

    selected = pd.concat(selected_parts, ignore_index=False) if selected_parts else eligible.iloc[0:0].copy()
    if len(selected) < total_min:
        counts = selected["project_code"].value_counts().to_dict()
        filler_rows = []
        for idx, row in eligible.drop(index=list(selected_indices), errors="ignore").sort_values(
            ["candidate_score", "nonsynonymous_count"], ascending=[False, False]
        ).iterrows():
            if len(selected) + len(filler_rows) >= total_max:
                break
            project = row["project_code"]
            if counts.get(project, 0) >= per_project_cap:
                continue
            filler_rows.append(row)
            counts[project] = counts.get(project, 0) + 1
            selected_indices.add(idx)
            if len(selected) + len(filler_rows) >= total_min:
                break
        if filler_rows:
            selected = pd.concat([selected, pd.DataFrame(filler_rows)], ignore_index=False)

    selected = selected.drop(columns=["_priority"], errors="ignore")
    selected = selected.sort_values(["project_code", "candidate_score"], ascending=[True, False]).reset_index(drop=True)
    return selected


def apply_pilot_status(candidates: pd.DataFrame, pilot: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    selected = set(pilot["sample_barcode"]) if not pilot.empty else set()
    out.loc[out["sample_barcode"].isin(selected), "inclusion_status"] = "selected_for_pilot"
    return out


def summarize_by_project(candidates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (project_id, project_code), group in candidates.groupby(["project_id", "project_code"], dropna=False):
        eligible = group[group["candidate_class"].isin(["copy_number_aware_clonal_phylogeny_candidate", "clonal_clustering_candidate_limited_cn"])]
        top = group.sort_values(["candidate_score", "nonsynonymous_count"], ascending=[False, False]).head(1)
        rows.append(
            {
                "project_id": project_id,
                "project_code": project_code,
                "n_mutation_samples": len(group),
                "n_samples_with_ref_alt_counts": int(group["ref_alt_counts_available"].sum()),
                "n_samples_with_purity_ploidy": int((group["purity_available"] & group["ploidy_available"]).sum()),
                "n_samples_with_copy_number": int(group["copy_number_available"].sum()),
                "n_samples_passing_minimum_filters": len(eligible),
                "median_mutation_count_proxy": group["nonsynonymous_count"].median(),
                "median_purity": group["purity"].median(),
                "median_ploidy": group["ploidy"].median(),
                "median_aneuploidy_score": group["aneuploidy_score"].median(),
                "top_candidate_sample": top["sample_barcode"].iloc[0] if not top.empty else pd.NA,
                "top_candidate_score": top["candidate_score"].iloc[0] if not top.empty else pd.NA,
                "notes": project_notes(group),
            }
        )
    return pd.DataFrame(rows).sort_values(["n_samples_passing_minimum_filters", "top_candidate_score"], ascending=[False, False])


def project_notes(group: pd.DataFrame) -> str:
    notes = []
    if group["has_segment_level_cn"].sum() == 0:
        notes.append("no_segment_level_cn_for_project")
    if group["candidate_class"].isin(["copy_number_aware_clonal_phylogeny_candidate", "clonal_clustering_candidate_limited_cn"]).sum() == 0:
        notes.append("no_samples_passing_minimum_filters")
    return ";".join(notes) if notes else "candidate_samples_available"


def exclusion_table(candidates: pd.DataFrame) -> pd.DataFrame:
    excluded = candidates[candidates["inclusion_status"] == "excluded"].copy()
    return excluded.sort_values(["project_code", "exclusion_reason", "sample_barcode"]).reset_index(drop=True)


def qc_summary(candidates: pd.DataFrame, pilot: pd.DataFrame) -> pd.DataFrame:
    rows = [
        ("global", "total_mutation_layer_samples", len(candidates)),
        ("global", "samples_with_ref_alt_counts", int(candidates["ref_alt_counts_available"].sum())),
        ("global", "samples_with_purity", int(candidates["purity_available"].sum())),
        ("global", "samples_with_ploidy", int(candidates["ploidy_available"].sum())),
        ("global", "samples_with_copy_number_or_aneuploidy", int(candidates["copy_number_available"].sum())),
        (
            "global",
            "samples_passing_minimum_filters",
            int(candidates["candidate_class"].isin(["copy_number_aware_clonal_phylogeny_candidate", "clonal_clustering_candidate_limited_cn"]).sum()),
        ),
        ("global", "samples_selected_for_pilot_cohort", len(pilot)),
        ("global", "projects_represented_in_pilot_cohort", ";".join(sorted(pilot["project_code"].unique())) if not pilot.empty else "none"),
        ("global", "segment_level_cn_candidates", int(candidates["has_segment_level_cn"].sum())),
        (
            "global",
            "limited_cn_clonal_clustering_candidates",
            int((candidates["candidate_class"] == "clonal_clustering_candidate_limited_cn").sum()),
        ),
        (
            "global",
            "limitations",
            "TCGA is mostly single bulk tumor; candidate selection supports clonal clustering prototype, while robust branching phylogeny needs serial, relapse, metastasis, or multi-region data; segment-level CN is absent for current candidates"
            if candidates["has_segment_level_cn"].sum() == 0
            else "TCGA is mostly single bulk tumor; branching phylogeny remains limited without serial or multi-region data",
        ),
    ]
    reason_counts: dict[str, int] = {}
    for value in candidates.loc[candidates["inclusion_status"] == "excluded", "exclusion_reason"].dropna():
        for reason in str(value).split(";"):
            if reason and reason != "none":
                reason_counts[reason] = reason_counts.get(reason, 0) + 1
    for reason, count in sorted(reason_counts.items()):
        rows.append(("excluded_by_reason", reason, count))
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def plot_mutation_vs_purity(candidates: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(10, 7))
    plot_data = candidates[candidates["purity"].notna()].copy()
    colors = plot_data["inclusion_status"].map(
        {"selected_for_pilot": "#1b9e77", "eligible_not_selected": "#377eb8", "excluded": "#bdbdbd"}
    ).fillna("#bdbdbd")
    ax.scatter(plot_data["purity"], plot_data["nonsynonymous_count"], c=colors, s=14, alpha=0.7, linewidths=0)
    ax.axvline(CONFIG["min_purity"], color="black", linestyle="--", linewidth=0.8)
    ax.axhline(CONFIG["min_nonsynonymous_solid"], color="black", linestyle=":", linewidth=0.8)
    ax.set_xlabel("Purity")
    ax.set_ylabel("Nonsynonymous mutation-count proxy")
    ax.set_title("Level 3 candidate mutation count vs purity")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_counts_by_project(summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = summary.sort_values("n_samples_passing_minimum_filters", ascending=False).head(33)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.bar(data["project_code"], data["n_samples_passing_minimum_filters"], color="#4daf4a")
    ax.set_ylabel("Samples passing minimum filters")
    ax.set_title("Level 3 candidate counts by TCGA project")
    ax.tick_params(axis="x", labelrotation=90)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_score_by_project(candidates: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = candidates[candidates["candidate_class"].isin(["copy_number_aware_clonal_phylogeny_candidate", "clonal_clustering_candidate_limited_cn"])]
    order = data.groupby("project_code")["candidate_score"].median().sort_values(ascending=False).index
    grouped = [data.loc[data["project_code"] == project, "candidate_score"].values for project in order]
    fig, ax = plt.subplots(figsize=(13, 7))
    try:
        ax.boxplot(grouped, tick_labels=order, showfliers=False)
    except TypeError:
        ax.boxplot(grouped, labels=order, showfliers=False)
    ax.set_ylabel("Candidate score")
    ax.set_title("Level 3 candidate score distribution by project")
    ax.tick_params(axis="x", labelrotation=90)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_pilot_overview(pilot: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(10, 6))
    if pilot.empty:
        ax.text(0.5, 0.5, "No pilot samples selected", ha="center", va="center")
        ax.axis("off")
    else:
        counts = pilot["project_code"].value_counts().sort_values(ascending=False)
        ax.bar(counts.index, counts.values, color="#984ea3")
        ax.set_ylabel("Pilot samples")
        ax.set_title("Level 3 pilot cohort composition")
        ax.tick_params(axis="x", labelrotation=45)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def run_candidate_selection(paths: Level3Paths, config: dict = CONFIG) -> dict[str, pd.DataFrame]:
    for path, description in [
        (paths.mutation_parquet, "standardized mutation parquet"),
        (paths.tmb_by_sample, "sample-level mutation-count proxy table"),
        (paths.purity_ploidy_by_sample, "sample-level purity/ploidy table"),
        (paths.aneuploidy_by_sample, "sample-level aneuploidy table"),
        (paths.projects, "TCGA project metadata"),
        (paths.cancer_group_map, "cancer group map"),
    ]:
        require_file(path, description)

    tmb = read_tsv(paths.tmb_by_sample)
    purity = read_tsv(paths.purity_ploidy_by_sample)
    aneuploidy = read_tsv(paths.aneuploidy_by_sample)
    projects = read_tsv(paths.projects)
    group_map = pd.read_csv(paths.cancer_group_map, dtype=str)
    ref_alt = mutation_ref_alt_summary(paths.mutation_parquet)
    segment_samples, arm_call_samples = load_optional_copy_number_sets(paths)

    candidates = build_candidate_table(tmb, purity, aneuploidy, projects, group_map, ref_alt, segment_samples, arm_call_samples, config)
    pilot = select_pilot_cohort(candidates, config)
    candidates = apply_pilot_status(candidates, pilot)
    pilot = candidates[candidates["inclusion_status"] == "selected_for_pilot"].copy().sort_values(
        ["project_code", "candidate_score"], ascending=[True, False]
    )
    summary = summarize_by_project(candidates)
    excluded = exclusion_table(candidates)
    qc = qc_summary(candidates, pilot)

    write_tsv(candidates, paths.candidate_samples)
    write_tsv(summary, paths.summary_by_project)
    write_tsv(pilot, paths.pilot_cohort)
    write_tsv(excluded, paths.excluded_samples)
    write_tsv(qc, paths.qc_summary)
    plot_mutation_vs_purity(candidates, paths.mutation_purity_figure)
    plot_counts_by_project(summary, paths.counts_by_project_figure)
    plot_score_by_project(candidates, paths.score_by_project_figure)
    plot_pilot_overview(pilot, paths.pilot_overview_figure)
    return {
        "candidates": candidates,
        "summary": summary,
        "pilot": pilot,
        "excluded": excluded,
        "qc": qc,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Select TCGA samples feasible for Level 3 clonal-analysis pilot work.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    paths = default_paths(root)
    logging.warning(
        "Level 3 candidate selection is not full clonal phylogeny. TCGA is mostly single bulk tumor; robust branching order usually needs serial, relapse, metastasis, or multi-region samples."
    )
    outputs = run_candidate_selection(paths, CONFIG)
    logging.info(
        "Level 3 candidate selection complete: %d samples evaluated, %d passing minimum filters, %d pilot samples",
        len(outputs["candidates"]),
        int(outputs["candidates"]["candidate_class"].isin(["copy_number_aware_clonal_phylogeny_candidate", "clonal_clustering_candidate_limited_cn"]).sum()),
        len(outputs["pilot"]),
    )


if __name__ == "__main__":
    main()

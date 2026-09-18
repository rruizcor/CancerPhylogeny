#!/usr/bin/env python

from __future__ import annotations

import argparse
import gzip
import logging
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcga_cancer_phylogeny_matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


CONFIG = {
    "allow_patient_level_segment_matching": False,
    "min_total_depth": 20,
    "candidate_nonsynonymous_only": True,
    "min_candidate_mutations_with_local_cn": 50,
    "min_annotation_rate": 0.80,
    "min_median_depth": 20,
    "gain_like_threshold": 0.2,
    "loss_like_threshold": -0.2,
    "scatter_max_points": 50000,
}

MUTATION_COLUMNS = [
    "Tumor_Sample_Barcode",
    "project_id",
    "project_code",
    "sample_barcode",
    "patient_barcode",
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "Variant_Type",
    "t_ref_count",
    "t_alt_count",
    "HGVSp_Short",
    "is_nonsynonymous",
]

ANNOTATED_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "chromosome",
    "position",
    "start_position",
    "end_position",
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
    "local_cn_segment_chromosome",
    "local_cn_segment_start",
    "local_cn_segment_end",
    "local_cn_num_probes",
    "local_segment_mean",
    "segment_mean_type",
    "local_cn_status",
    "local_cn_match_status",
    "local_cn_match_notes",
    "is_nonsynonymous",
    "is_clonal_input_candidate",
]

SUMMARY_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "n_total_mutations",
    "n_nonsynonymous_mutations",
    "n_mutations_with_ref_alt_counts",
    "n_mutations_with_local_cn",
    "pct_mutations_with_local_cn",
    "n_mutations_without_local_cn",
    "median_depth",
    "median_vaf",
    "median_local_segment_mean",
    "purity",
    "ploidy",
    "eligible_for_copy_number_aware_clonal_input",
    "exclusion_reason_or_warning",
]

PYCLONE_COLUMNS = [
    "sample_id",
    "mutation_id",
    "ref_counts",
    "var_counts",
    "normal_cn",
    "major_cn",
    "minor_cn",
    "local_segment_mean",
    "purity",
    "notes",
]


@dataclass
class MutationLocalCnPaths:
    mutations: Path
    segments: Path
    segment_summary: Path
    pilot_with_segments: Path
    candidates_with_segments: Path
    purity_ploidy: Path
    projects: Path
    best_files: Path
    download_log: Path
    segment_qc: Path
    annotated_mutations: Path
    summary_by_sample: Path
    pyclone_preview: Path
    qc_summary: Path
    coverage_figure: Path
    vaf_segment_figure: Path
    depth_figure: Path


def default_paths(root: Path) -> MutationLocalCnPaths:
    return MutationLocalCnPaths(
        mutations=project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root=root),
        segments=project_path("data", "processed", "copy_number", "segments_by_sample.tsv.gz", root=root),
        segment_summary=project_path("data", "processed", "copy_number", "segment_level_cn_summary_by_sample.tsv", root=root),
        pilot_with_segments=project_path("results", "tables", "level3_pilot_cohort_with_segments.tsv", root=root),
        candidates_with_segments=project_path("results", "tables", "level3_candidate_samples_with_segments.tsv", root=root),
        purity_ploidy=project_path("data", "processed", "features", "purity_ploidy_by_sample.tsv", root=root),
        projects=project_path("data", "interim", "tcga_projects.tsv", root=root),
        best_files=project_path("results", "tables", "gdc_segment_cn_best_file_per_pilot_sample.tsv", root=root),
        download_log=project_path("results", "tables", "gdc_segment_cn_download_log.tsv", root=root),
        segment_qc=project_path("results", "tables", "level3_segment_cn_qc_summary.tsv", root=root),
        annotated_mutations=project_path("data", "processed", "clonal", "level3_mutations_with_local_cn.tsv.gz", root=root),
        summary_by_sample=project_path("results", "tables", "level3_mutation_local_cn_summary_by_sample.tsv", root=root),
        pyclone_preview=project_path("data", "processed", "clonal", "level3_pyclone_input_preview.tsv", root=root),
        qc_summary=project_path("results", "tables", "level3_mutation_local_cn_annotation_qc_summary.tsv", root=root),
        coverage_figure=project_path(
            "results", "figures", "level3_mutation_cn_annotation_coverage_by_sample.pdf", root=root
        ),
        vaf_segment_figure=project_path(
            "results", "figures", "level3_mutation_vaf_vs_local_segment_mean.pdf", root=root
        ),
        depth_figure=project_path("results", "figures", "level3_mutation_depth_distribution.pdf", root=root),
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


def normalize_patient_barcode(barcode: str) -> str:
    text = str(barcode).strip()
    parts = text.split("-")
    if len(parts) >= 3 and parts[0] == "TCGA":
        return "-".join(parts[:3])
    return text[:12]


def normalize_sample_barcode(barcode: str) -> str:
    text = str(barcode).strip()
    parts = text.split("-")
    if len(parts) >= 4 and parts[0] == "TCGA":
        return "-".join(parts[:4])[:16]
    return text[:16]


def normalize_chromosome(value: Any) -> str:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na"}:
        return "NA"
    text = text.replace("chr", "").replace("CHR", "").strip()
    text = text.upper()
    if text in {"23"}:
        return "X"
    if text in {"24"}:
        return "Y"
    if text in {"25", "M", "MT", "MITO", "MITOCHONDRIA"}:
        return "MT"
    try:
        number = int(float(text))
        return str(number)
    except ValueError:
        return text


def is_autosome(chromosome: str) -> bool:
    return str(chromosome) in {str(i) for i in range(1, 23)}


def local_cn_status(segment_mean: Any, config: dict = CONFIG) -> str:
    value = pd.to_numeric(pd.Series([segment_mean]), errors="coerce").iloc[0]
    if pd.isna(value):
        return "unknown"
    if value >= float(config["gain_like_threshold"]):
        return "gain_like"
    if value <= float(config["loss_like_threshold"]):
        return "loss_like"
    return "neutral_like"


def load_pilot_mutations(mutation_parquet: Path, pilot_samples: set[str]) -> pd.DataFrame:
    require_file(mutation_parquet, "MC3 somatic mutation parquet")
    table = pq.read_table(mutation_parquet, columns=MUTATION_COLUMNS)
    mutations = table.to_pandas()
    mutations["sample_barcode"] = mutations["sample_barcode"].map(normalize_sample_barcode)
    mutations = mutations[mutations["sample_barcode"].isin(pilot_samples)].copy()
    mutations["patient_barcode"] = mutations["sample_barcode"].map(normalize_patient_barcode)
    mutations["chromosome"] = mutations["Chromosome"].map(normalize_chromosome)
    mutations["position"] = numeric(mutations["Start_Position"])
    mutations["start_position"] = numeric(mutations["Start_Position"])
    mutations["end_position"] = numeric(mutations["End_Position"])
    mutations["t_ref_count"] = numeric(mutations["t_ref_count"])
    mutations["t_alt_count"] = numeric(mutations["t_alt_count"])
    mutations["total_depth"] = mutations["t_ref_count"] + mutations["t_alt_count"]
    mutations["observed_vaf"] = np.where(
        mutations["total_depth"] > 0,
        mutations["t_alt_count"] / mutations["total_depth"],
        np.nan,
    )
    if "HGVSp_Short" not in mutations.columns:
        mutations["HGVSp_Short"] = pd.NA
    if "is_nonsynonymous" not in mutations.columns:
        mutations["is_nonsynonymous"] = mutations["Variant_Classification"].astype(str) != "Silent"
    return mutations


def load_segments(segments_path: Path, pilot_samples: set[str], config: dict = CONFIG) -> pd.DataFrame:
    require_file(segments_path, "processed segment-level CN table")
    segments = pd.read_csv(segments_path, sep="\t")
    segments["sample_barcode"] = segments["sample_barcode"].map(normalize_sample_barcode)
    if not bool(config.get("allow_patient_level_segment_matching", False)) and "match_level" in segments.columns:
        segments = segments[segments["match_level"].astype(str) == "sample"].copy()
    segments = segments[segments["sample_barcode"].isin(pilot_samples)].copy()
    segments["chromosome_norm"] = segments["chromosome"].map(normalize_chromosome)
    segments["start"] = numeric(segments["start"])
    segments["end"] = numeric(segments["end"])
    segments["num_probes"] = numeric(segments["num_probes"]) if "num_probes" in segments.columns else pd.NA
    segments["segment_mean"] = numeric(segments["segment_mean"])
    segments = segments[
        segments["chromosome_norm"].ne("MT")
        & segments["start"].notna()
        & segments["end"].notna()
        & segments["segment_mean"].notna()
    ].copy()
    return segments.sort_values(["sample_barcode", "chromosome_norm", "start", "end"]).reset_index(drop=True)


def annotate_mutations_with_segments(mutations: pd.DataFrame, segments: pd.DataFrame) -> pd.DataFrame:
    annotated = mutations.copy()
    for column in [
        "local_cn_segment_chromosome",
        "local_cn_segment_start",
        "local_cn_segment_end",
        "local_cn_num_probes",
        "local_segment_mean",
    ]:
        annotated[column] = pd.NA
    annotated["segment_mean_type"] = "log2_ratio_or_segment_mean_unverified"
    annotated["local_cn_status"] = "unknown"
    annotated["local_cn_match_status"] = "unmatched"
    annotated["local_cn_match_notes"] = "no_overlapping_sample_level_segment"

    if annotated.empty or segments.empty:
        return annotated

    for (sample, chrom), mut_idx in annotated.groupby(["sample_barcode", "chromosome"], dropna=False).groups.items():
        if chrom in {"NA", "MT"}:
            annotated.loc[mut_idx, "local_cn_match_notes"] = "chromosome_excluded_or_unrecognized"
            continue
        seg = segments[(segments["sample_barcode"] == sample) & (segments["chromosome_norm"] == chrom)]
        if seg.empty:
            continue
        seg = seg.sort_values(["start", "end"]).reset_index(drop=True)
        starts = seg["start"].to_numpy(dtype=float)
        ends = seg["end"].to_numpy(dtype=float)
        positions = annotated.loc[mut_idx, "position"].to_numpy(dtype=float)
        candidate_idx = np.searchsorted(starts, positions, side="right") - 1
        valid = (candidate_idx >= 0) & (positions <= ends[np.clip(candidate_idx, 0, len(ends) - 1)])
        idx_array = np.array(list(mut_idx))
        matched_mut_idx = idx_array[valid]
        matched_seg_idx = candidate_idx[valid]
        if len(matched_mut_idx) == 0:
            continue
        matched_segments = seg.iloc[matched_seg_idx].reset_index(drop=True)
        annotated.loc[matched_mut_idx, "local_cn_segment_chromosome"] = matched_segments["chromosome_norm"].to_numpy()
        annotated.loc[matched_mut_idx, "local_cn_segment_start"] = matched_segments["start"].to_numpy()
        annotated.loc[matched_mut_idx, "local_cn_segment_end"] = matched_segments["end"].to_numpy()
        annotated.loc[matched_mut_idx, "local_cn_num_probes"] = matched_segments["num_probes"].to_numpy()
        annotated.loc[matched_mut_idx, "local_segment_mean"] = matched_segments["segment_mean"].to_numpy()
        annotated.loc[matched_mut_idx, "local_cn_match_status"] = "matched"
        annotated.loc[matched_mut_idx, "local_cn_match_notes"] = "sample_level_segment_overlap"

    annotated["local_cn_status"] = annotated["local_segment_mean"].map(local_cn_status)
    return annotated


def attach_purity_ploidy_and_metadata(annotated: pd.DataFrame, purity_ploidy: pd.DataFrame, pilot: pd.DataFrame) -> pd.DataFrame:
    purity_cols = ["sample_barcode", "purity", "ploidy"]
    purity = purity_ploidy[purity_cols].drop_duplicates().copy()
    purity["sample_barcode"] = purity["sample_barcode"].map(normalize_sample_barcode)
    purity["purity"] = numeric(purity["purity"])
    purity["ploidy"] = numeric(purity["ploidy"])
    out = annotated.drop(columns=["purity", "ploidy"], errors="ignore").merge(purity, on="sample_barcode", how="left")

    pilot_meta_cols = ["sample_barcode", "patient_barcode", "project_id", "project_code"]
    pilot_meta = pilot[pilot_meta_cols].drop_duplicates().copy()
    pilot_meta["sample_barcode"] = pilot_meta["sample_barcode"].map(normalize_sample_barcode)
    out = out.drop(columns=["patient_barcode", "project_id", "project_code"], errors="ignore").merge(
        pilot_meta, on="sample_barcode", how="left"
    )
    return out


def flag_clonal_input_candidates(annotated: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    out = annotated.copy()
    has_counts = out["t_ref_count"].notna() & out["t_alt_count"].notna()
    candidate = (
        has_counts
        & (out["total_depth"] >= int(config["min_total_depth"]))
        & (out["observed_vaf"] > 0)
        & (out["local_cn_match_status"] == "matched")
    )
    if bool(config.get("candidate_nonsynonymous_only", True)):
        candidate &= out["is_nonsynonymous"].astype(bool)
    out["is_clonal_input_candidate"] = candidate
    return out


def build_sample_summary(
    annotated: pd.DataFrame,
    pilot: pd.DataFrame,
    segment_summary: pd.DataFrame,
    config: dict = CONFIG,
) -> pd.DataFrame:
    if annotated.empty:
        base = pilot[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates().copy()
        for column in SUMMARY_COLUMNS:
            if column not in base:
                base[column] = pd.NA
        base["eligible_for_copy_number_aware_clonal_input"] = False
        base["exclusion_reason_or_warning"] = "no_mutations_for_pilot_sample"
        return base[SUMMARY_COLUMNS]

    data = annotated.copy()
    data["has_ref_alt_counts"] = data["t_ref_count"].notna() & data["t_alt_count"].notna()
    data["has_local_cn"] = data["local_cn_match_status"] == "matched"
    grouped = data.groupby("sample_barcode", dropna=False)
    summary = grouped.agg(
        patient_barcode=("patient_barcode", "first"),
        project_id=("project_id", "first"),
        project_code=("project_code", "first"),
        n_total_mutations=("sample_barcode", "size"),
        n_nonsynonymous_mutations=("is_nonsynonymous", "sum"),
        n_mutations_with_ref_alt_counts=("has_ref_alt_counts", "sum"),
        n_mutations_with_local_cn=("has_local_cn", "sum"),
        median_depth=("total_depth", "median"),
        median_vaf=("observed_vaf", "median"),
        median_local_segment_mean=("local_segment_mean", lambda x: x.dropna().median() if x.notna().any() else np.nan),
        purity=("purity", "first"),
        ploidy=("ploidy", "first"),
        n_candidate_mutations=("is_clonal_input_candidate", "sum"),
    ).reset_index()
    summary["n_mutations_without_local_cn"] = summary["n_total_mutations"] - summary["n_mutations_with_local_cn"]
    summary["pct_mutations_with_local_cn"] = np.where(
        summary["n_total_mutations"] > 0,
        100 * summary["n_mutations_with_local_cn"] / summary["n_total_mutations"],
        0,
    )

    segment_samples = set(segment_summary.loc[as_bool(segment_summary["has_segment_level_cn"]), "sample_barcode"].astype(str))
    rows: list[str] = []
    eligible: list[bool] = []
    for row in summary.itertuples(index=False):
        reasons = []
        if pd.isna(row.purity):
            reasons.append("missing_purity")
        if pd.isna(row.ploidy):
            reasons.append("missing_ploidy")
        if row.sample_barcode not in segment_samples:
            reasons.append("missing_segment_level_cn")
        if int(row.n_candidate_mutations) < int(config["min_candidate_mutations_with_local_cn"]):
            reasons.append(f"too_few_candidate_mutations_lt_{config['min_candidate_mutations_with_local_cn']}")
        if pd.isna(row.median_depth) or float(row.median_depth) < float(config["min_median_depth"]):
            reasons.append(f"median_depth_lt_{config['min_median_depth']}")
        annotation_rate = float(row.n_mutations_with_local_cn) / float(row.n_total_mutations) if row.n_total_mutations else 0
        if annotation_rate < float(config["min_annotation_rate"]):
            reasons.append(f"local_cn_annotation_rate_lt_{config['min_annotation_rate']}")
        eligible.append(not reasons)
        rows.append("eligible" if not reasons else ";".join(reasons))
    summary["eligible_for_copy_number_aware_clonal_input"] = eligible
    summary["exclusion_reason_or_warning"] = rows

    pilot_meta = pilot[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates()
    summary = pilot_meta.merge(summary, on=["sample_barcode", "patient_barcode", "project_id", "project_code"], how="left")
    for column in [
        "n_total_mutations",
        "n_nonsynonymous_mutations",
        "n_mutations_with_ref_alt_counts",
        "n_mutations_with_local_cn",
        "n_mutations_without_local_cn",
    ]:
        summary[column] = summary[column].fillna(0).astype(int)
    summary["pct_mutations_with_local_cn"] = summary["pct_mutations_with_local_cn"].fillna(0)
    summary["eligible_for_copy_number_aware_clonal_input"] = summary[
        "eligible_for_copy_number_aware_clonal_input"
    ].fillna(False).astype(bool)
    summary["exclusion_reason_or_warning"] = summary["exclusion_reason_or_warning"].fillna("no_mutations_for_pilot_sample")
    return summary[SUMMARY_COLUMNS]


def build_pyclone_preview(annotated: pd.DataFrame) -> pd.DataFrame:
    if annotated.empty:
        return pd.DataFrame(columns=PYCLONE_COLUMNS)
    data = annotated[annotated["is_clonal_input_candidate"]].copy()
    if data.empty:
        return pd.DataFrame(columns=PYCLONE_COLUMNS)
    data["mutation_id"] = (
        data["sample_barcode"].astype(str)
        + ":"
        + data["chromosome"].astype(str)
        + ":"
        + data["position"].astype("Int64").astype(str)
        + ":"
        + data["Reference_Allele"].astype(str)
        + ">"
        + data["Tumor_Seq_Allele2"].astype(str)
    )
    preview = pd.DataFrame(
        {
            "sample_id": data["sample_barcode"],
            "mutation_id": data["mutation_id"],
            "ref_counts": data["t_ref_count"].astype("Int64"),
            "var_counts": data["t_alt_count"].astype("Int64"),
            "normal_cn": np.where(data["chromosome"].map(is_autosome), 2, pd.NA),
            "major_cn": pd.NA,
            "minor_cn": pd.NA,
            "local_segment_mean": data["local_segment_mean"],
            "purity": data["purity"],
            "notes": np.where(
                data["chromosome"].map(is_autosome),
                "preview_only_normal_cn_2_placeholder_for_autosome;major_minor_cn_unavailable;segment_mean_not_allele_specific",
                "preview_only_normal_cn_missing_for_non_autosome;major_minor_cn_unavailable;segment_mean_not_allele_specific",
            ),
        }
    )
    return preview[PYCLONE_COLUMNS]


def qc_summary(
    annotated: pd.DataFrame,
    sample_summary: pd.DataFrame,
    pilot: pd.DataFrame,
    segments: pd.DataFrame,
    config: dict = CONFIG,
) -> pd.DataFrame:
    total = len(annotated)
    annotated_n = int((annotated["local_cn_match_status"] == "matched").sum()) if not annotated.empty else 0
    pilot_with_segments = int(as_bool(pilot["has_segment_level_cn"]).sum()) if "has_segment_level_cn" in pilot else 0
    pilot_with_mutations = int((sample_summary["n_total_mutations"] > 0).sum()) if not sample_summary.empty else 0
    pilot_with_annotation = int((sample_summary["n_mutations_with_local_cn"] > 0).sum()) if not sample_summary.empty else 0
    eligible = int(sample_summary["eligible_for_copy_number_aware_clonal_input"].astype(bool).sum())
    excluded_reasons = (
        sample_summary.loc[~sample_summary["eligible_for_copy_number_aware_clonal_input"].astype(bool), "exclusion_reason_or_warning"]
        .value_counts()
        .to_dict()
    )
    mutation_chromosomes = set(annotated["chromosome"].astype(str)) if not annotated.empty else set()
    segment_chromosomes = set(segments["chromosome_norm"].astype(str)) if not segments.empty else set()
    chromosome_issues = sorted((mutation_chromosomes - segment_chromosomes) - {"MT", "NA"})
    unmatched = total - annotated_n
    rows = [
        ("global", "pilot_samples_evaluated", len(pilot)),
        ("global", "pilot_samples_with_segment_cn", pilot_with_segments),
        ("global", "pilot_samples_with_mutations", pilot_with_mutations),
        ("global", "pilot_samples_with_mutation_to_segment_annotation", pilot_with_annotation),
        ("global", "total_mutations_evaluated", total),
        ("global", "total_mutations_annotated_with_local_cn", annotated_n),
        ("global", "percentage_mutations_annotated", round(100 * annotated_n / total, 3) if total else 0),
        ("global", "samples_eligible_for_downstream_clonal_input", eligible),
        ("global", "samples_excluded_from_downstream_clonal_input", len(sample_summary) - eligible),
        ("global", "chromosome_naming_issues", ";".join(chromosome_issues) if chromosome_issues else "none"),
        (
            "global",
            "coordinate_mismatch_issues",
            f"{unmatched}_mutations_without_overlapping_sample_level_segment" if unmatched else "none",
        ),
        (
            "global",
            "non_allele_specific_cn_warning",
            "segment_mean is not allele-specific integer copy number; PyClone-style preview sets major_cn/minor_cn to NA and uses normal_cn=2 only as an autosomal placeholder.",
        ),
        ("global", "patient_level_segment_matching_allowed", bool(config.get("allow_patient_level_segment_matching", False))),
    ]
    for reason, count in excluded_reasons.items():
        rows.append(("sample_exclusion", str(reason), int(count)))
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def plot_coverage(summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = summary.sort_values("pct_mutations_with_local_cn", ascending=False)
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.bar(data["sample_barcode"], data["pct_mutations_with_local_cn"], color="#4c78a8")
    ax.set_ylabel("Mutations with local CN (%)")
    ax.set_title("Level 3 mutation-to-segment CN annotation coverage")
    ax.tick_params(axis="x", rotation=90, labelsize=6)
    ax.set_ylim(0, 105)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_vaf_vs_segment_mean(annotated: pd.DataFrame, path: Path, config: dict = CONFIG) -> None:
    ensure_dir(path.parent)
    data = annotated[
        (annotated["local_cn_match_status"] == "matched")
        & annotated["observed_vaf"].notna()
        & annotated["local_segment_mean"].notna()
    ].copy()
    if len(data) > int(config["scatter_max_points"]):
        data = data.sample(int(config["scatter_max_points"]), random_state=1)
    fig, ax = plt.subplots(figsize=(9, 6))
    if data.empty:
        ax.text(0.5, 0.5, "No mutation/local-CN overlaps available", ha="center", va="center")
        ax.axis("off")
    else:
        ax.scatter(data["local_segment_mean"], data["observed_vaf"], s=8, alpha=0.35, color="#4c78a8")
        ax.set_xlabel("Local segment mean")
        ax.set_ylabel("Observed VAF")
        ax.set_title("Pilot mutation VAF vs local segment mean")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_depth_distribution(annotated: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = annotated["total_depth"].dropna()
    fig, ax = plt.subplots(figsize=(9, 6))
    if data.empty:
        ax.text(0.5, 0.5, "No depth values available", ha="center", va="center")
        ax.axis("off")
    else:
        ax.hist(data.clip(upper=data.quantile(0.99)), bins=50, color="#54a24b")
        ax.set_xlabel("Total depth (clipped at 99th percentile)")
        ax.set_ylabel("Mutations")
        ax.set_title("Pilot mutation depth distribution")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def run_mutation_local_cn_annotation(paths: MutationLocalCnPaths, config: dict = CONFIG) -> dict[str, pd.DataFrame]:
    for path, description in [
        (paths.mutations, "MC3 somatic mutation parquet"),
        (paths.segments, "processed segment-level CN table"),
        (paths.segment_summary, "segment-level CN summary"),
        (paths.pilot_with_segments, "Level 3 pilot cohort with segments"),
        (paths.candidates_with_segments, "Level 3 candidates with segments"),
        (paths.purity_ploidy, "purity/ploidy by sample"),
        (paths.projects, "TCGA project metadata"),
    ]:
        require_file(path, description)

    pilot = read_tsv(paths.pilot_with_segments)
    pilot["sample_barcode"] = pilot["sample_barcode"].map(normalize_sample_barcode)
    pilot_samples = set(pilot["sample_barcode"].astype(str))
    purity_ploidy = read_tsv(paths.purity_ploidy)
    segment_summary_df = read_tsv(paths.segment_summary)
    if not segment_summary_df.empty:
        segment_summary_df["sample_barcode"] = segment_summary_df["sample_barcode"].map(normalize_sample_barcode)

    mutations = load_pilot_mutations(paths.mutations, pilot_samples)
    segments = load_segments(paths.segments, pilot_samples, config)
    annotated = annotate_mutations_with_segments(mutations, segments)
    annotated = attach_purity_ploidy_and_metadata(annotated, purity_ploidy, pilot)
    annotated = flag_clonal_input_candidates(annotated, config)
    annotated = annotated.rename(
        columns={
            "Chromosome": "chromosome_original",
            "Start_Position": "start_position_original",
            "End_Position": "end_position_original",
        }
    )
    for column in ANNOTATED_COLUMNS:
        if column not in annotated.columns:
            annotated[column] = pd.NA
    annotated_out = annotated[ANNOTATED_COLUMNS].copy()

    sample_summary = build_sample_summary(annotated_out, pilot, segment_summary_df, config)
    pyclone_preview = build_pyclone_preview(annotated_out)
    qc = qc_summary(annotated_out, sample_summary, pilot, segments, config)

    write_tsv_gz(annotated_out, paths.annotated_mutations)
    write_tsv(sample_summary, paths.summary_by_sample)
    write_tsv(pyclone_preview, paths.pyclone_preview)
    write_tsv(qc, paths.qc_summary)
    plot_coverage(sample_summary, paths.coverage_figure)
    plot_vaf_vs_segment_mean(annotated_out, paths.vaf_segment_figure, config)
    plot_depth_distribution(annotated_out, paths.depth_figure)

    return {
        "annotated": annotated_out,
        "summary": sample_summary,
        "pyclone_preview": pyclone_preview,
        "qc": qc,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Annotate Level 3 pilot mutations with local segment-level CN.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument("--min-total-depth", type=int, default=CONFIG["min_total_depth"])
    parser.add_argument(
        "--min-candidate-mutations-with-local-cn",
        type=int,
        default=CONFIG["min_candidate_mutations_with_local_cn"],
    )
    parser.add_argument("--min-annotation-rate", type=float, default=CONFIG["min_annotation_rate"])
    parser.add_argument("--include-silent", action="store_true", help="Include synonymous mutations in clonal preview.")
    parser.add_argument(
        "--allow-patient-level-segment-matching",
        action="store_true",
        help="Allow lower-confidence patient-level segment rows. Off by default.",
    )
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    config = dict(CONFIG)
    config["min_total_depth"] = args.min_total_depth
    config["min_candidate_mutations_with_local_cn"] = args.min_candidate_mutations_with_local_cn
    config["min_annotation_rate"] = args.min_annotation_rate
    config["candidate_nonsynonymous_only"] = not args.include_silent
    config["allow_patient_level_segment_matching"] = bool(args.allow_patient_level_segment_matching)
    outputs = run_mutation_local_cn_annotation(default_paths(root), config)
    logging.info(
        "Mutation local-CN annotation complete: %d mutations, %d annotated, %d eligible samples",
        len(outputs["annotated"]),
        int((outputs["annotated"]["local_cn_match_status"] == "matched").sum()),
        int(outputs["summary"]["eligible_for_copy_number_aware_clonal_input"].sum()),
    )


if __name__ == "__main__":
    main()

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
import pandas as pd
import pyarrow.parquet as pq

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


CONFIG = {
    "gain_like_threshold": 0.2,
    "loss_like_threshold": -0.2,
    "min_segments_for_usable_cn": 1,
    "attempt_remote_download": True,
}

SEGMENT_PATTERNS = [
    "*.seg",
    "*.seg.txt",
    "*.tsv",
    "*.txt",
    "*.maf.cnv",
    "copy_number_segments*.tsv",
    "segments_by_sample*.tsv",
    "segments_by_sample*.parquet",
]

SAMPLE_COLUMNS = ["sample_barcode", "Sample", "sample", "Sample_ID", "ID", "Tumor_Sample_Barcode", "aliquot_barcode"]
CHROMOSOME_COLUMNS = ["chromosome", "Chromosome", "chrom", "chr"]
START_COLUMNS = ["start", "Start", "loc.start", "chromStart", "Start_Position"]
END_COLUMNS = ["end", "End", "loc.end", "chromEnd", "End_Position"]
NUM_PROBES_COLUMNS = ["num_probes", "Num_Probes", "probes", "num.mark", "Num.Mark"]
SEGMENT_MEAN_COLUMNS = ["segment_mean", "Segment_Mean", "seg.mean", "log2_copy_ratio", "copy_ratio", "SegmentMean"]


@dataclass
class SegmentPaths:
    candidate_samples: Path
    pilot_cohort: Path
    mutation_parquet: Path
    purity_ploidy: Path
    aneuploidy: Path
    projects: Path
    processed_segments_tsv: Path
    processed_segments_gz: Path
    segment_summary: Path
    candidates_with_segments: Path
    pilot_with_segments: Path
    qc_summary: Path
    availability_by_project: Path
    availability_figure: Path
    count_distribution_figure: Path
    upgrade_status_figure: Path
    local_dirs: list[Path]


def default_paths(root: Path) -> SegmentPaths:
    return SegmentPaths(
        candidate_samples=project_path("results", "tables", "level3_candidate_samples.tsv", root=root),
        pilot_cohort=project_path("results", "tables", "level3_pilot_cohort.tsv", root=root),
        mutation_parquet=project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root=root),
        purity_ploidy=project_path("data", "processed", "features", "purity_ploidy_by_sample.tsv", root=root),
        aneuploidy=project_path("data", "processed", "features", "aneuploidy_by_sample.tsv", root=root),
        projects=project_path("data", "interim", "tcga_projects.tsv", root=root),
        processed_segments_tsv=project_path("data", "processed", "copy_number", "segments_by_sample.tsv", root=root),
        processed_segments_gz=project_path("data", "processed", "copy_number", "segments_by_sample.tsv.gz", root=root),
        segment_summary=project_path("data", "processed", "copy_number", "segment_level_cn_summary_by_sample.tsv", root=root),
        candidates_with_segments=project_path("results", "tables", "level3_candidate_samples_with_segments.tsv", root=root),
        pilot_with_segments=project_path("results", "tables", "level3_pilot_cohort_with_segments.tsv", root=root),
        qc_summary=project_path("results", "tables", "level3_segment_cn_qc_summary.tsv", root=root),
        availability_by_project=project_path("results", "tables", "level3_segment_cn_availability_by_project.tsv", root=root),
        availability_figure=project_path("results", "figures", "level3_segment_cn_availability_by_project.pdf", root=root),
        count_distribution_figure=project_path("results", "figures", "level3_segment_count_distribution.pdf", root=root),
        upgrade_status_figure=project_path("results", "figures", "level3_candidate_upgrade_status.pdf", root=root),
        local_dirs=[
            project_path("data", "raw", "copy_number", root=root),
            project_path("data", "raw", "copy_number_segments", root=root),
            project_path("data", "processed", "copy_number", root=root),
        ],
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


def first_present(columns: list[str], candidates: list[str]) -> str | None:
    lower_map = {column.lower(): column for column in columns}
    for candidate in candidates:
        if candidate in columns:
            return candidate
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    return None


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


def project_from_barcode(sample_barcode: str, project_lookup: dict[str, str] | None = None) -> tuple[str, str]:
    patient = normalize_patient_barcode(sample_barcode)
    if project_lookup and patient in project_lookup:
        code = project_lookup[patient]
        return f"TCGA-{code}", code
    return "TCGA-UNKNOWN", "UNKNOWN"


def find_local_segment_files(paths: SegmentPaths) -> list[Path]:
    found: set[Path] = set()
    excluded_names = {
        "arm_level_calls_by_sample.tsv",
        "segment_level_cn_summary_by_sample.tsv",
        "PANCAN_ArmCallsAndAneuploidyScore_092817.txt",
    }
    for directory in paths.local_dirs:
        if not directory.exists():
            continue
        for pattern in SEGMENT_PATTERNS:
            for path in directory.glob(pattern):
                if path.is_file():
                    lower_name = path.name.lower()
                    if path.name in excluded_names or "armcalls" in lower_name or "aneuploidy" in lower_name:
                        continue
                    if path.resolve() in {paths.processed_segments_tsv.resolve(), paths.processed_segments_gz.resolve()}:
                        continue
                    found.add(path)
    return sorted(found)


def read_segment_file(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pq.read_table(path).to_pandas()
    return pd.read_csv(path, sep=None, engine="python", comment="#", dtype=str)


def standardize_segment_table(raw: pd.DataFrame, source_file: Path, project_lookup: dict[str, str]) -> tuple[pd.DataFrame, list[str]]:
    warnings: list[str] = []
    if raw.empty:
        return empty_segments(), [f"{source_file.name}:empty_file"]
    columns = list(raw.columns)
    sample_col = first_present(columns, SAMPLE_COLUMNS)
    chrom_col = first_present(columns, CHROMOSOME_COLUMNS)
    start_col = first_present(columns, START_COLUMNS)
    end_col = first_present(columns, END_COLUMNS)
    num_probes_col = first_present(columns, NUM_PROBES_COLUMNS)
    segment_mean_col = first_present(columns, SEGMENT_MEAN_COLUMNS)
    required = {
        "sample": sample_col,
        "chromosome": chrom_col,
        "start": start_col,
        "end": end_col,
        "segment_mean": segment_mean_col,
    }
    missing = [name for name, column in required.items() if column is None]
    if missing:
        return empty_segments(), [f"{source_file.name}:missing_required_columns_{','.join(missing)}"]

    out = pd.DataFrame(
        {
            "sample_barcode": raw[sample_col].map(normalize_sample_barcode),
            "chromosome": raw[chrom_col].astype(str).str.replace("^chr", "", regex=True),
            "start": numeric(raw[start_col]),
            "end": numeric(raw[end_col]),
            "num_probes": numeric(raw[num_probes_col]) if num_probes_col else pd.NA,
            "segment_mean": numeric(raw[segment_mean_col]),
            "source_file": source_file.name,
            "source": "local_segment_file",
            "notes": "sample_level_match_preferred_from_segment_file",
        }
    )
    out["patient_barcode"] = out["sample_barcode"].map(normalize_patient_barcode)
    projects = out["patient_barcode"].map(lambda patient: project_from_barcode(patient, project_lookup))
    out["project_id"] = [item[0] for item in projects]
    out["project_code"] = [item[1] for item in projects]
    out = out[
        out["sample_barcode"].astype(str).str.startswith("TCGA-")
        & out["start"].notna()
        & out["end"].notna()
        & out["segment_mean"].notna()
    ].copy()
    out = out[out["end"] >= out["start"]]
    return out[segment_columns()], warnings


def empty_segments() -> pd.DataFrame:
    return pd.DataFrame(columns=segment_columns())


def segment_columns() -> list[str]:
    return [
        "sample_barcode",
        "patient_barcode",
        "project_id",
        "project_code",
        "chromosome",
        "start",
        "end",
        "num_probes",
        "segment_mean",
        "source_file",
        "source",
        "notes",
    ]


def project_lookup_from_candidates(candidates: pd.DataFrame) -> dict[str, str]:
    lookup = candidates[["patient_barcode", "project_code"]].dropna().drop_duplicates()
    return dict(zip(lookup["patient_barcode"], lookup["project_code"], strict=False))


def ingest_local_segments(paths: SegmentPaths, candidates: pd.DataFrame) -> tuple[pd.DataFrame, list[Path], list[str]]:
    files = find_local_segment_files(paths)
    project_lookup = project_lookup_from_candidates(candidates)
    frames: list[pd.DataFrame] = []
    warnings: list[str] = []
    for path in files:
        try:
            raw = read_segment_file(path)
            standardized, file_warnings = standardize_segment_table(raw, path, project_lookup)
            warnings.extend(file_warnings)
            if not standardized.empty:
                frames.append(standardized)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"{path.name}:parse_error:{exc}")
    if not frames:
        return empty_segments(), files, warnings
    segments = pd.concat(frames, ignore_index=True)
    segments = segments.drop_duplicates(
        subset=["sample_barcode", "chromosome", "start", "end", "segment_mean", "source_file"]
    )
    return segments[segment_columns()], files, warnings


def segment_summary(segments: pd.DataFrame, candidates: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    if segments.empty:
        return pd.DataFrame(
            columns=[
                "sample_barcode",
                "patient_barcode",
                "project_id",
                "project_code",
                "n_segments",
                "n_autosomal_segments",
                "median_segment_mean",
                "sd_segment_mean",
                "max_abs_segment_mean",
                "fraction_segments_gain_like",
                "fraction_segments_loss_like",
                "has_segment_level_cn",
                "segment_cn_quality_flag",
                "notes",
            ]
        )
    data = segments.copy()
    data["segment_mean"] = numeric(data["segment_mean"])
    data["chromosome_clean"] = data["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    data["is_autosomal"] = data["chromosome_clean"].isin([str(i) for i in range(1, 23)])
    grouped = data.groupby("sample_barcode", dropna=False)
    summary = grouped.agg(
        patient_barcode=("patient_barcode", "first"),
        project_id=("project_id", "first"),
        project_code=("project_code", "first"),
        n_segments=("sample_barcode", "size"),
        n_autosomal_segments=("is_autosomal", "sum"),
        median_segment_mean=("segment_mean", "median"),
        sd_segment_mean=("segment_mean", "std"),
        max_abs_segment_mean=("segment_mean", lambda x: x.abs().max()),
        fraction_segments_gain_like=("segment_mean", lambda x: (x >= config["gain_like_threshold"]).mean()),
        fraction_segments_loss_like=("segment_mean", lambda x: (x <= config["loss_like_threshold"]).mean()),
    ).reset_index()
    summary["has_segment_level_cn"] = summary["n_segments"] >= int(config["min_segments_for_usable_cn"])
    summary["segment_cn_quality_flag"] = summary.apply(
        lambda row: "usable_segment_cn" if row["has_segment_level_cn"] else "insufficient_segments",
        axis=1,
    )
    summary["notes"] = "segment_level_cn_available_sample_level_match"
    known = candidates[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates()
    summary = summary.drop(columns=["patient_barcode", "project_id", "project_code"]).merge(known, on="sample_barcode", how="left").merge(
        summary[["sample_barcode", "patient_barcode", "project_id", "project_code"]],
        on="sample_barcode",
        how="left",
        suffixes=("", "_from_segment"),
    )
    for column in ["patient_barcode", "project_id", "project_code"]:
        summary[column] = summary[column].fillna(summary[f"{column}_from_segment"])
        summary = summary.drop(columns=[f"{column}_from_segment"])
    ordered = [
        "sample_barcode",
        "patient_barcode",
        "project_id",
        "project_code",
        "n_segments",
        "n_autosomal_segments",
        "median_segment_mean",
        "sd_segment_mean",
        "max_abs_segment_mean",
        "fraction_segments_gain_like",
        "fraction_segments_loss_like",
        "has_segment_level_cn",
        "segment_cn_quality_flag",
        "notes",
    ]
    return summary[ordered]


def minimum_requirements_met(row: pd.Series) -> bool:
    return (
        bool(row.get("ref_alt_counts_available", False))
        and bool(row.get("purity_available", False))
        and bool(row.get("ploidy_available", False))
        and bool(row.get("copy_number_available", False))
        and str(row.get("exclusion_reason", "none")) == "none"
    )


def update_candidates_with_segments(candidates: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    merged = candidates.drop(columns=["has_segment_level_cn"], errors="ignore").merge(
        summary[
            [
                "sample_barcode",
                "n_segments",
                "n_autosomal_segments",
                "median_segment_mean",
                "sd_segment_mean",
                "max_abs_segment_mean",
                "fraction_segments_gain_like",
                "fraction_segments_loss_like",
                "has_segment_level_cn",
                "segment_cn_quality_flag",
            ]
        ],
        on="sample_barcode",
        how="left",
    )
    merged["has_segment_level_cn"] = merged["has_segment_level_cn"].where(merged["has_segment_level_cn"].notna(), False).astype(bool)
    upgrade_mask = merged["has_segment_level_cn"] & merged.apply(minimum_requirements_met, axis=1)
    limited_mask = (~merged["has_segment_level_cn"]) & merged["candidate_class"].isin(
        ["copy_number_aware_clonal_phylogeny_candidate", "clonal_clustering_candidate_limited_cn"]
    )
    merged.loc[upgrade_mask, "candidate_class"] = "copy_number_aware_clonal_phylogeny_candidate"
    merged.loc[upgrade_mask, "inclusion_status"] = merged.loc[upgrade_mask, "inclusion_status"].replace(
        {"eligible_not_selected": "eligible_not_selected", "selected_for_pilot": "selected_for_pilot"}
    )
    merged.loc[limited_mask, "candidate_class"] = "clonal_clustering_candidate_limited_cn"
    merged["notes"] = merged.apply(update_notes, axis=1)
    return merged


def update_notes(row: pd.Series) -> str:
    notes = [item for item in str(row.get("notes", "")).split(";") if item and item != "NA"]
    notes = [
        item
        for item in notes
        if item
        not in {
            "segment_level_cn_missing_suitable_for_clonal_clustering_prototype_not_full_cn_aware_phylogeny",
            "segment_level_cn_available",
        }
    ]
    if bool(row.get("has_segment_level_cn", False)):
        notes.append("segment_level_cn_available_candidate_upgraded_if_minimum_filters_pass")
    elif row.get("candidate_class") == "clonal_clustering_candidate_limited_cn":
        notes.append("segment_level_cn_missing_suitable_for_clonal_clustering_prototype_not_full_cn_aware_phylogeny")
    return ";".join(dict.fromkeys(notes))


def update_pilot_with_segments(pilot: pd.DataFrame, updated_candidates: pd.DataFrame) -> pd.DataFrame:
    selected = set(pilot["sample_barcode"])
    upgraded = updated_candidates[updated_candidates["sample_barcode"].isin(selected)].copy()
    if upgraded.empty:
        return pilot.copy()
    any_segments = updated_candidates["has_segment_level_cn"].sum() > 0
    if not any_segments:
        upgraded["notes"] = upgraded["notes"].astype(str) + ";pilot_preserved_no_segment_level_cn_available"
        return upgraded
    priority = updated_candidates[
        updated_candidates["candidate_class"] == "copy_number_aware_clonal_phylogeny_candidate"
    ].sort_values(["project_code", "candidate_score"], ascending=[True, False])
    if priority.empty:
        upgraded["notes"] = upgraded["notes"].astype(str) + ";pilot_preserved_no_candidate_overlap_with_segments"
        return upgraded
    counts = pilot["project_code"].value_counts().to_dict()
    replacement_parts = []
    for project, count in counts.items():
        project_priority = priority[priority["project_code"] == project].head(count)
        if len(project_priority) < count:
            project_priority = pd.concat(
                [
                    project_priority,
                    updated_candidates[
                        (updated_candidates["project_code"] == project)
                        & (~updated_candidates["sample_barcode"].isin(project_priority["sample_barcode"]))
                        & updated_candidates["candidate_class"].isin(
                            [
                                "copy_number_aware_clonal_phylogeny_candidate",
                                "clonal_clustering_candidate_limited_cn",
                            ]
                        )
                    ]
                    .sort_values("candidate_score", ascending=False)
                    .head(count - len(project_priority)),
                ]
            )
        replacement_parts.append(project_priority)
    out = pd.concat(replacement_parts, ignore_index=True) if replacement_parts else upgraded
    out["inclusion_status"] = "selected_for_pilot"
    return out


def availability_by_project(candidates: pd.DataFrame, pilot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (project_id, project_code), group in candidates.groupby(["project_id", "project_code"], dropna=False):
        pilot_group = pilot[pilot["project_code"] == project_code]
        n_candidates = len(group)
        n_with = int(group["has_segment_level_cn"].sum())
        n_pilot = len(pilot_group)
        n_pilot_with = int(pilot_group["has_segment_level_cn"].sum()) if "has_segment_level_cn" in pilot_group else 0
        rows.append(
            {
                "project_id": project_id,
                "project_code": project_code,
                "n_level3_candidates": n_candidates,
                "n_candidates_with_segments": n_with,
                "pct_candidates_with_segments": round(100 * n_with / n_candidates, 3) if n_candidates else 0,
                "n_pilot_samples": n_pilot,
                "n_pilot_samples_with_segments": n_pilot_with,
                "pct_pilot_with_segments": round(100 * n_pilot_with / n_pilot, 3) if n_pilot else 0,
                "notes": "segment_cn_available" if n_with else "no_segment_level_cn_overlap",
            }
        )
    return pd.DataFrame(rows).sort_values(["n_candidates_with_segments", "n_level3_candidates"], ascending=[False, False])


def manual_download_instructions() -> str:
    return (
        "No parseable segment-level CN rows were found locally. Manual option: use the GDC Data Portal/API; "
        "filter Program=TCGA, Data Category=Copy Number Variation, Data Type=Masked Copy Number Segment; "
        "download segment files and place them under data/raw/copy_number_segments/ or data/raw/copy_number/. "
        "Supported patterns include *.seg, *.seg.txt, *.tsv, *.txt, *.maf.cnv, copy_number_segments*.tsv, "
        "segments_by_sample*.tsv, and segments_by_sample*.parquet."
    )


def qc_summary(
    local_files: list[Path],
    remote_attempted: bool,
    remote_succeeded: bool,
    segments: pd.DataFrame,
    summary: pd.DataFrame,
    candidates: pd.DataFrame,
    updated_candidates: pd.DataFrame,
    pilot: pd.DataFrame,
    updated_pilot: pd.DataFrame,
    mutation_samples: set[str],
    warnings: list[str],
) -> pd.DataFrame:
    candidate_overlap = int(updated_candidates["has_segment_level_cn"].sum())
    pilot_overlap = int(updated_pilot["has_segment_level_cn"].sum()) if "has_segment_level_cn" in updated_pilot else 0
    upgraded = int((updated_candidates["candidate_class"] == "copy_number_aware_clonal_phylogeny_candidate").sum())
    remaining_limited = int((updated_candidates["candidate_class"] == "clonal_clustering_candidate_limited_cn").sum())
    projects = sorted(summary["project_code"].dropna().unique()) if not summary.empty else []
    rows = [
        ("global", "local_segment_files_found", ";".join(str(path) for path in local_files) if local_files else "none"),
        ("global", "remote_download_attempted", remote_attempted),
        ("global", "remote_download_succeeded", remote_succeeded),
        ("global", "total_segment_rows_processed", len(segments)),
        ("global", "samples_with_segment_level_cn", len(summary)),
        ("global", "tcga_projects_represented_in_segment_level_cn", ";".join(projects) if projects else "none"),
        ("global", "overlap_with_mutation_layer_samples", len(set(summary["sample_barcode"]) & mutation_samples) if not summary.empty else 0),
        ("global", "overlap_with_level3_candidates", candidate_overlap),
        ("global", "overlap_with_level3_pilot_cohort", pilot_overlap),
        ("global", "candidates_upgraded_to_copy_number_aware", upgraded),
        ("global", "candidates_remaining_limited_cn", remaining_limited),
        ("global", "parsing_warnings", ";".join(warnings) if warnings else "none"),
        ("global", "limitations", "none" if len(summary) else manual_download_instructions()),
    ]
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def mutation_sample_set(mutation_parquet: Path) -> set[str]:
    table = pq.read_table(mutation_parquet, columns=["sample_barcode"])
    return set(table.to_pandas()["sample_barcode"].dropna().astype(str))


def plot_availability(availability: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = availability.sort_values("n_candidates_with_segments", ascending=False)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.bar(data["project_code"], data["n_candidates_with_segments"], color="#4daf4a")
    ax.set_ylabel("Candidates with segment-level CN")
    ax.set_title("Level 3 segment CN availability by project")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_segment_counts(summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(9, 6))
    if summary.empty:
        ax.text(0.5, 0.5, "No segment-level CN rows available", ha="center", va="center")
        ax.axis("off")
    else:
        ax.hist(summary["n_segments"], bins=30, color="#377eb8")
        ax.set_xlabel("Segments per sample")
        ax.set_ylabel("Samples")
        ax.set_title("Segment-level CN segment-count distribution")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_upgrade_status(updated_candidates: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    counts = updated_candidates["candidate_class"].value_counts()
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.bar(counts.index, counts.values, color="#984ea3")
    ax.set_ylabel("Samples")
    ax.set_title("Level 3 candidate class after segment CN join")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def run_segment_cn(paths: SegmentPaths, config: dict = CONFIG) -> dict[str, pd.DataFrame]:
    for path, description in [
        (paths.candidate_samples, "Level 3 candidate sample table"),
        (paths.pilot_cohort, "Level 3 pilot cohort table"),
        (paths.mutation_parquet, "standardized mutation parquet"),
        (paths.purity_ploidy, "purity/ploidy table"),
        (paths.aneuploidy, "aneuploidy table"),
        (paths.projects, "TCGA project metadata"),
    ]:
        require_file(path, description)
    candidates = pd.read_csv(paths.candidate_samples, sep="\t")
    pilot = pd.read_csv(paths.pilot_cohort, sep="\t")
    segments, local_files, warnings = ingest_local_segments(paths, candidates)
    remote_attempted = False
    remote_succeeded = False
    if segments.empty and config.get("attempt_remote_download", True):
        remote_attempted = True
        warnings.append("remote_download_not_performed_in_restricted_environment_manual_download_required")
    summary = segment_summary(segments, candidates, config)
    updated_candidates = update_candidates_with_segments(candidates, summary)
    updated_pilot = update_pilot_with_segments(pilot, updated_candidates)
    availability = availability_by_project(updated_candidates, updated_pilot)
    mutation_samples = mutation_sample_set(paths.mutation_parquet)
    qc = qc_summary(
        local_files,
        remote_attempted,
        remote_succeeded,
        segments,
        summary,
        candidates,
        updated_candidates,
        pilot,
        updated_pilot,
        mutation_samples,
        warnings,
    )

    write_tsv_gz(segments[segment_columns()], paths.processed_segments_gz)
    write_tsv(summary, paths.segment_summary)
    write_tsv(updated_candidates, paths.candidates_with_segments)
    write_tsv(updated_pilot, paths.pilot_with_segments)
    write_tsv(qc, paths.qc_summary)
    write_tsv(availability, paths.availability_by_project)
    plot_availability(availability, paths.availability_figure)
    plot_segment_counts(summary, paths.count_distribution_figure)
    plot_upgrade_status(updated_candidates, paths.upgrade_status_figure)
    return {
        "segments": segments,
        "summary": summary,
        "candidates": updated_candidates,
        "pilot": updated_pilot,
        "availability": availability,
        "qc": qc,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Ingest segment-level copy-number data and update Level 3 candidate classes.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    paths = default_paths(root)
    logging.info("Checking local segment-level copy-number files before attempting remote acquisition")
    outputs = run_segment_cn(paths, CONFIG)
    logging.info(
        "Segment CN join complete: %d segment rows, %d samples with segment CN, %d candidates upgraded",
        len(outputs["segments"]),
        len(outputs["summary"]),
        int((outputs["candidates"]["candidate_class"] == "copy_number_aware_clonal_phylogeny_candidate").sum()),
    )


if __name__ == "__main__":
    main()

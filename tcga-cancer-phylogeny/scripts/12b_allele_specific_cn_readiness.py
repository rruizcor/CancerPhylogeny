#!/usr/bin/env python

from __future__ import annotations

import argparse
import fnmatch
import gzip
import logging
import os
import re
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

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


CONFIG = {
    "min_total_depth": 20,
    "min_pyclone_eligible_mutations": 50,
    "min_pyclone_mutation_coverage_fraction": 0.80,
    "min_phylowgs_eligible_mutations": 100,
    "min_phylowgs_mutation_coverage_fraction": 0.90,
    "min_phylowgs_expected_clusters": 2,
    "integer_cn_tolerance": 0.05,
    "recommended_external_tool": "FACETS",
}

CANDIDATE_FILE_PATTERNS = [
    "*facets*.tsv",
    "*facets*.txt",
    "*facets*.tsv.gz",
    "*facets*.txt.gz",
    "*cncf*.tsv",
    "*cncf*.txt",
    "*cncf*.tsv.gz",
    "*cncf*.txt.gz",
    "*sequenza*.tsv",
    "*sequenza*.txt",
    "*sequenza*.tsv.gz",
    "*sequenza*.txt.gz",
    "*segments.txt",
    "*segments.txt.gz",
    "*ascat*.tsv",
    "*ascat*.txt",
    "*ascat*.tsv.gz",
    "*ascat*.txt.gz",
    "*absolute*.tsv",
    "*absolute*.txt",
    "*absolute*.tsv.gz",
    "*absolute*.txt.gz",
    "*allele_specific*.tsv",
    "*allele_specific*.txt",
    "*allele_specific*.tsv.gz",
    "*allele_specific*.txt.gz",
    "*major_minor_cn*.tsv",
    "*major_minor_cn*.txt",
    "*major_minor_cn*.tsv.gz",
    "*major_minor_cn*.txt.gz",
]

SAMPLE_COLUMNS = [
    "sample",
    "sample_id",
    "sample_barcode",
    "Tumor_Sample_Barcode",
    "Sample",
    "Sample_ID",
    "tumor_sample",
    "tumor_barcode",
]
CHROMOSOME_COLUMNS = ["chromosome", "chrom", "Chromosome", "chr"]
START_COLUMNS = ["start", "Start", "loc.start", "chromStart", "start_position", "Start_Position"]
END_COLUMNS = ["end", "End", "loc.end", "chromEnd", "end_position", "End_Position"]
TOTAL_CN_COLUMNS = [
    "total_cn",
    "tcn",
    "total_copy_number",
    "CNt",
    "tcn.em",
    "TCN",
    "nTotal",
    "modal_total_cn",
    "Modal_Total_CN",
]
MAJOR_CN_COLUMNS = [
    "major_cn",
    "major",
    "major_copy_number",
    "A",
    "nMajor",
    "nA",
    "major.allele",
    "Modal_HSCN_1",
]
MINOR_CN_COLUMNS = [
    "minor_cn",
    "minor",
    "minor_copy_number",
    "B",
    "lcn",
    "lcn.em",
    "nMinor",
    "nB",
    "minor.allele",
    "Modal_HSCN_2",
]
LOH_COLUMNS = ["LOH", "loh", "loh_status", "LOH_status", "loss_of_heterozygosity"]
MAFR_COLUMNS = ["mafR", "mafr", "minor_allele_fraction", "BAF", "baf"]
CELLULAR_FRACTION_COLUMNS = ["cf", "cellular_fraction", "cellularity", "cf.em", "cancer_cell_fraction"]

STANDARDIZED_SEGMENT_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "chromosome",
    "start",
    "end",
    "total_cn",
    "major_cn",
    "minor_cn",
    "loh_status",
    "cellular_fraction",
    "source_file",
    "source_tool",
    "notes",
]

AVAILABILITY_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "in_pilot_cohort",
    "in_limited_vaf_pilot",
    "has_segment_mean_cn",
    "has_allele_specific_cn",
    "allele_specific_cn_source",
    "n_allele_specific_segments",
    "has_total_cn",
    "has_major_cn",
    "has_minor_cn",
    "has_loh_status",
    "can_prepare_cn_aware_pyclone_input",
    "reason_not_ready",
    "recommended_next_action",
    "notes",
]

MUTATION_READINESS_COLUMNS = [
    "sample_barcode",
    "n_mutations_with_segment_mean_cn",
    "n_mutations_with_allele_specific_cn",
    "pct_mutations_with_allele_specific_cn",
    "n_eligible_mutations_with_allele_specific_cn",
    "has_purity",
    "has_ref_alt_counts",
    "can_build_pyclone_vi_input",
    "can_build_phylowgs_input",
    "reason_not_ready",
]

MANIFEST_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "recommended_tool",
    "tumor_bam_required",
    "normal_bam_required",
    "tumor_vcf_required",
    "normal_vcf_required",
    "snp_pileup_required",
    "matched_normal_available_unknown",
    "current_available_inputs",
    "missing_required_inputs",
    "recommended_command_template",
    "notes",
]


@dataclass
class AlleleSpecificCnPaths:
    pilot_with_segments: Path
    limited_vaf_pilot: Path
    mutations_with_local_cn: Path
    segment_mean_cn: Path
    segment_mean_summary: Path
    purity_ploidy: Path
    clonal_input_readiness: Path
    limited_vaf_complexity: Path
    standardized_segments: Path
    availability_by_sample: Path
    mutation_readiness: Path
    external_manifest: Path
    qc_summary: Path
    planning_report: Path
    readiness_figure: Path
    missing_inputs_figure: Path
    local_dirs: list[Path]
    raw_data_root: Path


def default_paths(root: Path, local_dirs: list[Path] | None = None) -> AlleleSpecificCnPaths:
    return AlleleSpecificCnPaths(
        pilot_with_segments=project_path(
            "results", "tables", "level3_pilot_cohort_with_segments.tsv", root=root
        ),
        limited_vaf_pilot=project_path(
            "results", "tables", "level3_limited_vaf_clustering_pilot_samples.tsv", root=root
        ),
        mutations_with_local_cn=project_path(
            "data", "processed", "clonal", "level3_mutations_with_local_cn.tsv.gz", root=root
        ),
        segment_mean_cn=project_path(
            "data", "processed", "copy_number", "segments_by_sample.tsv.gz", root=root
        ),
        segment_mean_summary=project_path(
            "data", "processed", "copy_number", "segment_level_cn_summary_by_sample.tsv", root=root
        ),
        purity_ploidy=project_path(
            "data", "processed", "features", "purity_ploidy_by_sample.tsv", root=root
        ),
        clonal_input_readiness=project_path(
            "results", "tables", "level3_clonal_input_readiness_by_sample.tsv", root=root
        ),
        limited_vaf_complexity=project_path(
            "results", "tables", "level3_limited_vaf_clonal_complexity_by_sample.tsv", root=root
        ),
        standardized_segments=project_path(
            "data",
            "processed",
            "allele_specific_cn",
            "allele_specific_segments_by_sample.tsv.gz",
            root=root,
        ),
        availability_by_sample=project_path(
            "results", "tables", "level3_allele_specific_cn_availability_by_sample.tsv", root=root
        ),
        mutation_readiness=project_path(
            "results", "tables", "level3_mutation_to_allele_specific_cn_readiness.tsv", root=root
        ),
        external_manifest=project_path(
            "results", "tables", "level3_allele_specific_cn_external_inference_manifest.tsv", root=root
        ),
        qc_summary=project_path(
            "results", "tables", "level3_allele_specific_cn_readiness_qc_summary.tsv", root=root
        ),
        planning_report=project_path(
            "results", "reports", "level3_allele_specific_cn_and_pyclone_plan.md", root=root
        ),
        readiness_figure=project_path(
            "results", "figures", "level3_allele_specific_cn_readiness_by_project.pdf", root=root
        ),
        missing_inputs_figure=project_path(
            "results", "figures", "level3_allele_specific_cn_missing_inputs.pdf", root=root
        ),
        local_dirs=local_dirs
        or [
            project_path("data", "raw", "allele_specific_cn", root=root),
            project_path("data", "processed", "allele_specific_cn", root=root),
            project_path("data", "raw", "facets", root=root),
            project_path("data", "raw", "sequenza", root=root),
            project_path("data", "raw", "ascat", root=root),
            project_path("data", "raw", "absolute", root=root),
        ],
        raw_data_root=project_path("data", "raw", root=root),
    )


def read_tsv(path: Path, label: str) -> pd.DataFrame:
    require_file(path, label)
    return pd.read_csv(path, sep="\t", keep_default_na=True)


def write_tsv(data: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data.to_csv(path, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_tsv_gz(data: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    with gzip.open(path, "wt") as handle:
        data.to_csv(handle, sep="\t", index=False, na_rep="NA")
    logging.info("Wrote %s", path)


def write_text(text: str, path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")
    logging.info("Wrote %s", path)


def require_columns(data: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise ValueError(f"{label} lacks required columns: {missing}")


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def normalize_patient_barcode(barcode: Any) -> str:
    text = str(barcode).strip().upper()
    parts = text.split("-")
    if len(parts) >= 3 and parts[0] == "TCGA":
        return "-".join(parts[:3])
    return text[:12]


def normalize_sample_barcode(barcode: Any) -> str:
    text = str(barcode).strip().upper()
    parts = text.split("-")
    if len(parts) >= 4 and parts[0] == "TCGA":
        return "-".join(parts[:4])[:16]
    return text[:16]


def normalize_chromosome(value: Any) -> str:
    text = str(value).strip().upper().replace("CHR", "")
    if text in {"", "NAN", "NONE", "NA"}:
        return "NA"
    if text == "23":
        return "X"
    if text == "24":
        return "Y"
    if text in {"25", "M", "MT"}:
        return "MT"
    try:
        return str(int(float(text)))
    except ValueError:
        return text


def normalized_column(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def first_present(columns: list[str], candidates: list[str]) -> str | None:
    lookup: dict[str, str] = {}
    for column in columns:
        lookup.setdefault(normalized_column(column), column)
    for candidate in candidates:
        if normalized_column(candidate) in lookup:
            return lookup[normalized_column(candidate)]
    return None


def empty_standardized_segments() -> pd.DataFrame:
    return pd.DataFrame(columns=STANDARDIZED_SEGMENT_COLUMNS)


def is_candidate_allele_specific_file(path: Path) -> bool:
    lower_name = path.name.lower()
    return any(fnmatch.fnmatch(lower_name, pattern.lower()) for pattern in CANDIDATE_FILE_PATTERNS)


def find_local_allele_specific_files(paths: AlleleSpecificCnPaths) -> list[Path]:
    found: set[Path] = set()
    output_resolved = paths.standardized_segments.resolve()
    for directory in paths.local_dirs:
        if not directory.exists():
            continue
        for path in directory.rglob("*"):
            if not path.is_file() or path.resolve() == output_resolved:
                continue
            if is_candidate_allele_specific_file(path):
                found.add(path)
    return sorted(found)


def read_candidate_file(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep=None,
        engine="python",
        comment="#",
        dtype=str,
        compression="infer",
    )


def infer_source_tool(path: Path, columns: list[str]) -> str:
    name = path.name.lower()
    normalized = {normalized_column(column) for column in columns}
    if "facets" in name or "cncf" in name or {"tcnem", "lcnem"} & normalized:
        return "FACETS"
    if "sequenza" in name or {"cnt", "depthratio"} & normalized:
        return "Sequenza"
    if "ascat" in name or {"nmajor", "nminor"} & normalized:
        return "ASCAT"
    if "absolute" in name or {"modaltotalcn", "modalhscn1", "modalhscn2"} & normalized:
        return "ABSOLUTE-style"
    return "generic_allele_specific_cn"


def barcode_from_path(path: Path) -> str | None:
    match = re.search(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]?", str(path).upper())
    return normalize_sample_barcode(match.group(0)) if match else None


def pilot_metadata_maps(pilot: pd.DataFrame) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    meta = pilot[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates().copy()
    meta["sample_barcode"] = meta["sample_barcode"].map(normalize_sample_barcode)
    meta["patient_barcode"] = meta["patient_barcode"].map(normalize_patient_barcode)
    sample_map = meta.set_index("sample_barcode").to_dict(orient="index")
    patient_counts = meta.groupby("patient_barcode")["sample_barcode"].nunique()
    unique_patients = set(patient_counts[patient_counts == 1].index)
    patient_map = dict(
        zip(
            meta.loc[meta["patient_barcode"].isin(unique_patients), "patient_barcode"],
            meta.loc[meta["patient_barcode"].isin(unique_patients), "sample_barcode"],
            strict=False,
        )
    )
    return sample_map, patient_map


def resolve_pilot_sample(value: Any, sample_map: dict[str, Any], patient_map: dict[str, str]) -> str | None:
    sample = normalize_sample_barcode(value)
    if sample in sample_map:
        return sample
    patient = normalize_patient_barcode(value)
    return patient_map.get(patient)


def integer_cn(series: pd.Series, tolerance: float) -> pd.Series:
    values = numeric(series)
    rounded = values.round()
    valid = values.notna() & (values >= 0) & ((values - rounded).abs() <= tolerance)
    return rounded.where(valid)


def normalize_loh(value: Any) -> str | pd.NA:
    if pd.isna(value) or str(value).strip() == "":
        return pd.NA
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "loh", "loss", "loss_of_heterozygosity"}:
        return "LOH"
    if text in {"0", "false", "no", "n", "retained", "heterozygous", "no_loh"}:
        return "no_LOH"
    return f"reported:{str(value).strip()}"


def standardize_allele_specific_table(
    raw: pd.DataFrame,
    source_file: Path,
    pilot: pd.DataFrame,
    config: dict[str, Any] = CONFIG,
) -> tuple[pd.DataFrame, list[str]]:
    warnings: list[str] = []
    if raw.empty:
        return empty_standardized_segments(), [f"{source_file.name}:empty_file"]
    columns = list(raw.columns)
    sample_col = first_present(columns, SAMPLE_COLUMNS)
    chrom_col = first_present(columns, CHROMOSOME_COLUMNS)
    start_col = first_present(columns, START_COLUMNS)
    end_col = first_present(columns, END_COLUMNS)
    total_col = first_present(columns, TOTAL_CN_COLUMNS)
    major_col = first_present(columns, MAJOR_CN_COLUMNS)
    minor_col = first_present(columns, MINOR_CN_COLUMNS)
    loh_col = first_present(columns, LOH_COLUMNS)
    mafr_col = first_present(columns, MAFR_COLUMNS)
    cellular_fraction_col = first_present(columns, CELLULAR_FRACTION_COLUMNS)

    missing_coordinates = [
        name
        for name, column in {"chromosome": chrom_col, "start": start_col, "end": end_col}.items()
        if column is None
    ]
    if missing_coordinates:
        return empty_standardized_segments(), [
            f"{source_file.name}:missing_required_columns_{','.join(missing_coordinates)}"
        ]
    if total_col is None and major_col is None and minor_col is None:
        return empty_standardized_segments(), [
            f"{source_file.name}:no_integer_total_major_or_minor_cn_fields;segment_mean_not_allele_specific"
        ]
    fallback_sample = barcode_from_path(source_file)
    if sample_col is None and fallback_sample is None:
        return empty_standardized_segments(), [f"{source_file.name}:missing_sample_column_and_filename_barcode"]

    sample_map, patient_map = pilot_metadata_maps(pilot)
    sample_values = raw[sample_col] if sample_col else pd.Series(fallback_sample, index=raw.index)
    resolved_samples = sample_values.map(lambda value: resolve_pilot_sample(value, sample_map, patient_map))
    unmatched_count = int(resolved_samples.isna().sum())
    if unmatched_count:
        warnings.append(f"{source_file.name}:{unmatched_count}_rows_not_matched_to_pilot_samples")

    tolerance = float(config["integer_cn_tolerance"])
    total = integer_cn(raw[total_col], tolerance) if total_col else pd.Series(np.nan, index=raw.index)
    major = integer_cn(raw[major_col], tolerance) if major_col else pd.Series(np.nan, index=raw.index)
    minor = integer_cn(raw[minor_col], tolerance) if minor_col else pd.Series(np.nan, index=raw.index)
    notes = pd.Series("integer_allele_specific_cn_ingested_without_segment_mean_inference", index=raw.index, dtype=object)

    both_alleles = major.notna() & minor.notna()
    if both_alleles.any():
        allele_max = pd.concat([major, minor], axis=1).max(axis=1)
        allele_min = pd.concat([major, minor], axis=1).min(axis=1)
        reordered = both_alleles & ((allele_max != major) | (allele_min != minor))
        major = major.where(~both_alleles, allele_max)
        minor = minor.where(~both_alleles, allele_min)
        notes.loc[reordered] += ";allele_channels_ordered_to_major_minor"

    derive_major = major.isna() & total.notna() & minor.notna()
    major.loc[derive_major] = total.loc[derive_major] - minor.loc[derive_major]
    notes.loc[derive_major] += ";major_cn_derived_from_total_minus_minor"
    derive_minor = minor.isna() & total.notna() & major.notna()
    minor.loc[derive_minor] = total.loc[derive_minor] - major.loc[derive_minor]
    notes.loc[derive_minor] += ";minor_cn_derived_from_total_minus_major"
    derive_total = total.isna() & major.notna() & minor.notna()
    total.loc[derive_total] = major.loc[derive_total] + minor.loc[derive_total]
    notes.loc[derive_total] += ";total_cn_derived_from_major_plus_minor"

    both_alleles = major.notna() & minor.notna()
    if both_alleles.any():
        allele_max = pd.concat([major, minor], axis=1).max(axis=1)
        allele_min = pd.concat([major, minor], axis=1).min(axis=1)
        major = major.where(~both_alleles, allele_max)
        minor = minor.where(~both_alleles, allele_min)
    invalid_negative = (major.notna() & (major < 0)) | (minor.notna() & (minor < 0))
    inconsistent = total.notna() & major.notna() & minor.notna() & ((total - major - minor).abs() > tolerance)
    invalid = invalid_negative | inconsistent
    if invalid.any():
        warnings.append(f"{source_file.name}:{int(invalid.sum())}_rows_with_inconsistent_integer_cn")
        total.loc[invalid] = np.nan
        major.loc[invalid] = np.nan
        minor.loc[invalid] = np.nan
        notes.loc[invalid] += ";invalid_or_inconsistent_cn_excluded_from_readiness"

    if loh_col:
        loh_status = raw[loh_col].map(normalize_loh)
    else:
        loh_status = pd.Series(pd.NA, index=raw.index, dtype=object)
    inferred_loh = loh_status.isna() & major.notna() & minor.notna()
    loh_status.loc[inferred_loh & total.eq(0)] = "complete_copy_loss"
    loh_status.loc[inferred_loh & total.ne(0) & minor.eq(0)] = "LOH"
    loh_status.loc[inferred_loh & minor.gt(0)] = "no_LOH"
    notes.loc[inferred_loh] += ";loh_status_derived_from_integer_minor_cn"
    if mafr_col:
        notes += ";mafR_or_BAF_present_not_used_to_infer_integer_cn"

    if cellular_fraction_col:
        cellular_fraction = numeric(raw[cellular_fraction_col])
        percent_mask = cellular_fraction.gt(1) & cellular_fraction.le(100)
        cellular_fraction.loc[percent_mask] = cellular_fraction.loc[percent_mask] / 100.0
        notes.loc[percent_mask] += ";cellular_fraction_converted_from_percent"
        cellular_fraction = cellular_fraction.where(cellular_fraction.between(0, 1, inclusive="both"))
    else:
        cellular_fraction = pd.Series(np.nan, index=raw.index)

    out = pd.DataFrame(
        {
            "sample_barcode": resolved_samples,
            "chromosome": raw[chrom_col].map(normalize_chromosome),
            "start": numeric(raw[start_col]),
            "end": numeric(raw[end_col]),
            "total_cn": total,
            "major_cn": major,
            "minor_cn": minor,
            "loh_status": loh_status,
            "cellular_fraction": cellular_fraction,
            "source_file": str(source_file),
            "source_tool": infer_source_tool(source_file, columns),
            "notes": notes,
        }
    )
    out = out[
        out["sample_barcode"].notna()
        & out["chromosome"].ne("NA")
        & out["chromosome"].ne("MT")
        & out["start"].notna()
        & out["end"].notna()
        & (out["end"] >= out["start"])
        & (out[["total_cn", "major_cn", "minor_cn"]].notna().any(axis=1))
    ].copy()
    if out.empty:
        return empty_standardized_segments(), warnings + [f"{source_file.name}:no_usable_pilot_rows"]
    metadata = pd.DataFrame.from_dict(sample_map, orient="index").reset_index(names="sample_barcode")
    out = out.merge(metadata, on="sample_barcode", how="left")
    out = out[STANDARDIZED_SEGMENT_COLUMNS]
    return out, warnings


def ingest_local_allele_specific_cn(
    paths: AlleleSpecificCnPaths,
    pilot: pd.DataFrame,
    config: dict[str, Any] = CONFIG,
) -> tuple[pd.DataFrame, list[Path], list[Path], list[str]]:
    files = find_local_allele_specific_files(paths)
    parsed_files: list[Path] = []
    warnings: list[str] = []
    frames: list[pd.DataFrame] = []
    for path in files:
        try:
            raw = read_candidate_file(path)
            standardized, file_warnings = standardize_allele_specific_table(raw, path, pilot, config)
            warnings.extend(file_warnings)
            if not standardized.empty:
                frames.append(standardized)
                parsed_files.append(path)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"{path.name}:parse_error:{exc}")
    if not frames:
        return empty_standardized_segments(), files, parsed_files, warnings
    segments = pd.concat(frames, ignore_index=True)
    segments = segments.drop_duplicates(
        subset=[
            "sample_barcode",
            "chromosome",
            "start",
            "end",
            "total_cn",
            "major_cn",
            "minor_cn",
            "source_file",
        ]
    )
    segments = segments.sort_values(
        ["sample_barcode", "chromosome", "start", "end", "source_tool"], kind="stable"
    ).reset_index(drop=True)
    return segments[STANDARDIZED_SEGMENT_COLUMNS], files, parsed_files, warnings


def mark_mutation_overlaps(mutations: pd.DataFrame, segments: pd.DataFrame) -> pd.Series:
    overlap = pd.Series(False, index=mutations.index)
    usable = segments[segments["major_cn"].notna() & segments["minor_cn"].notna()].copy()
    if mutations.empty or usable.empty:
        return overlap
    usable["chromosome_norm"] = usable["chromosome"].map(normalize_chromosome)
    usable["start"] = numeric(usable["start"])
    usable["end"] = numeric(usable["end"])
    mutation_chromosomes = mutations["chromosome"].map(normalize_chromosome)
    mutation_positions = numeric(mutations["position"])
    for (sample, chromosome), indices in mutations.groupby(
        [mutations["sample_barcode"].map(normalize_sample_barcode), mutation_chromosomes], dropna=False
    ).groups.items():
        if chromosome in {"NA", "MT"}:
            continue
        sample_segments = usable[
            (usable["sample_barcode"] == sample) & (usable["chromosome_norm"] == chromosome)
        ].sort_values(["start", "end"], kind="stable")
        if sample_segments.empty:
            continue
        starts = sample_segments["start"].to_numpy(dtype=float)
        ends = sample_segments["end"].to_numpy(dtype=float)
        positions = mutation_positions.loc[indices].to_numpy(dtype=float)
        candidate = np.searchsorted(starts, positions, side="right") - 1
        clipped = np.clip(candidate, 0, len(ends) - 1)
        valid = (candidate >= 0) & np.isfinite(positions) & (positions <= ends[clipped])
        overlap.loc[np.asarray(list(indices))[valid]] = True
    return overlap


def build_mutation_readiness(
    pilot: pd.DataFrame,
    limited_vaf_pilot: pd.DataFrame,
    mutations: pd.DataFrame,
    segments: pd.DataFrame,
    purity_ploidy: pd.DataFrame,
    complexity: pd.DataFrame,
    config: dict[str, Any] = CONFIG,
) -> pd.DataFrame:
    data = mutations.copy()
    data["sample_barcode"] = data["sample_barcode"].map(normalize_sample_barcode)
    for column in ["t_ref_count", "t_alt_count", "total_depth", "observed_vaf", "local_segment_mean"]:
        data[column] = numeric(data[column]) if column in data.columns else np.nan
    segment_mean_match = data["local_segment_mean"].notna()
    if "local_cn_match_status" in data.columns:
        segment_mean_match &= data["local_cn_match_status"].astype(str).eq("matched")
    ascn_overlap = mark_mutation_overlaps(data, segments)
    has_counts = data["t_ref_count"].notna() & data["t_alt_count"].notna()
    base_eligible = has_counts & data["total_depth"].ge(int(config["min_total_depth"]))
    if "observed_vaf" in data.columns:
        base_eligible &= data["observed_vaf"].gt(0)
    if "is_clonal_input_candidate" in data.columns:
        base_eligible &= as_bool(data["is_clonal_input_candidate"])
    eligible_ascn = base_eligible & ascn_overlap
    data["segment_mean_match"] = segment_mean_match
    data["ascn_overlap"] = ascn_overlap
    data["eligible_ascn"] = eligible_ascn
    data["has_counts"] = has_counts

    grouped = data.groupby("sample_barcode", dropna=False).agg(
        n_mutations_with_segment_mean_cn=("segment_mean_match", "sum"),
        n_mutations_with_allele_specific_cn=("ascn_overlap", "sum"),
        n_eligible_mutations_with_allele_specific_cn=("eligible_ascn", "sum"),
        n_mutations_with_ref_alt_counts=("has_counts", "sum"),
    ).reset_index()
    metadata = pilot[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates().copy()
    metadata["sample_barcode"] = metadata["sample_barcode"].map(normalize_sample_barcode)
    readiness = metadata.merge(grouped, on="sample_barcode", how="left")
    count_columns = [
        "n_mutations_with_segment_mean_cn",
        "n_mutations_with_allele_specific_cn",
        "n_eligible_mutations_with_allele_specific_cn",
        "n_mutations_with_ref_alt_counts",
    ]
    readiness[count_columns] = readiness[count_columns].fillna(0).astype(int)
    denominator = readiness["n_mutations_with_segment_mean_cn"]
    readiness["pct_mutations_with_allele_specific_cn"] = np.where(
        denominator > 0,
        100.0 * readiness["n_mutations_with_allele_specific_cn"] / denominator,
        0.0,
    )

    purity = purity_ploidy[["sample_barcode", "purity"]].drop_duplicates().copy()
    purity["sample_barcode"] = purity["sample_barcode"].map(normalize_sample_barcode)
    purity["purity"] = numeric(purity["purity"])
    readiness = readiness.merge(purity, on="sample_barcode", how="left")
    readiness["has_purity"] = readiness["purity"].notna()
    readiness["has_ref_alt_counts"] = readiness["n_mutations_with_ref_alt_counts"] > 0

    segment_flags = segments.copy()
    if segment_flags.empty:
        segment_by_sample = pd.DataFrame(
            columns=["sample_barcode", "has_major_cn", "has_minor_cn"]
        )
    else:
        segment_by_sample = segment_flags.groupby("sample_barcode", as_index=False).agg(
            has_major_cn=("major_cn", lambda values: values.notna().any()),
            has_minor_cn=("minor_cn", lambda values: values.notna().any()),
        )
    readiness = readiness.merge(segment_by_sample, on="sample_barcode", how="left")
    readiness["has_major_cn"] = readiness["has_major_cn"].eq(True)
    readiness["has_minor_cn"] = readiness["has_minor_cn"].eq(True)
    coverage_fraction = readiness["pct_mutations_with_allele_specific_cn"] / 100.0
    readiness["can_build_pyclone_vi_input"] = (
        readiness["has_major_cn"]
        & readiness["has_minor_cn"]
        & readiness["has_purity"]
        & readiness["has_ref_alt_counts"]
        & (
            readiness["n_eligible_mutations_with_allele_specific_cn"]
            >= int(config["min_pyclone_eligible_mutations"])
        )
        & (coverage_fraction >= float(config["min_pyclone_mutation_coverage_fraction"]))
    )

    limited_samples = set(limited_vaf_pilot["sample_barcode"].map(normalize_sample_barcode))
    complexity_lookup = complexity[["sample_barcode", "n_input_mutations", "n_clusters"]].drop_duplicates().copy()
    complexity_lookup["sample_barcode"] = complexity_lookup["sample_barcode"].map(normalize_sample_barcode)
    readiness = readiness.merge(complexity_lookup, on="sample_barcode", how="left")
    readiness["in_limited_vaf_pilot"] = readiness["sample_barcode"].isin(limited_samples)
    readiness["can_build_phylowgs_input"] = (
        readiness["can_build_pyclone_vi_input"]
        & readiness["in_limited_vaf_pilot"]
        & (
            readiness["n_eligible_mutations_with_allele_specific_cn"]
            >= int(config["min_phylowgs_eligible_mutations"])
        )
        & (coverage_fraction >= float(config["min_phylowgs_mutation_coverage_fraction"]))
        & readiness["n_clusters"].ge(int(config["min_phylowgs_expected_clusters"]))
    )

    reasons: list[str] = []
    for row in readiness.itertuples(index=False):
        sample_reasons: list[str] = []
        if not bool(row.has_major_cn):
            sample_reasons.append("major_cn_missing")
        if not bool(row.has_minor_cn):
            sample_reasons.append("minor_cn_missing")
        if not bool(row.has_purity):
            sample_reasons.append("purity_missing")
        if not bool(row.has_ref_alt_counts):
            sample_reasons.append("ref_alt_counts_missing")
        if float(row.pct_mutations_with_allele_specific_cn) < 100 * float(
            config["min_pyclone_mutation_coverage_fraction"]
        ):
            sample_reasons.append(
                f"allele_specific_mutation_coverage_lt_{100 * float(config['min_pyclone_mutation_coverage_fraction']):.0f}pct"
            )
        if int(row.n_eligible_mutations_with_allele_specific_cn) < int(
            config["min_pyclone_eligible_mutations"]
        ):
            sample_reasons.append(
                f"eligible_allele_specific_mutations_lt_{int(config['min_pyclone_eligible_mutations'])}"
            )
        if bool(row.can_build_pyclone_vi_input) and not bool(row.can_build_phylowgs_input):
            if not bool(row.in_limited_vaf_pilot):
                sample_reasons.append("phylowgs_not_prioritized_outside_limited_vaf_pilot")
            if int(row.n_eligible_mutations_with_allele_specific_cn) < int(
                config["min_phylowgs_eligible_mutations"]
            ):
                sample_reasons.append(
                    f"phylowgs_eligible_mutations_lt_{int(config['min_phylowgs_eligible_mutations'])}"
                )
            if float(row.pct_mutations_with_allele_specific_cn) < 100 * float(
                config["min_phylowgs_mutation_coverage_fraction"]
            ):
                sample_reasons.append(
                    f"phylowgs_cn_coverage_lt_{100 * float(config['min_phylowgs_mutation_coverage_fraction']):.0f}pct"
                )
            if pd.isna(row.n_clusters) or float(row.n_clusters) < int(config["min_phylowgs_expected_clusters"]):
                sample_reasons.append("phylowgs_expected_cluster_evidence_missing_or_insufficient")
        if bool(row.can_build_phylowgs_input):
            sample_reasons.append("ready_for_input_preparation_only_single_bulk_phylogeny_remains_cautious")
        elif bool(row.can_build_pyclone_vi_input) and not sample_reasons:
            sample_reasons.append("ready_for_pyclone_vi_input_preparation_only")
        reasons.append(";".join(sample_reasons) if sample_reasons else "ready_for_input_preparation_review")
    readiness["reason_not_ready"] = reasons
    return readiness[MUTATION_READINESS_COLUMNS]


def build_availability_table(
    pilot: pd.DataFrame,
    limited_vaf_pilot: pd.DataFrame,
    segment_mean_summary: pd.DataFrame,
    segments: pd.DataFrame,
    mutation_readiness: pd.DataFrame,
) -> pd.DataFrame:
    availability = pilot[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates().copy()
    availability["sample_barcode"] = availability["sample_barcode"].map(normalize_sample_barcode)
    availability["in_pilot_cohort"] = True
    limited_samples = set(limited_vaf_pilot["sample_barcode"].map(normalize_sample_barcode))
    availability["in_limited_vaf_pilot"] = availability["sample_barcode"].isin(limited_samples)
    segment_mean = segment_mean_summary[["sample_barcode", "has_segment_level_cn"]].drop_duplicates().copy()
    segment_mean["sample_barcode"] = segment_mean["sample_barcode"].map(normalize_sample_barcode)
    segment_mean["has_segment_level_cn"] = as_bool(segment_mean["has_segment_level_cn"])
    availability = availability.merge(segment_mean, on="sample_barcode", how="left")
    availability["has_segment_mean_cn"] = availability.pop("has_segment_level_cn").fillna(False).astype(bool)

    if segments.empty:
        segment_summary = pd.DataFrame(
            columns=[
                "sample_barcode",
                "allele_specific_cn_source",
                "n_allele_specific_segments",
                "has_total_cn",
                "has_major_cn",
                "has_minor_cn",
                "has_loh_status",
            ]
        )
    else:
        segment_data = segments.copy()
        segment_data["valid_ascn"] = segment_data["major_cn"].notna() & segment_data["minor_cn"].notna()
        segment_data["valid_loh"] = segment_data["loh_status"].notna() & ~segment_data[
            "loh_status"
        ].astype(str).str.lower().isin({"na", "nan", "unknown"})
        segment_summary = segment_data.groupby("sample_barcode", as_index=False).agg(
            allele_specific_cn_source=(
                "source_tool",
                lambda values: ";".join(sorted(set(values.dropna().astype(str)))),
            ),
            n_allele_specific_segments=("valid_ascn", "sum"),
            has_total_cn=("total_cn", lambda values: values.notna().any()),
            has_major_cn=("major_cn", lambda values: values.notna().any()),
            has_minor_cn=("minor_cn", lambda values: values.notna().any()),
            has_loh_status=("valid_loh", "any"),
        )
    availability = availability.merge(segment_summary, on="sample_barcode", how="left")
    availability["allele_specific_cn_source"] = availability["allele_specific_cn_source"].fillna(
        "none_found_locally"
    )
    availability["n_allele_specific_segments"] = pd.to_numeric(
        availability["n_allele_specific_segments"], errors="coerce"
    ).fillna(0).astype(int)
    for column in ["has_total_cn", "has_major_cn", "has_minor_cn", "has_loh_status"]:
        availability[column] = availability[column].eq(True)
    availability["has_allele_specific_cn"] = availability["n_allele_specific_segments"] > 0
    availability = availability.merge(
        mutation_readiness[["sample_barcode", "can_build_pyclone_vi_input", "reason_not_ready"]],
        on="sample_barcode",
        how="left",
    )
    availability["can_prepare_cn_aware_pyclone_input"] = availability.pop(
        "can_build_pyclone_vi_input"
    ).fillna(False).astype(bool)

    actions: list[str] = []
    notes: list[str] = []
    for row in availability.itertuples(index=False):
        if bool(row.can_prepare_cn_aware_pyclone_input):
            actions.append(
                "review_standardized_segments_and_prepare_validated_pyclone_vi_input_without_running_in_this_step"
            )
        elif bool(row.has_allele_specific_cn):
            actions.append("review_or_regenerate_allele_specific_cn_and_improve_mutation_overlap")
        else:
            actions.append("obtain_tumor_normal_inputs_and_run_one_consistent_allele_specific_cn_method")
        notes.append(
            "segment_mean_is_not_allele_specific_integer_copy_number;single_bulk_sample_limits_branching_inference"
        )
    availability["recommended_next_action"] = actions
    availability["notes"] = notes
    return availability[AVAILABILITY_COLUMNS]


def discover_external_assets(raw_data_root: Path) -> list[Path]:
    if not raw_data_root.exists():
        return []
    assets: list[Path] = []
    for path in raw_data_root.rglob("*"):
        if not path.is_file():
            continue
        lower = path.name.lower()
        if (
            lower.endswith((".bam", ".cram", ".vcf", ".vcf.gz", ".seqz", ".seqz.gz", ".pileup"))
            or "snp_pileup" in lower
            or "snp-pileup" in lower
        ):
            assets.append(path)
    return sorted(assets)


def sample_assets(sample_barcode: str, patient_barcode: str, assets: list[Path]) -> list[Path]:
    sample = sample_barcode.upper()
    patient = patient_barcode.upper()
    return [path for path in assets if sample in str(path).upper() or patient in str(path).upper()]


def normal_bam_like(path: Path) -> bool:
    text = str(path).lower()
    if path.suffix.lower() not in {".bam", ".cram"}:
        return False
    if any(token in text for token in ["normal", "blood", "germline"]):
        return True
    match = re.search(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2})[A-Z]?", str(path).upper())
    return bool(match and 10 <= int(match.group(1)) <= 14)


def build_external_manifest(
    availability: pd.DataFrame,
    mutation_readiness: pd.DataFrame,
    assets: list[Path],
    config: dict[str, Any] = CONFIG,
) -> pd.DataFrame:
    readiness_lookup = mutation_readiness.set_index("sample_barcode")
    rows: list[dict[str, Any]] = []
    for row in availability.itertuples(index=False):
        local_assets = sample_assets(str(row.sample_barcode), str(row.patient_barcode), assets)
        tumor_bams = [path for path in local_assets if path.suffix.lower() in {".bam", ".cram"} and not normal_bam_like(path)]
        normal_bams = [path for path in local_assets if normal_bam_like(path)]
        vcfs = [path for path in local_assets if path.name.lower().endswith((".vcf", ".vcf.gz"))]
        pileups = [path for path in local_assets if "pileup" in path.name.lower() or ".seqz" in path.name.lower()]
        current = ["segment_mean_cn", "somatic_mutation_ref_alt_counts"]
        readiness_row = readiness_lookup.loc[row.sample_barcode]
        if bool(readiness_row["has_purity"]):
            current.append("purity_ploidy")
        if bool(row.has_allele_specific_cn):
            current.append("local_allele_specific_cn_segments")
        if tumor_bams:
            current.append("local_tumor_bam_or_cram")
        if normal_bams:
            current.append("local_matched_normal_bam_or_cram")
        if vcfs:
            current.append("local_vcf")
        if pileups:
            current.append("local_snp_pileup_or_seqz")

        if bool(row.can_prepare_cn_aware_pyclone_input):
            recommended_tool = "existing_local_ASCN_QC"
            missing: list[str] = []
            command = "not_required_existing_allele_specific_cn_requires_manual_qc_and_input_preparation"
        else:
            recommended_tool = str(config["recommended_external_tool"])
            missing = []
            if not tumor_bams and not pileups:
                missing.append("tumor_bam_or_precomputed_tumor_snp_counts")
            if not normal_bams and not pileups:
                missing.append("matched_normal_bam_or_precomputed_normal_snp_counts")
            if not pileups:
                missing.append("snp_pileup_or_equivalent_allele_count_input")
            command = (
                "<snp-pileup> <common_snp_vcf.gz> <sample_counts.csv.gz> <normal.bam> <tumor.bam>; "
                "Rscript <facets_wrapper.R> <sample_counts.csv.gz> <output_dir>"
            )
        rows.append(
            {
                "sample_barcode": row.sample_barcode,
                "patient_barcode": row.patient_barcode,
                "project_id": row.project_id,
                "project_code": row.project_code,
                "recommended_tool": recommended_tool,
                "tumor_bam_required": not bool(row.can_prepare_cn_aware_pyclone_input),
                "normal_bam_required": not bool(row.can_prepare_cn_aware_pyclone_input),
                "tumor_vcf_required": False,
                "normal_vcf_required": False,
                "snp_pileup_required": not bool(row.can_prepare_cn_aware_pyclone_input),
                "matched_normal_available_unknown": not bool(normal_bams),
                "current_available_inputs": ";".join(current),
                "missing_required_inputs": ";".join(missing) if missing else "none",
                "recommended_command_template": command,
                "notes": (
                    "placeholder_command_requires_tool_version_reference_build_common_snp_resource_and_manual_qc;"
                    "do_not_assume_bams_exist;do_not_run_pyclone_vi_or_phylowgs_in_this_step"
                ),
            }
        )
    return pd.DataFrame(rows, columns=MANIFEST_COLUMNS)


def missing_input_counts(manifest: pd.DataFrame) -> pd.Series:
    counts: dict[str, int] = {}
    for value in manifest["missing_required_inputs"].fillna("none").astype(str):
        for item in value.split(";"):
            item = item.strip()
            if item and item != "none":
                counts[item] = counts.get(item, 0) + 1
    return pd.Series(counts, dtype=int).sort_values(ascending=False)


def build_qc_summary(
    availability: pd.DataFrame,
    mutation_readiness: pd.DataFrame,
    standardized_segments: pd.DataFrame,
    candidate_files: list[Path],
    parsed_files: list[Path],
    parse_warnings: list[str],
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    missing_counts = missing_input_counts(manifest)
    missing_text = ";".join(f"{name}:{int(count)}" for name, count in missing_counts.items()) or "none"
    rows = [
        ("global", "pilot_samples_evaluated", len(availability)),
        ("global", "limited_vaf_pilot_samples_evaluated", int(availability["in_limited_vaf_pilot"].sum())),
        ("global", "candidate_local_allele_specific_cn_files_found", len(candidate_files)),
        ("global", "candidate_local_files_parsed_successfully", len(parsed_files)),
        ("global", "standardized_allele_specific_segments", len(standardized_segments)),
        ("global", "samples_with_existing_allele_specific_cn", int(availability["has_allele_specific_cn"].sum())),
        ("global", "samples_with_total_cn", int(availability["has_total_cn"].sum())),
        ("global", "samples_with_major_cn", int(availability["has_major_cn"].sum())),
        ("global", "samples_with_minor_cn", int(availability["has_minor_cn"].sum())),
        ("global", "samples_with_loh_status", int(availability["has_loh_status"].sum())),
        ("global", "samples_ready_for_pyclone_vi", int(mutation_readiness["can_build_pyclone_vi_input"].sum())),
        ("global", "samples_ready_for_phylowgs", int(mutation_readiness["can_build_phylowgs_input"].sum())),
        (
            "global",
            "samples_requiring_external_allele_specific_cn_inference",
            int((~availability["can_prepare_cn_aware_pyclone_input"]).sum()),
        ),
        ("global", "likely_missing_input_types", missing_text),
        ("global", "parse_warnings", ";".join(parse_warnings) if parse_warnings else "none"),
        (
            "global",
            "limitations",
            "segment_mean_is_not_allele_specific_integer_copy_number;major_minor_cn_not_inferred_from_segment_mean;"
            "no_pyclone_vi_or_phylowgs_run;single_bulk_samples_limit_branching_inference",
        ),
    ]
    for project_code, group in availability.groupby("project_code", sort=True):
        project_readiness = mutation_readiness[
            mutation_readiness["sample_barcode"].isin(group["sample_barcode"])
        ]
        rows.extend(
            [
                ("project", f"{project_code}_pilot_samples", len(group)),
                ("project", f"{project_code}_samples_with_allele_specific_cn", int(group["has_allele_specific_cn"].sum())),
                ("project", f"{project_code}_samples_ready_for_pyclone_vi", int(project_readiness["can_build_pyclone_vi_input"].sum())),
            ]
        )
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def markdown_table(data: pd.DataFrame, columns: list[tuple[str, str]]) -> str:
    headers = [label for _, label in columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in data.iterrows():
        values: list[str] = []
        for column, _ in columns:
            value = row[column]
            if isinstance(value, (float, np.floating)) and np.isfinite(value):
                values.append(f"{float(value):.3f}")
            else:
                values.append(str(value).replace("|", "\\|").replace("\n", " "))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_planning_report(
    availability: pd.DataFrame,
    mutation_readiness: pd.DataFrame,
    manifest: pd.DataFrame,
    candidate_files: list[Path],
    parsed_files: list[Path],
    parse_warnings: list[str],
    config: dict[str, Any] = CONFIG,
) -> str:
    n_pilot = len(availability)
    n_limited = int(availability["in_limited_vaf_pilot"].sum())
    n_existing = int(availability["has_allele_specific_cn"].sum())
    n_major_minor = int((availability["has_major_cn"] & availability["has_minor_cn"]).sum())
    n_pyclone = int(mutation_readiness["can_build_pyclone_vi_input"].sum())
    n_phylowgs = int(mutation_readiness["can_build_phylowgs_input"].sum())
    project_summary = availability.groupby("project_code", as_index=False).agg(
        pilot_samples=("sample_barcode", "size"),
        limited_vaf_samples=("in_limited_vaf_pilot", "sum"),
        samples_with_ascn=("has_allele_specific_cn", "sum"),
        pyclone_ready=("can_prepare_cn_aware_pyclone_input", "sum"),
    )
    project_table = markdown_table(
        project_summary,
        [
            ("project_code", "Project"),
            ("pilot_samples", "Pilot samples"),
            ("limited_vaf_samples", "Limited-VAF pilot"),
            ("samples_with_ascn", "Existing ASCN"),
            ("pyclone_ready", "PyClone-VI ready"),
        ],
    )
    missing_counts = missing_input_counts(manifest)
    missing_table = markdown_table(
        missing_counts.rename_axis("missing_input").reset_index(name="samples"),
        [("missing_input", "Likely missing input"), ("samples", "Samples")],
    ) if not missing_counts.empty else "No external input gaps were recorded."
    file_lines = (
        "\n".join(f"- `{path}`" for path in candidate_files)
        if candidate_files
        else "- No candidate local allele-specific CN files were found."
    )
    warning_lines = (
        "\n".join(f"- {warning}" for warning in parse_warnings)
        if parse_warnings
        else "- No parser warnings."
    )
    return f"""# Level 3 Allele-Specific CN Readiness and PyClone-VI/PhyloWGS Plan

Generated by `scripts/12b_allele_specific_cn_readiness.py`. This step assesses and plans allele-specific copy-number acquisition; it does not run PyClone-VI or PhyloWGS.

## A. Current Status

- Pilot samples evaluated: {n_pilot}
- Limited-VAF pilot samples evaluated: {n_limited}
- Candidate local allele-specific CN files found: {len(candidate_files)}
- Candidate files parsed with usable CN rows: {len(parsed_files)}
- Samples with existing allele-specific major/minor CN: {n_existing}
- Samples with both major and minor CN: {n_major_minor}
- Samples meeting the current PyClone-VI input-readiness rules: {n_pyclone}
- Samples meeting the current cautious PhyloWGS input-readiness rules: {n_phylowgs}

The project currently has sample-level segment-mean/log2-like CN for the pilot cohort, mutation ref/alt counts, purity, ploidy, and limited VAF clustering results. However, `segment_mean is not allele-specific integer copy number`. It does not directly provide `total_cn`, `major_cn`, `minor_cn`, or a validated LOH state. The existing segment-mean mutation annotation must therefore remain separate from true copy-number-aware clonal input.

{project_table}

Candidate local files:

{file_lines}

Parser notes:

{warning_lines}

## B. Why Allele-Specific CN Is Required

Observed VAF depends on tumor purity, local total copy number, the multiplicity of the mutated allele, and contamination by normal cells. Total CN describes the number of tumor copies at a locus. Major and minor CN describe the allele-specific split, which distinguishes balanced gains from copy-neutral or deletion-associated LOH. LOH status helps determine whether one parental allele has been lost. Without this information, the same VAF can correspond to different cancer-cell fractions and mutation multiplicities.

Segment means are useful for locating broad gain-like, loss-like, or near-neutral regions, but they are continuous assay signals rather than validated integer allele states. The script never derives major/minor CN from segment mean.

## C. Tool Options

### FACETS

FACETS is a practical primary option for paired tumor/normal sequencing. It generally requires common-SNP pileups or equivalent allele counts generated from tumor and matched-normal BAM/CRAM files, followed by segmentation and integer total/minor CN estimation. A placeholder workflow is:

`<snp-pileup> <common_snp_vcf.gz> <sample_counts.csv.gz> <normal.bam> <tumor.bam>; Rscript <facets_wrapper.R> <sample_counts.csv.gz> <output_dir>`

Exact flags, genome build, SNP resource, normal/tumor ordering, and QC thresholds must be fixed for one validated tool version before execution.

### Sequenza

Sequenza can estimate purity, ploidy, and allele-specific CN from paired tumor/normal depth and B-allele-frequency-style inputs. A placeholder preparation path is `sequenza-utils bam2seqz` followed by the Sequenza R workflow. Reference FASTA, GC correction resources, matched-normal BAM, and tumor BAM are required unless compatible seqz inputs already exist.

### ASCAT

ASCAT uses logR and BAF measurements from SNP arrays or sequencing-derived allele counts. It is suitable when harmonized SNP-array or sequencing inputs are available, but preprocessing and platform assumptions must be consistent across the pilot.

### ABSOLUTE-Style Outputs

Existing ABSOLUTE-style purity/ploidy summaries are already available, but the current project file does not contain sample-level allele-specific segment states. If a public or locally generated ABSOLUTE-style segment table with integer total and homolog-specific CN is obtained, it can be standardized here after source and coordinate-build review.

### Existing Public or Local Resources

The supported local directories are searched recursively for FACETS, Sequenza, ASCAT, ABSOLUTE-style, allele-specific, and major/minor segment tables. Any discovered resource must be checked for genome build, sample/aliquot identity, integer CN semantics, and whether values are truly allele-specific rather than segment means.

## D. Recommended Path for This Project

Recommended next action: prioritize the 30 limited-VAF pilot samples, obtain paired tumor/normal BAMs or validated precomputed SNP allele-count inputs, and run one allele-specific CN method consistently across that subset.

1. Obtain tumor/normal BAMs or appropriate allele-specific CN inputs for the 30 limited-VAF pilot samples first.
2. Run one allele-specific CN method consistently, with FACETS as the default planning option for paired WES if suitable inputs can be obtained.
3. Standardize reviewed outputs into `data/processed/allele_specific_cn/allele_specific_segments_by_sample.tsv.gz`.
4. Re-annotate eligible mutations with integer `total_cn`, `major_cn`, `minor_cn`, and LOH status.
5. Prepare and validate PyClone-VI input only for samples with at least {int(config['min_pyclone_eligible_mutations'])} eligible mutations and at least {100 * float(config['min_pyclone_mutation_coverage_fraction']):.0f}% allele-specific CN overlap.
6. Run PyClone-VI only after the standardized inputs pass per-sample QC.
7. Consider PhyloWGS only on a high-confidence subset with at least {int(config['min_phylowgs_eligible_mutations'])} eligible mutations, at least {100 * float(config['min_phylowgs_mutation_coverage_fraction']):.0f}% CN overlap, interpretable cluster complexity, and explicit single-bulk-sample caveats.

Likely missing external inputs:

{missing_table}

The external inference manifest records requirements conservatively. It does not assume BAMs, matched normals, SNP pileups, or VCFs are available merely because TCGA-derived mutation and segment tables exist locally.

## E. Explicit Caveat

Do not run PyClone-VI or PhyloWGS until allele-specific CN is available and mutation overlap, purity, ref/alt counts, coordinate build, and integer-CN semantics have been reviewed. A copy-neutral exploratory analysis may continue only if it remains clearly labeled as limited VAF-based clustering rather than copy-number-aware clonal inference or a definitive branching phylogeny.

Single bulk TCGA samples remain a fundamental limitation even after allele-specific CN is added. PyClone-VI may support more defensible clonal clustering, but PhyloWGS branching should be restricted to carefully reviewed high-confidence cases and interpreted cautiously.
"""


def plot_readiness_by_project(availability: pd.DataFrame, mutation_readiness: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    merged = availability.merge(
        mutation_readiness[["sample_barcode", "can_build_pyclone_vi_input", "can_build_phylowgs_input"]],
        on="sample_barcode",
        how="left",
    )
    summary = merged.groupby("project_code", as_index=False).agg(
        pilot_samples=("sample_barcode", "size"),
        existing_ascn=("has_allele_specific_cn", "sum"),
        pyclone_ready=("can_build_pyclone_vi_input", "sum"),
        phylowgs_ready=("can_build_phylowgs_input", "sum"),
    )
    x = np.arange(len(summary))
    width = 0.20
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.bar(x - 1.5 * width, summary["pilot_samples"], width, label="Pilot samples", color="#4B6A88")
    ax.bar(x - 0.5 * width, summary["existing_ascn"], width, label="Existing ASCN", color="#4D7A55")
    ax.bar(x + 0.5 * width, summary["pyclone_ready"], width, label="PyClone-VI ready", color="#B06A3C")
    ax.bar(x + 1.5 * width, summary["phylowgs_ready"], width, label="PhyloWGS ready", color="#7A5A8C")
    ax.set_xticks(x, summary["project_code"], rotation=45, ha="right")
    ax.set_ylabel("Samples")
    ax.set_title("Level 3 allele-specific CN readiness by project")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    fig.text(
        0.5,
        0.01,
        "ASCN means integer major/minor copy number; segment-mean CN alone does not count.",
        ha="center",
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_missing_inputs(manifest: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    counts = missing_input_counts(manifest)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    if counts.empty:
        ax.text(0.5, 0.5, "No missing external ASCN inputs recorded", ha="center", va="center", fontsize=13)
        ax.set_axis_off()
    else:
        plotted = counts.sort_values(ascending=True)
        labels = [label.replace("_", " ") for label in plotted.index]
        bars = ax.barh(labels, plotted.values, color="#B04A5A")
        ax.bar_label(bars, labels=[str(int(value)) for value in plotted.values], padding=3, fontsize=8)
        ax.set_xlabel("Pilot samples")
        ax.set_title("Likely missing inputs for external allele-specific CN inference")
        ax.set_xlim(0, max(float(plotted.max()) * 1.15, 1))
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    fig.text(
        0.5,
        0.01,
        "Availability is based on local file discovery only; remote controlled-access assets were not queried.",
        ha="center",
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", path)


def run_allele_specific_cn_readiness(
    paths: AlleleSpecificCnPaths,
    config: dict[str, Any] = CONFIG,
) -> dict[str, Any]:
    pilot = read_tsv(paths.pilot_with_segments, "Level 3 pilot cohort with segments")
    limited_vaf_pilot = read_tsv(paths.limited_vaf_pilot, "limited-VAF pilot samples")
    mutations = read_tsv(paths.mutations_with_local_cn, "mutations with local segment CN")
    segment_mean_cn = read_tsv(paths.segment_mean_cn, "segment-mean CN table")
    segment_mean_summary = read_tsv(paths.segment_mean_summary, "segment-mean CN summary")
    purity_ploidy = read_tsv(paths.purity_ploidy, "purity/ploidy table")
    clonal_readiness = read_tsv(paths.clonal_input_readiness, "Level 3 clonal input readiness")
    complexity = read_tsv(paths.limited_vaf_complexity, "limited-VAF clonal complexity")

    require_columns(
        pilot,
        ["sample_barcode", "patient_barcode", "project_id", "project_code", "has_segment_level_cn"],
        "Level 3 pilot cohort with segments",
    )
    require_columns(limited_vaf_pilot, ["sample_barcode", "project_code"], "limited-VAF pilot samples")
    require_columns(
        mutations,
        [
            "sample_barcode",
            "chromosome",
            "position",
            "t_ref_count",
            "t_alt_count",
            "total_depth",
            "observed_vaf",
            "local_segment_mean",
            "local_cn_match_status",
        ],
        "mutations with local segment CN",
    )
    require_columns(segment_mean_cn, ["sample_barcode", "chromosome", "start", "end", "segment_mean"], "segment-mean CN table")
    require_columns(
        segment_mean_summary,
        ["sample_barcode", "has_segment_level_cn"],
        "segment-mean CN summary",
    )
    require_columns(purity_ploidy, ["sample_barcode", "purity", "ploidy"], "purity/ploidy table")
    require_columns(
        clonal_readiness,
        ["sample_barcode", "eligible_for_copy_number_aware_clustering", "reason_not_copy_number_aware"],
        "Level 3 clonal input readiness",
    )
    require_columns(
        complexity,
        ["sample_barcode", "n_input_mutations", "n_clusters", "interpretation_class"],
        "limited-VAF clonal complexity",
    )

    pilot_samples = set(pilot["sample_barcode"].map(normalize_sample_barcode))
    if len(pilot_samples) != len(pilot):
        raise ValueError("Pilot cohort contains duplicate normalized sample barcodes")
    if not set(limited_vaf_pilot["sample_barcode"].map(normalize_sample_barcode)).issubset(pilot_samples):
        raise ValueError("Limited-VAF pilot contains samples outside the Level 3 pilot cohort")
    if not set(mutations["sample_barcode"].map(normalize_sample_barcode)).issubset(pilot_samples):
        raise ValueError("Mutation table contains samples outside the Level 3 pilot cohort")
    if any(column in segment_mean_cn.columns for column in ["major_cn", "minor_cn"]):
        logging.warning(
            "Segment-mean input unexpectedly contains major/minor columns; it remains excluded from ASCN discovery"
        )

    standardized, candidate_files, parsed_files, parse_warnings = ingest_local_allele_specific_cn(
        paths, pilot, config
    )
    mutation_readiness = build_mutation_readiness(
        pilot,
        limited_vaf_pilot,
        mutations,
        standardized,
        purity_ploidy,
        complexity,
        config,
    )
    availability = build_availability_table(
        pilot,
        limited_vaf_pilot,
        segment_mean_summary,
        standardized,
        mutation_readiness,
    )
    assets = discover_external_assets(paths.raw_data_root)
    manifest = build_external_manifest(availability, mutation_readiness, assets, config)
    qc = build_qc_summary(
        availability,
        mutation_readiness,
        standardized,
        candidate_files,
        parsed_files,
        parse_warnings,
        manifest,
    )
    report = build_planning_report(
        availability,
        mutation_readiness,
        manifest,
        candidate_files,
        parsed_files,
        parse_warnings,
        config,
    )

    write_tsv_gz(standardized, paths.standardized_segments)
    write_tsv(availability, paths.availability_by_sample)
    write_tsv(mutation_readiness, paths.mutation_readiness)
    write_tsv(manifest, paths.external_manifest)
    write_tsv(qc, paths.qc_summary)
    write_text(report, paths.planning_report)
    plot_readiness_by_project(availability, mutation_readiness, paths.readiness_figure)
    plot_missing_inputs(manifest, paths.missing_inputs_figure)
    return {
        "standardized_segments": standardized,
        "availability": availability,
        "mutation_readiness": mutation_readiness,
        "external_manifest": manifest,
        "qc": qc,
        "candidate_files": candidate_files,
        "parsed_files": parsed_files,
        "parse_warnings": parse_warnings,
        "external_assets": assets,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Assess local allele-specific integer CN availability and plan external inference before "
            "PyClone-VI or PhyloWGS."
        )
    )
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    paths = default_paths(root)
    outputs = run_allele_specific_cn_readiness(paths)
    logging.info(
        "Allele-specific CN readiness complete: pilot=%d files_found=%d samples_with_ascn=%d "
        "pyclone_ready=%d phylowgs_ready=%d",
        len(outputs["availability"]),
        len(outputs["candidate_files"]),
        int(outputs["availability"]["has_allele_specific_cn"].sum()),
        int(outputs["mutation_readiness"]["can_build_pyclone_vi_input"].sum()),
        int(outputs["mutation_readiness"]["can_build_phylowgs_input"].sum()),
    )


if __name__ == "__main__":
    main()

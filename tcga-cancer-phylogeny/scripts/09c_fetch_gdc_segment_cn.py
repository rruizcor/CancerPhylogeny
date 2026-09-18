#!/usr/bin/env python

from __future__ import annotations

import argparse
import gzip
import json
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
import pandas as pd
import requests

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file


GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"

SEGMENT_DATA_TYPES = [
    "Masked Copy Number Segment",
    "Copy Number Segment",
    "Filtered Copy Number Segment",
]

DOWNLOAD_MODE = "pilot_full"  # none, pilot_capped, pilot_full, candidates_capped
MAX_DOWNLOADS = 10
PREFER_SAMPLE_LEVEL_MATCHES = True
ALLOW_PATIENT_LEVEL_MATCHES = False
OVERWRITE_EXISTING = False

CONFIG = {
    "download_mode": DOWNLOAD_MODE,
    "request_timeout_seconds": 90,
    "download_timeout_seconds": 30,
    "query_page_size": 5000,
    "max_query_records": 50000,
    "max_downloads": MAX_DOWNLOADS,
    "max_download_files": MAX_DOWNLOADS,
    "prefer_sample_level_matches": PREFER_SAMPLE_LEVEL_MATCHES,
    "allow_patient_level_matches": ALLOW_PATIENT_LEVEL_MATCHES,
    "overwrite_existing": OVERWRITE_EXISTING,
    "gain_like_threshold": 0.2,
    "loss_like_threshold": -0.2,
    "min_segments_for_usable_cn": 1,
}

MANIFEST_COLUMNS = [
    "gdc_file_id",
    "file_name",
    "data_category",
    "data_type",
    "experimental_strategy",
    "access",
    "file_size",
    "project_id",
    "case_submitter_ids",
    "sample_submitter_ids",
    "aliquot_barcodes",
    "matched_sample_barcode",
    "matched_patient_barcode",
    "match_level",
    "download_status",
    "local_path",
    "parse_status",
    "notes",
]

BEST_FILE_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "selected_gdc_file_id",
    "selected_file_name",
    "selected_data_type",
    "match_level",
    "file_size",
    "n_candidate_files",
    "selection_reason",
    "selection_rank",
    "notes",
]

DOWNLOAD_LOG_COLUMNS = [
    "sample_barcode",
    "gdc_file_id",
    "file_name",
    "download_attempted",
    "download_status",
    "local_path",
    "error_message",
    "file_size",
]

SEGMENT_COLUMNS = [
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
    "gdc_file_id",
    "match_level",
    "notes",
]

SAMPLE_COLUMNS = [
    "sample_barcode",
    "Sample",
    "sample",
    "Sample_ID",
    "ID",
    "Tumor_Sample_Barcode",
    "aliquot_barcode",
    "Aliquot",
    "GDC_Aliquot",
]
CHROMOSOME_COLUMNS = ["chromosome", "Chromosome", "chrom", "chr"]
START_COLUMNS = ["start", "Start", "loc.start", "chromStart", "Start_Position"]
END_COLUMNS = ["end", "End", "loc.end", "chromEnd", "End_Position"]
NUM_PROBES_COLUMNS = ["num_probes", "Num_Probes", "probes", "num.mark", "Num.Mark", "Num_Probes_Segment"]
SEGMENT_MEAN_COLUMNS = [
    "segment_mean",
    "Segment_Mean",
    "seg.mean",
    "log2_copy_ratio",
    "Log2_Copy_Ratio",
    "copy_ratio",
    "SegmentMean",
    "mean",
]


@dataclass
class GdcSegmentPaths:
    candidate_samples: Path
    pilot_cohort: Path
    manifest: Path
    best_file_per_pilot_sample: Path
    download_log: Path
    raw_download_dir: Path
    processed_segments_gz: Path
    segment_summary: Path
    candidates_with_segments: Path
    pilot_with_segments: Path
    qc_summary: Path
    availability_by_project: Path
    availability_figure: Path
    pilot_coverage_figure: Path
    upgrade_by_project_figure: Path
    count_distribution_figure: Path
    upgrade_status_figure: Path


def default_paths(root: Path) -> GdcSegmentPaths:
    return GdcSegmentPaths(
        candidate_samples=project_path("results", "tables", "level3_candidate_samples.tsv", root=root),
        pilot_cohort=project_path("results", "tables", "level3_pilot_cohort.tsv", root=root),
        manifest=project_path("results", "tables", "gdc_segment_cn_manifest_pilot.tsv", root=root),
        best_file_per_pilot_sample=project_path(
            "results", "tables", "gdc_segment_cn_best_file_per_pilot_sample.tsv", root=root
        ),
        download_log=project_path("results", "tables", "gdc_segment_cn_download_log.tsv", root=root),
        raw_download_dir=project_path("data", "raw", "copy_number_segments", "gdc_pilot", root=root),
        processed_segments_gz=project_path("data", "processed", "copy_number", "segments_by_sample.tsv.gz", root=root),
        segment_summary=project_path("data", "processed", "copy_number", "segment_level_cn_summary_by_sample.tsv", root=root),
        candidates_with_segments=project_path("results", "tables", "level3_candidate_samples_with_segments.tsv", root=root),
        pilot_with_segments=project_path("results", "tables", "level3_pilot_cohort_with_segments.tsv", root=root),
        qc_summary=project_path("results", "tables", "level3_segment_cn_qc_summary.tsv", root=root),
        availability_by_project=project_path(
            "results", "tables", "level3_segment_cn_availability_by_project.tsv", root=root
        ),
        availability_figure=project_path("results", "figures", "level3_segment_cn_availability_by_project.pdf", root=root),
        pilot_coverage_figure=project_path("results", "figures", "level3_segment_cn_pilot_coverage.pdf", root=root),
        upgrade_by_project_figure=project_path("results", "figures", "level3_segment_cn_upgrade_by_project.pdf", root=root),
        count_distribution_figure=project_path("results", "figures", "level3_segment_count_distribution.pdf", root=root),
        upgrade_status_figure=project_path("results", "figures", "level3_candidate_upgrade_status.pdf", root=root),
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


def project_code_from_project_id(project_id: str) -> str:
    text = str(project_id)
    return text.replace("TCGA-", "") if text.startswith("TCGA-") else text


def empty_manifest() -> pd.DataFrame:
    return pd.DataFrame(columns=MANIFEST_COLUMNS)


def empty_segments() -> pd.DataFrame:
    return pd.DataFrame(columns=SEGMENT_COLUMNS)


def safe_join(values: list[str] | set[str]) -> str:
    clean = sorted({str(value) for value in values if str(value) and str(value) != "NA"})
    return ";".join(clean) if clean else "NA"


def build_gdc_filters(project_ids: list[str], *, use_file_prefix: bool = True, include_data_category: bool = True) -> dict:
    data_category_field = "files.data_category" if use_file_prefix else "data_category"
    data_type_field = "files.data_type" if use_file_prefix else "data_type"
    access_field = "files.access" if use_file_prefix else "access"
    content: list[dict] = [
        {"op": "in", "content": {"field": "cases.project.program.name", "value": ["TCGA"]}},
        {"op": "in", "content": {"field": "cases.project.project_id", "value": sorted(project_ids)}},
        {"op": "in", "content": {"field": data_type_field, "value": SEGMENT_DATA_TYPES}},
        {"op": "in", "content": {"field": access_field, "value": ["open"]}},
    ]
    if include_data_category:
        content.insert(
            2,
            {
                "op": "in",
                "content": {"field": data_category_field, "value": ["Copy Number Variation"]},
            },
        )
    return {"op": "and", "content": content}


def candidate_filter_sets(project_ids: list[str]) -> list[dict[str, Any]]:
    return [
        {
            "name": "files_prefixed_copy_number_variation",
            "filters": build_gdc_filters(project_ids, use_file_prefix=True, include_data_category=True),
        },
        {
            "name": "unprefixed_copy_number_variation",
            "filters": build_gdc_filters(project_ids, use_file_prefix=False, include_data_category=True),
        },
        {
            "name": "unprefixed_data_type_only",
            "filters": build_gdc_filters(project_ids, use_file_prefix=False, include_data_category=False),
        },
    ]


def gdc_fields() -> str:
    return ",".join(
        [
            "file_id",
            "file_name",
            "data_category",
            "data_type",
            "experimental_strategy",
            "access",
            "file_size",
            "cases.project.project_id",
            "cases.submitter_id",
            "cases.samples.submitter_id",
        ]
    )


def post_gdc_files(filters: dict, offset: int, config: dict = CONFIG) -> dict:
    payload = {
        "filters": filters,
        "fields": gdc_fields(),
        "format": "JSON",
        "size": int(config["query_page_size"]),
        "from": offset,
        "sort": "file_size:asc",
    }
    response = requests.post(GDC_FILES_ENDPOINT, json=payload, timeout=int(config["request_timeout_seconds"]))
    response.raise_for_status()
    return response.json()


def query_gdc_files(project_ids: list[str], config: dict = CONFIG) -> tuple[list[dict], dict, list[str]]:
    warnings: list[str] = []
    attempted_filters = candidate_filter_sets(project_ids)
    last_filters = attempted_filters[-1]["filters"]
    for filter_spec in attempted_filters:
        filters = filter_spec["filters"]
        last_filters = filters
        hits: list[dict] = []
        try:
            total = None
            offset = 0
            while True:
                payload = post_gdc_files(filters, offset, config)
                data = payload.get("data", {})
                page_hits = data.get("hits", [])
                hits.extend(page_hits)
                pagination = data.get("pagination", {})
                total = int(pagination.get("total", len(hits))) if total is None else total
                if len(hits) >= total:
                    break
                offset += int(config["query_page_size"])
                if offset >= int(config["max_query_records"]):
                    warnings.append(f"gdc_query_truncated_at_{config['max_query_records']}_records")
                    break
            logging.info("GDC query '%s' returned %d files", filter_spec["name"], len(hits))
            if hits:
                return hits, filters, warnings
        except requests.RequestException as exc:
            warnings.append(f"gdc_query_failed_{filter_spec['name']}:{exc}")
            logging.warning("GDC query '%s' failed: %s", filter_spec["name"], exc)
    return [], last_filters, warnings


def file_name_barcodes(file_name: str) -> list[str]:
    pattern = r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-[A-Z0-9]{3})?(?:-[A-Z0-9]{3})?(?:-[A-Z0-9]{4})?(?:-[A-Z0-9]{2})?"
    return re.findall(pattern, str(file_name).upper())


def extract_hit_ids(hit: dict) -> tuple[set[str], set[str], set[str]]:
    project_ids: set[str] = set()
    case_ids: set[str] = set()
    sample_ids: set[str] = set()
    for case in hit.get("cases", []) or []:
        project = case.get("project") or {}
        if project.get("project_id"):
            project_ids.add(str(project["project_id"]))
        if case.get("submitter_id"):
            case_ids.add(str(case["submitter_id"]))
        for sample in case.get("samples", []) or []:
            submitter_id = sample.get("submitter_id")
            if submitter_id:
                sample_ids.add(str(submitter_id))
    for barcode in file_name_barcodes(hit.get("file_name", "")):
        if len(barcode) >= 16:
            sample_ids.add(barcode)
        else:
            case_ids.add(barcode)
    return project_ids, case_ids, sample_ids


def aliquot_barcodes(sample_ids: set[str]) -> set[str]:
    return {sample for sample in sample_ids if len(sample) > 16}


def manifest_from_hits(hits: list[dict], pilot: pd.DataFrame) -> pd.DataFrame:
    if not hits:
        return empty_manifest()
    pilot_samples = set(pilot["sample_barcode"].astype(str))
    pilot_patients = set(pilot["patient_barcode"].astype(str))
    patient_to_samples = (
        pilot.groupby("patient_barcode")["sample_barcode"].apply(lambda x: sorted(set(x.astype(str)))).to_dict()
    )
    rows: list[dict] = []
    for hit in hits:
        gdc_file_id = str(hit.get("file_id") or hit.get("id") or "")
        project_ids, case_ids, sample_ids = extract_hit_ids(hit)
        normalized_samples = {normalize_sample_barcode(sample) for sample in sample_ids if sample}
        normalized_patients = {normalize_patient_barcode(item) for item in (sample_ids | case_ids) if item}
        sample_matches = sorted(normalized_samples & pilot_samples)
        patient_matches = sorted(normalized_patients & pilot_patients)
        common = {
            "gdc_file_id": gdc_file_id,
            "file_name": str(hit.get("file_name", "")),
            "data_category": str(hit.get("data_category", "")),
            "data_type": str(hit.get("data_type", "")),
            "experimental_strategy": str(hit.get("experimental_strategy", "")),
            "access": str(hit.get("access", "")),
            "file_size": str(hit.get("file_size", "")),
            "project_id": safe_join(project_ids),
            "case_submitter_ids": safe_join(case_ids),
            "sample_submitter_ids": safe_join(sample_ids),
            "aliquot_barcodes": safe_join(aliquot_barcodes(sample_ids)),
            "download_status": "not_matched",
            "local_path": "NA",
            "parse_status": "not_parsed",
        }
        if sample_matches:
            for sample in sample_matches:
                rows.append(
                    {
                        **common,
                        "matched_sample_barcode": sample,
                        "matched_patient_barcode": normalize_patient_barcode(sample),
                        "match_level": "sample",
                        "notes": "sample_level_match_from_gdc_sample_submitter_id_or_file_name",
                    }
                )
        elif patient_matches:
            for patient in patient_matches:
                matched_samples = patient_to_samples.get(patient, [])
                rows.append(
                    {
                        **common,
                        "matched_sample_barcode": matched_samples[0] if len(matched_samples) == 1 else safe_join(set(matched_samples)),
                        "matched_patient_barcode": patient,
                        "match_level": "patient",
                        "notes": "patient_level_match_lower_confidence_validate_before_clonal_inference",
                    }
                )
        else:
            rows.append(
                {
                    **common,
                    "matched_sample_barcode": "NA",
                    "matched_patient_barcode": "NA",
                    "match_level": "unmatched",
                    "notes": "no_overlap_with_level3_pilot_cohort",
                }
            )
    manifest = pd.DataFrame(rows, columns=MANIFEST_COLUMNS)
    return manifest.drop_duplicates().reset_index(drop=True)


def sanitize_file_name(file_name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9._-]+", "_", str(file_name).strip())
    return clean or "gdc_segment_file"


def expected_local_path(download_dir: Path, gdc_file_id: str, file_name: str) -> Path:
    return download_dir / f"{gdc_file_id}_{sanitize_file_name(file_name)}"


def find_existing_local_path(download_dir: Path, gdc_file_id: str, file_name: str) -> Path:
    expected = expected_local_path(download_dir, gdc_file_id, file_name)
    if expected.exists() and expected.stat().st_size > 0:
        return expected
    matches = sorted(download_dir.glob(f"{gdc_file_id}_*")) if download_dir.exists() else []
    for match in matches:
        if match.is_file() and match.stat().st_size > 0:
            return match
    return expected


def empty_best_files() -> pd.DataFrame:
    return pd.DataFrame(columns=BEST_FILE_COLUMNS)


def empty_download_log() -> pd.DataFrame:
    return pd.DataFrame(columns=DOWNLOAD_LOG_COLUMNS)


def data_type_priority(data_type: str) -> int:
    priorities = {
        "Masked Copy Number Segment": 0,
        "Copy Number Segment": 1,
        "Filtered Copy Number Segment": 2,
    }
    return priorities.get(str(data_type), 9)


def lightweight_segment_check(
    manifest_row: pd.Series,
    download_dir: Path,
    candidates: pd.DataFrame,
) -> tuple[bool, int | None, str]:
    local_path = find_existing_local_path(
        download_dir,
        str(manifest_row.get("gdc_file_id", "")),
        str(manifest_row.get("file_name", "")),
    )
    if not local_path.exists() or local_path.stat().st_size == 0:
        name = str(manifest_row.get("file_name", "")).lower()
        looks_parseable = any(token in name for token in [".seg", "segment", "copy_number", "copynumber"])
        return looks_parseable, None, "format_inferred_from_metadata_no_local_parse"
    metadata = {
        "gdc_file_id": str(manifest_row.get("gdc_file_id", "")),
        "file_name": str(manifest_row.get("file_name", "")),
        "project_id": str(manifest_row.get("project_id", "")),
        "match_level": str(manifest_row.get("match_level", "")),
        "matched_samples": [str(manifest_row.get("matched_sample_barcode", ""))],
        "matched_patients": [str(manifest_row.get("matched_patient_barcode", ""))],
    }
    try:
        raw = read_segment_file(local_path)
        parsed, warnings = standardize_segment_table(raw, local_path, metadata, patient_to_project_lookup(candidates))
        if parsed.empty:
            return False, 0, "local_parse_check_no_rows:" + safe_join(set(warnings))
        return True, len(parsed), "local_parse_check_passed"
    except Exception as exc:  # noqa: BLE001
        return False, 0, f"local_parse_check_failed:{exc}"


def select_best_files_per_pilot_sample(
    manifest: pd.DataFrame,
    pilot: pd.DataFrame,
    candidates: pd.DataFrame,
    download_dir: Path,
    config: dict = CONFIG,
) -> pd.DataFrame:
    if manifest.empty:
        return empty_best_files()
    matched = manifest[manifest["match_level"].isin(["sample", "patient"])].copy()
    if matched.empty:
        return empty_best_files()
    if not bool(config.get("allow_patient_level_matches", False)):
        matched = matched[matched["match_level"] == "sample"].copy()
    elif bool(config.get("prefer_sample_level_matches", True)):
        samples_with_sample_match = set(matched.loc[matched["match_level"] == "sample", "matched_sample_barcode"].astype(str))
        patient_mask = (matched["match_level"] == "patient") & matched["matched_sample_barcode"].astype(str).isin(
            samples_with_sample_match
        )
        matched = matched[~patient_mask].copy()
    if matched.empty:
        return empty_best_files()

    pilot_lookup = pilot.set_index("sample_barcode")[["patient_barcode", "project_id", "project_code"]].to_dict("index")
    rows: list[dict[str, Any]] = []
    for sample, group in matched.groupby("matched_sample_barcode", dropna=False):
        sample = str(sample)
        if not sample or sample == "NA" or ";" in sample or sample not in pilot_lookup:
            continue
        scored = group.copy()
        checks = scored.apply(lambda row: lightweight_segment_check(row, download_dir, candidates), axis=1)
        scored["_parseable"] = [item[0] for item in checks]
        scored["_n_parse_rows"] = [item[1] for item in checks]
        scored["_check_note"] = [item[2] for item in checks]
        scored["_match_priority"] = scored["match_level"].map({"sample": 0, "patient": 1}).fillna(9)
        scored["_access_priority"] = scored["access"].astype(str).str.lower().map({"open": 0}).fillna(9)
        scored["_data_type_priority"] = scored["data_type"].map(data_type_priority)
        scored["_parse_priority"] = (~scored["_parseable"].astype(bool)).astype(int)
        scored["_n_parse_rank"] = pd.to_numeric(scored["_n_parse_rows"], errors="coerce").fillna(-1).astype(int) * -1
        scored["_harmonized_priority"] = scored["file_name"].astype(str).str.lower().str.contains("grch38|hg38").map(
            {True: 0, False: 1}
        )
        scored["_file_size_numeric"] = pd.to_numeric(scored["file_size"], errors="coerce").fillna(10**18)
        sort_cols = [
            "_match_priority",
            "_access_priority",
            "_data_type_priority",
            "_parse_priority",
            "_n_parse_rank",
            "_harmonized_priority",
            "_file_size_numeric",
            "gdc_file_id",
        ]
        ranked = scored.sort_values(sort_cols).reset_index(drop=True)
        selected = ranked.iloc[0]
        rank_tie = False
        if len(ranked) > 1:
            tie_cols = sort_cols[:-2]
            rank_tie = bool((ranked.loc[1, tie_cols].astype(str).to_list() == selected[tie_cols].astype(str).to_list()))
        lookup = pilot_lookup[sample]
        reason = [
            f"{selected['match_level']}_level_match",
            str(selected["data_type"]),
            str(selected["_check_note"]),
        ]
        if pd.notna(selected["_n_parse_rows"]):
            reason.append(f"lightweight_parse_rows={int(selected['_n_parse_rows'])}")
        rows.append(
            {
                "sample_barcode": sample,
                "patient_barcode": lookup["patient_barcode"],
                "project_id": lookup["project_id"],
                "project_code": lookup["project_code"],
                "selected_gdc_file_id": str(selected["gdc_file_id"]),
                "selected_file_name": str(selected["file_name"]),
                "selected_data_type": str(selected["data_type"]),
                "match_level": str(selected["match_level"]),
                "file_size": str(selected["file_size"]),
                "n_candidate_files": int(group["gdc_file_id"].nunique()),
                "selection_reason": ";".join(reason),
                "selection_rank": 1,
                "notes": "tie_resolved_by_file_size_or_file_id" if rank_tie else "best_ranked_file_selected",
            }
        )
    best = pd.DataFrame(rows, columns=BEST_FILE_COLUMNS)
    return best.sort_values(["project_code", "sample_barcode"]).reset_index(drop=True) if not best.empty else best


def apply_download_mode(best_files: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    mode = str(config.get("download_mode", "pilot_full"))
    if best_files.empty or mode in {"none", "pilot_full"}:
        return best_files.copy()
    if mode in {"pilot_capped", "candidates_capped"}:
        max_downloads = int(config.get("max_downloads", config.get("max_download_files", MAX_DOWNLOADS)))
        return best_files.head(max_downloads).copy()
    raise ValueError(f"Unsupported DOWNLOAD_MODE: {mode}")


def manifest_for_selected_files(manifest: pd.DataFrame, selected_files: pd.DataFrame) -> pd.DataFrame:
    if manifest.empty or selected_files.empty:
        return empty_manifest()
    selected_pairs = set(
        zip(
            selected_files["selected_gdc_file_id"].astype(str),
            selected_files["sample_barcode"].astype(str),
            strict=False,
        )
    )
    selected_ids = set(selected_files["selected_gdc_file_id"].astype(str))
    selected = manifest[manifest["gdc_file_id"].astype(str).isin(selected_ids)].copy()
    if "matched_sample_barcode" in selected:
        selected = selected[
            selected.apply(
                lambda row: (
                    str(row["gdc_file_id"]),
                    str(row["matched_sample_barcode"]),
                )
                in selected_pairs,
                axis=1,
            )
        ].copy()
    return selected.reset_index(drop=True)


def download_selected_files(
    manifest: pd.DataFrame,
    selected_files: pd.DataFrame,
    download_dir: Path,
    config: dict = CONFIG,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    ensure_dir(download_dir)
    out = manifest.copy()
    warnings: list[str] = []
    log_rows: list[dict[str, Any]] = []
    if selected_files.empty:
        return out, empty_download_log(), warnings
    mode = str(config.get("download_mode", "pilot_full"))
    overwrite = bool(config.get("overwrite_existing", False))
    selected_ids = set(selected_files["selected_gdc_file_id"].astype(str))
    out.loc[out["match_level"].isin(["sample", "patient"]), "download_status"] = "download_not_selected"
    out.loc[out["gdc_file_id"].astype(str).isin(selected_ids), "download_status"] = "selected_not_downloaded"

    for row in selected_files.itertuples(index=False):
        file_id = str(row.selected_gdc_file_id)
        file_name = str(row.selected_file_name)
        local_path = find_existing_local_path(download_dir, file_id, file_name)
        attempted = False
        status = "download_not_attempted"
        error_message = "NA"
        if local_path.exists() and local_path.stat().st_size > 0 and not overwrite:
            status = "already_present"
        elif mode == "none":
            status = "file_already_absent" if not local_path.exists() else "already_present"
        else:
            attempted = True
            local_path = expected_local_path(download_dir, file_id, file_name)
            try:
                url = f"{GDC_DATA_ENDPOINT}/{file_id}"
                with requests.get(url, stream=True, timeout=int(config["download_timeout_seconds"])) as response:
                    response.raise_for_status()
                    with local_path.open("wb") as handle:
                        for chunk in response.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                handle.write(chunk)
                status = "downloaded" if local_path.stat().st_size > 0 else "download_failed"
            except requests.RequestException as exc:
                status = "download_failed"
                error_message = str(exc)
                warnings.append(f"download_failed_{file_id}:{exc}")
                logging.warning("Failed to download %s: %s", file_id, exc)
        out.loc[out["gdc_file_id"].astype(str) == file_id, "download_status"] = status
        out.loc[out["gdc_file_id"].astype(str) == file_id, "local_path"] = str(local_path) if local_path.exists() else "NA"
        log_rows.append(
            {
                "sample_barcode": str(row.sample_barcode),
                "gdc_file_id": file_id,
                "file_name": file_name,
                "download_attempted": attempted,
                "download_status": status,
                "local_path": str(local_path) if local_path.exists() else "NA",
                "error_message": error_message,
                "file_size": str(row.file_size),
            }
        )
    log = pd.DataFrame(log_rows, columns=DOWNLOAD_LOG_COLUMNS)
    return out, log, warnings


def download_matched_files(manifest: pd.DataFrame, download_dir: Path, config: dict = CONFIG) -> tuple[pd.DataFrame, list[str]]:
    ensure_dir(download_dir)
    warnings: list[str] = []
    out = manifest.copy()
    if out.empty:
        return out, warnings
    matched = out[out["match_level"].isin(["sample", "patient"])].copy()
    matched["_match_priority"] = matched["match_level"].map({"sample": 0, "patient": 1}).fillna(9)
    matched["_file_size_numeric"] = pd.to_numeric(matched["file_size"], errors="coerce")
    unique_files = (
        matched.sort_values(["_match_priority", "_file_size_numeric", "gdc_file_id"])
        [["gdc_file_id", "file_name"]]
        .drop_duplicates()
    )
    if len(unique_files) > int(config["max_download_files"]):
        allowed = set(unique_files.head(int(config["max_download_files"]))["gdc_file_id"])
        skipped = set(unique_files["gdc_file_id"]) - allowed
        out.loc[out["gdc_file_id"].isin(skipped), "download_status"] = "download_skipped_limit"
        warnings.append(f"download_limited_to_{config['max_download_files']}_files_from_{len(unique_files)}_matched_files")
        unique_files = unique_files[unique_files["gdc_file_id"].isin(allowed)]

    for row in unique_files.itertuples(index=False):
        file_id = str(row.gdc_file_id)
        file_name = sanitize_file_name(str(row.file_name))
        local_path = download_dir / f"{file_id}_{file_name}"
        if local_path.exists() and local_path.stat().st_size > 0:
            status = "already_present"
        else:
            try:
                url = f"{GDC_DATA_ENDPOINT}/{file_id}"
                with requests.get(url, stream=True, timeout=int(config["download_timeout_seconds"])) as response:
                    response.raise_for_status()
                    with local_path.open("wb") as handle:
                        for chunk in response.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                handle.write(chunk)
                status = "downloaded" if local_path.stat().st_size > 0 else "download_failed"
            except requests.RequestException as exc:
                status = "download_failed"
                warnings.append(f"download_failed_{file_id}:{exc}")
                logging.warning("Failed to download %s: %s", file_id, exc)
        out.loc[out["gdc_file_id"] == file_id, "download_status"] = status
        out.loc[out["gdc_file_id"] == file_id, "local_path"] = str(local_path) if local_path.exists() else "NA"
    return out, warnings


def metadata_by_file(manifest: pd.DataFrame) -> dict[str, dict[str, Any]]:
    metadata: dict[str, dict[str, Any]] = {}
    matched = manifest[manifest["match_level"].isin(["sample", "patient"])].copy()
    for file_id, group in matched.groupby("gdc_file_id", dropna=False):
        sample_rows = group[group["match_level"] == "sample"]
        if not sample_rows.empty:
            match_level = "sample"
            matched_samples = sorted(set(sample_rows["matched_sample_barcode"].astype(str)))
            matched_patients = sorted(set(sample_rows["matched_patient_barcode"].astype(str)))
        else:
            match_level = "patient"
            patient_rows = group[group["match_level"] == "patient"]
            matched_samples = sorted(
                sample
                for sample in set(patient_rows["matched_sample_barcode"].astype(str))
                if sample and sample != "NA" and ";" not in sample
            )
            matched_patients = sorted(set(patient_rows["matched_patient_barcode"].astype(str)))
        metadata[str(file_id)] = {
            "gdc_file_id": str(file_id),
            "file_name": str(group["file_name"].iloc[0]),
            "project_id": str(group["project_id"].iloc[0]),
            "match_level": match_level,
            "matched_samples": matched_samples,
            "matched_patients": matched_patients,
        }
    return metadata


def read_segment_file(path: Path) -> pd.DataFrame:
    try:
        raw = pd.read_csv(path, sep="\t", comment="#", dtype=str, compression="infer")
        if len(raw.columns) > 1:
            return raw
    except Exception:
        pass
    return pd.read_csv(path, sep=None, engine="python", comment="#", dtype=str, compression="infer")


def patient_to_project_lookup(candidates: pd.DataFrame) -> dict[str, tuple[str, str]]:
    lookup = candidates[["patient_barcode", "project_id", "project_code"]].drop_duplicates()
    return {
        str(row.patient_barcode): (str(row.project_id), str(row.project_code))
        for row in lookup.itertuples(index=False)
    }


def standardize_segment_table(
    raw: pd.DataFrame,
    source_file: Path,
    metadata: dict[str, Any],
    patient_project_lookup: dict[str, tuple[str, str]],
) -> tuple[pd.DataFrame, list[str]]:
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
    missing = [
        name
        for name, column in {
            "chromosome": chrom_col,
            "start": start_col,
            "end": end_col,
            "segment_mean": segment_mean_col,
        }.items()
        if column is None
    ]
    if missing:
        return empty_segments(), [f"{source_file.name}:missing_required_columns_{','.join(missing)}"]

    matched_samples = set(metadata.get("matched_samples", []))
    matched_patients = set(metadata.get("matched_patients", []))
    match_level = str(metadata.get("match_level", "unmatched"))
    if sample_col and raw[sample_col].astype(str).str.startswith("TCGA-").any():
        raw_sample = raw[sample_col].map(normalize_sample_barcode)
        raw_patient = raw[sample_col].map(normalize_patient_barcode)
    elif len(matched_samples) == 1:
        sample = next(iter(matched_samples))
        raw_sample = pd.Series([sample] * len(raw), index=raw.index)
        raw_patient = pd.Series([normalize_patient_barcode(sample)] * len(raw), index=raw.index)
        warnings.append(f"{source_file.name}:sample_barcode_assigned_from_gdc_manifest")
    else:
        return empty_segments(), [f"{source_file.name}:missing_sample_column_and_no_unique_manifest_sample_match"]

    data = raw.copy()
    data["_raw_sample_barcode"] = raw_sample
    data["_raw_patient_barcode"] = raw_patient
    if match_level == "sample" and matched_samples:
        data = data[data["_raw_sample_barcode"].isin(matched_samples)].copy()
    elif match_level == "patient" and matched_patients:
        data = data[data["_raw_patient_barcode"].isin(matched_patients)].copy()
        patient_to_matched_sample = {
            normalize_patient_barcode(sample): sample for sample in matched_samples if sample and sample != "NA"
        }
        if patient_to_matched_sample:
            data["_raw_sample_barcode"] = data["_raw_patient_barcode"].map(patient_to_matched_sample).fillna(
                data["_raw_sample_barcode"]
            )
    if data.empty:
        return empty_segments(), [f"{source_file.name}:no_rows_overlap_manifest_match"]

    out = pd.DataFrame(
        {
            "sample_barcode": data["_raw_sample_barcode"],
            "patient_barcode": data["_raw_patient_barcode"].map(normalize_patient_barcode),
            "chromosome": data[chrom_col].astype(str).str.replace("^chr", "", regex=True),
            "start": numeric(data[start_col]),
            "end": numeric(data[end_col]),
            "num_probes": numeric(data[num_probes_col]) if num_probes_col else pd.NA,
            "segment_mean": numeric(data[segment_mean_col]),
            "source_file": source_file.name,
            "gdc_file_id": str(metadata.get("gdc_file_id", "")),
            "match_level": match_level,
            "notes": (
                "gdc_segment_cn_sample_level_match"
                if match_level == "sample"
                else "gdc_segment_cn_patient_level_lower_confidence_match_validate_before_clonal_inference"
            ),
        }
    )
    project_pairs = out["patient_barcode"].map(
        lambda patient: patient_project_lookup.get(patient, (str(metadata.get("project_id", "TCGA-UNKNOWN")), "UNKNOWN"))
    )
    out["project_id"] = [item[0] for item in project_pairs]
    out["project_code"] = [
        item[1] if item[1] != "UNKNOWN" else project_code_from_project_id(item[0]) for item in project_pairs
    ]
    out = out[
        out["sample_barcode"].astype(str).str.startswith("TCGA-")
        & out["start"].notna()
        & out["end"].notna()
        & out["segment_mean"].notna()
    ].copy()
    out = out[out["end"] >= out["start"]]
    return out[SEGMENT_COLUMNS], warnings


def parse_downloaded_segments(manifest: pd.DataFrame, candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    out_manifest = manifest.copy()
    warnings: list[str] = []
    frames: list[pd.DataFrame] = []
    metadata = metadata_by_file(out_manifest)
    project_lookup = patient_to_project_lookup(candidates)
    downloaded = out_manifest[out_manifest["download_status"].isin(["downloaded", "already_present"])]
    for file_id, group in downloaded.groupby("gdc_file_id", dropna=False):
        local_paths = [Path(path) for path in group["local_path"].astype(str).unique() if path and path != "NA"]
        if not local_paths:
            out_manifest.loc[out_manifest["gdc_file_id"] == file_id, "parse_status"] = "not_downloaded"
            continue
        local_path = local_paths[0]
        try:
            raw = read_segment_file(local_path)
            parsed, file_warnings = standardize_segment_table(raw, local_path, metadata[str(file_id)], project_lookup)
            warnings.extend(file_warnings)
            if parsed.empty:
                out_manifest.loc[out_manifest["gdc_file_id"] == file_id, "parse_status"] = "parse_no_rows"
            else:
                frames.append(parsed)
                out_manifest.loc[out_manifest["gdc_file_id"] == file_id, "parse_status"] = "parsed"
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"{local_path.name}:parse_error:{exc}")
            out_manifest.loc[out_manifest["gdc_file_id"] == file_id, "parse_status"] = "parse_error"
    if not frames:
        return empty_segments(), out_manifest, warnings
    segments = pd.concat(frames, ignore_index=True)
    segments = segments.drop_duplicates(
        subset=["sample_barcode", "chromosome", "start", "end", "segment_mean", "gdc_file_id"]
    )
    return segments[SEGMENT_COLUMNS], out_manifest, warnings


def best_match_level(values: pd.Series) -> str:
    levels = set(values.dropna().astype(str))
    if "sample" in levels:
        return "sample"
    if "patient" in levels:
        return "patient"
    return "unmatched"


def segment_summary(segments: pd.DataFrame, candidates: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    columns = [
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
        "segment_cn_match_level",
        "segment_cn_quality_flag",
        "source_files",
        "gdc_file_ids",
        "notes",
    ]
    if segments.empty:
        return pd.DataFrame(columns=columns)
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
        segment_cn_match_level=("match_level", best_match_level),
        source_files=("source_file", lambda x: safe_join(set(x.astype(str)))),
        gdc_file_ids=("gdc_file_id", lambda x: safe_join(set(x.astype(str)))),
    ).reset_index()
    summary["has_segment_level_cn"] = summary["n_segments"] >= int(config["min_segments_for_usable_cn"])
    summary["segment_cn_quality_flag"] = summary.apply(
        lambda row: (
            "usable_segment_cn_sample_level_match"
            if row["has_segment_level_cn"] and row["segment_cn_match_level"] == "sample"
            else (
                "usable_segment_cn_patient_level_lower_confidence"
                if row["has_segment_level_cn"]
                else "insufficient_segments"
            )
        ),
        axis=1,
    )
    summary["notes"] = summary["segment_cn_quality_flag"]
    known = candidates[["sample_barcode", "patient_barcode", "project_id", "project_code"]].drop_duplicates()
    summary = summary.drop(columns=["patient_barcode", "project_id", "project_code"]).merge(
        known, on="sample_barcode", how="left"
    ).merge(
        summary[["sample_barcode", "patient_barcode", "project_id", "project_code"]],
        on="sample_barcode",
        how="left",
        suffixes=("", "_from_segment"),
    )
    for column in ["patient_barcode", "project_id", "project_code"]:
        summary[column] = summary[column].fillna(summary[f"{column}_from_segment"])
        summary = summary.drop(columns=[f"{column}_from_segment"])
    return summary[columns]


def minimum_requirements_met(row: pd.Series) -> bool:
    return (
        bool(row.get("ref_alt_counts_available", False))
        and bool(row.get("purity_available", False))
        and bool(row.get("ploidy_available", False))
        and bool(row.get("copy_number_available", False))
        and str(row.get("exclusion_reason", "none")) == "none"
    )


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
        match_level = str(row.get("segment_cn_match_level", "unknown"))
        notes.append(f"segment_level_cn_available_from_gdc_{match_level}_match")
        if match_level == "patient":
            notes.append("patient_level_segment_cn_lower_confidence_validate_before_clonal_inference")
    elif row.get("candidate_class") == "clonal_clustering_candidate_limited_cn":
        notes.append("segment_level_cn_missing_suitable_for_clonal_clustering_prototype_not_full_cn_aware_phylogeny")
    return ";".join(dict.fromkeys(notes))


def update_candidates_with_segments(candidates: pd.DataFrame, summary: pd.DataFrame, config: dict = CONFIG) -> pd.DataFrame:
    merge_cols = [
        "sample_barcode",
        "n_segments",
        "n_autosomal_segments",
        "median_segment_mean",
        "sd_segment_mean",
        "max_abs_segment_mean",
        "fraction_segments_gain_like",
        "fraction_segments_loss_like",
        "has_segment_level_cn",
        "segment_cn_match_level",
        "segment_cn_quality_flag",
        "source_files",
        "gdc_file_ids",
    ]
    available_cols = [column for column in merge_cols if column in summary.columns]
    merged = candidates.drop(
        columns=[
            "has_segment_level_cn",
            "n_segments",
            "n_autosomal_segments",
            "median_segment_mean",
            "sd_segment_mean",
            "max_abs_segment_mean",
            "fraction_segments_gain_like",
            "fraction_segments_loss_like",
            "segment_cn_match_level",
            "segment_cn_quality_flag",
            "source_files",
            "gdc_file_ids",
        ],
        errors="ignore",
    ).merge(summary[available_cols], on="sample_barcode", how="left")
    merged["has_segment_level_cn"] = merged["has_segment_level_cn"].where(merged["has_segment_level_cn"].notna(), False).astype(bool)
    for column in ["segment_cn_match_level", "segment_cn_quality_flag", "source_files", "gdc_file_ids"]:
        if column not in merged.columns:
            merged[column] = "NA"
        merged[column] = merged[column].fillna("NA")
    sample_level_segment_mask = merged["segment_cn_match_level"].astype(str) == "sample"
    patient_level_allowed_mask = (
        bool(config.get("allow_patient_level_matches", False))
        & (merged["segment_cn_match_level"].astype(str) == "patient")
    )
    upgrade_mask = (
        merged["has_segment_level_cn"]
        & (sample_level_segment_mask | patient_level_allowed_mask)
        & merged.apply(minimum_requirements_met, axis=1)
    )
    limited_mask = (~merged["has_segment_level_cn"]) & merged["candidate_class"].isin(
        ["copy_number_aware_clonal_phylogeny_candidate", "clonal_clustering_candidate_limited_cn"]
    )
    merged.loc[upgrade_mask, "candidate_class"] = "copy_number_aware_clonal_phylogeny_candidate"
    merged.loc[limited_mask, "candidate_class"] = "clonal_clustering_candidate_limited_cn"
    merged["notes"] = merged.apply(update_notes, axis=1)
    return merged


def update_pilot_with_segments(pilot: pd.DataFrame, updated_candidates: pd.DataFrame) -> pd.DataFrame:
    selected = set(pilot["sample_barcode"])
    updated = updated_candidates[updated_candidates["sample_barcode"].isin(selected)].copy()
    if updated.empty:
        return pilot.copy()
    updated["inclusion_status"] = "selected_for_pilot"
    if not updated["has_segment_level_cn"].any():
        updated["notes"] = updated["notes"].astype(str) + ";pilot_preserved_no_gdc_segment_level_cn_available"
    return updated


def availability_by_project(candidates: pd.DataFrame, pilot: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    manifest_counts = pd.DataFrame(
        columns=["project_code", "n_gdc_files_found_for_project", "n_gdc_files_matched_to_pilot"]
    )
    if not manifest.empty:
        rows = []
        for project_id in sorted(set(";".join(manifest["project_id"].astype(str)).split(";"))):
            if not project_id or project_id == "NA":
                continue
            project_manifest = manifest[manifest["project_id"].astype(str).str.contains(project_id, regex=False, na=False)]
            rows.append(
                {
                    "project_code": project_code_from_project_id(project_id),
                    "n_gdc_files_found_for_project": project_manifest["gdc_file_id"].nunique(),
                    "n_gdc_files_matched_to_pilot": project_manifest[
                        project_manifest["match_level"].isin(["sample", "patient"])
                    ]["gdc_file_id"].nunique(),
                }
            )
        manifest_counts = pd.DataFrame(rows)
    rows = []
    for (project_id, project_code), group in candidates.groupby(["project_id", "project_code"], dropna=False):
        pilot_group = pilot[pilot["project_code"] == project_code]
        n_candidates = len(group)
        n_with = int(group["has_segment_level_cn"].sum())
        n_pilot = len(pilot_group)
        n_pilot_with = int(pilot_group["has_segment_level_cn"].sum()) if "has_segment_level_cn" in pilot_group else 0
        sample_matches = int((pilot_group.get("segment_cn_match_level", pd.Series(dtype=str)) == "sample").sum())
        patient_matches = int((pilot_group.get("segment_cn_match_level", pd.Series(dtype=str)) == "patient").sum())
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
                "n_pilot_sample_level_matches": sample_matches,
                "n_pilot_patient_level_matches": patient_matches,
                "notes": "segment_cn_available" if n_with else "no_segment_level_cn_overlap",
            }
        )
    availability = pd.DataFrame(rows)
    if not manifest_counts.empty:
        availability = availability.merge(manifest_counts, on="project_code", how="left")
    for column in ["n_gdc_files_found_for_project", "n_gdc_files_matched_to_pilot"]:
        if column not in availability.columns:
            availability[column] = 0
        availability[column] = availability[column].fillna(0).astype(int)
    return availability.sort_values(["n_candidates_with_segments", "n_level3_candidates"], ascending=[False, False])


def manual_download_instructions(filters: dict) -> str:
    return (
        "No usable GDC pilot segment-level CN rows were parsed. Manual option: use the GDC API/Data Portal with "
        f"filters={json.dumps(filters, sort_keys=True)}; download open-access TCGA Copy Number Variation files with "
        "data_type in Masked Copy Number Segment, Copy Number Segment, or Filtered Copy Number Segment for pilot "
        "samples, then place them under data/raw/copy_number_segments/gdc_pilot/ or data/raw/copy_number_segments/."
    )


def pilot_missing_reasons(
    pilot: pd.DataFrame,
    manifest: pd.DataFrame,
    best_files: pd.DataFrame,
    download_log: pd.DataFrame,
    updated_pilot: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    sample_matches = (
        set(manifest.loc[manifest["match_level"] == "sample", "matched_sample_barcode"].astype(str))
        if not manifest.empty
        else set()
    )
    patient_matches = (
        set(manifest.loc[manifest["match_level"] == "patient", "matched_sample_barcode"].astype(str))
        if not manifest.empty
        else set()
    )
    selected = set(best_files["sample_barcode"].astype(str)) if not best_files.empty else set()
    downloaded_status = (
        dict(zip(download_log["sample_barcode"].astype(str), download_log["download_status"].astype(str), strict=False))
        if not download_log.empty
        else {}
    )
    parsed_samples = (
        set(updated_pilot.loc[updated_pilot["has_segment_level_cn"].astype(bool), "sample_barcode"].astype(str))
        if "has_segment_level_cn" in updated_pilot
        else set()
    )
    for sample in pilot["sample_barcode"].astype(str):
        if sample in parsed_samples:
            continue
        status = downloaded_status.get(sample, "")
        if sample not in sample_matches and sample in patient_matches:
            reason = "patient-level only"
        elif sample not in sample_matches:
            reason = "no sample-level match"
        elif sample not in selected:
            reason = "other"
        elif status == "download_failed":
            reason = "download failed"
        elif status == "file_already_absent":
            reason = "file already absent"
        elif status in {"downloaded", "already_present"}:
            reason = "parse failed"
        else:
            reason = "other"
        rows.append({"sample_barcode": sample, "missing_reason": reason})
    return pd.DataFrame(rows, columns=["sample_barcode", "missing_reason"])


def qc_summary(
    manifest: pd.DataFrame,
    best_files: pd.DataFrame,
    download_log: pd.DataFrame,
    segments: pd.DataFrame,
    summary: pd.DataFrame,
    candidates: pd.DataFrame,
    updated_candidates: pd.DataFrame,
    pilot: pd.DataFrame,
    updated_pilot: pd.DataFrame,
    filters: dict,
    warnings: list[str],
    config: dict = CONFIG,
) -> pd.DataFrame:
    matched = manifest[manifest["match_level"].isin(["sample", "patient"])] if not manifest.empty else manifest
    downloaded = manifest[manifest["download_status"] == "downloaded"] if not manifest.empty else manifest
    present = manifest[manifest["download_status"] == "already_present"] if not manifest.empty else manifest
    parsed = manifest[manifest["parse_status"] == "parsed"] if not manifest.empty else manifest
    sample_matches = int((matched["match_level"] == "sample").sum()) if not matched.empty else 0
    patient_matches = int((matched["match_level"] == "patient").sum()) if not matched.empty else 0
    sample_match_samples = (
        int(manifest.loc[manifest["match_level"] == "sample", "matched_sample_barcode"].nunique())
        if not manifest.empty
        else 0
    )
    any_match_samples = (
        int(matched["matched_sample_barcode"].nunique()) if not matched.empty and "matched_sample_barcode" in matched else 0
    )
    selected_samples = int(best_files["sample_barcode"].nunique()) if not best_files.empty else 0
    downloaded_samples = (
        int(download_log[download_log["download_status"].isin(["downloaded", "already_present"])]["sample_barcode"].nunique())
        if not download_log.empty
        else 0
    )
    pilot_projects = set(pilot["project_code"].astype(str))
    matched_projects: set[str] = set()
    if not matched.empty:
        for value in matched["project_id"].dropna().astype(str):
            matched_projects.update(
                project_code_from_project_id(project_id)
                for project_id in value.split(";")
                if project_id and project_id != "NA"
            )
    missing_projects = sorted(pilot_projects - matched_projects)
    upgrade_table = candidates[["sample_barcode", "candidate_class"]].merge(
        updated_candidates[["sample_barcode", "candidate_class"]],
        on="sample_barcode",
        how="right",
        suffixes=("_before", "_after"),
    )
    upgraded = int(
        (
            (upgrade_table["candidate_class_after"] == "copy_number_aware_clonal_phylogeny_candidate")
            & (upgrade_table["candidate_class_before"] != "copy_number_aware_clonal_phylogeny_candidate")
        ).sum()
    )
    upgraded_pilot = int(
        updated_pilot["candidate_class"].astype(str).eq("copy_number_aware_clonal_phylogeny_candidate").sum()
    ) if "candidate_class" in updated_pilot else 0
    missing = pilot_missing_reasons(pilot, manifest, best_files, download_log, updated_pilot)
    missing_counts = missing["missing_reason"].value_counts().to_dict() if not missing.empty else {}
    rows = [
        ("global", "gdc_query_endpoint", GDC_FILES_ENDPOINT),
        ("global", "gdc_query_filters", json.dumps(filters, sort_keys=True)),
        ("global", "download_mode", str(config.get("download_mode", "pilot_full"))),
        ("global", "gdc_segment_files_found", manifest["gdc_file_id"].nunique() if not manifest.empty else 0),
        ("global", "gdc_segment_files_matched_to_pilot", matched["gdc_file_id"].nunique() if not matched.empty else 0),
        ("global", "gdc_manifest_match_rows", len(matched) if not matched.empty else 0),
        ("global", "pilot_samples_total", len(pilot)),
        ("global", "pilot_samples_with_any_gdc_match", any_match_samples),
        ("global", "pilot_samples_with_sample_level_match", sample_match_samples),
        ("global", "pilot_samples_with_selected_best_file", selected_samples),
        ("global", "selected_best_files", best_files["selected_gdc_file_id"].nunique() if not best_files.empty else 0),
        ("global", "gdc_unique_files_downloaded", downloaded["gdc_file_id"].nunique() if not downloaded.empty else 0),
        ("global", "gdc_unique_files_already_present", present["gdc_file_id"].nunique() if not present.empty else 0),
        ("global", "gdc_unique_files_parsed_successfully", parsed["gdc_file_id"].nunique() if not parsed.empty else 0),
        ("global", "pilot_samples_downloaded", downloaded_samples),
        ("global", "total_segment_rows_processed", len(segments)),
        ("global", "samples_with_segment_level_cn", len(summary)),
        (
            "global",
            "pilot_samples_with_segment_level_cn",
            int(updated_pilot["has_segment_level_cn"].sum()) if "has_segment_level_cn" in updated_pilot else 0,
        ),
        (
            "global",
            "pilot_samples_parsed",
            int(updated_pilot["has_segment_level_cn"].sum()) if "has_segment_level_cn" in updated_pilot else 0,
        ),
        ("global", "overlap_with_level3_candidates", int(updated_candidates["has_segment_level_cn"].sum())),
        ("global", "candidates_upgraded_to_copy_number_aware", upgraded),
        ("global", "pilot_samples_upgraded", upgraded_pilot),
        (
            "global",
            "pilot_samples_still_lacking_segment_cn",
            len(pilot) - int(updated_pilot["has_segment_level_cn"].sum()) if "has_segment_level_cn" in updated_pilot else len(pilot),
        ),
        ("global", "sample_level_manifest_matches", sample_matches),
        ("global", "patient_level_manifest_matches_lower_confidence", patient_matches),
        ("global", "pilot_projects_with_no_matching_segment_files", ";".join(missing_projects) if missing_projects else "none"),
        ("global", "parsing_and_download_warnings", ";".join(warnings) if warnings else "none"),
        ("global", "limitations", "none" if len(summary) else manual_download_instructions(filters)),
    ]
    for reason in [
        "no sample-level match",
        "download failed",
        "parse failed",
        "file already absent",
        "patient-level only",
        "other",
    ]:
        rows.append(("missing_reason", f"pilot_missing_reason_{reason.replace(' ', '_').replace('-', '_')}", missing_counts.get(reason, 0)))
    if patient_matches:
        rows.append(
            (
                "global",
                "patient_level_match_caveat",
                "Patient-level segment matches are lower confidence and should be manually validated before full copy-number-aware clonal inference.",
            )
        )
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def plot_availability(availability: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = availability.sort_values("n_pilot_samples_with_segments", ascending=False)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.bar(data["project_code"], data["n_pilot_samples_with_segments"], color="#4daf4a")
    ax.set_ylabel("Pilot samples with segment-level CN")
    ax.set_title("GDC pilot segment-level CN availability by project")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def qc_metric_int(qc: pd.DataFrame, metric: str) -> int:
    values = qc.loc[qc["metric"] == metric, "value"]
    if values.empty:
        return 0
    return int(float(values.iloc[0]))


def plot_pilot_coverage(qc: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    metrics = [
        ("Total", "pilot_samples_total"),
        ("Any GDC match", "pilot_samples_with_any_gdc_match"),
        ("Sample match", "pilot_samples_with_sample_level_match"),
        ("Selected", "pilot_samples_with_selected_best_file"),
        ("Downloaded/present", "pilot_samples_downloaded"),
        ("Parsed", "pilot_samples_parsed"),
        ("Upgraded", "pilot_samples_upgraded"),
    ]
    labels = [item[0] for item in metrics]
    values = [qc_metric_int(qc, item[1]) for item in metrics]
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(labels, values, color=["#4c78a8", "#72b7b2", "#54a24b", "#eeca3b", "#f58518", "#b279a2", "#e45756"])
    ax.set_ylabel("Pilot samples")
    ax.set_title("Level 3 pilot segment-CN coverage")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_upgrade_by_project(updated_pilot: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = updated_pilot.copy()
    if "candidate_class" not in data:
        data["candidate_class"] = "NA"
    grouped = (
        data.assign(
            upgraded=data["candidate_class"].astype(str).eq("copy_number_aware_clonal_phylogeny_candidate")
        )
        .groupby("project_code", dropna=False)
        .agg(n_pilot=("sample_barcode", "size"), n_upgraded=("upgraded", "sum"))
        .reset_index()
        .sort_values(["n_upgraded", "n_pilot"], ascending=[False, False])
    )
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.bar(grouped["project_code"], grouped["n_upgraded"], color="#54a24b")
    ax.set_ylabel("Pilot samples upgraded")
    ax.set_title("Copy-number-aware Level 3 pilot upgrades by project")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_segment_counts(summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(9, 6))
    if summary.empty:
        ax.text(0.5, 0.5, "No GDC segment-level CN rows parsed", ha="center", va="center")
        ax.axis("off")
    else:
        ax.hist(summary["n_segments"], bins=30, color="#377eb8")
        ax.set_xlabel("Segments per pilot sample")
        ax.set_ylabel("Samples")
        ax.set_title("GDC segment-level CN segment-count distribution")
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
    ax.set_title("Level 3 candidate class after GDC segment CN join")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def run_gdc_segment_fetch(paths: GdcSegmentPaths, config: dict = CONFIG) -> dict[str, pd.DataFrame]:
    require_file(paths.candidate_samples, "Level 3 candidate sample table")
    require_file(paths.pilot_cohort, "Level 3 pilot cohort table")
    candidates = pd.read_csv(paths.candidate_samples, sep="\t")
    pilot = pd.read_csv(paths.pilot_cohort, sep="\t")
    project_ids = sorted(set(pilot["project_id"].dropna().astype(str)))
    logging.info("Querying GDC for pilot projects: %s", ", ".join(project_ids))
    hits, filters, query_warnings = query_gdc_files(project_ids, config)
    if not hits:
        logging.warning("GDC query returned no files. Exact filters used: %s", json.dumps(filters, sort_keys=True))
    if hits:
        manifest = manifest_from_hits(hits, pilot)
    elif paths.manifest.exists() and paths.manifest.stat().st_size > 0:
        manifest = pd.read_csv(paths.manifest, sep="\t")
        query_warnings.append("gdc_query_returned_no_files_reused_existing_manifest")
    else:
        manifest = empty_manifest()
    best_files = select_best_files_per_pilot_sample(manifest, pilot, candidates, paths.raw_download_dir, config)
    selected_for_download = apply_download_mode(best_files, config)
    if len(selected_for_download) < len(best_files):
        query_warnings.append(
            f"download_limited_to_{len(selected_for_download)}_selected_files_from_{len(best_files)}_best_files"
        )
    manifest, download_log, download_warnings = download_selected_files(
        manifest, selected_for_download, paths.raw_download_dir, config
    )
    selected_manifest = manifest_for_selected_files(manifest, selected_for_download)
    segments, selected_manifest, parse_warnings = parse_downloaded_segments(selected_manifest, candidates)
    if not selected_manifest.empty:
        status_map = selected_manifest.set_index("gdc_file_id")["parse_status"].to_dict()
        status_ids = set(status_map)
        manifest.loc[manifest["gdc_file_id"].isin(status_ids), "parse_status"] = manifest.loc[
            manifest["gdc_file_id"].isin(status_ids), "gdc_file_id"
        ].map(status_map)
    summary = segment_summary(segments, candidates, config)
    updated_candidates = update_candidates_with_segments(candidates, summary, config)
    updated_pilot = update_pilot_with_segments(pilot, updated_candidates)
    availability = availability_by_project(updated_candidates, updated_pilot, manifest)
    warnings = query_warnings + download_warnings + parse_warnings
    qc = qc_summary(
        manifest,
        best_files,
        download_log,
        segments,
        summary,
        candidates,
        updated_candidates,
        pilot,
        updated_pilot,
        filters,
        warnings,
        config,
    )

    write_tsv(manifest if not manifest.empty else empty_manifest(), paths.manifest)
    write_tsv(best_files if not best_files.empty else empty_best_files(), paths.best_file_per_pilot_sample)
    write_tsv(download_log if not download_log.empty else empty_download_log(), paths.download_log)
    write_tsv_gz(segments[SEGMENT_COLUMNS], paths.processed_segments_gz)
    write_tsv(summary, paths.segment_summary)
    write_tsv(updated_candidates, paths.candidates_with_segments)
    write_tsv(updated_pilot, paths.pilot_with_segments)
    write_tsv(qc, paths.qc_summary)
    write_tsv(availability, paths.availability_by_project)
    plot_availability(availability, paths.availability_figure)
    plot_pilot_coverage(qc, paths.pilot_coverage_figure)
    plot_upgrade_by_project(updated_pilot, paths.upgrade_by_project_figure)
    plot_segment_counts(summary, paths.count_distribution_figure)
    plot_upgrade_status(updated_candidates, paths.upgrade_status_figure)
    return {
        "manifest": manifest,
        "best_files": best_files,
        "download_log": download_log,
        "segments": segments,
        "summary": summary,
        "candidates": updated_candidates,
        "pilot": updated_pilot,
        "availability": availability,
        "qc": qc,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch pilot-cohort TCGA segment-level copy-number files from the GDC.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument(
        "--max-download-files",
        type=int,
        default=CONFIG["max_downloads"],
        help="Maximum selected GDC segment files to download in capped modes.",
    )
    parser.add_argument(
        "--download-mode",
        choices=["none", "pilot_capped", "pilot_full", "candidates_capped"],
        default=CONFIG["download_mode"],
        help="Download behavior for selected GDC segment files.",
    )
    parser.add_argument(
        "--allow-patient-level-matches",
        action="store_true",
        help="Allow lower-confidence patient-level segment matches to be selected and upgraded.",
    )
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Redownload selected files even when a non-empty local file is already present.",
    )
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    config = dict(CONFIG)
    config["download_mode"] = args.download_mode
    config["max_downloads"] = args.max_download_files
    config["max_download_files"] = args.max_download_files
    config["allow_patient_level_matches"] = bool(args.allow_patient_level_matches) or bool(
        CONFIG["allow_patient_level_matches"]
    )
    config["overwrite_existing"] = bool(args.overwrite_existing) or bool(CONFIG["overwrite_existing"])
    paths = default_paths(root)
    outputs = run_gdc_segment_fetch(paths, config)
    logging.info(
        "GDC segment CN pilot fetch complete: %d files found, %d matched files, %d parsed files, %d pilot samples with segments",
        outputs["manifest"]["gdc_file_id"].nunique() if not outputs["manifest"].empty else 0,
        outputs["manifest"][outputs["manifest"]["match_level"].isin(["sample", "patient"])]["gdc_file_id"].nunique()
        if not outputs["manifest"].empty
        else 0,
        outputs["manifest"][outputs["manifest"]["parse_status"] == "parsed"]["gdc_file_id"].nunique()
        if not outputs["manifest"].empty
        else 0,
        int(outputs["pilot"]["has_segment_level_cn"].sum()) if "has_segment_level_cn" in outputs["pilot"] else 0,
    )


if __name__ == "__main__":
    main()

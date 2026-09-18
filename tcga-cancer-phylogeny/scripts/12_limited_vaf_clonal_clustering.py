#!/usr/bin/env python

from __future__ import annotations

import argparse
import gzip
import logging
import math
import os
import tempfile
import warnings
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
from matplotlib.backends.backend_pdf import PdfPages

from lib.common import configure_logging, ensure_dir, find_project_root, project_path, require_file

try:
    from sklearn.mixture import GaussianMixture
except Exception:  # noqa: BLE001
    GaussianMixture = None


CONFIG = {
    "min_total_depth": 20,
    "min_observed_vaf": 0.0,
    "max_observed_vaf": 1.0,
    "neutral_segment_mean_min": -0.15,
    "neutral_segment_mean_max": 0.15,
    "min_input_mutations": 50,
    "max_clusters": 6,
    "min_cluster_size": 10,
    "random_state": 1,
    "bic_simpler_model_tolerance": 10.0,
    "expected_clonal_vaf_tolerance": 0.08,
    "dominant_cluster_fraction_threshold": 0.75,
    "high_entropy_threshold": 0.75,
}

ASSIGNMENT_COLUMNS = [
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
    "cluster_id",
    "cluster_mean_vaf",
    "cluster_median_vaf",
    "cluster_size",
    "cluster_interpretation",
    "notes",
]

SAMPLE_SUMMARY_COLUMNS = [
    "sample_barcode",
    "patient_barcode",
    "project_id",
    "project_code",
    "n_input_mutations",
    "n_clusters",
    "dominant_cluster_id",
    "dominant_cluster_mean_vaf",
    "dominant_cluster_size",
    "dominant_cluster_fraction",
    "subclonal_cluster_count",
    "subclonal_mutation_fraction",
    "vaf_entropy",
    "vaf_dispersion",
    "purity",
    "ploidy",
    "interpretation_class",
    "warnings",
]

PROJECT_SUMMARY_COLUMNS = [
    "project_code",
    "n_samples",
    "median_n_input_mutations",
    "median_n_clusters",
    "median_dominant_cluster_fraction",
    "median_subclonal_mutation_fraction",
    "median_vaf_entropy",
    "interpretation_summary",
]


@dataclass
class LimitedVafPaths:
    copy_neutral_input: Path
    selected_pilot: Path
    readiness: Path
    purity_ploidy: Path
    assignments: Path
    sample_summary: Path
    project_summary: Path
    qc_summary: Path
    density_figure: Path
    cluster_counts_figure: Path
    complexity_heatmap: Path
    example_samples_figure: Path


def default_paths(root: Path) -> LimitedVafPaths:
    return LimitedVafPaths(
        copy_neutral_input=project_path(
            "data", "processed", "clonal", "level3_copy_neutral_vaf_clustering_input.tsv.gz", root=root
        ),
        selected_pilot=project_path("results", "tables", "level3_limited_vaf_clustering_pilot_samples.tsv", root=root),
        readiness=project_path("results", "tables", "level3_clonal_input_readiness_by_sample.tsv", root=root),
        purity_ploidy=project_path("data", "processed", "features", "purity_ploidy_by_sample.tsv", root=root),
        assignments=project_path(
            "results", "clonal", "limited_vaf", "level3_limited_vaf_cluster_assignments.tsv.gz", root=root
        ),
        sample_summary=project_path("results", "tables", "level3_limited_vaf_clonal_complexity_by_sample.tsv", root=root),
        project_summary=project_path("results", "tables", "level3_limited_vaf_clonal_complexity_by_project.tsv", root=root),
        qc_summary=project_path("results", "tables", "level3_limited_vaf_clustering_qc_summary.tsv", root=root),
        density_figure=project_path(
            "results", "figures", "level3_limited_vaf_cluster_vaf_density_by_sample.pdf", root=root
        ),
        cluster_counts_figure=project_path(
            "results", "figures", "level3_limited_vaf_cluster_counts_by_project.pdf", root=root
        ),
        complexity_heatmap=project_path("results", "figures", "level3_limited_vaf_complexity_heatmap.pdf", root=root),
        example_samples_figure=project_path("results", "figures", "level3_limited_vaf_example_samples.pdf", root=root),
    )


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


def read_inputs(paths: LimitedVafPaths) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path, description in [
        (paths.copy_neutral_input, "copy-neutral VAF clustering input"),
        (paths.selected_pilot, "selected limited-VAF pilot table"),
        (paths.readiness, "clonal input readiness table"),
        (paths.purity_ploidy, "purity/ploidy table"),
    ]:
        require_file(path, description)
    copy_neutral = pd.read_csv(paths.copy_neutral_input, sep="\t")
    selected = pd.read_csv(paths.selected_pilot, sep="\t")
    readiness = pd.read_csv(paths.readiness, sep="\t")
    purity_ploidy = pd.read_csv(paths.purity_ploidy, sep="\t")
    return copy_neutral, selected, readiness, purity_ploidy


def filter_sample_mutations(data: pd.DataFrame, selected_samples: set[str], config: dict = CONFIG) -> pd.DataFrame:
    out = data[data["sample_barcode"].astype(str).isin(selected_samples)].copy()
    for col in ["total_depth", "observed_vaf", "local_segment_mean", "purity", "ploidy", "t_ref_count", "t_alt_count"]:
        out[col] = numeric(out[col])
    keep = (
        (out["total_depth"] >= int(config["min_total_depth"]))
        & (out["observed_vaf"] > float(config["min_observed_vaf"]))
        & (out["observed_vaf"] <= float(config["max_observed_vaf"]))
        & out["local_segment_mean"].between(
            float(config["neutral_segment_mean_min"]),
            float(config["neutral_segment_mean_max"]),
            inclusive="both",
        )
    )
    return out[keep].copy()


def relabel_clusters_by_mean(labels: np.ndarray, values: np.ndarray) -> np.ndarray:
    means = {label: values[labels == label].mean() for label in sorted(set(labels))}
    ordered = sorted(means, key=lambda label: means[label], reverse=True)
    mapping = {old: i + 1 for i, old in enumerate(ordered)}
    return np.array([mapping[label] for label in labels])


def quantile_fallback(values: np.ndarray, config: dict = CONFIG) -> tuple[np.ndarray, str]:
    n = len(values)
    min_size = int(config["min_cluster_size"])
    max_k = max(1, min(int(config["max_clusters"]), n // min_size))
    if max_k <= 1:
        return np.ones(n, dtype=int), "fallback_single_cluster_min_size"
    # Use broad quantile bins and merge duplicate edges if the distribution is discrete.
    for k in range(max_k, 1, -1):
        try:
            labels = np.asarray(pd.qcut(values, q=k, labels=False, duplicates="drop")) + 1
            counts = pd.Series(labels).value_counts()
            if counts.min() >= min_size and counts.nunique() > 0:
                labels = relabel_clusters_by_mean(labels, values)
                return labels, f"fallback_quantile_binning_k_{len(set(labels))}"
        except ValueError:
            continue
    return np.ones(n, dtype=int), "fallback_single_cluster_quantile_failed"


def fit_gmm_labels(values: np.ndarray, config: dict = CONFIG) -> tuple[np.ndarray, str, dict[str, Any]]:
    n = len(values)
    if n < int(config["min_input_mutations"]):
        return np.array([], dtype=int), "skipped_too_few_input_mutations", {}
    if GaussianMixture is None:
        labels, method = quantile_fallback(values, config)
        return labels, method, {"bic_by_k": "NA", "warning": "sklearn_unavailable"}

    x = values.reshape(-1, 1)
    max_k = max(1, min(int(config["max_clusters"]), n // int(config["min_cluster_size"])))
    candidates: list[tuple[float, int, np.ndarray, GaussianMixture, np.ndarray]] = []
    gmm_warnings: list[str] = []
    for k in range(1, max_k + 1):
        try:
            model = GaussianMixture(n_components=k, covariance_type="full", random_state=int(config["random_state"]))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                model.fit(x)
                labels = model.predict(x)
            counts = pd.Series(labels).value_counts().sort_index()
            bic = float(model.bic(x))
            candidates.append((bic, k, labels, model, counts.to_numpy()))
        except Exception as exc:  # noqa: BLE001
            gmm_warnings.append(f"gmm_k_{k}_failed:{exc}")
    if not candidates:
        labels, method = quantile_fallback(values, config)
        return labels, method, {"bic_by_k": "NA", "warning": ";".join(gmm_warnings) if gmm_warnings else "gmm_failed"}

    bic_by_k = ",".join(f"{k}:{bic:.3f}" for bic, k, *_ in candidates)
    valid_candidates = [item for item in candidates if item[4].min() >= int(config["min_cluster_size"])]
    if valid_candidates:
        best_bic = min(item[0] for item in valid_candidates)
        tolerance = float(config.get("bic_simpler_model_tolerance", 0))
        simple_valid = [item for item in valid_candidates if item[0] <= best_bic + tolerance]
        bic, k, labels, _model, counts = sorted(simple_valid, key=lambda item: (item[1], item[0]))[0]
        if counts.min() >= int(config["min_cluster_size"]):
            labels = relabel_clusters_by_mean(labels, values)
            warning = ";".join(gmm_warnings) if gmm_warnings else "none"
            return labels, f"GaussianMixture_BIC_k_{k}", {"bic_by_k": bic_by_k, "warning": warning}

    labels, method = quantile_fallback(values, config)
    return labels, f"{method}_after_small_gmm_clusters", {
        "bic_by_k": bic_by_k,
        "warning": "all_gmm_solutions_had_small_clusters" + (f";{';'.join(gmm_warnings)}" if gmm_warnings else ""),
    }


def cluster_interpretation(mean_vaf: float, purity: float, config: dict = CONFIG) -> str:
    if pd.isna(mean_vaf) or pd.isna(purity):
        return "uncertain"
    expected = float(purity) / 2.0
    tol = float(config["expected_clonal_vaf_tolerance"])
    if abs(float(mean_vaf) - expected) <= tol:
        return "clonal_like"
    if float(mean_vaf) < expected - tol:
        return "subclonal_like"
    return "uncertain"


def entropy_from_counts(counts: pd.Series) -> float:
    total = counts.sum()
    if total <= 0 or len(counts) <= 1:
        return 0.0
    probs = counts / total
    return float(-(probs * np.log(probs)).sum() / math.log(len(probs)))


def interpretation_class(
    n_input: int,
    n_clusters: int,
    dominant_fraction: float,
    subclonal_cluster_count: int,
    entropy: float,
    warning: str,
    config: dict = CONFIG,
) -> str:
    if n_input < int(config["min_input_mutations"]) or str(warning).startswith("skipped"):
        return "low_confidence"
    if dominant_fraction >= float(config["dominant_cluster_fraction_threshold"]) and subclonal_cluster_count <= 1 and n_clusters <= 2:
        return "predominantly_clonal_like"
    if n_clusters <= 3 and entropy < float(config["high_entropy_threshold"]):
        return "oligoclonal_like"
    return "multicluster_subclonal_like"


def cluster_one_sample(sample_data: pd.DataFrame, config: dict = CONFIG) -> tuple[pd.DataFrame, dict[str, Any]]:
    sample = sample_data["sample_barcode"].iloc[0]
    values = sample_data["observed_vaf"].to_numpy(dtype=float)
    labels, method, metadata = fit_gmm_labels(values, config)
    if len(labels) == 0:
        return pd.DataFrame(columns=ASSIGNMENT_COLUMNS), {
            "sample_barcode": sample,
            "warning": method,
            "method": method,
            "n_input_mutations": len(sample_data),
        }

    assigned = sample_data.copy()
    assigned["cluster_id"] = [f"cluster_{label}" for label in labels]
    stats = assigned.groupby("cluster_id").agg(
        cluster_mean_vaf=("observed_vaf", "mean"),
        cluster_median_vaf=("observed_vaf", "median"),
        cluster_size=("mutation_id", "size"),
    )
    assigned = assigned.merge(stats, on="cluster_id", how="left")
    purity = numeric(assigned["purity"]).dropna()
    purity_value = float(purity.iloc[0]) if not purity.empty else np.nan
    assigned["cluster_interpretation"] = assigned["cluster_mean_vaf"].map(
        lambda mean: cluster_interpretation(float(mean), purity_value, config)
    )
    assigned["notes"] = (
        "limited_vaf_based_clonal_clustering;"
        "copy_neutral_near_diploid_only;"
        "not_copy_number_aware_phylogeny;"
        f"method={method};"
        "branching_order_not_inferred"
    )

    counts = assigned["cluster_id"].value_counts()
    dominant_cluster = counts.idxmax()
    dominant_size = int(counts.max())
    dominant_fraction = dominant_size / len(assigned)
    cluster_interps = assigned.groupby("cluster_id")["cluster_interpretation"].first()
    subclonal_clusters = set(cluster_interps[cluster_interps == "subclonal_like"].index)
    subclonal_mutations = assigned[assigned["cluster_id"].isin(subclonal_clusters)]
    entropy = entropy_from_counts(counts)
    warning = str(metadata.get("warning", "none"))
    summary = {
        "sample_barcode": sample,
        "patient_barcode": assigned["patient_barcode"].iloc[0],
        "project_id": assigned["project_id"].iloc[0],
        "project_code": assigned["project_code"].iloc[0],
        "n_input_mutations": len(assigned),
        "n_clusters": assigned["cluster_id"].nunique(),
        "dominant_cluster_id": dominant_cluster,
        "dominant_cluster_mean_vaf": float(stats.loc[dominant_cluster, "cluster_mean_vaf"]),
        "dominant_cluster_size": dominant_size,
        "dominant_cluster_fraction": dominant_fraction,
        "subclonal_cluster_count": len(subclonal_clusters),
        "subclonal_mutation_fraction": len(subclonal_mutations) / len(assigned),
        "vaf_entropy": entropy,
        "vaf_dispersion": float(assigned["observed_vaf"].std(ddof=0)),
        "purity": assigned["purity"].iloc[0],
        "ploidy": assigned["ploidy"].iloc[0],
        "interpretation_class": interpretation_class(
            len(assigned), assigned["cluster_id"].nunique(), dominant_fraction, len(subclonal_clusters), entropy, warning, config
        ),
        "warnings": (
            f"{method};bic_by_k={metadata.get('bic_by_k', 'NA')};"
            f"{warning};limited_vaf_clustering_not_cn_aware_no_phylogeny"
        ),
    }
    for col in ASSIGNMENT_COLUMNS:
        if col not in assigned.columns:
            assigned[col] = pd.NA
    return assigned[ASSIGNMENT_COLUMNS], summary


def run_clustering(copy_neutral: pd.DataFrame, selected_pilot: pd.DataFrame, config: dict = CONFIG) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected_samples = set(selected_pilot["sample_barcode"].astype(str))
    filtered = filter_sample_mutations(copy_neutral, selected_samples, config)
    assignment_frames: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    for sample in selected_pilot["sample_barcode"].astype(str):
        sample_data = filtered[filtered["sample_barcode"].astype(str) == sample].copy()
        if len(sample_data) < int(config["min_input_mutations"]):
            pilot_row = selected_pilot[selected_pilot["sample_barcode"].astype(str) == sample].iloc[0]
            summary_rows.append(
                {
                    "sample_barcode": sample,
                    "patient_barcode": pd.NA,
                    "project_id": f"TCGA-{pilot_row['project_code']}",
                    "project_code": pilot_row["project_code"],
                    "n_input_mutations": len(sample_data),
                    "n_clusters": 0,
                    "dominant_cluster_id": "NA",
                    "dominant_cluster_mean_vaf": pd.NA,
                    "dominant_cluster_size": 0,
                    "dominant_cluster_fraction": 0,
                    "subclonal_cluster_count": 0,
                    "subclonal_mutation_fraction": 0,
                    "vaf_entropy": 0,
                    "vaf_dispersion": pd.NA,
                    "purity": pilot_row.get("purity", pd.NA),
                    "ploidy": pilot_row.get("ploidy", pd.NA),
                    "interpretation_class": "low_confidence",
                    "warnings": f"skipped_too_few_input_mutations_lt_{config['min_input_mutations']}",
                }
            )
            continue
        assignments, summary = cluster_one_sample(sample_data, config)
        if not assignments.empty:
            assignment_frames.append(assignments)
        summary_rows.append(summary)
    assignments = pd.concat(assignment_frames, ignore_index=True) if assignment_frames else pd.DataFrame(columns=ASSIGNMENT_COLUMNS)
    sample_summary = pd.DataFrame(summary_rows)
    for col in SAMPLE_SUMMARY_COLUMNS:
        if col not in sample_summary.columns:
            sample_summary[col] = pd.NA
    return assignments[ASSIGNMENT_COLUMNS], sample_summary[SAMPLE_SUMMARY_COLUMNS]


def build_project_summary(sample_summary: pd.DataFrame) -> pd.DataFrame:
    if sample_summary.empty:
        return pd.DataFrame(columns=PROJECT_SUMMARY_COLUMNS)
    grouped = sample_summary.groupby("project_code", dropna=False)
    out = grouped.agg(
        n_samples=("sample_barcode", "size"),
        median_n_input_mutations=("n_input_mutations", "median"),
        median_n_clusters=("n_clusters", "median"),
        median_dominant_cluster_fraction=("dominant_cluster_fraction", "median"),
        median_subclonal_mutation_fraction=("subclonal_mutation_fraction", "median"),
        median_vaf_entropy=("vaf_entropy", "median"),
    ).reset_index()
    class_summary = grouped["interpretation_class"].agg(lambda x: ";".join(f"{k}:{v}" for k, v in x.value_counts().items()))
    out = out.merge(class_summary.rename("interpretation_summary"), on="project_code", how="left")
    return out[PROJECT_SUMMARY_COLUMNS]


def build_qc(
    selected_pilot: pd.DataFrame,
    assignments: pd.DataFrame,
    sample_summary: pd.DataFrame,
    config: dict = CONFIG,
) -> pd.DataFrame:
    skipped = sample_summary[sample_summary["interpretation_class"] == "low_confidence"]
    rows = [
        ("global", "total_selected_pilot_samples", len(selected_pilot)),
        ("global", "samples_successfully_clustered", int((sample_summary["n_clusters"] > 0).sum()) if not sample_summary.empty else 0),
        ("global", "samples_skipped", len(skipped)),
        ("global", "samples_skipped_reasons", ";".join(skipped["warnings"].astype(str)) if not skipped.empty else "none"),
        ("global", "total_input_mutations", int(sample_summary["n_input_mutations"].sum()) if not sample_summary.empty else 0),
        ("global", "total_clustered_mutations", len(assignments)),
        ("global", "clustering_method_used", "GaussianMixture_BIC_per_sample_with_quantile_fallback"),
        (
            "global",
            "parameter_settings",
            ";".join(f"{key}={value}" for key, value in config.items() if key not in {"priority_projects"}),
        ),
        (
            "global",
            "limitations",
            "limited_vaf_based_clustering_in_copy_neutral_regions_only;not_full_copy_number_aware_phylogeny;"
            "branching_order_not_inferred;requires_allele_specific_integer_cn_for_pyClone_vi_or_phylowgs;"
            "single_bulk_tcga_samples_limit_phylogenetic_inference",
        ),
    ]
    if not sample_summary.empty:
        for cls, count in sample_summary["interpretation_class"].value_counts().items():
            rows.append(("interpretation_class", str(cls), int(count)))
    return pd.DataFrame(rows, columns=["qc_section", "metric", "value"])


def plot_density_by_sample(assignments: pd.DataFrame, sample_summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    samples = list(sample_summary["sample_barcode"].astype(str))
    with PdfPages(path) as pdf:
        for start in range(0, len(samples), 6):
            page_samples = samples[start : start + 6]
            fig, axes = plt.subplots(3, 2, figsize=(11, 12), squeeze=False)
            for ax, sample in zip(axes.ravel(), page_samples, strict=False):
                data = assignments[assignments["sample_barcode"].astype(str) == sample]
                if data.empty:
                    ax.text(0.5, 0.5, "Skipped", ha="center", va="center")
                    ax.axis("off")
                    continue
                for cluster, group in data.groupby("cluster_id"):
                    ax.hist(group["observed_vaf"], bins=30, alpha=0.55, label=cluster)
                ax.set_title(sample, fontsize=9)
                ax.set_xlabel("Observed VAF")
                ax.set_ylabel("Mutations")
                ax.legend(fontsize=6)
            for ax in axes.ravel()[len(page_samples) :]:
                ax.axis("off")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    logging.info("Wrote %s", path)


def plot_cluster_counts_by_project(sample_summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    data = sample_summary.groupby("project_code", dropna=False)["n_clusters"].median().reset_index()
    data = data.sort_values("n_clusters", ascending=False)
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(data["project_code"], data["n_clusters"], color="#4c78a8")
    ax.set_ylabel("Median clusters per sample")
    ax.set_title("Limited VAF cluster counts by project")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_complexity_heatmap(sample_summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    metrics = [
        "n_clusters",
        "dominant_cluster_fraction",
        "subclonal_mutation_fraction",
        "vaf_entropy",
        "vaf_dispersion",
    ]
    data = sample_summary.sort_values(["project_code", "sample_barcode"]).set_index("sample_barcode")
    matrix = data[metrics].apply(pd.to_numeric, errors="coerce").fillna(0)
    scaled = matrix.copy()
    for col in scaled.columns:
        span = scaled[col].max() - scaled[col].min()
        scaled[col] = 0 if span == 0 else (scaled[col] - scaled[col].min()) / span
    fig, ax = plt.subplots(figsize=(8, max(8, 0.22 * len(scaled))))
    im = ax.imshow(scaled.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(metrics)), metrics, rotation=45, ha="right")
    ax.set_yticks(range(len(scaled)), scaled.index, fontsize=5)
    ax.set_title("Limited VAF clonal complexity metrics")
    fig.colorbar(im, ax=ax, label="Scaled metric")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def plot_example_samples(assignments: pd.DataFrame, sample_summary: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    examples = (
        sample_summary.sort_values(["n_clusters", "n_input_mutations"], ascending=[False, False])
        .head(6)["sample_barcode"]
        .astype(str)
        .tolist()
    )
    fig, axes = plt.subplots(3, 2, figsize=(11, 12), squeeze=False)
    for ax, sample in zip(axes.ravel(), examples, strict=False):
        data = assignments[assignments["sample_barcode"].astype(str) == sample]
        if data.empty:
            ax.axis("off")
            continue
        ax.scatter(data["observed_vaf"], np.zeros(len(data)), c=pd.factorize(data["cluster_id"])[0], s=8, alpha=0.5)
        for cluster, group in data.groupby("cluster_id"):
            ax.axvline(group["cluster_mean_vaf"].iloc[0], linestyle="--", linewidth=1)
        ax.set_title(f"{sample}: {data['cluster_id'].nunique()} clusters", fontsize=9)
        ax.set_xlabel("Observed VAF")
        ax.set_yticks([])
    for ax in axes.ravel()[len(examples) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logging.info("Wrote %s", path)


def run_limited_vaf_clustering(paths: LimitedVafPaths, config: dict = CONFIG) -> dict[str, pd.DataFrame]:
    copy_neutral, selected_pilot, _readiness, _purity_ploidy = read_inputs(paths)
    assignments, sample_summary = run_clustering(copy_neutral, selected_pilot, config)
    project_summary = build_project_summary(sample_summary)
    qc = build_qc(selected_pilot, assignments, sample_summary, config)

    write_tsv_gz(assignments, paths.assignments)
    write_tsv(sample_summary, paths.sample_summary)
    write_tsv(project_summary, paths.project_summary)
    write_tsv(qc, paths.qc_summary)
    plot_density_by_sample(assignments, sample_summary, paths.density_figure)
    plot_cluster_counts_by_project(sample_summary, paths.cluster_counts_figure)
    plot_complexity_heatmap(sample_summary, paths.complexity_heatmap)
    plot_example_samples(assignments, sample_summary, paths.example_samples_figure)
    return {
        "assignments": assignments,
        "sample_summary": sample_summary,
        "project_summary": project_summary,
        "qc": qc,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run limited VAF-based clonal clustering prototype.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument("--max-clusters", type=int, default=CONFIG["max_clusters"])
    parser.add_argument("--min-cluster-size", type=int, default=CONFIG["min_cluster_size"])
    parser.add_argument("--min-input-mutations", type=int, default=CONFIG["min_input_mutations"])
    parser.add_argument("--expected-clonal-vaf-tolerance", type=float, default=CONFIG["expected_clonal_vaf_tolerance"])
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())
    config = dict(CONFIG)
    config["max_clusters"] = args.max_clusters
    config["min_cluster_size"] = args.min_cluster_size
    config["min_input_mutations"] = args.min_input_mutations
    config["expected_clonal_vaf_tolerance"] = args.expected_clonal_vaf_tolerance
    outputs = run_limited_vaf_clustering(default_paths(root), config)
    logging.info(
        "Limited VAF clustering complete: %d samples, %d clustered mutations",
        len(outputs["sample_summary"]),
        len(outputs["assignments"]),
    )


if __name__ == "__main__":
    main()

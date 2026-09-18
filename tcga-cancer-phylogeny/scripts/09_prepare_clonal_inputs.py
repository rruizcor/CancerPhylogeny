#!/usr/bin/env python

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from lib.common import configure_logging, ensure_dir, find_project_root, not_implemented, project_path, require_file


def select_high_information_samples(
    mutation_path: Path,
    purity_ploidy_path: Path,
    copy_number_dir: Path,
    output_path: Path,
) -> None:
    not_implemented(
        "High-information clonal pilot sample selection",
        "select high-TMB, high-purity samples with adequate depth, ref/alt counts, and local copy-number data",
    )


def prepare_pyclone_vi_input(sample_id: str, mutation_path: Path, copy_number_dir: Path, output_dir: Path) -> None:
    not_implemented(
        f"PyClone-VI input preparation for {sample_id}",
        "match variants to local copy-number segments and write per-mutation total/read-count/copy-number tables",
    )


def prepare_phylowgs_input(sample_id: str, mutation_path: Path, copy_number_dir: Path, output_dir: Path) -> None:
    not_implemented(
        f"PhyloWGS input preparation for {sample_id}",
        "write mutation and copy-number inputs only for samples passing QC thresholds",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare Level 3 clonal phylogeny pilot inputs.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())

    mutation_path = project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root=root)
    purity_ploidy_path = project_path("data", "processed", "features", "purity_ploidy_by_sample.tsv", root=root)
    copy_number_dir = project_path("data", "processed", "copy_number_segments", root=root)
    output_dir = ensure_dir(project_path("results", "clonal", "prepared_inputs", root=root))

    logging.warning(
        "Most TCGA cases have single bulk tumor samples; clonal clustering is usually more defensible than exact branching order."
    )
    require_file(mutation_path, "standardized mutation table")
    require_file(purity_ploidy_path, "purity/ploidy feature table")
    select_high_information_samples(mutation_path, purity_ploidy_path, copy_number_dir, output_dir / "pilot_samples.tsv")


if __name__ == "__main__":
    main()

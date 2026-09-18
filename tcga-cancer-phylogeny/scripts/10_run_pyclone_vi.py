#!/usr/bin/env python

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from lib.common import configure_logging, ensure_dir, find_project_root, not_implemented, project_path


def run_pyclone_vi_for_sample(sample_id: str, input_path: Path, output_dir: Path) -> None:
    not_implemented(
        f"PyClone-VI run for {sample_id}",
        "call pyclone-vi on prepared per-sample inputs and write cluster assignments plus cellular prevalence estimates",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run PyClone-VI for selected TCGA clonal pilot samples.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument("--sample-id", default=None, help="Optional single sample barcode to process.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())

    input_dir = project_path("results", "clonal", "prepared_inputs", "pyclone_vi", root=root)
    output_dir = ensure_dir(project_path("results", "clonal", "sample_cluster_assignments", root=root))

    logging.info("Planned PyClone-VI input directory: %s", input_dir)
    logging.info("Planned PyClone-VI output directory: %s", output_dir)
    sample_id = args.sample_id or "SAMPLE_PLACEHOLDER"
    run_pyclone_vi_for_sample(sample_id, input_dir / f"{sample_id}.tsv", output_dir)


if __name__ == "__main__":
    main()

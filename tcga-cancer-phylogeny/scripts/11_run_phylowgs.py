#!/usr/bin/env python

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from lib.common import configure_logging, ensure_dir, find_project_root, not_implemented, project_path


def run_phylowgs_for_sample(sample_id: str, input_dir: Path, output_dir: Path) -> None:
    not_implemented(
        f"PhyloWGS run for {sample_id}",
        "call a documented PhyloWGS installation only for samples with adequate mutation and copy-number inputs",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Optionally run PhyloWGS for selected TCGA clonal pilot samples.")
    parser.add_argument("--root", type=Path, default=None, help="Project root. Defaults to auto-detection.")
    parser.add_argument("--sample-id", default=None, help="Optional single sample barcode to process.")
    return parser.parse_args()


def main() -> None:
    configure_logging()
    args = parse_args()
    root = find_project_root(args.root) if args.root else find_project_root(Path(__file__).resolve())

    input_dir = project_path("results", "clonal", "prepared_inputs", "phylowgs", root=root)
    output_dir = ensure_dir(project_path("results", "clonal", "sample_tree_outputs", root=root))

    logging.warning(
        "PhyloWGS trees from single bulk TCGA samples should be reported as candidate structures, not definitive branching histories."
    )
    logging.info("Planned PhyloWGS input directory: %s", input_dir)
    logging.info("Planned PhyloWGS output directory: %s", output_dir)
    sample_id = args.sample_id or "SAMPLE_PLACEHOLDER"
    run_phylowgs_for_sample(sample_id, input_dir, output_dir)


if __name__ == "__main__":
    main()

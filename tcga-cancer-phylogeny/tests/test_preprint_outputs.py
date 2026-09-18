#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise AssertionError(message)


preprint = ROOT / "results" / "reports" / "tcga_cancer_phylogeny_preprint.md"
index = ROOT / "results" / "reports" / "tcga_cancer_phylogeny_preprint_figures_tables_index.tsv"

missing_or_empty = [str(path) for path in [preprint, index] if not path.exists() or path.stat().st_size == 0]
if missing_or_empty:
    fail(f"Missing or empty preprint outputs: {', '.join(missing_or_empty)}")

text = preprint.read_text(encoding="utf-8")

required_sections = [
    "## Introduction",
    "## Methods",
    "## Results",
    "## Discussion",
    "## Conclusion",
]
for section in required_sections:
    if section not in text:
        fail(f"Preprint manuscript is missing required section: {section}")

required_language = [
    "Level 1",
    "Level 2",
    "Level 3",
    "molecular similarity tree",
    "molecular trait dendrogram",
    "limited VAF-based clonal clustering",
    "97-feature",
    "Clinical WES/WTS",
    "CUP Interpretation Prototype",
    "descriptive molecular similarity weights, not tissue-of-origin probabilities",
]
for phrase in required_language:
    if phrase not in text:
        fail(f"Preprint manuscript is missing required language: {phrase}")

required_caveats = [
    "Level 1 and Level 2 are not literal species-like phylogenies",
    "Level 3 is not PyClone-VI/PhyloWGS and is not allele-specific CN-aware",
    "Segment means are not allele-specific integer copy number",
    "PyClone-VI and PhyloWGS remain paused",
    "not a clinically validated diagnostic classifier",
]
for phrase in required_caveats:
    if phrase not in text:
        fail(f"Preprint manuscript is missing required caveat: {phrase}")

if "70 retained features" in text or "70-feature" in text:
    fail("Preprint manuscript contains an outdated 70-feature atlas claim")

table_separator_count = text.count("|---")
if table_separator_count < 5:
    fail("Preprint manuscript should contain at least five Markdown tables")

if "## Figure Legends" not in text:
    fail("Preprint manuscript is missing figure legends")
for figure_number in range(1, 15):
    if f"Figure {figure_number}." not in text:
        fail(f"Preprint manuscript is missing Figure {figure_number} legend/reference")

if "## Table Legends" not in text:
    fail("Preprint manuscript is missing table legends")
for table_number in range(1, 9):
    if f"Table {table_number}." not in text:
        fail(f"Preprint manuscript is missing Table {table_number} legend/reference")

index_text = index.read_text(encoding="utf-8")
required_index_header = "item_type\titem_number\ttitle\tsource_file\tincluded_in_preprint\tnotes"
if not index_text.startswith(required_index_header):
    fail("Figures/tables index lacks the required header")

if "figure\t1\tConceptual overview of the multilevel framework\tpending\tTRUE" not in index_text:
    fail("Figures/tables index should record the pending conceptual figure placeholder")
if "table\t2\tTCGA reference atlas feature composition" not in index_text:
    fail("Figures/tables index should record the reference-atlas composition table")
if "table\t3\tSynthetic clinical projection example" not in index_text:
    fail("Figures/tables index should record the clinical projection table")
if "table\t4\tSynthetic CUP interpretation example" not in index_text:
    fail("Figures/tables index should record the CUP interpretation table")
if "table\t5\tPyClone-VI and PhyloWGS readiness" not in index_text:
    fail("Figures/tables index should record the clonal-readiness table")

print("Preprint manuscript output checks passed")

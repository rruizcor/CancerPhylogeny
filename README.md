# CancerPhylogeny

TCGA cancer molecular similarity, clonal analysis, and reference-atlas research project.

The complete working project is in [`tcga-cancer-phylogeny/`](tcga-cancer-phylogeny/). See its [detailed README](tcga-cancer-phylogeny/README.md) for methods, environment setup, workflow commands, completed analyses, and scientific limitations. Run project commands from that directory.

## Included work

- Analysis scripts, shared libraries, configuration, Snakemake workflow, notebooks, and environment specifications.
- Tests and available processed/intermediate datasets.
- Generated trees, figures, tables, clonal results, reference atlas, synthetic clinical example, validation and calibration results.
- Integrated research reports and preprint manuscript draft.
- Available source data files smaller than 100 MiB.

## Large source downloads

Four large source downloads (MC3 mutation files and the PanCanAtlas expression matrix) exceed GitHub's normal per-file limit and remain local. Their exact paths, byte sizes, and SHA-256 checksums are recorded in [`SOURCE_DATA_INVENTORY.tsv`](tcga-cancer-phylogeny/SOURCE_DATA_INVENTORY.tsv). The repository includes the resulting processed data and analysis outputs, but is not a byte-for-byte backup of those source downloads.

To restore source inputs, follow the detailed README's data acquisition instructions and the resource/download logic in `scripts/02_download_mc3_maf.R` and `scripts/04_download_expression.R`, inside the project directory. MC3 acquisition may require the documented manual source-file placement. The inventory preserves both root-level MC3 downloads and the cached raw MC3 file as distinct local files.

This is a snapshot of the available project files. Local caches and operating-system metadata are excluded. Research interpretation and validation limitations are documented in the project README and reports.

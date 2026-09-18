# tcga-cancer-phylogeny

Reproducible research project for analyzing TCGA cancer data with a three-level phylogenetic and evolutionary framework.

## Scientific framing

This project uses the language of phylogenetics carefully. Levels 1 and 2 are not literal species-like evolutionary trees and should not be interpreted as showing that one broad cancer type evolved from another. They are molecular similarity trees, trait dendrograms, or lineage-informed molecular phylograms that summarize shared molecular phenotypes across TCGA cancer types or cancer groups. Level 3 is closer to true tumor evolution because the unit of analysis is subclonal structure within individual tumors.

The guiding analogy is macroevolutionary trait mapping: build a tree-like representation of relatedness, then map malignant traits onto it. In this project, mapped traits may include tumor mutational burden, chromosomal instability, aneuploidy, whole-genome doubling, driver mutation prevalence, mismatch repair or MSI-related features, HRD/LOH-associated features, immune/stromal phenotypes, EMT/invasion programs, proliferation, lineage expression programs, clonal complexity, and possible therapy resistance features.

## Three analysis levels

### Level 1: Pan-cancer molecular similarity tree

Unit of analysis: TCGA cancer type/project.

Goal: build pan-cancer molecular similarity trees across available TCGA cancer types. Rows are TCGA projects and columns are aggregated molecular features. Candidate features include driver mutation prevalence, median or mean TMB, copy-number burden, aneuploidy, fraction genome altered, purity/ploidy summaries, pathway expression scores, immune/stromal features, proliferation, EMT, MSI/HRD-related features where available, and major histogenetic group.

Outputs include Gower-distance and correlation-distance trees, Newick files, tree figures, feature heatmaps, and cluster stability summaries.

### Level 2: Within-group molecular trees

Unit of analysis: cancer types, subtypes, or sample groups within broader cancer categories.

Separate molecular similarity trees are planned for carcinomas, sarcomas, hematolymphoid tumors, CNS/glial tumors, melanocytic tumors, and other biologically meaningful groups where feasible. Carcinoma analyses will test whether molecular clustering supports pan-squamous, pan-GI, pan-kidney, pan-gynecologic, and other lineage/histology-associated groupings. Sarcoma analyses will split TCGA-SARC by histologic subtype where metadata permits. Hematolymphoid analyses are expected to be limited because TCGA representation is mainly LAML and DLBC.

### Level 3: Within-patient clonal phylogeny pilot

Unit of analysis: inferred subclones within individual tumor samples.

This level is closest to a true tumor phylogeny. It uses somatic SNVs/indels, variant allele frequencies, ref/alt counts, local copy number, purity, and ploidy where available. Because most TCGA cases have a single bulk tumor sample, clonal clustering may be feasible, but exact branching order is often uncertain. The initial goal is a pilot set of high-information cases, not a full-TCGA clonal analysis.

Candidate pilot cases include high-TMB SKCM, high-TMB LUAD/LUSC, MSI-high or POLE-mutated UCEC/COAD if identifiable, and high-purity tumors with enough mutations and copy-number data.

## Repository layout

```text
config/                 Project lists, feature definitions, cancer groups, driver genes
data/raw/               Downloaded source data, ignored by git
data/interim/           Intermediate files, ignored by git
data/processed/         Analysis-ready derived files, ignored by git
scripts/                Workflow scripts and shared helpers
notebooks/              Quarto analysis notebooks
results/tables/         Final tables, ignored by git
results/trees/          Newick tree outputs, ignored by git
results/figures/        Figures, ignored by git
results/reports/        Rendered reports, ignored by git
results/clonal/         Level 3 pilot outputs, ignored by git
tests/                  Lightweight sanity checks
```

## Current status

Last updated: 2026-07-15.

This section is intended to be updated as the project progresses so a new contributor can quickly understand what is implemented, what has been validated, and what still needs attention.

Implemented and validated:

- Repository skeleton, configuration files, README, Snakefile, and starter notebooks.
- TCGA project metadata retrieval with all 33 configured TCGA projects represented.
- MC3 mutation layer with 3,592,393 mutation records, 10,201 tumor samples, 10,130 patients, all 33 projects represented, ref/alt counts present, `HGVSp_Short` present, and no configured driver genes absent from the MAF.
- Copy-number, purity, ploidy, and aneuploidy layer using PanCanAtlas/GDC ABSOLUTE and arm-call resources. Current validated coverage is 9,516 samples with purity/ploidy, 9,667 samples with aneuploidy/arm-level calls, all 33 projects represented, and no parsing warnings.
- RNA expression/pathway feature layer using the GDC PanCanAtlas `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv` matrix. Current validated coverage is 9,565 mutation-matched expression samples, 67 configured scoring genes, all 33 projects represented, and no parsing warnings.
- Level 1 project-level feature matrix assembly. The current `results/tables/level1_feature_matrix.tsv` has 33 project rows and 52 retained biologic feature columns spanning mutation-count proxies, driver mutation prevalence, purity/ploidy, aneuploidy/arm-level copy-number burden, and expression pathway scores. Two features, `median_fraction_genome_altered` and `median_copy_number_complexity_score`, are dropped for exceeding the configured 30% missingness threshold.
- Level 1 pan-cancer molecular similarity dendrograms from the validated feature matrix. The current tree step includes all 33 TCGA projects and all 52 retained biologic features, drops no additional features, imputes zero values, writes scaled features, Gower and correlation distance matrices, average-linkage Newick trees, a neighbor-joining sensitivity tree, tree/PCA/trait-heatmap figures, bootstrap co-clustering stability, and a tree QC summary.
- Level 1 interpretation and sensitivity analysis. `scripts/07b_interpret_level1_tree.R` cuts Gower and correlation trees at k = 3 through k = 8, writes cluster assignment and trait-summary tables, ranks feature associations with the Gower k = 6 clustering, runs feature-class and feature-weighting sensitivity analyses, interprets bootstrap co-clustering stability, generates summary figures, and writes `results/reports/level1_interpretation_summary.md`.
- Level 2 group-specific molecular similarity trees. `scripts/08_level2_within_group_trees.R` defines biologically meaningful project groups, filters features for within-group variability, scales features within group, builds Gower and correlation average-linkage trees where project counts are sufficient, records underpowered groups explicitly, writes group-specific matrices/distances/trees/figures/bootstrap tables, and writes `results/reports/level2_interpretation_summary.md`.
- Level 3 clonal-candidate feasibility and pilot selection. `scripts/09_select_level3_clonal_candidates.py` evaluates 10,201 mutation-layer tumor samples for ref/alt count availability, mutation-count proxy, purity, ploidy, aneuploidy/copy-number availability, and segment-level copy-number availability. The pre-GDC selection output had 3,017 samples passing minimum filters, 84 pilot-cohort samples, 11 pilot projects, zero segment-level copy-number candidates, and 3,017 `clonal_clustering_candidate_limited_cn` candidates.
- Level 3 segment-level copy-number acquisition/join step. `scripts/09b_download_segment_level_cn.py` checks local segment-level CN sources, parses supported TCGA/GDC segment formats when present, writes processed segment and segment-summary tables, joins segment availability back to Level 3 candidates, and updates candidate classes. Current validated run found no parseable local segment-level CN rows, attempted no remote download in the restricted shell, upgraded zero candidates, and preserved all 3,017 passing candidates as limited-CN clonal clustering candidates.
- Level 3 GDC pilot segment-level copy-number downloader. `scripts/09c_fetch_gdc_segment_cn.py` queries the GDC files API for open-access TCGA Copy Number Variation files with data types `Masked Copy Number Segment`, `Copy Number Segment`, and `Filtered Copy Number Segment`, matches files to the Level 3 pilot cohort, selects one best segment file per pilot sample, downloads the full pilot sample-level match set by default, parses GDC segment tables, and updates the same Level 3 segment-CN outputs. Current validated run found 20,274 GDC segment-CN files across the pilot projects, matched 434 files to pilot samples or patients, selected 84 best sample-level files, downloaded 74 new files while reusing 10 existing files, parsed 84 files, processed 34,028 segment rows, marked all 84 pilot samples as having segment-level CN, and upgraded 84 candidates to `copy_number_aware_clonal_phylogeny_candidate`.
- Level 3 mutation local copy-number annotation. `scripts/10_annotate_mutations_with_local_cn.py` joins somatic mutations from the 84 copy-number-aware pilot samples to overlapping sample-level segment-CN intervals, computes depth and observed VAF, records conservative local segment status, writes a mutation-level annotated table, sample-level annotation summary, non-final PyClone-style preview, QC table, and figures. Current validated run evaluated 220,124 pilot mutations, annotated 213,112 with local segment CN (96.815%), produced 131,398 preview rows passing default mutation/depth/local-CN filters, and flagged all 84 pilot samples as eligible for downstream clonal-input preparation. Major/minor allele-specific CN is not inferred; preview rows keep `major_cn` and `minor_cn` as `NA`.
- Level 3 clonal input stratification and readiness assessment. `scripts/11_prepare_clonal_input_sets.py` stratifies the mutation/local-CN annotations into analytically honest input sets: a copy-neutral/near-diploid limited-VAF clustering candidate set, a broader segment-annotated set that retains local segment mean without claiming allele-specific CN, sample-level readiness classes, a 30-sample limited-VAF pilot selection, QC, and figures. Current validated run produced 106,337 copy-neutral limited-VAF candidate mutations, 131,398 segment-annotated candidate mutations, 79 samples ready for limited VAF clustering, zero samples eligible for fully copy-number-aware clustering, 84 samples requiring allele-specific CN for copy-number-aware inputs, and a 30-sample selected prototype pilot spanning all 11 pilot projects.
- Level 3 limited VAF-based clonal clustering prototype. `scripts/12_limited_vaf_clonal_clustering.py` clusters copy-neutral candidate mutation VAFs independently within the 30 selected prototype samples using per-sample Gaussian mixture models selected by BIC, with a quantile fallback if needed. This is explicitly a limited VAF clustering prototype, not PyClone-VI, PhyloWGS, a full copy-number-aware analysis, or definitive phylogenetic branching. Current validated run clustered all 30 selected samples and 81,154 mutations, skipped zero samples, and classified samples as 5 predominantly clonal-like, 11 oligoclonal-like, and 14 multicluster subclonal-like.
- Level 3 allele-specific CN readiness and acquisition planning. `scripts/12b_allele_specific_cn_readiness.py` searches supported local FACETS, Sequenza, ASCAT, ABSOLUTE-style, and generic major/minor CN outputs; rejects segment-mean-only files as allele-specific evidence; standardizes valid integer total/major/minor states; measures mutation overlap; and writes an external-inference manifest without running PyClone-VI or PhyloWGS. The current validated run evaluated all 84 segment-CN pilot samples and all 30 limited-VAF pilot samples, found zero local allele-specific CN files, found zero samples with major/minor integer CN, and marked zero samples ready for PyClone-VI or PhyloWGS. All 84 currently require external allele-specific CN inference. The recommended next action is to obtain paired tumor/normal BAMs or validated precomputed SNP allele-count inputs for the 30 limited-VAF pilot samples and run one method consistently, with FACETS as the default paired-WES planning option.
- Hardened TCGA sample-level molecular reference atlas. `scripts/20_build_tcga_sample_reference_atlas.py` locks deterministic tumor-sample and feature ordering; preserves an observed pre-imputation matrix; records feature missingness, definitions, medians, means, standard deviations, and contract hashes; builds PCA coordinates/model/variance outputs; and saves Euclidean and cosine nearest-neighbor models with Euclidean recorded as the default. The validated atlas includes 10,201 TCGA tumor samples from 10,130 patients, all 33 TCGA projects, 97 raw molecular features, 97 retained features, zero dropped features at the 30% missingness threshold, and 23,222 median-imputed values in the modeling matrix. The retained set contains 2 mutation-count proxies, 29 driver-gene mutation indicators, 29 driver-gene mutation counts, 2 purity/ploidy features, 4 aneuploidy/arm-level copy-number features, 11 expression pathway scores, and 20 expression PCs. Four atlas QC figures and reusable PCA/nearest-neighbor metadata are generated. This unsupervised reference is consumed by the current clinical projection prototype and can support future harmonized WES/WTS and CUP research; it is not a clinical classifier. The separate PyClone-VI/PhyloWGS branch remains paused until validated allele-specific integer copy number is available.
- Clinical WES/WTS query projection prototype. `scripts/21_project_clinical_case_to_tcga.py` ingests already-computed clinical feature tables, matches exact locked atlas feature names, applies saved TCGA medians/means/standard deviations, reuses the saved PCA and nearest-neighbor models when compatible, and falls back to direct distance calculation from the scaled atlas when needed. The artificial `example_cup_case` supplies 77 of 97 locked features (79.381% coverage), median-imputes 20 expression-PC features, and returns 50 Euclidean neighbors. Its current descriptive project summary contains COAD (17 neighbors), PAAD (16), READ (11), STAD (3), and LUAD (3); all 50 neighbors are in the broad Carcinoma group. These artificial results demonstrate workflow behavior only and are not a tissue-of-origin prediction or diagnosis.
- CUP molecular interpretation prototype. `scripts/22_cup_nearest_neighbor_classifier.py` consumes the locked clinical projection, validates project and major-group summaries against the top-k neighbor rows, aggregates all 33 TCGA projects into 18 cautious lineage groups, computes overall and modality-specific ambiguity metrics and transparent confidence tiers, evaluates supplied-only heuristic feature evidence, constrains review to a submitted differential diagnosis when present, and writes a heuristic molecular-pathologic concordance table. For the artificial `example_cup_case`, lower_GI receives 0.5701 of normalized top-k similarity weight and pancreatobiliary/hepatobiliary receives 0.3092; the result is `moderate_ambiguity` with `moderate_support` for the top lineage. The workflow reports competing APC/KRAS and KRAS/SMAD4 heuristic patterns rather than resolving them as a diagnosis. This is a research interpretation framework, not a clinically validated tissue-of-origin classifier.
- Known-primary internal validation. `scripts/23_validate_known_primary_projection.py` treats all 10,201 atlas samples from 33 projects as TCGA pseudo-unknowns, excludes each query from its own reference neighbors, and summarizes 50 Euclidean neighbors using the same project, lineage, major-group, and ambiguity logic as the CUP prototype. The full leave-one-out run recovered the known project at top 1 for 59.455% of samples, top 3 for 82.364%, and top 5 for 89.550%; top-1 configured-lineage recovery was 66.503%, top-3 lineage recovery was 88.668%, and broad major-group recovery was 88.246%. Low-ambiguity samples had 82.998% top-1 project recovery versus 34.766% in high-ambiguity samples, supporting the descriptive utility of the current abstention signal without establishing probability calibration. This is internal TCGA validation with shared cohort and preprocessing, not external clinical validation.
- Calibration and feature-ablation validation. `scripts/24_calibrate_projection_confidence_and_feature_ablation.py` derives transparent project-, lineage-, and major-group confidence thresholds on a deterministic 5,092-sample within-project calibration split and evaluates the resulting tiers on the remaining 5,109 internal samples. Held-out top-1 recovery for low-, moderate-, and high-confidence tiers was 35.762%, 55.763%, and 84.837% at project level; 42.322%, 61.189%, and 92.719% at lineage level; and 65.482%, 94.875%, and 99.759% at major-group level. Leave-one-out ablation across 13 prespecified feature sets found 69.111% top-1 project recovery for the 31-feature WTS-like set, 38.535% for the 66-feature WES-like set, and 59.455% for the 97-feature combined set. These tiers are validation-informed empirical reliability labels, not probabilities or externally validated clinical cutoffs; expression dependence and incomplete-modality coverage must be reported explicitly.
- Refreshed integrated project reporting layer. `scripts/13_generate_integrated_project_report.py` now summarizes the validated Level 1, Level 2, Level 3, allele-specific-CN readiness, hardened 97-feature atlas, synthetic clinical projection, and CUP interpretation outputs without rerunning primary analyses. The main report, executive summary, methods, limitations/next-steps document, optional HTML, machine-readable summary, and output index all use the current atlas and clinical/CUP QC tables.
- Refreshed preprint-style manuscript drafting layer. `scripts/14_generate_preprint_manuscript.py` now incorporates the 97-feature atlas, synthetic WES/WTS projection, uncertainty-aware CUP interpretation, and current PyClone-VI/PhyloWGS readiness boundary into the abstract, methods, results, discussion, limitations, tables, and figure legends. The draft remains explicit that clinical similarity weights are not diagnostic probabilities and that the example case is synthetic and unvalidated.

Current final proof-of-concept status: the project has a complete validated path from TCGA data acquisition and feature engineering through Level 1 pan-cancer molecular similarity trees, Level 2 group-specific molecular similarity trees, Level 3 limited VAF-based clonal-structure prototyping, an explicit allele-specific CN readiness gate before copy-number-aware clonal inference, a hardened TCGA sample-level molecular reference atlas, a research-only clinical WES/WTS query projection prototype, a cautious CUP molecular interpretation layer, internal known-primary pseudo-unknown validation, validation-informed confidence calibration and feature-ablation benchmarking, integrated reporting, and a first preprint-style manuscript draft.

Current Level 1 interpretation highlights:

- Gower k = 6 separates a broad lower-mutation or lineage-distinct cluster including kidney, CNS, endocrine/adrenal, liver/biliary, mesothelioma, prostate, thyroid, thymoma, uveal melanoma, and sarcoma projects; a broad carcinoma-enriched cluster including BLCA, BRCA, CESC, COAD/READ, ESCA/HNSC, LUAD/LUSC, OV, and STAD; and singleton DLBC, SKCM, UCEC, and UCS clusters.
- Correlation k = 6 differs because it emphasizes relative feature-profile shape rather than absolute scaled trait levels. It groups ACC/CHOL/kidney/CNS/liver/mesothelioma/endocrine/prostate/uveal melanoma together, groups several GI/gynecologic carcinoma projects together, and creates separate BRCA/PAAD/SARC, DLBC/LAML/THYM, ESCA/HNSC/LUSC/OV/TGCT/UCS, and GBM/LUAD/SKCM clusters.
- The strongest Gower k = 6 feature associations are currently mutation-heavy, including BRCA1/BRCA2 mutation prevalence, mutation-count proxies, ALK, ERBB2, and related driver-prevalence features. These are descriptive cluster associations, not causal claims.
- Feature-class sensitivity suggests mutation features dominate the all-feature Gower structure more than copy-number-only or expression-only runs. Mutation plus expression is closest to the full tree among non-reference subsets; equal class weighting changes more project assignments, while copy-number and expression down-weighting preserve much of the reference structure.
- Bootstrap co-clustering supports stable project pairs and groups such as COAD/READ and a low-mutation cluster containing combinations of MESO, PCPG, PRAD, THCA, and UVM. Projects with fewer stable partners are flagged in `results/tables/level1_bootstrap_project_stability.tsv`.

Current Level 2 interpretation highlights:

- Project-level trees were generated for `carcinoma_all` (20 projects), `pan_squamous` (5), `gi_pancancreatobiliary` (7), `kidney` (3), and `gynecologic_breast` (5).
- `cns_glial` (GBM/LGG), `melanocytic` (SKCM/UVM), `hematolymphoid` (LAML/DLBC), and `sarcoma` (SARC only) are recorded as underpowered for robust project-level tree construction. No SARC subtype-level molecular feature matrix is available in the current workflow.
- The carcinoma-all Gower tree groups a broad carcinoma cluster containing BLCA, BRCA, CESC, ESCA, HNSC, LUAD, LUSC, OV, and STAD; a CHOL/KICH/KIRP/LIHC/PRAD/THCA cluster; a COAD/READ pair; and separate KIRC, PAAD, and UCEC branches at the current k = 6 summary.
- Expected subgroup signals appear cautiously: HNSC/ESCA cluster within the squamous-enriched set; COAD/READ/STAD and CHOL/LIHC cluster within the GI/hepatobiliary set; KIRC/KIRP cluster separately from KICH in kidney; and BRCA/OV cluster within the gynecologic/breast set, with UCEC, CESC, and UCS separate.
- Gower and correlation trees broadly agree on some close relationships such as COAD/READ and CHOL/KICH/KIRP/LIHC, but differ in exact placement because correlation emphasizes relative feature-profile shape rather than absolute scaled trait levels.

Current Level 3 candidate-selection highlights:

- Total samples evaluated: 10,201. Ref/alt counts are available for all 10,201 mutation-layer samples; purity and ploidy are available for 9,516; copy-number or aneuploidy features are available for 9,667.
- Samples passing minimum filters: 3,017. The GDC full-pilot segment-CN fetch currently upgrades 84 pilot samples to copy-number-aware clonal phylogeny candidates; the remaining passing candidates outside the pilot set remain limited-CN clonal clustering candidates until segment CN is acquired for them.
- Pilot cohort size: 84 samples across BLCA, COAD, DLBC, HNSC, LUAD, LUSC, READ, SARC, SKCM, UCEC, and UCS.
- Top projects by number of candidates include LUSC, SKCM, LUAD, BLCA, HNSC, STAD, COAD, UCEC, LIHC, BRCA, CESC, OV, ESCA, READ, KIRP, GBM, DLBC, and SARC.
- Main exclusion reasons are too few mutations, low purity, missing purity/ploidy, and missing copy-number data. The current mutation burden field remains a mutation-count proxy, not true TMB.
- Segment-level CN join result: the local-only scan found no parseable local source, but the GDC pilot fetcher now selects and parses one sample-level segment file for every pilot sample with a sample-level GDC match. `data/processed/copy_number/segments_by_sample.tsv.gz` currently contains 34,028 segment rows for all 84 pilot samples. `results/tables/gdc_segment_cn_manifest_pilot.tsv` records 20,274 files found and 434 pilot-matched rows; `results/tables/gdc_segment_cn_best_file_per_pilot_sample.tsv` records 84 selected best files; `results/tables/gdc_segment_cn_download_log.tsv` records 74 new downloads and 10 reused local files. No pilot samples currently lack segment-level CN.
- Mutation local-CN annotation result: all 84 pilot samples have mutation-to-segment annotations. `data/processed/clonal/level3_mutations_with_local_cn.tsv.gz` contains 220,124 pilot mutations, 213,112 local-CN matches, depth/VAF fields, purity/ploidy, segment coordinates, segment means, conservative gain/loss/neutral status, and match notes. `results/tables/level3_mutation_local_cn_summary_by_sample.tsv` marks all 84 samples eligible under the current thresholds. `data/processed/clonal/level3_pyclone_input_preview.tsv` contains 131,398 candidate preview rows, but it is not final PyClone-VI input because the segment data are not allele-specific integer copy number.
- Clonal input preparation result: `data/processed/clonal/level3_copy_neutral_vaf_clustering_input.tsv.gz` contains 106,337 candidate mutations in near-neutral segments for limited VAF-based clustering prototypes. `data/processed/clonal/level3_segment_annotated_clonal_input.tsv.gz` contains 131,398 broader segment-annotated candidate mutations for visualization and future CN-aware preparation. `results/tables/level3_clonal_input_readiness_by_sample.tsv` classifies 79 samples as `ready_for_limited_vaf_clustering` and 5 as `insufficient_neutral_mutations`; all 84 remain not copy-number-aware because allele-specific major/minor CN is unavailable. `results/tables/level3_limited_vaf_clustering_pilot_samples.tsv` selects 30 prototype samples with a default cap of 3 per project.
- Limited VAF clustering result: `results/clonal/limited_vaf/level3_limited_vaf_cluster_assignments.tsv.gz` contains 81,154 clustered copy-neutral mutations from the 30 selected prototype samples. `results/tables/level3_limited_vaf_clonal_complexity_by_sample.tsv` summarizes per-sample cluster counts, dominant cluster fraction, subclonal fraction, entropy, dispersion, and cautious interpretation class. `results/tables/level3_limited_vaf_clonal_complexity_by_project.tsv` summarizes project-level complexity. `results/tables/level3_limited_vaf_clustering_qc_summary.tsv` records the method, parameters, skipped samples, interpretation class counts, and limitations.

What worked:

- Public GDC/PanCanAtlas files for copy-number and expression were retrievable and are cached under `data/raw/`.
- Local-file fallbacks and explicit missing-feature placeholders are implemented for mutation, copy-number, and expression layers.
- Barcode harmonization now aligns expression, mutation, purity/ploidy, and aneuploidy sample outputs to the mutation-layer sample universe.
- Existing layer tests pass for mutation, copy-number, and expression outputs.
- The Level 1 feature matrix keeps sample-support counts in `results/tables/level1_feature_support.tsv`, separate from biologic features, so sample availability does not silently become a clustering feature.
- The Level 1 tree tests pass for synthetic feature handling and generated output contracts:
  `Rscript tests/test_level1_tree_synthetic.R` and `Rscript tests/test_level1_tree_outputs.R`.
- The Level 1 interpretation tests pass for synthetic sensitivity/ARI checks and generated output contracts:
  `Rscript tests/test_level1_interpretation_synthetic.R` and `Rscript tests/test_level1_interpretation_outputs.R`.
- The Level 2 tests pass for synthetic group/filtering/underpowered checks and generated output contracts:
  `Rscript tests/test_level2_tree_synthetic.R` and `Rscript tests/test_level2_tree_outputs.R`.
- The Level 3 candidate-selection tests pass for synthetic scoring/filtering/pilot-cohort behavior and generated output contracts:
  `python tests/test_level3_candidate_selection_synthetic.py` and `python tests/test_level3_candidate_selection_outputs.py`.
- The Level 3 segment-CN tests pass for synthetic segment parser/barcode/upgrade behavior and generated output contracts:
  `python tests/test_level3_segment_cn_synthetic.py` and `python tests/test_level3_segment_cn_outputs.py`.
- The GDC pilot segment-CN tests pass for synthetic GDC matching/parser behavior and generated output contracts:
  `python tests/test_gdc_segment_cn_fetch_synthetic.py` and `python tests/test_gdc_segment_cn_outputs.py`.
- The Level 3 mutation local-CN annotation tests pass for synthetic interval/barcode/filtering behavior and generated output contracts:
  `python tests/test_mutation_local_cn_annotation_synthetic.py` and `python tests/test_mutation_local_cn_annotation_outputs.py`.
- The Level 3 clonal input preparation tests pass for synthetic filtering/readiness/pilot-selection behavior and generated output contracts:
  `python tests/test_clonal_input_preparation_synthetic.py` and `python tests/test_clonal_input_preparation_outputs.py`.
- The Level 3 limited VAF clustering tests pass for synthetic one-/two-cluster GMM behavior, low-confidence handling, and generated output contracts:
  `python tests/test_limited_vaf_clonal_clustering_synthetic.py` and `python tests/test_limited_vaf_clonal_clustering_outputs.py`.
- The TCGA sample-level reference atlas tests pass for synthetic join/preprocessing/artifact behavior and generated output contracts:
  `python tests/test_tcga_sample_reference_atlas_synthetic.py` and `python tests/test_tcga_sample_reference_atlas_outputs.py`.
- The known-primary projection validation tests pass for self/same-patient exclusion, project ranking, lineage mapping, ambiguity summaries, confusion matrices, and generated output contracts:
  `python tests/test_known_primary_validation_synthetic.py` and `python tests/test_known_primary_validation_outputs.py`.
- The calibration and feature-ablation tests pass for feature grouping, exact self exclusion, absent-class ranking, threshold orientation, confidence-tier separation, small-matrix ablation, and generated output contracts:
  `python tests/test_calibration_feature_ablation_synthetic.py` and `python tests/test_calibration_feature_ablation_outputs.py`.
- The integrated report tests pass for generated report contracts, required Level 1/2/3 sections, output-index existence flags, and avoiding false full copy-number-aware phylogeny claims:
  `python tests/test_integrated_project_report_outputs.py`.

Issues encountered and fixes:

- `config/driver_genes.csv` originally had an unquoted comma in the APC rationale field. The row was quoted correctly and driver-gene CSV validation was strengthened.
- R `arrow` was not available in the original mutation-layer environment. The mutation layer uses a Python `pyarrow` fallback for Parquet output, and later layers avoid strict R `arrow` dependency.
- The first copy-number layer implementation produced missing-feature placeholders when no local files were present. The script was extended to retrieve the real PanCanAtlas ABSOLUTE purity/ploidy and aneuploidy/arm-call files.
- The first expression download attempt failed under restricted network access. A network-approved run downloaded the 1.8 GB GDC PanCanAtlas expression matrix successfully, and future runs reuse the cached local file.
- Initial expression outputs had a small set of barcode mismatches caused by vial suffix differences such as expression `03A` versus MC3 `03B`. The harmonizer now remaps by unique 15-character TCGA sample prefix when unambiguous; final expression-to-mutation sample mismatch count is zero.

Known limitations:

- `tmb_by_sample.tsv` is currently a mutation-count proxy, not true mutations/Mb, because callable territory has not been added.
- Expression pathway scores are simple average z-scored gene-set signatures and should be treated as first-pass phenotype summaries.
- Expression PCA is computed from the configured expression-gene subset, not the full transcriptome.
- Snakemake rules are wired, but `snakemake` was not installed in the active shell used for validation, so workflow dry-runs could not be executed there.
- Level 1 tree inputs are now scaled in `scripts/07_level1_pan_cancer_tree.R`; the original `level1_feature_matrix.tsv` remains raw.
- Level 1 trees are pan-cancer molecular similarity dendrograms or molecular trait trees, not literal species-like evolutionary phylogenies.
- Level 1 interpretation is sensitive to feature balance: mutation features currently have the strongest cluster associations, and equal class weighting changes several cluster assignments. Interpret exact branch placement and singleton clusters cautiously.
- Level 2 group-specific trees remain project-level molecular similarity analyses. Small groups, especially kidney with three projects and the underpowered two-project groups, should be interpreted as descriptive comparisons rather than robust within-lineage trees.
- Level 3 candidate selection, pilot segment-CN acquisition, and mutation-to-local-CN annotation are implemented, but full clonal phylogeny is not. The 84 pilot samples now have sample-level segment CN and mutation overlap annotations, but final PyClone-VI/PhyloWGS input generation and tree inference have not been run. Serial, relapse, metastasis, or multi-region samples remain preferred for robust branching inference.
- Segment-level CN acquisition now has a GDC-backed full-pilot downloader. The default `DOWNLOAD_MODE = "pilot_full"` selects one best sample-level file per pilot sample, skips existing files unless `OVERWRITE_EXISTING = True`, and leaves patient-level matches as lower-confidence evidence that does not auto-upgrade candidates by default.
- Mutation local-CN annotation is an interval overlap/preparation step, not clonal inference. It uses segment mean/log2-like values conservatively and does not transform them into allele-specific integer major/minor copy number. The PyClone-style preview is explicitly non-final until copy-number model assumptions or allele-specific CN sources are added.
- Clonal input preparation deliberately separates limited VAF-based clustering candidates from segment-annotated, not-copy-number-aware candidates. The current readiness assessment reports zero fully copy-number-aware samples because major/minor allele-specific CN is absent.
- Limited VAF clustering is a prototype summary of VAF cluster structure in copy-neutral regions only. It is not a definitive clonal phylogeny, does not infer branching order, and does not replace copy-number-aware methods that require allele-specific integer CN. Single bulk TCGA samples further limit phylogenetic interpretation.
- The sample-level reference atlas is a translational reference resource, not a supervised classifier and not a clinical case ingestion pipeline. Future clinical WES/WTS cases must be transformed using `tcga_reference_feature_definitions.yaml` and `tcga_reference_scaler_parameters.json`; otherwise nearest-neighbor distances or PCA projections will not be comparable to the TCGA reference.
- Confidence calibration and feature ablation remain internal to TCGA. The deterministic calibration/evaluation split reduces direct threshold-fit reporting bias but shares cohort selection, assays, feature engineering, global imputation, and scaling. WTS-like outperformance does not establish external clinical performance, and the ablation does not directly simulate incomplete WTS, platform effects, or clinical missingness mechanisms.
- Running Python plotting in the current sandbox emits fontconfig cache warnings, but the Level 3 PDF figures are created successfully.

Current next step:

- Validate the new confidence thresholds and modality benchmarks against independent known-primary WES/WTS cases, simulate realistic partial-modality coverage, and review assay harmonization before any CUP classifier training. For Level 3, add allele-specific CN resources before final PyClone-VI/PhyloWGS-style copy-number-aware inputs.

## Planned data products

Key planned outputs include:

- `data/interim/tcga_projects.tsv`
- `data/processed/mutations/mc3_somatic_mutations.parquet`
- `data/processed/features/tmb_by_sample.tsv`
- `data/processed/features/mutation_prevalence_by_project.tsv`
- `data/processed/features/aneuploidy_by_sample.tsv`
- `data/processed/features/aneuploidy_by_project.tsv`
- `data/processed/features/purity_ploidy_by_sample.tsv`
- `data/processed/features/expression_pathway_scores_by_sample.tsv`
- `data/processed/features/expression_pathway_scores_by_project.tsv`
- `data/processed/features/expression_pca_by_sample.tsv`
- `data/processed/features/expression_pca_variance.tsv`
- `data/processed/expression/expression_matrix_by_sample.tsv.gz`
- `results/tables/level1_feature_matrix.tsv`
- `results/tables/level1_feature_support.tsv`
- `results/tables/level1_feature_matrix_qc_summary.tsv`
- `results/tables/level1_feature_matrix_scaled.tsv`
- `results/tables/level1_gower_distance_matrix.tsv`
- `results/tables/level1_correlation_distance_matrix.tsv`
- `results/tables/level1_bootstrap_cluster_stability.tsv`
- `results/tables/level1_tree_qc_summary.tsv`
- `results/trees/level1_pan_cancer_gower_hclust_tree.nwk`
- `results/trees/level1_pan_cancer_correlation_hclust_tree.nwk`
- `results/trees/level1_pan_cancer_neighbor_joining_tree.nwk`
- `results/figures/level1_pan_cancer_gower_tree.pdf`
- `results/figures/level1_pan_cancer_correlation_tree.pdf`
- `results/figures/level1_pan_cancer_tree_with_trait_heatmap.pdf`
- `results/figures/level1_pan_cancer_feature_pca.pdf`
- `results/tables/level1_gower_cluster_assignments.tsv`
- `results/tables/level1_correlation_cluster_assignments.tsv`
- `results/tables/level1_gower_cluster_summaries.tsv`
- `results/tables/level1_correlation_cluster_summaries.tsv`
- `results/tables/level1_feature_pca_loadings.tsv`
- `results/tables/level1_feature_variance_summary.tsv`
- `results/tables/level1_feature_cluster_association.tsv`
- `results/tables/level1_feature_class_sensitivity_summary.tsv`
- `results/tables/level1_feature_class_cluster_assignments.tsv`
- `results/tables/level1_feature_weighting_sensitivity_summary.tsv`
- `results/tables/level1_weighted_cluster_assignments.tsv`
- `results/tables/level1_stable_project_pairs.tsv`
- `results/tables/level1_bootstrap_project_stability.tsv`
- `results/reports/level1_interpretation_summary.md`
- `results/figures/level1_cluster_trait_summary_heatmap.pdf`
- `results/figures/level1_feature_driver_barplot.pdf`
- `results/figures/level1_bootstrap_coclustering_heatmap.pdf`
- `results/figures/level1_feature_class_sensitivity_heatmap.pdf`
- `results/figures/level1_weighting_sensitivity_pca_or_mds.pdf`
- `results/tables/level2_group_definitions.tsv`
- `results/tables/level2_tree_qc_summary.tsv`
- `results/reports/level2_interpretation_summary.md`
- `results/tables/level2_<group_name>_feature_matrix.tsv`
- `results/tables/level2_<group_name>_feature_matrix_scaled.tsv`
- `results/tables/level2_<group_name>_gower_distance_matrix.tsv`
- `results/tables/level2_<group_name>_correlation_distance_matrix.tsv`
- `results/tables/level2_<group_name>_bootstrap_cluster_stability.tsv`
- `results/trees/level2_<group_name>_gower_hclust_tree.nwk`
- `results/trees/level2_<group_name>_correlation_hclust_tree.nwk`
- `results/figures/level2_<group_name>_gower_tree.pdf`
- `results/figures/level2_<group_name>_correlation_tree.pdf`
- `results/figures/level2_<group_name>_trait_heatmap.pdf`
- `results/figures/level2_<group_name>_feature_pca.pdf`
- `results/tables/level3_candidate_samples.tsv`
- `results/tables/level3_candidate_summary_by_project.tsv`
- `results/tables/level3_pilot_cohort.tsv`
- `results/tables/level3_excluded_samples.tsv`
- `results/tables/level3_candidate_selection_qc_summary.tsv`
- `results/figures/level3_candidate_mutation_count_vs_purity.pdf`
- `results/figures/level3_candidate_counts_by_project.pdf`
- `results/figures/level3_candidate_score_by_project.pdf`
- `results/figures/level3_pilot_cohort_overview.pdf`
- `data/processed/copy_number/segments_by_sample.tsv.gz`
- `data/processed/copy_number/segment_level_cn_summary_by_sample.tsv`
- `results/tables/gdc_segment_cn_manifest_pilot.tsv`
- `results/tables/gdc_segment_cn_best_file_per_pilot_sample.tsv`
- `results/tables/gdc_segment_cn_download_log.tsv`
- `results/tables/level3_candidate_samples_with_segments.tsv`
- `results/tables/level3_pilot_cohort_with_segments.tsv`
- `results/tables/level3_segment_cn_qc_summary.tsv`
- `results/tables/level3_segment_cn_availability_by_project.tsv`
- `results/figures/level3_segment_cn_availability_by_project.pdf`
- `results/figures/level3_segment_cn_pilot_coverage.pdf`
- `results/figures/level3_segment_cn_upgrade_by_project.pdf`
- `results/figures/level3_segment_count_distribution.pdf`
- `results/figures/level3_candidate_upgrade_status.pdf`
- `data/processed/clonal/level3_mutations_with_local_cn.tsv.gz`
- `data/processed/clonal/level3_pyclone_input_preview.tsv`
- `results/tables/level3_mutation_local_cn_summary_by_sample.tsv`
- `results/tables/level3_mutation_local_cn_annotation_qc_summary.tsv`
- `results/figures/level3_mutation_cn_annotation_coverage_by_sample.pdf`
- `results/figures/level3_mutation_vaf_vs_local_segment_mean.pdf`
- `results/figures/level3_mutation_depth_distribution.pdf`
- `data/processed/clonal/level3_copy_neutral_vaf_clustering_input.tsv.gz`
- `data/processed/clonal/level3_segment_annotated_clonal_input.tsv.gz`
- `results/tables/level3_clonal_input_readiness_by_sample.tsv`
- `results/tables/level3_limited_vaf_clustering_pilot_samples.tsv`
- `results/tables/level3_clonal_input_preparation_qc_summary.tsv`
- `results/figures/level3_copy_neutral_candidate_counts_by_sample.pdf`
- `results/figures/level3_copy_neutral_candidate_counts_by_project.pdf`
- `results/figures/level3_vaf_distribution_copy_neutral_candidates.pdf`
- `results/figures/level3_readiness_class_by_project.pdf`
- `results/clonal/limited_vaf/level3_limited_vaf_cluster_assignments.tsv.gz`
- `results/tables/level3_limited_vaf_clonal_complexity_by_sample.tsv`
- `results/tables/level3_limited_vaf_clonal_complexity_by_project.tsv`
- `results/tables/level3_limited_vaf_clustering_qc_summary.tsv`
- `results/figures/level3_limited_vaf_cluster_vaf_density_by_sample.pdf`
- `results/figures/level3_limited_vaf_cluster_counts_by_project.pdf`
- `results/figures/level3_limited_vaf_complexity_heatmap.pdf`
- `results/figures/level3_limited_vaf_example_samples.pdf`
- `results/reference_atlas/tcga_sample_reference_matrix.tsv.gz`
- `results/reference_atlas/tcga_sample_reference_metadata.tsv`
- `results/reference_atlas/tcga_sample_reference_scaled_matrix.tsv.gz`
- `results/reference_atlas/tcga_reference_scaler_parameters.json`
- `results/reference_atlas/tcga_reference_feature_definitions.yaml`
- `results/reference_atlas/tcga_reference_pca_model.pkl`
- `results/reference_atlas/tcga_reference_pca_coordinates.tsv`
- `results/reference_atlas/tcga_reference_nearest_neighbor_index.pkl`
- `results/tables/tcga_sample_reference_atlas_qc_summary.tsv`
- `results/figures/tcga_sample_reference_pca_by_project.pdf`
- `results/figures/tcga_sample_reference_pca_by_major_group.pdf`
- `results/reports/integrated_tcga_cancer_phylogeny_report.md`
- `results/reports/integrated_tcga_cancer_phylogeny_executive_summary.md`
- `results/reports/integrated_tcga_cancer_phylogeny_methods.md`
- `results/reports/integrated_tcga_cancer_phylogeny_limitations_next_steps.md`
- `results/reports/integrated_tcga_cancer_phylogeny_report.html`
- `results/tables/integrated_project_summary.tsv`
- `results/tables/integrated_project_outputs_index.tsv`
- `results/clonal/sample_cluster_assignments/`
- `results/clonal/sample_tree_outputs/`
- `results/tables/clonal_summary_by_sample.tsv`

## Reproducible workflow

The project uses R for GDC/TCGA querying, feature engineering, statistics, tree building, and visualization; Python for clonal-phylogeny preprocessing and optional PyClone-VI/PhyloWGS integration; Snakemake for workflow orchestration; and Quarto for notebooks/reports.

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate tcga-cancer-phylogeny
```

Run configuration sanity checks:

```bash
Rscript tests/test_config.R
```

Dry-run the workflow:

```bash
snakemake -n
```

Fetch GDC project metadata:

```bash
Rscript scripts/00_fetch_tcga_projects.R
```

Build the MC3 mutation data layer:

```bash
Rscript scripts/02_download_mc3_maf.R
```

To force a fresh MC3 retrieval attempt with `TCGAbiolinks::getMC3MAF()`, use:

```bash
Rscript scripts/02_download_mc3_maf.R --refresh
```

The mutation layer writes:

- `data/raw/mc3/mc3.maf.gz` when MC3 is retrieved through TCGAbiolinks
- `data/processed/mutations/mc3_somatic_mutations.parquet`
- `data/processed/features/tmb_by_sample.tsv`
- `data/processed/features/mutation_prevalence_by_project.tsv`
- `data/processed/features/mutation_prevalence_by_project_driver_only.tsv`
- `results/tables/mutation_layer_qc_summary.tsv`

The current `tmb_by_sample.tsv` output is labeled `mutation_count_proxy_no_callable_territory`. It contains nonsynonymous and total mutation counts per sample, not true mutations per megabase, until a reliable callable territory denominator is added.

Project-level mutation prevalence is appropriate as a molecular feature for Level 1 and Level 2 similarity trees. It is not a direct clonal phylogeny input; Level 3 requires sample-level variant read counts, purity, and local copy-number context.

Build the copy-number, purity, ploidy, and aneuploidy layer:

```bash
Rscript scripts/03_download_copy_number.R
Rscript tests/test_copy_number_layer_synthetic.R
Rscript tests/test_copy_number_layer_outputs.R
```

This layer first checks for two public PanCanAtlas/GDC files at exact local paths, then attempts to download them from the GDC API if they are absent:

- `data/raw/purity_ploidy/TCGA_mastercalls.abs_tables_JSedit.fixed.txt`
- `data/raw/copy_number/PANCAN_ArmCallsAndAneuploidyScore_092817.txt`

If automatic download fails, manually download the files from the GDC PanCanAtlas publication pages and place them at those exact paths. The ABSOLUTE purity/ploidy file is listed on the TCGA PanCanAtlas publication page, and the aneuploidy score/arm-call file is listed on the PanCanAtlas aneuploidy publication page.

After the two primary files are checked, the layer also scans `data/raw/copy_number/` and `data/raw/purity_ploidy/` for additional local `.tsv`, `.tsv.gz`, `.csv`, `.csv.gz`, `.txt`, `.txt.gz`, and `.parquet` files. The script recognizes common column names for sample barcodes, purity, ploidy, aneuploidy score, fraction genome altered, arm-level calls, and segment-level copy-number values. It does not require R `arrow`; TSV outputs are used for feature tables, and Parquet input falls back to Python `pyarrow` when needed.

Expected source columns can include:

- ABSOLUTE purity/ploidy: `sample`, `purity`, `ploidy`, `Genome doublings`, and optionally `Subclonal genome fraction`.
- PanCanAtlas aneuploidy/arm calls: `Sample`, `Type`, `Aneuploidy Score`, followed by wide chromosome-arm columns such as `1p`, `1q`, and `13 (13q)` encoded as `-1`, `0`, or `1`.
- Other purity/ploidy files: `Tumor_Sample_Barcode` or `sample_barcode`, plus `purity` and/or `ploidy`.
- Other direct copy-number traits: `aneuploidy_score`, `fraction_genome_altered` or `fga`, arm gain/loss counts, or copy-number complexity/segment counts.
- Arm calls: sample barcode, `arm` or `chromosome_arm`, and `call` values such as `gain`, `loss`, or `neutral`.
- Segments: sample barcode, chromosome, start, end, and a segment value such as `segment_mean`, `log2`, or `copy_number`.

The layer writes:

- `data/processed/features/purity_ploidy_by_sample.tsv`
- `data/processed/features/aneuploidy_by_sample.tsv`
- `data/processed/features/aneuploidy_by_project.tsv`
- `data/processed/copy_number/arm_level_calls_by_sample.tsv`
- `data/processed/copy_number/segments_by_sample.tsv`
- `results/tables/copy_number_layer_qc_summary.tsv`
- `results/figures/copy_number_layer_qc_overview.pdf`

Feature meaning:

- `purity`: estimated tumor cellularity from the PanCanAtlas ABSOLUTE purity/ploidy table when available.
- `ploidy`: estimated tumor ploidy from the PanCanAtlas ABSOLUTE purity/ploidy table when available.
- `whole_genome_doubling_status`: derived from the ABSOLUTE `Genome doublings` field as `WGD` for one or more genome doublings and `no_WGD` for zero.
- `aneuploidy_score`: directly sourced from `PANCAN_ArmCallsAndAneuploidyScore_092817.txt` when available.
- `fraction_genome_altered`: directly sourced FGA when available, or an approximate segment-length-weighted altered fraction when derived from segment means. The primary PanCanAtlas arm-call file does not provide FGA, so this field may remain missing.
- `arm_gain_count`, `arm_loss_count`, `total_arm_alteration_count`: derived from PanCanAtlas arm-level `-1/0/1` calls or compatible local arm-call files.
- `copy_number_complexity_score`: currently a conservative segment-count proxy when derived from segment-level data.

If the PanCanAtlas files cannot be downloaded and no local copy-number resources are available, the script still writes mutation-sample-aligned tables with missing feature values and explicit `no_source_available` notes. These copy-number and purity/ploidy features are intended for Level 1/2 molecular trait mapping and as covariates for Level 3 clonal analysis; segment-derived summaries are approximate and should not be treated as a substitute for fully harmonized ABSOLUTE or allele-specific copy-number calls.

Build the RNA expression and pathway-feature layer:

```bash
Rscript scripts/04_download_expression.R
Rscript tests/test_expression_layer_synthetic.R
Rscript tests/test_expression_layer_outputs.R
```

This layer first checks for the PanCanAtlas expression matrix at:

- `data/raw/expression/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv`

If the file is absent, the script attempts to download it from the GDC PanCanAtlas publication resource endpoint. If automatic download fails, manually download `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv` from the GDC PanCanAtlas publication page and place it at the exact path above. If no source is available, the script writes mutation-sample-aligned missing-feature tables with `no_source_available` notes so downstream workflow contracts remain explicit.

The expression layer writes:

- `data/processed/expression/expression_matrix_by_sample.tsv.gz`
- `data/processed/features/expression_pathway_scores_by_sample.tsv`
- `data/processed/features/expression_pathway_scores_by_project.tsv`
- `data/processed/features/expression_pca_by_sample.tsv`
- `data/processed/features/expression_pca_variance.tsv`
- `results/tables/expression_layer_qc_summary.tsv`
- `results/figures/expression_layer_qc_overview.pdf`

For the default PanCanAtlas source, expression values are treated as EBPlusPlus-adjusted RNASeqV2 gene expression values and used as provided. Pathway scores are calculated as the average z-scored expression of genes in each configured gene set from `config/expression_gene_sets.yaml`. If a local file name indicates raw counts, the script applies library-size normalization followed by log2 CPM plus one. If a local file name indicates FPKM or TPM, it applies log2(value + 1).

Expression feature meaning:

- `proliferation_score`: proliferation marker score using genes such as `MKI67`, `TOP2A`, and MCM family members.
- `immune_inflammatory_score`: leukocyte, T-cell, cytotoxic, and inflammatory chemokine marker score.
- `interferon_gamma_score`: interferon-gamma response and antigen-presentation-associated marker score.
- `cytotoxic_t_cell_score`: cytotoxic lymphocyte effector score.
- `stromal_score`: fibroblast, extracellular matrix, and activated stromal marker score.
- `epithelial_score`: epithelial differentiation and keratin marker score.
- `EMT_score`: epithelial-mesenchymal transition and invasion-associated marker score.
- `hypoxia_score`: hypoxia and glycolytic adaptation marker score.
- `cell_cycle_score`: cell-cycle and DNA replication marker score.
- `DNA_repair_score`: DNA damage response and homologous recombination repair marker score.
- `angiogenesis_score`: angiogenesis and endothelial activation marker score.
- `lineage_or_tissue_PC1` and `lineage_or_tissue_PC2`: first two PCs from the selected configured expression-gene matrix.

These expression features are intended for Level 1 and Level 2 molecular similarity trees and trait mapping. They are phenotype summaries, not direct clonal phylogeny inputs. The starter gene sets are deliberately simple and transparent; they should be treated as first-pass signatures rather than definitive pathway activity estimates. The default PCA is computed from the configured expression-gene subset to keep the layer memory-conscious, not from the full transcriptome.

Build the Level 1 project-level feature matrix:

```bash
Rscript scripts/06_build_pan_cancer_features.R
Rscript tests/test_level1_feature_matrix_synthetic.R
Rscript tests/test_level1_feature_matrix_outputs.R
```

This integration step joins project-level mutation-count proxies, driver mutation prevalence, copy-number/purity/ploidy summaries, aneuploidy/arm-level burden, expression pathway scores, and manual cancer group annotations. It writes:

- `results/tables/level1_feature_matrix.tsv`
- `results/tables/level1_feature_support.tsv`
- `results/tables/level1_feature_matrix_qc_summary.tsv`

The feature matrix contains project annotations plus retained biologic features. Sample-support counts are written to the support table rather than included as clustering features. Features with missingness greater than `missingness.max_feature_missing_fraction` in `config/feature_sets.yaml` are dropped from the matrix and recorded in the QC summary. Current feature values are raw and unscaled; scaling, imputation decisions, and distance-specific preprocessing should be handled explicitly in the Level 1 tree script.

Build the Level 1 pan-cancer molecular similarity dendrograms:

```bash
Rscript scripts/07_level1_pan_cancer_tree.R
Rscript tests/test_level1_tree_synthetic.R
Rscript tests/test_level1_tree_outputs.R
```

This step reads `results/tables/level1_feature_matrix.tsv`, excludes annotation columns from distance calculations, median-imputes any remaining numeric missing values, scales retained biologic features, and writes:

- `results/tables/level1_feature_matrix_scaled.tsv`
- `results/tables/level1_gower_distance_matrix.tsv`
- `results/tables/level1_correlation_distance_matrix.tsv`
- `results/tables/level1_bootstrap_cluster_stability.tsv`
- `results/tables/level1_tree_qc_summary.tsv`
- `results/trees/level1_pan_cancer_gower_hclust_tree.nwk`
- `results/trees/level1_pan_cancer_correlation_hclust_tree.nwk`
- `results/trees/level1_pan_cancer_neighbor_joining_tree.nwk`
- `results/figures/level1_pan_cancer_gower_tree.pdf`
- `results/figures/level1_pan_cancer_correlation_tree.pdf`
- `results/figures/level1_pan_cancer_tree_with_trait_heatmap.pdf`
- `results/figures/level1_pan_cancer_feature_pca.pdf`

Interpretation:

- The Gower tree is the primary average-linkage pan-cancer molecular similarity dendrogram on scaled project-level molecular traits. Because all tree features are numeric after preprocessing, it summarizes scaled absolute trait differences across mutation, copy-number, purity/ploidy, and expression features.
- The correlation tree is a sensitivity view based on similarity of each TCGA project's feature profile shape across traits. It can group projects with similar relative molecular profiles even when absolute scaled levels differ.
- The neighbor-joining Newick tree is written as a sensitivity representation when the Gower distance matrix is valid. It should not be read as a literal species-like cancer evolutionary phylogeny.
- The bootstrap table is a first-pass feature-bootstrap co-clustering summary: features are sampled with replacement, average-linkage Gower clustering is rebuilt, and project-pair co-clustering frequencies are reported at a practical fixed `cutree` resolution.

Scientific caveat: these Level 1 outputs are a pan-cancer molecular similarity dendrogram, molecular trait tree, or lineage-informed molecular similarity tree. They do not imply that one TCGA cancer type evolved from another. Serial or multi-region samples are not required for this project-level trait tree; those samples matter later for Level 3 within-patient clonal phylogeny.

Interpret the Level 1 tree and run sensitivity analyses:

```bash
Rscript scripts/07b_interpret_level1_tree.R
Rscript tests/test_level1_interpretation_synthetic.R
Rscript tests/test_level1_interpretation_outputs.R
```

This step reads the Level 1 feature matrix, scaled matrix, Gower and correlation distances, Newick trees, and bootstrap co-clustering table. It writes:

- `results/tables/level1_gower_cluster_assignments.tsv`
- `results/tables/level1_correlation_cluster_assignments.tsv`
- `results/tables/level1_gower_cluster_summaries.tsv`
- `results/tables/level1_correlation_cluster_summaries.tsv`
- `results/tables/level1_feature_pca_loadings.tsv`
- `results/tables/level1_feature_variance_summary.tsv`
- `results/tables/level1_feature_cluster_association.tsv`
- `results/tables/level1_feature_class_sensitivity_summary.tsv`
- `results/tables/level1_feature_class_cluster_assignments.tsv`
- `results/tables/level1_feature_weighting_sensitivity_summary.tsv`
- `results/tables/level1_weighted_cluster_assignments.tsv`
- `results/tables/level1_stable_project_pairs.tsv`
- `results/tables/level1_bootstrap_project_stability.tsv`
- `results/reports/level1_interpretation_summary.md`
- `results/figures/level1_cluster_trait_summary_heatmap.pdf`
- `results/figures/level1_feature_driver_barplot.pdf`
- `results/figures/level1_bootstrap_coclustering_heatmap.pdf`
- `results/figures/level1_feature_class_sensitivity_heatmap.pdf`
- `results/figures/level1_weighting_sensitivity_pca_or_mds.pdf`

The interpretation step is descriptive and sensitivity-focused. It should be used to identify which project groupings are robust across distance, feature subset, weighting, and bootstrap views, not to make directional evolutionary claims between cancer types.

Build Level 2 within-group molecular similarity trees:

```bash
Rscript scripts/08_level2_within_group_trees.R
Rscript tests/test_level2_tree_synthetic.R
Rscript tests/test_level2_tree_outputs.R
```

This step reads `results/tables/level1_feature_matrix.tsv`, `results/tables/level1_feature_matrix_scaled.tsv`, `config/cancer_group_map.csv`, and `data/interim/tcga_projects.tsv`. It defines carcinoma-all and selected biologically motivated subgroups, drops features that lack within-group variability, scales features within each group, and writes:

- `results/tables/level2_group_definitions.tsv`
- `results/tables/level2_tree_qc_summary.tsv`
- `results/reports/level2_interpretation_summary.md`
- group-specific feature matrices and scaled matrices under `results/tables/level2_<group_name>_*`
- group-specific Gower/correlation distance matrices where tree construction is attempted
- group-specific Gower/correlation Newick trees under `results/trees/`
- group-specific tree, trait-heatmap, and PCA figures under `results/figures/`
- bootstrap co-clustering tables for groups with at least five projects

Groups with fewer than three projects are not forced into trees. They are retained in the group-definition and QC outputs with explicit underpowered reasons.

Select Level 3 clonal-analysis candidates:

```bash
python scripts/09_select_level3_clonal_candidates.py
python tests/test_level3_candidate_selection_synthetic.py
python tests/test_level3_candidate_selection_outputs.py
```

This step is feasibility and pilot-cohort selection, not full clonal phylogeny. It reads:

- `data/processed/mutations/mc3_somatic_mutations.parquet`
- `data/processed/features/tmb_by_sample.tsv`
- `data/processed/features/purity_ploidy_by_sample.tsv`
- `data/processed/features/aneuploidy_by_sample.tsv`
- `data/interim/tcga_projects.tsv`
- `config/cancer_group_map.csv`
- optional copy-number segment and arm-call files under `data/processed/copy_number/`

It writes:

- `results/tables/level3_candidate_samples.tsv`
- `results/tables/level3_candidate_summary_by_project.tsv`
- `results/tables/level3_pilot_cohort.tsv`
- `results/tables/level3_excluded_samples.tsv`
- `results/tables/level3_candidate_selection_qc_summary.tsv`
- `results/figures/level3_candidate_mutation_count_vs_purity.pdf`
- `results/figures/level3_candidate_counts_by_project.pdf`
- `results/figures/level3_candidate_score_by_project.pdf`
- `results/figures/level3_pilot_cohort_overview.pdf`

How Level 3 differs from Levels 1 and 2:

- Levels 1 and 2 are project-level molecular similarity trees built from aggregated molecular traits. They do not require serial samples and are not within-patient clonal phylogenies.
- Level 3 works within individual tumors. The unit of analysis becomes mutations/subclones within a sample or patient, so ref/alt counts, purity, ploidy, and local copy-number context matter directly.
- Serial, relapse, metastasis, or multi-region samples are preferred for robust branching order. Most TCGA cases are single bulk samples, so clonal clustering is usually more defensible than exact branching phylogeny.
- PyClone-VI-style inputs generally require variant read counts, tumor purity, and copy-number context for each mutation. PhyloWGS-style workflows require mutation read counts plus copy-number information and careful sample-level QC.

Candidate scoring is transparent and intentionally simple. The current 0-100 score combines up to 40 points for log-scaled nonsynonymous mutation-count proxy, 25 points for purity, 25 points for availability of ref/alt counts plus purity/ploidy/copy-number summaries, and 10 points for segment-level copy-number availability. Penalties are applied for missing segment-level copy number, low mutation count, low purity, missing ref/alt counts, missing purity/ploidy, and missing copy-number data. Thresholds are defined at the top of `scripts/09_select_level3_clonal_candidates.py`; current defaults require ref/alt counts, purity, ploidy, nonsynonymous mutation count at least 100 for solid tumors or 50 for hematolymphoid tumors, and purity at least 0.30.

Candidate classes:

- `copy_number_aware_clonal_phylogeny_candidate`: passes filters and has segment-level copy number.
- `clonal_clustering_candidate_limited_cn`: passes filters but lacks segment-level copy number.
- `descriptive_only`: has partial data but is insufficient for reliable clonal inference.
- `excluded`: fails required minimum filters.

Acquire or join segment-level copy-number data:

```bash
python scripts/09c_fetch_gdc_segment_cn.py
python tests/test_gdc_segment_cn_fetch_synthetic.py
python tests/test_gdc_segment_cn_outputs.py
```

This GDC-backed step queries the GDC files API for open-access TCGA Copy Number Variation files, restricted to the current Level 3 pilot projects. It searches for `Masked Copy Number Segment`, `Copy Number Segment`, and `Filtered Copy Number Segment` files, matches them to pilot samples using TCGA sample barcodes when possible, records patient-level matches only as lower-confidence evidence, selects one best sample-level file per pilot sample, downloads to `data/raw/copy_number_segments/gdc_pilot/`, parses segment tables, and updates the Level 3 candidate and pilot cohort tables. The default top-of-script settings are `DOWNLOAD_MODE = "pilot_full"`, `ALLOW_PATIENT_LEVEL_MATCHES = False`, and `OVERWRITE_EXISTING = False`; capped modes use `MAX_DOWNLOADS` or `--max-download-files`.

The older local-file scanner remains useful after manually placing segment files:

```bash
python scripts/09b_download_segment_level_cn.py
python tests/test_level3_segment_cn_synthetic.py
python tests/test_level3_segment_cn_outputs.py
```

Both segment-CN paths use or write local segment-level CN files under:

- `data/raw/copy_number/`
- `data/raw/copy_number_segments/`
- `data/processed/copy_number/`

Supported local file patterns include `*.seg`, `*.seg.txt`, `*.tsv`, `*.txt`, `*.maf.cnv`, `copy_number_segments*.tsv`, `segments_by_sample*.tsv`, and `segments_by_sample*.parquet`. The parser recognizes common segment columns such as sample ID, chromosome, start, end, number of probes, and segment mean/log2 copy ratio. If no parseable segment rows are found, the script writes explicit no-source QC and manual GDC instructions rather than failing.

Segment-level CN matters for clonal inference because local copy number changes the expected variant allele fraction for a mutation. Arm-level aneuploidy says whether broad chromosome arms are gained or lost, but it does not tell whether a specific mutation sits in a locally amplified, deleted, or neutral segment. Copy-number-aware tools such as PyClone-VI or PhyloWGS-style workflows need local copy-number context to estimate cancer-cell fraction more defensibly from observed VAF. Candidates without sample-level segment CN can still support limited clonal clustering prototypes, but they should not be presented as full copy-number-aware clonal phylogenies. Patient-level segment matches are recorded as lower-confidence matches and should be manually validated before clonal inference; sample-level matches are required for automatic upgrade to `copy_number_aware_clonal_phylogeny_candidate`.

The segment-CN step writes:

- `results/tables/gdc_segment_cn_manifest_pilot.tsv`
- `results/tables/gdc_segment_cn_best_file_per_pilot_sample.tsv`
- `results/tables/gdc_segment_cn_download_log.tsv`
- `data/processed/copy_number/segments_by_sample.tsv.gz`
- `data/processed/copy_number/segment_level_cn_summary_by_sample.tsv`
- `results/tables/level3_candidate_samples_with_segments.tsv`
- `results/tables/level3_pilot_cohort_with_segments.tsv`
- `results/tables/level3_segment_cn_qc_summary.tsv`
- `results/tables/level3_segment_cn_availability_by_project.tsv`
- `results/figures/level3_segment_cn_availability_by_project.pdf`
- `results/figures/level3_segment_cn_pilot_coverage.pdf`
- `results/figures/level3_segment_cn_upgrade_by_project.pdf`
- `results/figures/level3_segment_count_distribution.pdf`
- `results/figures/level3_candidate_upgrade_status.pdf`

Annotate pilot mutations with local segment-level copy-number context:

```bash
python scripts/10_annotate_mutations_with_local_cn.py
python tests/test_mutation_local_cn_annotation_synthetic.py
python tests/test_mutation_local_cn_annotation_outputs.py
```

This step overlaps each somatic mutation from the copy-number-aware Level 3 pilot cohort with the sample-level segment interval covering the mutation coordinate. It harmonizes TCGA sample barcodes and chromosome names, computes total depth and observed VAF, records local segment coordinates and segment mean, labels conservative gain-like/loss-like/neutral-like status, and writes a non-final PyClone-style preview. Patient-level segment matching is off by default. The preview does not infer allele-specific major/minor copy number: `major_cn` and `minor_cn` remain `NA`, and autosomal `normal_cn = 2` is labeled as a placeholder.

The mutation local-CN annotation step writes:

- `data/processed/clonal/level3_mutations_with_local_cn.tsv.gz`
- `results/tables/level3_mutation_local_cn_summary_by_sample.tsv`
- `data/processed/clonal/level3_pyclone_input_preview.tsv`
- `results/tables/level3_mutation_local_cn_annotation_qc_summary.tsv`
- `results/figures/level3_mutation_cn_annotation_coverage_by_sample.pdf`
- `results/figures/level3_mutation_vaf_vs_local_segment_mean.pdf`
- `results/figures/level3_mutation_depth_distribution.pdf`

Prepare stratified clonal input sets and readiness tables:

```bash
python scripts/11_prepare_clonal_input_sets.py
python tests/test_clonal_input_preparation_synthetic.py
python tests/test_clonal_input_preparation_outputs.py
```

This step does not run a clonal clustering tool. It stratifies mutation/local-CN annotations into a copy-neutral near-diploid set suitable for limited VAF-based clonal clustering prototypes, a broader segment-annotated set for visualization and future CN-aware preparation, and a sample-level readiness table. It explicitly labels current outputs as `not_copy_number_aware` or `requires_allele_specific_cn` when major/minor CN is unavailable.

The clonal input preparation step writes:

- `data/processed/clonal/level3_copy_neutral_vaf_clustering_input.tsv.gz`
- `data/processed/clonal/level3_segment_annotated_clonal_input.tsv.gz`
- `results/tables/level3_clonal_input_readiness_by_sample.tsv`
- `results/tables/level3_limited_vaf_clustering_pilot_samples.tsv`
- `results/tables/level3_clonal_input_preparation_qc_summary.tsv`
- `results/figures/level3_copy_neutral_candidate_counts_by_sample.pdf`
- `results/figures/level3_copy_neutral_candidate_counts_by_project.pdf`
- `results/figures/level3_vaf_distribution_copy_neutral_candidates.pdf`
- `results/figures/level3_readiness_class_by_project.pdf`

Run the limited VAF-based clonal clustering prototype:

```bash
python scripts/12_limited_vaf_clonal_clustering.py
python tests/test_limited_vaf_clonal_clustering_synthetic.py
python tests/test_limited_vaf_clonal_clustering_outputs.py
```

This prototype clusters observed mutation VAFs independently within the selected 30 copy-neutral pilot samples. It uses one-dimensional Gaussian mixture models selected by BIC, keeps the simplest model within a conservative BIC tolerance, and falls back to quantile bins if mixture fitting is not available or fails quality checks. Cluster interpretations are cautious descriptors of VAF structure only. This step does not run PyClone-VI or PhyloWGS, is not copy-number-aware beyond the upstream near-neutral segment filter, does not infer allele-specific integer copy number, and does not infer branching order.

The limited VAF clustering step writes:

- `results/clonal/limited_vaf/level3_limited_vaf_cluster_assignments.tsv.gz`
- `results/tables/level3_limited_vaf_clonal_complexity_by_sample.tsv`
- `results/tables/level3_limited_vaf_clonal_complexity_by_project.tsv`
- `results/tables/level3_limited_vaf_clustering_qc_summary.tsv`
- `results/figures/level3_limited_vaf_cluster_vaf_density_by_sample.pdf`
- `results/figures/level3_limited_vaf_cluster_counts_by_project.pdf`
- `results/figures/level3_limited_vaf_complexity_heatmap.pdf`
- `results/figures/level3_limited_vaf_example_samples.pdf`

### Path to PyClone-VI and PhyloWGS

```bash
python scripts/12b_allele_specific_cn_readiness.py
python tests/test_allele_specific_cn_readiness_synthetic.py
python tests/test_allele_specific_cn_readiness_outputs.py
```

The current GDC segment files contain continuous segment means or log2-like signals. `segment_mean` is useful for gain-like, loss-like, and near-neutral filtering, but it is not allele-specific integer copy number. It does not directly encode `total_cn`, `major_cn`, `minor_cn`, mutation multiplicity, or a validated LOH state. The existing 30-sample limited-VAF clustering therefore remains a pre-PyClone exploratory prototype restricted to copy-neutral/near-diploid regions.

The readiness step searches these local locations recursively:

- `data/raw/allele_specific_cn/`
- `data/processed/allele_specific_cn/`
- `data/raw/facets/`
- `data/raw/sequenza/`
- `data/raw/ascat/`
- `data/raw/absolute/`

It supports flexible FACETS-, Sequenza-, ASCAT-, ABSOLUTE-style, and generic allele-specific segment columns. Deterministic total/major/minor relationships may be standardized when true integer allele states are present, such as FACETS total and lesser-copy-number fields. Major or minor CN is never inferred from `segment_mean`. Candidate rows with noninteger, negative, or inconsistent total/major/minor states do not contribute to readiness.

For the current pilot, no supported local allele-specific CN files, paired tumor/normal BAMs, SNP pileups, seqz files, or comparable allele-count inputs were found. The readiness outputs therefore record:

- pilot samples evaluated: 84;
- limited-VAF pilot samples evaluated: 30;
- samples with existing allele-specific CN: 0;
- samples with major and minor integer CN: 0;
- samples ready for PyClone-VI input preparation: 0;
- samples ready for cautious PhyloWGS input preparation: 0; and
- samples requiring external allele-specific CN inference: 84.

PyClone-VI input preparation can begin only after a sample has reviewed integer major/minor CN, purity, ref/alt counts, at least 50 eligible mutation overlaps, and at least 80% allele-specific CN coverage under the current conservative defaults. The script does not run PyClone-VI. After PyClone-VI inputs and clustering pass per-sample QC, PhyloWGS can be considered only for a smaller high-confidence subset with at least 100 eligible mutations, at least 90% CN overlap, interpretable cluster complexity, and explicit single-bulk-sample limitations. The script does not run PhyloWGS.

Recommended acquisition path:

1. Prioritize the 30 limited-VAF pilot samples.
2. Obtain paired tumor/normal BAMs or validated precomputed common-SNP allele counts.
3. Run one allele-specific CN method consistently, with FACETS as the default planning option for paired WES when appropriate.
4. Review genome build, sample/aliquot identity, segmentation, integer-CN semantics, purity/ploidy, and QC.
5. Standardize accepted outputs into `data/processed/allele_specific_cn/allele_specific_segments_by_sample.tsv.gz`.
6. Re-annotate mutations with total, major, and minor CN before preparing PyClone-VI input.
7. Consider PhyloWGS only after PyClone-VI and allele-specific CN review identify a high-confidence subset.

The readiness step writes:

- `results/tables/level3_allele_specific_cn_availability_by_sample.tsv`
- `data/processed/allele_specific_cn/allele_specific_segments_by_sample.tsv.gz`
- `results/tables/level3_mutation_to_allele_specific_cn_readiness.tsv`
- `results/tables/level3_allele_specific_cn_external_inference_manifest.tsv`
- `results/tables/level3_allele_specific_cn_readiness_qc_summary.tsv`
- `results/reports/level3_allele_specific_cn_and_pyclone_plan.md`
- `results/figures/level3_allele_specific_cn_readiness_by_project.pdf`
- `results/figures/level3_allele_specific_cn_missing_inputs.pdf`

Do not run PyClone-VI or PhyloWGS until allele-specific CN and mutation-overlap QC are available. Copy-neutral exploratory analysis may continue only under the existing limited-VAF label and must not be described as copy-number-aware clonal phylogeny.

Build the TCGA sample-level molecular reference atlas:

```bash
python scripts/20_build_tcga_sample_reference_atlas.py
python tests/test_tcga_sample_reference_atlas_synthetic.py
python tests/test_tcga_sample_reference_atlas_outputs.py
```

This translational atlas is different from the Level 1 and Level 2 trees. Levels 1 and 2 aggregate molecular features to the TCGA project level and build project-level molecular similarity trees. The reference atlas preserves sample-level variation across individual TCGA tumors so it can later support nearest-neighbor search, PCA projection, UMAP-style projection, classifier training, and cancer-of-unknown-primary molecular comparison for clinical WES/WTS cases.

The atlas uses sample-level mutation-count proxies, configured driver-gene mutation indicators and mutation counts, purity, ploidy, aneuploidy and arm-level alteration counts, expression pathway scores, and expression PCs when available. The raw matrix preserves observed pre-imputation values for all candidate features. The scaled matrix contains only retained features after TCGA-reference median imputation and standard scaling. Future clinical samples must be transformed with the saved artifacts:

- Use `results/reference_atlas/tcga_reference_feature_definitions.yaml` to compute exactly the retained WES/WTS-compatible features.
- Use `results/reference_atlas/tcga_reference_scaler_parameters.json` to apply TCGA reference medians, means, and standard deviations.
- Use `results/reference_atlas/tcga_reference_nearest_neighbor_metadata.json` to select the recorded Euclidean default or cosine alternative and map neighbor rows back to sample metadata.
- Then project or compare against `results/reference_atlas/tcga_reference_pca_model.pkl` and `results/reference_atlas/tcga_reference_nearest_neighbor_index.pkl`.

The atlas step writes:

- `results/reference_atlas/tcga_sample_reference_matrix.tsv.gz`
- `results/reference_atlas/tcga_sample_reference_metadata.tsv`
- `results/reference_atlas/tcga_sample_reference_scaled_matrix.tsv.gz`
- `results/reference_atlas/tcga_reference_feature_definitions.yaml`
- `results/reference_atlas/tcga_reference_scaler_parameters.json`
- `results/reference_atlas/tcga_reference_feature_missingness.tsv`
- `results/reference_atlas/tcga_reference_pca_coordinates.tsv`
- `results/reference_atlas/tcga_reference_pca_model.pkl`
- `results/reference_atlas/tcga_reference_pca_variance.tsv`
- `results/reference_atlas/tcga_reference_nearest_neighbor_index.pkl`
- `results/reference_atlas/tcga_reference_nearest_neighbor_metadata.json`
- `results/tables/tcga_sample_reference_atlas_qc_summary.tsv`
- `results/figures/tcga_sample_reference_pca_by_project.pdf`
- `results/figures/tcga_sample_reference_pca_by_major_group.pdf`
- `results/figures/tcga_sample_reference_feature_missingness.pdf`
- `results/figures/tcga_sample_reference_samples_by_project.pdf`

The locked atlas is consumed by the clinical projection prototype below. The atlas builder itself remains separate from clinical-query ingestion.

### Clinical WES/WTS query projection prototype

```bash
python scripts/21_project_clinical_case_to_tcga.py --case-id example_cup_case
python tests/test_clinical_projection_synthetic.py
python tests/test_clinical_projection_outputs.py
```

This prototype ingests already-computed clinical molecular features from `data/clinical_queries/<case_id>/`. It does not process raw FASTQ, BAM, VCF, or RNA-seq files. A query directory contains required `clinical_metadata.tsv` plus `wes_features.tsv`, `wts_features.tsv`, or both. Reference-compatible feature columns must use exact names from `results/reference_atlas/tcga_reference_feature_definitions.yaml`, and projection uses the locked 97-feature TCGA reference contract.

The projection step:

- records exact supplied, missing, median-imputed, and extra ignored features;
- applies the saved TCGA reference medians, means, and standard deviations without fitting on the clinical case;
- queries the compatible saved Euclidean or cosine nearest-neighbor model, with direct scaled-matrix distances as a fallback;
- projects the query with the saved TCGA PCA model without refitting PCA;
- summarizes top neighbors by TCGA project and major cancer group; and
- writes a research-use Markdown report with explicit clinical caveats.

The artificial example is stored under `data/clinical_queries/example_cup_case/`. It contains no real patient data and intentionally omits expression PCs to exercise the missing-feature audit and TCGA-median imputation path.

Outputs are written under `results/clinical_projection/<case_id>/`:

- `clinical_feature_vector.tsv`
- `clinical_feature_vector_scaled.tsv`
- `clinical_feature_missingness.tsv`
- `tcga_nearest_neighbors.tsv`
- `tcga_nearest_project_summary.tsv`
- `tcga_nearest_major_group_summary.tsv`
- `clinical_tcga_pca_projection.tsv`
- `clinical_projection_qc_summary.tsv`
- `clinical_tcga_pca_projection_by_project.pdf`
- `clinical_tcga_pca_projection_by_major_group.pdf`
- `clinical_nearest_neighbor_barplot.pdf`
- `clinical_feature_missingness.pdf`
- `clinical_tcga_projection_report.md`

This module is a research/prototype molecular projection framework. It is not clinically validated, does not establish tissue of origin, and does not train or apply a supervised CUP classifier. The separate PyClone-VI/PhyloWGS branch remains paused until validated allele-specific integer copy number is available.

### CUP molecular interpretation prototype

```bash
python scripts/22_cup_nearest_neighbor_classifier.py --case-id example_cup_case
python tests/test_cup_classifier_synthetic.py
python tests/test_cup_classifier_outputs.py
```

This downstream workflow interprets an existing TCGA nearest-neighbor projection against the current locked 97-feature TCGA atlas. It does not process raw clinical data, refit the atlas, or train a supervised CUP classifier. The script:

- validates the top-k neighbor rows against the project and major cancer-group summaries;
- maps all 33 TCGA projects to 18 explicit lineage groups from `config/cup_lineage_groups.yaml`;
- calculates top-lineage weight, top-1/top-2 margin, normalized entropy, meaningful-project count, overall/WES/WTS feature coverage, broad-group agreement, and neighbor concentration;
- assigns transparent `high_support`, `moderate_support`, `weak_support`, or `indeterminate` tiers that are not clinically calibrated;
- evaluates conservative rules from `config/cup_feature_evidence_rules.yaml` only when their locked features were supplied rather than median-imputed;
- summarizes nearest-neighbor, major-group, expression/pathway, driver, copy-number, immune/stromal/EMT, contradictory, and missing evidence;
- highlights candidates named in the submitted differential diagnosis without changing the unconstrained molecular ranking; and
- writes a transparent heuristic molecular-pathologic concordance table and provides an actionability placeholder without performing clinical actionability interpretation.

Outputs are written under `results/clinical_projection/<case_id>/cup_interpretation/`:

- `cup_ranked_lineage_interpretation.tsv`
- `cup_ambiguity_metrics.tsv`
- `cup_molecular_evidence_table.tsv`
- `cup_differential_diagnosis_constrained_summary.tsv`
- `cup_molecular_pathologic_concordance.tsv`
- `cup_molecular_interpretation_report.md`
- `cup_top_project_similarity_barplot.pdf`
- `cup_major_group_similarity_barplot.pdf`
- `cup_ambiguity_summary.pdf`
- `cup_evidence_heatmap.pdf`

For the artificial example, the unconstrained top lineage is `lower_GI` with 0.5701 of normalized top-k similarity weight, followed by `pancreatobiliary_hepatobiliary` with 0.3092. The normalized entropy is 0.7372, the top-1/top-2 margin is 0.2609, and the overall class is `moderate_ambiguity`. Within the submitted synthetic gastrointestinal/pancreatobiliary/lung differential, the broad gastrointestinal candidate has the greatest relative support. These values are workflow demonstrations, not tissue-of-origin probabilities or a clinical diagnosis.

The report uses this required interpretation boundary: the analysis compares a query tumor with TCGA molecular reference profiles, is not a clinically validated tissue-of-origin assay, and must be integrated with morphology, immunophenotype, imaging, and clinical findings. Full OncoKB/CIViC actionability is intentionally not implemented. The separate PyClone-VI/PhyloWGS branch remains paused until validated allele-specific integer copy number is available.

### Known-primary internal validation

```bash
python scripts/23_validate_known_primary_projection.py --top-k 50 --metric euclidean
python tests/test_known_primary_validation_synthetic.py
python tests/test_known_primary_validation_outputs.py
```

This workflow treats each locked-atlas TCGA sample as a pseudo-unknown query and searches the remaining atlas samples using the existing 97-feature scaled space. The query sample itself is always excluded; `--exclude-same-patient` can additionally exclude other aliquots from the same patient. `--max-samples 1000` provides a deterministic project-stratified quick run, while the default `--max-samples all` validates the complete cohort.

The full default run evaluated 10,201 samples from all 33 projects and 18 configured CUP lineage groups. Known-project recovery was 59.455% at top 1, 82.364% at top 3, and 89.550% at top 5. Top-1 and top-3 lineage recovery were 66.503% and 88.668%, respectively, and top-1 major-group recovery was 88.246%. The highest project-level top-1 recovery was observed for THCA (93.800%), UVM (93.750%), PCPG (90.761%), and KIRC (88.649%); the strongest configured lineages were thyroid (93.600%), CNS/glial (85.899%), breast/gynecologic/Mullerian (84.308%), prostate (81.325%), and kidney (80.362%). These are descriptive internal estimates, not external clinical performance.

The current ambiguity classes separated easier from harder internal cases: low-, moderate-, and high-ambiguity strata contained 4,370, 2,736, and 3,095 samples, with top-1 project recovery of 82.998%, 49.781%, and 34.766%, respectively. Prominent exact-project confusions included READ to COAD, BRCA to OV, BLCA to BRCA, STAD to BRCA, LUAD to OV, UCEC to OV, HNSC to LUSC/CESC, and LGG to GBM. These patterns identify failure modes and support future abstention calibration; current similarity weights and ambiguity tiers are not diagnostic probabilities.

Outputs are written under `results/validation/known_primary_projection/` and include the validation cohort, per-sample predictions, project/lineage/major-group performance tables, overall metrics, three confusion matrices, ambiguity summaries, seven PDF figures, and `known_primary_validation_report.md`.

This evaluation reuses the same TCGA cohort, assays, feature engineering, global imputation values, and scaling parameters for query and reference samples. Self exclusion prevents identity matching but does not remove cohort or preprocessing leakage. Independent known-primary WES/WTS cohorts, assay harmonization, patient-level exclusion sensitivity, missing-modality analysis, and external uncertainty calibration remain required before clinical use or any tissue-of-origin performance claim.

### Calibration and feature-ablation validation

```bash
python scripts/24_calibrate_projection_confidence_and_feature_ablation.py --top-k 50 --metric euclidean
python tests/test_calibration_feature_ablation_synthetic.py
python tests/test_calibration_feature_ablation_outputs.py
```

This workflow converts the known-primary similarity metrics into transparent, validation-informed confidence tiers and benchmarks 13 prespecified atlas feature sets. Thresholds are derived on a deterministic within-project half split of 5,092 samples and evaluated without refitting on the remaining 5,109 internal samples. Margin, top-1 similarity-weight share, neighbor concentration, and entropy each contribute an equal 0/1/2 vote separately for project, configured lineage, and broad major-group calls. A combined score of at least 6 of 8 defines `high_confidence`, at most 2 defines `low_confidence`, and intermediate scores define `moderate_confidence`; only the high tier receives the research `reliable` flag.

Held-out project accuracy was 35.762%, 55.763%, and 84.837% in low-, moderate-, and high-confidence tiers. Corresponding lineage accuracy was 42.322%, 61.189%, and 92.719%; broad major-group accuracy was 65.482%, 94.875%, and 99.759%. The exact metric-specific gates are versioned in `calibrated_confidence_thresholds.tsv`; entropy is lower-is-better, so its numeric high-confidence threshold is lower than its low-confidence threshold. These labels describe empirical internal reliability and must not be reported as probabilities or clinical certainty.

Feature ablation recomputes exact self-excluded 50-neighbor Euclidean projections while preserving the locked TCGA scaling. The 31-feature WTS-like set had the highest top-1 project recovery at 69.111%, compared with 38.535% for the 66-feature WES-like set and 59.455% for the full 97-feature WES-plus-WTS set. Expression dependence was especially pronounced for LAML, SARC, CESC, UVM, PCPG, THYM, KIRC, and PRAD and for the hematolymphoid, sarcoma/mesenchymal, endocrine/adrenal/paraganglioma, prostate, melanoma, and kidney lineages. The combined unweighted space performing below expression alone is a property of this equal-weight TCGA distance model, not evidence that WES features are uninformative.

Outputs are written under `results/validation/calibration_feature_ablation/` and include calibrated predictions and thresholds, tier-performance summaries, feature definitions, overall/project/lineage ablation tables, long-form confusion matrices, seven PDF figures, and `calibration_feature_ablation_report.md`. This remains internal TCGA validation: the same cohort, assays, preprocessing, median imputation, and scaling are shared across partitions. Independent clinical WES/WTS validation, assay harmonization, realistic partial-WTS simulations, and prospective abstention review are still required before these gates can support clinical reporting.

Generate the integrated project report:

```bash
python scripts/13_generate_integrated_project_report.py
python tests/test_integrated_project_report_outputs.py
```

This reporting layer does not perform new primary analysis. It summarizes existing Level 1, Level 2, Level 3, allele-specific-CN readiness, 97-feature TCGA reference-atlas, synthetic clinical projection, and CUP interpretation outputs into a publication-style integrated report, collaborator-facing executive summary, methods summary, limitations/next-steps document, machine-readable project summary, and expanded output index. The report explicitly preserves the core interpretation boundaries: Levels 1 and 2 are molecular similarity/trait trees, Level 3 is limited VAF clustering rather than a full copy-number-aware clonal phylogeny, and clinical/CUP scores are descriptive similarity weights rather than diagnostic probabilities.

The integrated reporting step writes:

- `results/reports/integrated_tcga_cancer_phylogeny_report.md`
- `results/reports/integrated_tcga_cancer_phylogeny_executive_summary.md`
- `results/reports/integrated_tcga_cancer_phylogeny_methods.md`
- `results/reports/integrated_tcga_cancer_phylogeny_limitations_next_steps.md`
- `results/reports/integrated_tcga_cancer_phylogeny_report.html`
- `results/tables/integrated_project_summary.tsv`
- `results/tables/integrated_project_outputs_index.tsv`

Generate the preprint-style manuscript draft:

```bash
python scripts/14_generate_preprint_manuscript.py
python tests/test_preprint_outputs.py
```

The preprint manuscript step writes:

- `results/reports/tcga_cancer_phylogeny_preprint.md`
- `results/reports/tcga_cancer_phylogeny_preprint_figures_tables_index.tsv`

The refreshed manuscript includes the locked 97-feature atlas, the artificial WES/WTS projection example, the uncertainty-aware CUP interpretation, and the current zero-sample PyClone-VI/PhyloWGS readiness status. It contains eight Markdown tables and references existing Level 1, Level 2, Level 3, atlas, clinical projection, and CUP figures. The conceptual framework figure and formal references remain pending manuscript tasks.

## Current limitations

- Levels 1 and 2 are molecular similarity and trait-mapping analyses, not literal evolutionary histories of cancer types.
- TCGA cohorts differ in tumor purity, assay availability, sample size, and clinical annotation depth.
- Current mutation-count/TMB summaries are mutation-count proxies unless callable territory is added.
- Copy-number and expression features depend on availability and harmonization of external PanCanAtlas/GDC-style resources. Without source files, the layers produce explicit missing-value tables rather than inferred values.
- Expression pathway scores are simple average z-score signatures and do not replace model-based pathway activity, immune deconvolution, or batch-aware transcriptomic analyses.
- Many Level 3 TCGA analyses will be single-bulk-sample analyses, which support clonal clustering more reliably than confident branching order.
- Limited VAF clustering is a prototype summary of VAF cluster structure in copy-neutral regions only. It is not a definitive clonal phylogeny, does not infer branching order, and does not replace copy-number-aware methods that require allele-specific integer copy number. Single bulk TCGA samples further limit phylogenetic interpretation.
- The allele-specific CN readiness step currently finds zero local major/minor CN resources and zero PyClone-VI/PhyloWGS-ready samples. Its command templates are planning placeholders that require fixed tool versions, genome-build resources, tumor/normal identity review, and per-sample QC before external execution.
- MSI, HRD, aneuploidy, purity, ploidy, and whole-genome doubling features depend on availability and harmonization of external PanCanAtlas or related resources.
- PyClone-VI and PhyloWGS integration may require additional installation and per-sample quality control beyond the base workflow.
- The TCGA sample-level reference atlas and clinical query projection module are not diagnostic classifiers. Internal TCGA pseudo-unknown validation and a within-TCGA held-out confidence evaluation are complete, but query and reference samples share cohort selection, assays, feature engineering, global imputation, and scaling. Independent known-primary WES/WTS validation, external threshold evaluation, assay harmonization, and clinical reporting safeguards remain required.
- The CUP interpretation workflow adds transparent lineage grouping, ambiguity metrics, differential matching, and heuristic evidence rules. Validation-informed threshold artifacts now exist, but they have not been independently validated as clinical cutoffs. TCGA primary-tumor bias, incomplete lineage-marker coverage, imputed expression PCs, keyword-based differential matching, shared molecular features across lineages, and absent independent known-primary clinical validation limit interpretation. Actionability integration remains a placeholder.

## Next executable steps

1. Run the existing validation scripts for the completed layers if the environment changes.
2. Install or activate Snakemake, then run `snakemake -n level1_pan_cancer_tree`, `snakemake -n level1_interpretation`, `snakemake -n level2_within_group_trees`, `snakemake -n level3_candidate_selection`, `snakemake -n level3_segment_cn`, `snakemake -n level3_mutation_local_cn_annotation`, `snakemake -n level3_clonal_input_preparation`, `snakemake -n level3_limited_vaf_clonal_clustering`, `snakemake -n level3_allele_specific_cn_readiness`, `snakemake -n tcga_sample_reference_atlas`, `snakemake -n clinical_case_projection_example`, `snakemake -n cup_interpretation_example`, `snakemake -n known_primary_projection_validation`, `snakemake -n calibration_feature_ablation`, `snakemake -n integrated_project_report`, and `snakemake -n preprint_manuscript` to validate workflow wiring outside the current shell.
3. Review the locked atlas feature contract and `results/clinical_projection/example_cup_case/clinical_tcga_projection_report.md` before adding additional clinical research queries.
4. If adding manually downloaded segment files, place them under `data/raw/copy_number_segments/` or `data/raw/copy_number/` and rerun `python scripts/09b_download_segment_level_cn.py`.
5. Extend the completed internal TCGA pseudo-unknown validation and within-TCGA confidence calibration to independent known-primary clinical WES/WTS cases, including assay-stratified feature coverage, realistic partial-WTS simulations, top-k and distance-metric sensitivity, patient-level exclusion, blinded molecular-pathologic concordance review, and external evaluation of abstention thresholds. Do not train a supervised CUP classifier until these checks are complete.
6. Build a separately validated actionability integration layer after molecular-context performance is benchmarked; OncoKB/CIViC/AMP interpretation is not implemented in the current prototype.
7. Use `results/tables/level3_allele_specific_cn_external_inference_manifest.tsv` to obtain paired tumor/normal or precomputed allele-count inputs for the 30 limited-VAF pilot samples first. Defer PyClone-VI and PhyloWGS execution until standardized integer major/minor CN, mutation overlap, coordinate build, purity, and per-sample input schemas pass review.

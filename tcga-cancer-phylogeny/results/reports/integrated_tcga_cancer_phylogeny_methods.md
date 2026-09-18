# Integrated Methods Summary

## Data Acquisition

The workflow retrieves or reads TCGA project metadata, MC3 somatic mutation calls, PanCanAtlas/GDC purity and ploidy resources, aneuploidy and arm-level copy-number calls, GDC segment-level copy-number files for Level 3 pilot samples, and PanCanAtlas expression data. Local source-file fallbacks are used where possible, and missing-source states are recorded explicitly in QC tables.

## Feature Engineering and Reference Atlas

Mutation features include mutation-count proxies and configured driver-gene mutation indicators and counts. Copy-number features include purity, ploidy, aneuploidy score, and arm gain/loss burden. Expression features include transparent first-pass pathway scores and expression principal components. Level 1 uses project-level aggregated traits. The locked sample-level atlas retains 97 features, uses TCGA-reference median imputation and standard scaling, and saves feature definitions, transformation parameters, PCA, and Euclidean/cosine nearest-neighbor artifacts.

## Level 1 Tree Construction

Level 1 reads the project-level feature matrix, removes annotation fields from tree calculations, median-imputes remaining numeric missing values, scales biologic features, and generates Gower and correlation average-linkage molecular similarity dendrograms. Bootstrap feature resampling summarizes project-pair co-clustering stability.

## Level 2 Group Definitions

Level 2 defines biologically motivated project groups, including carcinoma_all, pan_squamous, gi_pancancreatobiliary, kidney, and gynecologic_breast. It filters features for within-group variability, scales features within group, and builds Gower and correlation trees where project counts are sufficient. Underpowered groups are retained in QC tables rather than forced into trees.

## Level 3 Candidate Selection

Level 3 candidate selection evaluates mutation counts, ref/alt count availability, purity, ploidy, copy-number availability, and candidate project priorities. It selects a high-information pilot cohort for clonal-structure prototyping.

## Segment CN Processing

The GDC-backed segment-CN step queries open-access TCGA copy-number segment files for pilot projects, matches files to pilot samples where possible, selects one best sample-level file per pilot sample, parses segment intervals, and updates candidate/pilot copy-number availability tables.

## Mutation-to-Segment Annotation

The mutation local-CN annotation step overlaps somatic mutations with sample-level segment intervals, computes depth and observed VAF, records local segment coordinates and segment mean, and writes a non-final PyClone-style preview. Segment mean is treated conservatively and is not allele-specific integer copy number.

## Limited VAF Clustering

The limited VAF clustering prototype uses copy-neutral candidate mutations from selected pilot samples. It fits one-dimensional Gaussian mixture models to observed VAFs per sample, selects cluster counts by BIC with a conservative simpler-model tolerance, and falls back to quantile binning when needed. Cluster labels are descriptive VAF-structure summaries.

## Clinical WES/WTS Projection

The clinical projection module ingests already-computed harmonized feature tables rather than raw FASTQ, BAM, VCF, or RNA-seq data. Exact locked feature names are matched, absent features are filled with saved TCGA medians, and saved TCGA means and standard deviations are applied without refitting. The saved PCA model is reused, and nearest neighbors are queried from the compatible saved index with direct scaled-matrix distances as a fallback.

## CUP Molecular Interpretation

The CUP layer validates the top-k project and broad-group summaries, aggregates project weights into predefined cautious lineage groups, and computes top-lineage share, top-1/top-2 margin, normalized entropy, neighbor concentration, feature coverage, and broad-group agreement. Configurable heuristic rules summarize project, major-group, expression, driver, copy-number, microenvironment, contradictory, and missing evidence. Differential diagnosis and molecular-pathologic comparisons are keyword-based aids, not diagnostic inference.

## Allele-Specific CN Readiness

The readiness workflow searches supported local allele-specific CN outputs, validates integer total/major/minor copy-number fields, and measures mutation overlap. It does not infer major/minor CN from segment mean. No PyClone-VI or PhyloWGS execution is performed; those tools remain gated on reviewed allele-specific integer CN, purity, mutation read counts, and sufficient mutation overlap.

## Statistical Caveats

Mutation burden remains a count proxy rather than true mutations/Mb. Expression pathway scores are simple gene-set summaries. Levels 1 and 2 are project-level molecular similarity trees rather than literal evolutionary phylogenies. Level 3 does not infer allele-specific integer CN or branching order, and single bulk samples limit clonal-evolution inference. Clinical and CUP scores are descriptive similarity summaries, not probabilities, and have not been clinically calibrated.

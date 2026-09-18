# Limitations and Next Steps

## Key Limitations

- TCGA is mostly primary tumors and does not fully represent metastatic, relapse, small-biopsy, or cytology clinical specimens.
- Levels 1 and 2 summarize project-level molecular traits, not within-patient phylogeny.
- Level 1 and Level 2 trees are molecular similarity trees and should not be interpreted as literal species-like evolutionary phylogenies.
- Mutation burden is currently a mutation-count proxy without callable-territory normalization.
- Segment mean is not allele-specific major/minor integer copy number.
- Level 3 lacks allele-specific CN and therefore does not produce CN-aware PyClone-VI or PhyloWGS inputs.
- Limited VAF clustering cannot infer definitive branching order.
- Single bulk TCGA samples limit robust clonal-evolution inference.
- The example clinical query is synthetic, supplies only 11 of 31 WTS features, and median-imputes 20 expression PCs.
- Clinical projection and CUP scores are descriptive similarity weights, not tissue-of-origin probabilities.
- CUP confidence and evidence rules are heuristic and uncalibrated, and no supervised CUP classifier was trained.
- Clinical actionability is not implemented.
- Clinical use requires independent known-primary validation, assay harmonization, uncertainty calibration, and prospective reporting safeguards.

## Next Recommended Steps

1. Validate projection and CUP interpretation on independent known-primary WES/WTS cohorts, including metastatic and limited-specimen subsets.
2. Calibrate top-k, distance-metric, ambiguity, and abstention behavior before any supervised classifier development.
3. Add a separately validated actionability layer after the molecular-context workflow is benchmarked.
4. Obtain paired tumor/normal or validated allele-count inputs for the 30 limited-VAF pilot samples and run one allele-specific CN method consistently.
5. Prepare PyClone-VI inputs only after reviewed integer major/minor CN and mutation-overlap QC are available; consider PhyloWGS later on a high-confidence subset.
6. Create the pending conceptual manuscript figure, format references, and obtain external scientific and clinical review.

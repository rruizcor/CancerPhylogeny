# TCGA Multilevel Molecular Similarity, Clonal-Structure, and Clinical Projection Framework: Executive Summary

This project provides an end-to-end computational proof of concept spanning TCGA pan-cancer molecular similarity, group-specific molecular trait dendrograms, limited within-tumor VAF clustering, a locked sample-level reference atlas, and research-only clinical WES/WTS and CUP projection.

Level 1 builds pan-cancer molecular similarity trees across 33 TCGA cancer types using 52 biologic features. The Gower and correlation trees provide complementary views: Gower emphasizes absolute scaled trait differences, while correlation emphasizes relative feature-profile shape.

Level 2 builds group-specific molecular similarity trees. Trees were generated for carcinoma_all, pan_squamous, gi_pancancreatobiliary, kidney, gynecologic_breast. CNS/glial, melanocytic, hematolymphoid, and sarcoma project-level analyses remain underpowered in the current project-level design.

Level 3 evaluates within-tumor clonal structure in a limited prototype. Candidate selection screened 10201 tumor samples, selected 84 pilot samples, prepared copy-neutral mutation sets, and clustered 81154 mutations from 30 selected samples. Interpretation classes are 5 predominantly clonal-like, 11 oligoclonal-like, and 14 multicluster subclonal-like.

The main scientific caveat is that Levels 1 and 2 are molecular similarity and trait-mapping trees, not literal organismal phylogenies. Level 3 is closer to tumor evolution, but the current implementation is limited VAF-based clustering in copy-neutral candidate regions. It does not run PyClone-VI or PhyloWGS, does not infer allele-specific integer copy number, and does not infer definitive branching order.

The hardened sample-level atlas contains 10201 tumors from 10130 patients across 33 projects. It retains all 97 features, records 23222 median-imputed modeling values, and provides locked PCA and Euclidean/cosine nearest-neighbor artifacts.

The artificial `example_cup_case` supplied 77 of 97 features (79.381% coverage). Its top CUP lineage was `lower_GI` with score 0.5701 and `moderate_ambiguity`. These are descriptive similarity weights from a synthetic workflow demonstration, not diagnostic probabilities.

PyClone-VI and PhyloWGS remain paused: 0 and 0 samples, respectively, are ready because allele-specific integer major/minor CN is unavailable. The next priorities are known-primary clinical validation, later actionability integration, and external allele-specific CN acquisition for the Level 3 pilot.

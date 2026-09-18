# Level 2 Within-Group Molecular Similarity Summary

Level 2 builds group-specific molecular similarity trees from the same project-level molecular traits used in Level 1. Features are filtered for within-group variability, median-imputed within group when needed, scaled within group, and clustered with Gower and correlation distances using average linkage.

Groups generating project-level trees were: carcinoma_all, pan_squamous, gi_pancancreatobiliary, kidney, gynecologic_breast. Underpowered groups were: cns_glial (underpowered_project_count_2_lt_3); melanocytic (underpowered_project_count_2_lt_3); hematolymphoid (underpowered_project_count_2_lt_3); sarcoma (underpowered_project_count_1_lt_3;insufficient_variable_features_0_lt_2;no_clinical_or_sarc_subtype_feature_metadata_found).

The carcinoma-all Gower tree shows the following k-level structure: cluster 1: BLCA,BRCA,CESC,ESCA,HNSC,LUAD,LUSC,OV,STAD; cluster 2: CHOL,KICH,KIRP,LIHC,PRAD,THCA; cluster 3: COAD,READ; cluster 4: KIRC; cluster 5: PAAD; cluster 6: UCEC. Close project pairs include PRAD-THCA, COAD-READ, KIRP-LIHC, KICH-PRAD, KICH-THCA, BRCA-OV, LUAD-LUSC, HNSC-STAD. These results indicate similar molecular trait profiles within the carcinoma set and should not be interpreted as directional ancestry.

The carcinoma-all correlation tree differs as expected because it emphasizes relative feature-profile shape rather than absolute scaled trait levels. Its k-level structure is: cluster 1: BLCA,CESC,ESCA,HNSC,LUAD,LUSC,OV,STAD; cluster 2: BRCA; cluster 3: CHOL,KICH,KIRP,LIHC; cluster 4: COAD,READ; cluster 5: KIRC,PAAD,PRAD,THCA; cluster 6: UCEC.

Pan-squamous / squamous-enriched produced a group-specific molecular trait tree with Gower clusters: cluster 1: LUSC; cluster 2: HNSC,ESCA; cluster 3: CESC; cluster 4: BLCA. GI / pancreatobiliary adenocarcinoma-enriched produced a group-specific molecular trait tree with Gower clusters: cluster 1: COAD,READ,STAD; cluster 2: ESCA; cluster 3: PAAD; cluster 4: CHOL,LIHC. Kidney produced a group-specific molecular trait tree with Gower clusters: cluster 1: KIRC,KIRP; cluster 2: KICH. Gynecologic / breast produced a group-specific molecular trait tree with Gower clusters: cluster 1: BRCA,OV; cluster 2: UCEC; cluster 3: CESC; cluster 4: UCS.

Across Level 2 trees, clustering appears influenced by a mixture of mutation-count proxies, driver-gene prevalence, aneuploidy/copy-number burden, purity/ploidy, and expression programs. Because Level 1 sensitivity showed mutation-heavy features had strong cluster associations, exact Level 2 topology should be read together with feature-subset and weighting sensitivity rather than treated as a single definitive structure.

CNS/glial, melanocytic, hematolymphoid, and sarcoma project-level analyses are limited by having one or two TCGA projects in the relevant group. The SARC project lacks a project-internal subtype molecular feature matrix in the current workflow, so a subtype-level sarcoma tree is deferred.

Scientific caveat: Level 2 outputs are molecular similarity dendrograms or molecular trait trees. They are not literal organismal phylogenies and should not be used to infer directional evolutionary relationships between cancer types.


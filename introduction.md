---
title: 'Introduction'
teaching: 10
exercises: 2
---

# scRNA-seq Integration and Differential Expression

[**Slides available here.**](files/Seurate-DE-slides-NEW.pdf)

**Author:** Manveer Chauhan, Clark Lab, The University of Melbourne.\
**Contributors:** Vini Salazar, Melbourne Bioinformatics.

Last updated October 2024.

## Overview

**Topic**

-   [ ] Genomics
-   [x] Transcriptomics
-   [ ] Proteomics
-   [ ] Metabolomics
-   [ ] Statistics and visualisation
-   [ ] Structural Modelling
-   [ ] Basic skills

**Skill level**

-   [ ] Beginner\
-   [x] Intermediate\
-   [ ] Advanced

**Data:** IFNB-Stimulated and Control PBMCs.

**Tools:** R \>=4.4.0 and associated packages:

-   Seurat\
-   SeuratData\
-   tidyverse\
-   DESeq2\
-   patchwork\
-   pheatmap\
-   grid\
-   metap

**Pipeline:**\
*Section 1:* Setup, Quality Control and Sample Integration.\
*Section 2:* Differential Gene Expression when dealing with two
treatment conditions.\
*Section 3:* Differential Expression using a pseudobulk approach and
DESeq2.

**Learning objectives:**

-   Gain more familiarity with standard scRNA-seq Quality Control (QC)
    steps
-   Understand and get comfortable using various integration strategies
    (harmoy and seuratCCA)
-   Understand when and how to use all of the differential expression
    functions offered by Seurat: FindMarkers(), FindConservedMarkers(),
    and FindAllMarkers()
-   Learn how to use differential expression tools meant for bulk data,
    like DESeq2, for single-cell 'pseudobulk' data and understand why
    you might choose this approach.
-   Learn different ways to visualise both in-built Seurat functions and
    external packages like pheatmap.


:::: callout

This tutorial is partially based on existing
material from:

* https://satijalab.org/seurat/articles/seurat5_integration
* https://satijalab.org/seurat/articles/de_vignette
* https://hbctraining.github.io/scRNA-seq_online/
* https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/

::::

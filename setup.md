---
title: Setup
---


## Data Sets


:::: prereq

Some knowledge of R and scRNA-seq analysis is assumed. We recommend reviewing the following materials before the starting the materials:
- [Introduction to R](https://mbite.org/intro-to-r/)
- [10X single-cell RNA-seq analysis in R](https://mbite.org/tutorials/singlecell/)

This lesson assumes you have the latest versions of R and RStudio installed on your computer.


::::


## Software Setup

::::::::::::::::::::::::::::::::::::::: discussion

If you don't have R and RStudio already installed, please download them here:

- [Download and install the latest version of R using the UniMelb mirror](https://cran.ms.unimelb.edu.au/).
- [Download and install RStudio](https://posit.co/download/rstudio-desktop/#download).


:::::::::::::::::::::::::::::::::::::::::::::::::::


**Run the code block below to install the packages needed for this
workshop.**

To check if installed properly, load each package in one at a time using
the `library()` function.


``` r
install.packages('Seurat')

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("satijalab/seurat-data", quiet = TRUE)

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")
BiocManager::install("multtest") # dependency commonly missing
BiocManager::install(c("SingleR", "celldex"))
install.packages("harmony") # dependency needed for harmony analysis

install.packages("tidyverse")
install.packages("pheatmap")
install.packages("metap")
install.packages("ggplot2")


```


---
title: Setup
---


## Data Sets


/data/ifnb.RData


:::: prereq

Some knowledge of R and scRNA-seq analysis is assumed.

This lesson assumes you have R and RStudio installed on your computer.


::::


## Software Setup

::::::::::::::::::::::::::::::::::::::: discussion

If you don't have R and RStudio already installed, please download them here:

[Download and install the latest version of R using the UniMelb mirror](https://cran.ms.unimelb.edu.au/).
[Download and install RStudio](https://posit.co/download/rstudio-desktop/#download).


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
BiocManager::install("DropletUtils")

install.packages("harmony") # dependency needed for harmony analysis

install.packages("tidyverse")
install.packages("pheatmap")
install.packages("metap")
install.packages("ggplot2")
```


---
title: 'Section 2: Differential Gene Expression when dealing with two treatment conditions'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do conserved markers help us label clusters reliably across conditions?
- What exactly do avg_log2FC, pct.1, pct.2, and p_val_adj mean in FindMarkers?
- Why must DE be run within a cell type (e.g., CD16 Mono_STIM vs CD16 Mono_CTRL) rather than “all cells”?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Use FindConservedMarkers() to pick markers and label clusters.
- Set identities to annotations and create compound identities (celltype.and.stim) for clean contrasts.
- Run FindMarkers() to get DEGs between conditions within a cell type and interpret key columns.
- Visualize DEGs (FeaturePlot with split.by, DoHeatmap / pheatmap) and export results for downstream use.
- Recognize pitfalls (composition effects, inappropriate contrasts, overly lenient thresholds).

::::::::::::::::::::::::::::::::::::::::::::::::


``` error
Error in library(DESeq2): there is no package called 'DESeq2'
```



``` r
DimPlot(ifnb.filtered, reduction = "umap.cca", label = T)
```

``` error
Error: object 'ifnb.filtered' not found
```



``` r
DimPlot(ifnb.filtered, reduction = "umap.cca", group.by = "stim")
```

``` error
Error: object 'ifnb.filtered' not found
```


## Step 1. Find Conserved Markers to label our celltypes


``` r
## Let's look at conserved markers in cluster 4 across our two conditions (compared to all other clusters)
markers.cluster.4 <- FindConservedMarkers(ifnb.filtered, ident.1 = 4,
                     grouping.var = 'stim')
```

``` error
Error: Please install the metap package to use FindConservedMarkers.
This can be accomplished with the following commands: 
----------------------------------------
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
----------------------------------------
```

``` r
head(markers.cluster.4)
```

``` error
Error: object 'markers.cluster.4' not found
```

Let's visualise the top upregulated, conserved between condition,
marker genes identified in cluster 4 using the `FeaturePlot()` function.


::::: discussion

Try running the function in the code block below without
defining a min.cutoff, or changing the value of the min.cutoff
parameter.

:::::



``` r
# Try looking up some of these markers here: https://www.proteinatlas.org/
FeaturePlot(ifnb.filtered, reduction = "umap.cca", 
            features = c('VMO1', 'FCGR3A', 'MS4A7', 'CXCL16'), min.cutoff = 'q10')
```

``` error
Error: object 'ifnb.filtered' not found
```




``` r
# I happen to know that the cells in cluster 4 are CD16 monocytes - lets rename this cluster
# Idents(ifnb.filtered) # Let's look at the identities of our cells at the moment
```



``` r
ifnb.filtered <- RenameIdents(ifnb.filtered, '4' = 'CD16 Mono') # Let's rename cells in cluster 4 with a new cell type label
```

``` error
Error: object 'ifnb.filtered' not found
```

``` r
# Idents(ifnb.filtered) # we can take a look at the cell identities again
```



``` r
DimPlot(ifnb.filtered, reduction = "umap.cca", label = T) +
  ggtitle("After changing the identity of cluster 4")
```

``` error
Error: object 'ifnb.filtered' not found
```


## Step 2: Set the identity of our clusters to the annotations provided


``` r
Idents(ifnb.filtered) <- ifnb.filtered@meta.data$seurat_annotations
```

``` error
Error: object 'ifnb.filtered' not found
```

``` r
# Idents(ifnb.filtered)
```



``` r
DimPlot(ifnb.filtered, reduction = "umap.cca", label = T)
```

``` error
Error: object 'ifnb.filtered' not found
```


:::: callout

If you want to perform cell-type identification on your own data
when you don't have a ground-truth, using automatic cell type annotation
methods can be a good starting point. This approach is discussed in more
detail in the Intro to scRNA-seq workshop material.

::::

::::::::::::::::::::::::::::::::::::: challenge 
Automated Cell Type Annotation

:::::::::::::::::::::::: solution 


``` r
# Load reference data 
# Blood & immune lineages
ref.set <- celldex::BlueprintEncodeData()
```

``` error
Error in loadNamespace(x): there is no package called 'celldex'
```

``` r
ifnb.v4 <- JoinLayers(ifnb.filtered)
```

``` error
Error: object 'ifnb.filtered' not found
```

``` r
sce.ifnb.filtered <- as.SingleCellExperiment(ifnb.v4)
```

``` error
Error: object 'ifnb.v4' not found
```

``` r
sce.ifnb.filtered <- logNormCounts(sce.ifnb.filtered)
```

``` error
Error in logNormCounts(sce.ifnb.filtered): could not find function "logNormCounts"
```

``` r
pred.cnts <- SingleR(
  test = sce.ifnb.filtered,
  ref = ref.set,
  labels = ref.set$label.main
)
```

``` error
Error in SingleR(test = sce.ifnb.filtered, ref = ref.set, labels = ref.set$label.main): could not find function "SingleR"
```

``` r
lbls.keep <- table(pred.cnts$labels)>10
```

``` error
Error: object 'pred.cnts' not found
```

``` r
# Add SingleR labels to Seurat metadata
ifnb.filtered$SingleR.labels <- sce.ifnb.filtered$SingleR.labels
```

``` error
Error: object 'sce.ifnb.filtered' not found
```

``` r
ifnb.filtered$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')
```

``` error
Error: object 'lbls.keep' not found
```

``` r
# Run UMAP (based on PCA)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20)
```

``` error
Error: object 'ifnb.filtered' not found
```

``` r
DimPlot(ifnb.filtered, reduction='umap.cca', group.by='SingleR.labels',  label = TRUE, label.size = 3 )
```

``` error
Error: object 'ifnb.filtered' not found
```

:::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::

## Step 3: Find differentially expressed genes (DEGs) between our two conditions, using CD16 Mono cells as an example


``` r
# Make another column in metadata showing what cells belong to each treatment group (This will make more sense in a bit)
ifnb.filtered$celltype.and.stim <- paste0(ifnb.filtered$seurat_annotations, '_', ifnb.filtered$stim)
```

``` error
Error: object 'ifnb.filtered' not found
```

``` r
# (ifnb.filtered@meta.data)

Idents(ifnb.filtered) <- ifnb.filtered$celltype.and.stim
```

``` error
Error: object 'ifnb.filtered' not found
```

``` r
DimPlot(ifnb.filtered, reduction = "umap.cca", label = T) # each cluster is now made up of two labels (control or stimulated)
```

``` error
Error: object 'ifnb.filtered' not found
```




``` r
DimPlot(ifnb.filtered, reduction = "umap.cca", 
        label = T, split.by = "stim") # Lets separate by condition to see what we've done a bit more clearly
```

``` error
Error: object 'ifnb.filtered' not found
```

We'll now leverage these new identities to compare DEGs between our
treatment groups


``` r
treatment.response.CD16 <- FindMarkers(ifnb.filtered, ident.1 = 'CD16 Mono_STIM', 
                                       ident.2 = 'CD16 Mono_CTRL')
```

``` error
Error: object 'ifnb.filtered' not found
```

``` r
head(treatment.response.CD16) # These are the genes that are upregulated in the stimulated versus control group
```

``` error
Error: object 'treatment.response.CD16' not found
```

## Step 4: Lets plot conserved features vs DEGs between conditions


``` r
FeaturePlot(ifnb.filtered, reduction = 'umap.cca', 
            features = c('VMO1', 'FCGR3A', 'IFIT1', 'ISG15'),
            split.by = 'stim', min.cutoff = 'q10')
```

``` error
Error: object 'ifnb.filtered' not found
```


## Step 5: Create a Heatmap to visualise DEGs between our two conditions + cell types


``` r
# Find upregulated genes in each group (cell type and condition)
ifnb.treatVsCtrl.markers <- FindAllMarkers(ifnb.filtered,
                                          only.pos = TRUE)
```

``` error
Error: object 'ifnb.filtered' not found
```


``` r
saveRDS(ifnb.treatVsCtrl.markers, "ifnb_stimVsCtrl_markers.rds")
```



If the top code block takes too long to run - you can download the rds
file of the output using the code below:

Seurat's in-built heatmap function can be quite messy and hard to
interpret sometimes (we'll learn how to make better and clearer custom
heatmaps using the pheatmap package from our Seurat expression data
later on).


``` r
ifnb.treatVsCtrl.markers <- readRDS(url("https://github.com/manveerchauhan/Seurat_DE_Workshop/raw/refs/heads/main/ifnb_stimVsCtrl_markers.rds"))

top5 <- ifnb.treatVsCtrl.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

DEG.heatmap <- DoHeatmap(ifnb.filtered, features = top5$gene,
          label = FALSE)
```

``` error
Error: object 'ifnb.filtered' not found
```

``` r
DEG.heatmap
```

``` error
Error: object 'DEG.heatmap' not found
```

::::::::::::::::::::::::::::::::::::: keypoints 
- QC filtering removes low-quality cells (e.g., low gene count or high mitochondrial %).
- Integration corrects sample-to-sample variation so cells group by biology, not by batch.
- Harmony and CCA both align shared cell states but use different mathematical strategies.

::::::::::::::::::::::::::::::::::::::::::::::::

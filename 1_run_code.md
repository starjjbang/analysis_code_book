Untitled
================

git commit -am ‘date’ git push origin +master

## merge vs integration

<https://github.com/satijalab/seurat/issues/1787>

-   merge just concatenates the counts or data table from two object
    together.

-   if either object has unique genes, it will be added in the merged
    objects.

-   IntegrateData is a step in the integration, integrating two objects
    by anchor(matched biological state (‘anchors’)) cells.

-   The detail you could find in the paper, here.

-   You should only use merge for technical replicates, and in theory
    for a group of samples with a low batch effect.

-   Integration in Seurat (and related) was developed because there
    tends to be a relatively strong batch in the manifolds.

-   By this I mean that even if two cell populations are the same
    between two samples,

-   they will appear as two partially or fully separate clusters.

-   Integration tries to “smooth out” the differences in batches so that
    cells that are likely similar will cluster together.

-   As with anything single-cell related I would suggest exploring and
    validating the results from both merging and integration  
    to get a better feel for your data.

## recomand pipeline

Merge samples into one Seurat object Make QC plots, split by sample.
Merged object allows to easily examine different features on one plot
Filter cells and genes, depending on batch structure (see below). Split
merged object into samples again Integrate UMAP, clustering, etc.

1.  Create all seurat objects.
2.  Perform the quality-check and filtering for each one of them.
3.  Calculate percent.mito.
4.  Normalize each dataset separately with SCTransform. Regress out
    percent.mito if desired.
5.  Calculate cell cycle scores.
6.  Prepare for integration with SelectIntegrationFeatures() and
    PrepSCTIntegration().
7.  Integrate datasets.
8.  ScaleData, regressing out cell cycle score if desired
9.  Run PCA, UMAP, FindClusters, FindNeighbors on integrated assay.
10. Switch to “SCT” assay and continue with the DE analysis.

\##merge() conserved the original raw data of the objects; and
IntegrateData() normalized the raw data?

## path & package

    q_dir = '/home/starjjbang/project/1_prof_kim_lab_deg/1_single_cell_project/'
    fig_dir = paste0(q_dir, 'plot/7_0725_fig/',in_fod_dir)
    rds_dir = paste0(q_dir, 'rds/4_220725/',in_fod_dir)
    ref_dir = paste0(q_dir, 'result/5_0725_analysis/',in_fod_dir)
    system(paste0('mkdir ',fig_dir))
    system(paste0('mkdir ',rds_dir))
    system(paste0('mkdir ',ref_dir))
    source('/home/starjjbang/project/1_prof_kim_lab_deg/1_single_cell_project/code/single_function_ver02_0412.R')

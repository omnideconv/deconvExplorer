# Data Upload

**All files can be provided in csv, tsv and rds format. Vectors can be
uploaded as txt as well.**

#### Deconvolution Data

-   Bulk RNA-seq data
    -   *Genes* x *Samples* matrix
    -   Rownames (gene names) are provided in the same format as in the
        sc RNA-seq data, for instance HGNC symbols
-   Single cell RNA-seq data
    -   *Genes* x *Cells* matrix
    -   Counts are *not log-transformed*
    -   Rownames (gene names) are provided in the same format as in the
        bulk RNA-seq data, for instance HGNC symbols
-   Cell type annotations
    -   Vector containing cell type annotations
    -   Annotations are in the same order as the columns of the single
        cell matrix
-   Batch ids
    -   Vector containing batch ids, so sample or patient ids
    -   Ids are in the same order as the columns of the single cell
        matrix
    -   This is only necessary for Bisque, MuSiC and SCDC
-   (Marker genes)
    -   Vector containing gene names
    -   This is only necessary for BSeq-sc

#### SimBu Simulation

If no ground truth data is available for your bulk dataset you can
benchmark by simulating a pseudo-bulk sample with known cell type
fraction using SimBu. To transfer your simulation to DeconvExplorer save
it in rds format and upload it to retrieve the simulated bulk file and
the corresponding ground truth.

    simulation <- SimBu::simulate_bulk(...)
    saveRDS(simulation, 'filepath.rds') # upload this file

#### Custom Signature and Ground Truth Reference

Please make sure the first column contains gene identifiers matching the
bulk sample.

# Deconvolution

Choose a deconvolution method to run a deconvolution with omnideconv. In
total omnideconv offers 11 deconvolution algorithms. Each method builds
a custom gene expression signature internally or, if applicable, you can
choose to provide a custom one. The results can be visualized below and
you can select multiple deconvolutions to be plotted at the same time.
In addition you can inspect the results and available signatures in
table form further below.

# Signature Exploration

This module can visualize and compare the available signatures. Multiple
metrics and plots are available. The prominently placed heatmap is
interactive and allows you to select fragments of a signature for a more
in depth analysis in the plot and the table below. To compare the gene
sets of multiple signatures the UpSet Plot at the end of this module can
be utilized.

# Signature Refinement

If required this module offers four functions to further subselect genes
and enhance signature quality. Before a refinement a signature should be
loaded on the right side and needs to be saved by a new name afterwards.
For benchmarking purposes it is also possible to rename cell types.

1.  Remove Sequencing Artifacts Too many zeros in expression values of a
    gene across cell types can indicate a sequencing error. For this
    refinement you can choose the maximum percentage of zeros to be
    allowed for each gene. All other genes get discarded. It is
    recommended to use this function early in the refinement process as
    genes expressed in this manner get scored extraordinary high be the
    genes scores applied in the third function.

2.  Remove genes expressed not specific enough This function works by
    dividing the expression range of each gene in three equal parts and
    assigning the expression values to one of the bins “high”, “medium”
    and “low”. With the parameter (default=1) you can select the maximum
    number of cell types a gene is allowed to be expressed in the “high”
    bin. For example, when running with the parameter set to 1 each gene
    is only allowed to be expressed in the highest bin in only 1 cell
    type. All other cell types need to fall into “medium” or “low”
    categories for this gene. All genes not fulfilling this requirement
    are removed from the signature.

3.  Select best scored genes for each cell type For each cell type
    select a fixed amount of best scored genes. As scores Entropy and
    Gini Index can be applied. Both scores can be seen as heatmap
    annotation in the signature exploration and refinement module. The
    scores are used to quantify how cell type specific the expression of
    a gene is. The expression of a gene mainly in one cell type will
    result in low entropy values and a high Gini index.

# Benchmark

To test deconvolution performance and compare methods or signatures the
benchmarking module offers three different plots. For benchmarking a
ground truth of the deconvoluted bulk sample is necessary and can be
obtained from a SimBu simulation or by directly uploading a reference
file in the Data Upload module. While only one reference file can be
chosen multiple deconvolution results can be benchmarked at the same
time.

This module offers a scatterplot visualizing the true vs estimated cell
fractions. The result is better if it appears closer to the black line
(identity line). In addition the correlation between the datasets and a
Root Mean Square Error (RMSE) heatmap or boxplot can be displayed.

# Further information

The code for omnideconv, SimBu and DeconvExplorer can be found on github
and is linked above.

An interactive tour can be started in the menu bar at the top.

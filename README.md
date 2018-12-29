
eda - Exploratory Data Analysis for R
=====================================

The `eda` package for R provides general purposes, as well as domain-specific classes for performing basic exploratory data analysis in R. It is built on the [R6 class system](https://cran.r-project.org/web/packages/R6/index.html), as well as a number of other powerful statistics, machine learning, and visualization packages.

The primary goals of `eda` are to provide a simple interface for exploring one or more related datasets, with an emphasis on:

-   Quick and easy summary statistics and visualization
-   Integration of multiple related datasets
-   Streamlined data transformation and visualization using function chaining
-   Generalized base classes that can be extended for domain-specific application

This package is still in the early stages of development. While a reasonable overall class hierarchy has been implemented, significant work is still required with respect to normalization of function calls (especially those relating to plotting), and documentation.

Installation
============

You can install `eda` using Bioconductor with:

``` r
BiocManager::install('github.com/khughitt/eda')
```

Example Usage
=============

I still need to write a fair bit of documentation for how to actually use this package. That should happen soon (probably Jan 2019).

In the meantime, here are some example use cases to give you a sense for how this package can be used.

The below examples make use of the high-throughput biology (transcriptomic) data made available through the [recount2](http://bioconductor.org/packages/release/bioc/html/recount.html) package.

``` r
library(eda)
library(recount)

# download TCGA RNA-Seq data from ReCount
download_study('TCGA')
```

    ## Loading eda

``` r
# make output reproducible
set.seed(1)

# load TCGA RangedSummarizedExperiment
load(file.path('TCGA', 'rse_gene.Rdata'))

# to speed things up, let's grab a random subsample of the data;
# this can also be done within eda, but there is currently an issue
# relating to transposition of subsampled data, so for now, we will
# handle this externally
SAMPLE_NROWS <- 1000
SAMPLE_NCOLS <- 100

rse_gene <- rse_gene[sample(nrow(rse_gene), SAMPLE_NROWS), sample(ncol(rse_gene), SAMPLE_NCOLS)]

# convert RSE object to an BioDataSet instance
bdat <- BioDataSet$new(rse_gene)

# the expression data and gene / sample metadata are eached stored as
# separate dataset in the BioDataSet object
bdat
```

    ## =========================================
    ## =
    ## = BioDataSet (n=3)
    ## =
    ## =   assays : matrix (1000 x 100)
    ## =  rowdata : DataFrame (1000 x 3)
    ## =  coldata : DataFrame (100 x 864)
    ## =
    ## =========================================

``` r
# each of these can be accessed through the "datasets" property
bdat$datasets$assays[1:3, 1:3]
```

    ##                   AEC4FDAA-8F60-43CE-8EBB-5629FA662F16
    ## ENSG00000211680.2                                    0
    ## ENSG00000235760.3                                    0
    ## ENSG00000273129.1                                  904
    ##                   655B32BE-AC22-4E6C-97D2-CC6B7F06B55E
    ## ENSG00000211680.2                                    0
    ## ENSG00000235760.3                                   29
    ## ENSG00000273129.1                                    0
    ##                   BDBF3ACE-7DC0-47C1-9B03-E2EDCF0B6CE6
    ## ENSG00000211680.2                                    0
    ## ENSG00000235760.3                                  304
    ## ENSG00000273129.1                                   24

``` r
bdat$datasets$coldata[1:3, 1:3]
```

    ## DataFrame with 3 rows and 3 columns
    ##                                          project    sample experiment
    ##                                      <character> <logical>  <logical>
    ## AEC4FDAA-8F60-43CE-8EBB-5629FA662F16        TCGA        NA         NA
    ## 655B32BE-AC22-4E6C-97D2-CC6B7F06B55E        TCGA        NA         NA
    ## BDBF3ACE-7DC0-47C1-9B03-E2EDCF0B6CE6        TCGA        NA         NA

``` r
head(bdat$datasets$rowdata, 3)
```

    ## DataFrame with 3 rows and 3 columns
    ##                             gene_id bp_length          symbol
    ##                         <character> <integer> <CharacterList>
    ## ENSG00000211680.2 ENSG00000211680.2        33              NA
    ## ENSG00000235760.3 ENSG00000235760.3       762      HCG2040054
    ## ENSG00000273129.1 ENSG00000273129.1       825          PACERR

``` r
# subsample data, log-transform, and plot a heatmap
#bdat$subsample(row_n = 500, col_n = 100)$log1p()$plot_heatmap(interactive = FALSE)

# generate a sample PCA, t-SNE, and UMAP plots, colored by cancer 
# TODO: currently an issue with RSE ingestion / handling of colors.. need to look
# into this when I have more time..
#bdat$t()$plot_pca(color_var = 'gdc_cases.project.primary_site', color_key = 'coldata')
```

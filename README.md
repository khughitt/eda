eda - Exploratory Data Analysis for R
=====================================

[![Build Status](https://travis-ci.org/khughitt/eda.svg?branch=master)](https://travis-ci.org/khughitt/eda)

**Note Feb 26, 2018**: This library is currently in the early stages of development
and is likely to change significantly over time. Feedback and suggestions are
welcome!

Overview
--------

`eda` is a package containing helper classes to assist with exploratory data 
analysis work. Object instances are constructed from datasets, including 
optional row and column metadata, and methods are provided for computing useful
summary statistics and visualizations. The overall focus of the package is on
simplicity, with the aim of allowing users to quickly get a feel for a new 
dataset with as few lines of code as possible. As such, some assumptions are
made about the types of visualizations to be made, and customization of plots
is limited compared to those created by hand. In order to support as wide a
range of datasets as possible, a class hierarchy has been created, including
both generic data classes (`EDAMatrix` and `EDADataFrame`), as well as 
domain-specific classes (e.g. `BioExprSet`). Users also have the ability to
inherit from any of the provided classes (which are built on top of the 
[R6 class framework](https://cran.r-project.org/web/packages/R6/vignettes/Introduction.html))
in order to construct new domain-specific classes.

The basic assumptions of `eda` are that you have a dataset stored in something
like a matrix or dataframe, and optionally, additional metadata about the rows
and columns of the dataset, each stored as separate dataframes.

Data is assumed to be oriented such that:

- **Rows** = "observations"
- **Columns** = "variables" / "features"

If your data is ordered in a different manner, you can simply transpose it before
passing into an `eda` class constructor.

For domain-specific cases, however, this may differ. For instance, biological
expression data is usually stored with the "features" (usually genes or 
transcripts) as rows, and the "observations" (arrays or RNA-Seq samples) as
columns. As such, `BioExprSet` expects data to be passed in in this order.

For example:

```r
library('eda')

edm <- EDADataFrame$new(mtcars, color=cyl, shape=am)
```


Installation
------------

The easiest way to install `eda` is using the `install_github()` function from
devtools:

```r
devtools::install_github('khughitt/eda')
```

Usage Examples
--------------

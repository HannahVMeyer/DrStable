
<!-- INDEX.md is generated from INDEX.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/meyer-lab-cshl/drStable.svg?branch=master)](https://travis-ci.org/meyer-lab-cshl/drStable)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## <i class="fas fa-shapes" aria-hidden="true"></i> drStable

The analysis of high-dimensional datasets often requires extracting
meaningful variables from the data or compress the data into a more
tractable number of features. Feature extraction can be achieved by
using dimensionality methods.

**drStable** provides an interface for easy access to 12 different
dimensionality reduction methods: DiffusionMap, DRR, ICA, LLE, Isomap,
LaplacianEigenmap, MDS, PCA, kPCA, nMDS, tSNE and UMAP.

**drStable** introduces a novel stability criterion which provides the
user with a tool of estimating the number of low-dimensional features
that can be reliably recovered and thus determine the stable dimensions
of the low-dimensional space.

**drStable** will be useful in any application that requires the
selection of a low-dimensional feature representation of the original
data. For instance, in the filed of genetic association studies, these
features can be used as the response variable to find genetic variants
associated with the low-dimensional phenotype representation. Stable
features for association testing harbor two advantages: they increase
the power for genetic discovery and guarantee reliable association
results.

## <i class="fa fa-rocket" aria-hidden="true"></i> Installation

Currently, **drStable** is available in alpha-version via:

``` r
library(devtools)
install_github("meyer-lab-cshl/drStable")
```

We are working on the release of a stable version of the package,
including improving test coverage, extension of the documentation and
use-case vignettes.

A log of version changes can be found
[here](https://github.com/meyer-lab-cshl/drStable/blob/master/NEWS.md).

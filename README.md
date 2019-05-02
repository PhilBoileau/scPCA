
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`scPCA`

[![Travis-CI Build
Status](https://travis-ci.org/PhilBoileau/scPCA.svg?branch=master)](https://travis-ci.org/PhilBoileau/scPCA)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/PhilBoileau/scPCA?branch=master&svg=true)](https://ci.appveyor.com/project/PhilBoileau/scPCA/)
[![Coverage
Status](https://img.shields.io/codecov/c/github/PhilBoileau/scPCA/master.svg)](https://codecov.io/github/PhilBoileau/scPCA?branch=master)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/scPCA.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/scPCA)
[![Bioc
Time](http://bioconductor.org/shields/years-in-bioc/scPCA.svg)](https://bioconductor.org/packages/release/bioc/html/scPCA.html)
[![Bioc
Downloads](http://bioconductor.org/shields/downloads/scPCA.svg)](https://bioconductor.org/packages/release/bioc/html/scPCA.html)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Sparse and contrastive principal components analysis for computational
> biology

**Author:** [Philippe Boileau](https://github.com/PhilBoileau)

-----

## What’s `scPCA`?

The `scPCA` R package facilitates Abid et al. (2018)

-----

## Installation

For standard use, install from
[Bioconductor](https://bioconductor.org/packages/scPCA) using
[`BiocManager`](https://CRAN.R-project.org/package=BiocManager):

``` r
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("scPCA")
```

To contribute, install the bleeding-edge *development version* from
GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/):

``` r
devtools::install_github("PhilBoileau/scPCA")
```

Current and prior [Bioconductor](https://bioconductor.org) releases are
available under branches with numbers prefixed by “RELEASE\_”. For
example, to install the version of this package available via
Bioconductor 4.0, use

``` r
devtools::install_github("PhilBoileau/scPCA", ref = "RELEASE_4_0")
```

-----

## Example

For details on how to best use the `scPCA` R package, please consult the
most recent [package
vignette](https://bioconductor.org/packages/release/bioc/vignettes/scPCA/inst/doc/scpca_intro.html)
available through the [Bioconductor
project](https://bioconductor.org/packages/scPCA).

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/PhilBoileau/scPCA/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/PhilBoileau/scPCA/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## License

© 2019 [Philippe Boileau](https://github.com/PhilBoileau)

The contents of this repository are distributed under the MIT license.
See file `LICENSE` for details.

-----

## References

<div id="refs" class="references">

<div id="ref-abid2018exploring">

Abid, Abubakar, Martin J Zhang, Vivek K Bagaria, and James Zou. 2018.
“Exploring Patterns Enriched in a Dataset with Contrastive Principal
Component Analysis.” *Nature Communications* 9 (1). Nature Publishing
Group: 2134.

</div>

</div>

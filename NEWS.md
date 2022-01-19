# Changes in version 0.99.0 (2019-09-13)

+ Submitted to Bioconductor

# Changes in version 1.1.2 (2020-01-08)

+ Added hierarchical clustering options for clustering based cross-validation

# Changes in version 1.1.5 (2020-01-18)

+ Fixed issue with matrix normalization
+ Misc. bug fixes
+ Improvements to code coverage

# Changes in version 1.1.11 (2020-02-02)

+ Added more SPCA algorithm options
  - SPCA via variable projection
  - Randomized SPCA via variable projection
+ New vignette section comparing performance of SPCA algorithms
+ Improvements to code coverage

# Changes in version 1.1.12 (2020-04-21)

+ Updated citations
+ Fixed typos in documentation

# Changes in version 1.1.14 (2020-04-26)

+ Fixing broken link in an internal function documentation page.

# Changes in version 1.1.15 (2020-06-02)

+ Fixing Travis CI settings

# Changes in version 1.3.2 (2020-08-05)

+ Updated `scPCA()` function documentation
+ Corrected spelling mistakes

# Changes in version 1.3.3 (2020-08-08)

+ The `n_centers` argument no longer matters when When the contrasts argument is of length 1 and the penalty term is set to 0.
+ Users can now pass in their own cluster labels

# Changes in version 1.3.4 (2020-08-12)

+ Replaced calls to `base::eigen()` by `RSpectra::eigs_sym()` to speed up eigendecompositions of contrastive covariance matrices. cPCA is now performed much more quickly when only wishing to compute a handful of leading contrastive principal components.
+ Replaced calls to `stats::cov()` by `coop::covar()` to speed up computation of large sample covariance matrices.
+ In future updates, we'd like to explore using the `DelayedArray` framework to support the analysis of larger datasets.

# Changes in version 1.3.5 (2020-08-18)

+ Fixed citations in docs
+ Provided more detailed warning when `RSpectra::eigs_sym()` fails to converge
+ Included arguments in `scPCA()` to control `RSpectra::eigs_sym()` convergence: error tolerance and max number of iterations

# Changes in version 1.3.6 (2020-08-30)

+ Fixed issue where `n_centers` was required when only one penalty and contrast term were provided
+ Users can now pass factors and character vectors to the clusters argument.

# Changes in version 1.3.7 and 1.3.8 (2020-09-01)

+ Minor bug fixes

# Changes in version 1.3.9 (2020-10-12)

+ `scPCA()` now accepts `DelayedMatrix` objects as target and background datasets.

# Changes in version 1.3.10 (2020-10-16)

+ Implementing suggested improvements from Aaron Lun.

# Changes in version 1.5.1 (2020-12-17)

+ `scPCA()` and other internal functions may now take advantage of the
  `ScaledMatrix` object class. This allows more computationally efficient
  contrastive covariance matrix estimation when analyzing large datasets.
+ `safeColScale()` now used `MatrixGenerics` to handle feature standardization.

# Changes in version 1.5.2 (2020-12-21)

+ Adding `LTLA/ScaledMatrix` to "Remotes" section of `DESCRIPTION`.

# Changes in version 1.5.3 (2021-03-14)

+ Updating plotting issue in vignette: comparison of cPCA and scPCA loadings.
+ Adding `pkgdown` site.
+ Moving `ScaledMatrix` to "imports" section of `DESCRIPTION`.

# Changes in version 1.7.2 (2021-09-15)

+ Removing tests checking that sequential and parallel calls to `scPCA()`
  produce identical outputs when `BiocParallel`'s `SerialParam()` is used. This
  due to new handing of random number generation in `BiocParallel` version 1.28.
  
# Changes in version 1.7.3 (2021-10-07)

+ Removing more tests attempting to verify that parallelized outputs perfectly
  match their serial counterparts.

# Changes in version 1.9.1 (2022-01-19)

+ Updating copyright years
+ Updating citation information

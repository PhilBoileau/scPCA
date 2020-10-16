# Changes in version 0.99.0 (2019-09-13)
+ Submitted to Bioconductor

# Changes in version 1.1.2 (2020-01-08)
+ Added hierarchical clustering options for clustering based cross-validation

# Changes in version 1.1.5 (2020-01-18)
+ Fixed issue with matrix normalization
+ Misc. big fixes
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
+ Fixing travis CI settings

# Changes in version 1.3.2 (2020-08-05)
+ Updated scPCA function documentation
+ Corrected pelling mistakes

# Changes in version 1.3.3 (2020-08-08)
+ The n_centers argument no longer matters when When the contrasts argument is of length 1 and the penalty term is set to 0.
+ Users can now pass in their own cluster labels

# Changes in version 1.3.4 (2020-08-12)
+ Replaced calls to base::eigen by RSpectra::eigs_sym to speed up eigendecompositions of contrastive covariance matrices. cPCA is now performed much more quickly when only whishing to compute a handful of leading contrastive principal components.
+ Replaced calls to stats::cov by coop::covar to speed up computation of large sample covariance matrices.
+ In future updates, we'd like to explore using the DelayedArray framework to support the analysis of larger datasets.

# Changes in version 1.3.5 (2020-08-18)
+ Fixed citations in docs
+ Provided more detailed warning when RSpectra::eigs_sym fails to converge
+ Included argumetns in scPCA to control RSpectra::eigs_sym convergence: error tolerance and max number of iterations

# Changes in version 1.3.6 (2020-08-30)
+ Fixed issue where n_centers was required when only one penalty and contrast term were provided
+ Users can now pass factors and character vectors to the clusters argument.

# Changes in version 1.3.7 and 1.3.8 (2020-09-01)
+ Minor bug fixes

# Changes in version 1.3.9 (2020-10-12)
+ scPCA now accepts DelayedMatrix objects as target and background datasets.

# Changes in version 1.3.10 (2020-10-16)
+ Implementing suggested improvements from Aaron Lun
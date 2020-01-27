---
title: "`scPCA`: A toolbox for (sparse) contrastive principal component analysis in `R`"
tags:
  - R
  - dimensionality reduction
  - computational biology
  - unwanted variation
  - sparsity
authors:
  - name: Philippe Boileau
    orcid: 0000-0002-4850-2507
    affiliation: 1
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1, 2
  - name: Sandrine Dudoit
    orcid: 0000-0002-6069-8629
    affiliation: 2, 3, 4
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Center for Computational Biology, University of California, Berkeley
    index: 2
  - name: Department of Statistics, University of California, Berkeley
    index: 3
  - name: Division of Epidemiology and Biostatistics, School of Public Health, University of California, Berkeley
    index: 4
date: 27 January 2020
bibliography: paper.bib
---

# Summary

Data pre-processing and exploratory data analysis are crucial steps in the data
science life-cycle, often relying on dimensionality reduction techniques to
extract pertinent signal. As the collection of large datasets becomes the norm,
the need for methods that can successfully glean relevant information, from
among increasingly sophisticated technical artifacts, is greater than in years
past. What's more, as collected datasets increase in complexity, many of the
most reliable classical methods have been shown to provided degraded
performance, or even fail outright, in reducing the dimensionality of large and
noisy datasets.

Principal component analysis (PCA) is one such method. Although popular for its
interpretable results and ease of implementation, PCA’s performance on
high-dimensional data often leaves much to be desired. Its performance has been
characterized as unstable on large datasets [@Johnstone2009], and it has been
shown to often erroneously identify technical variation in the presence of
noise.

Consequently, modifications of PCA have been developed to remedy these issues.
Namely, sparse PCA (sPCA) [@Zou2006] was created to increase the stability of
the principal component loadings and variable scores in high dimensions, while
constrastive PCA (cPCA) [@Abid2018] has been proposed to leverage control data
to capture biologically relevant information in high-dimensional settings.

The `scPCA` `R` package implements sparse constrastive PCA (scPCA) [@Boileau], a
combination of these methods, drawing on cPCA to remove technical effects and on
sPCA to sparsify the principal component loadings. In both simulation studies
and data analysis, @Boileau provided practical demonstrations of scPCA's ability
to extract stable, interpretable, and uncontaminated signal from
high-dimensional biological data. Such demonstrations included the re-analysis
of several publicly available protein expression, microarray gene expression,
and single-cell transcriptome sequencing datasets [@Boileau].

As the `scPCA` software package was specially designed for use in disentangling
biological signal from technical noise in high-throughput sequencing data,
a free and open-source software implementation has been made available via the
Bioconductor Project [@gentleman2004bioconductor; @gentleman2006bioinformatics;
@huber2015orchestrating] for the `R` language and environment for statistical
computing [@R]. The `scPCA` package implements both cPCA, previously unavailable
in the `R` language, in two flavors: (1) the semi-automated version of @Abid2018
and (2) the automated version formulated by @Boileau. In order to interface
seamlessly with data structures common in computational biology, the `scPCA`
package integrates fully with the `SingleCellExperiment` container class
[@lun2018singlecellexperiment], utilizing the class to store the cPCA and scPCA
representations generated via the `reducedDims` accessor method. Finally, to
facilitate parallel computation, the `scPCA` package contains parallelized
versions of each of its core subroutines, making use of the infrastructure
provided by the [`BiocParallel`
package](https://bioconductor.org/packages/BiocParallel) package. In order to
effectively utilize paralleization, one need only set `parallel = TRUE` in
a call to the `scPCA` package, after having registered a particular
parallelization backend, as per the `BiocParallel` documentation.

# References


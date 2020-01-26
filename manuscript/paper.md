---
title: 'scPCA: A toolbox for (sparse) contrastive principal component analysis'
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
	- name: Nima Hejazi
		orcid: 0000-0002-7127-2789
		affiliation: 1, 2
	- name: Sandrine Dudoit
		orcid: 0000-0002-6069-862
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
date: 22 January 2020
bibliography: paper.bib
---

# Summary

Data pre-processing and exploratory data analysis are important steps in the data science life-cycle. These steps often rely on dimensionality reduction techniques to isolate pertinent patterns in data. As data sets become larger and the signal they possess weaker, these methods' ability to glean relevant information become increasingly important. However, many of the most commonly-used methods fail to reduce the dimensions of these large and noisy data sets successfully.

Principal Component Analysis (PCA) is one such method. Although popular for its interpretable results and ease of implementation, PCA’s performance on high-dimensional data often leaves much to be desired. Its results on these large datasets have been found unstable [@Johnstone2009], and it often identifies variation that is inconsequential.

Modifications of PCA have been developed to remedy these issues. Namely, sparse PCA (sPCA) [@Zou2006] was created to increase the stability of the principal component loadings and variable scores in high dimensions, and constrastive PCA (cPCA) [@Abid2018] was proposed as a method for capturing relevant information in the high-dimensional data with the help of control data.

Although sPCA and cPCA have proven useful in resolving individual shortcomings of PCA, neither is capable of tackling the issues of stability and relevance simultaneously. The `scPCA` package implements a combination of these methods, dubbed sparse constrastive PCA (scPCA) [@Boileau], which draws on cPCA to remove technical effects and on SPCA for sparsification of the loadings, thereby extracting stable, interpretable, and uncontaminated signal from high-dimensional data.

`scPCA` was specially designed for researchers wishing to disentangle biological signal from technical noise in their high-throughput sequencing data by harnessing the information in their control samples. Its utility has already been demonstrated on several publicly available protein expression, microarray gene expression, and single cell transcriptome sequencing datasets [@Boileau]. A free and open-source software implementation is made available via the Bioconductor Project [@gentleman2004bioconductor;gentleman2006bioinformatics;huber2015orchestrating] and GitHub. cPCA, previously unavailable to `R` users, is also implemented in two flavors: (1) the semi-automated version of @Abid2018, and (2) the automated version presented in @Boileau.

# References

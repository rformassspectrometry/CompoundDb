[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/EuracBiomedicalResearch/CompoundDb/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/EuracBiomedicalResearch/CompoundDb/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](https://codecov.io/github/EuracBiomedicalResearch/CompoundDb/coverage.svg?branch=master)](https://codecov.io/github/EuracBiomedicalResearch/CompoundDb?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

# Installation and requirements

The package can be installed with

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("EuracBiomedicalResearch/CompoundDb")
```

In case there are missing package dependencies, run the code below first:

```r

#' Install also Jan Stanstrup's commonMZ package
devtools::install_github("stanstrup/commonMZ")
```


# Creating and using (chemical) compound databases

This package provides functionality to create and use compound databases
generated from (mostly publicly) available resources such as
[HMDB](http://www.hmdb.ca), [ChEBI](https://www.ebi.ac.uk/chebi/) and [PubChem](https://pubchem.ncbi.nlm.nih.gov).

# How to contribute

Contributions are welcome, but should follow certain guidelines:
1) Open an issue.
2) Create a new branch (internal collaborator) or fork the repository (external
contributor).
3) Add your code (following loosly Bioconductor's [coding
style](http://bioconductor.org/developers/how-to/coding-style/).
4) Ensure the package passes `R CMD build` and `R CMD check`.
5) Make a pull request.

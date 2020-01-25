[![Build Status](https://travis-ci.org/EuracBiomedicalResearch/CompoundDb.svg?branch=master)](https://travis-ci.org/EuracBiomedicalResearch/CompoundDb)
[![codecov.io](https://codecov.io/github/EuracBiomedicalResearch/CompoundDb/coverage.svg?branch=master)](https://codecov.io/github/EuracBiomedicalResearch/CompoundDb?branch=master)

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

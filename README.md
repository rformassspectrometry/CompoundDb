[![Build Status](https://travis-ci.org/EuracBiomedicalResearch/CompoundDb.svg?branch=master)](https://travis-ci.org/EuracBiomedicalResearch/CompoundDb)
[![codecov.io](https://codecov.io/github/EuracBiomedicalResearch/CompoundDb/coverage.svg?branch=master)](https://codecov.io/github/EuracBiomedicalResearch/CompoundDb?branch=master)

# Requirements

The package can be installed with

```r
devtools::install_github("EuracBiomedicalResearch/CompoundDb")
```

In case there are missing package dependencies, run the code below first:

```r
install.packages("BiocManager")
BiocManager::install(c("devtools", "AnnotationFilter", "S4Vectors",
                       "MSnbase", "ChemmineR", "tibble", "jsonlite",
                       "dbplyr", "RSQLite", "Biobase", "BiocGenerics",
                       "ProtGenerics", "xml2"))

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


# project

to create general purpose and flexible compound databases (german: *eierlegende
Wollmilchsau*) facilitating compound annotation in untargeted metabolics data
processing. Databases should be generated from a variety of input sources (HMDB,
MassBank etc) and providing thus also different extent of annotations (only
compound information (LipidBlast), compound information and MS2 spectrum data
(HMDB, MassBank)). Databases should be self contained, offline SQLite databases
with the possibility to store also larger data to MySQL database
servers. Standalone SQLite databases could be distributed over AnnotationHub. In
addition to the functionality to generate databases and provide access to the
data, the package should also provide functionality to annotate features
e.g. based on m/z by selecting the adducts. Spectrum data should be exported as
MS2 `Spectra` objects hence allowing spectrum data combination, definition of
consensus spectra and spectrum similarity calculations (as provided by
`MSnbase`).

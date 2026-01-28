# Usage of Annotation Resources with the CompoundDb Package

**Authors**: Jan Stanstrup \[aut\] (ORCID:
<https://orcid.org/0000-0003-0541-7369>), Johannes Rainer \[aut, cre\]
(ORCID: <https://orcid.org/0000-0002-6977-7147>), Josep M. Badia \[ctb\]
(ORCID: <https://orcid.org/0000-0002-5704-1124>), Roger Gine \[aut\]
(ORCID: <https://orcid.org/0000-0003-0288-9619>), Andrea Vicini \[aut\]
(ORCID: <https://orcid.org/0000-0001-9438-6909>), Prateek Arora \[ctb\]
(ORCID: <https://orcid.org/0000-0003-0822-9240>)  
**Last modified:** 2026-01-28 07:26:44.279722  
**Compiled**: Wed Jan 28 07:31:21 2026

## Introduction

The *[CompoundDb](https://bioconductor.org/packages/3.23/CompoundDb)*
package provides the functionality to create chemical *compound*
databases from a variety of sources and to use such annotation databases
(`CompDb`) (Rainer et al. 2022). A detailed description on the creation
of annotation resources is given in the *Creating CompoundDb annotation
resources* vignette. This vignette focuses on how annotations can be
search for and retrieved.

## Installation

The package (including dependencies) can be installed with the code
below:

``` r

install.packages("BiocManager")
BiocManager::install("CompoundDb")
```

## General usage

In this vignette we use a small `CompDb` database containing annotations
for a small number of metabolites build using
[MassBank](https://massbank.eu/MassBank/) release *2020.09*. The
respective `CompDb` database which is loaded below contains in addition
to general compound annotations also MS/MS spectra for these compounds.

``` r

library(CompoundDb)
cdb <- CompDb(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"))
cdb
```

    ## class: CompDb 
    ##  data source: MassBank 
    ##  version: 2020.09 
    ##  organism: NA 
    ##  compound count: 70 
    ##  MS/MS spectra count: 70

General information about the database can be accessed with the
`metadata` function.

``` r

metadata(cdb)
```

    ##                 name                         value
    ## 1             source                      MassBank
    ## 2                url https://massbank.eu/MassBank/
    ## 3     source_version                       2020.09
    ## 4        source_date                    1603272565
    ## 5           organism                          <NA>
    ## 6   db_creation_date      Thu Oct 22 08:45:31 2020
    ## 7 supporting_package                    CompoundDb
    ## 8  supporting_object                        CompDb

### Querying compound annotations

The `CompoundDb` package is designed to provide annotation resources for
small molecules, such as metabolites, that are characterized by an exact
mass and additional information such as their IUPAC International
Chemical Identifier
[InChI](https://en.wikipedia.org/wiki/International_Chemical_Identifier)
or their chemical formula. The available annotations (*variables*) for
compounds can differ between databases. The
[`compoundVariables()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function can be used to retrieve a list of all available compound
annotations for a specific `CompDb` database.

``` r

compoundVariables(cdb)
```

    ## [1] "formula"   "exactmass" "smiles"    "inchi"     "inchikey"  "cas"      
    ## [7] "pubchem"   "name"

The actual compound annotations can then be extracted with the
[`compounds()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function which returns by default all columns listed by
[`compoundVariables()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md).
We can also define specific columns we want to extract with the
`columns` parameter.

``` r

head(compounds(cdb, columns = c("name", "formula", "exactmass")))
```

    ##     formula exactmass         name
    ## 1  C10H10O3  178.0630      Mellein
    ## 2 C25H47NO9  505.3251 AAL toxin TB
    ## 3  C17H12O6  312.0634 Aflatoxin B1
    ## 4  C17H14O6  314.0790 Aflatoxin B2
    ## 5  C17H12O7  328.0583 Aflatoxin G1
    ## 6  C17H14O7  330.0739 Aflatoxin G2

As a technical detail, `CompDb` databases follow a very simple database
layout with only few constraints to allow data import and representation
for a variety of sources (e.g. MassBank, HMDB, MoNa, ChEBI). For the
present database, which is based on MassBank, the mapping between
entries in the *ms_compound* database table and MS/MS spectra is for
example 1:1 and the *ms_compound* table contains thus highly redundant
information. Thus, if we would include the column `"compound_id"` in the
query we would end up with redundant values:

``` r

head(compounds(cdb, columns = c("compound_id", "name", "formula")))
```

    ##   compound_id   formula         name
    ## 1           1  C10H10O3      Mellein
    ## 2           2  C10H10O3      Mellein
    ## 3           3  C10H10O3      Mellein
    ## 4           4  C10H10O3      Mellein
    ## 5           5  C10H10O3      Mellein
    ## 6           6 C25H47NO9 AAL toxin TB

By default,
[`compounds()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
extracts the data for **all** compounds stored in the database. The
function supports however also *filters* to get values for specific
entries only. These can be defined as *filter expressions* which are
similar to the way how e.g. a `data.frame` would be subsetted in R. In
the example below we extract the compound ID, name and chemical formula
for a compound *Mellein*.

``` r

compounds(cdb, columns = c("compound_id", "name", "formula"),
          filter = ~ name == "Mellein")
```

    ##   compound_id  formula    name
    ## 1           1 C10H10O3 Mellein
    ## 2           2 C10H10O3 Mellein
    ## 3           3 C10H10O3 Mellein
    ## 4           4 C10H10O3 Mellein
    ## 5           5 C10H10O3 Mellein

Note that a filter expression always has to start with `~` followed by
the *variable* on which the data should be subsetted and the condition
to select the entries of interest. An overview of available filters for
a `CompDb` can be retrieved with the `supportedFilter()` function which
returns the name of the filter and the database column on which the
filter selects the values:

``` r

supportedFilters(cdb)
```

    ##                 filter             field
    ## 1     CompoundIdFilter       compound_id
    ## 2      ExactmassFilter         exactmass
    ## 3        FormulaFilter           formula
    ## 4          InchiFilter             inchi
    ## 5       InchikeyFilter          inchikey
    ## 8 MsmsMzRangeMaxFilter msms_mz_range_max
    ## 7 MsmsMzRangeMinFilter msms_mz_range_min
    ## 6           NameFilter              name
    ## 9     SpectrumIdFilter       spectrum_id

Also, filters can be combined to create more specific filters in the
same manner this would be done in R, i.e. using `&` for *and*, `|` for
*or* and `!` for *not*. To illustrate this we extract below all compound
entries from the table for compounds with the name *Mellein* and that
have a `"compound_id"` which is either 1 or 5.

``` r

compounds(cdb, columns = c("compound_id", "name", "formula"),
          filter = ~ name == "Mellein" & compound_id %in% c(1, 5))
```

    ##   compound_id  formula    name
    ## 1           1 C10H10O3 Mellein
    ## 2           5 C10H10O3 Mellein

Similarly, we can define a filter expression to retrieve compounds with
an exact mass between 310 and 320.

``` r

compounds(cdb, columns = c("name", "exactmass"),
          filter = ~ exactmass > 310 & exactmass < 320)
```

    ##   exactmass         name
    ## 1  312.0634 Aflatoxin B1
    ## 2  314.0790 Aflatoxin B2

In addition to *filter expressions*, we can also define and combine
filters using the actual filter classes. This provides additional
conditions that would not be possible with regular filter expressions.
Below we fetch for examples only compounds from the database that
contain a *H14* in their formula. To this end we use a `FormulaFilter`
with the condition `"contains"`. Note that all filters that base on
character matching (i.e. `FormulaFilter`, `InchiFilter`,
`InchikeyFilter`, `NameFilter`) support as conditions also `"contains"`,
`"startsWith"` and `"endsWith"` in addition to `"="` and `"!="`.

``` r

compounds(cdb, columns = c("name", "formula", "exactmass"),
          filter = FormulaFilter("H14", "contains"))
```

    ##    formula exactmass         name
    ## 1 C17H14O6  314.0790 Aflatoxin B2
    ## 2 C17H14O7  330.0739 Aflatoxin G2

It is also possible to combine filters if they are defined that way,
even if it is a little less straight forward than with the filter
expressions. Below we combine the `FormulaFilter` with the
`ExactmassFilter` to retrieve only compounds with an `"H14"` in their
formula and an exact mass between 310 and 320.

``` r

filters <- AnnotationFilterList(
    FormulaFilter("H14", "contains"),
    ExactmassFilter(310, ">"),
    ExactmassFilter(320, "<"),
    logicOp = c("&", "&"))
compounds(cdb, columns = c("name", "formula", "exactmass"),
          filter = filters)
```

    ##    formula exactmass         name
    ## 1 C17H14O6   314.079 Aflatoxin B2

### Additional functionality for `CompDb` databases

*CompoundDb* defines additional functions to work with `CompDb`
databases. One of them is the
[`mass2mz()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function that allows to directly calculate ion (adduct) m/z values for
exact (monoisotopic) masses of compounds in a database. Below we use
this function to calculate `[M+H]+` and `[M+Na]+` ions for all unique
chemical formulas in our example `CompDb` database.

``` r

mass2mz(cdb, adduct = c("[M+H]+", "[M+Na]+"))
```

    ##              [M+H]+  [M+Na]+
    ## C10H10O3   179.0703 201.0522
    ## C25H47NO9  506.3324 528.3143
    ## C17H12O6   313.0706 335.0526
    ## C17H14O6   315.0863 337.0682
    ## C17H12O7   329.0656 351.0475
    ## C17H14O7   331.0812 353.0632
    ## C20H20N2O3 337.1547 359.1366
    ## C15H16O6   293.1020 315.0839
    ## C14H10O5   259.0601 281.0420
    ## C15H12O5   273.0757 295.0577
    ## C16H16O8   337.0918 359.0737

To get a `matrix` with adduct m/z values for discrete compounds
(identified by their InChIKey) we specify `name = "inchikey"`.

``` r

mass2mz(cdb, adduct = c("[M+H]+", "[M+Na]+"), name = "inchikey")
```

    ##                               [M+H]+  [M+Na]+
    ## KWILGNNWGSNMPA-UHFFFAOYSA-N 179.0703 201.0522
    ## CTXQVLLVFBNZKL-YVEDVMJTSA-N 506.3324 528.3143
    ## OQIQSTLJSLGHID-WNWIJWBNSA-N 313.0706 335.0526
    ## WWSYXEZEXMQWHT-WNWIJWBNSA-N 315.0863 337.0682
    ## XWIYFDMXXLINPU-WNWIJWBNSA-N 329.0656 351.0475
    ## WPCVRWVBBXIRMA-WNWIJWBNSA-N 331.0812 353.0632
    ## MJBWDEQAUQTVKK-IAGOWNOFSA-N 329.0656 351.0475
    ## SZINUGQCTHLQAZ-DQYPLSBCSA-N 337.1547 359.1366
    ## MMHTXEATDNFMMY-WBIUFABUSA-N 293.1020 315.0839
    ## CEBXXEKPIIDJHL-UHFFFAOYSA-N 259.0601 281.0420
    ## LCSDQFNUYFTXMT-UHFFFAOYSA-N 273.0757 295.0577
    ## VSMBLBOUQJNJIL-JJXSEGSLSA-N 337.0918 359.0737

Alternatively we could also use `name = "compound_id"` to get a value
for each row in the compound database table, but for this example
database this would result in highly redundant information.

[`mass2mz()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
bases on the
[`MetaboCoreUtils::mass2mz`](https://rdrr.io/pkg/MetaboCoreUtils/man/mass2mz.html)
function and thus supports all pre-defined adducts from that function.
These are (for positive polarity):

``` r

MetaboCoreUtils::adductNames()
```

    ##  [1] "[M+3H]3+"          "[M+2H+Na]3+"       "[M+H+Na2]3+"      
    ##  [4] "[M+Na3]3+"         "[M+2H]2+"          "[M+H+NH4]2+"      
    ##  [7] "[M+H+K]2+"         "[M+H+Na]2+"        "[M+C2H3N+2H]2+"   
    ## [10] "[M+2Na]2+"         "[M+C4H6N2+2H]2+"   "[M+C6H9N3+2H]2+"  
    ## [13] "[M+H]+"            "[M+Li]+"           "[M+2Li-H]+"       
    ## [16] "[M+NH4]+"          "[M+H2O+H]+"        "[M+Na]+"          
    ## [19] "[M+CH4O+H]+"       "[M+K]+"            "[M+C2H3N+H]+"     
    ## [22] "[M+2Na-H]+"        "[M+C3H8O+H]+"      "[M+C2H3N+Na]+"    
    ## [25] "[M+2K-H]+"         "[M+C2H6OS+H]+"     "[M+C4H6N2+H]+"    
    ## [28] "[2M+H]+"           "[2M+NH4]+"         "[2M+Na]+"         
    ## [31] "[2M+K]+"           "[2M+C2H3N+H]+"     "[2M+C2H3N+Na]+"   
    ## [34] "[3M+H]+"           "[M+H-NH3]+"        "[M+H-H2O]+"       
    ## [37] "[M+H-Hexose-H2O]+" "[M+H-H4O2]+"       "[M+H-CH2O2]+"     
    ## [40] "[M]+"

and for negative polarity:

``` r

MetaboCoreUtils::adductNames(polarity = "negative")
```

    ##  [1] "[M-3H]3-"      "[M-2H]2-"      "[M-H]-"        "[M+Na-2H]-"   
    ##  [5] "[M+Cl]-"       "[M+K-2H]-"     "[M+C2H3N-H]-"  "[M+CHO2]-"    
    ##  [9] "[M+C2H3O2]-"   "[M+Br]-"       "[M+C2F3O2]-"   "[2M-H]-"      
    ## [13] "[2M+CHO2]-"    "[2M+C2H3O2]-"  "[3M-H]-"       "[M-H+HCOONa]-"
    ## [17] "[M]-"

In addition, user-supplied adduct definitions are also supported (see
the help of
[`mass2mz()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
in the
*[MetaboCoreUtils](https://bioconductor.org/packages/3.23/MetaboCoreUtils)*
package for details).

### Accessing and using MS/MS data

`CompDb` database can also store and provide MS/MS spectral data. These
can be accessed *via* a `Spectra` object from the
*[Spectra](https://bioconductor.org/packages/3.23/Spectra)*
Bioconductor. Such a `Spectra` object for a `CompDb` can be created with
the [`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) function
as in the example below.

``` r

sps <- Spectra(cdb)
sps
```

    ## MSn data (Spectra) with 70 spectra in a MsBackendCompDb backend:
    ##       msLevel precursorMz  polarity
    ##     <integer>   <numeric> <integer>
    ## 1           2      179.07         1
    ## 2           2      179.07         1
    ## 3           2      179.07         1
    ## 4           2      179.07         1
    ## 5           2      179.07         1
    ## ...       ...         ...       ...
    ## 66          2     337.091         1
    ## 67          2     337.091         1
    ## 68          2     337.091         1
    ## 69          2     337.091         1
    ## 70          2     337.091         1
    ##  ... 47 more variables/columns.
    ##  Use  'spectraVariables' to list all of them.
    ##  data source: MassBank 
    ##  version: 2020.09 
    ##  organism: NA

This `Spectra` object uses a `MsBackendCompDb` to *represent* the MS
data of the `CompDb` database. In fact, only the compound identifiers
and the precursor m/z values from all spectra are stored in memory while
all other data is retrieved on-the-fly from the database when needed.

The
[`spectraVariables()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function lists all available annotations for a spectrum from the
database, which includes also annotations of the associated compounds.

``` r

spectraVariables(sps)
```

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "dataStorage"             "dataOrigin"             
    ##  [7] "centroided"              "smoothed"               
    ##  [9] "polarity"                "precScanNum"            
    ## [11] "precursorMz"             "precursorIntensity"     
    ## [13] "precursorCharge"         "collisionEnergy"        
    ## [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [17] "isolationWindowUpperMz"  "compound_id"            
    ## [19] "formula"                 "exactmass"              
    ## [21] "smiles"                  "inchi"                  
    ## [23] "inchikey"                "cas"                    
    ## [25] "pubchem"                 "name"                   
    ## [27] "accession"               "spectrum_name"          
    ## [29] "date"                    "authors"                
    ## [31] "license"                 "copyright"              
    ## [33] "publication"             "splash"                 
    ## [35] "precursor_intensity"     "adduct"                 
    ## [37] "ionization"              "ionization_voltage"     
    ## [39] "fragmentation_mode"      "collisionEnergy_text"   
    ## [41] "instrument"              "instrument_type"        
    ## [43] "precursorMz_text"        "spectrum_id"            
    ## [45] "predicted"               "msms_mz_range_min"      
    ## [47] "msms_mz_range_max"       "synonym"

Individual variables can then be accessed with `$` and the variable
name:

``` r

head(sps$adduct)
```

    ## [1] "[M+H]+" "[M+H]+" "[M+H]+" "[M+H]+" "[M+H]+" "[M+H]+"

For more information on how to use `Spectra` objects in your analysis
have also a look at the package
[vignette](https://rformassspectrometry.github.io/Spectra/articles/Spectra.html)
or a [tutorial](https://jorainer.github.io/SpectraTutorials/) on how to
perform MS/MS spectra matching with `Spectra`.

Similar to the
[`compounds()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function, a call to
[`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) will give
access to **all** spectra in the database. Using the same filtering
framework it is however also possible to *extract* only specific spectra
from the database. Below we are for example accessing only the MS/MS
spectra of the compound *Mellein*. Using the `filter` in the
[`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) call can be
substantially faster than first initializing a `Spectra` with the full
data and then subsetting that to selected spectra.

``` r

mellein <- Spectra(cdb, filter = ~ name == "Mellein")
mellein
```

    ## MSn data (Spectra) with 5 spectra in a MsBackendCompDb backend:
    ##     msLevel precursorMz  polarity
    ##   <integer>   <numeric> <integer>
    ## 1         2      179.07         1
    ## 2         2      179.07         1
    ## 3         2      179.07         1
    ## 4         2      179.07         1
    ## 5         2      179.07         1
    ##  ... 47 more variables/columns.
    ##  Use  'spectraVariables' to list all of them.
    ##  data source: MassBank 
    ##  version: 2020.09 
    ##  organism: NA

Instead of all spectra we extracted now only a subset of 5 spectra from
the database.

As a simple toy example we perform next pairwise spectra comparison
between the 5 spectra from *Mellein* with all the MS/MS spectra in the
database.

``` r

library(Spectra)
cormat <- compareSpectra(mellein, sps, ppm = 40)
```

Note that the `MsBackendCompDb` does not support parallel processing,
thus, while
[`compareSpectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
would in general support parallel processing, it gets automatically be
disabled if a `Spectra` with a `MsBackendCompDb` is used.

``` r

cormat <- compareSpectra(mellein, sps, ppm = 40, BPPARAM = MulticoreParam(2))
```

## Ion databases

The `CompDb` database layout is designed to provide compound
annotations, but in mass spectrometry (MS) ions are measured. These ions
are generated e.g. by electro spray ionization (ESI) from the original
compounds in a sample. They are characterized by their specific
mass-to-charge ratio (m/z) which is measured by the MS instrument.
Eventually, also a retention time is available. Also, for the same
compound several different ions (adducts) can be formed and measured,
all with a different m/z. This type of data can be represented by an
`IonDb` database, which extends the `CompDb` and hence inherits all of
its properties but adds additional database tables to support also ion
annotations. Also, `IonDb` objects provide functionality to add new ion
annotations to an existing database. Thus, this type of database can be
used to build lab-internal annotation resources containing ions, m/z and
retention times for pure standards measured on a specific e.g. LC-MS
setup.

`CompDb` databases, such as the `cdb` from this example, are however by
default *read-only*, thus, we below create a new database connection,
copy the content of the `cdb` to that database and convert the `CompDb`
to an `IonDb`.

``` r

library(RSQLite)
## Create a temporary database
con <- dbConnect(SQLite(), tempfile())
## Create an IonDb copying the content of cdb to the new database
idb <- IonDb(con, cdb)
idb
```

    ## class: IonDb 
    ##  data source: MassBank 
    ##  version: 2020.09 
    ##  organism: NA 
    ##  compound count: 70 
    ##  MS/MS spectra count: 70 
    ##  ion count: 0

The `IonDb` defines an additional function `ions` that allows to
retrieve ion information from the database.

``` r

ions(idb)
```

    ## [1] compound_id ion_adduct  ion_mz      ion_rt     
    ## <0 rows> (or 0-length row.names)

The present database does not yet contain any ion information. Below we
define a data frame with ion annotations and add that to the database
with the
[`insertIon()`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
function. The column `"compound_id"` needs to contain the identifiers of
the compounds to which the ion should be related to. In the present
example we add 2 different ions for the compound with the ID 1
(*Mellein*). Note that the specified m/z values as well as the retention
times are completely arbitrary.

``` r

ion <- data.frame(compound_id = c(1, 1),
                  ion_adduct = c("[M+H]+", "[M+Na]+"),
                  ion_mz = c(123.34, 125.34),
                  ion_rt = c(196, 196))
idb <- insertIon(idb, ion)
```

These ions have now be added to the database.

``` r

ions(idb)
```

    ##   compound_id ion_adduct ion_mz ion_rt
    ## 1           1     [M+H]+ 123.34    196
    ## 2           1    [M+Na]+ 125.34    196

Ions can also be deleted from a database with the `deleteIon` function
(see the respective help page for more information).

Note that we can also retrieve compound annotation information for the
ions. Below we extract the associated compound name and its exact mass.

``` r

ions(idb, columns = c("ion_adduct", "name", "exactmass"))
```

    ##   ion_adduct    name exactmass
    ## 1     [M+H]+ Mellein   178.063
    ## 2    [M+Na]+ Mellein   178.063

## Session information

``` r

sessionInfo()
```

    ## R Under development (unstable) (2026-01-25 r89330)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] RSQLite_2.4.5           Spectra_1.21.1          BiocParallel_1.45.0    
    ## [4] CompoundDb_1.15.2       S4Vectors_0.49.0        BiocGenerics_0.57.0    
    ## [7] generics_0.1.4          AnnotationFilter_1.35.0 BiocStyle_2.39.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6           rjson_0.2.23           xfun_0.56             
    ##  [4] bslib_0.10.0           ggplot2_4.0.1          htmlwidgets_1.6.4     
    ##  [7] Biobase_2.71.0         vctrs_0.7.1            tools_4.6.0           
    ## [10] bitops_1.0-9           parallel_4.6.0         tibble_3.3.1          
    ## [13] blob_1.3.0             cluster_2.1.8.1        pkgconfig_2.0.3       
    ## [16] dbplyr_2.5.1           RColorBrewer_1.1-3     S7_0.2.1              
    ## [19] desc_1.4.3             lifecycle_1.0.5        compiler_4.6.0        
    ## [22] farver_2.1.2           textshaping_1.0.4      Seqinfo_1.1.0         
    ## [25] codetools_0.2-20       clue_0.3-66            htmltools_0.5.9       
    ## [28] sass_0.4.10            RCurl_1.98-1.17        yaml_2.3.12           
    ## [31] lazyeval_0.2.2         pkgdown_2.2.0.9000     pillar_1.11.1         
    ## [34] jquerylib_0.1.4        MASS_7.3-65            DT_0.34.0             
    ## [37] cachem_1.1.0           MetaboCoreUtils_1.19.1 tidyselect_1.2.1      
    ## [40] digest_0.6.39          stringi_1.8.7          dplyr_1.1.4           
    ## [43] bookdown_0.46          rsvg_2.7.0             fastmap_1.2.0         
    ## [46] grid_4.6.0             cli_3.6.5              magrittr_2.0.4        
    ## [49] base64enc_0.1-3        ChemmineR_3.63.0       scales_1.4.0          
    ## [52] bit64_4.6.0-1          rmarkdown_2.30         bit_4.6.0             
    ## [55] otel_0.2.0             gridExtra_2.3          ragg_1.5.0            
    ## [58] png_0.1-8              memoise_2.0.1          evaluate_1.0.5        
    ## [61] knitr_1.51             GenomicRanges_1.63.1   IRanges_2.45.0        
    ## [64] rlang_1.1.7            Rcpp_1.1.1             glue_1.8.0            
    ## [67] DBI_1.2.3              xml2_1.5.2             BiocManager_1.30.27   
    ## [70] jsonlite_2.0.0         R6_2.6.1               ProtGenerics_1.39.2   
    ## [73] systemfonts_1.3.1      fs_1.6.6               MsCoreUtils_1.23.2

## References

Rainer, Johannes, Andrea Vicini, Liesa Salzer, et al. 2022. “A Modular
and Expandable Ecosystem for Metabolomics Data Annotation in R.”
*Metabolites* 12 (2): 173. <https://doi.org/10.3390/metabo12020173>.

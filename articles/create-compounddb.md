# Creating CompoundDb annotation resources

**Authors**: Johannes Rainer  
**Modified**: 2026-01-28 07:26:44.279722  
**Compiled**: Wed Jan 28 07:31:30 2026

## Introduction

Chemical compound annotation and information can be retrieved from a
variety of sources including [HMDB](http://www.hmdb.ca),
[LipidMaps](http://www.lipidmaps.org) or
[ChEBI](https://www.ebi.ac.uk/chebi). The
*[CompoundDb](https://bioconductor.org/packages/3.23/CompoundDb)*
package provides functionality to extract data relevant for
(chromatographic) peak annotations in metabolomics/lipidomics
experiments from these sources and to store it into a common format
(i.e. an `CompDb` object/database). This vignette describes how such
`CompDb` objects can be created exemplified with package-internal test
files that represent data subsets from some annotation resources.

The R object to represent the compound annotation is the `CompDb`
object. Each object (respectively its database) is supposed to contain
and provide annotations from a single source (e.g. HMDB or LipidMaps)
but it is also possible to create cross-source databases too.

## Creating `CompDb` databases

`CompDb` databases can be created from existing data resources such as
the Human Metabolome Database (HMDB) by importing all of their data or
can alternatively be *build* by sequentially adding data and information
to the database. This section explains first how data can be loaded from
existing resources to create a `CompDb` database and in the last
subsection how an empty `CompDb` can be sequentially and manually filled
with annotation data @ref(sec:fill).

The *CompoundDb* package provides the
[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
and the
[`compound_tbl_lipidblast()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_lipidblast.md)
functions to extract relevant compound annotation from files in SDF
(structure-data file) format or in the json files from LipidBlast
(<http://mona.fiehnlab.ucdavis.edu/downloads>). *CompoundDb* allows to
process SDF files from:

- Human Metabolome Database (HMDB), <http://www.hmdb.ca>
- Chemical Entities of Biological Interest (ChEBI):
  <https://www.ebi.ac.uk/chebi>
- LIPID MAPS Structure Database (LMSD): <http://www.lipidmaps.org>
- PubChem: <https://pubchem.ncbi.nlm.nih.gov>
- MoNa (Massbank of North America): <http://mona.fiehnlab.ucdavis.edu>
  (for MoNa import see the next section).

Note however that it is also possible to define such a table manually
and use that to create the database. As simple example for this is
provided in the section *`CompDb` from custom input* @ref(sec:custom)
below or the help page of
[`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
for more details on that.

### `CompDb` from HMDB data

Below we use the
[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
to extract compound annotations from a SDF file representing a very
small subset of the HMDB database. To generate a database for the full
HMDB we would have to download the *structures.sdf* file containing all
metabolites and load that file instead.

``` r

library(CompoundDb)

## Locate the file
hmdb_file <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
## Extract the data
cmps <- compound_tbl_sdf(hmdb_file)
```

The function returns by default a (`data.frame`-equivalent) `tibble`
(from the *tidyverse*’s *tibble* package).

``` r

cmps
```

    ## # A tibble: 9 × 8
    ##   compound_id name              inchi inchikey formula exactmass synonyms smiles
    ##   <chr>       <chr>             <chr> <chr>    <chr>       <dbl> <named > <chr> 
    ## 1 HMDB0000001 1-Methylhistidine InCh… BRMWTNU… C7H11N…     169.  <chr>    "CN1C…
    ## 2 HMDB0000002 1,3-Diaminopropa… InCh… XFNJVJP… C3H10N2      74.1 <chr>    "NCCC…
    ## 3 HMDB0000005 2-Ketobutyric ac… InCh… TYEYBOS… C4H6O3      102.  <chr>    "CCC(…
    ## 4 HMDB0000008 2-Hydroxybutyric… InCh… AFENDNX… C4H8O3      104.  <chr>    "CCC(…
    ## 5 HMDB0000010 2-Methoxyestrone  InCh… WHEUWNK… C19H24…     300.  <chr>    "[H][…
    ## 6 HMDB0000011 (R)-3-Hydroxybut… InCh… WHBMMWS… C4H8O3      104.  <chr>    "C[C@…
    ## 7 HMDB0000012 Deoxyuridine      InCh… MXHRCPN… C9H12N…     228.  <chr>    "OC[C…
    ## 8 HMDB0004370 N-Methyltryptami… InCh… NCIKQJB… C11H14…     174.  <chr>    "CNCC…
    ## 9 HMDB0006719 5,6-trans-Vitami… InCh… QYSXJUF… C27H44O     384.  <chr>    "CC(C…

The `tibble` contains columns

- `compound_id`: the resource-specific ID of the compound. Can be an
  `integer` or a `character`.
- `name`: the name of the compound, mostly a generic or common name.
- `inchi`: the compound’s inchi.
- `inchikey`: the INCHI key.
- `formula`: the chemical formula of the compound.
- `exactmass`: the compounds (monoisotopic) mass.
- `synonyms`: a `list` of aliases/synonyms for the compound.
- `smiles`: the SMILES of the compound.

To create a simple compound database, we could pass this `tibble` along
with additional required metadata information to the
[`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
function. In the present example we want to add however also MS/MS
spectrum data to the database. We thus load below the MS/MS spectra for
some of the compounds from the respective xml files downloaded from
HMDB. To this end we pass the path to the folder in which the files are
located to the
[`msms_spectra_hmdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_hmdb.md)
function. The function identifies the xml files containing MS/MS spectra
based on their their file name and loads the respective spectrum data.
The folder can therefore also contain other files, but the xml files
from HMDB should not be renamed or the function will not recognice them.
Note also that at present only MS/MS spectrum xml files from HMDB are
supported (one xml file per spectrum); these could be downloaded from
HMDB with the *hmdb_all_spectra.zip* file.

``` r

## Locate the folder with the xml files
xml_path <- system.file("xml", package = "CompoundDb")
spctra <- msms_spectra_hmdb(xml_path)
```

Also here, spectra information can be manually provided by adhering to
the expected structure of the `data.frame` (see
[`?createCompDb`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
for details).

At last we have to create the metadata for the resource. The metadata
information for a `CompDb` resource is crucial as it defines the origin
of the annotations and its version. This information should thus be
carefully defined by the user. Below we use the
[`make_metadata()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
helper function to create a `data.frame` in the expected format. The
organism should be provided in the format e.g. `"Hsapiens"` for human or
`"Mmusculus"` for mouse, i.e. capital first letter followed by lower
case characters without whitespaces.

``` r

metad <- make_metadata(source = "HMDB", url = "http://www.hmdb.ca",
                       source_version = "4.0", source_date = "2017-09",
                       organism = "Hsapiens")
```

With all the required data ready we create the SQLite database for the
HMDB subset. With `path` we specify the path to the directory in which
we want to save the database. This defaults to the current working
directory, but for this example we save the database into a temporary
folder.

``` r

db_file <- createCompDb(cmps, metadata = metad, msms_spectra = spctra,
                        path = tempdir())
```

The variable `db_file` is now the file name of the SQLite database. We
can pass this file name to the
[`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function to get the `CompDb` objects acting as the interface to the
database.

``` r

cmpdb <- CompDb(db_file)
cmpdb
```

    ## class: CompDb 
    ##  data source: HMDB 
    ##  version: 4.0 
    ##  organism: Hsapiens 
    ##  compound count: 9 
    ##  MS/MS spectra count: 4

To extract all compounds from the database we can use the
[`compounds()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function. The parameter `columns` allows to choose the database columns
to return. Any columns for any of the database tables are supported. To
get an overview of available database tables and their columns, the
[`tables()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function can be used:

``` r

tables(cmpdb)
```

    ## $ms_compound
    ## [1] "compound_id" "name"        "inchi"       "inchikey"    "formula"    
    ## [6] "exactmass"   "smiles"     
    ## 
    ## $msms_spectrum
    ##  [1] "original_spectrum_id" "compound_id"          "polarity"            
    ##  [4] "collision_energy"     "predicted"            "splash"              
    ##  [7] "instrument_type"      "instrument"           "precursor_mz"        
    ## [10] "spectrum_id"          "msms_mz_range_min"    "msms_mz_range_max"   
    ## 
    ## $msms_spectrum_peak
    ## [1] "spectrum_id" "mz"          "intensity"   "peak_id"    
    ## 
    ## $synonym
    ## [1] "compound_id" "synonym"

Below we extract only selected columns from the *compounds* table.

``` r

compounds(cmpdb, columns = c("name", "formula", "exactmass"))
```

    ##                        name   formula exactmass
    ## 1 (R)-3-Hydroxybutyric acid    C4H8O3  104.0473
    ## 2        1,3-Diaminopropane   C3H10N2   74.0844
    ## 3         1-Methylhistidine C7H11N3O2  169.0851
    ## 4     2-Hydroxybutyric acid    C4H8O3  104.0473
    ## 5        2-Ketobutyric acid    C4H6O3  102.0317
    ## 6          2-Methoxyestrone  C19H24O3  300.1725
    ## 7      5,6-trans-Vitamin D3   C27H44O  384.3392
    ## 8              Deoxyuridine C9H12N2O5  228.0746
    ## 9        N-Methyltryptamine  C11H14N2  174.1157

Analogously we can use the
[`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) function to
extract spectrum data from the database. The function returns by default
a `Spectra` object from the
*[Spectra](https://bioconductor.org/packages/3.23/Spectra)* package with
all spectra metadata available as *spectra variables*.

``` r

library(Spectra)
sps <- Spectra(cmpdb)
sps
```

    ## MSn data (Spectra) with 4 spectra in a MsBackendCompDb backend:
    ##     msLevel precursorMz  polarity
    ##   <integer>   <numeric> <integer>
    ## 1        NA          NA         1
    ## 2        NA          NA         1
    ## 3        NA          NA         1
    ## 4        NA          NA         0
    ##  ... 32 more variables/columns.
    ##  Use  'spectraVariables' to list all of them.
    ##  data source: HMDB 
    ##  version: 4.0 
    ##  organism: Hsapiens

The available *spectra variables* for the `Spectra` object can be
retrieved with
[`spectraVariables()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):

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
    ## [19] "name"                    "inchi"                  
    ## [21] "inchikey"                "formula"                
    ## [23] "exactmass"               "smiles"                 
    ## [25] "original_spectrum_id"    "predicted"              
    ## [27] "splash"                  "instrument_type"        
    ## [29] "instrument"              "spectrum_id"            
    ## [31] "msms_mz_range_min"       "msms_mz_range_max"      
    ## [33] "synonym"

Individual spectra variables can be accessed with the `$` operator:

``` r

sps$collisionEnergy
```

    ## [1] 10 25 NA 20

And the actual m/z and intensity values with
[`mz()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html) and
[`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):

``` r

mz(sps)
```

    ## NumericList of length 4
    ## [[1]] 109.2 124.2 124.5 170.16 170.52
    ## [[2]] 83.1 96.12 97.14 109.14 124.08 125.1 170.16
    ## [[3]] 44.1 57.9 61.4 71.2 73.8 78.3 78.8 ... 142.9 144.1 157.6 158 175.2 193.2
    ## [[4]] 111.0815386 249.2587746 273.2587746 ... 367.3006394 383.3319396

``` r

## m/z of the 2nd spectrum
mz(sps)[[2]]
```

    ## [1]  83.10  96.12  97.14 109.14 124.08 125.10 170.16

Note that it is also possible to retrieve specific spectra, e.g. for a
provided compound, or add compound annotations to the `Spectra` object.
Below we use the filter expression `~ compound_id == "HMDB0000001"`to
get only MS/MS spectra for the specified compound. In addition we ask
for the `"name"` and `"inchikey"` of the compound.

``` r

sps <- Spectra(cmpdb, filter = ~ compound_id == "HMDB0000001",
               columns = c(tables(cmpdb)$msms_spectrum, "name",
                           "inchikey"))
sps
```

    ## MSn data (Spectra) with 2 spectra in a MsBackendCompDb backend:
    ##     msLevel precursorMz  polarity
    ##   <integer>   <numeric> <integer>
    ## 1        NA          NA         1
    ## 2        NA          NA         1
    ##  ... 32 more variables/columns.
    ##  Use  'spectraVariables' to list all of them.
    ##  data source: HMDB 
    ##  version: 4.0 
    ##  organism: Hsapiens

The available spectra variables:

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
    ## [19] "name"                    "inchi"                  
    ## [21] "inchikey"                "formula"                
    ## [23] "exactmass"               "smiles"                 
    ## [25] "original_spectrum_id"    "predicted"              
    ## [27] "splash"                  "instrument_type"        
    ## [29] "instrument"              "spectrum_id"            
    ## [31] "msms_mz_range_min"       "msms_mz_range_max"      
    ## [33] "synonym"

The compound’s name and INCHI key have thus also been added as spectra
variables:

``` r

sps$inchikey
```

    ## [1] "BRMWTNUJHUMWMS-LURJTMIESA-N" "BRMWTNUJHUMWMS-LURJTMIESA-N"

To share or archive the such created `CompDb` database, we can also
create a dedicated R package containing the annotation. To enable
reproducible research, each `CompDb` package should contain the version
of the originating data source in its file name (which is by default
extracted from the metadata of the resource). Below we create a `CompDb`
package from the generated database file. Required additional
information we have to provide to the function are the package
creator/maintainer and its version.

``` r

createCompDbPackage(
    db_file, version = "0.0.1", author = "J Rainer", path = tempdir(),
    maintainer = "Johannes Rainer <johannes.rainer@eurac.edu>")
```

    ## Creating package in /tmp/RtmpcR5ugD/CompDb.Hsapiens.HMDB.4.0

The function creates a folder (in our case in a temporary directory)
that can be build and installed with `R CMD build` and `R CMD INSTALL`.

Special care should also be put on the license of the package that can
be passed with the `license` parameter. The license of the package and
how and if the package can be distributed will depend also on the
license of the originating resource.

### `CompDb` from custom data

A `CompDb` database can also be created from custom, manually defined
annotations. To illustrate this we create below first a `data.frame`
with some arbitrary compound annotations. According to the
[`?createCompDb`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
help page, the data frame needs to have columns `"compound_id"`,
`"name"`, `"inchi"`, `"inchikey"`, `"formula"`, `"exactmass"`,
`"synonyms"`. All columns except `"compound_id"` can also contain
missing values. It is also possible to define additional columns. Below
we thus create a `data.frame` with some compound annotations as well as
additional columns. Note that all these annotations in this example are
for illustration purposes only and are by no means *real*. Also, we
don’t provide any information for columns `"inchi"`, `"inchikey"` and
`"formula"` setting all values for these to `NA`.

``` r

cmps <- data.frame(
    compound_id = c("CP_0001", "CP_0002", "CP_0003", "CP_0004"),
    name = c("A", "B", "C", "D"),
    inchi = NA_character_,
    inchikey = NA_character_,
    formula = NA_character_,
    exactmass = c(123.4, 234.5, 345.6, 456.7),
    compound_group = c("G01", "G02", "G01", "G03")
)
```

Next we add also *synonyms* for each compound. This columns supports
multiple values for each row.

``` r

cmps$synonyms <- list(
    c("a", "AA", "aaa"),
    c(),
    c("C", "c"),
    ("d")
)
```

We also need to define the *metadata* for our database, which we do with
the
[`make_metadata()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
function. With this information we can already create a first
rudimentary `CompDb` database that contains only compound annotations.
We thus create below our custom `CompDb` database in a temporary
directory. We also manually specify the name of our database with the
`dbFile` parameter - if not provided, the name of the database will be
constructed based on information from the `metadata` parameter. In a
real-case scenario, `path` and `dbFile` should be changed to something
more meaningful.

``` r

metad <- make_metadata(source = "manually defined", url = "",
                       source_version = "1.0.0", source_date = "2022-03-01",
                       organism = NA_character_)

db_file <- createCompDb(cmps, metadata = metad, path = tempdir(),
                        dbFile = "CompDb.test.sqlite")
```

We can now load this toy database using the
[`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function providing the full path to the database file. Note that we load
the database in read-write mode by specifying
`flags = RSQLite::SQLITE_RW` - by default
[`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
will load databases in read-only mode hence ensuring that the data
within the database can not be compromised. In our case we would however
like to add more information to this database later and hence we load it
in read-write mode.

``` r

cdb <- CompDb(db_file, flags = RSQLite::SQLITE_RW)
cdb
```

    ## class: CompDb 
    ##  data source: manually defined 
    ##  version: 1.0.0 
    ##  organism: NA 
    ##  compound count: 4

We can now retrieve annotations from the database with the `compound()`
function.

``` r

compounds(cdb)
```

    ##   name inchi inchikey formula exactmass compound_group
    ## 1    A  <NA>     <NA>    <NA>     123.4            G01
    ## 2    B  <NA>     <NA>    <NA>     234.5            G02
    ## 3    C  <NA>     <NA>    <NA>     345.6            G01
    ## 4    D  <NA>     <NA>    <NA>     456.7            G03

Or also search and filter the annotations.

``` r

compounds(cdb, filter = ~ name %in% c("B", "A"))
```

    ##   name inchi inchikey formula exactmass compound_group
    ## 1    A  <NA>     <NA>    <NA>     123.4            G01
    ## 2    B  <NA>     <NA>    <NA>     234.5            G02

Next we would like to add also MS2 spectra data to the database. This
could be either done directly in the
[`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
call with parameter `msms_spectra`, or with the
[`insertSpectra()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function that allows to add MS2 spectra data to an existing `CompDb`
which can be provided as a `Spectra` object. We thus below manually
create a `Spectra` object with some arbitrary MS2 spectra -
alternatively, `Spectra` can be imported from a variety of input
sources, including MGF or MSP files using e.g. the
*[MsBackendMgf](https://bioconductor.org/packages/3.23/MsBackendMgf)* or
*[MsBackendMsp](https://bioconductor.org/packages/3.23/MsBackendMsp)*
packages.

``` r

#' Define basic spectra variables
df <- DataFrame(msLevel = 2L, precursorMz = c(124.4, 124.4, 235.5))
#' Add m/z and intensity information for each spectrum
df$mz <- list(
    c(3, 20, 59.1),
    c(2, 10, 30, 59.1),
    c(100, 206, 321.1))
df$intensity <- list(
    c(10, 13, 45),
    c(5, 8, 9, 43),
    c(10, 20, 400))
#' Create the Spectra object
sps <- Spectra(df)
```

The `Spectra` object needs also to have a variable (column) called
`"compound_id"` which provides the information with which existing
compound in the database the spectrum is associated.

``` r

compounds(cdb, "compound_id")
```

    ##   compound_id
    ## 1     CP_0001
    ## 2     CP_0002
    ## 3     CP_0003
    ## 4     CP_0004

``` r

sps$compound_id <- c("CP_0001", "CP_0001", "CP_0002")
```

We can also add additional information to the spectra, such as the
instrument.

``` r

sps$instrument <- "AB Sciex TripleTOF 5600+"
```

And we can now add these spectra to our existing toy `CompDb`. Parameter
`columns` allows to specify which of the *spectra variables* should be
stored into the database.

``` r

cdb <- insertSpectra(cdb, spectra = sps,
                     columns = c("compound_id", "msLevel",
                                 "precursorMz", "instrument"))
cdb
```

    ## class: CompDb 
    ##  data source: manually defined 
    ##  version: 1.0.0 
    ##  organism: NA 
    ##  compound count: 4 
    ##  MS/MS spectra count: 3

We have thus now a `CompDb` database with compound annotations and 3 MS2
spectra. We could for example also retrieve the MS2 spectra for the
compound with the name *A* from the database with:

``` r

Spectra(cdb, filter = ~ name == "A")
```

    ## MSn data (Spectra) with 2 spectra in a MsBackendCompDb backend:
    ##     msLevel precursorMz  polarity
    ##   <integer>   <numeric> <integer>
    ## 1         2       124.4        NA
    ## 2         2       124.4        NA
    ##  ... 31 more variables/columns.
    ##  Use  'spectraVariables' to list all of them.
    ##  data source: manually defined 
    ##  version: 1.0.0 
    ##  organism: NA

### `CompDb` from MoNA data

MoNa (Massbank of North America) provides a large SDF file with all
spectra which can be used to create a `CompDb` object with compound
information and MS/MS spectra. Note however that MoNa is organized by
spectra and the annotation of the compounds is not consistent and
normalized. Spectra from the same compound can have their own compound
identified and data that e.g. can differ in their chemical formula,
precision of their exact mass or other fields.

Similar to the example above, compound annotations can be imported with
the
[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
function while spectrum data can be imported with
[`msms_spectra_mona()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_mona.md).
In the example below we use however the
[`import_mona_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/import_mona_sdf.md)
that wraps both functions to reads both compound and spectrum data from
a SDF file without having to import the file twice. As an example we use
a small subset from a MoNa SDF file that contains only 7 spectra.

``` r

mona_sub <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
                        package = "CompoundDb")
mona_data <- import_mona_sdf(mona_sub)
```

    ## Warning: MoNa data can currently not be normalized and the compound table
    ## contains thus highly redundant data.

As a result we get a `list` with a data.frame each for compound and
spectrum information. These can be passed along to the
[`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
function to create the database (see below).

``` r

metad <- make_metadata(source = "MoNa",
                       url = "http://mona.fiehnlab.ucdavis.edu/",
                       source_version = "2018.11", source_date = "2018-11",
                       organism = "Unspecified")
mona_db_file <- createCompDb(mona_data$compound, metadata = metad,
                             msms_spectra = mona_data$msms_spectrum,
                             path = tempdir())
```

We can now load and use this database, e.g. by extracting all compounds
as shown below.

``` r

mona <- CompDb(mona_db_file)
compounds(mona)
```

    ##                   name
    ## 1 Sulfachlorpyridazine
    ## 2         Sulfaclozine
    ## 3        Sulfadimidine
    ## 4       Sulfamethazine
    ## 5       Sulfamethazine
    ##                                                                                                      inchi
    ## 1             InChI=1S/C10H9ClN4O2S/c11-9-5-6-10(14-13-9)15-18(16,17)8-3-1-7(12)2-4-8/h1-6H,12H2,(H,14,15)
    ## 2             InChI=1S/C10H9ClN4O2S/c11-9-5-13-6-10(14-9)15-18(16,17)8-3-1-7(12)2-4-8/h1-6H,12H2,(H,14,15)
    ## 3 InChI=1S/C12H14N4O2S/c1-8-7-9(2)15-12(14-8)16-19(17,18)11-5-3-10(13)4-6-11/h3-7H,13H2,1-2H3,(H,14,15,16)
    ## 4 InChI=1S/C12H14N4O2S/c1-8-7-9(2)15-12(14-8)16-19(17,18)11-5-3-10(13)4-6-11/h3-7H,13H2,1-2H3,(H,14,15,16)
    ## 5 InChI=1S/C12H14N4O2S/c1-8-7-9(2)15-12(14-8)16-19(17,18)11-5-3-10(13)4-6-11/h3-7H,13H2,1-2H3,(H,14,15,16)
    ##                      inchikey      formula exactmass smiles
    ## 1 XOXHILFPRYWFOD-UHFFFAOYSA-N C10H9ClN4O2S  284.0135   <NA>
    ## 2 QKLPUVXBJHRFQZ-UHFFFAOYSA-N C10H9ClN4O2S  284.0135   <NA>
    ## 3 ASWVTGNCAZCNNR-UHFFFAOYSA-N  C12H14N4O2S  278.0837   <NA>
    ## 4 ASWVTGNCAZCNNR-UHFFFAOYSA-N  C12H14N4O2S  278.0837   <NA>
    ## 5 ASWVTGNCAZCNNR-UHFFFAOYSA-N  C12H14N4O2S  278.0837   <NA>

As stated in the introduction of this section the `compound` information
contains redundant information and the table has essentially one row per
spectrum. Feedback on how to reduce the redundancy in the ms_compound
table is highly appreciated.

### `CompDb` by sequentially filling with data

As an alternative to creating a full database from an existing resource
it is also possible to create an empty `CompDb` database and
*sequentially* filling it with data. This could for example be used to
create a laboratory specific annotation library with compound, ion and
fragment spectra of pure standards measured on a certain LC-MS setup.
Below we create an empty `CompDb` database providing the file name of
the database. In the example we store the database to a temporary file
but in a real use case a meaningful file name and file path should be
used instead.

``` r

dbfile <- tempfile()
mydb <- emptyCompDb(dbfile)
mydb
```

    ## class: CompDb 
    ##  data source: NA 
    ##  version: NA 
    ##  organism: NA 
    ##  compound count: 0

We next define some first compound annotation we want to add to the
database. For compound annotations, fields `"compound_id"` (an arbitrary
ID of the compound), `"name"` (the compound name), `"inchi"`,
`"inchikey"`, `"formula"` (the chemical formula) and `"exactmass"` (the
monoisotopic mass) are expected, but, except of `"compound_id"`, they
can also contain missing values or be completely omitted. Below we
define a `data.frame` with annotations for some compounds and add this
annotation to the database using the
[`insertCompound()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function.

``` r

cmp <- data.frame(compound_id = c("1", "2"),
                  name = c("Caffeine", "Glucose"),
                  formula = c("C8H10N4O2", "C6H12O6"),
                  exactmass = c(194.080375584, 180.063388116))
mydb <- insertCompound(mydb, cmp)
mydb
```

    ## class: CompDb 
    ##  data source: NA 
    ##  version: NA 
    ##  organism: NA 
    ##  compound count: 2

We next add fragment spectra for the compounds. These could for example
represent MS2 spectra measured for the pure standard of the compound and
could be extracted for example from
*[xcms](https://bioconductor.org/packages/3.23/xcms)* result objects or
other sources. Below we load some fragment spectra for caffeine from an
MGF file distributed with this package. We use the
*[MsBackendMgf](https://bioconductor.org/packages/3.23/MsBackendMgf)*
package to import that data into a `Spectra` object.

``` r

library(MsBackendMgf)
caf_ms2 <- Spectra(system.file("mgf", "caffeine.mgf", package = "CompoundDb"),
                   source = MsBackendMgf())
caf_ms2
```

    ## MSn data (Spectra) with 2 spectra in a MsBackendMgf backend:
    ##     msLevel     rtime scanIndex
    ##   <integer> <numeric> <integer>
    ## 1         2        NA        NA
    ## 2         2        NA        NA
    ##  ... 25 more variables/columns.

We can evaluate what spectra variables are available in the imported
data.

``` r

spectraVariables(caf_ms2)
```

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "dataStorage"             "dataOrigin"             
    ##  [7] "centroided"              "smoothed"               
    ##  [9] "polarity"                "precScanNum"            
    ## [11] "precursorMz"             "precursorIntensity"     
    ## [13] "precursorCharge"         "collisionEnergy"        
    ## [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [17] "isolationWindowUpperMz"  "TITLE"                  
    ## [19] "original_spectrum_id"    "compound_id"            
    ## [21] "collision_energy"        "predicted"              
    ## [23] "splash"                  "spectrum_id"            
    ## [25] "msms_mz_range_min"       "msms_mz_range_max"

``` r

caf_ms2$rtime
```

    ## [1] NA NA

There are many variables available, but most of them, like for example
the retention time, or not defined as this information was not provided
in the MGF file. In order to associate these fragment spectra to the
caffeine compound we just added to the database, we need to assign them
the ID of the compound (in our case `"1"`).

``` r

caf_ms2$compound_id <- "1"
```

We can then add the spectra to the database using the
[`insertSpectra()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
function. With parameter `columns` we specify which of the spectra
variables we actually want to store in the database.

``` r

mydb <- insertSpectra(mydb, caf_ms2,
                      columns = c("compound_id", "msLevel", "splash",
                                  "precursorMz", "collisionEnergy"))
mydb
```

    ## class: CompDb 
    ##  data source: NA 
    ##  version: NA 
    ##  organism: NA 
    ##  compound count: 2 
    ##  MS/MS spectra count: 2

We thus have now 2 compounds in the database and 2 fragment spectra:

``` r

compounds(mydb)
```

    ##       name inchi inchikey   formula exactmass
    ## 1 Caffeine  <NA>     <NA> C8H10N4O2  194.0804
    ## 2  Glucose  <NA>     <NA>   C6H12O6  180.0634

``` r

sps <- Spectra(mydb)
sps$name
```

    ## [1] "Caffeine" "Caffeine"

With `insertCommpound()` and
[`insertSpectra()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
further compounds and fragment spectra could be added to the database.
Note that both functions support also to add additional columns
(*fields* or *variables*) to the database. As an example we define below
a compound with an arbitrary additional column and add this to the
database using parameter `addColumns = TRUE`.

``` r

cmps <- data.frame(compound_id = "3", name = "X003",
                   formula = "C5H2P3O", extra_field = "artificial compound")
mydb <- insertCompound(mydb, cmps, addColumns = TRUE)
compounds(mydb)
```

    ##       name inchi inchikey   formula exactmass         extra_field
    ## 1 Caffeine  <NA>     <NA> C8H10N4O2  194.0804                <NA>
    ## 2  Glucose  <NA>     <NA>   C6H12O6  180.0634                <NA>
    ## 3     X003  <NA>     <NA>   C5H2P3O        NA artificial compound

The additional column is now available in the database. Existing entries
in a `CompDb` can also be deleted using the
[`deleteCompound()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
or
[`deleteSpectra()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
functions. Both require as additional input the IDs of the compound(s)
(or spectra) to delete. Below we extract the IDs and names of the
compounds from our database.

``` r

compounds(mydb, columns = c("compound_id", "name"))
```

    ##   compound_id     name
    ## 1           1 Caffeine
    ## 2           2  Glucose
    ## 3           3     X003

We can now delete the compound `"X003"` with
[`deleteCompound()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
and the ID of this compound.

``` r

mydb <- deleteCompound(mydb, ids = "3")
compounds(mydb)
```

    ##       name inchi inchikey   formula exactmass extra_field
    ## 1 Caffeine  <NA>     <NA> C8H10N4O2  194.0804        <NA>
    ## 2  Glucose  <NA>     <NA>   C6H12O6  180.0634        <NA>

Note that deleting a compound with associated spectra (or ions) will
result in an error, thus it would not be possible to delete caffeine
from the database, because it contains also MS2 spectra for that
compound. Using parameter `recursive = TRUE` in the
[`deleteCompound()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
call would however allow to delete the compound **and all** associated
spectra (and/or ions) along with it. Below we delete thus caffeine and
the associated MS2 spectra which leaves us a `CompDb` with a single
compound and no more MS2 spectra.

``` r

mydb <- deleteCompound(mydb, ids = "1", recursive = TRUE)
compounds(mydb)
```

    ##      name inchi inchikey formula exactmass extra_field
    ## 1 Glucose  <NA>     <NA> C6H12O6  180.0634        <NA>

``` r

Spectra(mydb)
```

    ## MSn data (Spectra) with 0 spectra in a MsBackendCompDb backend:

Note that these functions can also be used to add or remove annotations
to/from any `CompDb` database, as long as the database is *writeable*
(i.e. the database is loaded by specifying `flags = RSQLite::SQLITE_RW`
as additional parameter to the `CompDb` call to load the database).

## Extending `CompDb` databases

The `CompDb` database layout is very flexible supporting extension with
new database tables. The content of such additional database tables can
be automatically integrated into the information returned by the
[`compounds()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
or
[`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
functions if the relationships between the new and existing database
tables is defined. Elements of the new tables need to be linked
(related) to entries in at least one of the *core* `CompDb` database
tables. The respective relationship information can then be added to a
`CompDb` object using the
[`addJoinDefinition()`](https://rformassspectrometry.github.io/CompoundDb/reference/addJoinDefinition.md)
function.

As an example we add a new database table *experiment* to a `CompDb`
database. This new database table could contain information on the
individual experiments in which certain compounds stored in the
`CompDb`’s *ms_compound* table were measured.

We add this new table to the test database of the *CompoundDb*. Below we
first copy this database to a temporary file and connect through it
using the *RSQLite* R package.

``` r

#' Copy the test database to a temporary file
tf <- tempfile()
file.copy(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"), tf)
```

    ## [1] TRUE

``` r

#' Connect to this new database
library(RSQLite)
con <- dbConnect(SQLite(), tf)
```

The tables in this database are:

``` r

dbListTables(con)
```

    ## [1] "metadata"           "ms_compound"        "msms_spectrum"     
    ## [4] "msms_spectrum_peak" "synonym"

As detailed above, we want to add a new database table *experiment* with
some additional (experiment-related) information on the
experiments/measurement runs in which certain compounds were measured.
Such a table could for example be defined as follows:

``` r

exp <- data.frame(
    experiment_id = 1:4,
    experiment_name = c("exp_a", "exp_b", "exp_c", "exp_d"),
    experiment_date = c("2022-11-24", "2023-06-12", "2023-09-21", "2024-12-07"),
    experiment_operator = c("AA", "VV", "AA", "AA"))
```

We next need to *link* entries in this table with rows in the
*ms_compound* database table. We thus load the content of this table
from the database and add a column `"experiment_id"` to link it’s
entries to the new *experiment* table defined above.

``` r

cmp <- dbGetQuery(con, "select * from ms_compound")
nrow(cmp)
```

    ## [1] 70

We add the new column and fill it with values. Each row should contain
one of the values of the `"experiment_id"` column of `exp`. In our
example we fill the column randomly with values between 1 and 4 (the
primary keys of the elements in the `exp` data frame). In a real use
case this reference should of course represent a real relationship
between elements in the two tables.

``` r

cmp$experiment_id <- sample(1:4, nrow(cmp), replace = TRUE)
```

We next write both tables to the database, replacing the existing
*ms_compound* table and adding a new *experiment* table.

``` r

dbWriteTable(con, name = "ms_compound", cmp, overwrite = TRUE)
dbWriteTable(con, name = "experiment", exp)
dbListTables(con)
```

    ## [1] "experiment"         "metadata"           "ms_compound"       
    ## [4] "msms_spectrum"      "msms_spectrum_peak" "synonym"

``` r

dbDisconnect(con)
```

We can now load this database as a `CompDb` object.

``` r

cdb <- CompDb(tf)
```

    ## Note: found one or more database tables without defined relationship to other database tables. These are: "experiment". You might consider to use the 'addJoinDefinition()' to specify how that table(s) could be joined with the other database tables.

The message indicates that a database table was found which is not
related to any other database table. We could use the `cdb` variable as
any other `CompDb` object, the content of the new *experiment* would
however not (yet) be included in the results returned from this
database. The full list of available database tables (including their
columns) is:

``` r

tables(cdb)
```

    ## $experiment
    ## [1] "experiment_id"       "experiment_name"     "experiment_date"    
    ## [4] "experiment_operator"
    ## 
    ## $ms_compound
    ##  [1] "compound_id"   "formula"       "exactmass"     "smiles"       
    ##  [5] "inchi"         "inchikey"      "cas"           "pubchem"      
    ##  [9] "name"          "experiment_id"
    ## 
    ## $msms_spectrum
    ##  [1] "accession"             "spectrum_name"         "date"                 
    ##  [4] "authors"               "license"               "copyright"            
    ##  [7] "publication"           "ms_level"              "polarity"             
    ## [10] "splash"                "compound_id"           "precursor_intensity"  
    ## [13] "precursor_mz"          "adduct"                "ionization"           
    ## [16] "ionization_voltage"    "fragmentation_mode"    "collision_energy_text"
    ## [19] "instrument"            "instrument_type"       "precursor_mz_text"    
    ## [22] "spectrum_id"           "collision_energy"      "predicted"            
    ## [25] "msms_mz_range_min"     "msms_mz_range_max"    
    ## 
    ## $msms_spectrum_peak
    ## [1] "spectrum_id" "mz"          "intensity"   "peak_id"    
    ## 
    ## $synonym
    ## [1] "compound_id" "synonym"

In order to allow extraction of content of the new database table, we
need to add the relationship information between the new table
*experiment* and the other existing database tables in the `CompDb`
object. We can add this information using the
[`addJoinDefinition()`](https://rformassspectrometry.github.io/CompoundDb/reference/addJoinDefinition.md)
function. With parameters `table_a` and `table_b` we provide the names
of the database tables that are related to each other and with
`column_a` and `column_b` we specify the names of the columns in the two
tables that define the relationship (i.e. contain the primary and
foreign keys of the relationship).

``` r

cdb <- addJoinDefinition(
    cdb, table_a = "ms_compound", table_b = "experiment",
    column_a = "experiment_id", column_b = "experiment_id")
```

We can now also request information from the *experiment* database
table. Below we show the first 6 rows returned by joining the
*ms_compound* and *experiment* tables.

``` r

#' Show the first 6 rows
compounds(cdb, columns = c("compound_id", "formula", "experiment_name",
                           "experiment_date")) |>
    head()
```

    ##   compound_id   formula experiment_name experiment_date
    ## 1           1  C10H10O3           exp_d      2024-12-07
    ## 2           2  C10H10O3           exp_d      2024-12-07
    ## 3           3  C10H10O3           exp_a      2022-11-24
    ## 4           4  C10H10O3           exp_d      2024-12-07
    ## 5           5  C10H10O3           exp_b      2023-06-12
    ## 6           6 C25H47NO9           exp_c      2023-09-21

Information on the *experiment* table is also available through the
`Spectra` interface. Below we return the full content from the database
through the
[`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function.

``` r

s <- Spectra(cdb)
spectraData(s)
```

    ## DataFrame with 70 rows and 52 columns
    ##       msLevel     rtime acquisitionNum scanIndex dataStorage  dataOrigin
    ##     <integer> <numeric>      <integer> <integer> <character> <character>
    ## 1           2        NA             NA        NA        <db>          NA
    ## 2           2        NA             NA        NA        <db>          NA
    ## 3           2        NA             NA        NA        <db>          NA
    ## 4           2        NA             NA        NA        <db>          NA
    ## 5           2        NA             NA        NA        <db>          NA
    ## ...       ...       ...            ...       ...         ...         ...
    ## 66          2        NA             NA        NA        <db>          NA
    ## 67          2        NA             NA        NA        <db>          NA
    ## 68          2        NA             NA        NA        <db>          NA
    ## 69          2        NA             NA        NA        <db>          NA
    ## 70          2        NA             NA        NA        <db>          NA
    ##     centroided  smoothed  polarity precScanNum precursorMz precursorIntensity
    ##      <logical> <logical> <integer>   <integer>   <numeric>          <numeric>
    ## 1           NA        NA         1          NA      179.07                 NA
    ## 2           NA        NA         1          NA      179.07                 NA
    ## 3           NA        NA         1          NA      179.07                 NA
    ## 4           NA        NA         1          NA      179.07                 NA
    ## 5           NA        NA         1          NA      179.07                 NA
    ## ...        ...       ...       ...         ...         ...                ...
    ## 66          NA        NA         1          NA     337.091                 NA
    ## 67          NA        NA         1          NA     337.091                 NA
    ## 68          NA        NA         1          NA     337.091                 NA
    ## 69          NA        NA         1          NA     337.091                 NA
    ## 70          NA        NA         1          NA     337.091                 NA
    ##     precursorCharge collisionEnergy isolationWindowLowerMz
    ##           <integer>       <numeric>              <numeric>
    ## 1                NA              NA                     NA
    ## 2                NA              NA                     NA
    ## 3                NA              NA                     NA
    ## 4                NA              NA                     NA
    ## 5                NA              NA                     NA
    ## ...             ...             ...                    ...
    ## 66               NA              NA                     NA
    ## 67               NA              NA                     NA
    ## 68               NA              NA                     NA
    ## 69               NA              NA                     NA
    ## 70               NA              NA                     NA
    ##     isolationWindowTargetMz isolationWindowUpperMz compound_id experiment_id
    ##                   <numeric>              <numeric> <character>     <integer>
    ## 1                        NA                     NA           1             4
    ## 2                        NA                     NA           2             4
    ## 3                        NA                     NA           3             1
    ## 4                        NA                     NA           4             4
    ## 5                        NA                     NA           5             2
    ## ...                     ...                    ...         ...           ...
    ## 66                       NA                     NA          66             2
    ## 67                       NA                     NA          67             4
    ## 68                       NA                     NA          68             3
    ## 69                       NA                     NA          69             4
    ## 70                       NA                     NA          70             1
    ##     experiment_name experiment_date experiment_operator     formula exactmass
    ##         <character>     <character>         <character> <character> <numeric>
    ## 1             exp_d      2024-12-07                  AA    C10H10O3   178.063
    ## 2             exp_d      2024-12-07                  AA    C10H10O3   178.063
    ## 3             exp_a      2022-11-24                  AA    C10H10O3   178.063
    ## 4             exp_d      2024-12-07                  AA    C10H10O3   178.063
    ## 5             exp_b      2023-06-12                  VV    C10H10O3   178.063
    ## ...             ...             ...                 ...         ...       ...
    ## 66            exp_b      2023-06-12                  VV    C16H16O8   336.084
    ## 67            exp_d      2024-12-07                  AA    C16H16O8   336.084
    ## 68            exp_c      2023-09-21                  AA    C16H16O8   336.084
    ## 69            exp_d      2024-12-07                  AA    C16H16O8   336.084
    ## 70            exp_a      2022-11-24                  AA    C16H16O8   336.084
    ##                     smiles                  inchi               inchikey
    ##                <character>            <character>            <character>
    ## 1   CC1CC2=C(C(=CC=C2)O).. InChI=1S/C10H10O3/c1.. KWILGNNWGSNMPA-UHFFF..
    ## 2   CC1CC2=C(C(=CC=C2)O).. InChI=1S/C10H10O3/c1.. KWILGNNWGSNMPA-UHFFF..
    ## 3   CC1CC2=C(C(=CC=C2)O).. InChI=1S/C10H10O3/c1.. KWILGNNWGSNMPA-UHFFF..
    ## 4   CC1CC2=C(C(=CC=C2)O).. InChI=1S/C10H10O3/c1.. KWILGNNWGSNMPA-UHFFF..
    ## 5   CC1CC2=C(C(=CC=C2)O).. InChI=1S/C10H10O3/c1.. KWILGNNWGSNMPA-UHFFF..
    ## ...                    ...                    ...                    ...
    ## 66  C[C@]1([C@@H]([C@H](.. InChI=1S/C16H16O8/c1.. VSMBLBOUQJNJIL-JJXSE..
    ## 67  C[C@]1([C@@H]([C@H](.. InChI=1S/C16H16O8/c1.. VSMBLBOUQJNJIL-JJXSE..
    ## 68  C[C@]1([C@@H]([C@H](.. InChI=1S/C16H16O8/c1.. VSMBLBOUQJNJIL-JJXSE..
    ## 69  C[C@]1([C@@H]([C@H](.. InChI=1S/C16H16O8/c1.. VSMBLBOUQJNJIL-JJXSE..
    ## 70  C[C@]1([C@@H]([C@H](.. InChI=1S/C16H16O8/c1.. VSMBLBOUQJNJIL-JJXSE..
    ##             cas     pubchem           name   accession          spectrum_name
    ##     <character> <character>    <character> <character>            <character>
    ## 1    17397-85-2   CID:28516        Mellein    AC000001 Mellein; LC-ESI-ITFT..
    ## 2    17397-85-2   CID:28516        Mellein    AC000002 Mellein; LC-ESI-ITFT..
    ## 3    17397-85-2   CID:28516        Mellein    AC000003 Mellein; LC-ESI-ITFT..
    ## 4    17397-85-2   CID:28516        Mellein    AC000004 Mellein; LC-ESI-ITFT..
    ## 5    17397-85-2   CID:28516        Mellein    AC000005 Mellein; LC-ESI-ITFT..
    ## ...         ...         ...            ...         ...                    ...
    ## 66   22268-16-2   CID:89644 Altersolanol A    AC000066 Altersolanol A; LC-E..
    ## 67   22268-16-2   CID:89644 Altersolanol A    AC000067 Altersolanol A; LC-E..
    ## 68   22268-16-2   CID:89644 Altersolanol A    AC000068 Altersolanol A; LC-E..
    ## 69   22268-16-2   CID:89644 Altersolanol A    AC000069 Altersolanol A; LC-E..
    ## 70   22268-16-2   CID:89644 Altersolanol A    AC000070 Altersolanol A; LC-E..
    ##            date                authors     license          copyright
    ##     <character>            <character> <character>        <character>
    ## 1    2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## 2    2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## 3    2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## 4    2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## 5    2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## ...         ...                    ...         ...                ...
    ## 66   2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## 67   2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## 68   2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## 69   2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ## 70   2017.07.07 Justin B. Renaud, Ma..    CC BY-SA Copyright (C) 2017
    ##                publication                 splash precursor_intensity
    ##                <character>            <character>           <numeric>
    ## 1   Renaud, J. B.; Sumar.. splash10-03fr-090000..             161.059
    ## 2   Renaud, J. B.; Sumar.. splash10-03fr-090000..             161.059
    ## 3   Renaud, J. B.; Sumar.. splash10-03fr-090000..             161.059
    ## 4   Renaud, J. B.; Sumar.. splash10-03di-090000..             161.059
    ## 5   Renaud, J. B.; Sumar.. splash10-01q9-090000..             133.064
    ## ...                    ...                    ...                 ...
    ## 66  Renaud, J. B.; Sumar.. splash10-0uk9-008900..             301.070
    ## 67  Renaud, J. B.; Sumar.. splash10-0uk9-009500..             301.070
    ## 68  Renaud, J. B.; Sumar.. splash10-0pi1-009200..             273.075
    ## 69  Renaud, J. B.; Sumar.. splash10-0pi1-009100..             245.080
    ## 70  Renaud, J. B.; Sumar.. splash10-0kus-029000..             245.080
    ##          adduct  ionization ionization_voltage fragmentation_mode
    ##     <character> <character>        <character>        <character>
    ## 1        [M+H]+         ESI             3.9 kV                HCD
    ## 2        [M+H]+         ESI             3.9 kV                HCD
    ## 3        [M+H]+         ESI             3.9 kV                HCD
    ## 4        [M+H]+         ESI             3.9 kV                HCD
    ## 5        [M+H]+         ESI             3.9 kV                HCD
    ## ...         ...         ...                ...                ...
    ## 66       [M+H]+         ESI             3.9 kV                HCD
    ## 67       [M+H]+         ESI             3.9 kV                HCD
    ## 68       [M+H]+         ESI             3.9 kV                HCD
    ## 69       [M+H]+         ESI             3.9 kV                HCD
    ## 70       [M+H]+         ESI             3.9 kV                HCD
    ##     collisionEnergy_text             instrument instrument_type
    ##              <character>            <character>     <character>
    ## 1                10(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## 2                20(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## 3                30(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## 4                35(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## 5                50(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## ...                  ...                    ...             ...
    ## 66               10(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## 67               20(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## 68               30(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## 69               35(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ## 70               50(NCE) Q-Exactive Orbitrap ..     LC-ESI-ITFT
    ##     precursorMz_text spectrum_id predicted msms_mz_range_min msms_mz_range_max
    ##          <character>   <integer> <logical>         <numeric>         <numeric>
    ## 1           179.0697           1        NA           133.065            179.07
    ## 2           179.0697           2        NA           133.065            179.07
    ## 3           179.0697           3        NA           105.070            179.07
    ## 4           179.0697           4        NA           105.070            179.07
    ## 5           179.0697           5        NA           105.070            179.07
    ## ...              ...         ...       ...               ...               ...
    ## 66          337.0912          66        NA          245.0808           337.092
    ## 67          337.0912          67        NA          245.0808           319.081
    ## 68          337.0912          68        NA          151.0390           301.071
    ## 69          337.0912          69        NA          151.0390           301.071
    ## 70          337.0912          70        NA           95.0491           301.071
    ##                                     synonym
    ##                             <CharacterList>
    ## 1   Mellein,Ochracin,8-hydroxy-3-methyl-3..
    ## 2   Mellein,Ochracin,8-hydroxy-3-methyl-3..
    ## 3   Mellein,Ochracin,8-hydroxy-3-methyl-3..
    ## 4   Mellein,Ochracin,8-hydroxy-3-methyl-3..
    ## 5   Mellein,Ochracin,8-hydroxy-3-methyl-3..
    ## ...                                     ...
    ## 66    Altersolanol A,(1S,2R,3S,4R)-1,2,3,..
    ## 67    Altersolanol A,(1S,2R,3S,4R)-1,2,3,..
    ## 68    Altersolanol A,(1S,2R,3S,4R)-1,2,3,..
    ## 69    Altersolanol A,(1S,2R,3S,4R)-1,2,3,..
    ## 70    Altersolanol A,(1S,2R,3S,4R)-1,2,3,..

This `DataFrame` includes thus also the columns `"experiment_id"`,
`"experiment_name"`, `"experiment_date"` and `"experiment_operator"`
from the new *experiment* table. Similarly, columns from this table can
also be extracted using `$`, thus returning the values from the
respective columns associated to the individual spectra.

``` r

s$experiment_name
```

    ##  [1] "exp_d" "exp_d" "exp_a" "exp_d" "exp_b" "exp_c" "exp_a" "exp_a" "exp_d"
    ## [10] "exp_b" "exp_c" "exp_b" "exp_c" "exp_b" "exp_b" "exp_b" "exp_c" "exp_b"
    ## [19] "exp_a" "exp_b" "exp_c" "exp_b" "exp_a" "exp_c" "exp_b" "exp_c" "exp_d"
    ## [28] "exp_c" "exp_b" "exp_d" "exp_b" "exp_d" "exp_c" "exp_c" "exp_d" "exp_a"
    ## [37] "exp_a" "exp_a" "exp_a" "exp_c" "exp_d" "exp_b" "exp_b" "exp_c" "exp_b"
    ## [46] "exp_a" "exp_a" "exp_a" "exp_b" "exp_a" "exp_c" "exp_b" "exp_a" "exp_d"
    ## [55] "exp_a" "exp_c" "exp_d" "exp_c" "exp_c" "exp_d" "exp_d" "exp_a" "exp_c"
    ## [64] "exp_b" "exp_b" "exp_b" "exp_d" "exp_c" "exp_d" "exp_a"

``` r

file.remove(tf)
```

    ## [1] TRUE

## Session information

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
    ##  [1] RSQLite_2.4.5           MsBackendMgf_1.19.0     Spectra_1.21.1         
    ##  [4] BiocParallel_1.45.0     CompoundDb_1.15.2       S4Vectors_0.49.0       
    ##  [7] BiocGenerics_0.57.0     generics_0.1.4          AnnotationFilter_1.35.0
    ## [10] BiocStyle_2.39.0       
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
    ## [49] base64enc_0.1-3        utf8_1.2.6             ChemmineR_3.63.0      
    ## [52] scales_1.4.0           bit64_4.6.0-1          rmarkdown_2.30        
    ## [55] bit_4.6.0              otel_0.2.0             gridExtra_2.3         
    ## [58] ragg_1.5.0             png_0.1-8              memoise_2.0.1         
    ## [61] evaluate_1.0.5         knitr_1.51             GenomicRanges_1.63.1  
    ## [64] IRanges_2.45.0         rlang_1.1.7            Rcpp_1.1.1            
    ## [67] glue_1.8.0             DBI_1.2.3              xml2_1.5.2            
    ## [70] BiocManager_1.30.27    jsonlite_2.0.0         R6_2.6.1              
    ## [73] ProtGenerics_1.39.2    systemfonts_1.3.1      fs_1.6.6              
    ## [76] MsCoreUtils_1.23.2

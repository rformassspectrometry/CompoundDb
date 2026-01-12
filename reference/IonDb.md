# IonDb: compound database with additional ion information

`IonDb` objects extends `CompDb` by allowing to store also information
about measured ions to a
[`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
database. This information includes the type (adduct) of the ion, it's
measured (or expected) retention time for a certain LC-MS setup and its
mass-to-charge ratio.

As suggested use case, users might create (or download) a `CompDb`
(SQLite) database e.g. containing compound (and eventually MS/MS
spectra) annotations from public databases such as the Human Metabolome
Database (HMDB) or MassBank. To store now measured ions (e.g. of
lab-internal standards) for a certain LC-MS setup, such a `CompDb` can
then be converted to an `IonDb` using the `IonDb()` constructor
function. Ions can be subsequently added using the `insertIon()`
function. In general, it is suggested to create one `IonDb` database for
one specific LC-MS setup. Such an `IonDb` database can then be used to
match experimental m/z and retention times against ions defined in the
database (using the functionality of the
[MetaboAnnotation](https://rformassspectrometry.github.io/MetaboAnnotation)
package).

## Usage

``` r
# S4 method for class 'IonDb'
ionVariables(object, includeId = FALSE, ...)

# S4 method for class 'IonDb'
ions(
  object,
  columns = ionVariables(object),
  filter,
  return.type = c("data.frame", "tibble"),
  ...
)

# S4 method for class 'IonDb'
insertIon(object, ions, addColumns = FALSE)

# S4 method for class 'IonDb'
deleteIon(object, ids = integer(0), ...)

# S4 method for class 'missing,missing'
IonDb(x, cdb, flags = SQLITE_RWC, ...)

# S4 method for class 'CompDb,missing'
IonDb(x, cdb, ions = data.frame(), ...)

# S4 method for class 'character,missing'
IonDb(x, cdb, flags = SQLITE_RW, ...)

# S4 method for class 'DBIConnection,missing'
IonDb(
  x,
  cdb,
  ions = data.frame(),
  flags = SQLITE_RW,
  ...,
  .DBNAME = character()
)

# S4 method for class 'character,CompDb'
IonDb(x, cdb, ions = data.frame(), flags = SQLITE_RW, ...)

# S4 method for class 'DBIConnection,CompDb'
IonDb(
  x,
  cdb,
  ions = data.frame(),
  flags = SQLITE_RW,
  ...,
  .DBNAME = character()
)
```

## Arguments

- object:

  For all methods: a `IonDb` object.

- includeId:

  For `ionVariables()`: `logical(1)` whether the ion ID (column
  `"ion_id"`) should be included in the result. The default is
  `includeId = FALSE`.

- ...:

  additional arguments. Currently not used.

- columns:

  For `ions()`: `character` with the names of the database columns that
  should be retrieved. Use `ionVariables` for a list of available column
  names.

- filter:

  For `ions()`: filter expression or
  [`AnnotationFilter::AnnotationFilter()`](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html)
  defining a filter to be used to retrieve specific elements from the
  database.

- return.type:

  For `ions()`: either `"data.frame"` or `"tibble"` to return the result
  as a [`data.frame()`](https://rdrr.io/r/base/data.frame.html) or
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html),
  respectively.

- ions:

  for `insertIon()` and `IonDb()`: `data.frame` with ion definitions to
  be added to the `IonDb` database. Columns `"compound_id"`
  ([`character()`](https://rdrr.io/r/base/character.html)),
  `"ion_adduct"`
  ([`character()`](https://rdrr.io/r/base/character.html)), `"ion_mz"`
  ([`numeric()`](https://rdrr.io/r/base/numeric.html)) and `"ion_rt"`
  ([`numeric()`](https://rdrr.io/r/base/numeric.html)) are mandatory
  (but, with the exception of `"compound_id"`, can contain `NA`).

- addColumns:

  For `insertIons()`: `logical(1)` whether columns being present in the
  submitted `data.frame` but not in the database table should be added
  to the database's ion table.

- ids:

  For `deleteIon()`:
  [`character()`](https://rdrr.io/r/base/character.html) or
  (alternatively [`integer()`](https://rdrr.io/r/base/integer.html))
  specifying the IDs of the ions to delete. IDs in `ids` that are not
  associated to any ion in the `IonDb` object are ignored.

- x:

  For `IonDb()`: database connection or `character(1)` with the file
  name of the SQLite database where the `IonDb` data will be stored or a
  [`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  object that should be converted into an `IonDb` object.

      For all other methods: an `IonDb` object.

- cdb:

  For `IonDb()`: `CompDb` object from which data should be transferred
  to the `IonDb` database.

- flags:

  For `IonDb()`: optional `integer(1)` defining the flags for the SQLite
  database connection. Only used if `x` is a
  [`character()`](https://rdrr.io/r/base/character.html).

- .DBNAME:

  `character(1)` defining the SQLite database file. This is an internal
  parameter not intended to be used/provided by the user.

## Value

See description of the respective function.

## Creation of `IonDb` objects/databases

- A new `IonDb` database can be created and initialized with data from
  an existing `CompDb` database by passing either the database
  connection (e.g. an `SQLiteConnection`) or the file path of a (to be
  created) SQLite database with parameter `x` to the `IonDb()` function
  and the `CompDb` object with parameter `cdb`. Optional parameter
  `ions` allows insert in addition ion definitions (which can also be
  added later using `insertIon()` function calls).

- An existing `CompDb` can be converted to an `IonDb` by passing the
  [`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  object with parameter `x` to the `IonDb` function. Optional parameter
  `ions` allows to provide a `data.frame` with ion definitions to be
  inserted in to the database (which can also be added later using
  `insertIon()` function calls). Note that this fails if the database
  connection for the `CompDb` is read-only.

- Previously created `IonDb` databases can be loaded by passing either
  the database connection (e.g. an `SQLiteConnection`) or the file path
  of the (SQLite) database with parameter `x` to the `IonDb()` function.

## Retrieve annotations and ion information from the database

Annotations/compound informations can be retrieved from a `IonDb` in the
same way as thay are extracted from a `CompDb`. In addition, the
function `ions()` allows to retrieve the specific ion information from
the database. It returns the actual data as a `data.frame` (if
`return.type = "data.frame"`) or a
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
(if `return.type = "tibble"`). An `ions()` call will always return all
elements from the *ms_ion* table (unless a `filter` is used).

## General functions (beside those inherited from `CompDb`)

- `IonDb()`: connect to or create a compound/ion database.

- `ionVariables()`: returns all available columns/database fields for
  ions.

## Adding and removing data from a database

`IonDb` inherits the
[`insertCompound()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md),
[`insertSpectra()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md),
[`deleteCompound()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
and
[`deleteSpectra()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
functions from
[`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md).
In addition, `IonDb` defines the functions:

- `insertIon()`: adds ions to the `IonDb` object. Note that
  `insertIon()` always adds all the ions specified through the `ions`
  parameter and does not check if they are already in the database. To
  add columns present in the submitted `data.frame` to the database
  table set `addColumns = TRUE` (default is `addColumns = FALSE`).

- `deleteIon()`: deletes ions from the `IonDb` object by specifying
  their IDs.

## Filtering the database

Like
[`compounds()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
and [`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) also
`ions()` allows to filter the results using specific filter classes and
expressions. Filtering uses the concepts from Bioconductor's
*AnnotationFilter* package. All information for a certain compound with
the ID `"1"` can for example be retrieved by passing the filter
expression `filter = ~ ion_id == 1` to the `ions()` function.

Use the
[`AnnotationFilter::supportedFilters()`](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html)
function on the `IonDb` object to get a list of all supported filters.
See also examples below or the usage vignette for details.

## Author

Andrea Vicini, Johannes Rainer

## Examples

``` r

# We load a small compound test database based on MassBank which is
# distributed with this package.
cdb <- CompDb(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"))
cdb
#> class: CompDb 
#>  data source: MassBank 
#>  version: 2020.09 
#>  organism: NA 
#>  compound count: 70 
#>  MS/MS spectra count: 70 

# We next want to convert this CompDb into an IonDb, but the original CompDb
# database is read only, thus we have to provide the name (or connection) of
# an other database to transfer all the data from the CompDb to that.
idb <- IonDb(paste0(tempdir(), "/idb_ex.db"), cdb)
idb
#> class: IonDb 
#>  data source: MassBank 
#>  version: 2020.09 
#>  organism: NA 
#>  compound count: 70 
#>  MS/MS spectra count: 70 
#>  ion count: 0 

# It is also possible to load a previously created IonDb passing only the
# connection to the database.
idb2 <- IonDb(paste0(tempdir(), "/idb_ex.db"))

# Ion definitions can be added to the database with the `insertIon` function
# providing a `data.frame` with ion definition. This `data.frame` is expected
# to provide the IDs of the compounds, an adduct name/definition and the
# (experimentally determined) m/z and retention time of the ion. To list
# compound IDs from the CompDb database:
head(compounds(cdb, c("compound_id", "name")))
#>   compound_id         name
#> 1           1      Mellein
#> 2           2      Mellein
#> 3           3      Mellein
#> 4           4      Mellein
#> 5           5      Mellein
#> 6           6 AAL toxin TB

ions = data.frame(compound_id = c("1", "1", "2", "3", "6", "35"),
                  ion_adduct = c("[M+H]+", "[M+Na]+", "[M+Na]+",
                                 "[M+Na]+", "[M+2H]2+", "[M+H-NH3]+"),
                  ion_mz = c(179.0703, 201.0522, 201.0522,
                             201.0522, 253.66982, 312.0390),
                  ion_rt = 1:6)

# Inserting ion definitions.
idb <- insertIon(idb, ions)
idb
#> class: IonDb 
#>  data source: MassBank 
#>  version: 2020.09 
#>  organism: NA 
#>  compound count: 70 
#>  MS/MS spectra count: 70 
#>  ion count: 6 

ions(idb, columns = c("name", "formula", "ion_adduct", "ion_mz", "ion_rt"))
#>   ion_adduct   ion_mz ion_rt         name   formula
#> 1     [M+H]+ 179.0703      1      Mellein  C10H10O3
#> 2    [M+Na]+ 201.0522      2      Mellein  C10H10O3
#> 3    [M+Na]+ 201.0522      3      Mellein  C10H10O3
#> 4    [M+Na]+ 201.0522      4      Mellein  C10H10O3
#> 5   [M+2H]2+ 253.6698      5 AAL toxin TB C25H47NO9
#> 6 [M+H-NH3]+ 312.0390      6 Aflatoxin G1  C17H12O7

## List all available ion variables
ionVariables(idb)
#> [1] "compound_id" "ion_adduct"  "ion_mz"      "ion_rt"     

## Extract a data.frame with ion variables for all ions
ions(idb)
#>   compound_id ion_adduct   ion_mz ion_rt
#> 1           1     [M+H]+ 179.0703      1
#> 2           1    [M+Na]+ 201.0522      2
#> 3           2    [M+Na]+ 201.0522      3
#> 4           3    [M+Na]+ 201.0522      4
#> 5          35 [M+H-NH3]+ 312.0390      6
#> 6           6   [M+2H]2+ 253.6698      5

## List all database tables and their columns
tables(idb)
#> $ms_compound
#> [1] "compound_id" "formula"     "exactmass"   "smiles"      "inchi"      
#> [6] "inchikey"    "cas"         "pubchem"     "name"       
#> 
#> $ms_ion
#> [1] "ion_id"      "compound_id" "ion_adduct"  "ion_mz"      "ion_rt"     
#> 
#> $msms_spectrum
#>  [1] "accession"             "spectrum_name"         "date"                 
#>  [4] "authors"               "license"               "copyright"            
#>  [7] "publication"           "ms_level"              "polarity"             
#> [10] "splash"                "compound_id"           "precursor_intensity"  
#> [13] "precursor_mz"          "adduct"                "ionization"           
#> [16] "ionization_voltage"    "fragmentation_mode"    "collision_energy_text"
#> [19] "instrument"            "instrument_type"       "precursor_mz_text"    
#> [22] "spectrum_id"           "collision_energy"      "predicted"            
#> [25] "msms_mz_range_min"     "msms_mz_range_max"    
#> 
#> $msms_spectrum_peak
#> [1] "spectrum_id" "mz"          "intensity"   "peak_id"    
#> 
#> $synonym
#> [1] "compound_id" "synonym"    
#> 

## Filtering the database
##
## Get all ions with an m/z between 200 and 300
res <- ions(idb, filter = ~ ion_mz > 200 & ion_mz < 300)
res
#>   compound_id ion_adduct   ion_mz ion_rt
#> 1           1    [M+Na]+ 201.0522      2
#> 2           2    [M+Na]+ 201.0522      3
#> 3           3    [M+Na]+ 201.0522      4
#> 4           6   [M+2H]2+ 253.6698      5

## Get all ions that have a H in their adduct definition.
res <- ions(idb, filter = IonAdductFilter("H", "contains"))
res
#>   compound_id ion_adduct   ion_mz ion_rt
#> 1           1     [M+H]+ 179.0703      1
#> 2          35 [M+H-NH3]+ 312.0390      6
#> 3           6   [M+2H]2+ 253.6698      5
```

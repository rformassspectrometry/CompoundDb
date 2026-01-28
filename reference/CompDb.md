# Simple compound (metabolite) databases

`CompDb` objects provide access to general (metabolite) compound
annotations along with *metadata* information such as the annotation's
source, date and release version. The data is stored internally in a
database (usually an SQLite database).

`hasMsMsSpectra` returns `TRUE` if MS/MS spectrum data is available in
the database and `FALSE` otherwise.

## Usage

``` r
CompDb(x, flags = SQLITE_RO)

hasMsMsSpectra(x)

src_compdb(x)

tables(x)

copyCompDb(x, y)

# S4 method for class 'CompDb'
dbconn(x)

# S4 method for class 'CompDb'
Spectra(object, filter, ...)

# S4 method for class 'CompDb'
supportedFilters(object)

# S4 method for class 'CompDb'
metadata(x, ...)

# S4 method for class 'CompDb'
spectraVariables(object, ...)

# S4 method for class 'CompDb'
compoundVariables(object, includeId = FALSE, ...)

# S4 method for class 'CompDb'
compounds(
  object,
  columns = compoundVariables(object),
  filter,
  return.type = c("data.frame", "tibble"),
  ...
)

# S4 method for class 'CompDb,Spectra'
insertSpectra(object, spectra, columns = spectraVariables(spectra), ...)

# S4 method for class 'CompDb'
deleteSpectra(object, ids = integer(0), ...)

# S4 method for class 'CompDb'
mass2mz(x, adduct = c("[M+H]+"), name = "formula")

# S4 method for class 'CompDb'
insertCompound(object, compounds = data.frame(), addColumns = FALSE)

# S4 method for class 'CompDb'
deleteCompound(object, ids = character(), recursive = FALSE, ...)
```

## Arguments

- x:

  For `CompDb()`: `character(1)` with the file name of the SQLite
  compound database. Alternatively it is possible to provide the
  connection to the database with parameter `x`. For `copyCompDb()`:
  either a `CompDb` or a database connection.

      For all other methods: a `CompDb` object.

- flags:

  flags passed to the SQLite database connection. See
  [`RSQLite::SQLite()`](https://rsqlite.r-dbi.org/reference/SQLite.html).
  Defaults to read-only, i.e.
  [`RSQLite::SQLITE_RO`](https://rsqlite.r-dbi.org/reference/SQLite.html).

- y:

  For `copyCompDb()`: connection to a database to which the content
  should be copied.

- object:

  For all methods: a `CompDb` object.

- filter:

  For `compounds()` and
  [`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html): filter
  expression or
  [`AnnotationFilter::AnnotationFilter()`](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html)
  defining a filter to be used to retrieve specific elements from the
  database.

- ...:

  additional arguments. Currently not used.

- includeId:

  for `compoundVariables()`: `logical(1)` whether the comound ID (column
  `"compound_id"`) should be included in the result. The default is
  `includeIds = FALSE`.

- columns:

  For `compounds()`, `Spectra`: `character` with the names of the
  database columns that should be retrieved. Use `compoundVariables()`
  and/or
  [`spectraVariables()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  for a list of available column names. For `insertSpectra()`: columns
  (spectra variables) that should be inserted into the database (to
  avoid inserting all variables).

- return.type:

  For `compounds()`: either `"data.frame"` or `"tibble"` to return the
  result as a [`data.frame()`](https://rdrr.io/r/base/data.frame.html)
  or
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html),
  respectively.

- spectra:

  For `insertSpectra()`: `Spectra` object containing the spectra to be
  added to the `IonDb` database.

- ids:

  For `deleteSpectra()`:
  [`integer()`](https://rdrr.io/r/base/integer.html) specifying the IDs
  of the spectra to delete. IDs in `ids` that are not associated to any
  spectra in the `CompDb` object are ignored. For `deleteCompound`:
  [`character()`](https://rdrr.io/r/base/character.html) with the
  compound IDs to be deleted.

- adduct:

  either a `character` specifying the name(s) of the adduct(s) for which
  the m/z should be calculated or a `data.frame` with the adduct
  definition. See
  [`adductNames()`](https://rdrr.io/pkg/MetaboCoreUtils/man/adductNames.html)
  for supported adduct names and the description for more information on
  the expected format if a `data.frame` is provided.

- name:

  For `mass2mz()`: `character(1)`. Defines the `CompDb` column that will
  be used to name/identify the returned m/z values. By default
  (`name = "formula"`) m/z values for all unique molecular formulas are
  calculated and these are used as `rownames` for the returned `matrix`.
  With `name = "compound_id"` the adduct m/z for all compounds (even
  those with equal formulas) are calculated and returned.

- compounds:

  For `insertCompound()`: `data.frame` with compound data to be inserted
  into a `CompDb` database. See function description for details.

- addColumns:

  For `insertCompound()`: `logical(1)` whether all (extra) columns in
  parameter `compounds` should be stored also in the database table. The
  default is `addColumns = FALSE`.

- recursive:

  For `deleteCompound()`: `logical(1)` whether also MS2 spectra
  associated with the compounds should be deleted.

## Value

See description of the respective function.

## Details

`CompDb` objects should be created using the constructor function
`CompDb()` providing the name of the (SQLite) database file providing
the compound annotation data.

## Retrieve annotations from the database

Annotations/compound informations can be retrieved from a `CompDb`
database with the `compounds()` and
[`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) functions:

- `compounds()` extracts compound data from the `CompDb` object. In
  contrast to `src_compdb` it returns the actual data as a `data.frame`
  (if `return.type = "data.frame"`) or a
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  (if `return.type = "tibble"`). A `compounds()` call will always return
  all elements from the *ms_compound* table (unless a `filter` is used).

- [`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) extract
  spectra from the database and returns them as a
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object from the *Spectra* package. Additional annotations requested
  with the `columns` parameter are added as additional spectra
  variables.

## General functions

- `CompDb()`: connect to a compound database.

- `compoundVariables()`: returns all available columns/database fields
  for compounds.

- `copyCompDb()`: allows to copy the content from a CompDb to another
  database. Parameter `x` is supposed to be either a `CompDb` or a
  database connection from which the data should be copied and `y` a
  connection to a database to which it should be copied.

- `dbconn()`: returns the connection (of type `DBIConnection`) to the
  database.

- [`metadata()`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html):
  returns general meta data of the compound database.

- [`spectraVariables()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  returns all spectra variables (i.e. columns) available in the
  `CompDb`.

- `src_compdb()` provides access to the `CompDb`'s database *via* the
  functionality from the `dplyr`/`dbplyr` package.

- `supportedFilters()`: provides an overview of the filters that can be
  applied on a `CompDb` object to extract only specific data from the
  database.

- `tables()`: returns a named `list` (names being table names) with the
  fields/columns from each table in the database.

- `mass2mz()`: calculates a table of the m/z values for each compound
  based on the provided set of adduct(s). Adduct definitions can be
  provided with parameter `adduct`. See
  [`MetaboCoreUtils::mass2mz()`](https://rdrr.io/pkg/MetaboCoreUtils/man/mass2mz.html)
  for more details. Parameter `name` defines the database table column
  that should be used as `rownames` of the returned `matrix`. By default
  `name = "formula"`, m/z values are calculated for each unique formula
  in the `CompDb` `x`.

## Adding and removing data from a database

Note that inserting and deleting data requires read-write access to the
database. Databases returned by `CompDb` are by default *read-only*. To
get write access `CompDb` should be called with parameter
`flags = RSQLite::SQLITE_RW`.

- `insertCompound()`: adds additional compound(s) to a `CompDb`. The
  compound(s) to be added can be specified with parameter `compounds`
  that is expected to be a `data.frame` with columns `"compound_id"`,
  `"name"`, `"inchi"`, `"inchikey"`, `"formula"`, `"exactmass"`. Column
  `"exactmass"` is expected to contain numeric values, all other columns
  `character`. Missing values are allowed for all columns except
  `"compound_id"`. An optional column `"synonyms"` can be used to
  provide alternative names for the compound. This column can contain a
  single `character` by row, or a `list` with multiple `character`
  (names) per row/compound (see examples below for details). By setting
  parameter `addColumns = TRUE` any additional columns in `compound`
  will be added to the database table. The default is
  `addColumns = FALSE`. The function returns the `CompDb` with the
  compounds added. See also
  [`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
  for more information and details on expected compound data and the
  examples below for general usage.

- `deleteCompound()`: removes specified compounds from the `CompDb`
  database. The IDs of the compounds that should be deleted need to be
  provided with parameter `ids`. To include compound IDs in the output
  of a `compounds()` call `"compound_id"` should be added to the
  `columns` parameter. By default an error is thrown if for some of the
  specified compounds also MS2 spectra are present in the database. To
  force deletion of the compounds along with all associated MS2 spectra
  use `recursive = TRUE`. See examples below for details. The function
  returns the updated `CompDb` database.

- `insertSpectra()`: adds further spectra to the database. The method
  always adds all the spectra specified through the `spectra` parameter
  and does not check if they are already in the database. Note that the
  input spectra must have the variable `compound_id` and only `Spectra`
  whose `compound_id` values are also in
  `compounds(object, "compound_id")` can be added. Parameter `columns`
  defines which spectra variables from the `spectra` should be inserted
  into the database. By default, all spectra variables are added but it
  is strongly suggested to specifically select (meaningful) spectra
  variables that should be stored in the database. Note that a spectra
  variable `"compound_id"` is mandatory. If needed, the function adds
  additional columns to the `msms_spectrum` database table. The function
  returns the updated `CompDb` object.

- `deleteSpectra()`: deletes specified spectra from the database. The
  IDs of the spectra to be deleted need to be provided with parameter
  `ids`.

Note that it is also possible to add new database tables and include
them in the data retrieval queries. See
[`addJoinDefinition()`](https://rformassspectrometry.github.io/CompoundDb/reference/addJoinDefinition.md)
for more information.

## Filtering the database

Data access methods such as `compounds()` and `Spectra` allow to filter
the results using specific filter classes and expressions. Filtering
uses the concepts from Bioconductor's `AnnotationFilter` package. All
information for a certain compound with the ID `"HMDB0000001"` can for
example be retrieved by passing the filter expression
`filter = ~ compound_id == "HMDB0000001"` to the `compounds` function.

Use the
[`AnnotationFilter::supportedFilters()`](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html)
function on the CompDb object to get a list of all supported filters.
See also examples below or the usage vignette for details.

## See also

- [`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
  for the function to create a SQLite compound database.

- [`CompoundIdFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  for filters that can be used on the `CompDb` database.

- [`addJoinDefinition()`](https://rformassspectrometry.github.io/CompoundDb/reference/addJoinDefinition.md)
  to expand a *CompDb* with additional, related, database tables.

## Author

Johannes Rainer

## Examples

``` r

## We load a small compound test database based on MassBank which is
## distributed with this package.
cdb <- CompDb(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"))
cdb
#> class: CompDb 
#>  data source: MassBank 
#>  version: 2020.09 
#>  organism: NA 
#>  compound count: 70 
#>  MS/MS spectra count: 70 

## Get general metadata information from the database, such as originating
## source and version:
metadata(cdb)
#>                 name                         value
#> 1             source                      MassBank
#> 2                url https://massbank.eu/MassBank/
#> 3     source_version                       2020.09
#> 4        source_date                    1603272565
#> 5           organism                          <NA>
#> 6   db_creation_date      Thu Oct 22 08:45:31 2020
#> 7 supporting_package                    CompoundDb
#> 8  supporting_object                        CompDb

## List all available compound annotations/fields
compoundVariables(cdb)
#> [1] "formula"   "exactmass" "smiles"    "inchi"     "inchikey"  "cas"      
#> [7] "pubchem"   "name"     

## Extract a data.frame with these annotations for all compounds
compounds(cdb)
#>       formula exactmass
#> 1    C10H10O3  178.0630
#> 2   C25H47NO9  505.3251
#> 3    C17H12O6  312.0634
#> 4    C17H14O6  314.0790
#> 5    C17H12O7  328.0583
#> 6    C17H14O7  330.0739
#> 7    C17H12O7  328.0583
#> 8  C20H20N2O3  336.1474
#> 9    C15H16O6  292.0947
#> 10   C14H10O5  258.0528
#> 11   C15H12O5  272.0685
#> 12   C16H16O8  336.0845
#>                                                                                   smiles
#> 1                                                            CC1CC2=C(C(=CC=C2)O)C(=O)O1
#> 2  CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 3                               COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 4                                COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 5                              COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 6                               COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 7                            COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4C(=C1)O[C@@H]5[C@]4(C=CO5)O
#> 8                    CC(=O)C1=C([C@@H]2[C@@H]3[C@@H](CC4=C5C3=CNC5=CC=C4)C(N2C1=O)(C)C)O
#> 9                                 C[C@]12C[C@@H]([C@H](C=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O)O
#> 10                                              CC1=CC(=CC2=C1C3=CC(=CC(=C3C(=O)O2)O)O)O
#> 11                                             CC1=CC(=CC2=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O
#> 12                    C[C@]1([C@@H]([C@H](C2=C([C@H]1O)C(=O)C3=CC(=CC(=C3C2=O)O)OC)O)O)O
#>                                                                                                                                                                                                            inchi
#> 1                                                                                                                                        InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-6/h2-4,6,11H,5H2,1H3
#> 2  InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 3                                                                                 InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 4                                                                                   InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 5                                                                             InChI=1S/C17H12O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h3,5-6,8,17H,2,4H2,1H3/t8-,17+/m0/s1
#> 6                                                                                 InChI=1S/C17H14O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 7                                                                            InChI=1S/C17H12O7/c1-21-9-6-10-13(17(20)4-5-22-16(17)23-10)14-12(9)7-2-3-8(18)11(7)15(19)24-14/h4-6,16,20H,2-3H2,1H3/t16-,17-/m1/s1
#> 8                                                      InChI=1S/C20H20N2O3/c1-9(23)14-18(24)17-16-11-8-21-13-6-4-5-10(15(11)13)7-12(16)20(2,3)22(17)19(14)25/h4-6,8,12,16-17,21,24H,7H2,1-3H3/t12-,16+,17+/m1/s1
#> 9                                                                                InChI=1S/C15H16O6/c1-15-6-12(18)10(16)5-9(15)8-3-7(20-2)4-11(17)13(8)14(19)21-15/h3-5,10,12,16-18H,6H2,1-2H3/t10-,12-,15-/m0/s1
#> 10                                                                                                                     InChI=1S/C14H10O5/c1-6-2-7(15)5-11-12(6)9-3-8(16)4-10(17)13(9)14(18)19-11/h2-5,15-17H,1H3
#> 11                                                                                                               InChI=1S/C15H12O5/c1-7-3-8(16)4-12-13(7)10-5-9(19-2)6-11(17)14(10)15(18)20-12/h3-6,16-17H,1-2H3
#> 12                                                                   InChI=1S/C16H16O8/c1-16(23)14(21)10-9(13(20)15(16)22)12(19)8-6(11(10)18)3-5(24-2)4-7(8)17/h3-4,13-15,17,20-23H,1-2H3/t13-,14+,15+,16-/m0/s1
#>                       inchikey         cas       pubchem
#> 1  KWILGNNWGSNMPA-UHFFFAOYSA-N  17397-85-2     CID:28516
#> 2  CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 3  OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 4  WWSYXEZEXMQWHT-WNWIJWBNSA-N   7220-81-7   CID:2724360
#> 5  XWIYFDMXXLINPU-WNWIJWBNSA-N   1165-39-5   CID:2724361
#> 6  WPCVRWVBBXIRMA-WNWIJWBNSA-N   7241-98-7   CID:2724362
#> 7  MJBWDEQAUQTVKK-IAGOWNOFSA-N   6795-23-9  CID:15558498
#> 8  SZINUGQCTHLQAZ-DQYPLSBCSA-N  18172-33-3  CID:54682463
#> 9  MMHTXEATDNFMMY-WBIUFABUSA-N  29752-43-0     CID:34687
#> 10 CEBXXEKPIIDJHL-UHFFFAOYSA-N    641-38-3   CID:5359485
#> 11 LCSDQFNUYFTXMT-UHFFFAOYSA-N  23452-05-3   CID:5360741
#> 12 VSMBLBOUQJNJIL-JJXSEGSLSA-N  22268-16-2     CID:89644
#>                        name
#> 1                   Mellein
#> 2              AAL toxin TB
#> 3              Aflatoxin B1
#> 4              Aflatoxin B2
#> 5              Aflatoxin G1
#> 6              Aflatoxin G2
#> 7              Aflatoxin M1
#> 8  alpha-Cyclopiazonic acid
#> 9                 Altenuene
#> 10              Alternariol
#> 11 Alternariol methyl ether
#> 12           Altersolanol A

## Note that the `compounds` function will by default always return a
## data frame of **unique** entries for the specified columns. Including
## also the `"compound_id"` to the requested columns will ensure that all
## data is returned from the tables.
compounds(cdb, columns = c("compound_id", compoundVariables(cdb)))
#>    compound_id    formula exactmass
#> 1            1   C10H10O3  178.0630
#> 2            2   C10H10O3  178.0630
#> 3            3   C10H10O3  178.0630
#> 4            4   C10H10O3  178.0630
#> 5            5   C10H10O3  178.0630
#> 6            6  C25H47NO9  505.3251
#> 7            7  C25H47NO9  505.3251
#> 8            8  C25H47NO9  505.3251
#> 9            9  C25H47NO9  505.3251
#> 10          10  C25H47NO9  505.3251
#> 11          11  C25H47NO9  505.3251
#> 12          12  C25H47NO9  505.3251
#> 13          13  C25H47NO9  505.3251
#> 14          14  C25H47NO9  505.3251
#> 15          15  C25H47NO9  505.3251
#> 16          16   C17H12O6  312.0634
#> 17          17   C17H12O6  312.0634
#> 18          18   C17H12O6  312.0634
#> 19          19   C17H12O6  312.0634
#> 20          20   C17H12O6  312.0634
#> 21          21   C17H12O6  312.0634
#> 22          22   C17H12O6  312.0634
#> 23          23   C17H12O6  312.0634
#> 24          24   C17H12O6  312.0634
#> 25          25   C17H12O6  312.0634
#> 26          26   C17H14O6  314.0790
#> 27          27   C17H14O6  314.0790
#> 28          28   C17H14O6  314.0790
#> 29          29   C17H14O6  314.0790
#> 30          30   C17H14O6  314.0790
#> 31          31   C17H14O6  314.0790
#> 32          32   C17H12O7  328.0583
#> 33          33   C17H12O7  328.0583
#> 34          34   C17H12O7  328.0583
#> 35          35   C17H12O7  328.0583
#> 36          36   C17H12O7  328.0583
#> 37          37   C17H12O7  328.0583
#> 38          38   C17H14O7  330.0739
#> 39          39   C17H14O7  330.0739
#> 40          40   C17H14O7  330.0739
#> 41          41   C17H14O7  330.0739
#> 42          42   C17H14O7  330.0739
#> 43          43   C17H14O7  330.0739
#> 44          44   C17H12O7  328.0583
#> 45          45   C17H12O7  328.0583
#> 46          46   C17H12O7  328.0583
#> 47          47   C17H12O7  328.0583
#> 48          48   C17H12O7  328.0583
#> 49          49   C17H12O7  328.0583
#> 50          50 C20H20N2O3  336.1474
#> 51          51 C20H20N2O3  336.1474
#> 52          52 C20H20N2O3  336.1474
#> 53          53 C20H20N2O3  336.1474
#> 54          54 C20H20N2O3  336.1474
#> 55          55   C15H16O6  292.0947
#> 56          56   C15H16O6  292.0947
#> 57          57   C15H16O6  292.0947
#> 58          58   C15H16O6  292.0947
#> 59          59   C15H16O6  292.0947
#> 60          60   C14H10O5  258.0528
#> 61          61   C14H10O5  258.0528
#> 62          62   C14H10O5  258.0528
#> 63          63   C15H12O5  272.0685
#> 64          64   C15H12O5  272.0685
#> 65          65   C15H12O5  272.0685
#> 66          66   C16H16O8  336.0845
#> 67          67   C16H16O8  336.0845
#> 68          68   C16H16O8  336.0845
#> 69          69   C16H16O8  336.0845
#> 70          70   C16H16O8  336.0845
#>                                                                                   smiles
#> 1                                                            CC1CC2=C(C(=CC=C2)O)C(=O)O1
#> 2                                                            CC1CC2=C(C(=CC=C2)O)C(=O)O1
#> 3                                                            CC1CC2=C(C(=CC=C2)O)C(=O)O1
#> 4                                                            CC1CC2=C(C(=CC=C2)O)C(=O)O1
#> 5                                                            CC1CC2=C(C(=CC=C2)O)C(=O)O1
#> 6  CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 7  CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 8  CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 9  CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 10 CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 11 CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 12 CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 13 CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 14 CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 15 CC[C@@H](C)[C@H]([C@H](C[C@@H](C)CCCCCC[C@H](C[C@@H](CN)O)O)OC(=O)CC(CC(=O)O)C(=O)O)O
#> 16                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 17                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 18                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 19                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 20                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 21                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 22                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 23                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 24                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 25                              COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 26                               COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 27                               COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 28                               COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 29                               COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 30                               COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 31                               COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 32                             COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 33                             COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 34                             COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 35                             COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 36                             COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 37                             COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 38                              COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 39                              COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 40                              COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 41                              COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 42                              COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 43                              COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 44                           COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4C(=C1)O[C@@H]5[C@]4(C=CO5)O
#> 45                           COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4C(=C1)O[C@@H]5[C@]4(C=CO5)O
#> 46                           COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4C(=C1)O[C@@H]5[C@]4(C=CO5)O
#> 47                           COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4C(=C1)O[C@@H]5[C@]4(C=CO5)O
#> 48                           COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4C(=C1)O[C@@H]5[C@]4(C=CO5)O
#> 49                           COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4C(=C1)O[C@@H]5[C@]4(C=CO5)O
#> 50                   CC(=O)C1=C([C@@H]2[C@@H]3[C@@H](CC4=C5C3=CNC5=CC=C4)C(N2C1=O)(C)C)O
#> 51                   CC(=O)C1=C([C@@H]2[C@@H]3[C@@H](CC4=C5C3=CNC5=CC=C4)C(N2C1=O)(C)C)O
#> 52                   CC(=O)C1=C([C@@H]2[C@@H]3[C@@H](CC4=C5C3=CNC5=CC=C4)C(N2C1=O)(C)C)O
#> 53                   CC(=O)C1=C([C@@H]2[C@@H]3[C@@H](CC4=C5C3=CNC5=CC=C4)C(N2C1=O)(C)C)O
#> 54                   CC(=O)C1=C([C@@H]2[C@@H]3[C@@H](CC4=C5C3=CNC5=CC=C4)C(N2C1=O)(C)C)O
#> 55                                C[C@]12C[C@@H]([C@H](C=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O)O
#> 56                                C[C@]12C[C@@H]([C@H](C=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O)O
#> 57                                C[C@]12C[C@@H]([C@H](C=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O)O
#> 58                                C[C@]12C[C@@H]([C@H](C=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O)O
#> 59                                C[C@]12C[C@@H]([C@H](C=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O)O
#> 60                                              CC1=CC(=CC2=C1C3=CC(=CC(=C3C(=O)O2)O)O)O
#> 61                                              CC1=CC(=CC2=C1C3=CC(=CC(=C3C(=O)O2)O)O)O
#> 62                                              CC1=CC(=CC2=C1C3=CC(=CC(=C3C(=O)O2)O)O)O
#> 63                                             CC1=CC(=CC2=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O
#> 64                                             CC1=CC(=CC2=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O
#> 65                                             CC1=CC(=CC2=C1C3=CC(=CC(=C3C(=O)O2)O)OC)O
#> 66                    C[C@]1([C@@H]([C@H](C2=C([C@H]1O)C(=O)C3=CC(=CC(=C3C2=O)O)OC)O)O)O
#> 67                    C[C@]1([C@@H]([C@H](C2=C([C@H]1O)C(=O)C3=CC(=CC(=C3C2=O)O)OC)O)O)O
#> 68                    C[C@]1([C@@H]([C@H](C2=C([C@H]1O)C(=O)C3=CC(=CC(=C3C2=O)O)OC)O)O)O
#> 69                    C[C@]1([C@@H]([C@H](C2=C([C@H]1O)C(=O)C3=CC(=CC(=C3C2=O)O)OC)O)O)O
#> 70                    C[C@]1([C@@H]([C@H](C2=C([C@H]1O)C(=O)C3=CC(=CC(=C3C2=O)O)OC)O)O)O
#>                                                                                                                                                                                                            inchi
#> 1                                                                                                                                        InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-6/h2-4,6,11H,5H2,1H3
#> 2                                                                                                                                        InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-6/h2-4,6,11H,5H2,1H3
#> 3                                                                                                                                        InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-6/h2-4,6,11H,5H2,1H3
#> 4                                                                                                                                        InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-6/h2-4,6,11H,5H2,1H3
#> 5                                                                                                                                        InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-6/h2-4,6,11H,5H2,1H3
#> 6  InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 7  InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 8  InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 9  InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 10 InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 11 InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 12 InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 13 InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 14 InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 15 InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-18(25(33)34)12-22(29)30)11-16(2)9-7-5-6-8-10-19(27)14-20(28)15-26/h16-21,24,27-28,32H,4-15,26H2,1-3H3,(H,29,30)(H,33,34)/t16-,17+,18?,19+,20-,21-,24+/m0/s1
#> 16                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 17                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 18                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 19                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 20                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 21                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 22                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 23                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 24                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 25                                                                                InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 26                                                                                  InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 27                                                                                  InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 28                                                                                  InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 29                                                                                  InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 30                                                                                  InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 31                                                                                  InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 32                                                                            InChI=1S/C17H12O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h3,5-6,8,17H,2,4H2,1H3/t8-,17+/m0/s1
#> 33                                                                            InChI=1S/C17H12O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h3,5-6,8,17H,2,4H2,1H3/t8-,17+/m0/s1
#> 34                                                                            InChI=1S/C17H12O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h3,5-6,8,17H,2,4H2,1H3/t8-,17+/m0/s1
#> 35                                                                            InChI=1S/C17H12O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h3,5-6,8,17H,2,4H2,1H3/t8-,17+/m0/s1
#> 36                                                                            InChI=1S/C17H12O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h3,5-6,8,17H,2,4H2,1H3/t8-,17+/m0/s1
#> 37                                                                            InChI=1S/C17H12O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h3,5-6,8,17H,2,4H2,1H3/t8-,17+/m0/s1
#> 38                                                                                InChI=1S/C17H14O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 39                                                                                InChI=1S/C17H14O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 40                                                                                InChI=1S/C17H14O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 41                                                                                InChI=1S/C17H14O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 42                                                                                InChI=1S/C17H14O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 43                                                                                InChI=1S/C17H14O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 44                                                                           InChI=1S/C17H12O7/c1-21-9-6-10-13(17(20)4-5-22-16(17)23-10)14-12(9)7-2-3-8(18)11(7)15(19)24-14/h4-6,16,20H,2-3H2,1H3/t16-,17-/m1/s1
#> 45                                                                           InChI=1S/C17H12O7/c1-21-9-6-10-13(17(20)4-5-22-16(17)23-10)14-12(9)7-2-3-8(18)11(7)15(19)24-14/h4-6,16,20H,2-3H2,1H3/t16-,17-/m1/s1
#> 46                                                                           InChI=1S/C17H12O7/c1-21-9-6-10-13(17(20)4-5-22-16(17)23-10)14-12(9)7-2-3-8(18)11(7)15(19)24-14/h4-6,16,20H,2-3H2,1H3/t16-,17-/m1/s1
#> 47                                                                           InChI=1S/C17H12O7/c1-21-9-6-10-13(17(20)4-5-22-16(17)23-10)14-12(9)7-2-3-8(18)11(7)15(19)24-14/h4-6,16,20H,2-3H2,1H3/t16-,17-/m1/s1
#> 48                                                                           InChI=1S/C17H12O7/c1-21-9-6-10-13(17(20)4-5-22-16(17)23-10)14-12(9)7-2-3-8(18)11(7)15(19)24-14/h4-6,16,20H,2-3H2,1H3/t16-,17-/m1/s1
#> 49                                                                           InChI=1S/C17H12O7/c1-21-9-6-10-13(17(20)4-5-22-16(17)23-10)14-12(9)7-2-3-8(18)11(7)15(19)24-14/h4-6,16,20H,2-3H2,1H3/t16-,17-/m1/s1
#> 50                                                     InChI=1S/C20H20N2O3/c1-9(23)14-18(24)17-16-11-8-21-13-6-4-5-10(15(11)13)7-12(16)20(2,3)22(17)19(14)25/h4-6,8,12,16-17,21,24H,7H2,1-3H3/t12-,16+,17+/m1/s1
#> 51                                                     InChI=1S/C20H20N2O3/c1-9(23)14-18(24)17-16-11-8-21-13-6-4-5-10(15(11)13)7-12(16)20(2,3)22(17)19(14)25/h4-6,8,12,16-17,21,24H,7H2,1-3H3/t12-,16+,17+/m1/s1
#> 52                                                     InChI=1S/C20H20N2O3/c1-9(23)14-18(24)17-16-11-8-21-13-6-4-5-10(15(11)13)7-12(16)20(2,3)22(17)19(14)25/h4-6,8,12,16-17,21,24H,7H2,1-3H3/t12-,16+,17+/m1/s1
#> 53                                                     InChI=1S/C20H20N2O3/c1-9(23)14-18(24)17-16-11-8-21-13-6-4-5-10(15(11)13)7-12(16)20(2,3)22(17)19(14)25/h4-6,8,12,16-17,21,24H,7H2,1-3H3/t12-,16+,17+/m1/s1
#> 54                                                     InChI=1S/C20H20N2O3/c1-9(23)14-18(24)17-16-11-8-21-13-6-4-5-10(15(11)13)7-12(16)20(2,3)22(17)19(14)25/h4-6,8,12,16-17,21,24H,7H2,1-3H3/t12-,16+,17+/m1/s1
#> 55                                                                               InChI=1S/C15H16O6/c1-15-6-12(18)10(16)5-9(15)8-3-7(20-2)4-11(17)13(8)14(19)21-15/h3-5,10,12,16-18H,6H2,1-2H3/t10-,12-,15-/m0/s1
#> 56                                                                               InChI=1S/C15H16O6/c1-15-6-12(18)10(16)5-9(15)8-3-7(20-2)4-11(17)13(8)14(19)21-15/h3-5,10,12,16-18H,6H2,1-2H3/t10-,12-,15-/m0/s1
#> 57                                                                               InChI=1S/C15H16O6/c1-15-6-12(18)10(16)5-9(15)8-3-7(20-2)4-11(17)13(8)14(19)21-15/h3-5,10,12,16-18H,6H2,1-2H3/t10-,12-,15-/m0/s1
#> 58                                                                               InChI=1S/C15H16O6/c1-15-6-12(18)10(16)5-9(15)8-3-7(20-2)4-11(17)13(8)14(19)21-15/h3-5,10,12,16-18H,6H2,1-2H3/t10-,12-,15-/m0/s1
#> 59                                                                               InChI=1S/C15H16O6/c1-15-6-12(18)10(16)5-9(15)8-3-7(20-2)4-11(17)13(8)14(19)21-15/h3-5,10,12,16-18H,6H2,1-2H3/t10-,12-,15-/m0/s1
#> 60                                                                                                                     InChI=1S/C14H10O5/c1-6-2-7(15)5-11-12(6)9-3-8(16)4-10(17)13(9)14(18)19-11/h2-5,15-17H,1H3
#> 61                                                                                                                     InChI=1S/C14H10O5/c1-6-2-7(15)5-11-12(6)9-3-8(16)4-10(17)13(9)14(18)19-11/h2-5,15-17H,1H3
#> 62                                                                                                                     InChI=1S/C14H10O5/c1-6-2-7(15)5-11-12(6)9-3-8(16)4-10(17)13(9)14(18)19-11/h2-5,15-17H,1H3
#> 63                                                                                                               InChI=1S/C15H12O5/c1-7-3-8(16)4-12-13(7)10-5-9(19-2)6-11(17)14(10)15(18)20-12/h3-6,16-17H,1-2H3
#> 64                                                                                                               InChI=1S/C15H12O5/c1-7-3-8(16)4-12-13(7)10-5-9(19-2)6-11(17)14(10)15(18)20-12/h3-6,16-17H,1-2H3
#> 65                                                                                                               InChI=1S/C15H12O5/c1-7-3-8(16)4-12-13(7)10-5-9(19-2)6-11(17)14(10)15(18)20-12/h3-6,16-17H,1-2H3
#> 66                                                                   InChI=1S/C16H16O8/c1-16(23)14(21)10-9(13(20)15(16)22)12(19)8-6(11(10)18)3-5(24-2)4-7(8)17/h3-4,13-15,17,20-23H,1-2H3/t13-,14+,15+,16-/m0/s1
#> 67                                                                   InChI=1S/C16H16O8/c1-16(23)14(21)10-9(13(20)15(16)22)12(19)8-6(11(10)18)3-5(24-2)4-7(8)17/h3-4,13-15,17,20-23H,1-2H3/t13-,14+,15+,16-/m0/s1
#> 68                                                                   InChI=1S/C16H16O8/c1-16(23)14(21)10-9(13(20)15(16)22)12(19)8-6(11(10)18)3-5(24-2)4-7(8)17/h3-4,13-15,17,20-23H,1-2H3/t13-,14+,15+,16-/m0/s1
#> 69                                                                   InChI=1S/C16H16O8/c1-16(23)14(21)10-9(13(20)15(16)22)12(19)8-6(11(10)18)3-5(24-2)4-7(8)17/h3-4,13-15,17,20-23H,1-2H3/t13-,14+,15+,16-/m0/s1
#> 70                                                                   InChI=1S/C16H16O8/c1-16(23)14(21)10-9(13(20)15(16)22)12(19)8-6(11(10)18)3-5(24-2)4-7(8)17/h3-4,13-15,17,20-23H,1-2H3/t13-,14+,15+,16-/m0/s1
#>                       inchikey         cas       pubchem
#> 1  KWILGNNWGSNMPA-UHFFFAOYSA-N  17397-85-2     CID:28516
#> 2  KWILGNNWGSNMPA-UHFFFAOYSA-N  17397-85-2     CID:28516
#> 3  KWILGNNWGSNMPA-UHFFFAOYSA-N  17397-85-2     CID:28516
#> 4  KWILGNNWGSNMPA-UHFFFAOYSA-N  17397-85-2     CID:28516
#> 5  KWILGNNWGSNMPA-UHFFFAOYSA-N  17397-85-2     CID:28516
#> 6  CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 7  CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 8  CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 9  CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 10 CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 11 CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 12 CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 13 CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 14 CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 15 CTXQVLLVFBNZKL-YVEDVMJTSA-N 149849-90-1 CID:102004382
#> 16 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 17 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 18 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 19 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 20 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 21 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 22 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 23 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 24 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 25 OQIQSTLJSLGHID-WNWIJWBNSA-N   1162-65-8    CID:186907
#> 26 WWSYXEZEXMQWHT-WNWIJWBNSA-N   7220-81-7   CID:2724360
#> 27 WWSYXEZEXMQWHT-WNWIJWBNSA-N   7220-81-7   CID:2724360
#> 28 WWSYXEZEXMQWHT-WNWIJWBNSA-N   7220-81-7   CID:2724360
#> 29 WWSYXEZEXMQWHT-WNWIJWBNSA-N   7220-81-7   CID:2724360
#> 30 WWSYXEZEXMQWHT-WNWIJWBNSA-N   7220-81-7   CID:2724360
#> 31 WWSYXEZEXMQWHT-WNWIJWBNSA-N   7220-81-7   CID:2724360
#> 32 XWIYFDMXXLINPU-WNWIJWBNSA-N   1165-39-5   CID:2724361
#> 33 XWIYFDMXXLINPU-WNWIJWBNSA-N   1165-39-5   CID:2724361
#> 34 XWIYFDMXXLINPU-WNWIJWBNSA-N   1165-39-5   CID:2724361
#> 35 XWIYFDMXXLINPU-WNWIJWBNSA-N   1165-39-5   CID:2724361
#> 36 XWIYFDMXXLINPU-WNWIJWBNSA-N   1165-39-5   CID:2724361
#> 37 XWIYFDMXXLINPU-WNWIJWBNSA-N   1165-39-5   CID:2724361
#> 38 WPCVRWVBBXIRMA-WNWIJWBNSA-N   7241-98-7   CID:2724362
#> 39 WPCVRWVBBXIRMA-WNWIJWBNSA-N   7241-98-7   CID:2724362
#> 40 WPCVRWVBBXIRMA-WNWIJWBNSA-N   7241-98-7   CID:2724362
#> 41 WPCVRWVBBXIRMA-WNWIJWBNSA-N   7241-98-7   CID:2724362
#> 42 WPCVRWVBBXIRMA-WNWIJWBNSA-N   7241-98-7   CID:2724362
#> 43 WPCVRWVBBXIRMA-WNWIJWBNSA-N   7241-98-7   CID:2724362
#> 44 MJBWDEQAUQTVKK-IAGOWNOFSA-N   6795-23-9  CID:15558498
#> 45 MJBWDEQAUQTVKK-IAGOWNOFSA-N   6795-23-9  CID:15558498
#> 46 MJBWDEQAUQTVKK-IAGOWNOFSA-N   6795-23-9  CID:15558498
#> 47 MJBWDEQAUQTVKK-IAGOWNOFSA-N   6795-23-9  CID:15558498
#> 48 MJBWDEQAUQTVKK-IAGOWNOFSA-N   6795-23-9  CID:15558498
#> 49 MJBWDEQAUQTVKK-IAGOWNOFSA-N   6795-23-9  CID:15558498
#> 50 SZINUGQCTHLQAZ-DQYPLSBCSA-N  18172-33-3  CID:54682463
#> 51 SZINUGQCTHLQAZ-DQYPLSBCSA-N  18172-33-3  CID:54682463
#> 52 SZINUGQCTHLQAZ-DQYPLSBCSA-N  18172-33-3  CID:54682463
#> 53 SZINUGQCTHLQAZ-DQYPLSBCSA-N  18172-33-3  CID:54682463
#> 54 SZINUGQCTHLQAZ-DQYPLSBCSA-N  18172-33-3  CID:54682463
#> 55 MMHTXEATDNFMMY-WBIUFABUSA-N  29752-43-0     CID:34687
#> 56 MMHTXEATDNFMMY-WBIUFABUSA-N  29752-43-0     CID:34687
#> 57 MMHTXEATDNFMMY-WBIUFABUSA-N  29752-43-0     CID:34687
#> 58 MMHTXEATDNFMMY-WBIUFABUSA-N  29752-43-0     CID:34687
#> 59 MMHTXEATDNFMMY-WBIUFABUSA-N  29752-43-0     CID:34687
#> 60 CEBXXEKPIIDJHL-UHFFFAOYSA-N    641-38-3   CID:5359485
#> 61 CEBXXEKPIIDJHL-UHFFFAOYSA-N    641-38-3   CID:5359485
#> 62 CEBXXEKPIIDJHL-UHFFFAOYSA-N    641-38-3   CID:5359485
#> 63 LCSDQFNUYFTXMT-UHFFFAOYSA-N  23452-05-3   CID:5360741
#> 64 LCSDQFNUYFTXMT-UHFFFAOYSA-N  23452-05-3   CID:5360741
#> 65 LCSDQFNUYFTXMT-UHFFFAOYSA-N  23452-05-3   CID:5360741
#> 66 VSMBLBOUQJNJIL-JJXSEGSLSA-N  22268-16-2     CID:89644
#> 67 VSMBLBOUQJNJIL-JJXSEGSLSA-N  22268-16-2     CID:89644
#> 68 VSMBLBOUQJNJIL-JJXSEGSLSA-N  22268-16-2     CID:89644
#> 69 VSMBLBOUQJNJIL-JJXSEGSLSA-N  22268-16-2     CID:89644
#> 70 VSMBLBOUQJNJIL-JJXSEGSLSA-N  22268-16-2     CID:89644
#>                        name
#> 1                   Mellein
#> 2                   Mellein
#> 3                   Mellein
#> 4                   Mellein
#> 5                   Mellein
#> 6              AAL toxin TB
#> 7              AAL toxin TB
#> 8              AAL toxin TB
#> 9              AAL toxin TB
#> 10             AAL toxin TB
#> 11             AAL toxin TB
#> 12             AAL toxin TB
#> 13             AAL toxin TB
#> 14             AAL toxin TB
#> 15             AAL toxin TB
#> 16             Aflatoxin B1
#> 17             Aflatoxin B1
#> 18             Aflatoxin B1
#> 19             Aflatoxin B1
#> 20             Aflatoxin B1
#> 21             Aflatoxin B1
#> 22             Aflatoxin B1
#> 23             Aflatoxin B1
#> 24             Aflatoxin B1
#> 25             Aflatoxin B1
#> 26             Aflatoxin B2
#> 27             Aflatoxin B2
#> 28             Aflatoxin B2
#> 29             Aflatoxin B2
#> 30             Aflatoxin B2
#> 31             Aflatoxin B2
#> 32             Aflatoxin G1
#> 33             Aflatoxin G1
#> 34             Aflatoxin G1
#> 35             Aflatoxin G1
#> 36             Aflatoxin G1
#> 37             Aflatoxin G1
#> 38             Aflatoxin G2
#> 39             Aflatoxin G2
#> 40             Aflatoxin G2
#> 41             Aflatoxin G2
#> 42             Aflatoxin G2
#> 43             Aflatoxin G2
#> 44             Aflatoxin M1
#> 45             Aflatoxin M1
#> 46             Aflatoxin M1
#> 47             Aflatoxin M1
#> 48             Aflatoxin M1
#> 49             Aflatoxin M1
#> 50 alpha-Cyclopiazonic acid
#> 51 alpha-Cyclopiazonic acid
#> 52 alpha-Cyclopiazonic acid
#> 53 alpha-Cyclopiazonic acid
#> 54 alpha-Cyclopiazonic acid
#> 55                Altenuene
#> 56                Altenuene
#> 57                Altenuene
#> 58                Altenuene
#> 59                Altenuene
#> 60              Alternariol
#> 61              Alternariol
#> 62              Alternariol
#> 63 Alternariol methyl ether
#> 64 Alternariol methyl ether
#> 65 Alternariol methyl ether
#> 66           Altersolanol A
#> 67           Altersolanol A
#> 68           Altersolanol A
#> 69           Altersolanol A
#> 70           Altersolanol A

## Add also the synonyms (aliases) for the compounds. This will cause the
## tables compound and synonym to be joined. The elements of the compound_id
## and name are now no longer unique
res <- compounds(cdb, columns = c("name", "synonym"))
head(res)
#>           name
#> 1      Mellein
#> 2      Mellein
#> 3      Mellein
#> 4 AAL toxin TB
#> 5 AAL toxin TB
#> 6 Aflatoxin B1
#>                                                                                                               synonym
#> 1                                                                      8-hydroxy-3-methyl-3,4-dihydroisochromen-1-one
#> 2                                                                                                             Mellein
#> 3                                                                                                            Ochracin
#> 4 2-[2-[(3R,4R,5S,7S,14R,16S)-17-amino-4,14,16-trihydroxy-3,7-dimethylheptadecan-5-yl]oxy-2-oxoethyl]butanedioic acid
#> 5                                                                                                        AAL toxin TB
#> 6                  (6aR,9aS)-4-Methoxy-2,3,6a,9a-tetrahydrocyclopenta[c]furo[3',2':4,5]furo[2,3-h]chromene-1,11-dione

## List all database tables and their columns
tables(cdb)
#> $ms_compound
#> [1] "compound_id" "formula"     "exactmass"   "smiles"      "inchi"      
#> [6] "inchikey"    "cas"         "pubchem"     "name"       
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

## Any of these columns can be used in the `compounds` call to retrieve
## the specific annotations. The corresponding database tables will then be
## joined together
compounds(cdb, columns = c("formula", "publication"))
#>       formula
#> 1    C10H10O3
#> 2   C25H47NO9
#> 3    C17H12O6
#> 4    C17H14O6
#> 5    C17H12O7
#> 6    C17H14O7
#> 7  C20H20N2O3
#> 8    C15H16O6
#> 9    C14H10O5
#> 10   C15H12O5
#> 11   C16H16O8
#>                                                                                                                                                                                                                                                                              publication
#> 1  Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5
#> 2                  Renaud, J. B.; Kelman, M. J.; Qi, T. F.; Seifert, K. A.; Sumarah, M. W. Product Ion Filtering with Rapid Polarity Switching for the Detection of All Fumonisins and AAL-Toxins. Rapid Communications in Mass Spectrometry 2015, 29 (22), 2131–9. DOI:10.1002/rcm.7374
#> 3  Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5
#> 4  Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5
#> 5  Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5
#> 6  Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5
#> 7  Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5
#> 8                                                                                                                                                                                                                                                                                   <NA>
#> 9  Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5
#> 10 Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5
#> 11 Renaud, J. B.; Sumarah, M. W. Data Independent Acquisition-Digital Archiving Mass Spectrometry: Application to Single Kernel Mycotoxin Analysis of Fusarium Graminearum Infected Maize. Analytical and Bioanalytical Chemistry 2016, 408 (12), 3083–91. DOI:10.1007/s00216-016-9391-5

## Calculating m/z values for the exact masses of unique chemical formulas
## in the database:
mass2mz(cdb, adduct = c("[M+H]+", "[M+Na]+"))
#>              [M+H]+  [M+Na]+
#> C10H10O3   179.0703 201.0522
#> C25H47NO9  506.3324 528.3143
#> C17H12O6   313.0706 335.0526
#> C17H14O6   315.0863 337.0682
#> C17H12O7   329.0656 351.0475
#> C17H14O7   331.0812 353.0632
#> C20H20N2O3 337.1547 359.1366
#> C15H16O6   293.1020 315.0839
#> C14H10O5   259.0601 281.0420
#> C15H12O5   273.0757 295.0577
#> C16H16O8   337.0918 359.0737

## By using `name = "compound_id"` the calculation will be performed for
## each unique compound ID instead (resulting in potentially redundant
## results)
mass2mz(cdb, adduct = c("[M+H]+", "[M+Na]+"), name = "compound_id")
#>      [M+H]+  [M+Na]+
#> 1  179.0703 201.0522
#> 2  179.0703 201.0522
#> 3  179.0703 201.0522
#> 4  179.0703 201.0522
#> 5  179.0703 201.0522
#> 6  506.3324 528.3143
#> 7  506.3324 528.3143
#> 8  506.3324 528.3143
#> 9  506.3324 528.3143
#> 10 506.3324 528.3143
#> 11 506.3324 528.3143
#> 12 506.3324 528.3143
#> 13 506.3324 528.3143
#> 14 506.3324 528.3143
#> 15 506.3324 528.3143
#> 16 313.0706 335.0526
#> 17 313.0706 335.0526
#> 18 313.0706 335.0526
#> 19 313.0706 335.0526
#> 20 313.0706 335.0526
#> 21 313.0706 335.0526
#> 22 313.0706 335.0526
#> 23 313.0706 335.0526
#> 24 313.0706 335.0526
#> 25 313.0706 335.0526
#> 26 315.0863 337.0682
#> 27 315.0863 337.0682
#> 28 315.0863 337.0682
#> 29 315.0863 337.0682
#> 30 315.0863 337.0682
#> 31 315.0863 337.0682
#> 32 329.0656 351.0475
#> 33 329.0656 351.0475
#> 34 329.0656 351.0475
#> 35 329.0656 351.0475
#> 36 329.0656 351.0475
#> 37 329.0656 351.0475
#> 38 331.0812 353.0632
#> 39 331.0812 353.0632
#> 40 331.0812 353.0632
#> 41 331.0812 353.0632
#> 42 331.0812 353.0632
#> 43 331.0812 353.0632
#> 44 329.0656 351.0475
#> 45 329.0656 351.0475
#> 46 329.0656 351.0475
#> 47 329.0656 351.0475
#> 48 329.0656 351.0475
#> 49 329.0656 351.0475
#> 50 337.1547 359.1366
#> 51 337.1547 359.1366
#> 52 337.1547 359.1366
#> 53 337.1547 359.1366
#> 54 337.1547 359.1366
#> 55 293.1020 315.0839
#> 56 293.1020 315.0839
#> 57 293.1020 315.0839
#> 58 293.1020 315.0839
#> 59 293.1020 315.0839
#> 60 259.0601 281.0420
#> 61 259.0601 281.0420
#> 62 259.0601 281.0420
#> 63 273.0757 295.0577
#> 64 273.0757 295.0577
#> 65 273.0757 295.0577
#> 66 337.0918 359.0737
#> 67 337.0918 359.0737
#> 68 337.0918 359.0737
#> 69 337.0918 359.0737
#> 70 337.0918 359.0737

## Create a Spectra object with all MS/MS spectra from the database.
library(Spectra)
#> Loading required package: BiocParallel
sps <- Spectra(cdb)
sps
#> MSn data (Spectra) with 70 spectra in a MsBackendCompDb backend:
#>       msLevel precursorMz  polarity
#>     <integer>   <numeric> <integer>
#> 1           2      179.07         1
#> 2           2      179.07         1
#> 3           2      179.07         1
#> 4           2      179.07         1
#> 5           2      179.07         1
#> ...       ...         ...       ...
#> 66          2     337.091         1
#> 67          2     337.091         1
#> 68          2     337.091         1
#> 69          2     337.091         1
#> 70          2     337.091         1
#>  ... 47 more variables/columns.
#>  Use  'spectraVariables' to list all of them.
#>  data source: MassBank 
#>  version: 2020.09 
#>  organism: NA 

## Extract spectra for a specific compound.
sps <- Spectra(cdb, filter = ~ name == "Mellein")
sps
#> MSn data (Spectra) with 5 spectra in a MsBackendCompDb backend:
#>     msLevel precursorMz  polarity
#>   <integer>   <numeric> <integer>
#> 1         2      179.07         1
#> 2         2      179.07         1
#> 3         2      179.07         1
#> 4         2      179.07         1
#> 5         2      179.07         1
#>  ... 47 more variables/columns.
#>  Use  'spectraVariables' to list all of them.
#>  data source: MassBank 
#>  version: 2020.09 
#>  organism: NA 

## List all available annotations for MS/MS spectra
spectraVariables(sps)
#>  [1] "msLevel"                 "rtime"                  
#>  [3] "acquisitionNum"          "scanIndex"              
#>  [5] "dataStorage"             "dataOrigin"             
#>  [7] "centroided"              "smoothed"               
#>  [9] "polarity"                "precScanNum"            
#> [11] "precursorMz"             "precursorIntensity"     
#> [13] "precursorCharge"         "collisionEnergy"        
#> [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
#> [17] "isolationWindowUpperMz"  "compound_id"            
#> [19] "name"                    "formula"                
#> [21] "exactmass"               "smiles"                 
#> [23] "inchi"                   "inchikey"               
#> [25] "cas"                     "pubchem"                
#> [27] "accession"               "spectrum_name"          
#> [29] "date"                    "authors"                
#> [31] "license"                 "copyright"              
#> [33] "publication"             "splash"                 
#> [35] "precursor_intensity"     "adduct"                 
#> [37] "ionization"              "ionization_voltage"     
#> [39] "fragmentation_mode"      "collisionEnergy_text"   
#> [41] "instrument"              "instrument_type"        
#> [43] "precursorMz_text"        "spectrum_id"            
#> [45] "predicted"               "msms_mz_range_min"      
#> [47] "msms_mz_range_max"       "synonym"                

## Get access to the m/z values of these
mz(sps)
#> NumericList of length 5
#> [[1]] 133.0648 151.0754 155.9743 161.0597 179.0703
#> [[2]] 133.0648 151.0754 155.9745 161.0597 179.0703
#> [[3]] 105.0699 133.0648 151.0754 161.0597 179.0703
#> [[4]] 105.0699 133.0648 151.0754 161.0597 179.0703
#> [[5]] 105.0699 115.0542 133.0648 151.0754 161.0597 179.0703

library(Spectra)
## Plot the first spectrum
plotSpectra(sps[1])



#########
## Filtering the database
##
## Get all compounds with an exact mass between 310 and 320
res <- compounds(cdb, filter = ~ exactmass > 310 & exactmass < 320)
res
#>    formula exactmass                                                   smiles
#> 1 C17H12O6  312.0634 COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1
#> 2 C17H14O6  314.0790  COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#>                                                                                                                            inchi
#> 1 InChI=1S/C17H12O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h4-6,8,17H,2-3H2,1H3/t8-,17+/m0/s1
#> 2   InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#>                      inchikey       cas     pubchem         name
#> 1 OQIQSTLJSLGHID-WNWIJWBNSA-N 1162-65-8  CID:186907 Aflatoxin B1
#> 2 WWSYXEZEXMQWHT-WNWIJWBNSA-N 7220-81-7 CID:2724360 Aflatoxin B2

## Get all compounds that have an H14 in their formula.
res <- compounds(cdb, filter = FormulaFilter("H14", "contains"))
res
#>    formula exactmass                                                   smiles
#> 1 C17H14O6  314.0790  COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#> 2 C17H14O7  330.0739 COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1
#>                                                                                                                            inchi
#> 1   InChI=1S/C17H14O6/c1-20-10-6-11-14(8-4-5-21-17(8)22-11)15-13(10)7-2-3-9(18)12(7)16(19)23-15/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#> 2 InChI=1S/C17H14O7/c1-20-9-6-10-12(8-3-5-22-17(8)23-10)14-11(9)7-2-4-21-15(18)13(7)16(19)24-14/h6,8,17H,2-5H2,1H3/t8-,17+/m0/s1
#>                      inchikey       cas     pubchem         name
#> 1 WWSYXEZEXMQWHT-WNWIJWBNSA-N 7220-81-7 CID:2724360 Aflatoxin B2
#> 2 WPCVRWVBBXIRMA-WNWIJWBNSA-N 7241-98-7 CID:2724362 Aflatoxin G2

#########
## Using CompDb with the *tidyverse*
##
## Using return.type = "tibble" the result will be returned as a "tibble"
compounds(cdb, return.type = "tibble")
#> # A tibble: 12 × 8
#>    formula    exactmass smiles                inchi inchikey cas   pubchem name 
#>    <chr>          <dbl> <chr>                 <chr> <chr>    <chr> <chr>   <chr>
#>  1 C10H10O3        178. CC1CC2=C(C(=CC=C2)O)… InCh… KWILGNN… 1739… CID:28… Mell…
#>  2 C25H47NO9       505. CC[C@@H](C)[C@H]([C@… InCh… CTXQVLL… 1498… CID:10… AAL …
#>  3 C17H12O6        312. COC1=C2C3=C(C(=O)CC3… InCh… OQIQSTL… 1162… CID:18… Afla…
#>  4 C17H14O6        314. COC1=C2C3=C(C(=O)CC3… InCh… WWSYXEZ… 7220… CID:27… Afla…
#>  5 C17H12O7        328. COC1=C2C3=C(C(=O)OCC… InCh… XWIYFDM… 1165… CID:27… Afla…
#>  6 C17H14O7        330. COC1=C2C3=C(C(=O)OCC… InCh… WPCVRWV… 7241… CID:27… Afla…
#>  7 C17H12O7        328. COC1=C2C3=C(C(=O)CC3… InCh… MJBWDEQ… 6795… CID:15… Afla…
#>  8 C20H20N2O3      336. CC(=O)C1=C([C@@H]2[C… InCh… SZINUGQ… 1817… CID:54… alph…
#>  9 C15H16O6        292. C[C@]12C[C@@H]([C@H]… InCh… MMHTXEA… 2975… CID:34… Alte…
#> 10 C14H10O5        258. CC1=CC(=CC2=C1C3=CC(… InCh… CEBXXEK… 641-… CID:53… Alte…
#> 11 C15H12O5        272. CC1=CC(=CC2=C1C3=CC(… InCh… LCSDQFN… 2345… CID:53… Alte…
#> 12 C16H16O8        336. C[C@]1([C@@H]([C@H](… InCh… VSMBLBO… 2226… CID:89… Alte…

## Use the CompDb in a dplyr setup
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:S4Vectors’:
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from ‘package:BiocGenerics’:
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from ‘package:generics’:
#> 
#>     explain
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
src_cmp <- src_compdb(cdb)
src_cmp
#> src:  sqlite 3.51.1 [/__w/_temp/Library/CompoundDb/sql/CompDb.MassBank.sql]
#> tbls: metadata, ms_compound, msms_spectrum, msms_spectrum_peak, synonym

## Get a tbl for the ms_compound table
cmp_tbl <- tbl(src_cmp, "ms_compound")

## Extract the id, name and inchi
cmp_tbl %>% select(compound_id, name, inchi) %>% collect()
#> # A tibble: 70 × 3
#>    compound_id name         inchi                                               
#>    <chr>       <chr>        <chr>                                               
#>  1 1           Mellein      InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-…
#>  2 2           Mellein      InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-…
#>  3 3           Mellein      InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-…
#>  4 4           Mellein      InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-…
#>  5 5           Mellein      InChI=1S/C10H10O3/c1-6-5-7-3-2-4-8(11)9(7)10(12)13-…
#>  6 6           AAL toxin TB InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-1…
#>  7 7           AAL toxin TB InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-1…
#>  8 8           AAL toxin TB InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-1…
#>  9 9           AAL toxin TB InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-1…
#> 10 10          AAL toxin TB InChI=1S/C25H47NO9/c1-4-17(3)24(32)21(35-23(31)13-1…
#> # ℹ 60 more rows

########
## Creating an empty CompDb and sequentially adding content
##
## Create an empty CompDb and store the database in a temporary file
cdb <- emptyCompDb(tempfile())
cdb
#> class: CompDb 
#>  data source: NA 
#>  version: NA 
#>  organism: NA 
#>  compound count: 0 

## Define a data.frame with some compounds to add
cmp <- data.frame(
    compound_id = c(1, 2),
    name = c("Caffeine", "Glucose"),
    formula = c("C8H10N4O2", "C6H12O6"),
    exactmass = c(194.080375584, 180.063388116))

## We can also add multiple synonyms for each compound
cmp$synonyms <- list(c("Cafeina", "Koffein"), "D Glucose")
cmp
#>   compound_id     name   formula exactmass     synonyms
#> 1           1 Caffeine C8H10N4O2  194.0804 Cafeina,....
#> 2           2  Glucose   C6H12O6  180.0634    D Glucose

## These compounds can be added to the empty database with insertCompound
cdb <- insertCompound(cdb, compounds = cmp)
compounds(cdb)
#>       name inchi inchikey   formula exactmass
#> 1 Caffeine  <NA>     <NA> C8H10N4O2  194.0804
#> 2  Glucose  <NA>     <NA>   C6H12O6  180.0634

## insertCompound would also allow to add additional columns/annotations to
## the database. Below we define a new compound adding an additional column
## hmdb_id
cmp <- data.frame(
    compound_id = 3,
    name = "Alpha-Lactose",
    formula = "C12H22O11",
    exactmass = 342.116211546,
    hmdb_id = "HMDB0000186")

## To add additional columns we need to set addColumns = TRUE
cdb <- insertCompound(cdb, compounds = cmp, addColumns = TRUE)
cdb
#> class: CompDb 
#>  data source: NA 
#>  version: NA 
#>  organism: NA 
#>  compound count: 3 
compounds(cdb)
#>            name inchi inchikey   formula exactmass     hmdb_id
#> 1 Alpha-Lactose  <NA>     <NA> C12H22O11  342.1162 HMDB0000186
#> 2      Caffeine  <NA>     <NA> C8H10N4O2  194.0804        <NA>
#> 3       Glucose  <NA>     <NA>   C6H12O6  180.0634        <NA>

######
## Deleting selected compounds from a database
##
## Compounds can be deleted with the deleteCompound function providing the
## IDs of the compounds that should be deleted. IDs of compounds in the
## database can be retrieved by adding "compound_id" to the columns parameter
## of the compounds function:
compounds(cdb, columns = c("compound_id", "name"))
#>   compound_id          name
#> 1           1      Caffeine
#> 2           2       Glucose
#> 3           3 Alpha-Lactose

## Compounds can be deleted with the deleteCompound function. Below we delete
## the compounds with the IDs "1" and "3" from the database
cdb <- deleteCompound(cdb, ids = c("1", "3"))
compounds(cdb)
#>      name inchi inchikey formula exactmass hmdb_id
#> 1 Glucose  <NA>     <NA> C6H12O6  180.0634    <NA>

## If also MS2 spectra associated with any of these two compounds an error
## would be thrown. Setting the parameter `recursive = TRUE` in the
## `deleteCompound` call would delete the compounds along with their MS2
## spectra.
```

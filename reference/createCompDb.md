# Create a CompDb database

`CompDb` databases can be created with the `createCompDb()` or the
`emptyCompDb()` functions, the former creating and initializing
(filling) the database with existing data, the latter creating an empty
database that can be subsequently filled with
[`insertCompound()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
or
[`insertSpectra()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
calls.

`emptyCompDb()` requires only the file name of the database that should
be created as input and returns a `CompDb` representing the empty
database.

`createCompDb()` creates a `SQLite`-based
[`CompDb`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
object/database from a compound resource provided as a `data.frame` or
`tbl`. Alternatively, the name(s) of the file(s) from which the
annotation should be extracted can be provided. Supported are SDF files
(such as those from the *Human Metabolome Database* HMDB) that can be
read using the
[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
or LipidBlast files (see
[`compound_tbl_lipidblast()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_lipidblast.md).

An additional `data.frame` providing metadata information including the
data source, date, version and organism is mandatory. By default, the
function will define the name of the database based on the provided
metadata, but it is also possible to define this manually with the
`dbFile` parameter.

Optionally MS/MS (MS2) spectra for compounds can be also stored in the
database. Currently only MS/MS spectra from HMDB are supported. These
can be downloaded in XML format from HMDB (http://www.hmdb.ca), loaded
with the
[`msms_spectra_hmdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_hmdb.md)
or
[`msms_spectra_mona()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_mona.md)
function and passed to the function with the `msms_spectra` argument.
See
[`msms_spectra_hmdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_hmdb.md)
or
[`msms_spectra_mona()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_mona.md)
for information on the expected columns and format.

Required columns for the `data.frame` providing the compound information
( parameter `x`) are:

- `"compound_id"`: the ID of the compound. Can be an `integer` or
  `character`. Duplicated IDs are supported (for compatibility reasons),
  but not suggested. No missing values allowed.

- `"name"`: the compound's name.

- `"inchi"`: the InChI of the compound.

- `"inchikey"`: the InChI key.

- `"formula"`: the chemical formula.

- `"exactmass"`: the compound's (exact) mass.

- `"synonyms"`: additional synonyms/aliases for the compound. Should be
  either a single character or a list of values for each compound.

Any additional columns in the provided `data.frame` (such as e.g.
`"smiles"` providing the compound's SMILES) are also supported and will
be inserted into the database table.

See e.g.
[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
or
[`compound_tbl_lipidblast()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_lipidblast.md)
for functions creating such compound tables.

The table containing the MS2 spectra data should have the following
format and columns:

- `"spectrum_id"`: an arbitrary ID for the spectrum. Has to be an
  `integer`.

- `"compound_id"`: the ID of the compound to which the spectrum can be
  associated with. This has to be present in the `data.frame` defining
  the compounds.

- `"polarity"`: the polarity (as an `integer`, `0` for negative, `1` for
  positive, `NA` for not set).

- `"collision_energy"`: the collision energy.

- `"predicted"`: whether the spectrum was predicted or measured.

- `"splash"`: the SPLASH of the spectrum.

- `"instrument_type"`: the instrument type.

- `"instrument"`: the name of the instrument.

- `"precursor_mz"`: the precursor m/z (as a `numeric`).

- `"mz"`: the m/z values.

- `"intensity"`: the intensity values.

Only for columns `"spectrum_id"`, `"compound_id"`, `"mz"` and
`"intensity"` a value has to be provided in each row of the
`data.frame`. The others are optional. Note that the `data.frame` can be
either in the format as in the example below (i.e. each row being one
spectrum and columns `"mz"` and `"intensity"` being of type `list` each
element being the m/z or intensity values of one spectrum) or in a
*full* form, in which each row represents one *peak* and all columns
except `"mz"` and `"intensity"` containing redundant information of each
spectrum (hence columns `"mz"` and `"intensity"` being of type
`numeric`).

The metadata `data.frame` is supposed to have two columns named `"name"`
and `"value"` providing the following minimal information as key-value
pairs (see `make_metadata` for a utility function to create such a
`data.frame`):

- `"source"`: the source from which the data was retrieved (e.g.
  `"HMDB"`).

- `"url"`: the url from which the original data was retrieved.

- `"source_version"`: the version from the original data source (e.g.
  `"v4"`).

- `"source_date"`: the date when the original data source was generated.

- `"organism"`: the organism. Should be in the form `"Hsapiens"` or
  `"Mmusculus"`.

`createCompDbPackage` creates an R data package with the data from a
[`CompDb`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
object.

`make_metadata()` helps generating a metadata `data.frame` in the
correct format expected by the `createCompDb` function. The function
returns a `data.frame`.

## Usage

``` r
createCompDb(x, metadata, msms_spectra, path = ".", dbFile = character())

createCompDbPackage(
  x,
  version,
  maintainer,
  author,
  path = ".",
  license = "Artistic-2.0"
)

make_metadata(
  source = character(),
  url = character(),
  source_version = character(),
  source_date = character(),
  organism = NA_character_
)

emptyCompDb(dbFile = character())
```

## Arguments

- x:

  For `createCompDb()`: `data.frame` or `tbl` with the compound
  annotations or `character` with the file name(s) from which the
  compound annotations should be retrieved. See description for details.

      For `createCompDbPackage()`: `character(1)` with the file name of the
      `CompDb` SQLite file (created by `createCompDb`).

- metadata:

  For `createCompDb()`: `data.frame` with metadata information. See
  description for details.

- msms_spectra:

  For `createCompDb()`: `data.frame` with MS/MS spectrum data. See
  [`msms_spectra_hmdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_hmdb.md)
  for the expected format and a function to import such data from
  spectrum xml files from HMDB.

- path:

  `character(1)` with the path to the directory where the database file
  or package folder should be written. Defaults to the current
  directory.

- dbFile:

  `character(1)` to optionally provide the name of the SQLite database
  file. If not provided (the default) the database name is defined using
  information from the provided `metadata`.

- version:

  For `createCompDbPackage()`: `character(1)` with the version of the
  package (ideally in the format `"x.y.z"`).

- maintainer:

  For `createCompDbPackage()`: `character(1)` with the name and email
  address of the package maintainer (in the form
  `"First Last <first.last@provider.com>"`.

- author:

  For `createCompDbPackage()`: `character(1)` with the name of the
  package author.

- license:

  For `createCompDbPackage()`: `character(1)` with the license of the
  package respectively the originating provider.

- source:

  For `make_metadata()`: `character(1)` with the name of the resource
  that provided the compound annotation.

- url:

  For `make_metadata()`: `character(1)` with the url to the original
  resource.

- source_version:

  For `make_metadata()`: `character(1)` with the version of the original
  resource providing the annotation.

- source_date:

  For `make_metadata()`: `character(1)` with the date of the resource's
  release.

- organism:

  For `make_metadata()`: `character(1)` with the name of the organism.
  This should be in the format `"Hsapiens"` for human, `"Mmusculus"` for
  mouse etc. Leave to `NA` if not applicable.

## Value

For `createCompDb()`: a `character(1)` with the database name
(invisibly).

## Details

Metadata information is also used to create the file name for the
database file. The name starts with `"CompDb"`, followed by the
organism, the data source and its version. A compound database file for
HMDB version 4 with human metabolites will thus be named:
`"CompDb.Hsapiens.HMDB.v4"`.

A single `CompDb` database is created from multiple SDF files (e.g. for
*PubChem*) if all the file names are provided with parameter `x`.
Parallel processing is currently not enabled because SQLite does not
support it yet natively.

## See also

[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
and
[`compound_tbl_lipidblast()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_lipidblast.md)
for functions to extract compound annotations from files in SDF format,
or files from LipidBlast.

[`import_mona_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/import_mona_sdf.md)
to import both the compound and spectrum data from a SDF file from MoNa
(Massbank of North America) in one call.

[`msms_spectra_hmdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_hmdb.md)
and
[`msms_spectra_mona()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_mona.md)
for functions to import MS/MS spectrum data from xml files from HMDB or
an SDF file from MoNa.

[`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
for how to use a compound database.

## Author

Johannes Rainer

## Examples

``` r

## Read compounds for a HMDB subset
fl <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
cmps <- compound_tbl_sdf(fl)

## Create a metadata data.frame for the compounds.
metad <- data.frame(name = c("source", "url", "source_version",
    "source_date", "organism"), value = c("HMDB", "http://www.hmdb.ca",
    "v4", "2017-08-27", "Hsapiens"))

## Alternatively use the make_metadata helper function
metad <- make_metadata(source = "HMDB", source_version = "v4",
    source_date = "2017-08", organism = "Hsapiens",
    url = "http://www.hmdb.ca")
## Create a SQLite database in the temporary folder
db_f <- createCompDb(cmps, metadata = metad, path = tempdir())

## The database can be loaded and accessed with a CompDb object
db <- CompDb(db_f)
db
#> class: CompDb 
#>  data source: HMDB 
#>  version: v4 
#>  organism: Hsapiens 
#>  compound count: 9 

## Create a database for HMDB that includes also MS/MS spectrum data
metad2 <- make_metadata(source = "HMDB_with_spectra", source_version = "v4",
    source_date = "2017-08", organism = "Hsapiens",
    url = "http://www.hmdb.ca")

## Import spectrum information from selected MS/MS xml files from HMDB
## that are provided in the package
xml_path <- system.file("xml", package = "CompoundDb")
spctra <- msms_spectra_hmdb(xml_path)
#> Going to process 4 xml files.
#> Postprocessing data ... 
#> OK

## Create a SQLite database in the temporary folder
db_f2 <- createCompDb(cmps, metadata = metad2, msms_spectra = spctra,
    path = tempdir())

## The database can be loaded and accessed with a CompDb object
db2 <- CompDb(db_f2)
db2
#> class: CompDb 
#>  data source: HMDB_with_spectra 
#>  version: v4 
#>  organism: Hsapiens 
#>  compound count: 9 
#>  MS/MS spectra count: 4 

## Does the database contain MS/MS spectrum data?
hasMsMsSpectra(db2)
#> [1] TRUE

## Create a database for a ChEBI subset providing the file name of the
## corresponding SDF file
metad <- make_metadata(source = "ChEBI_sub", source_version = "2",
    source_date = NA, organism = "Hsapiens", url = "www.ebi.ac.uk/chebi")
db_f <- createCompDb(system.file("sdf/ChEBI_sub.sdf.gz",
    package = "CompoundDb"), metadata = metad, path = tempdir())
#> Import data from ChEBI_sub.sdf.gz ...
#> OK
db <- CompDb(db_f)
db
#> class: CompDb 
#>  data source: ChEBI_sub 
#>  version: 2 
#>  organism: Hsapiens 
#>  compound count: 6 

compounds(db)
#>                       name
#> 1          (-)-epicatechin
#> 2         (1S,4R)-fenchone
#> 3   1-alkyl-2-acylglycerol
#> 4   16alpha-hydroxyestrone
#> 5 2,6-dichlorobenzonitrile
#> 6    2-hydroxybutyric acid
#>                                                                                                                                                 inchi
#> 1                                    InChI=1S/C15H14O6/c16-8-4-11(18)9-6-13(20)15(21-14(9)5-8)7-1-2-10(17)12(19)3-7/h1-5,13,15-20H,6H2/t13-,15-/m1/s1
#> 2                                                                         InChI=1S/C10H16O/c1-9(2)7-4-5-10(3,6-7)8(9)11/h7H,4-6H2,1-3H3/t7-,10+/m1/s1
#> 3                                                                                                                                                <NA>
#> 4 InChI=1S/C18H22O3/c1-18-7-6-13-12-5-3-11(19)8-10(12)2-4-14(13)15(18)9-16(20)17(18)21/h3,5,8,13-16,19-20H,2,4,6-7,9H2,1H3/t13-,14-,15+,16-,18+/m1/s1
#> 5                                                                                                     InChI=1S/C7H3Cl2N/c8-6-2-1-3-7(9)5(6)4-10/h1-3H
#> 6                                                                                                InChI=1S/C4H8O3/c1-2-3(5)4(6)7/h3,5H,2H2,1H3,(H,6,7)
#>                      inchikey  formula exactmass
#> 1 PFTAWBLQPZVEMU-UKRRQHHQSA-N C15H14O6   290.079
#> 2 LHXDLQBQYFFVNW-XCBNKYQSSA-N  C10H16O   152.120
#> 3                        <NA> C4H6O4R2   118.027
#> 4 WPOCIZJTELRQMF-QFXBJFAPSA-N C18H22O3   286.157
#> 5 YOYAIZYFCNQIRF-UHFFFAOYSA-N C7H3Cl2N   170.964
#> 6 AFENDNXGAFYKQO-UHFFFAOYSA-N   C4H8O3   104.047
#>                                                                 smiles
#> 1                     [H][C@@]1(Oc2cc(O)cc(O)c2C[C@H]1O)c1ccc(O)c(O)c1
#> 2                                      CC1(C)[C@@H]2CC[C@@](C)(C2)C1=O
#> 3                                                  OCC(CO[*])OC([*])=O
#> 4 [H][C@]12CC[C@]3(C)C(=O)[C@H](O)C[C@@]3([H])[C@]1([H])CCc1cc(O)ccc21
#> 5                                                    Clc1cccc(Cl)c1C#N
#> 6                                                         CCC(O)C(O)=O

## connect to the database and query it's tables using RSQlite
library(RSQLite)
con <- dbConnect(dbDriver("SQLite"), db_f)

dbGetQuery(con, "select * from metadata")
#>                 name                    value
#> 1             source                ChEBI_sub
#> 2                url      www.ebi.ac.uk/chebi
#> 3     source_version                        2
#> 4        source_date                     <NA>
#> 5           organism                 Hsapiens
#> 6   db_creation_date Mon Jan 12 14:58:38 2026
#> 7 supporting_package               CompoundDb
#> 8  supporting_object                   CompDb
dbGetQuery(con, "select * from ms_compound")
#>   compound_id                     name
#> 1    CHEBI:90          (-)-epicatechin
#> 2   CHEBI:165         (1S,4R)-fenchone
#> 3   CHEBI:598   1-alkyl-2-acylglycerol
#> 4   CHEBI:776   16alpha-hydroxyestrone
#> 5   CHEBI:943 2,6-dichlorobenzonitrile
#> 6  CHEBI:1148    2-hydroxybutyric acid
#>                                                                                                                                                 inchi
#> 1                                    InChI=1S/C15H14O6/c16-8-4-11(18)9-6-13(20)15(21-14(9)5-8)7-1-2-10(17)12(19)3-7/h1-5,13,15-20H,6H2/t13-,15-/m1/s1
#> 2                                                                         InChI=1S/C10H16O/c1-9(2)7-4-5-10(3,6-7)8(9)11/h7H,4-6H2,1-3H3/t7-,10+/m1/s1
#> 3                                                                                                                                                <NA>
#> 4 InChI=1S/C18H22O3/c1-18-7-6-13-12-5-3-11(19)8-10(12)2-4-14(13)15(18)9-16(20)17(18)21/h3,5,8,13-16,19-20H,2,4,6-7,9H2,1H3/t13-,14-,15+,16-,18+/m1/s1
#> 5                                                                                                     InChI=1S/C7H3Cl2N/c8-6-2-1-3-7(9)5(6)4-10/h1-3H
#> 6                                                                                                InChI=1S/C4H8O3/c1-2-3(5)4(6)7/h3,5H,2H2,1H3,(H,6,7)
#>                      inchikey  formula exactmass
#> 1 PFTAWBLQPZVEMU-UKRRQHHQSA-N C15H14O6   290.079
#> 2 LHXDLQBQYFFVNW-XCBNKYQSSA-N  C10H16O   152.120
#> 3                        <NA> C4H6O4R2   118.027
#> 4 WPOCIZJTELRQMF-QFXBJFAPSA-N C18H22O3   286.157
#> 5 YOYAIZYFCNQIRF-UHFFFAOYSA-N C7H3Cl2N   170.964
#> 6 AFENDNXGAFYKQO-UHFFFAOYSA-N   C4H8O3   104.047
#>                                                                 smiles
#> 1                     [H][C@@]1(Oc2cc(O)cc(O)c2C[C@H]1O)c1ccc(O)c(O)c1
#> 2                                      CC1(C)[C@@H]2CC[C@@](C)(C2)C1=O
#> 3                                                  OCC(CO[*])OC([*])=O
#> 4 [H][C@]12CC[C@]3(C)C(=O)[C@H](O)C[C@@]3([H])[C@]1([H])CCc1cc(O)ccc21
#> 5                                                    Clc1cccc(Cl)c1C#N
#> 6                                                         CCC(O)C(O)=O

## To create a CompDb R-package we could simply use the
## createCompDbPackage function on the SQLite database file name.
```

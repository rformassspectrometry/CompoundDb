# Expand a CompDb database with additional, related tables

The `CompDb` object uses a simple relational database model that
consists of the following database tables, some of which are optional:

- *ms_compound*: annotation(s) of compounds.

- *metadata*: general metadata information on the database. This
  database table is **not** related to any other table in the database
  and its content is thus also not joined with other database tables.

- *synonym* (optional): database table containing optional additional
  synonym(s) for compounds in the *ms_compound* table. Rows in this
  table are linked to a row in *ms_compound* through the `"compound_id"`
  database table column.

- *msms_spectrum* (optional): database table with information on
  individual mass spectra (each row containing the metadata for one
  spectrum). Database table column `"compound_id"` links entries in this
  database table to a single row in the *ms_compound* table.

- *msms_spectrum_peak* (otional): database table containing mass peak
  data. Each row in this table is related to one row in the
  *msms_spectrum* table (through the `"spectrum_id"` column present in
  both tables).

In addition, the `CompDb` database layout can be extended by adding
additional tables. To make their content automatically available through
the built-in
[`compounds()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
or [`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
functions, the information on how to combine/join these tables with the
existing ones needs to be provided. This can be done using the
`addJoinDefinition()` function: the relationship of a new table with one
of the existing tables can be defined with this function providing the
names of the two database tables as well as the names of the columns
containing the primary/foreign keys defining the the relationship.

See the section *Extending CompDb databases* in the *Creating CompoundDb
annotation resources* package vignette for a detailed example.

## Usage

``` r
addJoinDefinition(
  x,
  table_a = character(),
  table_b = character(),
  column_a = character(),
  column_b = character(),
  join = "left outer join"
)
```

## Arguments

- x:

  `CompDb` to which the join definition should be added.

- table_a:

  `character(1)` with the name of one of the two tables that are related
  to each other (and can be joined).

- table_b:

  `character(1)` with the name of the second of the two tables that are
  related to each other (and can be joined).

- column_a:

  `character(1)` with the name of the column in `table_a` containing the
  keys for the relationship to table `table_b`.

- column_b:

  `character(1)` with the name of the column in `table_b` containing the
  keys for the relationship to table `table_a`.

- join:

  `character(1)` with the type of join. Defaults to
  `join = "left outer join"`.

## Value

The input `CompDb` with tha added information on how to join the
respective database tables.

## Author

Johannes Rainer

## Examples

``` r

## The pre-defined table join definitions:
CompoundDb:::.JOINS
#>      [,1]            [,2]                
#> [1,] "ms_compound"   "synonym"           
#> [2,] "ms_ion"        "ms_compound"       
#> [3,] "ms_compound"   "msms_spectrum"     
#> [4,] "ms_ion"        "msms_spectrum"     
#> [5,] "msms_spectrum" "synonym"           
#> [6,] "msms_spectrum" "msms_spectrum_peak"
#>      [,3]                                                           
#> [1,] "on (ms_compound.compound_id=synonym.compound_id)"             
#> [2,] "on (ms_ion.compound_id=ms_compound.compound_id)"              
#> [3,] "on (ms_compound.compound_id=msms_spectrum.compound_id)"       
#> [4,] "on (ms_ion.compound_id=msms_spectrum.compound_id)"            
#> [5,] "on (msms_spectrum.compound_id=synonym.compound_id)"           
#> [6,] "on (msms_spectrum.spectrum_id=msms_spectrum_peak.spectrum_id)"
#>      [,4]             
#> [1,] "left outer join"
#> [2,] "left outer join"
#> [3,] "left outer join"
#> [4,] "left outer join"
#> [5,] "left outer join"
#> [6,] "left outer join"

## See section "Extending CompDb databases" in the *Creating CompoundDb
## annotation resources* package vignette for examples
```

# CompDb-based MS spectrum backend

The `MsBackendCompDb` *represents* MS2 spectra data from a
[`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
object/database. The object keeps only the primary keys of the spectra,
the associated compound IDs and the precursor m/z values in memory and
has thus only a very low memory footprint. All spectra variables,
including m/z and intensity values are retrieved from the database
*on-demand*. By extending the
[`Spectra::MsBackendCached()`](https://rdrr.io/pkg/Spectra/man/MsBackendCached.html)
class directly, `MsBackendCompDb` supports adding/replacing spectra
variables. These values are however only cached within the object and
not propagated (written) to the database.

It is not intended that users create or use instances of this class
directly, the
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
call on
[`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
will return a `Spectra` object that uses this backend.

The `MsBackendCompDb` does not support parallel processing because the
database connection stored within the object can not be used across
multiple parallel processes. The
[`backendBpparam()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method for `MsBackendCompDb` thus returns always `SerialParam` and hence
any function that uses this method to check for parallel processing
capability of a `MsBackend` will by default disable parallel processing.

## Usage

``` r
MsBackendCompDb()

# S4 method for class 'MsBackendCompDb'
backendInitialize(object, x, filter, ...)

# S4 method for class 'MsBackendCompDb'
show(object)

# S4 method for class 'MsBackendCompDb'
peaksData(object, columns = c("mz", "intensity"))

# S4 method for class 'MsBackendCompDb'
peaksVariables(object)

# S4 method for class 'MsBackendCompDb'
dataStorage(object)

# S4 method for class 'MsBackendCompDb'
intensity(object) <- value

# S4 method for class 'MsBackendCompDb'
mz(object) <- value

# S4 method for class 'MsBackendCompDb'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendCompDb'
spectraNames(object)

# S4 method for class 'MsBackendCompDb'
spectraNames(object) <- value

# S4 method for class 'MsBackendCompDb,ANY'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MsBackendCompDb,ANY'
extractByIndex(object, i)

# S4 method for class 'MsBackendCompDb'
x$name <- value

# S4 method for class 'MsBackendCompDb'
tic(object, initial = TRUE)

# S4 method for class 'MsBackendCompDb'
backendBpparam(object, BPPARAM = bpparam())
```

## Arguments

- object:

  an `MsBackendCompDb` instance.

- x:

  an `MsBackendCompDb` instance.

- filter:

  for
  [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  optional filter expression to specify which elements to retrieve from
  the database.

- ...:

  ignored.

- columns:

  for
  [`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `character` with names of columns/spectra variables that should be
  returned. Defaults to `spectraVariables(object)`. Database columns
  `"ms_level"`, `"precursor_mz"`, `"precursor_intensity"`,
  `"precursor_charge"` are mapped to the core `Spectra` variables
  `msLevel`, `precursorMz`, `precursorIntensity` and `precursorCharge`,
  respectively. For `peaksData`: `character` with the names of the peaks
  columns to return. Use `peaksVariables` for supported values.

- value:

  for `$<-`: the replacement values.

- i:

  For `[`: `integer`, `logical` or `character` to subset the object.

- j:

  For `[`: not supported.

- drop:

  For `[`: not considered.

- name:

  for `$<-`: the name of the spectra variable to replace.

- initial:

  for [`tic()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `logical(1)` whether original total ion current values should be
  returned or if the values should be calculated based on the actual
  intensity values of each spectrum.

- BPPARAM:

  for
  [`backendBpparam()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  `BiocParallel` parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

## Value

See the description of the respective function.

## Note

For higher performance it is suggested to change the backend of the
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
object to an
[`Spectra::MsBackendMemory()`](https://rdrr.io/pkg/Spectra/man/MsBackend.html)
backend with the
[`Spectra::setBackend()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
method of `Spectra` objects.

## Methods implemented for `MsBackendCompDb`

The methods listed here are implemented for the `MsBackendCompDb`. All
other methods are inherited directly from the parent
[`Spectra::MsBackendCached()`](https://rdrr.io/pkg/Spectra/man/MsBackendCached.html)
class. See the help of
[`Spectra::MsBackend()`](https://rdrr.io/pkg/Spectra/man/MsBackend.html)
in the `Spectra` package for a complete listing of methods.

- [`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html):
  gets the full list of peak matrices. Returns a
  [`list()`](https://rdrr.io/r/base/list.html), length equal to the
  number of spectra and each element being a `matrix` with columns
  `"mz"` and `"intensity"` with the spectra's m/z and intensity values.

- [`peaksVariables()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html):
  lists the available peaks variables in the backend (database). These
  can be used for parameter `columns` of
  [`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html).

- `intensity<-`: not supported.

- `mz<-`: not supported.

- [`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  returns the complete spectrum data including m/z and intensity values
  as a
  [`S4Vectors::DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html).

- `$<-`: replace or add a spectrum variable. Note that `mz`, `intensity`
  and `spectrum_id` variables can not be replaced.

- [`spectraNames()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  returns values from *spectrum_id* database column.

## Author

Johannes Rainer

## Examples

``` r

## MsBackendCompDb are not expected to be created/instanciated by users
## directly. Users also almost never directly interact with this type of
## object, as it is intented as a pure data backend for the `Spectra` object.
## Users will thus access MS data through such `Spectra` object, which can
## be created for `CompDb` objects using the `Spectra` method (see help
## of the `CompDb` object for more information. This examples shows how
## a `MsBackendCompDb` could be created purely from an SQLite database
## with data from a CompoundDb database.

## Connect to the SQLite database of a `CompDb` distributed via this package
library(RSQLite)
library(Spectra)
cdb <- CompDb(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"))

be <- backendInitialize(MsBackendCompDb(), cdb)
be
#> MsBackendCompDb with 70 spectra
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

## Accessing m/z values
mz(be)
#> NumericList of length 70
#> [[1]] 133.0648 151.0754 155.9743 161.0597 179.0703
#> [[2]] 133.0648 151.0754 155.9745 161.0597 179.0703
#> [[3]] 105.0699 133.0648 151.0754 161.0597 179.0703
#> [[4]] 105.0699 133.0648 151.0754 161.0597 179.0703
#> [[5]] 105.0699 115.0542 133.0648 151.0754 161.0597 179.0703
#> [[6]] 506.3324
#> [[7]] 276.2686 294.2792 312.287 330.2976 470.3112 488.3218 506.3324
#> [[8]] 56.0502 57.0706 60.045 69.0699 ... 312.287 330.2976 488.3218 506.3324
#> [[9]] 56.0502 57.0706 60.045 67.0542 ... 242.2115 276.2686 294.2792 312.287
#> [[10]] 56.0502 57.0706 60.045 67.0542 ... 159.0288 224.2009 276.2686 294.2792
#> ...
#> <60 more elements>
```

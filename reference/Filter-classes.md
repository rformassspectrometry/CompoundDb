# Filters supported by CompDb and IonDb

A variety of different filters can be applied to the `CompDb` object to
retrieve only subsets of the data. These filters extend the
[AnnotationFilter::AnnotationFilter](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html)
class and support the filtering concepts introduced by Bioconductor's
`AnnotationFilter` package.

The supported filters are:

- `CompoundIdFilter`: filter based on the compound ID.

- `FormulaFilter`: filter based on the compound's formula.

- `InchiFilter`: filter based on the compound's InChI.

- `InchikeyFilter`: filter based on the compound's InChI key.

- `ExactmassFilter`: filter based on the compound's (exact) mass.

- `NameFilter`: filter based on the compound name.

- `MsmsMzRangeMinFilter`: retrieve entries based on the smallest m/z of
  all peaks of their MS/MS spectra. Requires that MS/MS spectra data are
  present (i.e. `hasMsMsSpectra(cmp_db)` returns `TRUE`).

- `MsmsMzRangeMaxFilter`: retrieve entries based on the largest m/z of
  all peaks of their MS/MS spectra. Requires that MS/MS spectra data are
  present (i.e. `hasMsMsSpectra(cmp_db)` returns `TRUE`).

- `SpectrumIdFilter`: retrieve entries associated with the provided IDs
  of MS/MS spectra.

In addition to the filters listed above, the following ones are
supported by a IonDb (but not by a CompDb):

- `IonIdFilter`: filter based on the ion ID.

- `IonAdductFilter`: filter based on the adduct.

- `IonMzFilter`: filter based on the mz of the ion.

- `IonRtFilter`: filter based on the rt of the ion.

## Usage

``` r
CompoundIdFilter(value, condition = "==")

SpectrumIdFilter(value, condition = "==")

NameFilter(value, condition = "==")

MsmsMzRangeMinFilter(value, condition = ">=")

MsmsMzRangeMaxFilter(value, condition = "<=")

ExactmassFilter(value, condition = "==")

FormulaFilter(value, condition = "==")

InchiFilter(value, condition = "==")

InchikeyFilter(value, condition = "==")

IonIdFilter(value, condition = "==")

IonAdductFilter(value, condition = "==")

IonMzFilter(value, condition = "==")

IonRtFilter(value, condition = "==")
```

## Arguments

- value:

  The value for the filter. For details see
  [`AnnotationFilter::AnnotationFilter()`](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html).

- condition:

  The condition for the filter. For details see
  [`AnnotationFilter::AnnotationFilter()`](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html).

## Value

Constructor functions return an instance of the respective class.

## See also

[`AnnotationFilter::supportedFilters()`](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html)
for the method to list all supported filters for a `CompDb` (or a IonDb)
object.

## Author

Johannes Rainer

## Examples

``` r
## Create a filter for the compound id
cf <- CompoundIdFilter("comp_a")
cf
#> class: CompoundIdFilter 
#> condition: == 
#> value: comp_a 

## Create a filter using a formula expression
AnnotationFilter(~ compound_id == "comp_b")
#> class: CompoundIdFilter 
#> condition: == 
#> value: comp_b 

## Combine filters
AnnotationFilterList(CompoundIdFilter("a"), NameFilter("b"))
#> AnnotationFilterList of length 2 
#> compound_id == 'a' & name == 'b'

## Using a formula expression
AnnotationFilter(~ compound_id == "a" | name != "b")
#> AnnotationFilterList of length 2 
#> compound_id == 'a' | name != 'b'
```

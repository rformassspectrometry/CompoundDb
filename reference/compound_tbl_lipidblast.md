# Extract compound data from LipidBlast

`compound_tbl_lipidblast()` extracts basic compound annotations from a
LipidBlast file in (json format) downloaded from
http://mona.fiehnlab.ucdavis.edu/downloads . Note that no mass spectra
data is extracted from the json file.

## Usage

``` r
compound_tbl_lipidblast(
  file,
  collapse = character(),
  n = -1L,
  verbose = FALSE,
  BPPARAM = bpparam()
)
```

## Arguments

- file:

  `character(1)` with the name of the file name.

- collapse:

  optional `character(1)` to be used to collapse multiple values in the
  columns `"synonyms"`. See examples for details.

- n:

  `integer(1)` defining the number of rows from the json file that
  should be read and processed at a time. By default (`n = -1L`) the
  complete file is imported and processed. For large json files it is
  suggested to set e.g. `n = 100000` to enable chunk-wise processing and
  hence reduce the memory demand.

- verbose:

  `logical(1)` whether some progress information should be provided.
  Defaults to `verbose = FALSE`, but for parsing very large files
  (specifically with chunk-wise processing enabled with `n` \> 0) it
  might be helpful to set to `verbose = TRUE`.

- BPPARAM:

  `BiocParallelParam` object to configure parallel processing. Defaults
  to [`bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with general compound information (one row per compound):

- `compound_id`: the ID of the compound.

- `name`: the compound's name.

- `inchi`: the InChI of the compound.

- `inchikey`: the InChI key.

- `smiles`: the SMILES representation of the compound.

- `formula`: the chemical formula.

- `exactmass`: the compound's mass.

- `compound_class`: the class of the compound.

- `ionization_mode`: the ionization mode.

- `precursor_mz`: the precursor m/z value.

- `precursor_type`: the precursor type.

- `retention_time`: the retention time.

- `ccs`: the collision cross-section.

- `spectrum`: the spectrum data (i.e. the mass peaks, as a concatenated
  character string).

- `synonyms`: the compound's synonyms (aliases). This type of this
  column is by default a `list` to support multiple aliases per
  compound, unless argument `collapse` is provided, in which case
  multiple synonyms are pasted into a single element separated by the
  value of `collapse`.

## See also

Other compound table creation functions:
[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)

## Author

Johannes Rainer, Jan Stanstrup and Prateek Arora

## Examples

``` r

## Read compound information from a subset of HMDB
fl <- system.file("json/MoNa-LipidBlast_sub.json", package = "CompoundDb")
cmps <- compound_tbl_lipidblast(fl, n = 50000, verbose = TRUE)
#> Processed 1 elements
cmps
#> # A tibble: 8 × 15
#>   compound_id      name   inchi inchikey smiles formula exactmass compound_class
#>   <chr>            <chr>  <chr> <chr>    <chr>  <chr>       <dbl> <chr>         
#> 1 LipidBlast000001 CerP … InCh… FIIUKKI… "CCCC… C24H50…      479. CerP          
#> 2 LipidBlast000002 CerP … InCh… LVZNXIM… "CCCC… C26H54…      507. CerP          
#> 3 LipidBlast000003 CerP … InCh… KKMUWYA… "CCCC… C28H58…      535. CerP          
#> 4 LipidBlast000004 CerP … InCh… LYSHUFF… "CCCC… C30H62…      563. CerP          
#> 5 LipidBlast000005 CerP … InCh… OKWWFSI… "CCCC… C32H66…      591. CerP          
#> 6 LipidBlast000006 CerP … InCh… MKJOCAJ… "CCCC… C32H64…      589. CerP          
#> 7 LipidBlast000007 CerP … InCh… KZLFXFL… "CCCC… C34H70…      619. CerP          
#> 8 LipidBlast000008 CerP … InCh… GTSHVCJ… "CCCC… C36H74…      648. CerP          
#> # ℹ 7 more variables: ionization_mode <chr>, precursor_mz <dbl>,
#> #   precursor_type <chr>, retention_time <dbl>, ccs <dbl>, spectrum <chr>,
#> #   synonyms <list>
```

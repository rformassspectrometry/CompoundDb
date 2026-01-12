# Import compound and spectrum information from MoNa

`import_mona_sdf()` allows to import compound and spectrum information
from an SDF file from MoNa (Massbank of North America
http://mona.fiehnlab.ucdavis.edu/). This function is a convenience
function using the
[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
and
[`msms_spectra_mona()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_mona.md)
functions for data import but avoiding to read the SDF files twice.

## Usage

``` r
import_mona_sdf(x, nonStop = TRUE)
```

## Arguments

- x:

  `character(1)` being the SDF file name.

- nonStop:

  `logical(1)` whether file content specific errors should only reported
  as warnings and not break the full import process. The value of this
  parameter is passed to parameter `skipErrors` of the
  [`ChemmineR::read.SDFset()`](https://rdrr.io/pkg/ChemmineR/man/read.SDFset.html)
  function.

## Value

A `list` with elements `"compound"` and `"msms_spectrum"` containing
data.frames with compound and MS/MS spectrum data, respectively.

## Note

MoNa SDF files organize the data by individual spectra (i.e. each
element is one spectrum) and individual compounds can not easily and
consistently defined (i.e. not all entries have an InChI ID or other
means to uniquely identify compounds). Thus, the function returns a
highly redundant compound table. Feedback on how to reduce this
redundancy would be highly welcome!

## See also

[`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
to read only the compound information.

[`msms_spectra_mona()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_mona.md)
to read only the spectrum data.

## Author

Johannes Rainer

## Examples

``` r

## Define the test file containing a small subset from MoNa
fl <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
    package = "CompoundDb")

## Import the data
res <- import_mona_sdf(fl)
#> Reading SDF file ... 
#> OK
#> Extracting compound information ... 
#> Warning: MoNa data can currently not be normalized and the compound table contains thus highly redundant data.
#> OK
#> Extracting spectrum information ... 
#> OK
```

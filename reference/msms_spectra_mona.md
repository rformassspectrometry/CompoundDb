# Import MS/MS spectra from MoNa

`msms_spectra_mona()` imports MS/MS spectra from a MoNa (Massbank of
North America, http://mona.fiehnlab.ucdavis.edu/downloads) SDF file and
returns the data as a `data.frame`.

Depending on the parameter `collapsed`, the returned `data.frame` is
either *collapsed*, meaning that each row represents data from one
spectrum, or *expanded* with one row for each m/z and intensity pair for
each spectrum. Columns `"mz"` and `"intensity"` are of type `list` for
`collapsed = TRUE` and `numeric` for `collapsed = FALSE`.

## Usage

``` r
msms_spectra_mona(x, collapsed = TRUE)
```

## Arguments

- x:

  `character(1)`: with the path to directory containing the xml files.

- collapsed:

  `logical(1)` whether the returned `data.frame` should be *collapsed*
  or *expanded*. See description for more details.

## Value

`data.frame` with as many rows as there are peaks and columns:

- spectrum_id (`integer`): an arbitrary, unique ID for each spectrum.

- original_spectrum_id (`character`): The ID from the spectrum as
  specified in the MoNa SDF.

- compound_id (`character`): the compound ID the spectrum is associated
  with.

- polarity (`integer`): 0 for negative, 1 for positive, `NA` for not
  set.

- collision_energy (`character`): collision energy voltage.

- predicted (`logical`): whether the spectrum is predicted or
  experimentally verified.

- splash (`character`): `NA` since SPLASH (SPectraL hASH) keys are not
  provided.

- instrument_type (`character`): the type of MS instrument on which the
  spectrum was measured.

- instrument (`character`): the MS instrument.

- precursor_mz (`numeric`): precursor m/z.

- adduct (`character`): ion formed from the precursor ion.

- ms_level (`integer`): stage of the sequential mass spectrometry (MSn).

- mz (`numeric` or `list` of `numeric`): m/z values of the spectrum.

- intensity (`numeric` or `list` of `numeric`): intensity of the
  spectrum.

## Note

The identifiers provided by MoNa are used as *original_spectrum_id*.
Note also that the MoNa data is not normalized in the sense that each
spectrum is associated to one compound and the compound data is
partially redundant. Also, MoNa does not provide a *splash* for a
spectrum, hence the corresponding column will only contain `NA`.

## See also

[`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
for the function to create a
[CompDb](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
database with compound annotation and spectrum data.

Other spectrum data import functions.:
[`msms_spectra_hmdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_hmdb.md)

## Author

Johannes Rainer

## Examples

``` r

## Define the test file containing the data
fl <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
    package = "CompoundDb")
## Import spectrum data from the SDF file with a subset of the MoNa data
msms_spectra_mona(fl)
#> Reading SDF file ... 
#> OK
#> Extracting spectrum information ... 
#> OK
#>      original_spectrum_id compound_id polarity  collision_energy predicted
#> CMP1             AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> CMP2             AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> CMP3             AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> CMP4             AU100802        CMP4        1             20 eV        NA
#> CMP5             AU100803        CMP5        1             30 eV        NA
#> CMP6             AU100804        CMP6        1             40 eV        NA
#> CMP7             AU100805        CMP7        1             50 eV        NA
#>      splash instrument_type          instrument precursor_mz adduct
#> CMP1   <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> CMP2   <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> CMP3   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> CMP4   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> CMP5   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> CMP6   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> CMP7   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#>      spectrum_type spectrum_id           mz    intensity ms_level
#> CMP1           MS2           1 53.0389,.... 0.594951....        2
#> CMP2           MS2           2 53.0389,.... 0.594951....        2
#> CMP3           MS2           3 53.0379,.... 0.894101....        2
#> CMP4           MS2           4 122.0703.... 0.766124....        2
#> CMP5           MS2           5 108.0441.... 1.285794....        2
#> CMP6           MS2           6 108.0445.... 1.153673....        2
#> CMP7           MS2           7 108.0453.... 1.770916....        2

## Import the data as an *expanded* data frame, i.e. with a row for each
## single m/z (intensity) value.
msms_spectra_mona(fl, collapsed = FALSE)
#> Reading SDF file ... 
#> OK
#> Extracting spectrum information ... 
#> OK
#>     original_spectrum_id compound_id polarity  collision_energy predicted
#> 1               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 2               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 3               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 4               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 5               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 6               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 7               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 8               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 9               AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 10              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 11              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 12              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 13              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 14              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 15              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 16              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 17              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 18              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 19              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 20              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 21              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 22              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 23              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 24              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 25              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 26              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 27              AU100601        CMP1        1 Ramp 21.1-31.6 eV        NA
#> 28              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 29              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 30              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 31              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 32              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 33              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 34              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 35              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 36              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 37              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 38              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 39              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 40              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 41              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 42              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 43              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 44              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 45              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 46              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 47              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 48              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 49              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 50              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 51              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 52              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 53              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 54              AU100701        CMP2        1 Ramp 21.1-31.6 eV        NA
#> 55              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 56              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 57              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 58              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 59              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 60              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 61              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 62              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 63              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 64              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 65              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 66              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 67              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 68              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 69              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 70              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 71              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 72              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 73              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 74              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 75              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 76              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 77              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 78              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 79              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 80              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 81              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 82              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 83              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 84              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 85              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 86              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 87              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 88              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 89              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 90              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 91              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 92              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 93              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 94              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 95              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 96              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 97              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 98              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 99              AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 100             AU100801        CMP3        1 Ramp 20.8-31.3 eV        NA
#> 101             AU100802        CMP4        1             20 eV        NA
#> 102             AU100802        CMP4        1             20 eV        NA
#> 103             AU100802        CMP4        1             20 eV        NA
#> 104             AU100802        CMP4        1             20 eV        NA
#> 105             AU100802        CMP4        1             20 eV        NA
#> 106             AU100802        CMP4        1             20 eV        NA
#> 107             AU100802        CMP4        1             20 eV        NA
#> 108             AU100802        CMP4        1             20 eV        NA
#> 109             AU100802        CMP4        1             20 eV        NA
#> 110             AU100802        CMP4        1             20 eV        NA
#> 111             AU100802        CMP4        1             20 eV        NA
#> 112             AU100802        CMP4        1             20 eV        NA
#> 113             AU100802        CMP4        1             20 eV        NA
#> 114             AU100802        CMP4        1             20 eV        NA
#> 115             AU100802        CMP4        1             20 eV        NA
#> 116             AU100802        CMP4        1             20 eV        NA
#> 117             AU100803        CMP5        1             30 eV        NA
#> 118             AU100803        CMP5        1             30 eV        NA
#> 119             AU100803        CMP5        1             30 eV        NA
#> 120             AU100803        CMP5        1             30 eV        NA
#> 121             AU100803        CMP5        1             30 eV        NA
#> 122             AU100803        CMP5        1             30 eV        NA
#> 123             AU100803        CMP5        1             30 eV        NA
#> 124             AU100803        CMP5        1             30 eV        NA
#> 125             AU100803        CMP5        1             30 eV        NA
#> 126             AU100803        CMP5        1             30 eV        NA
#> 127             AU100803        CMP5        1             30 eV        NA
#> 128             AU100803        CMP5        1             30 eV        NA
#> 129             AU100803        CMP5        1             30 eV        NA
#> 130             AU100803        CMP5        1             30 eV        NA
#> 131             AU100803        CMP5        1             30 eV        NA
#> 132             AU100803        CMP5        1             30 eV        NA
#> 133             AU100803        CMP5        1             30 eV        NA
#> 134             AU100804        CMP6        1             40 eV        NA
#> 135             AU100804        CMP6        1             40 eV        NA
#> 136             AU100804        CMP6        1             40 eV        NA
#> 137             AU100804        CMP6        1             40 eV        NA
#> 138             AU100804        CMP6        1             40 eV        NA
#> 139             AU100804        CMP6        1             40 eV        NA
#> 140             AU100804        CMP6        1             40 eV        NA
#> 141             AU100804        CMP6        1             40 eV        NA
#> 142             AU100804        CMP6        1             40 eV        NA
#> 143             AU100804        CMP6        1             40 eV        NA
#> 144             AU100804        CMP6        1             40 eV        NA
#> 145             AU100804        CMP6        1             40 eV        NA
#> 146             AU100804        CMP6        1             40 eV        NA
#> 147             AU100804        CMP6        1             40 eV        NA
#> 148             AU100804        CMP6        1             40 eV        NA
#> 149             AU100804        CMP6        1             40 eV        NA
#> 150             AU100804        CMP6        1             40 eV        NA
#> 151             AU100804        CMP6        1             40 eV        NA
#> 152             AU100804        CMP6        1             40 eV        NA
#> 153             AU100804        CMP6        1             40 eV        NA
#> 154             AU100804        CMP6        1             40 eV        NA
#> 155             AU100804        CMP6        1             40 eV        NA
#> 156             AU100805        CMP7        1             50 eV        NA
#> 157             AU100805        CMP7        1             50 eV        NA
#> 158             AU100805        CMP7        1             50 eV        NA
#> 159             AU100805        CMP7        1             50 eV        NA
#> 160             AU100805        CMP7        1             50 eV        NA
#> 161             AU100805        CMP7        1             50 eV        NA
#> 162             AU100805        CMP7        1             50 eV        NA
#> 163             AU100805        CMP7        1             50 eV        NA
#> 164             AU100805        CMP7        1             50 eV        NA
#> 165             AU100805        CMP7        1             50 eV        NA
#> 166             AU100805        CMP7        1             50 eV        NA
#> 167             AU100805        CMP7        1             50 eV        NA
#> 168             AU100805        CMP7        1             50 eV        NA
#> 169             AU100805        CMP7        1             50 eV        NA
#> 170             AU100805        CMP7        1             50 eV        NA
#> 171             AU100805        CMP7        1             50 eV        NA
#> 172             AU100805        CMP7        1             50 eV        NA
#> 173             AU100805        CMP7        1             50 eV        NA
#> 174             AU100805        CMP7        1             50 eV        NA
#> 175             AU100805        CMP7        1             50 eV        NA
#> 176             AU100805        CMP7        1             50 eV        NA
#> 177             AU100805        CMP7        1             50 eV        NA
#> 178             AU100805        CMP7        1             50 eV        NA
#> 179             AU100805        CMP7        1             50 eV        NA
#>     splash instrument_type          instrument precursor_mz adduct
#> 1     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 2     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 3     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 4     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 5     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 6     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 7     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 8     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 9     <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 10    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 11    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 12    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 13    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 14    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 15    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 16    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 17    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 18    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 19    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 20    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 21    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 22    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 23    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 24    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 25    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 26    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 27    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 28    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 29    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 30    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 31    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 32    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 33    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 34    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 35    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 36    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 37    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 38    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 39    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 40    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 41    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 42    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 43    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 44    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 45    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 46    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 47    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 48    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 49    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 50    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 51    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 52    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 53    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 54    <NA>     LC-ESI-QTOF Bruker maXis Impact     285.0208 [M+H]+
#> 55    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 56    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 57    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 58    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 59    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 60    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 61    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 62    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 63    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 64    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 65    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 66    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 67    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 68    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 69    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 70    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 71    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 72    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 73    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 74    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 75    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 76    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 77    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 78    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 79    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 80    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 81    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 82    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 83    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 84    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 85    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 86    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 87    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 88    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 89    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 90    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 91    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 92    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 93    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 94    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 95    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 96    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 97    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 98    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 99    <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 100   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 101   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 102   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 103   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 104   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 105   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 106   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 107   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 108   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 109   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 110   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 111   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 112   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 113   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 114   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 115   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 116   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 117   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 118   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 119   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 120   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 121   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 122   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 123   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 124   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 125   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 126   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 127   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 128   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 129   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 130   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 131   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 132   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 133   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 134   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 135   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 136   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 137   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 138   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 139   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 140   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 141   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 142   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 143   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 144   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 145   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 146   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 147   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 148   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 149   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 150   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 151   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 152   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 153   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 154   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 155   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 156   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 157   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 158   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 159   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 160   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 161   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 162   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 163   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 164   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 165   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 166   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 167   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 168   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 169   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 170   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 171   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 172   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 173   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 174   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 175   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 176   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 177   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 178   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#> 179   <NA>     LC-ESI-QTOF Bruker maXis Impact     279.0910 [M+H]+
#>     spectrum_type spectrum_id       mz  intensity ms_level
#> 1             MS2           1  53.0389   0.594951        2
#> 2             MS2           1  54.0333   0.566811        2
#> 3             MS2           1  55.0178   0.522592        2
#> 4             MS2           1  60.0552   0.542692        2
#> 5             MS2           1  65.0382   3.822962        2
#> 6             MS2           1  66.0423   0.506512        2
#> 7             MS2           1  68.0490   7.963499        2
#> 8             MS2           1  78.0333   0.727609        2
#> 9             MS2           1  79.0177   1.057244        2
#> 10            MS2           1  92.0498   7.702203        2
#> 11            MS2           1  93.0532   0.731629        2
#> 12            MS2           1  96.0443   0.623091        2
#> 13            MS2           1 108.0457  12.172375        2
#> 14            MS2           1 109.0483   1.181862        2
#> 15            MS2           1 110.0609   4.904325        2
#> 16            MS2           1 120.0562   3.095353        2
#> 17            MS2           1 130.0172   5.656054        2
#> 18            MS2           1 132.0138   1.515517        2
#> 19            MS2           1 156.0118 100.000000        2
#> 20            MS2           1 157.0150   8.884065        2
#> 21            MS2           1 158.0080   3.891301        2
#> 22            MS2           1 174.0228   0.751729        2
#> 23            MS2           1 184.0757   0.619071        2
#> 24            MS2           1 191.9647   0.590931        2
#> 25            MS2           1 219.0438   0.723589        2
#> 26            MS2           1 285.0221   3.694324        2
#> 27            MS2           1 287.0184   0.840167        2
#> 28            MS2           2  53.0389   0.594951        2
#> 29            MS2           2  54.0333   0.566811        2
#> 30            MS2           2  55.0178   0.522592        2
#> 31            MS2           2  60.0552   0.542692        2
#> 32            MS2           2  65.0382   3.822962        2
#> 33            MS2           2  66.0423   0.506512        2
#> 34            MS2           2  68.0490   7.963499        2
#> 35            MS2           2  78.0333   0.727609        2
#> 36            MS2           2  79.0177   1.057244        2
#> 37            MS2           2  92.0498   7.702203        2
#> 38            MS2           2  93.0532   0.731629        2
#> 39            MS2           2  96.0443   0.623091        2
#> 40            MS2           2 108.0457  12.172375        2
#> 41            MS2           2 109.0483   1.181862        2
#> 42            MS2           2 110.0609   4.904325        2
#> 43            MS2           2 120.0562   3.095353        2
#> 44            MS2           2 130.0172   5.656054        2
#> 45            MS2           2 132.0138   1.515517        2
#> 46            MS2           2 156.0118 100.000000        2
#> 47            MS2           2 157.0150   8.884065        2
#> 48            MS2           2 158.0080   3.891301        2
#> 49            MS2           2 174.0228   0.751729        2
#> 50            MS2           2 184.0757   0.619071        2
#> 51            MS2           2 191.9647   0.590931        2
#> 52            MS2           2 219.0438   0.723589        2
#> 53            MS2           2 285.0221   3.694324        2
#> 54            MS2           2 287.0184   0.840167        2
#> 55            MS2           3  53.0379   0.894101        2
#> 56            MS2           3  54.0335   0.661867        2
#> 57            MS2           3  55.0176   0.598003        2
#> 58            MS2           3  65.0379   8.717487        2
#> 59            MS2           3  68.0491  13.013818        2
#> 60            MS2           3  69.0329   1.640153        2
#> 61            MS2           3  78.0334   1.477589        2
#> 62            MS2           3  79.0178   2.261379        2
#> 63            MS2           3  80.0489   1.431143        2
#> 64            MS2           3  81.0444   1.950766        2
#> 65            MS2           3  82.0284   0.606712        2
#> 66            MS2           3  92.0499  30.585230        2
#> 67            MS2           3  93.0558   2.844868        2
#> 68            MS2           3  94.0647   1.686600        2
#> 69            MS2           3  95.0608   3.027752        2
#> 70            MS2           3  96.0443   1.300511        2
#> 71            MS2           3 108.0461  33.946818        2
#> 72            MS2           3 109.0497   2.360079        2
#> 73            MS2           3 110.0616   6.107757        2
#> 74            MS2           3 111.0651   0.519624        2
#> 75            MS2           3 120.0565   1.962378        2
#> 76            MS2           3 122.0716   6.078727        2
#> 77            MS2           3 123.0794   2.246865        2
#> 78            MS2           3 124.0872  71.211681        2
#> 79            MS2           3 125.0905   6.398049        2
#> 80            MS2           3 126.0663  17.911054        2
#> 81            MS2           3 127.0697   0.595100        2
#> 82            MS2           3 156.0117  82.855318        2
#> 83            MS2           3 157.0148   5.739085        2
#> 84            MS2           3 158.0072   1.544357        2
#> 85            MS2           3 174.0224   1.106015        2
#> 86            MS2           3 186.0334  11.263353        2
#> 87            MS2           3 187.0368   0.775081        2
#> 88            MS2           3 188.0128   1.637250        2
#> 89            MS2           3 188.0291   0.534138        2
#> 90            MS2           3 204.0445 100.000000        2
#> 91            MS2           3 205.0473   6.972829        2
#> 92            MS2           3 206.0406   3.358686        2
#> 93            MS2           3 213.1141  18.259405        2
#> 94            MS2           3 214.1167   2.241059        2
#> 95            MS2           3 215.0927   3.071296        2
#> 96            MS2           3 215.1291   1.320831        2
#> 97            MS2           3 279.0925  61.483976        2
#> 98            MS2           3 280.0953   8.438806        2
#> 99            MS2           3 281.0725   7.837901        2
#> 100           MS2           3 282.0742   1.222132        2
#> 101           MS2           4 122.0703   0.766124        2
#> 102           MS2           4 124.0861  36.693459        2
#> 103           MS2           4 125.0892   1.930893        2
#> 104           MS2           4 149.0227   0.828453        2
#> 105           MS2           4 156.0104  53.249536        2
#> 106           MS2           4 157.0129   2.999571        2
#> 107           MS2           4 158.0061   1.778967        2
#> 108           MS2           4 174.0209   0.627183        2
#> 109           MS2           4 186.0321  22.621444        2
#> 110           MS2           4 187.0346   1.719235        2
#> 111           MS2           4 188.0285   0.646661        2
#> 112           MS2           4 204.0431 100.000000        2
#> 113           MS2           4 213.1128   8.749399        2
#> 114           MS2           4 214.1159   1.407591        2
#> 115           MS2           4 215.1281   0.658348        2
#> 116           MS2           4 279.0909  80.894937        2
#> 117           MS2           5 108.0441   1.285794        2
#> 118           MS2           5 122.0704   6.630847        2
#> 119           MS2           5 123.0781   2.170942        2
#> 120           MS2           5 124.0861 100.000000        2
#> 121           MS2           5 125.0889   6.093546        2
#> 122           MS2           5 149.0221   1.388285        2
#> 123           MS2           5 156.0106  50.043481        2
#> 124           MS2           5 158.0064   1.615007        2
#> 125           MS2           5 186.0323  15.118951        2
#> 126           MS2           5 187.0355   1.323064        2
#> 127           MS2           5 196.0858   1.220573        2
#> 128           MS2           5 204.0429  70.964035        2
#> 129           MS2           5 205.0455   4.931983        2
#> 130           MS2           5 213.1123  22.610100        2
#> 131           MS2           5 214.1155   3.003292        2
#> 132           MS2           5 215.1283   0.804398        2
#> 133           MS2           5 279.0903   3.580968        2
#> 134           MS2           6 108.0445   1.153673        2
#> 135           MS2           6 122.0702   5.323878        2
#> 136           MS2           6 123.0772   2.202467        2
#> 137           MS2           6 124.0862 100.000000        2
#> 138           MS2           6 125.0890   6.847126        2
#> 139           MS2           6 134.0701   0.714179        2
#> 140           MS2           6 149.0224   1.747990        2
#> 141           MS2           6 154.0624   0.644259        2
#> 142           MS2           6 155.0685   0.624282        2
#> 143           MS2           6 156.0104  10.373071        2
#> 144           MS2           6 157.0126   0.933926        2
#> 145           MS2           6 172.0852   0.564351        2
#> 146           MS2           6 186.0324   3.845578        2
#> 147           MS2           6 196.0852   5.209010        2
#> 148           MS2           6 197.0903   1.378415        2
#> 149           MS2           6 198.0888   2.362283        2
#> 150           MS2           6 204.0427  15.422264        2
#> 151           MS2           6 205.0463   0.869001        2
#> 152           MS2           6 206.0375   0.759127        2
#> 153           MS2           6 212.1036   0.659242        2
#> 154           MS2           6 213.1121  18.109174        2
#> 155           MS2           6 214.1152   2.577036        2
#> 156           MS2           7 108.0453   1.770916        2
#> 157           MS2           7 122.0703   2.803951        2
#> 158           MS2           7 123.0780   2.792598        2
#> 159           MS2           7 124.0859 100.000000        2
#> 160           MS2           7 125.0891   7.901010        2
#> 161           MS2           7 149.0231   1.623340        2
#> 162           MS2           7 154.0639   2.111477        2
#> 163           MS2           7 155.0605   2.463390        2
#> 164           MS2           7 155.0714   2.690430        2
#> 165           MS2           7 156.0100   2.713134        2
#> 166           MS2           7 169.0745   1.475763        2
#> 167           MS2           7 171.0781   1.555228        2
#> 168           MS2           7 172.0869   1.271427        2
#> 169           MS2           7 181.0634   0.930866        2
#> 170           MS2           7 186.1022   1.033034        2
#> 171           MS2           7 195.0786   1.555228        2
#> 172           MS2           7 196.0859   7.628562        2
#> 173           MS2           7 197.0856   3.871041        2
#> 174           MS2           7 198.0886   5.903054        2
#> 175           MS2           7 199.0904   0.998978        2
#> 176           MS2           7 204.0438   2.622318        2
#> 177           MS2           7 212.1048   2.327165        2
#> 178           MS2           7 213.1122   9.342718        2
#> 179           MS2           7 214.1153   1.725508        2
```

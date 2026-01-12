# Import MS/MS spectra from HMDB xml files

`msms_spectra_hmdb()` imports MS/MS spectra from corresponding xml files
from HMDB (http://www.hmdb.ca) and returns the data as a `data.frame`.
HMDB stores MS/MS spectrum data in xml files, one file per spectrum.

Depending on the parameter `collapsed`, the returned `data.frame` is
either *collapsed*, meaning that each row represents data from one
spectrum xml file, or *expanded* with one row for each m/z and intensity
pair for each spectrum. Columns `"mz"` and `"intensity"` are of type
`list` for `collapsed = TRUE` and `numeric` for `collapsed = FALSE`.

## Usage

``` r
msms_spectra_hmdb(x, collapsed = TRUE)
```

## Arguments

- x:

  `character(1)`: with the path to directory containing the xml files.

- collapsed:

  `logical(1)` whether the returned `data.frame` should be *collapsed*
  or *expanded*. See description for more details.

## Value

`data.frame` with as many rows as there are peaks and columns:

- spectrum_id (`integer`): an arbitrary, unique ID identifying values
  from one xml file.

- original_spectrum_id (`character`): the HMDB-internal ID of the
  spectrum.

- compound_id (`character`): the HMDB compound ID the spectrum is
  associated with.

- polarity (`integer`): 0 for negative, 1 for positive, `NA` for not
  set.

- collision_energy (`numeric`): collision energy voltage.

- predicted (`logical`): whether the spectrum is predicted or
  experimentally verified.

- splash (`character`): the SPLASH (SPectraL hASH) key of the spectrum
  (Wohlgemuth 2016).

- instrument_type (`character`): the type of MS instrument on which the
  spectrum was measured.

- instrument (`character`): the MS instrument (not available for all
  spectra in HMDB).

- precursor_mz (`numeric`): not provided by HMDB and thus `NA`.

- mz (`numeric` or `list` of `numeric`): m/z values of the spectrum.

- intensity (`numeric` or `list` of `numeric`): intensity of the
  spectrum.

## Note

The HMDB xml files are supposed to be extracted from the downloaded zip
file into a folder and should not be renamed. The function identifies
xml files containing MS/MS spectra by their file name.

The same spectrum ID can be associated with multiple compounds. Thus,
the function assignes an arbitrary ID (column `"spectrum_id"`) to values
from each file. The original ID of the spectrum in HMDB is provided in
column `"original_spectrum_id"`.

## References

Wohlgemuth G, Mehta SS, Mejia RF, Neumann S, Pedrosa D, Pluskal T,
Schymanski EL, Willighagen EL, Wilson M, Wishart DS, Arita M, Dorrestein
PC, Bandeira N, Wang M, Schulze T, Selak RM, Steinbeck C, Nainala VC,
Mistrik R, Nishioka T, Fiehn O. SPLASH, A hashed identifier for mass
spectra. Nature Biotechnology 2016 34(11):1099-1101

## See also

[`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
for the function to create a
[CompDb](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
database with compound annotation and spectrum data.

Other spectrum data import functions.:
[`msms_spectra_mona()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_mona.md)

## Author

Johannes Rainer

## Examples

``` r

## Locate the folder within the package containing test xml files.
pth <- system.file("xml", package = "CompoundDb")

## List all files in that directory
dir(pth)
#> [1] "HMDB0000001_ms_ms_spectrum_1_experimental.xml"     
#> [2] "HMDB0000001_ms_ms_spectrum_2_experimental.xml"     
#> [3] "HMDB0004370_ms_ms_spectrum_446575_experimental.xml"
#> [4] "HMDB0006719_ms_ms_spectrum_370739_predicted.xml"   
#> [5] "fail"                                              

## Import spectrum data from HMDB MS/MS spectrum xml files in that directory
msms_spectra_hmdb(pth)
#> Going to process 4 xml files.
#> Postprocessing data ... 
#> OK
#>   original_spectrum_id compound_id polarity collision_energy predicted
#> 1                    1 HMDB0000001        1               10     FALSE
#> 2                    2 HMDB0000001        1               25     FALSE
#> 3               446575 HMDB0004370        1               NA     FALSE
#> 4               370739 HMDB0006719        0               20      TRUE
#>                                          splash instrument_type
#> 1 splash10-00di-0900000000-037d24a7d65676b7e356     Quattro_QQQ
#> 2 splash10-00di-0900000000-03e99316bd6c098f5d11     Quattro_QQQ
#> 3 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 4 splash10-001i-0009000000-c16445787a057993951a            <NA>
#>                    instrument precursor_mz           mz    intensity
#> 1                        <NA>           NA 109.2, 1.... 3.406999....
#> 2                        <NA>           NA 83.1, 96.... 6.685028....
#> 3 API3000, Applied Biosystems           NA 44.1, 57.... 0.051124....
#> 4                        <NA>           NA 111.0815.... 0.550973....
#>   spectrum_id
#> 1           1
#> 2           2
#> 3           3
#> 4           4

## Import the data as an *expanded* data frame, i.e. with a row for each
## single m/z (intensity) value.
msms_spectra_hmdb(pth, collapsed = FALSE)
#> Going to process 4 xml files.
#> Postprocessing data ... 
#> OK
#>    original_spectrum_id compound_id polarity collision_energy predicted
#> 1                     1 HMDB0000001        1               10     FALSE
#> 2                     1 HMDB0000001        1               10     FALSE
#> 3                     1 HMDB0000001        1               10     FALSE
#> 4                     1 HMDB0000001        1               10     FALSE
#> 5                     1 HMDB0000001        1               10     FALSE
#> 6                     2 HMDB0000001        1               25     FALSE
#> 7                     2 HMDB0000001        1               25     FALSE
#> 8                     2 HMDB0000001        1               25     FALSE
#> 9                     2 HMDB0000001        1               25     FALSE
#> 10                    2 HMDB0000001        1               25     FALSE
#> 11                    2 HMDB0000001        1               25     FALSE
#> 12                    2 HMDB0000001        1               25     FALSE
#> 13               446575 HMDB0004370        1               NA     FALSE
#> 14               446575 HMDB0004370        1               NA     FALSE
#> 15               446575 HMDB0004370        1               NA     FALSE
#> 16               446575 HMDB0004370        1               NA     FALSE
#> 17               446575 HMDB0004370        1               NA     FALSE
#> 18               446575 HMDB0004370        1               NA     FALSE
#> 19               446575 HMDB0004370        1               NA     FALSE
#> 20               446575 HMDB0004370        1               NA     FALSE
#> 21               446575 HMDB0004370        1               NA     FALSE
#> 22               446575 HMDB0004370        1               NA     FALSE
#> 23               446575 HMDB0004370        1               NA     FALSE
#> 24               446575 HMDB0004370        1               NA     FALSE
#> 25               446575 HMDB0004370        1               NA     FALSE
#> 26               446575 HMDB0004370        1               NA     FALSE
#> 27               446575 HMDB0004370        1               NA     FALSE
#> 28               446575 HMDB0004370        1               NA     FALSE
#> 29               446575 HMDB0004370        1               NA     FALSE
#> 30               446575 HMDB0004370        1               NA     FALSE
#> 31               446575 HMDB0004370        1               NA     FALSE
#> 32               446575 HMDB0004370        1               NA     FALSE
#> 33               446575 HMDB0004370        1               NA     FALSE
#> 34               446575 HMDB0004370        1               NA     FALSE
#> 35               446575 HMDB0004370        1               NA     FALSE
#> 36               446575 HMDB0004370        1               NA     FALSE
#> 37               446575 HMDB0004370        1               NA     FALSE
#> 38               446575 HMDB0004370        1               NA     FALSE
#> 39               446575 HMDB0004370        1               NA     FALSE
#> 40               446575 HMDB0004370        1               NA     FALSE
#> 41               446575 HMDB0004370        1               NA     FALSE
#> 42               446575 HMDB0004370        1               NA     FALSE
#> 43               446575 HMDB0004370        1               NA     FALSE
#> 44               446575 HMDB0004370        1               NA     FALSE
#> 45               446575 HMDB0004370        1               NA     FALSE
#> 46               370739 HMDB0006719        0               20      TRUE
#> 47               370739 HMDB0006719        0               20      TRUE
#> 48               370739 HMDB0006719        0               20      TRUE
#> 49               370739 HMDB0006719        0               20      TRUE
#> 50               370739 HMDB0006719        0               20      TRUE
#> 51               370739 HMDB0006719        0               20      TRUE
#> 52               370739 HMDB0006719        0               20      TRUE
#> 53               370739 HMDB0006719        0               20      TRUE
#> 54               370739 HMDB0006719        0               20      TRUE
#> 55               370739 HMDB0006719        0               20      TRUE
#>                                           splash instrument_type
#> 1  splash10-00di-0900000000-037d24a7d65676b7e356     Quattro_QQQ
#> 2  splash10-00di-0900000000-037d24a7d65676b7e356     Quattro_QQQ
#> 3  splash10-00di-0900000000-037d24a7d65676b7e356     Quattro_QQQ
#> 4  splash10-00di-0900000000-037d24a7d65676b7e356     Quattro_QQQ
#> 5  splash10-00di-0900000000-037d24a7d65676b7e356     Quattro_QQQ
#> 6  splash10-00di-0900000000-03e99316bd6c098f5d11     Quattro_QQQ
#> 7  splash10-00di-0900000000-03e99316bd6c098f5d11     Quattro_QQQ
#> 8  splash10-00di-0900000000-03e99316bd6c098f5d11     Quattro_QQQ
#> 9  splash10-00di-0900000000-03e99316bd6c098f5d11     Quattro_QQQ
#> 10 splash10-00di-0900000000-03e99316bd6c098f5d11     Quattro_QQQ
#> 11 splash10-00di-0900000000-03e99316bd6c098f5d11     Quattro_QQQ
#> 12 splash10-00di-0900000000-03e99316bd6c098f5d11     Quattro_QQQ
#> 13 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 14 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 15 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 16 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 17 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 18 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 19 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 20 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 21 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 22 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 23 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 24 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 25 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 26 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 27 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 28 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 29 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 30 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 31 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 32 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 33 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 34 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 35 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 36 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 37 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 38 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 39 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 40 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 41 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 42 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 43 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 44 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 45 splash10-004i-0900000000-9f3d3580f44ff890d948       LC-ESI-QQ
#> 46 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 47 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 48 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 49 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 50 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 51 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 52 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 53 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 54 splash10-001i-0009000000-c16445787a057993951a            <NA>
#> 55 splash10-001i-0009000000-c16445787a057993951a            <NA>
#>                     instrument precursor_mz       mz   intensity spectrum_id
#> 1                         <NA>           NA 109.2000   3.4069998           1
#> 2                         <NA>           NA 124.2000  47.4945731           1
#> 3                         <NA>           NA 124.5000   3.0943658           1
#> 4                         <NA>           NA 170.1600 100.0000000           1
#> 5                         <NA>           NA 170.5200  13.2396932           1
#> 6                         <NA>           NA  83.1000   6.6850283           2
#> 7                         <NA>           NA  96.1200   4.3812987           2
#> 8                         <NA>           NA  97.1400   3.0221439           2
#> 9                         <NA>           NA 109.1400  16.7082568           2
#> 10                        <NA>           NA 124.0800 100.0000000           2
#> 11                        <NA>           NA 125.1000   4.5651409           2
#> 12                        <NA>           NA 170.1600  40.6434125           2
#> 13 API3000, Applied Biosystems           NA  44.1000   0.0511240           3
#> 14 API3000, Applied Biosystems           NA  57.9000   0.0065970           3
#> 15 API3000, Applied Biosystems           NA  61.4000   0.0120940           3
#> 16 API3000, Applied Biosystems           NA  71.2000   0.0016490           3
#> 17 API3000, Applied Biosystems           NA  73.8000   0.0032980           3
#> 18 API3000, Applied Biosystems           NA  78.3000   0.0060470           3
#> 19 API3000, Applied Biosystems           NA  78.8000   0.0049480           3
#> 20 API3000, Applied Biosystems           NA  83.1000   0.0153920           3
#> 21 API3000, Applied Biosystems           NA  88.2000   0.0038480           3
#> 22 API3000, Applied Biosystems           NA  98.8000   0.0087960           3
#> 23 API3000, Applied Biosystems           NA 101.7000   0.0038480           3
#> 24 API3000, Applied Biosystems           NA 102.8000   0.0021990           3
#> 25 API3000, Applied Biosystems           NA 106.8000   0.0109940           3
#> 26 API3000, Applied Biosystems           NA 111.6000   0.0038480           3
#> 27 API3000, Applied Biosystems           NA 115.1000   0.1016990           3
#> 28 API3000, Applied Biosystems           NA 116.9000   0.0555220           3
#> 29 API3000, Applied Biosystems           NA 118.1000   0.4254850           3
#> 30 API3000, Applied Biosystems           NA 121.0000   0.0082460           3
#> 31 API3000, Applied Biosystems           NA 124.9000   0.0087960           3
#> 32 API3000, Applied Biosystems           NA 126.7000   0.0032980           3
#> 33 API3000, Applied Biosystems           NA 127.2000   0.0115440           3
#> 34 API3000, Applied Biosystems           NA 128.9000   0.0038480           3
#> 35 API3000, Applied Biosystems           NA 129.9000   0.0164920           3
#> 36 API3000, Applied Biosystems           NA 131.9000   8.9527790           3
#> 37 API3000, Applied Biosystems           NA 139.0000   0.0241880           3
#> 38 API3000, Applied Biosystems           NA 140.4000   0.0115440           3
#> 39 API3000, Applied Biosystems           NA 141.2000   0.0588200           3
#> 40 API3000, Applied Biosystems           NA 142.9000   0.0571710           3
#> 41 API3000, Applied Biosystems           NA 144.1000  36.3817270           3
#> 42 API3000, Applied Biosystems           NA 157.6000   0.0357320           3
#> 43 API3000, Applied Biosystems           NA 158.0000   0.3848060           3
#> 44 API3000, Applied Biosystems           NA 175.2000 100.0000000           3
#> 45 API3000, Applied Biosystems           NA 193.2000   0.0087960           3
#> 46                        <NA>           NA 111.0815   0.5509733           4
#> 47                        <NA>           NA 249.2588   0.7330056           4
#> 48                        <NA>           NA 273.2588   0.5171464           4
#> 49                        <NA>           NA 341.2850   0.6527944           4
#> 50                        <NA>           NA 353.3214   0.8764642           4
#> 51                        <NA>           NA 355.3006   0.7582467           4
#> 52                        <NA>           NA 355.3370   0.4985064           4
#> 53                        <NA>           NA 365.3214  21.2824831           4
#> 54                        <NA>           NA 367.3006   2.4354476           4
#> 55                        <NA>           NA 383.3319  52.1357554           4
```

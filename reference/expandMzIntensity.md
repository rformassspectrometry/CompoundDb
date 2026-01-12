# Expand m/z and intensity values in a data.frame

`expandMzIntensity()` *expands* a `data.frame` with m/z and/or intensity
values stored as a `list` in columns `"mz"` and `"intensity"`. The
resulting `data.frame` has the m/z and intensity values stored as
`numeric` in columns `"mz"` and `"intensity"`, one value per row, with
the content of other columns repeated as many times as there are m/z and
intensity values.

## Usage

``` r
expandMzIntensity(x)
```

## Arguments

- x:

  `data.frame` with *collapsed* m/z and intensity values in columns
  `"mz"` and `"intensity"`, such as returned by
  [`msms_spectra_hmdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_hmdb.md)
  with parameter `collapsed = TRUE`, or by `spectra` or `compounds`
  calls.

## Value

`data.frame` with `"mz"` and `"intensity"` columns *expanded*. See
description for details.

## Author

Johannes Rainer

## Examples

``` r

## Read a data.frame with collapsed columns mz and intensity columns
dr <- system.file("xml/", package = "CompoundDb")
msms_spctra <- msms_spectra_hmdb(dr)
#> Going to process 4 xml files.
#> Postprocessing data ... 
#> OK

msms_spctra
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

## Columns mz and intensity are "collased"
msms_spctra$mz
#> [[1]]
#> [1] 109.20 124.20 124.50 170.16 170.52
#> 
#> [[2]]
#> [1]  83.10  96.12  97.14 109.14 124.08 125.10 170.16
#> 
#> [[3]]
#>  [1]  44.1  57.9  61.4  71.2  73.8  78.3  78.8  83.1  88.2  98.8 101.7 102.8
#> [13] 106.8 111.6 115.1 116.9 118.1 121.0 124.9 126.7 127.2 128.9 129.9 131.9
#> [25] 139.0 140.4 141.2 142.9 144.1 157.6 158.0 175.2 193.2
#> 
#> [[4]]
#>  [1] 111.0815 249.2588 273.2588 341.2850 353.3214 355.3006 355.3370 365.3214
#>  [9] 367.3006 383.3319
#> 

## Expand the data.frame to get one row per m/z and intensity value
spctra_exp <- expandMzIntensity(msms_spctra)
spctra_exp
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

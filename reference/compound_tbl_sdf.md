# Extract compound data from a file in SDF format

`compound_tbl_sdf()` extracts basic compound annotations from a file in
SDF format (structure-data file). The function currently supports SDF
files from:

- HMDB (Human Metabolome Database): http://www.hmdb.ca

- ChEBI (Chemical Entities of Biological Interest):
  http://ebi.ac.uk/chebi

- LMSD (LIPID MAPS Structure Database): http://www.lipidmaps.org

- PubChem: https://pubchem.ncbi.nlm.nih.gov/

- MoNa: http://mona.fiehnlab.ucdavis.edu/ (see notes below!)

## Usage

``` r
compound_tbl_sdf(file, collapse, onlyValid = TRUE, nonStop = TRUE)
```

## Arguments

- file:

  `character(1)` with the name of the SDF file.

- collapse:

  optional `character(1)` to be used to collapse multiple values in the
  columns `"synonyms"`. See examples for details.

- onlyValid:

  `logical(1)` to import only valid or all elements (defaults to
  `onlyValid = TRUE`)

- nonStop:

  `logical(1)` whether file content specific errors should only reported
  as warnings and not break the full import process. The value of this
  parameter is passed to parameter `skipErrors` of the
  [`ChemmineR::read.SDFset()`](https://rdrr.io/pkg/ChemmineR/man/read.SDFset.html)
  function.

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with general compound information (one row per compound):

- `compound_id`: the ID of the compound.

- `name`: the compound's name.

- `inchi`: the InChI of the compound.

- `inchikey`: the InChI key.

- `formula`: the chemical formula.

- `exactmass`: the compound's (monoisotopic exact) mass.

- `synonyms`: the compound's synonyms (aliases). This type of this
  column is by default a `list` to support multiple aliases per
  compound, unless argument `collapse` is provided, in which case
  multiple synonyms are pasted into a single element separated by the
  value of `collapse`.

- `smiles`: the compound's SMILES (if provided).

## Details

Column `"name"` reports for HMDB files the `"GENERIC_NAME"`, for ChEBI
the `"ChEBI Name"`, for PubChem the `"PUBCHEM_IUPAC_TRADITIONAL_NAME"`,
and for Lipid Maps the `"COMMON_NAME"`, if that is not available, the
first of the compounds synonyms and, if that is also not provided, the
`"SYSTEMATIC_NAME"`.

## Note

`compound_tbl_sdf()` supports also to read/process gzipped files.

MoNa SDF files organize the data by individual spectra (i.e. each
element is one spectrum) and individual compounds can not easily and
consistently defined (i.e. not all entries have an InChI ID or other
means to uniquely identify compounds). Thus, the function returns a
highly redundant compound table. Feedback on how to reduce this
redundancy would be highly welcome!

LIPID MAPS was tested August 2020. Older SDF files might not work as the
field names were changed.

## See also

[`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
for a function to create a SQLite-based compound database.

Other compound table creation functions:
[`compound_tbl_lipidblast()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_lipidblast.md)

## Author

Johannes Rainer and Jan Stanstrup

## Examples

``` r

## Read compound information from a subset of HMDB
fl <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
cmps <- compound_tbl_sdf(fl)
cmps
#> # A tibble: 9 × 8
#>   compound_id name              inchi inchikey formula exactmass synonyms smiles
#>   <chr>       <chr>             <chr> <chr>    <chr>       <dbl> <named > <chr> 
#> 1 HMDB0000001 1-Methylhistidine InCh… BRMWTNU… C7H11N…     169.  <chr>    "CN1C…
#> 2 HMDB0000002 1,3-Diaminopropa… InCh… XFNJVJP… C3H10N2      74.1 <chr>    "NCCC…
#> 3 HMDB0000005 2-Ketobutyric ac… InCh… TYEYBOS… C4H6O3      102.  <chr>    "CCC(…
#> 4 HMDB0000008 2-Hydroxybutyric… InCh… AFENDNX… C4H8O3      104.  <chr>    "CCC(…
#> 5 HMDB0000010 2-Methoxyestrone  InCh… WHEUWNK… C19H24…     300.  <chr>    "[H][…
#> 6 HMDB0000011 (R)-3-Hydroxybut… InCh… WHBMMWS… C4H8O3      104.  <chr>    "C[C@…
#> 7 HMDB0000012 Deoxyuridine      InCh… MXHRCPN… C9H12N…     228.  <chr>    "OC[C…
#> 8 HMDB0004370 N-Methyltryptami… InCh… NCIKQJB… C11H14…     174.  <chr>    "CNCC…
#> 9 HMDB0006719 5,6-trans-Vitami… InCh… QYSXJUF… C27H44O     384.  <chr>    "CC(C…

## Column synonyms contains a list
cmps$synonyms
#> $CMP1
#> [1] "1 Methylhistidine"      "1-Methyl histidine"     "1-Methyl-Histidine"    
#> [4] "1-Methyl-L-histidine"   "1-MHis"                 "1-N-Methyl-L-histidine"
#> [7] "L-1-Methylhistidine"    "N1-Methyl-L-histidine"  "Pi-methylhistidine"    
#> 
#> $CMP2
#> [1] "1,3-Diamino-N-propane"   "1,3-Propanediamine"     
#> [3] "1,3-Propylenediamine"    "1,3-Trimethylenediamine"
#> [5] "3-Aminopropylamine"      "a,w-Propanediamine"     
#> [7] "Propane-1,3-diamine"     "Trimethylenediamine"    
#> 
#> $CMP3
#>  [1] "2-Ketobutanoate"           "2-Ketobutanoic acid"      
#>  [3] "2-Ketobutyrate"            "2-Oxo-Butanoate"          
#>  [5] "2-Oxo-Butanoic acid"       "2-Oxo-Butyrate"           
#>  [7] "2-Oxo-Butyric acid"        "2-Oxo-N-butyrate"         
#>  [9] "2-Oxo-N-butyric acid"      "2-Oxobutanoate"           
#> [11] "2-Oxobutanoic acid"        "2-Oxobutyrate"            
#> [13] "2-Oxobutyric acid"         "3-Methylpyruvate"         
#> [15] "3-Methylpyruvic acid"      "a-Keto-N-Butyrate"        
#> [17] "a-Keto-N-Butyric acid"     "a-Ketobutyrate"           
#> [19] "a-Ketobutyric acid"        "a-Oxo-N-butyrate"         
#> [21] "a-Oxo-N-butyric acid"      "a-Oxobutyrate"            
#> [23] "a-Oxobutyric acid"         "alpha-Keto-N-butyrate"    
#> [25] "alpha-Keto-N-butyric acid" "alpha-Ketobutric acid"    
#> [27] "alpha-Ketobutyrate"        "alpha-Ketobutyric acid"   
#> [29] "alpha-Oxo-N-butyrate"      "alpha-Oxo-N-butyric acid" 
#> [31] "alpha-Oxobutyrate"         "alpha-Oxobutyric acid"    
#> [33] "Methyl-Pyruvate"           "Methyl-Pyruvic acid"      
#> [35] "Propionyl-formate"         "Propionyl-formic acid"    
#> 
#> $CMP4
#>  [1] "(RS)-2-Hydroxybutyrate"       "(RS)-2-Hydroxybutyric acid"  
#>  [3] "2-Hydroxy-Butanoate"          "2-Hydroxy-Butanoic acid"     
#>  [5] "2-Hydroxy-DL-Butyrate"        "2-Hydroxy-DL-Butyric acid"   
#>  [7] "2-Hydroxy-N-butyrate"         "2-Hydroxy-N-butyric acid"    
#>  [9] "2-Hydroxybutanoate"           "2-Hydroxybutanoic acid"      
#> [11] "2-Hydroxybutyrate"            "a-Hydroxy-N-butyrate"        
#> [13] "a-Hydroxy-N-butyric acid"     "a-Hydroxybutanoate"          
#> [15] "a-Hydroxybutanoic acid"       "a-Hydroxybutyrate"           
#> [17] "a-Hydroxybutyric acid"        "alpha-Hydroxy-N-butyrate"    
#> [19] "alpha-Hydroxy-N-butyric acid" "alpha-Hydroxybutanoate"      
#> [21] "alpha-Hydroxybutanoic acid"   "alpha-Hydroxybutyrate"       
#> [23] "alpha-Hydroxybutyric acid"    "DL-2-Hydroxybutanoate"       
#> [25] "DL-2-Hydroxybutanoic acid"    "DL-a-Hydroxybutyrate"        
#> [27] "DL-a-Hydroxybutyric acid"     "DL-alpha-Hydroxybutyrate"    
#> [29] "DL-alpha-Hydroxybutyric acid"
#> 
#> $CMP5
#> [1] "2-(8S,9S,13S,14S)-3-Hydroxy-2-methoxy-13-methyl-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]phenanthren-17-one"
#> [2] "2-Hydroxyestrone 2-methyl ether"                                                                                  
#> [3] "2-Methoxy-17-oxoestra-1,3,5(10)-trien-3-ol"                                                                       
#> [4] "2-Methoxy-3-hydroxyestra-1,3,5(10)-trien-17-one"                                                                  
#> [5] "3-Hydroxy-2-methoxy-Estra-1,3,5(10)-trien-17-one"                                                                 
#> [6] "3-Hydroxy-2-methoxyestra-1,3,5(10)-trien-17-one"                                                                  
#> [7] "Methoxy-Estrone"                                                                                                  
#> 
#> $CMP6
#>  [1] "(R)-(-)-b-Hydroxybutyrate"        "(R)-(-)-b-Hydroxybutyric acid"   
#>  [3] "(R)-(-)-beta-Hydroxybutyrate"     "(R)-(-)-beta-Hydroxybutyric acid"
#>  [5] "(R)-3-Hydroxybutanoate"           "(R)-3-Hydroxybutanoic acid"      
#>  [7] "(R)-3-Hydroxybutyrate"            "3-D-Hydroxybutyrate"             
#>  [9] "3-D-Hydroxybutyric acid"          "3-delta-Hydroxybutyrate"         
#> [11] "3-delta-Hydroxybutyric acid"      "BHIB"                            
#> [13] "D-(-)-3-Hydroxybutyrate"          "D-3-Hydroxybutyrate"             
#> [15] "D-3-Hydroxybutyric acid"          "D-beta-Hydroxybutyrate"          
#> [17] "delta-(-)-3-Hydroxybutyrate"      "delta-3-Hydroxybutyrate"         
#> [19] "delta-3-Hydroxybutyric acid"      "delta-beta-Hydroxybutyrate"      
#> 
#> $CMP7
#>  [1] "1-(2-Deoxy-beta-D-erythro-pentofuranosyl)-2,4(1H,3H)-Pyrimidinedione"    
#>  [2] "1-(2-Deoxy-beta-D-ribofuranosyl)-2,4(1H,3H)-Pyrimidinedione"             
#>  [3] "1-(2-Deoxy-beta-delta-erythro-pentofuranosyl)-2,4(1H,3H)-Pyrimidinedione"
#>  [4] "1-(2-Deoxy-beta-delta-ribofuranosyl)-2,4(1H,3H)-Pyrimidinedione"         
#>  [5] "1-(2-Deoxy-D-erythro-pentofuranosyl)uracil"                              
#>  [6] "1-(2-Deoxy-delta-erythro-pentofuranosyl)uracil"                          
#>  [7] "2'-Deoxyuridine"                                                         
#>  [8] "2'-Desoxyuridine"                                                        
#>  [9] "Deoxyribose uracil"                                                      
#> [10] "Desoxyuridine"                                                           
#> [11] "Uracil deoxyriboside"                                                    
#> [12] "Uracil desoxyuridine"                                                    
#> 
#> $CMP8
#> [1] "1-Methyl-2-(3-indolyl)ethylamine"    
#> [2] "2-(1H-Indol-3-yl)-N-methylethanamine"
#> [3] "3-(2-(Methylamino)ethyl)indole"      
#> [4] "Dipterine"                           
#> [5] "Dl-Methyltryptamine"                 
#> [6] "Methyltryptamine"                    
#> [7] "N-Methylindoleethylamine"            
#> [8] "N-Monomethyltryptamine"              
#> [9] "N-Omega-methyltryptamine"            
#> 
#> $CMP9
#> [1] "3-[(2E)-2-[(1R,3aS,7aR)-1-[(1R)-1,5-dimethylhexyl]octahydro-7a-methyl-4H-inden-4-ylidene]ethylidene]-4-methylene-Cyclohexanol"
#> [2] "5,6-trans-Cholecalciferol"                                                                                                    
#> [3] "5,6-trans-Vitamin D3"                                                                                                         
#> [4] "trans-Vitamin D3"                                                                                                             
#> 

## If we provide the optional argument collapse, multiple entries will be
## collapsed.
cmps <- compound_tbl_sdf(fl, collapse = "|")
cmps
#> # A tibble: 9 × 8
#>   compound_id name              inchi inchikey formula exactmass synonyms smiles
#>   <chr>       <chr>             <chr> <chr>    <chr>       <dbl> <chr>    <chr> 
#> 1 HMDB0000001 1-Methylhistidine InCh… BRMWTNU… C7H11N…     169.  1 Methy… "CN1C…
#> 2 HMDB0000002 1,3-Diaminopropa… InCh… XFNJVJP… C3H10N2      74.1 1,3-Dia… "NCCC…
#> 3 HMDB0000005 2-Ketobutyric ac… InCh… TYEYBOS… C4H6O3      102.  2-Ketob… "CCC(…
#> 4 HMDB0000008 2-Hydroxybutyric… InCh… AFENDNX… C4H8O3      104.  (RS)-2-… "CCC(…
#> 5 HMDB0000010 2-Methoxyestrone  InCh… WHEUWNK… C19H24…     300.  2-(8S,9… "[H][…
#> 6 HMDB0000011 (R)-3-Hydroxybut… InCh… WHBMMWS… C4H8O3      104.  (R)-(-)… "C[C@…
#> 7 HMDB0000012 Deoxyuridine      InCh… MXHRCPN… C9H12N…     228.  1-(2-De… "OC[C…
#> 8 HMDB0004370 N-Methyltryptami… InCh… NCIKQJB… C11H14…     174.  1-Methy… "CNCC…
#> 9 HMDB0006719 5,6-trans-Vitami… InCh… QYSXJUF… C27H44O     384.  3-[(2E)… "CC(C…
cmps$synonyms
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP1 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "1 Methylhistidine|1-Methyl histidine|1-Methyl-Histidine|1-Methyl-L-histidine|1-MHis|1-N-Methyl-L-histidine|L-1-Methylhistidine|N1-Methyl-L-histidine|Pi-methylhistidine" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP2 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "1,3-Diamino-N-propane|1,3-Propanediamine|1,3-Propylenediamine|1,3-Trimethylenediamine|3-Aminopropylamine|a,w-Propanediamine|Propane-1,3-diamine|Trimethylenediamine" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP3 
#> "2-Ketobutanoate|2-Ketobutanoic acid|2-Ketobutyrate|2-Oxo-Butanoate|2-Oxo-Butanoic acid|2-Oxo-Butyrate|2-Oxo-Butyric acid|2-Oxo-N-butyrate|2-Oxo-N-butyric acid|2-Oxobutanoate|2-Oxobutanoic acid|2-Oxobutyrate|2-Oxobutyric acid|3-Methylpyruvate|3-Methylpyruvic acid|a-Keto-N-Butyrate|a-Keto-N-Butyric acid|a-Ketobutyrate|a-Ketobutyric acid|a-Oxo-N-butyrate|a-Oxo-N-butyric acid|a-Oxobutyrate|a-Oxobutyric acid|alpha-Keto-N-butyrate|alpha-Keto-N-butyric acid|alpha-Ketobutric acid|alpha-Ketobutyrate|alpha-Ketobutyric acid|alpha-Oxo-N-butyrate|alpha-Oxo-N-butyric acid|alpha-Oxobutyrate|alpha-Oxobutyric acid|Methyl-Pyruvate|Methyl-Pyruvic acid|Propionyl-formate|Propionyl-formic acid" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP4 
#>      "(RS)-2-Hydroxybutyrate|(RS)-2-Hydroxybutyric acid|2-Hydroxy-Butanoate|2-Hydroxy-Butanoic acid|2-Hydroxy-DL-Butyrate|2-Hydroxy-DL-Butyric acid|2-Hydroxy-N-butyrate|2-Hydroxy-N-butyric acid|2-Hydroxybutanoate|2-Hydroxybutanoic acid|2-Hydroxybutyrate|a-Hydroxy-N-butyrate|a-Hydroxy-N-butyric acid|a-Hydroxybutanoate|a-Hydroxybutanoic acid|a-Hydroxybutyrate|a-Hydroxybutyric acid|alpha-Hydroxy-N-butyrate|alpha-Hydroxy-N-butyric acid|alpha-Hydroxybutanoate|alpha-Hydroxybutanoic acid|alpha-Hydroxybutyrate|alpha-Hydroxybutyric acid|DL-2-Hydroxybutanoate|DL-2-Hydroxybutanoic acid|DL-a-Hydroxybutyrate|DL-a-Hydroxybutyric acid|DL-alpha-Hydroxybutyrate|DL-alpha-Hydroxybutyric acid" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP5 
#>                                                                                                                                                                                                                                                                                                                                            "2-(8S,9S,13S,14S)-3-Hydroxy-2-methoxy-13-methyl-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]phenanthren-17-one|2-Hydroxyestrone 2-methyl ether|2-Methoxy-17-oxoestra-1,3,5(10)-trien-3-ol|2-Methoxy-3-hydroxyestra-1,3,5(10)-trien-17-one|3-Hydroxy-2-methoxy-Estra-1,3,5(10)-trien-17-one|3-Hydroxy-2-methoxyestra-1,3,5(10)-trien-17-one|Methoxy-Estrone" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP6 
#>                                                                                                                                                                                                 "(R)-(-)-b-Hydroxybutyrate|(R)-(-)-b-Hydroxybutyric acid|(R)-(-)-beta-Hydroxybutyrate|(R)-(-)-beta-Hydroxybutyric acid|(R)-3-Hydroxybutanoate|(R)-3-Hydroxybutanoic acid|(R)-3-Hydroxybutyrate|3-D-Hydroxybutyrate|3-D-Hydroxybutyric acid|3-delta-Hydroxybutyrate|3-delta-Hydroxybutyric acid|BHIB|D-(-)-3-Hydroxybutyrate|D-3-Hydroxybutyrate|D-3-Hydroxybutyric acid|D-beta-Hydroxybutyrate|delta-(-)-3-Hydroxybutyrate|delta-3-Hydroxybutyrate|delta-3-Hydroxybutyric acid|delta-beta-Hydroxybutyrate" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP7 
#>                                                                                                                                                                                                                          "1-(2-Deoxy-beta-D-erythro-pentofuranosyl)-2,4(1H,3H)-Pyrimidinedione|1-(2-Deoxy-beta-D-ribofuranosyl)-2,4(1H,3H)-Pyrimidinedione|1-(2-Deoxy-beta-delta-erythro-pentofuranosyl)-2,4(1H,3H)-Pyrimidinedione|1-(2-Deoxy-beta-delta-ribofuranosyl)-2,4(1H,3H)-Pyrimidinedione|1-(2-Deoxy-D-erythro-pentofuranosyl)uracil|1-(2-Deoxy-delta-erythro-pentofuranosyl)uracil|2'-Deoxyuridine|2'-Desoxyuridine|Deoxyribose uracil|Desoxyuridine|Uracil deoxyriboside|Uracil desoxyuridine" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP8 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "1-Methyl-2-(3-indolyl)ethylamine|2-(1H-Indol-3-yl)-N-methylethanamine|3-(2-(Methylamino)ethyl)indole|Dipterine|Dl-Methyltryptamine|Methyltryptamine|N-Methylindoleethylamine|N-Monomethyltryptamine|N-Omega-methyltryptamine" 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CMP9 
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "3-[(2E)-2-[(1R,3aS,7aR)-1-[(1R)-1,5-dimethylhexyl]octahydro-7a-methyl-4H-inden-4-ylidene]ethylidene]-4-methylene-Cyclohexanol|5,6-trans-Cholecalciferol|5,6-trans-Vitamin D3|trans-Vitamin D3" 
```

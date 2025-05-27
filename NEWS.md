# CompoundDb version 1.13

## Changes in version 1.13.3

- Add support for additional database tables.

## Changes in version 1.13.2

- Add support for `precScanNum()`.

# CompoundDb version 1.11

## Changes in version 1.11.2

- Import `extractByIndex()` from ProtGenerics.

## Changes in version 1.11.1

- Complete unit test coverage.

# CompoundDb version 1.9

## Changes in version 1.9.5

- Add new `extractByIndex()` method.

## Changes in version 1.9.4

- `compound_tbl_lipidblast` supports now parallel processing and extracts more
  information from MoNA's JSON format (thanks to Prateek Arora for
  contribution).

## Changes in version 1.9.3

- `compound_tbl_lipidblast`: ensure *exactmass* is of type `numeric`.

## Changes in version 1.9.2

- `compound_tbl_lipidblast`: add parameter `n` to support reading and
  processing MoNA json files in sets (chunks) of lines at a time and hence
  reduce memory demand for very large files.

## Changes in version 1.9.1

- Allow `CompDb` to store that database name as alternative to an active
  database connection. This allows to serialize and load an object to/from disk
  (serializing an active database connection would not be possible) . Each call
  to extract data from the database will however open (and close) its own
  connection.

# CompoundDb version 1.7

## Changes in version 1.7.2

- Import method generics from `ProtGenerics`.

## Changes in version 1.7.1

- Adapt script to create CompDb from MassBank to new MassBank database format.

# CompoundDb version 1.3

## Changes in version 1.3.3

- Add `backendBpparam` to define (disable) parallel processing for the
  `MsBackendCompDb` backend.

## Changes in version 1.3.2

- Evaluate validity of the `MsBackendCompDb` using the full test suite from the
  `Spectra` package.

## Changes in version 1.3.2

- Add parameter `nonStop` to `compound_tbl_sdf` that is passed to parameter
  `skipErrors` of `ChemmineR::read.SDFset`.
  Issue [#110](https://github.com/rformassspectrometry/CompoundDb/issues/110)

## Changes in version 1.3.0

- Bioconductor 3.17 developmental version.

# CompoundDb version 1.1

## Changes in version 1.1.6

- `CompDb` tests also for `NA` input.

## Changes in version 1.1.5

- `MsBackendCompDb` always returns `collisionEnergy` as `numeric`.

## Changes in version 1.1.4

- Add script to create a `CompDb` from a MassBank database.

## Changes in version 1.1.3

- Expand vignette with examples to create `CompDb` databases from scratch.

## Changes in version 1.1.2

- Add `insertCompound` and `deleteCompound` functions to add or remove compounds
  from a `CompDb` or `IonDb`.

## Changes in version 1.1.1

- Fix wrong warning message in `deleteIon`.
- Change database data type for internal `ion_id` from `character` to `integer`.

# Version 0.99

## Changes in version 0.99.12

- Add `mass2mz` method for `CompDb` databases.

## Changes in version 0.99.11

- Add `peaksVariables` method.

## Changes in version 0.99.10

- Add parameter `columns` to `peaksData`.

## Changes in version 0.99.9

- Add parameter `dbFile` to `createCompDb` and add an example on how to create
  a `CompDb` database from custom input.

## Changes in version 0.99.8

- Add citation.

## Changes in version 0.99.7

- Add bug reports link to DESCRIPTION.

## Changes in version 0.99.6

- `MsBackendCompDb` extends `Spectra::MsBackendCached` instead of
  `Spectra::MsBackendDataFrame`.

## Changes in version 0.99.5

- No updates, just version bump to cause a new build.

## Changes in version 0.99.4

- Address more comments from @jianhong.

## Changes in version 0.99.3

- Fix `BiocCheck` warnings.

## Changes in version 0.99.2

- Fix `BiocCheck` warnings.

## Changes in version 0.99.1

- Address comments/change requests from @jianhong.

## Changes in version 0.99.0

- Preparing for Bioconductor submission.

# Version 0.9

## Changes in version 0.9.4

- Add `deleteIon` and `deleteSpectra` allowing to delete ions or spectra.

## Changes in version 0.9.3

- `insertIons` supports adding additional database columns.

## Changes in version 0.9.2

- Add `instertSpectra` method to add MS/MS spectra from a `Spectra` object to
  the database.

## Changes in version 0.9.1

- Add `IonDb` constructor methods.
- Expand documentation and examples.
- Add and fix unit tests.

## Changes in version 0.9.0

- Add `IonDb` class as extension of `CompDb` (to allow adding ion information
  to the database) and the functionalities to create such object.
- Add `insertIon` to allow adding new ions to an `IonDb` object
- Add `ionVariables`, `ions` functions to access the ions data in the database.
- Add filters: `IonIdFilter`, `IonAdductFilter`, `IonMzFilter`, `IonRtFilter`.

# Version 0.8

## Changes in version 0.8.1

- Import spectra type (MS level) and precursor type from MoNa.

## Changes in version 0.8.0

- Rename database table name *compound* into *ms_compound* [issue
  #74](https://github.com/EuracBiomedicalResearch/CompoundDb/issues/74).

# Version 0.7

## Changes in version 0.7.0

- Remove `mass2mz` and `mz2mass` function in favour of the functions
  implemented in `MetaboCoreUtils`.

# Version 0.6

## Changes in version 0.6.6

- Import `compounds` method from `ProtGenerics`.

## Changes in version 0.6.5

- Add parameter `onlyValid` to `compound_tbl_sdf` to allow importing of only
  valid elements
  [issue #69](https://github.com/EuracBiomedicalResearch/CompoundDb/issues/69).

## Changes in version 0.6.4

- Add additional filters: `MassFilter`, `FormulaFilter`, `InchiFilter` and
  `InchikeyFilter`.

## Changes in version 0.6.3

- Add `metadata`, `spectraVariables` and `compoundVariables` functions.

## Changes in version 0.6.2

- Support creation of databases without specifying the organism.
- Ensure database columns are mapped correctly to Spectra variable names.


## Changes in version 0.6.1

- Add `SpectrumIdFilter` to support filtering by spectrum IDs.

## Changes in version 0.6.0

- Rename column names: compound_name -> name, mass -> exactmass, inchi_key ->
  inchikey.


# Version 0.5

## Changes in version 0.5.0

- Replace `as.list` with `peaksData`.
- Replace `asDataFrame` with `spectraData`.


# Version 0.4

## Changes in version 0.4.3

- Updated to match new LIPID MAPS field names.


## Changes in version 0.4.2

- Fix bug in `as.list,MsBackendCompDb` which returned a `SimpleList` instead of
  a `list`.


## Changes in version 0.4.0

- Rename method `spectraData` for `MsBackendCompDb` into `asDataFrame`
  (adapting to the changes in `Spectra`).

# Version 0.3

## Changes in version 0.3.2

- Import also smiles from SDF files.


## Changes in version 0.3.1

- Move package Spectra from Depends to Imports


## Changes in version 0.3.0

- Change from MSnbase to RforMassSpectrometry packages (Spectra and
  MsCoreUtils).
- Store MS/MS spectra in two tables, msms_spectrum and msms_spectrum_peak.


# Version 0.2

## Changes in version 0.2.3

- Add instrument and precursor_mz spectra data columns (issue #32).


## Changes in version 0.2.2

- Add adduct information from Jan Stanstrup's commonMZ package.
- Add matchWithPpm function to match numeric values allowing for a small
  difference.
- Add adducts function to retrieve adduct definitions.
- Add mass2mz and mz2mass to convert between mass and m/z for provided adducts.
- Add annotateMz method to annotate m/z values.


## Changes in version 0.2.1

- Change field collision_energy to character to support values from
  MoNa (issue #31).
- Add functions import_mona_sdf and msms_spectra_mona functions to enable
  import of spectrum data from MoNa SDF files (issue #30).
- Add support for MoNa SDF files (issue #30).


## Changes in version 0.2.0

- Add hasMz,Spectrum and hasMz,Spectra methods to look for m/z values within
  spectra (issue #28).
- Add MsmsMzRangeMinFilter and MsmsMzRangeMaxFilter (issue #29).
- Re-use Spectra object from MSnbase.
- Add supportedFilters,CompDb method.


# Version 0.1

## Changes in version 0.1.1

- Add precursorMz, precursorCharge, precursorIntensity, acquisitionNum,
  scanIndex, peaksCount, msLevel, tic, ionCount, collisionEnergy, fromFile,
  polarity, smoothed, isEmpty, centroided and isCentroided methods for
  Spectrum2List.


## Changes in version 0.1.0

- Add expandMzIntensity function.
- Add spectra method to extract spectra from the CompDb database.
- Add functionality to store MS/MS spectra in a CompDb database (m/z and
  intensity values stored as BLOB).
- Add functionality to load MS/MS spectra from HMDB xml files.


# Version 0.0

## Changes in version 0.0.3

- Add CompoundIdFilter and CompoundNameFilter classes and filtering framework.


## Changes in version 0.0.2

- Define CompDb class and all functionality to create CompDb databases.
- createCompDb supports file names as input and create a database including
  annotations from all files.
- Add create-compounddb vignette.

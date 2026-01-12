# Package index

## All functions

- [`CompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`hasMsMsSpectra()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`src_compdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`tables()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`copyCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`dbconn(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`Spectra(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`supportedFilters(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`metadata(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`spectraVariables(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`compoundVariables(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`compounds(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`insertSpectra(`*`<CompDb>`*`,`*`<Spectra>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`deleteSpectra(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`mass2mz(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`insertCompound(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  [`deleteCompound(`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/CompDb.md)
  : Simple compound (metabolite) databases
- [`CompoundIdFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`SpectrumIdFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`NameFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`MsmsMzRangeMinFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`MsmsMzRangeMaxFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`ExactmassFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`FormulaFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`InchiFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`InchikeyFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`IonIdFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`IonAdductFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`IonMzFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  [`IonRtFilter()`](https://rformassspectrometry.github.io/CompoundDb/reference/Filter-classes.md)
  : Filters supported by CompDb and IonDb
- [`ionVariables(`*`<IonDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`ions(`*`<IonDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`insertIon(`*`<IonDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`deleteIon(`*`<IonDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`IonDb(`*`<missing>`*`,`*`<missing>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`IonDb(`*`<CompDb>`*`,`*`<missing>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`IonDb(`*`<character>`*`,`*`<missing>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`IonDb(`*`<DBIConnection>`*`,`*`<missing>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`IonDb(`*`<character>`*`,`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  [`IonDb(`*`<DBIConnection>`*`,`*`<CompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/IonDb.md)
  : IonDb: compound database with additional ion information
- [`MsBackendCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`backendInitialize(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`show(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`peaksData(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`peaksVariables(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`dataStorage(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`` `intensity<-`( ``*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`` `mz<-`( ``*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`spectraData(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`spectraNames(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`` `spectraNames<-`( ``*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`` `[`( ``*`<MsBackendCompDb>`*`,`*`<ANY>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`extractByIndex(`*`<MsBackendCompDb>`*`,`*`<ANY>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`` `$<-`( ``*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`tic(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  [`backendBpparam(`*`<MsBackendCompDb>`*`)`](https://rformassspectrometry.github.io/CompoundDb/reference/MsBackendCompDb.md)
  : CompDb-based MS spectrum backend
- [`addJoinDefinition()`](https://rformassspectrometry.github.io/CompoundDb/reference/addJoinDefinition.md)
  : Expand a CompDb database with additional, related tables
- [`compound_tbl_lipidblast()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_lipidblast.md)
  : Extract compound data from LipidBlast
- [`compound_tbl_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/compound_tbl_sdf.md)
  : Extract compound data from a file in SDF format
- [`createCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
  [`createCompDbPackage()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
  [`make_metadata()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
  [`emptyCompDb()`](https://rformassspectrometry.github.io/CompoundDb/reference/createCompDb.md)
  : Create a CompDb database
- [`expandMzIntensity()`](https://rformassspectrometry.github.io/CompoundDb/reference/expandMzIntensity.md)
  : Expand m/z and intensity values in a data.frame
- [`import_mona_sdf()`](https://rformassspectrometry.github.io/CompoundDb/reference/import_mona_sdf.md)
  : Import compound and spectrum information from MoNa
- [`msms_spectra_hmdb()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_hmdb.md)
  : Import MS/MS spectra from HMDB xml files
- [`msms_spectra_mona()`](https://rformassspectrometry.github.io/CompoundDb/reference/msms_spectra_mona.md)
  : Import MS/MS spectra from MoNa

library("testthat")
library("CompoundDb")
library(RSQLite)

hmdb <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
cmps <- compound_tbl_sdf(hmdb)
metad <- data.frame(name = c("source", "url", "source_version",
                             "source_date", "organism"),
                    value = c("HMDB_test", "http://www.hmdb.ca",
                              "v4", "2017-08-27", "Hsapiens"))
db_file <- createCompDb(cmps, metadata = metad, path = tempdir())
cmp_db <- CompDb(db_file)

dr <- system.file("xml/", package = "CompoundDb")
msms_spctra <- msms_spectra_hmdb(dr)
## spl_ <- as(msms_spctra, "Spectra")

metad2 <- data.frame(name = c("source", "url", "source_version",
                              "source_date", "organism"),
                     value = c("HMDB_spctra", "http://www.hmdb.ca",
                               "v4", "2017-08-27", "Hsapiens"))
db_spctra_file <- createCompDb(cmps, metadata = metad2, path = tempdir(),
                               msms_spectra = msms_spctra)
cmp_spctra_db <- CompDb(db_spctra_file)

cdb <- CompDb(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"))

ions <- data.frame(compound_id = paste0("HMDB000000",
                                        c("1", "1", "2", "2", "5")),
                   ion_adduct = c("A", "B", "B", "C", "D"),
                   ion_mz = c(100, 110, 150, 170, 200),
                   ion_rt = c(50, 60, 100, 110, 90))
ion_db <- IonDb(paste0(tempdir(), "/ion_db.db"), cmp_db, ions)

ion_spctra_db <- IonDb(paste0(tempdir(), "/ion_spctra_db.db"),
                       cmp_spctra_db, ions)

test_check("CompoundDb")

library(Spectra)
be <- backendInitialize(MsBackendCompDb(), cmp_spctra_db)
test_suite <- system.file("test_backends", "test_MsBackend",
                          package = "Spectra")
test_dir(test_suite, stop_on_failure = TRUE)

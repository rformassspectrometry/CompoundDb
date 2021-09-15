library("testthat")
library("CompoundDb")

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

ions <- data.frame(compound_id = paste0("HMDB000000", c("1", "1", "2", "2", "5")),
                  ion_adduct = c("A", "B", "B", "C", "D"),
                  ion_mz = c(100, 110, 150, 170, 200),
                  ion_rt = c(50, 60, 100, 110, 90))
con <- DBI::dbConnect(RSQLite::SQLite(), paste0(tempdir(), "/ion_db.sqlite"))
ion_db <- IonDb(x = con, cdb = cmp_db, ions)

con_spctra <- DBI::dbConnect(RSQLite::SQLite(),
                             paste0(tempdir(), "/ion_spctra_db.sqlite"))
ion_spctra_db <- IonDb(x = con_spctra, cdb = cmp_spctra_db, ions)

test_check("CompoundDb")

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

test_check("CompoundDb")

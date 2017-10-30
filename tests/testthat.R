library("testthat")
library("CompoundDb")

hmdb <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
cmps <- compound_tbl_sdf(hmdb)
metad <- data.frame(name = c("source", "url", "source_version",
    "source_date", "organism"), value = c("HMDB_test", "http://www.hmdb.ca",
                                          "v4", "2017-08-27", "Hsapiens"))
db_file <- createCompDb(cmps, metadata = metad, path = tempdir())
cmp_db <- CompDb(db_file)

test_check("CompoundDb")

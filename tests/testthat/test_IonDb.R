cdb <- CompDb(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"))
test_that("CompDb constructor and low level functions", {
    expect_error(IonDb())
    expect_error(IonDb(3))
    idb <- new("IonDb")

    ions = data.frame(compound_id = c(1, 1, 2, 3, 6, 35),
                      ion_adduct = c("[M+H]+", "[M+Na]+", "[M+Na]+",
                                     "[M+Na]+", "[M+2H]2+", "[M+H-NH3]+"),
                      ion_mz = c(179.0703, 201.0522, 201.0522,
                                 201.0522, 253.66982, 312.0390),
                      ion_rt = 1:6)
    icon <- dbConnect(RSQLite::SQLite(), paste0(tempdir(), "/idb.sqlite"))
    idb <- IonDb(x = icon, cdb = cdb, ions)
    
    # expect_true(!is.null(.dbconn(icon)))
    # expect_true(.validIonDb(dbconn(icon)))
    expect_equal(tables(idb)$ms_ion, c("ion_id", "compound_id", "ion_adduct",
                                       "ion_mz", "ion_rt"))
    expect_equal(tables(idb)[-2], tables(cdb))

    # Connection to an already existing database with ion information
    idb2 <- IonDb(x = icon)
    expect_equal(tables(idb2), tables(idb))
    dbDisconnect(icon)
})

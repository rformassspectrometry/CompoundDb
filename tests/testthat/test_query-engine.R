test_that(".reduce_tables works", {
    tabs <- list(gene = c("gene_name", "gene_id", "redundant_field"),
                 tx = c("tx_id", "gene_id", "tx_name", "redundant_field",
                        "tx_start"),
                 exon = c("exon_id", "redundant_field", "tx_id"))

    res <- .reduce_tables(tabs, columns = c("tx_id", "gene_id"))
    expect_equal(length(res), 1)
    expect_equal(res, list(tx = c("tx_id", "gene_id")))

    res <- .reduce_tables(tabs, columns = c("gene_id", "gene_name"))
    expect_equal(length(res), 1)
    expect_equal(res, list(gene = c("gene_id", "gene_name")))

    res <- .reduce_tables(tabs, columns = c("gene_name", "exon_id",
                                            "redundant_field",
                                            "tx_id"))
    expect_equal(length(res), 2)
    expect_equal(res, list(exon = c("exon_id", "redundant_field", "tx_id"),
                           gene = "gene_name"))
})

test_that(".prefix_columns works", {
    res <- .prefix_columns(list(a = c("b", "c")))
    expect_equal(res, list(a = c("a.b", "a.c")))
    res <- .prefix_columns(list(a = c("b", "c"),
                                b = c("d", "e")))
    expect_equal(res, list(a = c("a.b", "a.c"),
                           b = c("b.d", "b.e")))
})

test_that(".add_join_tables works", {
    expect_equal(.add_join_tables(c("a", "b")), c("a", "b"))
    expect_equal(.add_join_tables(c("msms_spectrum_peak")),
                 c("msms_spectrum_peak"))
    expect_equal(.add_join_tables(c("msms_spectrum", "compound")),
                 c("msms_spectrum", "compound"))
    ## expect_equal(.add_join_tables(c("msms_spectrum_peak", "compound")),
    ##              c("msms_spectrum_peak", "compound", "msms_spectrum_metadata"))
})

test_that(".where works", {
    expect_equal(.where(), NULL)
    expect_equal(.where(CompoundIdFilter("a")), " where compound_id = 'a'")
    expect_equal(.where(CompoundIdFilter("a"), list(table = "compound_id")),
                 " where table.compound_id = 'a'")
})

test_that(".order works", {
    expect_equal(.order(), NULL)
    expect_equal(.order("a"), " order by a")
})

test_that(".from works", {
    expect_equal(.from("compound"), " from compound")
    expect_error(.from(c("a", "b")))
    expect_equal(.from(c("compound", "synonym")),
                 paste0(" from compound left outer join synonym on (compound.",
                        "compound_id=synonym.compound_id)"))
    expect_equal(.from(c("synonym", "compound")),
                 paste0(" from synonym left outer join compound on (compound.",
                        "compound_id=synonym.compound_id)"))
})

test_that(".select works", {
    expect_equal(.select(c("a", "b"), distinct = FALSE), "select a,b")
    expect_equal(.select(c("a", "b")), "select distinct a,b")
    expect_error(.select())
})

test_that(".build_query_CompDb works", {
    expect_error(.build_query_CompDb())
    expect_error(.build_query_CompDb(columns = c("compound_id", "inchi")))
    expect_error(.build_query_CompDb(cmp_db))
    res <- .build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"))
    expect_equal(res, paste0("select distinct compound.compound_id,compound.",
                             "inchi from compound"))
    expect_error(.build_query_CompDb(
        cmp_db, columns = c("od", "inchi")))
    res <- .build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"), order = "something")
    expect_equal(res, paste0("select distinct compound.compound_id,compound.",
                             "inchi from compound order by something"))
    res <- .build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"),
        filter = ~ compound_id == "a")
    expect_equal(res,
                 paste0("select distinct compound.compound_id,compound.",
                        "inchi from compound where compound.compound_id = 'a'"))
    res <- .build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"),
        filter = ~ compound_id == "a" | compound_name != "b")
    expect_equal(res, paste0("select distinct compound.compound_id,compound.",
                             "inchi,compound.compound_name from compound ",
                             "where (compound.compound_id = 'a' or ",
                             "compound.compound_name != 'b')"))
    expect_error(.build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"),
        filter = ~ compound_id == "a" | gene_id != "b"))
    ## Include columns from msms_spectrum_*
    expect_error(.build_query_CompDb(cmp_db,
                                     columns = c("compound_id", "intensity")))
    ## compound, msms_spectrum
    expect_error(.build_query_CompDb(cmp_spctra_db, start_from = "compound",
                                     columns = c("compound_id", "intensity")),
                 "can not be joined")
    res <- .build_query_CompDb(
        cmp_spctra_db, columns = c("compound_id", "intensity", "splash"))
    expect_equal(
        res, paste0("select distinct msms_spectrum.compound_id,msms_spectrum.s",
                    "plash,msms_spectrum_peak.intensity from msms_spectrum lef",
                    "t outer join msms_spectrum_peak on (msms_spectrum.spectru",
                    "m_id=msms_spectrum_peak.spectrum_id)"))
    res <- .build_query_CompDb(cmp_spctra_db,
                               columns = c("splash", "inchi"))
    expect_equal(
        res, paste0("select distinct compound.inchi,msms_spectrum.splash from ",
                    "compound left outer join msms_spectrum on (compound.compo",
                    "und_id=msms_spectrum.compound_id)"))
    res <- .build_query_CompDb(cmp_spctra_db, start_from = "msms_spectrum",
                               columns = c("splash", "inchi"),
                               filter = ~ compound_id == "a")
    expect_equal(
        res, paste0("select distinct msms_spectrum.compound_id,msms_spectrum.s",
                    "plash,compound.inchi from msms_spectrum left outer join c",
                    "ompound on (compound.compound_id=msms_spectrum.compound_i",
                    "d) where msms_spectrum.compound_id = 'a'"))
    ## msms_spectrum, synonym
    expect_error(.build_query_CompDb(
        cmp_spctra_db, columns = c("synonym", "mz")), "can not be")
    res <- .build_query_CompDb(
        cmp_spctra_db, columns = c("synonym", "mz", "polarity"))
    expect_equal(
        res, paste0("select distinct msms_spectrum.polarity,msms_spectrum_peak",
                    ".mz,synonym.synonym from msms_spectrum left outer join sy",
                    "nonym on (msms_spectrum.compound_id=synonym.compound_id) ",
                    "left outer join msms_spectrum_peak on (msms_spectrum.spec",
                    "trum_id=msms_spectrum_peak.spectrum_id)"))
})

test_that(".join_tables works", {
    res <- .join_tables(c("compound", "synonym"))
    expect_equal(res, paste0("compound left outer join synonym on ",
                             "(compound.compound_id=synonym.compound_id)"))
    expect_equal(.join_tables("compound"), "compound")
    expect_error(.join_tables(c("cmps", "synonym")))
    res <- .join_tables(c("compound", "msms_spectrum"))
    expect_equal(
        res, paste0("compound left outer join msms_spectrum on (compound.compo",
                    "und_id=msms_spectrum.compound_id)"))
    res <- .join_tables(c("msms_spectrum", "compound"))
    expect_equal(
        res, paste0("msms_spectrum left outer join compound on (compound.",
                    "compound_id=msms_spectrum.compound_id)"))
    res <- .join_tables(c("synonym", "msms_spectrum"))
    expect_equal(
        res, paste0("synonym left outer join msms_spectrum on (msms_spectrum.c",
                    "ompound_id=synonym.compound_id)"))
    ## res <- .join_tables(c("compound", "msms_spectrum_metadata"))
    ## expect_equal(res, paste0("compound left outer join msms_spectrum_metadata ",
    ##                          "on (compound.compound_id=msms_spectrum_metadata.",
    ##                          "compound_id)"))
    ## ## Same query but starting from other table.
    ## res <- .join_tables(c("msms_spectrum_metadata", "compound"))
    ## expect_equal(res, paste0("msms_spectrum_metadata left outer join compound ",
    ##                          "on (compound.compound_id=msms_spectrum_metadata.",
    ##                          "compound_id)"))
    ## res <- .join_tables(c("msms_spectrum_peak", "compound"))
    ## expect_equal(res, paste0("msms_spectrum_peak left outer join msms_spectrum",
    ##                          "_metadata on (msms_spectrum_metadata.spectrum_id",
    ##                          "=msms_spectrum_peak.spectrum_id) left outer join",
    ##                          " compound on (compound.compound_id=msms_spectrum",
    ##                          "_metadata.compound_id)"))
    ## res <- .join_tables(c("synonym", "msms_spectrum_peak"))
    ## expect_equal(res, paste0("synonym left outer join msms_spectrum_metadata",
    ##                          " on (msms_spectrum_metadata.compound_id=synonym",
    ##                          ".compound_id) left outer join msms_spectrum_",
    ##                          "peak on (msms_spectrum_metadata.spectrum_id=",
    ##                          "msms_spectrum_peak.spectrum_id)"))
})

test_that(".reduce_tables_start_from works", {
    tabs <- list(
        compound = c("compound_id", "compound_name", "red_field"),
        spectrum = c("spectrum_id", "compound_id"),
        other_tab = c("compound_id", "red_field"))
    res <- .reduce_tables_start_from(tabs, c("compound_id"))
    expect_equal(res, list(compound = "compound_id"))
    res <- .reduce_tables_start_from(tabs, c("compound_id", "red_field"))
    expect_equal(res, list(compound = c("compound_id", "red_field")))
    res <- .reduce_tables_start_from(tabs, c("compound_id", "red_field"),
                                     start_from = "other_tab")
    expect_equal(res, list(other_tab = c("compound_id", "red_field")))
    res <- .reduce_tables_start_from(tabs, c("compound_id", "red_field"),
                                     start_from = "spectrum")
    expect_equal(res, list(spectrum = "compound_id", compound = "red_field"))
    expect_warning(res <- .reduce_tables_start_from(
                       tabs, c("compound_name", "red_field"),
                       start_from = "spectrum"))
    expect_equal(res, list(compound = c("compound_name", "red_field")))
    expect_error(.reduce_tables_start_from(tabs,
                                           c("compound_name", "red_field"),
                                           start_from = "spectrum_bla") )
})

test_that(".deserialize_mz_intensity works", {
    df <- data.frame(cola = 1:3, colb = c("a", "b", "c"),
                     stringsAsFactors = FALSE)
    df$mz <- list(serialize(1:3, NULL), serialize(5:10, NULL),
                  serialize(3:9, NULL))
    df_ds <- .deserialize_mz_intensity(df)
    expect_equal(df[, 1:2], df_ds[, 1:2])
    expect_equal(df_ds$mz[[1]], 1:3)
    expect_equal(df_ds$mz[[2]], 5:10)
    expect_equal(df_ds$mz[[3]], 3:9)

    ## ## Now from the database...
    ## library(RSQLite)
    ## res <- dbGetQuery(dbconn(cmp_spctra_db), "select * from msms_spectrum")
    ## res <- .deserialize_mz_intensity(res)
    ## expect_true(is.numeric(res$mz[[1]]))
    ## expect_true(is.numeric(res$mz[[2]]))
    ## expect_true(is.numeric(res$mz[[3]]))
    ## expect_true(is.numeric(res$intensity[[1]]))
    ## expect_true(is.numeric(res$intensity[[2]]))
    ## expect_true(is.numeric(res$intensity[[3]]))
})

test_that(".fetch_data works", {
    clmns <- c("compound_id", "compound_name", "inchi")
    res <- .fetch_data(cmp_db, clmns)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), clmns)

    ## With filter.
    res <- .fetch_data(cmp_db, clmns,
                       filter = ~ compound_id == "HMDB0000002")
    expect_equal(colnames(res), clmns)
    expect_true(all(res$compound_id == "HMDB0000002"))

    ## MS/MS spectra
    res <- .fetch_data(cmp_spctra_db,
                       columns = c("mz", "compound_name", "polarity"))
    expect_equal(colnames(res), c("polarity", "spectrum_id", "mz",
                                  "compound_name"))
    expect_true(is.numeric(res$mz[[1]]))
    res <- CompoundDb:::.fetch_data(cmp_spctra_db,
                       columns = c("compound_name", "spectrum_id",
                                   "compound_id"),
                       filter = ~ compound_id == "HMDB0000001")
    expect_equal(colnames(res), c("compound_name", "compound_id",
                                  "spectrum_id"))
    expect_true(nrow(res) == 2)
})

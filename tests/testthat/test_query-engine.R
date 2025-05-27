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
    expect_equal(.add_join_tables(c("msms_spectrum", "ms_compound")),
                 c("msms_spectrum", "ms_compound"))
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
    expect_equal(.from("ms_compound"), " from ms_compound")
    expect_error(.from(c("a", "b")))
    expect_equal(
        .from(c("ms_compound", "synonym")),
        paste0(" from ms_compound left outer join synonym on (ms_compound.",
               "compound_id=synonym.compound_id)"))
    expect_equal(
        .from(c("synonym", "ms_compound")),
        paste0(" from synonym left outer join ms_compound on (ms_compound.",
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
    expect_equal(res,
                 paste0("select distinct ms_compound.compound_id,ms_compound.",
                        "inchi from ms_compound"))
    expect_error(.build_query_CompDb(
        cmp_db, columns = c("od", "inchi")))
    res <- .build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"), order = "something")
    expect_equal(res,
                 paste0("select distinct ms_compound.compound_id,ms_compound.",
                        "inchi from ms_compound order by something"))
    res <- .build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"),
        filter = ~ compound_id == "a")
    expect_equal(
        res,
        paste0("select distinct ms_compound.compound_id,ms_compound.",
               "inchi from ms_compound where ms_compound.compound_id = 'a'"))
    res <- .build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"),
        filter = ~ compound_id == "a" | name != "b")
    expect_equal(res,
                 paste0("select distinct ms_compound.compound_id,ms_compound.",
                        "inchi,ms_compound.name from ms_compound ",
                        "where (ms_compound.compound_id = 'a' or ",
                        "ms_compound.name != 'b')"))
    expect_error(.build_query_CompDb(
        cmp_db, columns = c("compound_id", "inchi"),
        filter = ~ compound_id == "a" | gene_id != "b"))
    ## Include columns from msms_spectrum_*
    expect_error(.build_query_CompDb(cmp_db,
                                     columns = c("compound_id", "intensity")))
    ## compound, msms_spectrum
    expect_error(.build_query_CompDb(cmp_spctra_db, start_from = "ms_compound",
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
        res,
        paste0("select distinct ms_compound.inchi,msms_spectrum.splash from ",
               "ms_compound left outer join msms_spectrum on (ms_compound.",
               "compound_id=msms_spectrum.compound_id)"))
    res <- .build_query_CompDb(cmp_spctra_db, start_from = "msms_spectrum",
                               columns = c("splash", "inchi"),
                               filter = ~ compound_id == "a")
    expect_equal(
        res,
        paste0("select distinct msms_spectrum.compound_id,msms_spectrum.s",
               "plash,ms_compound.inchi from msms_spectrum left outer join ",
               "ms_compound on (ms_compound.compound_id=msms_spectrum.",
               "compound_id) where msms_spectrum.compound_id = 'a'"))
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
    res <- .join_tables(c("ms_compound", "synonym"))
    expect_equal(res, paste0("ms_compound left outer join synonym on ",
                             "(ms_compound.compound_id=synonym.compound_id)"))
    expect_equal(.join_tables("ms_compound"), "ms_compound")
    expect_error(.join_tables(c("cmps", "synonym")))
    res <- .join_tables(c("ms_compound", "msms_spectrum"))
    expect_equal(
        res,
        paste0("ms_compound left outer join msms_spectrum on (ms_compound.",
               "compound_id=msms_spectrum.compound_id)"))
    res <- .join_tables(c("msms_spectrum", "ms_compound"))
    expect_equal(
        res,
        paste0("msms_spectrum left outer join ms_compound on (ms_compound.",
               "compound_id=msms_spectrum.compound_id)"))
    res <- .join_tables(c("synonym", "msms_spectrum"))
    expect_equal(
        res, paste0("synonym left outer join msms_spectrum on (msms_spectrum.c",
                    "ompound_id=synonym.compound_id)"))
})

test_that(".reduce_tables_start_from works", {
    tabs <- list(
        ms_compound = c("compound_id", "name", "red_field"),
        spectrum = c("spectrum_id", "compound_id"),
        other_tab = c("compound_id", "red_field"))
    res <- .reduce_tables_start_from(tabs, c("compound_id"))
    expect_equal(res, list(ms_compound = "compound_id"))
    res <- .reduce_tables_start_from(tabs, c("compound_id", "red_field"))
    expect_equal(res, list(ms_compound = c("compound_id", "red_field")))
    res <- .reduce_tables_start_from(tabs, c("compound_id", "red_field"),
                                     start_from = "other_tab")
    expect_equal(res, list(other_tab = c("compound_id", "red_field")))
    res <- .reduce_tables_start_from(tabs, c("compound_id", "red_field"),
                                     start_from = "spectrum")
    expect_equal(res, list(spectrum = "compound_id", ms_compound = "red_field"))
    expect_warning(res <- .reduce_tables_start_from(
                       tabs, c("name", "red_field"),
                       start_from = "spectrum"))
    expect_equal(res, list(ms_compound = c("name", "red_field")))
    expect_error(.reduce_tables_start_from(tabs,
                                           c("name", "red_field"),
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

    df$intensity <- list(serialize(5:7, NULL), serialize(15:20, NULL),
                  serialize(13:19, NULL))
    df_ds <- .deserialize_mz_intensity(df)
    expect_equal(df[, 1:2], df_ds[, 1:2])
    expect_equal(df_ds$intensity[[1]], 5:7)
    expect_equal(df_ds$intensity[[2]], 15:20)
    expect_equal(df_ds$intensity[[3]], 13:19)

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
    clmns <- c("compound_id", "name", "inchi")
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
                       columns = c("mz", "name", "polarity"))
    expect_equal(colnames(res), c("polarity", "spectrum_id", "mz",
                                  "name"))
    expect_true(is.numeric(res$mz[[1]]))
    res <- .fetch_data(cmp_spctra_db,
                       columns = c("name", "spectrum_id",
                                   "compound_id"),
                       filter = ~ compound_id == "HMDB0000001")
    expect_equal(colnames(res), c("name", "compound_id",
                                  "spectrum_id"))
    expect_true(nrow(res) == 2)
})

test_that(".path_nodes works", {
    expect_equal(.path_nodes(), Inf)
    expect_equal(.path_nodes(1:3), 3)
    expect_equal(.path_nodes(c("a", "b")), 2)
})

test_that(".shortest_path works", {
    ## Construct graph
    g <- list(a = c("b", "c"),
              b = c("a", "e", "f", "g"),
              c = c("a", "d"),
              d = c("c"),
              e = c("b", "g"),
              f = c("b"),
              g = c("b", "e"))

    res <- .shortest_path(g, "a", "b")
    expect_equal(res, c("a", "b"))
    res <- .shortest_path(g, "a", "d")
    expect_equal(res, c("a", "c", "d"))
    res <- .shortest_path(g, "a", "g")
    expect_equal(res, c("a", "b", "g"))

    expect_equal(.shortest_path(g, "z"), NULL)

    g <- list(a = c("b", "c"),
              b = c("a", "c"),
              c = c("a", "b"),
              d = "e",
              e = "d")
    res <- .shortest_path(g, "a", "c")
    expect_equal(res, c("a", "c"))
    expect_equal(.shortest_path(g, "a", "e"), NULL)
})

test_that(".table_to_graph works", {
    m <- rbind(c("a", "b"),
               c("a", "c"),
               c("d", "a"),
               c("b", "e"),
               c("f", "g"))
    res <- .table_to_graph(m)
    expect_equal(names(res), unique(as.vector(m)))
    expect_equal(res[["a"]], c("d", "b", "c"))
    expect_equal(res[["d"]], c("a"))
    expect_equal(res[["b"]], c("a", "e"))
    expect_equal(res[["f"]], c("g"))
    expect_equal(res[["c"]], c("a"))
    expect_equal(res[["e"]], c("b"))
    expect_equal(res[["g"]], c("f"))

    ## res <- .table_to_graph(.JOINS[, 1:2])
    ## expect_equal(res[["ms_compound"]], c("ms_ion", "synonym", "msms_spectrum"))
    ## expect_equal(res[["ms_ion"]], c("ms_compound", "msms_spectrum"))
    ## expect_equal(res[["msms_spectrum"]], c("ms_compound", "ms_ion", "synonym",
    ##                                        "msms_spectrum_peak"))
    ## expect_equal(res[["synonym"]], c("ms_compound", "msms_spectrum"))
    ## expect_equal(res[["msms_spectrum_peak"]], c("msms_spectrum"))
})

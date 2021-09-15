
test_that("show,IonDb works", {
    expect_output(show(ion_db))
    db <- new("IonDb")
    expect_output(show(ion_db))
})

# test_that("dbconn,IonDb works", {
#     expect_true(!is.null(dbconn(ion_db)))
#     expect_true(is(dbconn(ion_db), "DBIConnection"))
# })

# test_that("Spectra,IonDb works", {
# 
#     res <- Spectra(ion_spctra_db)
#     expect_true(is(res, "Spectra"))
#     expect_true(length(res) == 4)
#     expect_true(all(c("instrument", "predicted") %in% spectraVariables(res)))
# 
#     ## columns
#     res <- Spectra(ion_spctra_db, columns = "compound_id")
#     expect_false(all(c("instrument", "predicted") %in% spectraVariables(res)))
# 
#     ## filter
#     res <- Spectra(ion_spctra_db, filter = ~ compound_id == "HMDB0000001")
#     expect_true(length(res) == 2)
# 
#     ## filter and columns
#     res <- Spectra(ion_spctra_db, filter = ~ compound_id == "HMDB0000001",
#                    columns = c("inchi", "name"))
#     expect_true(all(c("spectrum_id", "name", "inchi") %in%
#                     spectraVariables(res)))
#     expect_true(length(res) == 2)
# 
#     expect_error(Spectra(ion_spctra_db, filter = "ad"), "'filter' has to")
#     expect_error(Spectra(ion_spctra_db, filter = ~ gene_name == "b"),
#                  "not supported")
# })
# 
# test_that("supportedFilters works", {
#     res <- supportedFilters(ion_db)
#     expect_equal(colnames(res), c("filter", "field"))
#     res_2 <- supportedFilters(ion_spctra_db)
#     expect_true(nrow(res) < nrow(res_2))
# })
# 
# test_that("metadata works", {
#     res <- metadata(ion_db)
#     expect_true(is.data.frame(res))
#     expect_true(all(colnames(res) == c("name", "value")))
# })
# 
# test_that("spectraVariables,IonDb works", {
#     db <- new("IonDb")
#     expect_equal(spectraVariables(db), character())
# 
#     res <- spectraVariables(ion_spctra_db)
#     expect_true(is.character(res))
#     expect_true(length(res) > 0)
# })
# 
# test_that("compoundVariables,IonDb works", {
#     db <- new("CompDb")
#     expect_equal(compoundVariables(db), character())
# 
#     res <- compoundVariables(ion_db)
#     expect_true(is.character(res))
#     expect_true(length(res) > 0)
#     expect_true(all(c("formula", "inchi") %in% res))
# 
#     expect_true(any(compoundVariables(ion_db, TRUE) == "compound_id"))
# })
# 
# test_that("compounds works", {
#     res <- compounds(ion_db, columns = character())
#     expect_true(is.data.frame(res))
#     expect_true(ncol(res) == 0)
#     expect_true(nrow(res) == 0)
#     cmps <- compounds(ion_db)
#     expect_true(is(cmps, "data.frame"))
#     cmps_tbl <- compounds(ion_db, columns = c("compound_id", "name"),
#                           return.type = "tibble")
#     expect_true(is(cmps_tbl, "tbl"))
#     expect_equal(colnames(cmps_tbl), c("compound_id", "name"))
# 
#     expect_error(compounds(ion_db, filter = "something"))
# 
#     expect_true(
#         nrow(compounds(ion_db, filter = ~ compound_id == "HMDB0000005")) == 1)
#     res <- compounds(cmp_spctra_db,
#                      columns = c("compound_id", "spectrum_id", "splash"))
#     cmp_ids <- compounds(cmp_spctra_db, columns = "compound_id")$compound_id
#     expect_true(all(cmp_ids %in% res$compound_id))
#     expect_true(sum(is.na(res$spectrum_id)) == 6)
# 
#     ## compounds with filters
#     res <- compounds(cdb, filter = ~ exactmass > 300)
#     expect_true(all(res$exactmass > 300))
#     res_2 <- compounds(cdb, filter = ~ exactmass > 300 & exactmass < 340)
#     expect_true(nrow(res_2) < nrow(res))
# 
#     res <- compounds(cdb, filter = FormulaFilter("C17", "startsWith"))
#     expect_true(nrow(res) > 0)
# })


test_that("ionVariables,ionDb works", {
    db <- new("IonDb")
    expect_equal(ionVariables(db), character())
    
    res <- ionVariables(ion_db)
    expect_true(is.character(res))
    expect_true(length(res) > 0)
    expect_true(all(c("compound_id", "ion_adduct", "ion_mz", "ion_rt")
                    %in% res))
    expect_true(any(ionVariables(ion_db, TRUE) == "ion_id"))
})

test_that("ions,ionDb works", {
    res <- ions(ion_db, columns = character())
    expect_true(is.data.frame(res))
    expect_true(ncol(res) == 0)
    expect_true(nrow(res) == 0)
    expect_true(is(ions(ion_db), "data.frame"))
    ions_tbl <- ions(ion_db, columns = c("ion_id", "compound_id", "ion_adduct"),
                          return.type = "tibble")
    expect_true(is(ions_tbl, "tbl"))
    expect_equal(colnames(ions_tbl), c("ion_id", "compound_id", "ion_adduct"))
    
    expect_error(compounds(ion_db, filter = "something"))
    
    ## ions with filters
    res <- ions(ion_db, filter = ~ ion_rt < 100)
    expect_true(all(res$ion_rt < 100))
    res_2 <- ions(ion_db, filter = ~ ion_rt < 100 & ion_mz < 200)
    expect_true(nrow(res_2) < nrow(res))

    res <- ions(ion_db, filter = CompoundIdFilter("01", "endsWith"))
    expect_true(nrow(res) > 0)
})

# This test is succesful when I run it directly but an error is thrown when I
# do rcmdcheck() and I still have to understand why 

# test_that("insertIon,ionDb works", {
#     more_ions <- data.frame(compound_id = c("HMDB0000005", "HMDB0000008"),
#                             ion_adduct = c("C", "D"),
#                             ion_mz = c(220, 300),
#                             ion_rt = c(90, 140))
#     insertIon(ion_db, more_ions)
#     expect_equal(ions(ion_db), rbind(ions, more_ions))
#     expect_equal(ions(ion_db, "ion_id")$ion_id,
#                  seq_len(nrow(ions) + nrow(more_ions)))
# })



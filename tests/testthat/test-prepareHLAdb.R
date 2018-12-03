### Unit tests for functions in prepareHLAdb.R file

library(HLAClustRView)

directory <- system.file("extdata", package = "HLAClustRView")

### Tests parseHLADbAlignment() results

context("parseHLADbAlignment() results")

test_that("parseHLADbAlignment() must retun an error when directory is a number", {
    message <- "The hlaDbPath parameter must by a character string"
    expect_error(parseHLADbAlignment(hlaDbPath = 33, seqType = "nuc"), message)
})

test_that("parseHLADbAlignment() must retun an error when directory does not exist", {
    path_to_test <- "./3_21_33eddger_39kdf_23sk2-2221fvd_qqqqqqqqqq_2"
    if (!dir.exists(path_to_test)) {
        message <- "The hlaDbPath parameter must by a valid directory"
        expect_error(parseHLADbAlignment(hlaDbPath = path_to_test, seqType = "nuc"), message)
    }
})

test_that("parseHLADbAlignment() must retun an error when no file in directory", {
    message <- "There must be at least one alignment file in the hlaDbPath directory"
    expect_error(parseHLADbAlignment(hlaDbPath = directory, seqType = "nuc"), message)
})

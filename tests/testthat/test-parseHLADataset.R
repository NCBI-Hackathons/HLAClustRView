### Unit tests for functions in parseHLADb.R file

library(HLAClustRView)


data("demoHLADataset")


### Tests readHLADataset() function

context("Test for readHLADataset() function")

test_that("readHLADataset() must return an error when hlaFilePath is a number", {
    message <- "The hlaFilePath parameter must by a character string"
    expect_error(readHLADataset(hlaFilePath = 33), message)
})

test_that("readHLADataset() must return an error when hlaFilePath is a number", {
    message <- "The hlaFilePath parameter must by a character string"
    expect_error(readHLADataset(hlaFilePath = TRUE), message)
})

test_that("parseHLADbAlignment() must return an error when hlaFilePath does not exist", {
    path_to_test <- "./3_21_33eddger_39kdf_23sk2-2221fvd_qqqqqqqqqq_2.txt"
    if (!dir.exists(path_to_test)) {
        message <- paste0("The file '", path_to_test, "' must be a valid text file")
        expect_error(readHLADataset(hlaFilePath = path_to_test), message)
    }
})

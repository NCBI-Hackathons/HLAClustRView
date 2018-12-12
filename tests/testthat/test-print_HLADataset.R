### Unit tests for functions in prepareHLAdb.R file

library(HLAClustRView)


data("demoHLADataset")


### Tests print() for HLADataset class

context("print() HLADataset class")

test_that("print() for HLADataset object must return identical object", {

    result <- print(demoHLADataset)

    expect_equal(result, demoHLADataset)
})

test_that("print() for HLADataset object without data must return identical object", {

    demo <- demoHLADataset
    demo[["data"]] <- NULL

    result <- print(demo)

    expect_equal(result, demo)
})

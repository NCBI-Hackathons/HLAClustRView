### Unit tests for functions in print.HLAMetric.R file

library(HLAClustRView)


data("demoHLADataset")


### Tests print() for HLAMetric class

context("print() HLAMetric class")

test_that("print() for HLAMetric object must return identical object", {

    expeced <- calculateHamming(demoHLADataset)

    result <- print(expeced)

    expect_equal(result$dist, expeced$dist)
})

test_that("print() for HLAMetric object without distance must return identical object", {

    demo <- calculateHamming(demoHLADataset)
    demo[["dist"]] <- NULL

    result <- print(demo)

    expect_equal(result$dist, demo$dist)
})

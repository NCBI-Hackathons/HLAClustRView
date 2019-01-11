### Unit tests for functions in createCluster.R file

library(HLAClustRView)

### Tests hclustHLA() results

context("hclustHLA() results")

test_that("hclustHLA() must failed when wrong dist input 01", {
    expect_error(hclustHLA(dist = "hello"))
})

test_that("hclustHLA() must failed when wrong dist input 02", {
    expect_error(hclustHLA(dist = 33))
})

### Unit tests for functions in createCluster.R file

library(HLAClustRView)

### Tests createCluster() results

context("createCluster() results")

test_that("applyCNPmask() must failed when wrong dist input", {
    expect_error(createCluster(dist = "hello"))
})

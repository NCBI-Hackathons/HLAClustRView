### Unit tests for functions in splitTyping.R file

library(HLAClustRView)


context("Test splitTyping() function")

test_that("splitTyping() must return the good result 01", {
    input <- c("G*01:04:01:02", "DQA1*01:01:01:03", "DPA1*01:03:01:18Q",
               "DRB3*01:01:02:02")

    result <- HLAClustRView:::splitTyping(input)

    expected <- matrix(c("G", "01", "04", "01", "02", NA,
                         "DQA1", "01", "01", "01", "03", NA,
                         "DPA1", "01", "03", "01", "18", "Q",
                         "DRB3", "01", "01", "02", "02", NA),
                        nrow = 4, byrow = TRUE)

    expect_identical(result, expected)
})

test_that("splitTyping() must return the good result 02", {
    input <- c("DRA*01:01:02", "E*01:10")

    result <- HLAClustRView:::splitTyping(input)

    expected <- matrix(c("DRA", "01", "01", "02", NA, NA,
                         "E", "01", "10", NA, NA, NA),
                       nrow = 2, byrow = TRUE)

    expect_identical(result, expected)
})

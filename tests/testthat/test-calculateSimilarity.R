### Unit tests for functions in calculateSimilarity.R file

library(HLAClustRView)
library(dplyr)

data("example_sample_pair_data")
data("demoHLADataset")


context("Test hamming_distance_digit1() function")


demo1 <- data.frame(SampleName=c("s1", "s1", "s2", "s2"),
                    AlleleName=c(1, 2, 1, 2), AlleleGroup=c(1, 3, 1, 5))

demo2 <- data.frame(SampleName=c("s1", "s1", "s2", "s2"),
                    AlleleName=c(1, 2, 1, 2), AlleleGroup=c(1, 1, 1, 1))

# test function hamming_distance_digit1 -----------------------------------
test_that("correct Hamming Value", {
    d1 <- dplyr::as_tibble(demo1)
    d2 <- dplyr::as_tibble(demo2)
    expect_equal(HLAClustRView:::hamming_distance_digit1(d1)$distance, 1)
    expect_equal(HLAClustRView:::hamming_distance_digit1(d2)$distance, 0)
    ## returns NA when distance is the same for both allels
    expect_true(is.na(HLAClustRView:::hamming_distance_digit1(d2)$same_allele))
})

test_that("correct output format", {
    d1 <- dplyr::as_tibble(demo1)
    d2 <- dplyr::as_tibble(demo2)
    expect_equal(nrow(HLAClustRView:::hamming_distance_digit1(d2)), 1)
    expect_equal(ncol(HLAClustRView:::hamming_distance_digit1(d2)), 2)

})

# test function sample_pair_distance --------------------------------------

context("Test sample_pair_distance() function")

test_that("sample_pair_distance() should return good result", {
    res <- HLAClustRView:::sample_pair_distance(example_sample_pair_data)
    expect_equal(res$HammingDistance, 3)
    # keep only genes that are in both samples
    expect_true(all(c("A", "C") %in% tidyr::unnest(res)$GeneName))

    dataRes <- tibble(GeneName=c("A", "C"), same_allele=c(TRUE, NA))
    expect_equal(dataRes, res$data[[1]])
})


# test function parse_hla_dataset()  --------------------------------------

context("Test parse_hla_dataset() function")

test_that("parse_hla_dataset() should return good result 01", {
    demo <- data.frame(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                AlleleName=c(1, 2, 1, 2), GeneName=c("A", "A", "A", "A"),
                AlleleGroup=c("02", "02", "03", NA),
                Protein=c("01", "01", "01", "02"),
                SynSubst=c("01", "02", "01", "03"),
                NonCoding=c("01", "01", NA, NA),
                Suffix=c(NA, NA, "L", NA))

    res <- HLAClustRView:::parse_hla_dataset(demo)

    expected <- tibble(SampleName=c("DEMO1", "DEMO1"),
                       AlleleName=c("1", "2"), GeneName=c("A", "A"),
                       AlleleGroup=c("02", "02"), Protein=c("01", "01"),
                       SynSubst=c("01", "02"), NonCoding=c("01", "01"),
                       Suffix=as.character(c(NA, NA)))

    expect_equal(res, expected)
})


test_that("parse_hla_dataset() should return good result 02", {
    demo <- data.frame(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                       AlleleName=c(1, 2, 1, 2), GeneName=c("A", "A", "A", "A"),
                       AlleleGroup=c("02", "02", "03", "03"),
                       Protein=c("01", "01", "01", "02"),
                       SynSubst=c("01", "02", "01", "03"),
                       NonCoding=c("01", "01", NA, NA),
                       Suffix=c(NA, NA, "L", NA))

    res <- HLAClustRView:::parse_hla_dataset(demo)

    expected <- tibble(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                       AlleleName=as.character(c(1, 2, 1, 2)), GeneName=c("A", "A", "A", "A"),
                       AlleleGroup=c("02", "02", "03", "03"), Protein=c("01", "01", "01", "02"),
                       SynSubst=c("01", "02", "01", "03"), NonCoding=c("01", "01", NA, NA),
                       Suffix=as.character(c(NA, NA, "L", NA)))

    expect_equal(res, expected)
})

test_that("parse_hla_dataset() should return good result 03", {
    demo <- data.frame(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                       AlleleName=c(1, NA, 1, 2), GeneName=c("A", "A", "A", "A"),
                       AlleleGroup=c("02", "02", "03", "03"), Protein=c("01", "01", "01", "02"),
                       SynSubst=c("01", "02", "01", "03"), NonCoding=c("01", "01", NA, NA),
                       Suffix=c(NA, NA, "L", NA))

    res <- HLAClustRView:::parse_hla_dataset(demo)

    expected <- tibble(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                       AlleleName=as.character(c(1, NA, 1, 2)), GeneName=c("A", "A", "A", "A"),
                       AlleleGroup=c("02", "02", "03", "03"), Protein=c("01", "01", "01", "02"),
                       SynSubst=c("01", "02", "01", "03"), NonCoding=c("01", "01", NA, NA),
                       Suffix=as.character(c(NA, NA, "L", NA)))

    expect_equal(res, expected)
})


# test function calculateHamming()  --------------------------------------

context("Test calculateHamming() function")


test_that("calculateHamming() should return an error when numerical for hla_data", {
    message <- "hla_data must be of class \"HLADataset\""
    expect_error(calculateHamming(33), message)
})

test_that("calculateHamming() should return an error when string for hla_data", {
    message <- "hla_data must be of class \"HLADataset\""
    expect_error(calculateHamming("Canada"), message)
})

test_that("calculateHamming() should return an error when data.frame for hla_data", {
    message <- "hla_data must be of class \"HLADataset\""
    demo <- data.frame(SampleName=c("s1", "s2"), GeneName=c("A", "A"))
    expect_error(calculateHamming(demo), message)
})

test_that("calculateHamming() should return an error when only one sample", {
    tempData <- demoHLADataset
    infoSamples <- tempData$data
    demo <- subset(infoSamples, infoSamples$SampleName %in% c("s1"))

    tempData$data <- demo

    message <- "hla_data must contain information for at least 2 samples"
    expect_error(calculateHamming(tempData), message)
})

test_that("calculateHamming() should return an error when not data", {
    tempData <- demoHLADataset

    tempData$data <- NULL

    message <- "A entry called \"data\" is missing from hla_data"
    expect_error(calculateHamming(tempData), message)
})


test_that("calculateHamming() should return good result 01", {
    tempData <- demoHLADataset
    infoSamples <- tempData$data
    demo <- subset(infoSamples, infoSamples$SampleName %in% c("s1", "s3"))

    tempData$data <- demo

    res <- calculateHamming(tempData)


    expected <- list()
    expected$dist <- matrix(c(0, 0, 5, 0), nrow=2, byrow = TRUE)
    colnames(expected$dist) <- c("s1", "s3")
    rownames(expected$dist) <- colnames(expected$dist)
    expected$metric <- "Hamming Distance"
    class(expected) <- "HLAMetric"

    expect_equal(res$dist, expected$dist)
})

test_that("calculateHamming() should return good result 02", {
    tempData <- demoHLADataset
    infoSamples <- tempData$data
    demo <- subset(infoSamples, infoSamples$SampleName %in% c("s3", "s4", "s5"))

    tempData$data <- demo

    res <- calculateHamming(tempData)

    alleleInfo <- list()
    alleleInfo[[1]] <-  tibble(GeneName=c("A", "B", "C"), same_allele=c(NA, NA, FALSE))
    alleleInfo[[2]] <-  tibble(GeneName=c("B", "C"), same_allele=c(NA, NA))
    alleleInfo[[3]] <-  tibble(GeneName=c("B", "C"), same_allele=c(TRUE, NA))


    expected <- matrix(c(0, 0, 0, 5, 0, 0, 4, 3, 0), nrow=3, byrow = TRUE)
    colnames(expected) <- c("s3", "s4", "s5")
    rownames(expected) <- colnames(expected)


    expected <- list()
    expected$dist <- matrix(c(0, 0, 0, 5, 0, 0, 4, 3, 0), nrow=3, byrow = TRUE)
    colnames(expected$dist) <- c("s3", "s4", "s5")
    rownames(expected$dist) <- colnames(expected$dist)

    expected$alleleInfo <- tibble(SampleName1=c("s3", "s3", "s4"), SampleName2=c("s4", "s5", "s5"),
           HammingDistance=c(5,4,3),
           data=alleleInfo)

    expected$metric <- "Hamming Distance"
    class(expected) <- "HLAMetric"


    expect_equal(res$dist, expected$dist)

})



# test function makeDistanceMatrix()  --------------------------------------

context("Test makeDistanceMatrix() function")

test_that("makeDistanceMatrix() should return good result 01", {

    inputTibble <- tibble(SampleName1=c("ERR188053", "ERR188053", "ERR188465"),
                       SampleName2=c("ERR188465", "ERR188040", "ERR188040"),
                       HammingDistance = c(17, 18, 19))

    result <- HLAClustRView:::makeDistanceMatrix(inputTibble)

    expected <- matrix(c(rep(0, 3), 17, 0, 0, 18, 19, 0), byrow = TRUE,
                            nrow = 3)
    rownames(expected) <- c("ERR188053", "ERR188465", "ERR188040")
    colnames(expected) <- rownames(expected)

    expect_equal(result, expected)
})

test_that("makeDistanceMatrix() should return good result 02", {

    inputTibble <- tibble(SampleName1=c("A", "A", "A", "A", "A", "A",
                                        "B", "B", "B", "B", "B",
                                        "C", "C", "C", "C",
                                        "D", "D", "D",
                                        "E", "E",
                                        "F"),
                          SampleName2=c("B", "C", "D", "E", "F", "G",
                                        "C", "D", "E", "F", "G",
                                        "D", "E", "F", "G",
                                        "E", "F", "G",
                                        "F", "G",
                                        "G"),
                          HammingDistance = c(1:21))

    result <- HLAClustRView:::makeDistanceMatrix(inputTibble)

    expected <- matrix(c(rep(0, 7),
                         1, rep(0, 6),
                         2, 7, rep(0, 5),
                         3, 8, 12, rep(0, 4),
                         4, 9, 13, 16, rep(0, 3),
                         5, 10, 14, 17, 19, rep(0, 2),
                         6, 11, 15, 18, 20, 21, 0), byrow = TRUE, nrow = 7)
    rownames(expected) <- c("A", "B", "C", "D", "E", "F", "G")
    colnames(expected) <- rownames(expected)

    expect_equal(result, expected)
})


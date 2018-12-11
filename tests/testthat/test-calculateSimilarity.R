### Unit tests for functions in calculateSimilarity.R file

library(HLAClustRView)
library(tibble)

data("example_sample_pair_data")


context("Teste hamming_distance_digit1() function")

d1 <- tibble::tribble(
        ~SampleName, ~AlleleName, ~AlleleGroup,
        "s1", 1, 1,
        "s1", 2, 3,
        "s2", 1, 1,
        "s2", 2, 5
        )
d2 <- tibble::tribble(
    ~SampleName, ~AlleleName, ~AlleleGroup,
    "s1", 1, 1,
    "s1", 2, 1,
    "s2", 1, 1,
    "s2", 2, 1
)

# test function hamming_distance_digit1 -----------------------------------
test_that("correct Hamming Value", {
    expect_equal(HLAClustRView:::hamming_distance_digit1(d1)$distance, 1)
    expect_equal(HLAClustRView:::hamming_distance_digit1(d2)$distance, 0)
    ## returns NA when distance is the same for both allels
    expect_true(is.na(HLAClustRView:::hamming_distance_digit1(d2)$same_allele))
})

test_that("correct output format", {
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


# test function parse_hla_data()  --------------------------------------

context("Test parse_hla_data() function")

test_that("parse_hla_data() should return good result 01", {
    demo <- data.frame(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                AlleleName=c(1, 2, 1, 2), GeneName=c("A", "A", "A", "A"),
                AlleleGroup=c("02", "02", "03", NA), Protein=c("01", "01", "01", "02"),
                SynSubst=c("01", "02", "01", "03"), NonCoding=c("01", "01", NA, NA),
                Suffix=c(NA, NA, "L", NA))

    res <- HLAClustRView:::parse_hla_data(demo)

    expected <- tibble(SampleName=c("DEMO1", "DEMO1"),
                       AlleleName=c("1", "2"), GeneName=c("A", "A"),
                       AlleleGroup=c("02", "02"), Protein=c("01", "01"),
                       SynSubst=c("01", "02"), NonCoding=c("01", "01"),
                       Suffix=as.character(c(NA, NA)))

    expect_equal(res, expected)
})


test_that("parse_hla_data() should return good result 02", {
    demo <- data.frame(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                       AlleleName=c(1, 2, 1, 2), GeneName=c("A", "A", "A", "A"),
                       AlleleGroup=c("02", "02", "03", "03"), Protein=c("01", "01", "01", "02"),
                       SynSubst=c("01", "02", "01", "03"), NonCoding=c("01", "01", NA, NA),
                       Suffix=c(NA, NA, "L", NA))

    res <- HLAClustRView:::parse_hla_data(demo)

    expected <- tibble(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                       AlleleName=as.character(c(1, 2, 1, 2)), GeneName=c("A", "A", "A", "A"),
                       AlleleGroup=c("02", "02", "03", "03"), Protein=c("01", "01", "01", "02"),
                       SynSubst=c("01", "02", "01", "03"), NonCoding=c("01", "01", NA, NA),
                       Suffix=as.character(c(NA, NA, "L", NA)))

    expect_equal(res, expected)
})

test_that("parse_hla_data() should return good result 03", {
    demo <- data.frame(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                       AlleleName=c(1, NA, 1, 2), GeneName=c("A", "A", "A", "A"),
                       AlleleGroup=c("02", "02", "03", "03"), Protein=c("01", "01", "01", "02"),
                       SynSubst=c("01", "02", "01", "03"), NonCoding=c("01", "01", NA, NA),
                       Suffix=c(NA, NA, "L", NA))

    res <- HLAClustRView:::parse_hla_data(demo)

    expected <- tibble(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
                       AlleleName=as.character(c(1, NA, 1, 2)), GeneName=c("A", "A", "A", "A"),
                       AlleleGroup=c("02", "02", "03", "03"), Protein=c("01", "01", "01", "02"),
                       SynSubst=c("01", "02", "01", "03"), NonCoding=c("01", "01", NA, NA),
                       Suffix=as.character(c(NA, NA, "L", NA)))

    expect_equal(res, expected)
})


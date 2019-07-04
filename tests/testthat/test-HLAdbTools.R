### Unit tests for functions in HLAdbTools.R file

library(HLAClustRView)
library(dplyr)

directory <- system.file("extdata", package = "HLAClustRView")
HLA_INFO <- parseHLADbAlignment(hlaDbPath = directory, seqType = "prot")

### Tests getSeqCMP() function

context("Test for getSeqCMP() function")

test_that("getSeqCMP() must return good result", {

    #sample1 <- "DRA*01:01:01:03"
    #sample2 <- "DRA*01:01:02"
    inputTibble <- tibble(SampleName=c("ERR188021", "ERR188021"),
                          AlleleName=c("1","2"),
                          GeneName=c("DRA", "DRA"),
                          AlleleGroup=c("01", "01"),
                          Protein=c("01", "01"),
                          SynSubst=c("01", "02"),
                          NonCoding=c("03", NA),
                          Suffix=c(NA, NA)
                          )
    sample1 <- inputTibble[1,]
    sample2 <- inputTibble[2,]

    regions <- data.frame(start=c(160, 200, 240), end=c(180, 220, 260))

    result <- HLAClustRView:::getSeqCMP(HLAInfo = HLA_INFO, regionExt = regions,
                                        typeS1 = sample1, typeS2 = sample2)

    expected <- list(refSeq="VYDCRVEHWGLDEPLLKHWEFTVGLVGIIIGTIFIIKGVRKS",
                        seqS1="------------------------------------------",
                        seqS2="------------------------------------------")


    expect_equal(result, expected)
})


### Tests getSubSeq() function

context("Test for getSubSeq() function")

test_that("getSubSeq() must return good result", {

    sequence <- "    MGSGWVPWVVALLVNLTRLDSSMTQGTDSPEDFVIQAKADCYFTNGTEKVQFVVRFIFNLEEYVRFDSDVGM"

    position <- -30

    regions <- data.frame(start=c(15, 25), end=c(20, 30))


    result <- HLAClustRView:::getSubSeq(seq = sequence, posInit = position,
                                            regionExt = regions)

    expected <-  "CYFTNGQFVVRF"


    expect_equal(result, expected)
})

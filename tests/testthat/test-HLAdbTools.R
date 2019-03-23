### Unit tests for functions in HLAdbTools.R file

library(HLAClustRView)


directory <- system.file("extdata", package = "HLAClustRView")
HLA_INFO <- parseHLADbAlignment(hlaDbPath = directory, seqType = "prot")

### Tests getSeqCMP() function

context("Test for getSeqCMP() function")

test_that("getSeqCMP() must return good result", {

    sample1 <- "DRA*01:01:01:03"
    sample2 <- "DRA*01:01:02"

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

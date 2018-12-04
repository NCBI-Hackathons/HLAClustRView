### Unit tests for functions in prepareHLAdb.R file

library(HLAClustRView)

directory <- system.file("extdata", package = "HLAClustRView")

### Tests parseHLADbAlignment() results

context("parseHLADbAlignment() results")

test_that("parseHLADbAlignment() must return an error when directory is a number", {
    message <- "The hlaDbPath parameter must by a character string"
    expect_error(parseHLADbAlignment(hlaDbPath = 33, seqType = "nuc"), message)
})

test_that("parseHLADbAlignment() must return an error when no a valid sequence type", {
    message <- "Not validate sequence type parameter for parseHLADbAlignment: toto"
    expect_error(parseHLADbAlignment(hlaDbPath = directory, seqType = "toto"), message)
})


test_that("parseHLADbAlignment() must return an error when directory does not exist", {
    path_to_test <- "./3_21_33eddger_39kdf_23sk2-2221fvd_qqqqqqqqqq_2"
    if (!dir.exists(path_to_test)) {
        message <- "The hlaDbPath parameter must by a valid directory"
        expect_error(parseHLADbAlignment(hlaDbPath = path_to_test, seqType = "nuc"), message)
    }
})

test_that("parseHLADbAlignment() must return an error when no file in directory", {
    message <- "There must be at least one alignment file in the hlaDbPath directory"
    expect_error(parseHLADbAlignment(hlaDbPath = directory, seqType = "nuc"), message)
})

#
# test_that("parseHLADbAlignment() must return an error when directory does not exist", {
#      result <- parseHLADbAlignment(hlaDbPath = directory, seqType = "prot")
#
#      expected <- list()
#      expected[["refSeq"]] <- list()
#      expected[["refSeq"]][['DOB']] <- paste0("NLTRLDSSMTQGTDSPEDFVIQAKADCYFTNGTEKVQFVVRFIFN",
#                 "LEEYVRFDSDVGMFVALTKLGQPDAEQWNSRLDLLERFTVGRKVQPEVTVYPERTPLLHQHNLLHCSVTGFY",
#                 "PGDIKIKWFLNGQEERAGVMSTGPIRNGDWTFQTVVMLEMTPELGHVRAQSEYSWRKMLSGIAAFLLGLIFL",
#                 "LVGIVIQLRAQKGYVRTQMSGNEVSRAVLLPQS")
#      expected[["refSeq"]][['DRA']] <- paste0("IAVLMSAQESWAIKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFH",
#                 "VDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLVTVLTNSPVELREPNVLICFIDKFTPPVVNVTWLR",
#                 "NGKPVTTGVSETVFLPREDHLFRKFHYLPFLPSTEDVYDCRVEHWGLTTENVVCALGLTVGLVGIIIGTIFI",
#                 "IKGVRKSNAAERRGPL")
#      expected[["posInit"]] <- list()
#      expected[["posInit"]][['DOB']] <- 171
#      expected[["posInit"]][['DRA']] <- 171
#
#      expected[["HLAAlignment"]] <- list()
#
#
#      ## TODO
# })


context("extractTyping() results")

test_that("extractTyping() must return good result 01", {
    sequence <- " DOB*01:01:01:02       MGSGWV PWVVALLVNL TRLDSSMTQG"
    result <- HLAClustRView:::extractTyping(seq=sequence, endPos=20)

    expected <- "DOB*01:01:01:02"
    expect_equal(result, expected)
})

test_that("extractTyping() must return good result 02", {
    sequence <- "            MGSGWV PWVVALLVNL TRLDSSMTQG"
    result <- HLAClustRView:::extractTyping(seq=sequence, endPos=20)

    expected <- ""
    expect_equal(result, expected)
})

context("extractRef() results")

test_that("extractRef() must return good result 01", {
    sequence <- " DOB*01:01:01:01       MGSGWV PWVVALLVNL TRLDSSMTQG TDSPEDFVIQ AKADCYFTNG"
    result <- HLAClustRView:::extractRef(seq=sequence, startPos=20)

    expected <- list()
    expected[["refSeq"]] <- "    MGSGWVPWVVALLVNLTRLDSSMTQGTDSPEDFVIQAKADCYFTNG"
    expected[["seqDiff"]] <- "    ----------------------------------------------"

    expect_equal(result, expected)
})

test_that("extractRef() must return good result 02", {
    sequence <- " A*01:01:01:01     ATG GCC GTC ATG GCG CCC CGA ACC CTC CTC CTG CTA CTC TCG GGG  "
    result <- HLAClustRView:::extractRef(seq=sequence, startPos=21)

    expected <- list()
    expected[["refSeq"]] <- "TGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGG"
    expected[["seqDiff"]] <- "--------------------------------------------"

    expect_equal(result, expected)
})

context("extractSeq() results")

test_that("extractSeq() must return good result 01", {
    sequence <- " F*01:06           --- --- --- --- --- --- --- --- --- -|-- "
    result <- HLAClustRView:::extractSeq(seq=sequence, startPos=20)

    expected <- "------------------------------"
    expect_equal(result, expected)
})

test_that("extractSeq() must return good result 02", {
    sequence <- " G*01:17           ----H----- --******** ********** ****"
    result <- HLAClustRView:::extractSeq(seq=sequence, startPos=20)

    expected <- "----H-------**********************"
    expect_equal(result, expected)
})



context("getTypingPos() results")


test_that("getTypingPos() must return good result 01", {
    seqData <- data.table(GeneName=rep("DOB", 3), AlleleGroup=rep("01", 3),
                    Protein=rep("01", 3), SynSubst=rep("01", 3),
                    Noncoding=c("01", "02", "03"), Suffix=rep(NA, 3))

    typing <- matrix(data=c("DOB", "01", "01", "01", "03", NA), nrow=1)

    result <- HLAClustRView:::getTypingPos(seqProcess=seqData, curTyping=typing)

    expected <- 3
    expect_equal(result, expected)
})


test_that("extractSeq() must return good result 02", {
    seqData <- data.table(GeneName=c(rep("DOB", 3), "DOA"), AlleleGroup=rep("01", 4),
                          Protein=rep("01", 4), SynSubst=rep("01", 4),
                          Noncoding=c("01", "02", "03", "02"), Suffix=rep(NA, 4))

    typing <- matrix(data=c("DOA", "01", "01", "01", "02", NA), nrow=1)


    result <- HLAClustRView:::getTypingPos(seqProcess=seqData, curTyping=typing)

    expected <- 4
    expect_equal(result, expected)
})

test_that("extractSeq() must return NA when typing not present", {
    seqData <- data.table(GeneName=c(rep("DOB", 3), "DOA"), AlleleGroup=rep("01", 4),
                          Protein=rep("01", 4), SynSubst=rep("01", 4),
                          Noncoding=c("01", "02", "03", "02"), Suffix=rep(NA, 4))

    typing <- matrix(data=c("DOA", "01", "01", "01", "03", NA), nrow=1)


    result <- HLAClustRView:::getTypingPos(seqProcess=seqData, curTyping=typing)

    expected <- integer()
    expect_equal(result, expected)
})

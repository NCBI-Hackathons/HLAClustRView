library(HLAClustRView)

d1 <- tibble::tribble(
        ~sample, ~Allele, ~Digit1,
        "s1", 1, 1,
        "s1", 2, 3,
        "s2", 1, 1,
        "s2", 2, 5,
        )
d2 <- tibble::tribble(
    ~sample, ~Allele, ~Digit1,
    "s1", 1, 1,
    "s1", 2, 1,
    "s2", 1, 1,
    "s2", 2, 1,
)
# test function hamming_distance_digit1 -----------------------------------
test_that("correct Hamming Value", {
    expect_equal(hamming_distance_digit1(d1)$distance, 1)
    expect_equal(hamming_distance_digit1(d2)$distance, 0)
})

test_that("correct output format", {
    expect_equal(nrow(hamming_distance_digit1(d2)), 1)
    expect_equal(ncol(hamming_distance_digit1(d2)), 2)

})

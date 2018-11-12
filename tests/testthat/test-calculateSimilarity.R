d1 <- tibble::tribble(
        ~sample, ~Allele, ~Digit1,
        "s1", 1, 1,
        "s1", 2, 3,
        "s2", 1, 1,
        "s2", 2, 5
        )
d2 <- tibble::tribble(
    ~sample, ~Allele, ~Digit1,
    "s1", 1, 1,
    "s1", 2, 1,
    "s2", 1, 1,
    "s2", 2, 1
)
# test function hamming_distance_digit1 -----------------------------------
test_that("correct Hamming Value", {
    expect_equal(hamming_distance_digit1(d1)$distance, 1)
    expect_equal(hamming_distance_digit1(d2)$distance, 0)
    ## returns NA when distance is the same for both allels
    expect_true(is.na(hamming_distance_digit1(d2)$same_allele))
})

test_that("correct output format", {
    expect_equal(nrow(hamming_distance_digit1(d2)), 1)
    expect_equal(ncol(hamming_distance_digit1(d2)), 2)

})

# test function sample_pair_distance --------------------------------------
data("example_sample_pair_data")
res <- sample_pair_distance(example_sample_pair_data)
test_that("correct distance (agregated sum)", {
    expect_equal(res$HammingDistance, 3)
    # keep only genes that are in both samples
    expect_true(all(c("A", "C") %in% tidyr::unnest(res)$HLAgene))

})

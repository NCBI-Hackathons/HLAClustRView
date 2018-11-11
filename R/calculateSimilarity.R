#' Compute the hamming distance between two alleles
#'
#' This function computes the hamming distance for two samples and two alleles,
#'  the hamming distance is given by
#' \eqn{min(s1A1 != s2A1 + s1A2 != s2A2, s1A1 != s2A2 + s1A2 != s2A1)}
#' @param allele tibble with 3 columns: *sample*, *Allele*, *Digit1*
#' (see example)
#' @author Santiago Medina, Nissim Ranade
#' @return One row tibble with minimum hamming distance &
#' logical value column indicating whether same allele was used If the distance is the
#' same for both then NA is returned
#'
#' @examples
#' d <-
#' tibble::tribble(
#'    ~sample, ~Allele, ~Digit1,
#'    "s1", 1, 1,
#'    "s1", 2, 3,
#'    "s2", 1, 1,
#'    "s2", 2, 5,
#')
#' hamming_distance_digit1(d)
#' @importFrom dplyr %>% mutate group_by summarise filter slice
#' @importFrom purrr reduce
#' @export
hamming_distance_digit1 <- function(allele) {
    stopifnot(nrow(allele) == 4)
    hdist <-
        allele %>%
        split(.$sample) %>%
        reduce(.f = tidyr::crossing) %>% # biparty graph to compute simularities
        mutate(
            is_similar = Digit1 != Digit11,
            same_allele = Allele == Allele1
            ) %>%
        group_by(same_allele) %>%
        summarise(distance = sum(is_similar)) %>%
        filter(distance == min(distance))

    ## check if the minimum distance is the same
    if (nrow(hdist) == 2) {
        mutate(hdist, same_allele = NA) %>%
            slice(1:1)
    }
    else hdist
}

#' Sample Pair Distance
#'
#' Computes the Hamming Distance for a pair of samples and then aggregates the
#' sum over the genes keeping information of which allele was used for the distance
#' @param sample_pair_data a tible with columns HLAgene (genes), Allele (alleles),
#' Digit1 (digit 1) and sample (see example data)
#'
#' @author Santiago Medina, Nissim Ranade
#' @return  a tibble with one row and column HammingDistance & column data
#' corresponding to the same_allele information
#' @export
#' @importFrom tidyr nest unnest spread
#' @importFrom dplyr group_by %>% summarise filter pull select mutate
#' @importFrom utils data
#' @importFrom purrr map
#' @examples
#' data("example_sample_pair_data") # example data
#' sample_pair_distance(example_sample_pair_data)
sample_pair_distance <- function(sample_pair_data) {
    # make sure only two samples are given
    stopifnot(length(unique(sample_pair_data$sample)) == 2)

    # get complete cases, a gene should apear in both samples
    # to get the distance
    complete_genes <-
        sample_pair_data %>%
        group_by(HLAgene) %>%
        summarise(total = n()) %>%
        filter(total == 4) %>%
        pull(HLAgene)

    distances_per_gene <-
        sample_pair_data %>%
        filter(HLAgene %in% complete_genes) %>% # keep complete cases
        group_by(HLAgene) %>%
        nest() %>%
        mutate(x = map(data, hamming_distance_digit1)) %>%
        select(-data) %>%
        unnest(x)

    distances_per_gene %>%
        mutate(HammingDistance = sum(distance)) %>%
        select(-distance) %>%
        group_by(HammingDistance) %>%
        nest()
}


#' calculate hamming distance between samples
#'
#' @param hla_data data frame with allele data
#'
#' @importFrom dplyr filter mutate select group_by inner_join as_tibble rename %>%
#' @importFrom tidyr unnest nest
#' @importFrom utils combn
#' @importFrom purrr map_lgl
#' @return
#' @export
#'
#' @examples
calculateSimilarity <- function(hla_data) {
    ssample_pair_distance <- possibly(sample_pair_distance, otherwise = FALSE)

    map_data_to_distance <- function(s1, s2) {
        # get's the data for s1 and s2 and applys the distance function
        filter(hla_data, sample %in% c(s1, s2)) %>%
            ssample_pair_distance()

    }

    sample_pairs <-
        unique(hla_data$sample) %>%
        combn(m = 2) %>%
        t() %>%
        as_tibble() %>%
        rename(sample1 = V1, sample2 = V2)

    distances <-
        sample_pairs %>%
        mutate(x = map2(sample1, sample2, map_data_to_distance)) %>%
        mutate(is_error = map_lgl(x, is.logical)) %>%
        filter(!is_error) %>%
        select(-is_error) %>%
        unnest(x)

    distances %>%
        select(-HammingDistance) %>%
        unnest(data) %>%
        group_by(sample1, sample2) %>%
        nest() %>%
        inner_join(
            y = select(distances, sample1, sample2, HammingDistance),
            by = c("sample1", "sample2")
            )

}


#' Compute the hamming distance between two alleles
#'
#' This function computes the hamming distance for two samples and two alleles,
#'  the hamming distance is given by
#' \eqn{min(s1A1 != s2A1 + s1A2 != s2A2, s1A1 != s2A2 + s1A2 != s2A1)}
#' @param allele tibble with 3 columns: *sample*, *Allele*, *Digit1*
#' (see example)
#' @return One row tibble with minimum hamming distance &
#' logical value column indicating whether same allele was used
#' @export
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
hamming_distance_digit1 <- function(allele) {
    stopifnot(nrow(allele) == 4)
    allele %>%
        split(.$sample) %>%
        purrr::reduce(.f = tidyr::crossing) %>% # biparty graph to compute simularities
        dplyr::mutate(
            is_similar = Digit1 != Digit11,
            same_allele = Allele == Allele1
            ) %>%
        dplyr::group_by(same_allele) %>%
        dplyr::summarise(distance = sum(is_similar)) %>%
        dplyr::filter(distance == min(distance)) %>%
        dplyr::slice(1:1) # make sure you return one row in case of repeats
}

#' Compute the hamming distance between two alleles
#'
#' This function computes the hamming distance for two samples and two alleles,
#'  the hamming distance is given by
#' \eqn{min(s1A1 != s2A1 + s1A2 != s2A2, s1A1 != s2A2 + s1A2 != s2A1)}
#' @param allele tibble with 3 columns: *sample*, *Allele*, *Digit1*
#' (see example)
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
#' @importFrom dplyr %>% mutate group_by summarise filter
#' @export
hamming_distance_digit1 <- function(allele) {
    stopifnot(nrow(allele) == 4)
    hdist <-
        allele %>%
        split(.$sample) %>%
        purrr::reduce(.f = tidyr::crossing) %>% # biparty graph to compute simularities
        mutate(
            is_similar = Digit1 != Digit11,
            same_allele = Allele == Allele1
            ) %>%
        group_by(same_allele) %>%
        summarise(distance = sum(is_similar)) %>%
        dplyr::filter(distance == min(distance))

    ## check if the minimum distance is the same
    if (nrow(hdist) == 2) {
        mutate(hdist, same_allele = NA) %>%
            dplyr::slice(1:1)
    }
    else hdist
}

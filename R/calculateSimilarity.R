#' Compute the hamming distance between two alleles
#'
#' This function computes the hamming distance for two samples and two alleles,
#'  the hamming distance is given by
#' \eqn{min(s1A1 != s2A1 + s1A2 != s2A2, s1A1 != s2A2 + s1A2 != s2A1)}
#' @param allele tibble with 3 columns: SampleName, AlleleName, AlleleGroup,
#' (see example)
#' @author Santiago Medina, Nissim Ranade
#' @return One row tibble with minimum hamming distance &
#' logical value column indicating whether same allele was used If the distance is the
#' same for both then NA is returned
#'
#' @examples
#' d <-
#' tibble::tribble(
#'    ~SampleName, ~AlleleName, ~AlleleGroup,
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
        split(.$SampleName) %>%
        reduce(.f = tidyr::crossing) %>% # biparty graph to compute simularities
        mutate(
            is_similar = AlleleGroup != AlleleGroup1,
            same_allele = AlleleName == AlleleName1
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
#' @param sample_pair_data a tible with columns GeneName (genes), Allele (alleles),
#' Digit1 (digit 1) and sample (see example data)
#'
#' @return  a tibble with one row and column HammingDistance & column data
#' corresponding to the same_allele information
#'
#' @examples
#'
#' ## Load example dataset
#' data("example_sample_pair_data")
#'
#' ## Computes the Hamming distance
#' #sample_pair_distance(example_sample_pair_data)
#'
#' @author Santiago Medina, Nissim Ranade
#'
#' @importFrom tidyr nest unnest spread
#' @importFrom dplyr group_by %>% summarise filter pull select mutate
#' @importFrom utils data
#' @importFrom purrr map
#' @export
sample_pair_distance <- function(sample_pair_data) {
    # make sure only two samples are given
    stopifnot(length(unique(sample_pair_data$SampleName)) == 2)

    # get complete cases, a gene should apear in both samples
    # to get the distance
    complete_genes <-
        sample_pair_data %>%
        group_by(GeneName) %>%
        summarise(total = n()) %>%
        filter(total == 4) %>%
        pull(GeneName)

    distances_per_gene <-
        sample_pair_data %>%
        filter(GeneName %in% complete_genes) %>% # keep complete cases
        group_by(GeneName) %>%
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


#' @title calculate hamming distance between samples
#'
#' @description Takes the parsed hla database and computes the hamming distance
#' for each pair of samples
#'
#' @param hla_data data frame with allele data
#'
#' @return TODO
#'
#' @examples
#'
#' ## Load example dataset
#' data(example_calculateSimilarity)
#'
#' ## Calculate similarity metrics
#' #calculateSimilarity(example_calculateSimilarity)
#'
#' @importFrom dplyr filter mutate select group_by inner_join as_tibble rename %>%
#' @importFrom tidyr unnest nest
#' @importFrom utils combn
#' @importFrom purrr map_lgl possibly map2
#' @export
calculateSimilarity <- function(hla_data) {
    hla_data <- select(hla_data, SampleName, GeneName, AlleleName, AlleleGroup)
    ssample_pair_distance <- possibly(sample_pair_distance, otherwise = FALSE)

    map_data_to_distance <- function(s1, s2) {
        # get's the data for s1 and s2 and applys the distance function
        filter(hla_data, SampleName %in% c(s1, s2)) %>%
            ssample_pair_distance()

    }

    sample_pairs <-
        unique(hla_data$SampleName) %>%
        combn(m = 2) %>%
        t() %>%
        as_tibble() %>%
        rename(SampleName1 = V1, SampleName2 = V2)

    distances <-
        sample_pairs %>%
        mutate(x = map2(SampleName1, SampleName2, map_data_to_distance)) %>%
        mutate(is_error = map_lgl(x, is.logical)) %>%
        filter(!is_error) %>%
        select(-is_error) %>%
        unnest(x)

    distances %>%
        select(-HammingDistance) %>%
        unnest(data) %>%
        group_by(SampleName1, SampleName2) %>%
        nest() %>%
        inner_join(
            y = select(distances, SampleName1, SampleName2, HammingDistance),
            by = c("SampleName1", "SampleName2")
            ) %>%
        rename(AlleleName_info = data)

}

#' Parse hladb
#'
#' Converts object to tibble, removes samples with missing values
#' @param hladb hla database the output of function \code{parseHLADb}
#'
#' @return
#' @export
#' @importFrom purrr map reduce set_names
#' @importFrom dplyr mutate_all filter pull bind_cols as_tibble %>%
#' @examples
parse_hla_data <- function(hladb) {

    hla_data <-
        hladb %>%
        map(as_tibble) %>%
        reduce(bind_cols) %>%
        set_names(names(hladb)) %>%
        mutate_all(.funs = as.character)

    missing_samples <-
        hla_data %>%
        filter(is.na(AlleleGroup)) %>%
        pull(SampleName)

    hla_data %>%
        filter(!SampleName %in% missing_samples)


}

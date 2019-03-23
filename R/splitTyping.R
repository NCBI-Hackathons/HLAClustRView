#' @title splitTyping
#'
#' @description Converts a well-formed HLA allele name to a matrix of fields
#'
#' @param aChar a \code{character} array containing one or more well-formed
#' HLA alleles,
#' see http://hla.alleles.org/nomenclature/naming.html
#'
#' @return a n x 6 \code{matrix} of fields
#'
#' @examples
#'
#' ## Split one HLA typing into an array of entries
#' output <- HLAClustRView:::splitTyping("HLA-DRB1*13:01:01:02")
#'
#' print(output)
#'
#' ## Split one HLA typing into an array of entries
#' output <- HLAClustRView:::splitTyping(c("G*01:04:01:02", "DQA1*01:01:01:03"))
#'
#' print(output)
#'
#' @author Adewunmi Adelaja
#'
#' @importFrom stringr str_extract
#' @keywords internal
splitTyping <- function (aChar) {

    prefix <- str_extract(aChar, "(^.*[:digit:])|(^[:alpha:]+$)")
    splitted <- strsplit(prefix,"[.*:_]")
    splitted <- t(vapply(splitted,
                            FUN=function(x) {c(x, rep(NA, 5 - length(x)))},
                            FUN.VALUE=character(5)))
    suffix <- str_extract(aChar, "(?<=([:digit:]|[.*:_]))[:alpha:]+$")
    output <- cbind(splitted, as.array(suffix))

    return(output)
}

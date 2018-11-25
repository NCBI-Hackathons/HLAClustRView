#' @title
#'splitTyping
#' @description
#'Converts a well-formed HLA allele name to a matrix of fields
#' @param aChar
#'Character array containing well-formed HLA allele
#'see http://hla.alleles.org/nomenclature/naming.html
#' @return
#'n x 6 matrix of fields
#' @examples
#'output <- splitTyping("HLA-DRB1*13:01:01:02")
#'print(output)
#'#[,1] [,2] [,3] [,4] [,5] [,6]
#'#"HLA-DRB1" "13" "01" "01" "02" NA
#'
#' @author Adewunmi Adelaja
#'
#' @importFrom stringr str_extract
#' @export
splitTyping <- function (aChar)
{
    options('stringAsFactor'=FALSE)
    prefix <-str_extract(aChar, "(^.*[:digit:])|(^[:alpha:]+$)")
    splitted <- strsplit(prefix,"[.*:_]")
    splitted<-t(sapply(splitted, function(x) c(x,rep(NA, 5-length(x)))))
    suffix <-str_extract(aChar, "(?<=([:digit:]|[.*:_]))[:alpha:]+$")
    output <-cbind(splitted, as.array(suffix));

    return(output)

}

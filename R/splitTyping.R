#' @title TODO
#'
#' @description TODO
#'
#' @param aChar TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
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
    # output <-matrix(unlist(cbind(splitted, as.array(suffix))), ncol = 6, byrow = TRUE);
    return(output)

}

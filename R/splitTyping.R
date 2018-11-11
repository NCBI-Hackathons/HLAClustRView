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
    aChar<-"A*31:01:02:01"
    prefix <-str_extract(aChar, '(^.*[:digit:])|(^[:alpha:]+$)')
    splitted <- strsplit(prefix,'[.*:_]')
    suffix <-str_extract(aChar, '(?<=([:digit:]|[.*:_]))[:alpha:]+$')
    output <- append(unlist(splitted), suffix)
    return(output)

}

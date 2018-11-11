splitTyping <- function (aChar)
{
    library (stringr)
    options('stringAsFactor'=FALSE)
    prefix <-str_extract(aChar, "(^.*[:digit:])|(^[:alpha:]+$)")
    splitted <- strsplit(prefix,"[.*:_]")
    pad <-5-lengths(splitted)
    splitted <- append(splitted,rep(NA, pad) )
    suffix <-str_extract(aChar, "(?<=([:digit:]|[.*:_]))[:alpha:]+$")
    output <- c(unlist(splitted), suffix)

    return(output)

}

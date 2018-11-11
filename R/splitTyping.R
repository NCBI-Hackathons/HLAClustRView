splitTyping <- function (aChar)
{
    library (stringr)
    options('stringAsFactor'=FALSE)
    prefix <-str_extract(aChar, "(^.*[:digit:])|(^[:alpha:]+$)")
    splitted <- strsplit(prefix,"[.*:_]")

    splitted<-t(sapply(splitted, function(x) c(x,rep(NA, 5-length(x)))))

    suffix <-str_extract(aChar, "(?<=([:digit:]|[.*:_]))[:alpha:]+$")
    output <-cbind(splitted, as.array(suffix))
    return(output)

}

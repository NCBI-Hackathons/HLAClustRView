#' @title Get the sequence from 2 allele types from the same HLA gene.
#'
#' @description Extract from the \code{HLAdb} object, two sequences and
#' the reference for a specific region related to 2 specific allele types.
#'
#' @param HLAInfo An object of class \code{HLAdb} contain information from
#' HLA database.
#'
#' @param regionExt TODO
#'
#' @param typeS1 A line of the tibble data from the object data return readHLADataset
#'
#' @param typeS2 A line of the tibble data from the object data return readHLADataset
#' HLA allele.
#'
#' @return An object of class TODO
#'
#' @details TODO
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom data.table data.table rbindlist
#' @export
getSeqCMP <- function(HLAInfo, regionExt, typeS1, typeS2) {

    if (!("HLAdb" %in% class(HLAInfo))) {
        stop("HLAInfo must be of class \"HLAdb\"")
    }

    splitS1 <- unlist(typeS1[,c("GeneName", "AlleleGroup",
                                "Protein", "SynSubst",
                                "NonCoding", "Suffix")])
    splitS2 <- unlist(typeS2[,c("GeneName", "AlleleGroup",
                                "Protein", "SynSubst",
                                "NonCoding", "Suffix")]) #splitTyping(typeS2)

    posS1 <- getIncompleteTypingPos(HLAInfo$HLAAlignment, splitS1)
    posS1 <- reduceTypingPos(HLAInfo$HLAAlignment, posS1)
    posS2 <- getIncompleteTypingPos(HLAInfo$HLAAlignment, splitS2)
    posS2 <- reduceTypingPos(HLAInfo$HLAAlignment, posS2)

    if (splitS1[1] != splitS2[1]) {
        stop("Call get seq with type from 2 genes")
    }
    if (is.na(posS1) || is.na(posS2)) {
        stop(paste0("Typing without specific sequence ", typeS1, " ", typeS2))
    }

    refSeq <- HLAInfo$refSeq[[splitS1[1]]]
    posInit <- HLAInfo$posInit[[splitS1[1]]]

    seqCMP <- list(refSeq="", seqS1="", seqS2="")
    seqCMP$refSeq <- getSubSeq(refSeq, posInit, regionExt)
    if (!(is.na(posS1))) {
        seqCMP$seqS1 <- getSubSeq(HLAInfo$HLAAlignment[posS1]$SeqDiff,
                                    posInit, regionExt)
    }
    if (!(is.na(posS2))) {
        seqCMP$seqS2 <- getSubSeq(HLAInfo$HLAAlignment[posS2]$SeqDiff,
                                    posInit, regionExt)
    }

    return(seqCMP)
}

#' @title Get the subsequence
#'
#' @description Get from the object HLAdb two sequences and the reference
#' for a region
#'
#' @param seq TODO
#'
#' @param posInit TODO
#'
#' @param regionExt TODO
#'
#' @return An object of class TODO
#'
#' @details TODO
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @export
getSubSeq <- function(seq, posInit, regionExt) {

    subSeq <- ""

    for(i in seq_len(nrow(regionExt))) {
        subSeq <- paste0(subSeq, substr(seq, regionExt$start[i] - posInit,
                                            regionExt$end[i] - posInit))
    }

    return(subSeq)
}

#' @title load matrix of substitution
#'
#' @description Load matrix BLOSOM or PAM
#' from ftp://ftp.ncbi.nih.gov/blast/matrices/
#'
#' @param fileName TODO
#'
#' @return a matrix
#'
#' @details TODO
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom utils read.table
#' @export
parseSubMatrix <- function(fileName) {

    subMatrix <- read.table(fileName, header = TRUE)
    return(subMatrix)
}



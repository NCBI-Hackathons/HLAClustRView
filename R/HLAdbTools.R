#' @title Get the subsequences from 2 allele types from the same HLA gene.
#'
#' @description Extract from the \code{HLAdb} object, two sequences and
#' the reference for a group of specific regions. The final sequences are
#' the concatenation of the sequences of each region.
#'
#' @param HLAInfo An object of class \code{HLAdb} containing information from
#' HLA database.
#'
#' @param regionExt A \code{data.frame} containing 2 columns called
#' \code{Start} and \code{end}. Those columns must contain \code{integer}
#' entries corresponding to the regions to retain to create the subsequence.
#' When more than one entry is present, the final subsequence will be the
#' concatenation of all retained regions.
#'
#' @param typeS1 A \code{character} array containing one well-formed
#' HLA allele.
#'
#' @param typeS2 A \code{character} array containing one well-formed
#' HLA allele.
#'
#' @return A \code{list} containing:
#' \itemize{
#' \item \code{refSeq} a \code{character} string containing the reference.
#' \item \code{seqS1} a \code{character} string containing the sequence of the
#' first HLA allele.
#' \item \code{seqS2} a \code{character} string containing the sequence of the
#' second HLA allele.
#' }
#'
#' @details TODO
#'
#' @examples
#'
#' ## Get HLA protein database
#' data(hladb_protein_3.35.0)
#'
#' ## Two samples from same HLA gene
#' sample1 <- "DRA*01:01:01:03"
#' sample2 <- "DRA*01:01:02"
#'
#' ## Fix regions that will be extracted
#' regions <- data.frame(start=c(160, 200, 240), end=c(180, 220, 260))
#'
#' ## Extract subsequences
#' getSeqCMP(HLAInfo = hladb_protein_3.35.0, regionExt = regions,
#'     typeS1 = sample1, typeS2 = sample2)
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom data.table data.table rbindlist
#' @export
getSeqCMP <- function(HLAInfo, regionExt, typeS1, typeS2) {

    if (!("HLAdb" %in% class(HLAInfo))) {
        stop("HLAInfo must be of class \"HLAdb\"")
    }

    if (!is.data.frame(regionExt)) {
        stop("regionExt must a \"data.frame\"")
    }

    splitS1 <- splitTyping(typeS1)
    splitS2 <- splitTyping(typeS2)

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
#' @description Extract, from a HLAdb object, two sequences and the reference
#' for a region specified by user
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
#' @keywords internal
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



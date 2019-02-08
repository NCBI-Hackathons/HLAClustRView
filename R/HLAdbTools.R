#' @title Get the reference sequence and the sequences of two allele typing
#' for a specific region.
#'
#' @description From an \code{HLAdb} object that contains HLA aligment
#' information for all typing, extract the reference sequence of a
#' specific region of a HLA
#' gene and the sequences associated to two allele typing for the same region.
#' The two typing must be for the same HLA gene. This enable an aligment
#' comparaison of the two allele typing.
#'
#' @param HLAInfo an object of class \code{HLAdb} containing the global HLA
#' aligment for all typing.
#'
#' @param regionExt a \code{data.frame} containing a \code{start} and an
#' \code{end} column. The \code{start} and \code{end} positions are both
#' included in the extracted region.
#'
#' @param typeS1 a \code{character} string that represent the typing of the
#' first sample. The \code{typeS1} and \code{typeS2} typing must be for the
#' same HLA gene.
#'
#' @param typeS2 a \code{character} string that represent the typing of the
#' second sample. The \code{typeS1} and \code{typeS2} typing must be for the
#' same HLA gene.
#'
#' @return A \code{list} containing the following elements:
#' \itemize{
#' \item \code{data} a \code{tibble} object containing the HLA typing
#' information for all samples. The columns are:
#' \itemize{
#' \item \code{refSeq} a \code{character} string that represent the
#' reference sequence for the specified region of the HLA gene.
#' \item \code{seqS1} a \code{character} string that represent the
#' sequence of the first sample for the specified region of the allele typing.
#' \item \code{seqS2} a \code{character} string that represent the
#' sequence of the second sample for the specified region of the allele typing.
#' }
#' }
#'
#' @examples
#'
#' ## Get path where some HLA database files are stored
#' directory <- system.file("extdata", package = "HLAClustRView")
#'
#' ## Parse HLA database files of protein type
#' HLAInfo <- parseHLADbAlignment(hlaDbPath=directory, seqType="prot")
#'
#' ## Select two allele typing for the same HLA gene
#' sample1 <- "DRA*01:01:01:03"
#' sample2 <- "DRA*01:02:03"
#'
#' ## Select a region to extract
#' regions <- data.frame(start=c(160, 200), end=c(180, 220))
#'
#' ## Extract the aligment information for the specified region included the
#' ## reference sequence and the sequences for both allele typing
#' getSeqCMP(HLAInfo = HLAInfo, regionExt = regions,
#'     typeS1 = sample1, typeS2 = sample2)
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @export
getSeqCMP <- function(HLAInfo, regionExt, typeS1, typeS2) {

    splitS1 <- splitTyping(typeS1)
    splitS2 <- splitTyping(typeS2)

    posS1 <- getTypingPos(HLAInfo$HLAAlignment, splitS1)
    posS2 <- getTypingPos(HLAInfo$HLAAlignment, splitS2)

    if (splitS1[1] != splitS2[1]) {
        stop("The typing for the two alleles must be for the same HLA gene.")
    }

    refSeq <- HLAInfo$refSeq[[splitS1[1]]]
    posInit <- HLAInfo$posInit[[splitS1[1]]]

    seqCMP <- list(refSeq="", seqS1="", seqS2="")
    seqCMP$refSeq <- getSubSeq(refSeq, posInit, regionExt)
    seqCMP$seqS1 <- getSubSeq(HLAInfo$HLAAlignment[posS1]$SeqDiff,
                                posInit, regionExt)
    seqCMP$seqS2 <- getSubSeq(HLAInfo$HLAAlignment[posS2]$SeqDiff,
                                posInit, regionExt)

    return(seqCMP)
}

#' @title Get the subsequence from a selected sequence for a specific region.
#'
#' @description From a selected sequence, extract a subsequence using the list
#' of specified regions.
#'
#' @param seq a \code{character} string that represent the sequence used to
#' extract the subsequence.
#'
#' @param posInit a \code{integer} that represent the starting position of the
#' gene.
#'
#' @param regionExt a \code{data.frame} containing a \code{start}
#' and an \code{end} column. The \code{start} and \code{end} positions are both
#' included in the extracted region.
#'
#' @return a \code{character} string containing the extracted subsequence.
#'
#'
#' @examples
#'
#' ## Sequence used for extraction
#' sequence <- "    MGSGWVPWVVALLVNLTRLDSSMTQGTDSPEDFVIQAKADCYFTNGTEKVQFVVRFIF"
#'
#' ## The starting position of the gene
#' position <- -10
#'
#' ## The region to extact
#' regions <- data.frame(start=c(15, 25), end=c(20, 30))
#'
#' ## Extract subsequence
#' HLAClustRView:::getSubSeq(seq = sequence, posInit = position,
#'     regionExt = regions)
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal
getSubSeq <- function(seq, posInit, regionExt) {

    subSeq <- ""

    for(i in seq_len(nrow(regionExt))) {
        subSeq <- paste0(subSeq, substr(seq,
                            regionExt$start[i] - posInit,
                            regionExt$end[i] - posInit))
    }

    return(subSeq)
}


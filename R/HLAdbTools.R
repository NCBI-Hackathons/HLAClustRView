#' @title Get the sequence from the type
#'
#' @description Get from the object HLAdb two sequences and the reference
#' for a region
#'
#' @param HLAInfo TODO
#'
#' @param regionExt TODO
#'
#' @param typeS1 Sample 1 TODO
#'
#' @param typeS2 Sample 2 TODO
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

    splitS1 <- splitTyping(typeS1)
    splitS2 <- splitTyping(typeS2)

    posS1 <- getTypingPos(HLAInfo$HLAAlignment, splitS1)
    posS2 <- getTypingPos(HLAInfo$HLAAlignment, splitS2)

    if(splitS1[1] != splitS2[1]){
        stop("Call get seq with type from 2 genes")
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
#' @importFrom data.table data.table rbindlist
#' @export
getSubSeq <- function(seq, posInit, regionExt){

    subSeq <- ""

    for(i in seq_len(nrow(regionExt))){
        subSeq <- paste0(subSeq, substr(seq,
                            regionExt$start[i] - posInit,
                            regionExt$end[i] - posInit))
    }

    return(subSeq)
}


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
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom data.table data.table rbindlist
#' @export


getSeqCMP <- function(HLAInfo, regionExt, typeS1, typeS2){


    splitS1 <- splitTyping(typeS1)
    splitS2 <- splitTyping(typeS2)

    posS1 <- getIncompleteTypingPos(HLAInfo$HLAAlignment, splitS1)
    posS1 <- reduceTypingPos(HLAInfo$HLAAlignment, posS1)
    posS2 <- getIncompleteTypingPos(HLAInfo$HLAAlignment, splitS2)
    posS2 <- reduceTypingPos(HLAInfo$HLAAlignment, posS2)

    if(splitS1[1] != splitS2[1]){
        stop("Call get seq with type from 2 genes")
    }
    if(is.na(posS1) || is.na(posS2)){
        stop(paste0("Typing without specific sequence ", typeS1, " ", typeS2))
    }

    refSeq <- HLAInfo$refSeq[[splitS1[1]]]
    posInit <- HLAInfo$posInit[[splitS1[1]]]

    seqCMP <- list(refSeq="", seqS1="", seqS2="")
    seqCMP$refSeq <- getSubSeq(refSeq, posInit, regionExt)
    if(!(is.na(posS1))){
        seqCMP$seqS1 <- getSubSeq(HLAInfo$HLAAlignment[posS1]$SeqDiff, posInit, regionExt)
    }
    if(!(is.na(posS2))){
        seqCMP$seqS2 <- getSubSeq(HLAInfo$HLAAlignment[posS2]$SeqDiff, posInit, regionExt)
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
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom data.table data.table rbindlist
#' @export


getSubSeq <- function(seq, posInit, regionExt){

    subSeq <- ""

    for(i in seq_len(nrow(regionExt))){
        subSeq <- paste0(subSeq, substr(seq, regionExt$start[i] - posInit, regionExt$end[i] - posInit))
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
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom data.table data.table rbindlist
#' @export

parseSubMatrix <- function(fileName){

    subMatrix <- read.table(fileName, header = TRUE)

}

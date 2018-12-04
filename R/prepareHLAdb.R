#' @title Parse HLA Database alignment files for a specific alignment type
#'
#' @description Parse HLA Database alignment files that are present in a
#' directory, as specified by input, an
#' generate a class object that can be used in further analysis. The function
#' only parse one type of aligment file at the time. There is 3 types of
#' alignment files that can be parsed: CDS sequence, genomic and protein.
#'
#' Beware that the names of the alignment files should not be changed as
#' the name is used to identify the gene that is currently parsed.
#'
#' @param hlaDbPath A \code{character} string, the path to the directory of
#' the alignment files from
#' \url{http://hla.alleles.org/alleles/text_index.html}. The directory must
#' exist and contain at least one HLA alignment file.
#'
#' @param seqType A \code{character} string, the sequence type of the file to
#' parse. The choices are: "nuc" for CDS sequence, "gen" for genomic, and
#' "prot" for protein. Default: "nuc".
#'
#' @return An object of class \code{HLAdb} with 3 entries. The entries are:
#' \itemize{
#' \item \code{refSeq} A \code{list} of \code{character} string that
#' represent reference sequences; one sequence
#' per HLA gene.
#' \item \code{posInit} A \code{list} of \code{integer};
#' one starting position per HLA gene.
#' \item \code{HLAAlignment} A \code{data.table} containing the information
#' for each allele of each HLA gene.
#' }
#'
#' @details See \url{http://hla.alleles.org/alleles/text_index.html}
#'
#' @examples
#'
#' ## Get path where some HLA database files are stored
#' directory <- system.file("extdata", package = "HLAClustRView")
#'
#' ## Parse HLA database files of protein type
#' HLAInfo <- parseHLADbAlignment(hlaDbPath=directory, seqType="prot")
#'
#' ## Show reference sequences
#' HLAInfo$refSeq
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom data.table data.table rbindlist
#' @export
parseHLADbAlignment <- function(hlaDbPath, seqType=c("nuc", "gen", "prot")) {

    ## Validate that the sequence type is known
    if(! seqType %in% c("nuc", "gen", "prot")) {
        stop(paste0("Not validate sequence type parameter for ",
                        "parseHLADbAlignment: ", seqType))
    }

    ## Validate that the directory is a string
    if (!is.character(hlaDbPath)) {
        stop("The hlaDbPath parameter must by a character string")
    }

    ## Validate that the directory exists
    if (!dir.exists(hlaDbPath)) {
        stop("The hlaDbPath parameter must by a valid directory")
    }

    files <- dir(path = hlaDbPath, pattern = paste0("*_", seqType, ".txt"))

    ## Validate that there is at least one file in the directory
    if(length(files) == 0) {
        stop(paste0("There must be at least one alignment file in the ",
                    "hlaDbPath directory"))
    }

    refSeq <- list()
    posInit <- list()
    HLAAlignment <- list()

    for(fileName in files) {
        tmp <- parseAlignment(paste0(hlaDbPath, "/", fileName))
        geneName <- gsub(paste0("_", seqType, ".txt"), "", fileName)
        refSeq[[geneName]] <- tmp$refSeq
        posInit[[geneName]] <- tmp$posInit
        HLAAlignment[[geneName]] <- tmp$HLAalignment
    }

    ## Create object to return
    HLAdb <- list(refSeq=refSeq, posInit=posInit,
                  HLAAlignment=rbindlist(HLAAlignment))
    class(HLAdb) <- "HLAdb"
    return(HLAdb)
}

#' @title Pre-process HLA Database alignment file
#'
#' @description Extract information from one HLA alignment file.
#'
#' @param fileName A \code{character} string, the name of the alignment file
#' from \url{http://hla.alleles.org/alleles/text_index.html}.
#'
#' @return A \code{list} containing:
#' \itemize{
#' \item refSeq A \code{character} string that represent the reference
#' sequence.
#' \item posInit A \code{integer} that represent the starting position of the
#' gene.
#' \item HLAalignment A \code{data.table} containing the information for all
#' types associated to the HLA gene.
#' }
#'
#' @details See \url{http://hla.alleles.org/alleles/text_index.html}
#'
#' @examples
#'
#' ## Set the name of the aligment file to be processed
#' fileInfo <- paste0(system.file("extdata", package = "HLAClustRView"),
#'     "/DRA_prot.txt")
#'
#' ## Parse aligment file
#' HLAClustRView:::parseAlignment(fileName=fileInfo)
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom data.table data.table
#' @keywords internal
parseAlignment <- function(fileName) {
    # Position where the sequence start
    startPosFile <- 20
    maxTyping <- 50000

    allLines <- readLines(fileName)
    listPos <- NULL
    seqType <- NULL

    # Find which type of alignment Prot, cDNA or gDNA
    for (curType in c("Prot", "cDNA", "gDNA")) {

        listPos <- grep(curType, allLines)
        if (length(listPos) > 0) {
            seqType <- curType
            break;
        }
    }

    # Get the position of the reference
    offSet <- ifelse(seqType == "cDNA", 1,0)
    cpt <- listPos[1] + 2 + offSet

    while (cpt < maxTyping) {
        if (regexpr("[A-Z]", allLines[cpt], perl=TRUE)[1] == -1) {
            break
        }
        cpt <- cpt + 1
    }

    if (cpt == maxTyping) {
        stop(paste0("The program reach the max of typing in ", fileName), "\n")
    }

    nbType <- cpt - listPos[1] + 2 + offSet

    # The reference sequence
    refSeq <- ""

    # data.table of each type with a representation of the diffence
    # between the sequence and the reference.
    HLAalignment <- data.table(GeneName=character(nbType),
                               AlleleGroup=character(nbType),
                               Protein=character(nbType),
                               SynSubst=character(nbType),
                               Noncoding=character(nbType),
                               Suffix=character(nbType),
                               SeqDiff=character(nbType))

    # Loop on the position in allLine of the sequence type
    # before the reference sequence
    # The reference sequence position is + 2 + offSet
    for (startLine in listPos) {
        # Get the initial position of the sequence
        i <- 1
        s <- substr(allLines[startLine + offSet], startPosFile, startPosFile)
        posInit <- ""
        while (s != " " && s != "") {
            posInit <- paste0(posInit, s)
            s <- substr(allLines[startLine + offSet], startPosFile + i,
                        startPosFile + i)
            i<-i+1
        }
        posInit <- as.integer(posInit)

        # The reference sequence
        startSeq <- startLine + 2 + offSet
        parseRef <- extractRef(allLines[startSeq], startPosFile)
        refSeq <- paste0(refSeq, parseRef$refSeq)

        # Parse the type of the reference sequence
        if(startLine == listPos[1]){
            nameTyping <- extractTyping(allLines[startSeq], startPosFile)
            curTyping <- splitTyping(nameTyping)
            HLAalignment$GeneName[1] <- curTyping[1]
            HLAalignment$AlleleGroup[1] <- curTyping[2]
            HLAalignment$Protein[1] <- curTyping[3]
            HLAalignment$SynSubst[1] <- curTyping[4]
            HLAalignment$Noncoding[1] <- curTyping[5]
            HLAalignment$Suffix[1] <- curTyping[6]
        }

        HLAalignment$SeqDiff[1] <- paste0(HLAalignment$SeqDiff[1],
                                                    parseRef$seqDiff)

        cpt <- startSeq + 1

        # loop on all type for thee current portion of the
        # sequence
        while (cpt < maxTyping) {
            seqCur <- allLines[cpt]
            if (regexpr("[A-Z]", seqCur, perl=TRUE)[1] >= 1) {
                nameTyping <- extractTyping(seqCur, startPosFile)
                # parse the typing
                curTyping <- splitTyping(nameTyping)
                # get the position of the typing in
                # HLAalignment
                typingPos <- getTypingPos(HLAalignment, curTyping)

                if (length(typingPos) == 0) {
                    # If the typing is not initialise
                    HLAalignment$GeneName[cpt - startSeq + 1] <- curTyping[1]
                    HLAalignment$AlleleGroup[cpt - startSeq + 1] <- curTyping[2]
                    HLAalignment$Protein[cpt - startSeq + 1] <- curTyping[3]
                    HLAalignment$SynSubst[cpt - startSeq + 1] <- curTyping[4]
                    HLAalignment$Noncoding[cpt - startSeq + 1] <- curTyping[5]
                    HLAalignment$Suffix[cpt - startSeq + 1] <- curTyping[6]
                    HLAalignment$SeqDiff[cpt - startSeq + 1] <- paste0(
                        HLAalignment$SeqDiff[cpt - startSeq + 1],
                        extractSeq(seqCur, startPosFile))

                } else if (length(typingPos) == 1) {
                    # if the typing is initialize just
                    # add the diff seq
                    HLAalignment$SeqDiff[typingPos] <- paste0(
                        HLAalignment$SeqDiff[typingPos],
                        extractSeq(seqCur, startPosFile))
                }
            } else {
                # stop the loop if line without [A-Z]
                break
            }
            cpt <- cpt + 1
        }

        # Just to prevent infinite loop
        if(cpt == maxTyping){
            stop(paste0("The program reach the max of typing in ",
                        fileName, "\n"))
        }
    }

    HLAGene <- list(refSeq=refSeq, posInit=posInit, HLAalignment=HLAalignment)

    return(HLAGene)
}

#' @title Process the line contening the reference sequence
#'
#' @description Extract the information about the reference sequence.
#'
#' @param seq A \code{character} string contening the reference sequence.
#'
#' @param startPos A \code{integer} corresponding to the position in the string
#' where the sequence is starting.
#'
#' @return A \code{list} with two entries. The first is \code{refSeq} which is
#' containing the sequence of the reference. The second is seqDiff is a
#' \code{string} representing the alignment format for the current sequence.
#'
#' @examples
#'
#' ## String with the sequence information
#' sequence <- " DOB*01:01:01:01       MGSGWV PWVVALLVNL TRLDSSMTQG TDSPEDFVIQ"
#'
#' ## Extract sequence information
#' HLAClustRView:::extractRef(seq=sequence, startPos=20)
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal
extractRef <- function(seq, startPos) {

    tmpSeq <- substr(seq, startPos, nchar(seq))

    refSeq <- ""
    seqDiff <- ""
    flag <- TRUE

    for (i in seq_len(nchar(tmpSeq))) {
        b <- substr(seq, startPos - 1 + i, startPos - 1 + i)
        if (flag) {
            refSeq <- paste0(refSeq, b)

            if (b != " ") {
                seqDiff <- paste0(seqDiff, "-")
                flag <- FALSE
            } else {
                seqDiff <- paste0(seqDiff, " ")
            }
        } else {
            if (!(b %in% c(" ", "|"))) {
                refSeq <- paste0(refSeq, b)
                if (regexpr("[A-Z]", b, perl=TRUE)[1] >= 1) {
                    seqDiff <- paste0(seqDiff, "-")
                } else {
                    seqDiff <- paste0(seqDiff, b)
                }
            }
        }
    }
    return(list(refSeq=refSeq, seqDiff=seqDiff))
}

#' @title Process the line containing the alignment to the reference
#'
#' @description Extract the information about the alignment compared to the
#' reference sequence.
#'
#' @param seq A \code{character} string containing the alignment sequence.
#'
#' @param startPos A \code{integer} corresponding to the position in the
#' \code{character} string where the sequence is starting.
#'
#' @return A \code{character} string representing the aligment sequence
#' compared to the reference sequence.
#'
#' @examples
#'
#' ## Define a alignment sequence string
#' sequence <- " DOB*01:01:01:02       ------ ---------- ---------- "
#'
#' ## Extract the alignment section
#' HLAClustRView:::extractSeq(seq=sequence, startPos=20)
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal
extractSeq <- function(seq, startPos) {

    tmpSeq <- substr(seq, startPos, nchar(seq))

    seqDiff <- ""
    flag <- TRUE

    for (i in seq_len(nchar(tmpSeq))) {
        b <- substr(seq, startPos - 1 + i, startPos - 1 + i)
        if (flag) {
            seqDiff <- paste0(seqDiff, b)
            if (b != " ") {
                flag <- FALSE
            }
        } else {
            if (!(b %in% c(" ", "|"))) {
                seqDiff <- paste0(seqDiff, b)
            }
        }
    }

    return(seqDiff)
}


#' @title Extract the HLA typing string from one line of aligment
#'
#' @description Extract the HLA typing string from a character string
#' containing one line of alignment
#'
#' @param seq A \code{character} string contening the sequence.
#'
#' @param endPos A \code{integer} fixing the end of the introduction
#' section of the sequence.
#'
#' @return a \code{character} string containing the HLA type of a empty
#' \code{character} string
#'
#' @examples
#'
#' ## One line of aligment
#' sequence <- " DOB*01:01:01:01       MGSGWV PWVVALLVNL TRLDSSMTQG"
#'
#' ## Extract HLA type
#' HLAClustRView:::extractTyping(seq=sequence, endPos=20)
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal
extractTyping <- function(seq, endPos) {

    tmpSeq <- substr(seq, 2, endPos)

    nameTyping <- ""
    for (i in seq_len(nchar(tmpSeq))) {
        b <- substr(tmpSeq, i, i)
        if (b != " ") {
            nameTyping <- paste0(nameTyping, b)
        } else {
            break
        }
    }

    return(nameTyping)
}


#' @title Get the position in the HLA table of the corresponding HLA
#' typing
#'
#' @description Get the position in the HLA table of the corresponding HLA
#' typing. The HLA typing can use up to 6 fields for the identification.
#'
#' @param seq A \code{character} string containing the sequence.
#'
#' @param curTyping A \code{vector} of \code{character} string defining
#' the HLA type. The \code{vector} should have 6 fields.
#'
#' @return A \code{integer} representing the position, in the data.table, of
#' the HLA typing. If the HLA typing is not present, an empty \code{integer}
#' is returned.
#'
#' @examples
#'
#' ## Information about HLA alignment file
#' fileInfo <- paste0(system.file("extdata", package = "HLAClustRView"),
#'     "/DRA_prot.txt")
#'
#' ## Parse HLA alignment file
#' ## DRAInfo$HLAalignment contains the HLA table
#' DRAInfo <- HLAClustRView:::parseAlignment(fileName=fileInfo)
#'
#' ## Create a vector with the information about one HLA typing
#' typing <- matrix(data=c("DRA", "01", "01", "02", NA, NA), nrow=1)
#'
#' ## Get the position of the specific typing in the HLA table
#' HLAClustRView:::getTypingPos(seqProcess=DRAInfo$HLAalignment,
#'     curTyping=typing)
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal
getTypingPos <- function(seqProcess, curTyping) {

    curPos <- integer()

    if(is.na(curTyping[4])){
        curPos <- which(seqProcess$GeneName == curTyping[1] &
                        seqProcess$AlleleGroup == curTyping[2] &
                        seqProcess$Protein == curTyping[3] &
                        is.na(seqProcess$SynSubst))
    } else if(is.na(curTyping[5])){
        curPos <- which(seqProcess$GeneName == curTyping[1] &
                        seqProcess$AlleleGroup == curTyping[2] &
                        seqProcess$Protein == curTyping[3] &
                        seqProcess$SynSubst == curTyping[4] &
                        is.na(seqProcess$Noncoding) )
    } else if(is.na(curTyping[6])){
        curPos <- which(seqProcess$GeneName == curTyping[1] &
                            seqProcess$AlleleGroup == curTyping[2] &
                            seqProcess$Protein == curTyping[3] &
                            seqProcess$SynSubst == curTyping[4] &
                            seqProcess$Noncoding == curTyping[5] &
                            is.na(seqProcess$Suffix) )
    } else {
        curPos <- which(seqProcess$GeneName == curTyping[1] &
                        seqProcess$AlleleGroup == curTyping[2] &
                        seqProcess$Protein == curTyping[3] &
                        seqProcess$Noncoding == curTyping[4] &
                        seqProcess$Noncoding == curTyping[6] &
                        seqProcess$Suffix == curTyping[7] )
    }

    return(curPos)
}


#' @title Extract HLA typing information for all samples present in a text file
#'
#' @description Parses a text file containing HLA allele typing for one or
#' multiple samples and transforms it into an
#' HLADb object.
#'
#' @param hlaDbPath a \code{character} string, the filepath to a text file
#' containing the information about the samples typing. The name must
#' correspond to an existing text file.
#'
#' @return a \code{data.frame} containing the following column names:
#' \itemize{
#' \item \code{SampleName} The name of the sample.
#' \item \code{AllelName} The name of the allele (1 or 2).
#' \item \code{GeneName} The name of the HLA gene.
#' \item \code{AlleleGroup} The section identifying the subtype group.
#' \item \code{Protein} The section identifying the nucleotide substitutions
#' group.
#' \item \code{SynSubst} The section identifying the synonymous nucleotide
#' substitutions group.
#' \item \code{NonCoding} The section identifying the non-coding
#' polymorphisms group.
#' \item \code{Suffix} The suffix of the HLA typing or \code{NA}.
#' }
#'
#' @examples
#'
#' ## Get path where some HLA database files are stored
#' directory <- system.file("extdata", package = "HLAClustRView")
#' fileName <- paste0(directory, "/Samples_HLA_typing.txt")
#'
#' ## Parse file to extract HLA typing for all samples
#' ex <- readHLADataset(hlaDbPath=fileName)
#'
#' ## Show the first lines of the output dataset
#' head(ex)
#'
#'
#' @author Adewunmi Adelaja
#'
#' @importFrom stringr str_split
#' @importFrom utils read.table
#' @export
readHLADataset <-function (hlaDbPath)
{
    # GOAL: load information from HLA database and create a usable R
    #  object that contains it
    # Input validation:
    #     must be valid database path
    # minimum expected files must be present
    db <- data.frame()

     #check if RDS file is present
    #folderName <- dirname(hlaDbPath)
    #rdsFileName <-paste(folderName,'hladb.rds', sep = '/')
    #if (!file.exists(rdsFileName))
    #{

    ## Validate that the parameter is a character string
    if (!is.character(hlaDbPath)) {
        stop("The hlaDbPath parameter must by a character string")
    }

    ## Validate that the parameter is a valid file
    if (!file.exists(hlaDbPath)) {
        stop(paste0("The file '", hlaDbPath, "' must be a valid text file"))
    }

    #if (file.exists(hlaDbPath))
    #{
    data <- read.table(hlaDbPath, header = TRUE, stringsAsFactors = FALSE)
    rownames(data) <- data[,1]
    data <- data[,2:ncol(data)]

    tData <- t(data)
    rNames <- rownames(tData)
    #split alleles numbers from gene name
    nameList <- rep(str_split(rNames, '_'), times = ncol(tData))
    sampleNames <- rep(colnames(tData), each = nrow(tData))
    geneNames <-matrix(unlist(nameList),  ncol = 2, byrow=TRUE)
    colName <-c('SampleName', 'AlleleName', 'GeneName', 'AlleleGroup',
                    'Protein', 'SynSubst', 'NonCoding', 'Suffix')

     #split the hla digits
    typing <- matrix(unlist(lapply(tData, splitTyping)), ncol = 6, byrow = TRUE)
    db <- data.frame(cbind(matrix(sampleNames), geneNames[,2], typing),
                        stringsAsFactors = FALSE)
    colnames(db) <- colName

    #create HlADb class
    #hla <-setClass('HLADb',
    #              contains = c("data.frame"))
    #db<-hla(db)

    res<-list()
    res[[1]]<- db
    class(res) <- "HLADb"

    #saveRDS(db, file= rdsFileName)
    #}
    #    else{
    #        print("Can't find hlaDbPath.")
    #    }
    #
    #}else
    #{
    #    db <-readRDS(rdsFileName)
    #}

    #
    #remove samples with NA allele groups
    drop <- db[is.na(db$AlleleGroup), 'SampleName']
    db <- db[!db$SampleName %in% drop,]

    return(db)
}

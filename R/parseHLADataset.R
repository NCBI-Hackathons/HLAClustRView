#' @title Extract HLA typing information for all samples present in a text file
#'
#' @description Parses a text file containing HLA allele typing for one or
#' multiple samples and transforms it into an object.
#'
#' @param hlaFilePath a \code{character} string, the filepath to a text file
#' containing the information about the samples typing. The name must
#' correspond to an existing text file. The text file must have a header.
#'
#' @return a \code{list} of class \code{HLADataset}
#' containing the following elements:
#' \itemize{
#' \item \code{data} a \code{tibble} object containing the HLA typing
#' information for all samples. The columns are:
#' \itemize{
#' \item \code{SampleName} a \code{character} string that represent the
#' name of the sample.
#' \item \code{AllelName} a \code{character} string that represent the
#' name of the allele (1 or 2).
#' \item \code{GeneName} a \code{character} string that represent the
#' name of the HLA gene.
#' \item \code{AlleleGroup} a \code{character} string that represent the
#' section identifying the subtype group.
#' \item \code{Protein} a \code{character} string that represent the
#' section identifying the nucleotide substitutions group.
#' \item \code{SynSubst} a \code{character} string that represent the
#' section identifying the synonymous nucleotide substitutions group.
#' \item \code{NonCoding} a \code{character} string that represent the
#' section identifying the non-coding polymorphisms group.
#' \item \code{Suffix} a \code{character} string that represent the
#' suffix of the HLA typing or \code{NA}.
#' }
#' }
#'
#' @examples
#'
#' ## Get path where some HLA database files are stored
#' directory <- system.file("extdata", package = "HLAClustRView")
#' fileName <- paste0(directory, "/Samples_HLA_typing.txt")
#'
#' ## Parse file to extract HLA typing for all samples
#' ex <- readHLADataset(hlaFilePath=fileName)
#'
#' ## Show the output HLADataset object
#' print(ex)
#'
#'
#' @author Adewunmi Adelaja
#'
#' @importFrom stringr str_split
#' @importFrom utils read.table
#' @export
readHLADataset <- function(hlaFilePath)
{
    ## Validate that the parameter is a character string
    if (!is.character(hlaFilePath)) {
        stop("The hlaFilePath parameter must by a character string")
    }

    ## Validate that the parameter is a valid file
    if (!file.exists(hlaFilePath)) {
        stop(paste0("The file '", hlaFilePath, "' must be a valid text file"))
    }

    data <- read.table(hlaFilePath, header = TRUE, stringsAsFactors = FALSE)
    rownames(data) <- data[,1]
    data <- data[, 2:ncol(data)]

    tData <- t(data)
    rNames <- rownames(tData)
    #split alleles numbers from gene name
    nameList <- rep(str_split(rNames, '_'), times = ncol(tData))
    sampleNames <- rep(colnames(tData), each = nrow(tData))
    geneNames <- matrix(unlist(nameList),  ncol = 2, byrow=TRUE)
    colName <- c('SampleName', 'AlleleName', 'GeneName', 'AlleleGroup',
                    'Protein', 'SynSubst', 'NonCoding', 'Suffix')

    ## Split the hla digits
    typing <- matrix(unlist(lapply(tData, splitTyping)), ncol = 6, byrow = TRUE)
    db <- data.frame(cbind(matrix(sampleNames), geneNames[, 2], typing),
                        stringsAsFactors = FALSE)
    colnames(db) <- colName

    ## Remove samples with NA allele groups
    drop <- db[is.na(db$AlleleGroup), 'SampleName']
    db <- db[!db$SampleName %in% drop, ]

    ## Format data.frame into tibble
    db <- parse_hla_dataset(db)

    ## Prepare HLADataset object
    res <- list()
    res[["data"]] <- db
    class(res) <- "HLADataset"

    return(res)
}


#' @title Convert data.frame with HLA typing information to tibble object
#'
#' @description Converts a data.frame that contains the HLA typing informatiOn
#' to a tibble object. It also removes samples with missing allele group
#' information.
#'
#' @param hladata a \code{data.frame} with the HLA typing information from
#' all samples.
#'
#' @return a \code{tibble} object with the HLA information for each sample.
#'
#' @examples
#'
#' ## Create a data.frame as demo
#' demo <- data.frame(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
#'     AlleleName=c(1, 2, 1, 2), GeneName=c("A", "A", "A", "A"),
#'     AlleleGroup=c("02", "02", "03", "03"), Protein=c("01", "01", "01", "02"),
#'     SynSubst=c("01", "02", "01", "01"), NonCoding=c("01", "01", NA, NA),
#'     Suffix=c(NA, NA, NA, NA))
#'
#' ## Convert the data.frame to a tibble object
#' HLAClustRView:::parse_hla_dataset(demo)
#'
#' ## Create a data.frame with missing information
#' demoMissing <- data.frame(SampleName=c("DEMO1", "DEMO1", "DEMO2", "DEMO2"),
#'     AlleleName=c(1, 2, 1, 2), GeneName=c("A", "A", "A", "A"),
#'     AlleleGroup=c("02", "02", NA, "03"), Protein=c("01", "01", "01", "02"),
#'     SynSubst=c("01", "02", "01", "01"), NonCoding=c("01", "01", NA, NA),
#'     Suffix=c(NA, NA, NA, NA))
#'
#' ## Convert the data.frame to a tibble object
#' HLAClustRView:::parse_hla_dataset(demoMissing)
#'
#'
#' @author Santiago Medina
#'
#' @importFrom purrr map reduce set_names
#' @importFrom dplyr mutate_all filter pull bind_cols as_tibble %>%
#' @importFrom rlang .data
#' @keywords internal
parse_hla_dataset <- function(hladata) {

    hla_data <-
        hladata %>%
        map(as_tibble) %>%
        reduce(bind_cols) %>%
        set_names(names(hladata)) %>%
        mutate_all(.funs = as.character)

    missing_samples <-
        hla_data %>%
        filter(is.na(.data$AlleleGroup)) %>%
        pull(.data$SampleName)

    return(hla_data %>%
               filter(!.data$SampleName %in% missing_samples))
}

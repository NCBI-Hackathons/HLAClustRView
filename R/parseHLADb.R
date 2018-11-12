#' @title
#'parseHLADb
#' @description
#'Parses a text file containing HLA allele names into an HLADb object
#' @param hlaDbPath
#'Filepath to a text file
#' @return
#' HLADb object with the following column names:
#' "SampleName"  "AlleleName"  "GeneName"    "AlleleGroup" "Protein"     "SynSubst"
#' "NonCoding"   "Suffix"
#' @examples
#'ex <- parseHLADb("pathToFile.txt")
#'head(ex)
#'# SampleName AlleleName GeneName AlleleGroup Protein SynSubst NonCoding Suffix
#'# 1  ERR188053          1        A          31      01       02        01   <NA>
#'# 2  ERR188053          2        A          68      01       01        01   <NA>
#'# 3  ERR188053          1        B          27      05       02      <NA>   <NA>
#'# 4  ERR188053          2        B          27      05       02      <NA>   <NA>
#'# 5  ERR188053          1        C          02      02       02        01   <NA>
#'# 6  ERR188053          2        C          02      02       02        01   <NA>
#'
#'
#' @author Adewunmi Adelaja
#'
#' @importFrom stringr str_split
#' @importFrom utils read.table
#' @importFrom methods setClass
#' @export
parseHLADb <-function (hlaDbPath)
{    # GOAL: load information from HLA database and create a usable R object that contains it
    # Input validation:
    #     must be valid database path
    # minimum expected files must be present
    db <- data.frame
    options(stringsAsFactors = FALSE)

     #check if RDS file is present
    folderName <- dirname(hlaDbPath);
    rdsFileName <-paste(folderName,'hladb.rds', sep = '/')
    if (!file.exists(rdsFileName))
    {

    if (file.exists(hlaDbPath))
    {
    data <- read.table(hlaDbPath, header = TRUE, stringsAsFactors = FALSE)
    rownames(data) <- data[,1]
    data <- data[,2:ncol(data)]

    tData <- t(data)
    rNames <- rownames(tData)
    #split alleles numbers from gene name
    nameList <- rep(str_split(rNames, '_'), times = ncol(tData))
    sampleNames <- rep(colnames(tData), each = nrow(tData))
    geneNames <-matrix(unlist(nameList),  ncol = 2,byrow=T)
    colName <-c('SampleName', 'AlleleName', 'GeneName', 'AlleleGroup','Protein', 'SynSubst','NonCoding','Suffix')

     #split the hla digits
    typing <-matrix(unlist(lapply(tData, splitTyping)),ncol = 6, byrow =T)
    db <- data.frame(cbind(matrix(sampleNames), geneNames[,2],typing))
    colnames(db)<-colName

    #create HlADb class
    #hla <-setClass('HLADb',
    #              contains = c("data.frame"))
    #db<-hla(db)

    res<-list()
    res[[1]]<- db
    class(res) <- "HLADb"

    #saveRDS(db, file= rdsFileName)
    }
        else{
            print("Can't find hlaDbPath.")
        }

    }else
    {
        db <-readRDS(rdsFileName)
    }

    #
    #remove samples with NA allele groups
    drop <- db[is.na(db$AlleleGroup), 'SampleName']
    db <- db[db$SampleName != drop,]
    return (db)
}

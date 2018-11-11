parseHLADb <-function (hlaDbPath)
{
    options(stringsAsFactors = FALSE)
    library (stringr)

    #check if RDS file is present
    #if RDS file is not present, load file

    setwd('D:/GoogleDrive/grad_school/conferences/CSH/2018/hackathon/results/')
    fileName <- paste(getwd(), 'E_GEUVexp.txt',sep = '/')
    data <- read.table(fileName, header = TRUE, stringsAsFactors = FALSE)
    rownames(data) <- data[,1]
    data <- data[,2:ncol(data)]


    #Name pattern :
    #Add NA if digit is not prsent
    tData <- t(data)


    rNames <- rownames(tData)
    #split alleles numbers from gename
    nameList <- str_split(rNames, '_')
    sampleNames <- colnames(tData)
    geneName <- data.frame(matrix(unlist(nameList), nrow = length(geneNames), byrow=T),stringsAsFactors=FALSE)

    #split the hla digits
    apply(tData, 1, splitTyping)
     # prefix <- apply(tData,1, function (row) str_extract(row, '(^.*[:digit:])|(^[:alpha:]+$)'))
    # splitted <- sapply(prefix, function (row) str_split(row, '[.*:_]'))
    splitted <- apply(prefix, 1,function (row) str_split(row, '[.*:_]'))

    splitted <- apply(prefix,function (row) str_split(row, '[.*:_]'))

    suffix <- apply(tData,1, function (row) str_match(row, '(?<=([:digit:]|[.*:_]))[:alpha:]+$'))
    # splitted <- sapply(tData, function (row) str_split(row, '[.*:_]'))

    # pad to common length
    # splitMat <-sapply(splitted, matrix)
    # splitMat <-sapply(splitted, matrix)
    # splitMat <-matrix(splitted, )
    # elemSz <-apply(splitted,1,length)
    # # elemSz <-lapply(splitted,lengths)
    # maxSz <-max(elemSz)
    # pad <- apply(elemSz,2,function(x)  maxSz-x)
    colName <-c('SampleName', 'AlleleName', 'GeneName', 'AlleleGroup','Protein', 'SynSubst','Noncoding')
    db<-  data.frame(matrix(ncol = 8, nrow = nrow(prefix)*ncol(prefix)))

    db <-  data.frame(matrix(ncol = 8, nrow = nrow(prefix)*ncol(prefix)))
    db<-  data.frame(matrix(ncol = 8, nrow = nrow(prefix)))
    i <-1
    for (r in 1:nrow(tData))
    {
        for (c in 1:ncol(tData))
        {

            splitted <- strsplit(prefix[r,c],'[.*:_]')
            pad <-5-lengths(splitted)
            dataRow <-c(sampleNames[r],
                        geneName[c,2],
                       splitted,
                        rep("NA", pad),
                        array(suffix[r,c]))
            db[i,] <-dataRow
             i<-i+1
        }
    }





    # db<-  data.frame(matrix(ncol = 8, nrow = nrow(prefix)))
    # i <-1
    # for (r in 1:nrow(prefix))
    # {
    #     for (c in 1:ncol(prefix))
    #     {
    #         # db[i,] <-c(sampleNames[r], geneName[c,2], splitMat[[r,c]], rep("NA", pad[r,c]))
    #         # splitMat[r,c]<-c(sampleNames[r], geneNames[plitMat[[r,c]], rep("NA", pad[r,c]))
    #
    #         splitted <- strsplit(prefix[r,c],'[.*:_]')
    #         pad <-5-lengths(splitted)
    #         dataRow <-c(sampleNames[r],
    #                     geneName[c,2],
    #                    splitted,
    #                     rep("NA", pad),
    #                     array(suffix[r,c]))
    #         db[i,] <-dataRow
    #          i<-i+1
    #     }
    # }
    #
    #
    # nonSampleCols <- setdiff(colnames(data), 'Sample')
    # names<-str_split(data[,nonSampleCols], '[*.:]')


    #create HLADb class using setClass()
    setClass("HLADb", representation = ())

    # R object assigned HLADb class
    #
    #       List of sequence type
    #       Protein
    #       - List of HLA Genes
    #       - Gene A
    #       - Reference (string)
    #       - data.table with all alleles
    #       - Gene B
    #       - Reference (string)
    #       - data.table with all alleles
    #       DNA
    #       RNA
    # overwrite public function print(class)

    setMethod("print","HLADb", function (x){
        cat(paste9)
    } )
    HLADbObj <-

    return HLADbObj
}



# GOAL: load information from HLA database and create a usable R object that contains it
# Input validation:
#     must be valid database path
# minimum expected files must be present

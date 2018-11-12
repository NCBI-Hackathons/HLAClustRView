#' Convert distance table to distance matrix
#'
#' This function takes the output of \code{calculateSimilarity} function
#' and returns a distance matrix to use for clustering
#' @param outMet tibble with hamming distance for sample pairs \code{calculateSimilarity}
#'
#' @return
#' @export
#'
#' @examples
make_distance_matrix <- function(outMet) {
    nbCase <- length(unique(outMet$SampleName1)) + 1
    matDia <- matrix(unlist(sapply(seq_len(nbCase),
                                   FUN=function(x,outMetric, nbCase){
                                       l<-NULL
                                       if(x < nbCase){
                                           l <- c(rep(0,x),
                                                  outMetric$HammingDistance[
                                                      ((x-1) * (nbCase-1) - ((x-1)*x)/2 + 1):
                                                          (x * (nbCase-1) - (x*(x+1))/2 +1)])
                                       } else{
                                           l <- rep(0,x)
                                       }

                                       return(l)
                                   }, outMetric=outMet,nbCase=nbCase)), nc=nbCase )

    # add row names and column names to the distance matrix
    nameMat <- c(unique(unlist(outMet[,1])), unique(unlist(outMet[,2]))[length(unique(unlist(outMet[,2])))])
    rownames(matDia) <- nameMat
    colnames(matDia) <- nameMat
    matDia

}

#' Draw cluster
#'
#' @param matrixDistance distance matrix
#' @importFrom dplyr %>%
#' @importFrom graphics plot
#' @importFrom stats as.dist
#' @return
#' @export
#'
#' @examples
draw_cluster <- function(matrixDistance) {
    hc <- as.dist(matrixDistance) %>%
        hclust()
    plot(hc)
}

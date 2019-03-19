#' @title Hierarchical clustering using a HLA typing similarity metric
#'
#' @description Hierarchical cluster analysis on a set of similarity values
#' related to HLA typing.
#'
#' @param hlaMetric an object of class \code{HLAMetric}.
#'
#' @param method a string \code{character} string, the agglomeration method
#' to be used for clustering.
#'
#' @param \ldots arguments passed to or from other methods.
#'
#' @return an object of class \code{hclust} which describes the tree produced
#' by the clustering process. See \link[stats]{hclust}.
#'
#'
#' @examples
#'
#' ## Load example dataset
#' data(demoHLADataset)
#'
#' ## Calculate Hamming distance metric
#' hammingDistance <- calculateHamming(demoHLADataset)
#'
#' ## Run clustering
#' hclustHLA(hammingDistance)
#'
#' @author Astrid Deschenes
#'
#' @importFrom stats hclust
#' @export
hclustHLA <- function(hlaMetric, method="complete", ...) {


    distance <- as.dist(hlaMetric$dist)
    hclust(distance, method=method, ...)
}

#' @title Hierarchical clustering using a HLA typing similarity metric
#'
#' @description Hierarchical cluster analysis on a set of similarity values
#' related to HLA typing.
#'
#' @param dist an object of class \code{dist}. See \link[stats]{dist}.
#'
#' @param method a string \code{character} string, the agglomeration method
#' to be used for clustering.
#'
#' @return an object of class \code{hclust} which describes the tree produced
#' by the clustering process. See \link[stats]{hclust}.
#'
#' @param \ldots arguments passed to or from other methods.
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Astrid Deschenes
#'
#' @importFrom stats hclust
#' @export
createCluster <- function(dist, method="complete", ...) {
    hclust(dist, method=method, ...)
}

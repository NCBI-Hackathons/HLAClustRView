
#' @title Draw dendrogram
#'
#' @description This function creates and draws a dendrogram using the
#' metric distances from an \code{hlaMetric} object.
#'
#' @param hlaMetric an object of class \code{hlaMetric} that contain the
#' distance between samples.
#'
#' @param as.phylo a \code{logical} indicating if the dendrogram has to
#' be drawn using a phylogenetic tree display. Default: \code{FALSE}.
#'
#' @param \ldots arguments passed to the \code{plot} method.
#'
#' @return \code{NULL}
#'
#' @examples
#'
#' ## Load example dataset
#' data(demoHLADataset)
#'
#' ## Calculate Hamming distance metric
#' hammingMetric <- calculateHamming(demoHLADataset)
#'
#' ## Create dendrogram
#' draw_dendrogram(hammingMetric)
#'
#' @author Santiago Medina, Astrid Deschênes, Pascal Belleau
#'
#' @importFrom dplyr %>%
#' @importFrom graphics plot
#' @importFrom ape  as.phylo
#' @importFrom stats as.dendrogram
#' @export
draw_dendrogram <- function(hlaMetric, as.phylo=FALSE, ...) {

    myArgs <- match.call()

    if (!("HLAMetric" %in% class(hlaMetric))) {
        stop("hlaMetric must be of class \"HLAMetric\"")
    }

    if (as.phylo) {
        pp <- as.dist(hlaMetric$dist) %>%
                hclust() %>% as.phylo()
    } else {
        dd <- as.dist(hlaMetric$dist) %>%
                hclust() %>% as.dendrogram()
    }

    if (as.phylo) {
        plot(x=pp, ...)
    } else {
        plot(x=dd, ...)
    }
}


#' @title Draw sample-to-sample heatmap
#'
#' @description Create and draw a sample-to-sample heatmap using the
#' metric distances.
#'
#' @param hlaMetric an object of class \code{hlaMetric} that contain the
#' distance between samples.
#'
#' @param \ldots arguments passed to the \code{Heatmap} method
#'
#' @return A \code{Heatmap-class} object.
#'
#' @examples
#'
#' ## Load example dataset
#' data(demoHLADataset)
#'
#' ## Calculate Hamming distance metric
#' hammingMetric <- calculateHamming(demoHLADataset)
#'
#' ## Create sample-to-sample heatmap
#' draw_heatmap(hammingMetric)
#'
#' @author Astrid Deschênes, Pascal Belleau
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @export
draw_heatmap <- function(hlaMetric, ...) {

    myArgs <- match.call()

    if (!("HLAMetric" %in% class(hlaMetric))) {
        stop("hlaMetric must be of class \"HLAMetric\"")
    }

    distMat <- hlaMetric$dist

    newMat <- t(distMat)
    newMat[lower.tri(newMat)] <- distMat[lower.tri(distMat)]

    if ("name" %in% names(myArgs)) {
        if ("heatmap_legend_param" %in% names(myArgs)) {
            heatGraph <- Heatmap(newMat, ...)
        } else {
            heatGraph <- Heatmap(newMat,
                        heatmap_legend_param = list(direction="horizontal"),
                        ...)
        }
    } else {
        if ("heatmap_legend_param" %in% names(myArgs)) {
            heatGraph <- Heatmap(newMat, name = hlaMetric$metric, ...)
        } else {
            heatGraph <- Heatmap(newMat, name = hlaMetric$metric,
                        heatmap_legend_param = list(direction="horizontal"),
                        ...)
        }
    }


    draw(heatGraph, heatmap_legend_side = "bottom")

    return(invisible(heatGraph))
}



#' @title Draw dendrogram
#'
#' @description Create and draw a dendrogram using the
#' metric distances.
#'
#' @param hlaMetric an object of class \code{hlaMetric} that contain the
#' distance between samples.
#'
#' @return TODO
#'
#' @examples
#'
#' ## Load example dataset
#' data(demoHLADataset)
#'
#' ## Calculate Hamming distance metric
#' hammingMetric <- calculateHamming(demoHLADataset)
#'
#' ## Create dendogram
#' draw_dendogram(hammingMetric)
#'
#' @author Santiago Medina, Astrid Deschenes, Pascal Belleau
#'
#' @importFrom dplyr %>%
#' @importFrom graphics plot
#' @importFrom stats as.dist
#' @export
draw_dendogram <- function(hlaMetric) {

    if (!("HLAMetric" %in% class(hlaMetric))) {
        stop("hlaMetric must be of class \"HLAMetric\"")
    }


    hc <- as.dist(hlaMetric$dist) %>%
        hclust()
    plot(hc)
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
#' @author Astrid Deschenes, Pascal Belleau
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

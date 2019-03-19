#' @rdname HLAMetric
#'
#' @title  Print a \code{HLAMetric} object
#'
#' @method print HLAMetric
#'
#' @description Print a \code{HLAMetric} object
#'
#' @param x the output object from \code{calculateHamming}
#' function to be printed
#'
#' @param \ldots arguments passed to or from other methods
#'
#' @return an object of class
#' \code{HLAMetric}
#'
#' @examples
#'
#' ## Load example dataset
#' data(demoHLADataset)
#'
#' ## Calculate Hamming distance metric
#' metric <- calculateHamming(demoHLADataset)
#'
#' ## Print information about metric
#' print(metric)
#'
#' @export
print.HLAMetric <- function(x, ...) {
    # Print title before printing the content of the object
    cat("HLA Dataset\n\n")

    distTemp <- x[["dist"]]
    distType <- x[["metric"]]
    if (!is.null(distType)) {

        cat(paste0("Metric: ", distType, "\n"))
    }
    if (!is.null(distTemp)) {
        cat(paste0("Distance matrix: \n"))
        print(distTemp)
    } else {
        cat("\tNo HLA metric data.")
    }

    invisible(x)
}

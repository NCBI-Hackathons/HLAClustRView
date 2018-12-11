#' @rdname HLADataset
#'
#' @title  Print a \code{HLADataset} object
#'
#' @method print HLADataset
#'
#' @description Print a \code{HLADataset} object
#'
#' @param x the output object from \code{readHLADataset}
#' function to be printed
#'
#' @param \ldots arguments passed to or from other methods
#'
#' @return an object of class
#' \code{HLADataset}
#'
#' @examples
#'
#' ## Load dataset
#' ## TODO
#'
#' @export
print.HLADataset <- function(x, ...) {
    # Print title before printing the content of the object
    cat("HLA Dataset\n\n")

    dataTemp <- x[["data"]]
    if (!is.null(dataTemp)) {
        cat(paste0("Number of Samples: ", length(unique(dataTemp$SampleName)),
                    "\n"))
        cat(paste0("Number of Genes: ", length(unique(dataTemp$GeneName)),
                    "\n\n"))
        print.data.frame(dataTemp, max=100)
    } else {
        cat("\tNo HLA typing data.")
    }

    invisible(x)
}



#' @title Draw cluster
#'
#' @description TODO
#'
#' @param matrixDistance distance matrix
#'
#' @return TODO
#'
#' @examples
#'
#' ##TODO
#'
#' @author Santiago Medina, Pascal Belleau
#'
#' @importFrom dplyr %>%
#' @importFrom graphics plot
#' @importFrom stats as.dist
#' @export
draw_cluster <- function(matrixDistance) {
    hc <- as.dist(matrixDistance) %>%
        hclust()
    plot(hc)
}

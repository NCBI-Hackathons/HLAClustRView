#' HLAClustRView: HLA typing clustering and visualization based on specific similarity metrics
#'
#' This package implements specialized similarity metrics that
#' quantify the similarity between HLA typing from multiple samples. Using
#' these metrics, the package enables hierarchical cluster analysis and can
#' bind RNA-seq expression (when available) to the clusters.
#'
#' @docType package
#'
#' @name HLAClustRView-package
#'
#' @aliases HLAClustRView-package HLAClustRView
#'
#' @author  Pascal Belleau, Adewunmi Adelaja, Astrid Deschenes,
#' Santiago Medina, Nissim Ranade and Allissa Dillman
#'
#' Maintainer:
#' Pascal Belleau <belleau@@cshl.edu>
#'
#' @seealso

#' \itemize{

#'     \item \code{\link{readHLADataset}} {for extracting HLA typing
#'     information for all samples present in a text file}
#'     \item \code{\link{calculateHamming}} {for calculating the Hamming
#'     distance metric between all samples}
#'     \item \code{\link{parseHLADbAlignment} {for parsing HLA Database
#'     alignment files for a specific alignment type}}
#' }
#'
#' @keywords package
NULL

#' This dataset can be used TODO.
#'
#' @name hladb
#'
#' @aliases hladb
#'
#' @docType data
#'
#' @format TODO
#'
#' @return TODO
#'
#' @usage TODO
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## TODO
#'
NULL

#' The HLA typing information from multiple samples stored in a
#' \code{HLADataset} format (for demo purpose).
#'
#' This dataset can be used to test \code{calculateHamming} and
#' TODO(AUTRE METRIC) functions.
#'
#' @name demoHLADataset
#'
#' @aliases demoHLADataset
#'
#' @docType data
#'
#' @format TODO
#'
#' @return a \code{list} of class \code{HLADataset}
#' containing the following elements:
#' \itemize{
#' \item \code{data} a \code{tibble} object containing the HLA typing
#' information for all samples. The columns are:
#' \itemize{
#' \item \code{SampleName} a \code{character} string that represent the
#' name of the sample.
#' \item \code{AllelName} a \code{character} string that represent the
#' name of the allele (1 or 2).
#' \item \code{GeneName} a \code{character} string that represent the
#' name of the HLA gene.
#' \item \code{AlleleGroup} a \code{character} string that represent the
#' section identifying the subtype group.
#' \item \code{Protein} a \code{character} string that represent the
#' section identifying the nucleotide substitutions group.
#' \item \code{SynSubst} a \code{character} string that represent the
#' section identifying the synonymous nucleotide substitutions group.
#' \item \code{NonCoding} a \code{character} string that represent the
#' section identifying the non-coding polymorphisms group.
#' \item \code{Suffix} a \code{character} string that represent the
#' suffix of the HLA typing or \code{NA}.
#' }}
#'
#' @usage data(demoHLADataset)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Load example dataset
#' data(demoHLADataset)
#'
#' ## Calculate hamming distance metric
#' calculateHamming(demoHLADataset)
#'
#'
NULL

#' This dataset contains the allele and digit1 information for
#' HLA genes in one pair of samples.
#'
#' This dataset can be used to test the \code{sample_pair_distance} function.
#'
#' @name example_sample_pair_data
#'
#' @aliases example_sample_pair_data
#'
#' @docType data
#'
#' @format A \code{tibble} object with 10 rows and 4 variables:
#' \itemize{
#'     \item{SampleName}{sample id}
#'     \item{GeneName}{hla gene id}
#'     \item{AlleleName}{allele information, 1 or 2}
#'     \item{AlleleGroup}{allele group data}
#' }
#'
#' @return TODO
#'
#' @usage data(example_sample_pair_data)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## TODO
#'
NULL


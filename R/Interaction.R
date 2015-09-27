#' Interaction class to hold results from mode classification with the edgeR
#' limma pipelines
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{interaction.tbl}:}{Data.frame containing mode classification
#'    and other data from edgeR or limma}
#'    \item{\code{modes}:}{Factor vector of length n where each element is an
#'    interaction mode classification that can mapped to a gene.}
#'    \item{\code{outcome.vectors}:}{Matrix of dimensions i x 6 where each row
#'    represents a 6-dimensional outcome vector corresponding to the possible
#'    comparisons when gene expression outputs from 4 conditions are compared}
#'    \item{\code{outcome.vectors.s}:}{Character vector of length n where each
#'    element is an outcome vector as a string}
#'    \item{\code{voom}:}{EList object}
#'    \item{\code{d}:}{DGEList object}
#'    \item{\code{design}:}{list with design matrix and contrasts from
#'    makeContrasts()}
#'    }
#'  @exportClass Interaction
setClass("Interaction",
         representation(
           interaction.tbl = "data.frame",
           modes = "factor",
           outcome.vectors = "matrix",
           outcome.vectors.s = "character",
           voom = "EList",
           d = "list",
           design = "list")
)

#' initalize an Interaction object from a list generated with the limma or edge classifier
#' functions
#'
#' @param l list with elements seen in the function definition.
#' @return Interaction object
#' @export initInteraction
initInteraction <- function(l) {
  this_combo <- new("Interaction", interaction.tbl = l$interaction.tbl,
                    modes = l$modes, outcome.vectors = l$outcome.vectors,
                    outcome.vectors.s = l$outcome.vectors.s,
                    voom = l$voom, d = l$d, design = l$design)
  return(this_combo)
}

#' code pairwise comparison outcome vectors with limma
#'
#' @param eb limma::eBayes()
#' @param design design matrix initialized with model.matrix
#' @param DRUG_COMBINATION Boolean TRUE when conditions of interaction are vehicle, signal
#' x, signal y, and signal x + y respectively
#' @return matrix where matrix[,1:6] are the six elements of the outcome vectors and
#' matrix[,7:10] are the log fold changes from limma::topTable()
#' @examples
#' library(limma)

codeLimmaOutcomeVectors <- function(eb, design) {
  # differentially expressed genes ----------------------------------------------
  top_tables <- list()
  for (i in 1:ncol(design)) {
    top_tables[[i]] <- topTable(eb, coef = i, number = nrow(v$E), confint = ci,
                                p.value = 1, sort.by = "none", adjust.method = "fdr")
  }
  # Class Classification ----------------------------------------------------------------
  # log fold changes and confidence interval bounds for contrasts of interest
  # conditions - [Vehicle][Signal X][Signal Y][Signal X + Y]
  if (DRUG_COMBINATION == TRUE) {
    x <- top_tables[[which(colnames(design) == T1)]][,6:8]
    o <- top_tables[[which(colnames(design) == T0)]][,6:8]
    y <- top_tables[[which(colnames(design) == T2)]][,6:8]
    xy <- top_tables[[which(colnames(design) == T12)]][,6:8]
    fold_changes <- cbind(top_tables[[which(colnames(design) == T0)]]$logFC,
                          top_tables[[which(colnames(design) == T1)]]$logFC,
                          top_tables[[which(colnames(design) == T2)]]$logFC,
                          top_tables[[which(colnames(design) == T12)]]$logFC)
  }
  # conditions - [None:F][None:M][Signal X:F][Signal X:M]
  if (DRUG_COMBINATION == FALSE) {
    x <- top_tables[[which(colnames(design) == "None:F")]][,6:8]
    o <- top_tables[[which(colnames(design) == "None:M")]][,6:8]
    y <- top_tables[[which(colnames(design) == paste(T1, "F", sep = ":"))]][,6:8]
    xy <- top_tables[[which(colnames(design) == paste(T1, "M", sep = ":"))]][,6:8]
    fold_changes <- cbind(top_tables[[which(colnames(design) == "None:F")]]$logFC,
                          top_tables[[which(colnames(design) == "None:M")]]$logFC,
                          top_tables[[which(colnames(design) == paste(T1, "F", sep = ":"))]]$logFC,
                          top_tables[[which(colnames(design) == paste(T1, "M", sep = ":"))]]$logFC)
  }

  # code outcome vectors --------------------------------------------------------
  x_o <- codeOVElement(x, o)
  y_o <- codeOVElement(y, o)
  xy_o <- codeOVElement(xy, o)
  y_x <- codeOVElement(y, x)
  xy_x <- codeOVElement(xy, x)
  xy_y <- codeOVElement(xy, y)

  # combine elements of each comparison to form outcome vector
  limma_outcome_vectors <- cbind(x_o, y_o, xy_o, y_x, xy_x, xy_x)
  limma_outcome_vectors <- cbind(validateOutcomes(limma_outcome_vectors), fold_changes)
  return(limma_outcome_vectors)
}

#' Interaction mode classification with matrix of outcome vectors and log fold changes
#'
#' @param res matrix where columns 1:6 are the elements of the outcome vectors and columns
#' 7:length(res) are the log fold changes for interaction class classification
#' @param pipeline character vector of length n = 1 indicating pipeline to use for
#' interaction mode classification ("edgeR" or "limma")
#' @return character vector of length m = nrow(res) indicating interaction modees
#' @export modeClassification
#' @examples
#' \dontrun{
#' library(limma)
#' ov <- codeLimmaOutcomeVectors(eb, design)
#' limma_modes <- modeClassification(ov, pipeline = "limma")
#' }

modeClassification <- function(res, pipeline) {
  if (pipeline == "edgeR") {
    modes <- classifyByThOutcomeVector(res, "edgeR")
  }
  if (pipeline == "limma") {
    modes <- classifyByThOutcomeVector(res, "limma")
  }
  return(modes)
}


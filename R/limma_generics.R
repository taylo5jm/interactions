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


#' Code elements of pairwise outcome vectors with confidence intervals from limma
#'
#' @param condition1 matrix where conditions 1:3 are the log fold changes
#' @param condition2 ""
#' @return numeric vector of length n where each element represents the outcome of the ith
#' pairwise comparison
#' @export codeOVElement

codeOVElement <- function(condition1, condition2) {
  con <- cbind(condition1, condition2)
  assignOutcome <- function(v) {
    if (v[1] > v[5] & v[1] < v[6] & v[4] > v[2] & v[4] < v[3]) {
      return(0)
    }
    if (v[1] < v[5] & v[4] > v[3]) {
      return(-1)
    }
    else {
      return(1)
    }
  }
  apply(con, 1, assignOutcome)
}

#knnplot.R

library(GGally)

#' Make scatter matrix from PCA scores
#'
#' @param modes character vector of length n indicating the modes of genes 1:n
#' @param log.fc matrix of log fold changes where ncol(log.fc) == 4 if pipeline == "limma" and
#' ncol(log.fc) == 6 if pipeline == "edgeR"
#' @return PCA scatter matrix generated with GGally::ggpairs
#' @export pcaScatterMatrix
#' @examples
#' library(GGally)
#' \dontrun{
#' scatter_matrix <- pcaScatterMatrix(limma_modes, logfc)
#' }
pcaScatterMatrix <- function(modes, log.fc) {
    int_fold_changes <- fold_changes[which(limma_modes != "A" & limma_modes != "NI"),]
    int_modes <- limma_modes[which(limma_modes != "A" & limma_modes != "NI")]

    limma_pca <- princomp(int_fold_changes, cor = T)
    pca_res <- data.frame(limma_pca$scores, int_modes)
    if (pipeline == "limma") {
      colnames(pca_res) <- c(colnames(pca_res)[-5], "Mode")
      pca_cols <- 4
    }
    if (pipeline == "edgeR") {
    colnames(pca_res) <- c("XvC", "YvC", "XYvC", "YvX", "XYvX",
                           "XYvY", "Mode")
    pca_cols <- 6
    }
    return(ggpairs(pca_res, columns = 1:pca_cols, colour = 'Mode'))
}

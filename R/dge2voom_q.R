
#' calculate normalization factors with edgeR and voom with limma to model as
#' approximately normal
#'
#' @param counts numeric matrix of counts where rows are genes and columns are libraries
#' @param genes data.frame of gene metadata
#' @param null Boolean TRUE if columns of counts matrix should be permuted to
#'        simulate null distribution
#' @return \code{d} DGEList with normalization factors calculated and
#' \code{voom(d)} is returned locally
#' @export dge2voom_q
#' @examples
#' \dontrun{
#' # v <- dge2voom_q(counts, genes)
#' # d <- estimateCommonDisp
#' }

dge2voom_q <- function(counts, genes, null = FALSE) {
  
  # permute columns of counts matrix for null distribution
  if (null == TRUE) {
    counts <- counts[,sample(ncol(counts))]
  }
  d <- DGEList(counts, genes = genes, remove.zeros = TRUE)
  d <- massageCounts(d, cpm_min = 2, sample_min = 5)
  d <- calcNormFactors(d)
  output <- list(voom = voom(d),
                 dgelist = d)
  input <- list(raw.counts = counts,
                genes = genes)
  list(output = output,
       input = input)
}



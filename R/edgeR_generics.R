#' Fit generalized linear model and make pairwise comparisons based on contrasts
#'
#' @param d DGEList object with normalization factors calculated
#' @param design design matrix initialized with model.matrix
#' @return d DGEList with common, trended and tagwise dispersion estimated
#' @export estimateAllGLMDisp
#' @examples
#' # dge2voom_q returns d, a DGEList with normalization factors
#' # calculated to the global environment
#' \dontrun{v <- dge2voom_q(counts, genes)
#' d <- estimateAllGLMDisp(d, design)
#' }

estimateAllGLMDisp <- function(d, design) {
  if (is.list(design) & length(design) == 2) {
    this_design <- design[[1]]
    contrasts <- design[[2]]
  }
  d <- estimateGLMCommonDisp(d, this_design) %>%
       estimateGLMTrendedDisp(this_design) %>%
       estimateGLMTagwiseDisp(this_design)
  return(d)
}

#' make design matrix and contrasts based on whether signal_combination is 2 drugs or M/F/Dug
#'
#' @importFrom limma makeContrasts
#' @param metadata data.frame with metadata retrieved for 4 conditions
#' @param drug_combination Boolean TRUE when signal x, signal y, and signal x + y
#' are conditions 2, 3, and 4 respectively
#' @param pipeline character vector of length 1 indicating what pipeline is
#' being used to classify interaction modes for each gene
#' @return list where list[[1]] is the design matrix and list[[2]]
#' represent the contrasts created with limma::make.contrasts if pipeline == "edgeR"
#' if pipeline == "limma", then a design matrix is returned.
#' @export designMatrixAndContrasts

designMatrixAndContrasts <- function(metadata, pipeline,
                                     DRUG_COMBINATION = TRUE) {

  if (pipeline == "edgeR") {
    if (DRUG_COMBINATION == TRUE) {
      metadata$treatment_name <- factor(metadata$treatment_name)
      metadata$donor_id <- factor(metadata$donor_id)
      colnames(d$counts) <- paste0(metadata$treatment_name, metadata$donor_id)
      this_design <- model.matrix(~0 + treatment_name + donor_id, data = metadata)
      
      # nrows <- which(metadata$treatment_name == T1) %>% length
      
      nrows <- unique(metadata$treatment_name) / 4
      this_contrasts <- makeContrasts(X-C, Y-C, XY-C, Y-X, XY-X, XY-Y,
                                      levels = c("C", "X", "Y", "XY")) %>%
                        rbind(matrix(0, 
                              nrow = length(unique(metadata$donor_id)) - 1,
                              ncol = 6))
      this_list <- list(design_matrix = this_design,
                        contrasts = this_contrasts)
      return(this_list)
    }
    
    if (DRUG_COMBINATION == FALSE) {
      this_design <- model.matrix(~0 + treatment_name:donor_sex, data = metadata)
      colnames(this_design) <- c("X.F", "C.F", "X.M", "C.M")
      nrows <- which(metadata$treatment_name == T1) %>% length - 1
      this_contrasts <- makeContrasts(X.F-C.F, C.M-C.F, X.M-C.F, C.M-X.F, X.M-X.F, X.M-C.M,
                                      levels = c("C.F", "X.F", "C.M", "X.M"))
      this_list <- list(design_matrix = this_design,
                        contrasts = this_contrasts)
      return(this_list)
    }

  else {
    return("Failed to make design matrix and/or contrasts")
  }
}

  if (pipeline == "limma") {
    if (DRUG_COMBINATION == TRUE) {
      this_design <- model.matrix(~0 + treatment_name + donor_id, data = metadata)
      colnames(this_design) <- c(levels(factor(metadata$treatment_name)), levels(factor(donors[-1])))
      return(this_design)
    }

    if (DRUG_COMBINATION == FALSE) {
      this_design <- model.matrix(~0 + treatment_name:donor_sex, data = metadata)
      colnames(this_design) <- sub("treatment_name", "", colnames(this_design))
      colnames(this_design) <- sub("donor_sex", "", colnames(this_design))
      return(this_design)
    }
  }
}

#' generalized linear model and likelihood ratio tests with edgeR
#'
#' @param d DGEList with dispersions estimated
#' @param design list where design[[1]] is the design matrix constructed by model.matrix
#' and design[[2]] are the contrasts constructed by make.contrasts
#' @param contrasts optional constrasts constructed by make.contrasts
#' @return pairwise list of likelihood ratio test tables generated with
#' edgeR::topTags. 
#' @export glmFitAndLRT
#' @examples
#' \dontrun{
#' d <- dge2voom_q(counts, genes)
#' design <- designMatrixAndContrasts(metadata, drug_combination = TRUE,
#' pipeline = "edgeR")
#' my_lrt_res <- glmFitAndLRT(d, design)
#' }

glmFitAndLRT <- function(d, design, contrasts) {
  # designMatrixAndContrasts() can output a list with both
  if (is.list(design) & length(design) == 2) {
    this_design <- design$design_matrix
    contrasts <- design$contrasts
  }
  # fit generalized linear model
  fit <- glmFit(d, this_design)
  pairwise <- list()
  # likelihood ratio tests for each element in the outcome vector
  for (i in 1:6) {
    pairwise[[i]] <- glmLRT(fit, contrast = contrasts[,i])
  }
  pairwise <- lapply(pairwise, topTags, n = nrow(d$counts), 
                     adjust.method = "fdr", sort.by = "none")
  return(pairwise)
}


#' Code pairwise comparison outcome vectors with edgeR likelilood ratio tests p-values
#' and log fold changes
#'
#' @param pairwise list of results from glmFitAndLRT()
#' @param alpha numeric vector of length n = 1 indicating alpha
#' @return res matrix where res[,1:6] are the elements of the pairwise comparison outcome vectors
#' and res[,7:12] are the log-fold-change ratios from the 6 contrasts
#' @export codeEdgeROutcomeVectors
#' @examples
#' library(edgeR)
#' library(limma)


codeEdgeROutcomeVectors <- function(pairwise, alpha) {
  codeOutcomeVector <- function(lrt, alpha) {
    lrt_res <- lrt$table
    assignOutcome <- function(v, alpha) {
      if (v[4] < alpha) {
        return(sign(v[1]))
      }
      else {
        return(0)
      }
    }
    apply(lrt_res, 1, assignOutcome, alpha)
  }
  my_res <- lapply(pairwise, codeOutcomeVector, alpha)
  edgeR_outcome_vectors <- mapply(cbind, my_res)
  edgeR_pw_s <<- apply(edgeR_outcome_vectors, 1, stringCoerce)

  getLogFC <- function(lrt_list) {
    l <- lapply(lrt_list, function(x) (x$table$logFC))
    return(mapply(cbind, l))
  }

  logfc <- lapply(pairwise, function(x) (x$table$logFC))
  logfc <- mapply(cbind, logfc)
  res <- matrix(c(edgeR_outcome_vectors, logfc), nrow = nrow(logfc),
                ncol = 12)
  return(res)
}

#' check column names of design matrix
#'
#' @param design design matrix initialized with model.matrix
#' @return design matrix with validated column names
#' @export designColnames
#' @examples
#' \dontrun{
#' this_design <- designMatrixAndContrasts(metadata, drug_combination = TRUE,
#' pipeline = "edgeR")
#' this_design <- designColnames(design)
#' }
designColnames <- function(design) {
  treatment_cols <- grep("treatment_name", colnames(design))
  treatment_cols <- colnames(design)[treatment_cols]
  new_treatment_cols <- lapply(treatment_cols,
                               function(x) (strsplit(x, "treatment_name"))) %>% unlist
  new_treatment_cols <- new_treatment_cols[new_treatment_cols != ""]

  genericTreatmentColnames <- function(treatment_name) {
    if (treatment_name == T0) {
      return("C")
    }
    if (treatment_name == T1) {
      return("X")
    }
    if (treatment_name == T2) {
      return("Y")
    }
    if (treatment_name == T12) {
      return("XY")
    }
    else {
      return("ERROR")
    }
  }
  design_colnames <- sapply(new_treatment_cols, genericTreatmentColnames)
  colnames(design)[1:length(treatment_cols)] <- design_colnames
  return(design)
}




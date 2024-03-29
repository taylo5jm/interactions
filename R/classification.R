# classification.R

#' Interaction mode classification with matrix of counts and data.frames for genes and metadata
#'
#' @param counts numeric matrix of counts where rows are genes and columns are libraries
#' @param genes data.frame of gene metadata
#' @param metadata data.frame with metadata retrieved for 4 conditions
#' @param design design matrix constructed with model.matrix
#' @param drug_combination Boolean TRUE when signal x, signal y, and signal x + y
#' are conditions 2, 3, and 4 respectively
#' @param ci Confidence level for confidence intervals generated by limma
#' @param null Boolean TRUE if null distribution should be simulated by permuting columns
#' of counts matrix
#' @return list of interaction modes & classes, outcome vectors, EList (voom), lmFit() object,
#' eBayes() object, design matrix and summary table
#' @export limmaModeClassifier
limmaModeClassifier <- function(counts, genes, metadata,
                                design,
                                DRUG_COMBINATION = TRUE,
                                ci = 0.99, null = FALSE, ...) {
  # permute columns of counts matrix for random
  if (null == TRUE) {
    counts <- counts[,sample(ncol(counts))]
  }

  # initialize DGEList (edgeR) and EList (limma voom)
  v <- dge2voom_q(counts, genes)

  # guess design matrix with donor_id as covariates
  if (DRUG_COMBINATION == TRUE & missing(design)) {
    design <- model.matrix(~0 + treatment_name + donor_id, data = metadata)
    colnames(design) <- c(levels(factor(metadata$treatment_name)),
                          paste("Donor", levels(factor(metadata$donor_id))[-1], sep = "_"))
  }

  fit <- lmFit(v, design)
  eb <- eBayes(fit)
  top_tables <- list()

  for (i in 1:ncol(design)) {
      top_tables[[i]] <- limma::topTable(eb, coef = i, number = nrow(v$E), confint = ci,
                                  p.value = 1, sort.by = "none", adjust.method = "fdr")
    }

  # Class Classification ----------------------------------------------------------------
  # log fold changes and confidence interval bounds for contrasts of interest
  if (DRUG_COMBINATION == TRUE) {
    x <- top_tables[[which(colnames(design) == T1)]][,6:8]
    o <- top_tables[[which(colnames(design) == T0)]][,6:8]
    y <- top_tables[[which(colnames(design) == T2)]][,6:8]
    xy <- top_tables[[which(colnames(design) == T12)]][,6:8]
  }

  if (DRUG_COMBINATION == FALSE) {
    x <- top_tables[[which(colnames(design) == "None:F")]][,6:8]
    o <- top_tables[[which(colnames(design) == "None:M")]][,6:8]
    y <- top_tables[[which(colnames(design) == paste(T1, "F", sep = ":"))]][,6:8]
    xy <- top_tables[[which(colnames(design) == paste(T1, "M", sep = ":"))]][,6:8]
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
  rm(x_o, y_o, xy_o, y_x, xy_x, xy_y)

  # as.character for easy tabulation
  limma_pw_s <- apply(limma_outcome_vectors, 1, stringCoerce)

  # logFC from limma::topTable (control, signal x, signal y, signal xy)
  fold_changes <- cbind(o$logFC, x$logFC,
                        y$logFC, xy$logFC)

  limma_outcome_vectors <- cbind(limma_outcome_vectors, fold_changes)
  limma_classes <- apply(fold_changes, 1, getIntClass, "limma")
  limma_modes <- classifyByThOutcomeVector(limma_outcome_vectors, pipeline = "limma")

  limma_classes_exc <- limma_classes[which(limma_modes != "NI" &
                                             limma_modes != "A" & limma_modes != "UC")]

  # summary table ---------------------------------------------------------------
  limma_interaction_tbl <- sitt::interactionGeneTable(v, limma_modes, fold_changes,
                                                limma_pw_s, limma_classes_exc,
                                                DRUG_COMBINATION = TRUE)

  return(list(interaction.tbl = limma_interaction_tbl,
              classes = limma_classes,
              modes = limma_modes,
              outcome.vectors = limma_outcome_vectors,
              outcome.vectors.s = limma_pw_s,
              voom = v,
              fit = fit,
              eBayes = eb,
              design = design))
}

#' Interaction mode classification with edgeR likelihood ratio tests
#'
#' @param counts numeric matrix of counts where rows are genes and columns are 
#'  libraries
#' @param genes data.frame of gene metadata
#' @param metadata data.frame with metadata retrieved for 4 conditions
#' @param alpha numeric vector indicating what alpha should be used to
#' code pairwise comparison outcome vectors
#' @param null Boolean TRUE if null distribution should be simulated by 
#' permuting columns of counts matrix
#' @return list of interaction modes & classes, outcome vectors, EList (voom),
#'  lmFit() object, eBayes() object, design matrix and summary table
#' @export edgeRModeClassifier

edgeRModeClassifier <- function(counts, genes, metadata, design, alpha = 0.01,
                                DRUG_COMBINATION = TRUE) {
  
  pairwise <- glmFitAndLRT(d, design, contrasts)
  # pairwise <- lapply(pairwise, topTags, n = nrow(d$counts),
  #                  adjust.method = "fdr", sort.by = "none")
  edgeR_outcome_vectors <- codeEdgeROutcomeVectors(pairwise, alpha)
  edgeR_modes <- classifyByThOutcomeVector(edgeR_outcome_vectors, 
                                           pipeline = "edgeR")

  # interaction modes
  edgeR_interaction_tbl <- interactionGeneTable(v, edgeR_modes,
                                                edgeR_outcome_vectors[,7:12],
                                                edgeR_pw_s, T0 = T0,
                                                T1 = T1, T2 = T2,
                                                T12 = T12) # summary table
  return(list(metadata = metadata,
              counts = counts,
              genes = genes,
              interaction.tbl = edgeR_interaction_tbl,
              modes = edgeR_modes,
              outcome.vectors = edgeR_outcome_vectors,
              outcome.vectors.s = edgeR_pw_s,
              voom = v,
              DGEList = d,
              design = this_design))
}

#' Interaction mode classification with limma and edgeR
#'
#' @param counts numeric matrix of counts where rows are genes and columns are 
#'  libraries
#' @param genes data.frame of gene metadata
#' @param metadata data.frame with metadata retrieved for 4 conditions
#' @param alpha numeric vector indicating what alpha should be used to
#' code pairwise comparison outcome vectors
#' @param null Boolean TRUE if null distribution should be simulated by 
#' permuting columns of counts matrix
#' @return nested list
#' @export interactionModeClassifier
#'
interactionModeClassifier <- function(data, 
                                      DRUG_COMBINATION = TRUE,
                                      ci = 0.99, alpha = 0.01) {
  
  .limma <- limmaModeClassifier(counts, genes, metadata, design,
                                DRUG_COMBINATION = DRUG_COMBINATION,
                                ci = ci)
  
  .edgeR <- edgeRModeClassifier(counts, genes, metadata, design,
                                DRUG_COMBINATION = DRUG_COMBINATION,
                                alpha = alpha)
  
  
  .null_limma <- limmaModeClassifier(counts, genes, metadata, design,
                                     DRUG_COMBINATION = DRUG_COMBINATION,
                                     ci = ci)

  .null_edgeR <- edgeRModeClassifier(counts, genes, metadata, design,
                                     DRUG_COMBINATION = DRUG_COMBINATION,
                                     alpha = alpha)
  return(list(
            limma = .limma,
            null.limma = .null_limma,
            edgeR = .edgeR,
            null.edgeR = .null_edgeR))
}


#' Classify a gene into an interaction mode given an outcome vector and log-fold changes
#'
#' @param ov_logfc vector of length 1:6 + n, where elements 1:6 are the elements of the pairwise comparison outcome vector,
#' n is a positive integer that satisfies n %% 2 == 0 and elements 6 : 6 + n are log fold changes associated with
#' each gene
#' @param character vector indicating which pipeline should be used to code outcome vectors
#' @return character vector of length 1 indicating interaction mode or lack thereof (null vector)
#' @export classifyByThOutcomeVector

 classifyByThOutcomeVector <- function(ov_logfc, pipeline = "limma") {
   
  # enumerate theoretical outcome vectors ------------------------------------
  th.vectors <- outcomeVectorsByMode()
  
  # mathematically defined interaction modes ---------------------------------
  MODES <<- allModes()
  
  # classify modes -------------------------------------------------
  modeClassifier <- function(ov_logfc, th.vectors, pipeline) {
    ov <- ov_logfc[1:6]
    logfc <- ov_logfc[7:length(ov_logfc)]
    
    # additive vectors ----------------------------------------------
    if (isMode(ov, th.vectors$sym_left)) {
      return("Sym.Left")
    }

    if (isMode(ov, th.vectors$sym_right)) {
      return("Sym.Right")
    }

    if (isMode(ov, th.vectors$step_down)) {
      return("Step.Down")
    }

    if (isMode(ov, th.vectors$step_up)) {
      return("Step.Up")
    }

    # null vector ---------------------------------------------------
    if (isNi(ov)) {
      return("NI")
    }
    # anomalous vector ---------------------------------------------
    if (isAnomaly(ov, th.vectors)) {
      return("A")
    }

    # Negative Interaction, XY < X + Y ---------------------------------------
    if(getIntClass(logfc, pipeline) == "Negative") {

      if (isMode(ov, th.vectors$high_stab)) {
        return("High.Stab")
      }
      if (isMode(ov, th.vectors$x_inhibits_y)) {
        return("X.Inhibits.Y")
      }
      if (isMode(ov, th.vectors$y_inhibits_x)) {
        return("Y.Inhibits.X")
      }
      if (isMode(ov, th.vectors$neg_syn)) {
        return("Neg.Syn")
      }
      if (isMode(ov, th.vectors$emer_neg_syn)) {
        return("Emer.Neg.Syn")
      }
      else {
        return("UC")
      }
    }

    # Positive Interaction, XY > X + Y ----------------------------------------
    if(getIntClass(logfc, pipeline) == "Positive") {
      if (isMode(ov, th.vectors$low_stab)) {
        return("Low.Stab")
      }
      if (isMode(ov, th.vectors$x_restores_y)) {
        return("X.Restores.Y")
      }
      if (isMode(ov, th.vectors$y_restores_x)) {
        return("Y.Restores.X")
      }
      if (isMode(ov, th.vectors$pos_syn)) {
        return("Pos.Syn")
      }
      if (isMode(ov, th.vectors$emer_pos_syn)) {
        return("Emer.Pos.Syn")
      }
      
      # unclassified
      else {
        return("UC")
      }
    }

    else {
      return("NI")
    }
  }
  return(apply(ov_logfc, 1, modeClassifier, th.vectors, pipeline = pipeline) %>%
           factor(levels = MODES))
}

#' Evaluate model performance by initializing many analysis objects given
#' a vector of alphas for coding confidence intervals
#'
#' @param counts numeric matrix of counts where rows are genes and columns are libraries
#' @param genes data.frame of gene metadata
#' @param metadata data.frame with metadata retrieved for 4 conditions
#' @param alphas numeric vector of length n with alphas to code outcome vectors
#' with
#' @param DRUG_COMBINATION Boolean TRUE if drug combination
#' @param null Boolean TRUE if null distribution should be simulated by permuting
#' columns of count matrix
#' @return list of n = length(alphas), see output for \code{edgeRModeClassifier}
#' and \code{limmaModeClassifier}
#' @export tuneClassifier
#'
tuneClassifier <- function(counts, genes, metadata, alphas,
                           DRUG_COMBINATION = TRUE, null = FALSE) {
  res_list <- list()
  for (i in 1:length(alphas)) {
    res_list[[i]] <- edgeRModeClassifier(counts, genes, metadata,
                                         alpha = alphas[i],
                                         DRUG_COMBINATION = DRUG_COMBINATION,
                                         null = null)
  }
  return(res_list)
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


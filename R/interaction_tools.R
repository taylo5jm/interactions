#' Get character string classification for interaction class
#'
#' @param logfc Numeric vector of log fold changes or character vector indicating interaction mode
#' @param pipeline Character vector ("limma" or "edgeR")
#' @return The interaction class of \code{logfc}
#' @export getIntClass
#' @examples 
#' \dontrun{
#' logfc <- c(4, 1, 1, 1, 1, 2)
#' getIntClass(logfc, pipeline = "edgeR")
#' logfc <- "Pos.Syn"
#' getIntClass(logfc)
#' }

getIntClass <- function(logfc, pipeline = "limma") {
  # logfc is a vector where each element is a log fold change for a comparison or condition
  if (is.numeric(logfc)) {
    if (pipeline == "edgeR") {
      if (logfc[3] > (logfc[2] + logfc[1])) {
        return("Positive")
      }
      if (logfc[3] < (logfc[2] + logfc[1])) {
        return("Negative")
      }
      else {
        return("L")
      }
    }

    if (pipeline == "limma") {
      if (logfc[4] > (logfc[2] + logfc[3] - logfc[1])) {
        return("Positive")
      }
      if (logfc[4] < (logfc[2] + logfc[3] - logfc[1])) {
        return("Negative")
      }
      else {
        return("L")
      }
    }

  }
  # logfc is actually a character string indicating a mode
  if (is.character(logfc)) {
    findClass <- function(m) {
      if (m %in% c("Emer.Pos.Syn", "Low.Stab", "X.Restores.Y", "Y.Restores.X",
                 "Pos.Syn", "Pos")) {
      return("Positive")
    }

    if (m %in% c("Emer.Neg.Syn", "High.Stab", "X.Inhibits.Y", "Y.Inhibits.X",
                 "Neg.Syn", "Neg")) {
      return("Negative")
    }

    if (m %in% c("Sym.Right", "Sym.Left", "Step.Up", "Step.Down")) {
      return("No.Interaction")
    }

    else {
      return("NI")
    }
  }

  if (length(logfc) == 1) {
    findClass(logfc)
  }
  else {
    sapply(modes, findClass)
  }
}
  else {
    return(NaN)
  }
}


getInteractionData <- function(drug_combination = TRUE, local) {
  # get interaction metadata & counts -------------------------------------------
  #source('get_data.R')
  if (drug_combination == TRUE) {
    source('pipelines/functions/get_data/get_metadata_revised.R', local = local)
    print(dplyr::select(metadata, donor_id, treatment_label, treatment_name, cell_name))
    dim(counts)
    metadataCheck(metadata, cell_type)
  }

  else {
    source('pipelines/functions/get_data/get_mf_metadata.R', local = local)
    print(dplyr::select(metadata, donor_id, treatment_label, treatment_name, cell_name))
    dim(counts)
    #metadataCheck(metadata, cell_type)
  }
}

# print barplots for interactions in interaction table ------------------------
scanModeBarplots <- function(interaction_tbl, int_mode) {
  my_table <- filter(interaction_tbl, Mode == int_mode)
  for (i in 1:nrow(my_table)) {
    print(geom_bar_constructor(as.numeric(my_table[i,9:12]), int_mode, my_table[i, 2]))
    Sys.sleep(5)
  }
}


# get matrix of log fold changes from excised interaction_top_tables
getLogFC <- function(interaction_top_tables) {
  log_fcs <- lapply(interaction_top_tables, function(x) (x$logFC))
  logfc_matrix <- matrix(NaN, nrow = length(log_fcs[[1]]), ncol = 4)
  colnames(logfc_matrix) <- c("0.logFC", "X.logFC", "Y.logFC", "X+Y.logFC")
  for (i in 1:4) {
    logfc_matrix[,i] <- log_fcs[[i]]
  }
  logfc_matrix
}


# get confidence intervals from excised top tables -----------------------
# output = ci_matrix of confidence intervals, where
#           ncol(ci_matrix) = 8, in the order:
# [Signal1.LowerBound] [Signal1.UpperBound] , ..., [Signal2.UpperBound]
getLimmaCI <- function(interaction_top_tables, return.matrix = FALSE) {
  ci_l <- lapply(interaction_top_tables, function(x) (x$CI.L))
  ci_r <- lapply(interaction_top_tables, function(x) (x$CI.R))


  if (return.matrix == TRUE) {
    ci_matrix <- cbind(ci_l[[1]], ci_r[[1]], ci_l[[2]], ci_r[[2]],
                       ci_l[[3]], ci_r[[3]], ci_l[[4]], ci_r[[4]])
    colnames(ci_matrix) <- c("0.CI.L", "0.CI.R", "X.CI.L", "X.CI.R", "Y.CI.L",
                             "Y.CI.R", "XY.C+I.L", "X+Y.CI.R")
    return(ci_matrix)
  }

  if (return.matrix == FALSE) {
    return(list(mapply(cbind, ci_l), mapply(cbind, ci_r)))
  }

}


#' Generate table with genes classified into a mathematically defined interaction mode
#'
#' @param v EList constructed with \code{voom}
#' @param modes Character vector of length n indicating interaction modes of genes
#' @param log.fcs Numeric matrix of log fold changes generated from the limma or edgeR pipelines
#' @param classes Character vector of interaction classes for n genes. If missing, getIntClass()
#' is called to assign interaction classes to each interaction mode.
#' @param exclude.null.vector Boolean TRUE when genes with the null outcome vector
#' (0, 0, 0, 0, 0, 0) should be excluded from the table.
#' @param exclusivity.vectors Boolean FALSE when exclusivity vectors should be excluded
#' from the table.
#' @return \code{data.frame} with outcome vectors, interaction classes & modes, fold changes
#' for all genes classified into an interaction class
#' @export interactionGeneTable


interactionGeneTable <- function(v, modes, log.fcs, pw_s, classes,
                                 exclude.null.vector = TRUE,
                                 exclusivity.vectors = FALSE,
                                 DRUG_COMBINATION = TRUE, ...) {
  if (missing(classes)) {
    classes <- getIntClass(modes[which(
      modes != "NI" & modes != "A" & modes != "UC")])
  }

  # map from modes to outcome vectors, classess, etc.
  getInteractionIndices <- function(modes, exclude.null.vector = TRUE,
                                    exclusivity.vectors = FALSE) {
    if (exclude.null.vector == TRUE & exclusivity.vectors == TRUE) {
      return(which(modes != "NI" & modes != "A" & modes != "Step.Up" &
                     modes != "Step.Down" & modes != "UC"))
    }

    if (exclude.null.vector == TRUE & exclusivity.vectors == FALSE) {
      return(which(modes != "NI" & modes != "A" & modes != "UC"))
    }
  }

  interactions <- getInteractionIndices(modes,
    exclude.null.vector = exclude.null.vector,
                                        exclusivity.vectors = exclusivity.vectors)

  if (nrow(log.fcs) != length(interactions)) {
    log.fcs <- log.fcs[interactions,]
  }

  # initialize data.frame
  int_genes <- data.frame(v$genes[interactions,], classes, modes[interactions],
                          pw_s[interactions], log.fcs)
  # return empty df if data.frame does not fit dimensions for
  # edgeR and limma pipelines
  if (ncol(int_genes) != 14 & ncol(int_genes) != 12) {
    return(data.frame())
  }
  # 4 fold changes --------------------------------------------------------
  # make data.frame for limma pipeline results with 4 coef
  if (ncol(int_genes) == 12) {
    # make data.frame if signal combination satisfies
    # T0 = Condition 1, Vehicle
    # T1 = Condition 2, Drug X
    # T2 = Condition 3, Drug Y
    # T12 = Condition 4, Drug X + Y
    # T0 <- "None"; T1 <- "Sitagliptin". T2 <- "Metformin";
    # T12 <- "Sitagliptin / Metformin"
    if (DRUG_COMBINATION == TRUE) {
    colnames(int_genes) <- c(colnames(v$genes), "Class", "Mode", "Outcome.Vector",
                             "Control.logFC", paste0(T1, ".logFC"),
                             paste0(T2, ".logFC"), paste0(T12, ".logFC"))
    }
    # make data.frame if signal combination satisfies
    # T0 = Condition 1, Female.Vehicle
    # T1 = Condition 2, Male.Vehicle
    # T2 = Condition 3, Female.SignalX
    # T12 = Condition 4, Male.SignalY
    # T1 <- "Sitagliptin"
    if (DRUG_COMBINATION == FALSE) {
      colnames(int_genes) <- c(colnames(v$genes), "Class", "Mode", "Outcome.Vector",
                               "None:F.logFC", "None:M.logFC",
                               paste(T1, "F.logFC", sep = ":"),
                               paste(T1, "M.logFC", sep = ":"))
    }
  }

  # 6 fold changes from contrasts ---------------------------------------------
  if (ncol(int_genes) == 14) {
    if (DRUG_COMBINATION == FALSE) {
      T12 <- paste("Male", T1, sep = ":")
      T2 <- paste("Female", T1, sep = ":")
      T0 <- "Female:None"
      T1 <- "Male:None"
    }

      colnames(int_genes) <- c(colnames(v$genes), "Class", "Mode", "Outcome.Vector",
                               paste0(T1, "v", T0, ".logFC"),
                               paste0(T2, "v", T0, ".logFC"),
                               paste0(T12, "v", T0, ".logFC"),
                               paste0(T2, "v", T1, ".logFC"),
                               paste0(T12, "v", T1, ".logFC"),
                               paste0(T12, "v", T2, ".logFC"))

  }
  arrange(int_genes, Class, Mode)
}

#' Create data.frame for genes that were classified into modes by both limma and edgeR
#'
#' @param edgeR_interaction_tbl interaction table (data.frame) generated by edgeR
#' interaction mode classification
#' @param limma_interaction_tbl interaction table (data.frame) generated by limma
#' interaction mode classification
#' @return interaction table (data.frame) with genes classified into non-null modes
#' by both limma and edgeR
#' @export sharedGenesTable
#' @examples
#' library(limma)
#' library(edgeR)
#' \dontrun{
#' master_interaction_tbl <- sharedGenesTable(edgeR_interaction_tbl,
#'                                          limma_interaction_tbl)
#' }


sharedGenesTable <- function(edgeR_interaction_tbl, limma_interaction_tbl) {
  e_genes <- edgeR_interaction_tbl$hg19_gene_id
  l_genes <- limma_interaction_tbl$hg19_gene_id
  prop_shared_e <- length(which(e_genes %in% l_genes)) / length(e_genes)
  prop_shared_l <- length(which(l_genes %in% e_genes)) / length(l_genes)
  shared_gene_indices <- which(e_genes %in% l_genes)
  master_interaction_tbl <- edgeR_interaction_tbl[shared_gene_indices,]
  return(master_interaction_tbl)
}

#' Report number of genes classified into a non-null or non-anomalous interaction mode
#'
#' @param modes character vector of length n indicating the interaction modes classified
#' by a particular pipeline
#' @return number of genes classified into a non-null and non-anomalous interaction mode
#' @examples
#' \dontrun{
#' reportnInteractions(limma_modes)
#' }
#'
# Print statistics to console -------------------------------------------------
reportnInteractions <- function(modes) {
  print(paste0(length(which(modes != "A" &
                              modes != "NI")),
               " interactions were classified with confidence intervals
generated with limma."))
}

#' Report number of genes classified into an anomalous mode classifcation
#'
#' @param modes character vector of length n indicating the interaction mode classifications
#' of each gene in v$genes
#' @return number of genes coded into anomalous outcome vectors
#' @examples
#' \dontrun{
#' reportNAnomaly(edgeR_modes)
#' }
reportnAnomaly <- function(modes) {
  print(paste0(length(which(modes == "A")), " genes exhibited
               anomalous outcome vectors. These outcome vectors may have been able
               to be classified into an interaction class, but the outcome
               vectors are not included in the 75 unique vectors enumerated."))
}

#' Report number of null interaction vectors coded
#'
#' @param character vector of length n indicating the interaction mode classifications
#' of each gene in v$genes
#' @return number of null outcome vectors (0, 0, 0, 0, 0, 0) coded
#' @examples
#' \dontrun{
#' reportnNI(limma_modes)
#' }
reportnNI <- function(modes) {
  print(paste0(length(which(modes == "NI")), " genes were consistent
               with a non-interaction. Genes are classified as non-interactions
               if their class inequality confidence interval contains 0, or if
               the pairwise comparison outcome vector is known to exhibit no
               interactive effects"))
}

#' Print statistics about the interaction mode distribution and return a frequency distribution
#' table of interaction mode classifications
#'
#' @param modes character vector of length n indicating the interaction modes for each gene
#' @return table of
reportModeClassStats <- function(modes) {
  reportnInteractions(modes)
  reportnAnomaly(modes)
  reportnNI(modes)
  print(table(modes))
  return(table(modes))
}

### Interaction Mode Classification -------------------------------------------

#' Check to see if an outcome vector is consistent with a "non-interaction" mode
#'
#' @param v character vector of length n = 6 or numeric vector of length m = 1 indicating the outcome vector of interest
#' @return Boolean TRUE when vector is consistent with 1/5 "non-interaction" vector
#' @export niCheck
#' @examples
#' null_vector <- c(0, 0, 0, 0, 0, 0)
#' niCheck(null_vector)
#' high_stab_vector <- c(1, 1, 1, 0, 0, 0)
#' niCheck(high_stab_vector)
niCheck <- function(v) {

  if (is.numeric(v)) {
    v <- v[-7]
    if (all(v == rep(0, 6)) | all(v == c(1, 0, 1, -1, 0, 1)) |
          all(v == c(0, 1, 1, 1, 1, 0)) | all(v == c(-1, 0, -1, 1, 0, -1)) |
          all(v == c(0, -1, -1, -1, -1, 0))) {
      return(TRUE)
    }
    else {
      return(FALSE)
    }
  }

  if (is.character(v)) {
    if (!(v %in% all_vecs)) {
      return(TRUE)
    }
    else {
      return(FALSE)
    }
  }

}

#' apply expected counts filter to DGEList object
#' Aaron Mackey
#' @param dge DGEList
#' @param cpm_min minimum counts per million
#' @param sample_min minimum number of samples
#' @return filtered DGEList
#' @export massageCounts

massageCounts <- function(dge, cpm_min=2, sample_min=5) {
  # filter em (this can change in the future)
  # current method is using a cpm of greater than 2 in at least 10 samples
  keep <- rowSums(cpm(dge) > cpm_min) >= round(sample_min)
  dge <- dge[keep,]
  dge$samples$lib.size <- colSums(dge$counts)
  return(dge)
}

#' generate frequency distribution table for outcome vectors
#'
#' @param pw_s Character vector of length n with pairwise comparison outcome vectors
#' if pw_s is a matrix, then each row will be coerced to a character vector of
#' length 1
#' @return A data.frame with the frequency distribution of outcome vectors and
#' associated interaction classes
#' @export ovTable

ovTable <- function(pw_s, get.class = FALSE) {
  if (is.matrix(pw_s)) {
    pw_s <- apply(pw_s, 1, stringCoerce)
  }
  if (is.character(pw_s)) {
    pw_s <- pw_s %>% table %>% as.data.frame
    colnames(pw_s) <- c("Outcome.Vector", "Freq")
    pw_s$Log.Freq <- log(pw_s$Freq)
    pw_s$Prop <- pw_s$Freq / sum(pw_s$Freq)
    if (get.class == TRUE) {
       pw_s$Class <- lapply(pw_s$Outcome.Vector, classFromOV) %>% unlist
    }
    return(pw_s)
  }
}


summarizeOVTable <- function(ov.table, var) {
  vars <- unique(ov.table$var)
  ov_sum_table <- data.frame()
  for (i in 1:length(vars)) {
    this_tbl <- filter(ov.table, Class == var)
  }
}




#' Generate interaction table with genes classified into mathematically defined modes by both limma and edgeR
#'
#' @param edgeR_interaction_tbl interaction table with data generated with edgeR
#' @param limma_interaction_tbl interaction table with data generated with limma
#' @return master_interaction_tbl, interaction table with genes classified into modes by limma and edgeR
#' @export poolInteractionGenes

poolInteractionGenes <- function(edgeR_interaction_tbl, limma_interaction_tbl) {

  e_genes <- edgeR_interaction_tbl$hg19_gene_id
  l_genes <- limma_interaction_tbl$hg19_gene_id
  prop_shared_e <- length(which(e_genes %in% l_genes)) / length(e_genes)
  prop_shared_l <- length(which(l_genes %in% e_genes)) / length(l_genes)
  shared_gene_indices <- which(e_genes %in% l_genes)
  master_interaction_tbl <- edgeR_interaction_tbl[shared_gene_indices,]
  return(master_interaction_tbl)
}

#' Generate table with outcome vectors, interaction classes & modes for all genes
#'
#' @param outcome.vectors matrix where nrow(matrix) == n or character vector of length = n indicating outcome vectors
#' @param modes character vector of length n indicating interaction modes
#' @param classes character vector of length n indicating interaction classes
#' @return data.frame with outcome vectors, interaction modes, and classes for each gene
#' @export interactionDistTable
interactionDistTable <- function(outcome.vectors, modes, classes) {
  if (is.matrix(outcome.vectors)) {
    outcome.vectors <- apply(outcome.vectors, 1, stringCoerce)
  }
  this_df <- data.frame(outcome.vectors, modes, classes)
  colnames(this_df) <- c("Outcome.Vector", "Mode", "Class")
  return(this_df)
}

summarizeDistTable <- function(interaction.dist.table) {
  return(count(interaction.dist.table))
}

#' Generate character vector of all possible interaction mode classifications
#' implemented in classifyByThOutcomeVector()
#'
#' @return MODES character vector of interaction mode classifications. See
#' the vignettes on the implementation of mode classification.
#' @export allModes
#' @examples
#' MODES <- allModes()

allModes <- function() {
  MODES <- c("Low.Stab", "X.Restores.Y", "Y.Restores.X", "Pos.Syn", "Emer.Pos.Syn",
              "High.Stab", "X.Inhibits.Y", "Y.Inhibits.X", "Neg.Syn", "Emer.Neg.Syn",
              "Sym.Right", "Sym.Left", "Step.Up", "Step.Down", "NI", "A", "UC")
  return(MODES)
}


#' @export metadataTable
metadataTable <- function(counts, genes, metadata, alpha, drug.combo) {
  return(list(
    counts = counts,
    genes = genes,
    metadata = metadata,
    alpha = alpha,
    drug.combo = drug.combo))
  }










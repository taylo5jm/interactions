#' Principal Component Analysis scores using log-fold changes as feature vectors
#'
#' @param pooled_lfc matrix of log fold change ratios for pairwise comparisons.
#' \code{pooled_lfc} is intended to be a matrix of pooled log fold change 
#'  ratios that map to known interaction modes
#' @return matrix of PCA scores
#' @examples
#' \dontrun{
#' getPCAScores(pooled_lfc)
#' }

getPCAScores <- function(pooled_lfc) {
  my_pca <- stats::princomp(pooled_lfc)
  pca_scores <- my_pca$scores
  return(pca_scores)
}

#' Partition feature vectors (pca_scores) and interaction mode classifications
#'  into training and test sets for classification with KNN
#'
#' @importFrom caret createDataPartition
#' @param pooled_modes character vector representing interaction mode 
#'  classifications pooled from n signal combinations
#' @param pca_scores matrix of PCA scores generated from \code{stats::princomp}
#' @return list with training scores, training modes, testing scores, and 
#'  testing modes
#' @examples
#' \dontrun{
#' #FINISH LATER
#' }

partitionData <- function(pooled_modes, pca_scores) {

  main_indices <- which(pooled_modes != "A" & pooled_modes != "NI")
  my_train_modes <- factor(pooled_modes[main_indices])
  my_train_scores <- pca_scores[main_indices,]
  test_scores <- pca_scores[which(pooled_modes == 16),]
  test_modes <- pooled_modes[which(pooled_modes) == "A" | 
                               which(pooled_modes) == 16]
  # split data with known modes into training and test sets -------------------
  my_split <- createDataPartition(my_train_modes, times = 1)
  my_scores <- my_train_scores[my_split[[1]],]
  my_modes <- factor(my_train_modes[my_split[[1]]])
  
  return(list(
    training.scores = my_scores,
    training.modes = my_modes,
    testing.scores = test_scores,
    testing.modes = test_modes))
}

#' Train k-nearest neighbours model with 10-fold cross-validation using PCA 
#' scores as feature vectors and interaction modes as classifications
#'
#' @importFrom caret trainControl train
#' @param train_scores matrix where each row represents a feature vector (gene)
#' @param train_modes character vector of length n where each element indicates 
#' the interaction mode classification for a gene
#' @return caret k-nearest neighbours model
#' @examples
#' \dontrun{
#' knn_model <- trainKNN(my_pca_scores, my_modes)
#' }
#' 
trainKNN <- function(train_scores, train_modes) {
  # 10-fold cross-validation
  train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
  # train k-nearest neighbours model
  model <- train(train_scores, train_modes,
                  trControl = train_control,
                  method = "knn")
  return(model)
}

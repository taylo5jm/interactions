Anomaly Classification with k-nearest neighbours
========================================================

Some genes may be coded as an outcome vector that was not enumerated as one of the 75 theoretical outcome vectors. Recall there are 256 possible permutations with repetition when there are n = 4 conditions and k = (1, 2, 3, 4) possible ranks. We would expect to see only a subset of these n-tuples when classifying a matrix of counts, since some n-tuples would not be meaningful in terms of biology in our current framework.

This subset, however, can often have a greater number of outcome vectors than the number of outcome vectors that were enumerated with the closed-form recurrence relation solution and [maybe another example]. Outcome vectors that were coded but do not match an outcome vector in the 75 enumerated, as anomalous outcome vectors.

The authors of "combinatorial code governing cellular response to complex stimuli" coded outcome vectors with confidence intervals (see [LINK TO LIMMA PAGE]). If the profile could not be matched to a theoretically enumerated outcome vector due to "statistical inconsistencies", the nearest-neighbour algorithm with correlation distance was used to classify these anomalies. It is later mentioned, however, that the MATLAB function *knnclassify* was used to classify the profile into an interaction mode with k-nearest neighbors.

We will use the knn classifier function in caret to train a model, so that some of these anomalous outcome vectors can be classified into an interaction mode.


## Principal Component Analysis

We will use principal component analysis (*stats::princomp*) to convert our variables of interest to principal components. In this application, our possibly correlated variables are log fold changes for contrasts. `pooled_lfc` is a matrix where each row represents a gene and each column represents a condition or contrast. This variable is named `pooled_lfc`, because ideally we would have a large matrix of fold-changes generated from *multiple* signal combinations, not just our signal combination of interest.


```r
getPCAScores <- function(pooled_lfc) {
  my_pca <- princomp(pooled_lfc)
  pca_scores <- my_pca$scores
  return(pca_scores)
}

# ... run analysis
this_pca_scores <- getPCAScores(pooled_lfc)
```

### Plot Principal Components
Now that we have our PCA scores, we should verify that the PCA scores associated with each interaction are clustered by at least one of the components. If we were to graph all of the possible combinations of principal componenents, then a given interaction mode should be observed in an identifiable cluster in at least one of these graphs.

`pcaScatterMatrix` uses *GGally::ggpairs* to render a scatter plot matrix for this purpose. We will first plot the PCA scores of all non-additive interactions (interaction modes belonging to a positive or negative class). The null-vector is excluded, as it is a super majority class. The presence of this class in a knn model would cause many outcome vectors to be mistakenly classified as a null-interaction.



```r
library(GGally)
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

pcaScatterMatrix(this_modes, this_pca_scores)
```

## Interaction Mode Classification with caret
We see that most (if not all) [CHECK] of the interaction modes are separated from the rest by at least one principal componenent. Now, a model can be trained so that predictions can be made to classify anomalous pairwise comparison outcome vectors into defined interaction modes.

### Partition Data into Training and Testing Sets
The data is partitioned into two subsets, a training and test set, both of which include PCA scores and interaction modes. Our training set will be used to teach the model how to classify feature vectors into known interaction modes. Our test set will be used to evaluate the models performance, since we know the true classifications for the test set.


```r
library(caret)
partitionData <- function(pooled_modes, pca_scores) {

  main_indices <- which(pooled_modes != "A" & pooled_modes != "NI")
  my_train_modes <- factor(pooled_modes[main_indices])
  my_train_scores <- pca_scores[main_indices,]

  test_scores <- pca_scores[which(pooled_modes == 16),]
  test_modes <- pooled_modes[which(pooled_modes) == "A" | which(pooled_modes) == 16]
  # split data with known modes into training and test sets ---------------------

  my_split <- createDataPartition(my_train_modes, times = 1)
  my_scores <- my_train_scores[my_split[[1]],]
  my_modes <- factor(my_train_modes[my_split[[1]]])
  return(list(
    training.scores = my_scores,
    training.modes = my_modes,
    testing.scores = test_scores,
    testing.modes = test_modes))
}
```

### Train k-nearest neighbours model

*caret* will be used to train a knn model with our training data set, comprised of PCA scores and true classifications for the PCA scores. 10-fold repeated cross-validation will be used so that we are able to validate the predictive accuracy of the model.


```r
trainKNN <- function(my_train_scores, my_train_modes) {
  # 10-fold cross-validation
  train_control <- trainControl(method="repeatedcv", number=10, repeats=3)

  # train k-nearest neighbours model
  model <- train(my_train_scores, my_train_modes,
                  trControl = train_control,
                  method = "knn")
  return(model)
}

this_model <- trainKNN(training_scores, training_modes)
this_model
```

### Make predictions
Now that we have trained our knn model, we can make predictions using our testing PCA scores as feature vectors.

```r
predicted_modes <- predict(this_model, newdata = testing_scores)
table %>% predicted_modes
```

### Evaluate model performance
Finally, we can retrieve some statistics (such as kappa) to evaluate the performance of our *k*-nearest neighbours model.

```r
confusionMatrix(predicted_modes, test_modes)
```





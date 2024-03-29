Sample Workflows: Interaction Mode Classification with limma and edgeR
====================================================================
This report describes the steps involved in classifying genes into a mathematically and biologically defined interaction mode and class using the `limma` and `edgeR` pipelines. Please refer to function documentation for more information about a function defined or called in this report, as the documentation may contain details about a function that are omitted here for conceptual clarity.

## Retrieve Metadata
In order to classify gene expression profiles into interaction modes and classes, we will require 3 primary objects: `counts`, `genes`, and `metadata`.

1. `metadata` is a `data.frame` with information about the libraries retrieved for `counts`. This `data.frame` should include variables related to the libraries of interest.  batch(es), treatments, age, sex, etc.

```r
CON <- dbConnect(MySQL(), user = "eld", password = '5thStreet',
                 dbname = 'eld_v4')

# metadata --------------------------------------------------------------------
combo_exists <- sitt::isValidCombination(T1, T2,
                                   DRUG_COMBINATION = TRUE)
metadata <- getInteraction(CON, T1, T2,
                           T12, "EC") %>%
            sitt::filterMetadata(conc = "Mid", sex = "M", cell.type = "EC")
metadata <- getVehicles(CON, "EC", metadata)
```

2. `counts` matrix of expected counts where each row is a gene and each column is a library. We refer to these counts as *expected* in this report, since the counts retrieved were obtained by random sampling of the population of all transcripts present in the cell during sample collected, thus the counts obtained do not represent the absolute total number of transcripts present in the cell at a given time, *t*.


```r
getExpectedCounts <- function(metadata, level = "gene") {
 getCounts <- function(rnaseq_assay_id, level) {
   dbGetQuery(CON, paste0('SELECT * from ', paste0(level, '_rsem_count '),
                          'WHERE rnaseq_assay_id = ', rnaseq_assay_id))
 }
 counts <- sapply(metadata$rnaseq_assay_id, getCounts, level) %>%
   as.matrix()
 return(counts)
}

# strip extra variables from counts
asCountsMatrix <- function(counts) {
 counts <- counts[4,1:ncol(counts)] %>% sapply(unlist)
}
```

3. `genes` refers to a `data.frame` with gene metadata retrieved from the *hg19_gene* table in *eld_v4*.


```r
getGeneMetadata <- function(CON, counts, level = "gene") {
 gene_id <- counts[,2] %>% unlist()
 dbGetQuery(CON, paste0('SELECT * from hg19_', level)) %>%
   filter(hg19_gene_id %in% gene_id)
  }
```


```
## Error in match.fun(FUN): object 'dbDisconnect' not found
```

```
## Loading required package: DBI
## 
## Attaching package: 'dplyr'
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
## 
## Loading required package: limma
## 
## Attaching package: 'sitt'
## 
## The following objects are masked _by_ '.GlobalEnv':
## 
##     asCountsMatrix, getExpectedCounts, getGeneMetadata
```

```
## Error in getExpectedCounts(CON, metadata, level = "gene"): unused argument (metadata)
```

```
## Error in eval(expr, envir, enclos): object 'counts' not found
```

```
## Error in eval(expr, envir, enclos): object 'counts' not found
```

## Initalize DGEList, massage counts and calculate normalization factors
Suppose we have a matrix where rows are genes and columns are libraries. We also have a `data.frame` of gene metadata (`genes`). We will use the `DGEList` class from `edgeR` to hold our matrix of counts and gene metadata. Rows with all 0 counts are removed. We then use `massageCounts` to eliminate genes from our analysis below a threshold cpm. Then, we calculate normalization factors

This object is also returned to the global environment so that it can be used in the `edgeR` pipeline. Finally, we use `voom` from limma to estimate the mean-variance relationship and compute observational-level weights. After this function executes, `d`, a `DGEList` with normalization factors calculated and `v`, an `EList` object, can be input to the `edgeR` and `limma` pipelines respectively.


```r
dge2voom_q <- function(counts, genes) {
  d <<- DGEList(counts,
               genes = genes, remove.zeros = TRUE) %>%
       massageCounts(cpm_min=2, sample_min=5) %>%
       calcNormFactors
  voom(d)
}

v <- dge2voom_q(counts, genes)
```

## Design matrix and contrasts
Construct design matrix and contrasts for a given pipeline. We pass `designMatrixAndContrasts` a metadata data.frame, a character vector indicating which pipeline ("limma" or "edgeR we are using")

```r
designMatrixAndContrasts <- function(metadata, pipeline,
                                     DRUG_COMBINATION = TRUE) {

  if (pipeline == "edgeR") {
    if (DRUG_COMBINATION == TRUE) {
      metadata$treatment_name <- factor(metadata$treatment_name)
      metadata$donor_id <- factor(metadata$donor_id)

      colnames(d$counts) <- paste0(metadata$treatment_name, metadata$donor_id)

      this_design <- model.matrix(~0 + treatment_name + donor_id, data = metadata)
      print(this_design)
      #colnames(this_design) <- c("X", "C", "Y", "XY", factor(donors[-1]))
      this_design <- designColnames(this_design)
      nrows <- which(metadata$treatment_name == T1) %>% length
      this_contrasts <- makeContrasts(X-C, Y-C, XY-C, Y-X, XY-X, XY-Y,
                                      levels = c("C", "X", "Y", "XY")) %>%
        rbind(matrix(0, nrow = length(unique(metadata$donor_id)) - 1, ncol = 6))
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
```

# edgeR

To begin interaction mode classification with the `edgeR` pipeline, we will need `d`, a DGEList object, `design` a list where the first element is a design matrix and the second element are contrasts as well as `metadata`, a data.frame with sample metadata. With `edgeR`, we can use the negative binomial model for differential expression analysis.

First, we estimate the common, trended and tagwise dispersion of our libraries and their respective gene expression outputs. Next, we fit a generalized linear model and make contrasts. Our six contrasts of interest will directly map to the six elements in the pairwise comparison outcome vectors.

With our model and contrasts, we can run likelihood ratio tests to estimate the model parameter for each contrast. With the results from our likelihood ratio tests, we can use the p-values and fold-change ratios to code our outcome vectors. With outcome vectors and fold-change ratios, genes can be classified into an interaction class and interaction mode.

### Estimate common, trended and tagwise dispersion

Recall that `d` has already had its normalization factors calculated. Now, we proceed to estimate the common, trended and tagwise dispersions of `d` as a function of `design`. `estimateAllGLMDisp` below returns `d` with all dispersions estimated.


```r
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
```


### Fit generalized linear model and run likelihood ratio tests for each contrast
Now that we have estimated the dispersions for our DGEList object, a generalized linear model can be fit and likelihood ratio tests can be run for each of the 6 contrasts. `glmFitAndLRT` returns a list where each element holds the results from a likelihood ratio test.


```r
glmFitAndLRT <- function(d, design, contrasts) {

  if (is.list(design) & length(design) == 2) {
    this_design <- design[[1]]
    contrasts <- design[[2]]
  }

  fit <- glmFit(d, this_design)
  pairwise <- list()
  for (i in 1:6) {
    pairwise[[i]] <- glmLRT(fit, contrast = contrasts[,i])
  }
  return(pairwise)
}
```

### Code pairwise comparison outcome vectors for each gene

We now have log-fold-change ratios and p-values generated from likelihood ratio tests (`edgeR::glmLRT`). With these two sets of values, elements of the pairwise comparison outcome vectors can be coded with the following algorithm.

1. if p-value for likelihood ratio test < alpha
   outcome vector element = sign(log-fold-change ratio)
2. if p-value for likelihood ratio test > alpha
   outcome vector element = 0

Recall that log-fold-change-ratios for each contrast are needed to classify each gene into a mode, along with the pairwise comparison outcome vectors. We construct a matrix where the first 6 columns represent the pairwise comparison outcome vectors and columns 7:12 represent the log fold change ratios for each contrast.

```r
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
  my_res <- lapply(pairwise, codeOutcomeVector, p_val)
  edgeR_outcome_vectors <- mapply(cbind, my_res)

  # mode classification ---------------------------------------------------------
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
```

### Sample Workflow: Interaction Mode Classification with edgeR

The following is a sample workflow for interaction mode classification with edgeR. The code below assumes that counts, genes and metadata have been retrieved

```r
# get counts, genes and metadata
# ...

v <- dge2voom_q(counts, genes) #init DGEList
this_design <- designMatrixAndContrasts(metadata, DRUG_COMBINATION, pipeline = "edgeR")
# estimate dispersion ---------------------------------------------------------
d <- estimateAllGLMDisp(d, this_design) #estimate all glm dispersions
signif(sqrt(d$common.dispersion))       # Biological coeff. of variation

pairwise <- glmFitAndLRT(d, this_design, contrasts)
# outcome vectors and fold changes
edgeR_outcome_vectors <- codeEdgeROutcomeVectors(pairwise, p_val)
edgeR_pw_s <- apply(edgeR_outcome_vectors, 1, stringCoerce)
edgeR_modes <- modeClassification(edgeR_outcome_vectors, pipeline = "edgeR")
# interaction modes
edgeR_interaction_tbl <- interactionGeneTable(edgeR_modes,
                                              edgeR_outcome_vectors[,7:12],
                                              edgeR_pw_s) # summary table
edgeR_interaction_tbl
```


# limma

We will assume that a `counts` matrix and `genes` data.frame exists to initalize an EList object with `voom`.

```r
v <- dge2voom_q(counts, genes)
```

We then construct our design matrix.

```r
design <- designMatrixAndContrasts(arb.args, pipeline == "limma")
```

With out EList object and design matrix, we can fit a linear model and use empirical Bayes techniques to squeeze the standard errors. `eBayes` will output t-statistics, moderated F-statistics and log-odds of differential expression

```r
fit <- lmFit(v, design)
eb <- eBayes(fit)
```

A table of differentially expressed genes can then be produced for each coefficient associated with a treatment in the linear model. Note that the parameters of the function below are set so that all genes appear unsorted in the data.frame returned by topTable, since we will need confidence intervals for each gene to code the pairwise comparison outcome vectors.

When passing `ci` to the `confint` parameter in `limma::topTable`, we are telling limma to make confidence intervals at `ci` confidence. Finally, `adjust.method` will adjust the p-values with false discovery rate to adjust for multiple testing


```r
# differentially expressed genes ----------------------------------------------
ci <- 0.99
top_tables <- list()
for (i in 1:ncol(design)) {
  top_tables[[i]] <- topTable(eb, coef = i, number = nrow(v$E), confint = ci,
                              p.value = 1, sort.by = "none", adjust.method = "fdr")
}
```

When constructing our design matrix, the coefficients may not be in the same column as we would expect them to be. This occurrence is not necessarily undesirable, but it is convenient to order the coefficients in a way that the results from one signal combination can be compared to another signal combination.

For a drug combination, the conditions will be ordered as such:
  1. Vehicle
  2. Signal X             (ex. Sitagliptin)
  3. Signal Y             (ex. Metformin)
  4. Signal X + Y         (ex. Sitagliptin / Metformin)

  note: current code will assume that all libraries were retrieved from one sex

For comparing females and males treated with one signal (e.g. drug):
  1. Vehicle:Female       (ex. None:Female)
  2. Vehicle:Male         (ex. None:Male)
  3. Signal X:Female      (ex. Signal X:Female)
  4. Signal X + Y:Male    (ex. Signal X:Male)

We subset a topTable for each condition so that we have means and confidence intervals for each condition.


```r
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
```

### Code outcome vector elements

Outcome vectors with limma are coded using the confidence intervals generated from topTable. The steps in the algorithm are listed below:

1. If the logFC of condition 1 is greater than the lower bound and less than the upper bound of the confidence intervals for condition 2 (i.e. if sample mean logFC of condition 1 is in the confidence interval of condition 2)
  assign 0 to element of outcome vector
2. If the logFC of condition 1 is less than the lower bound CI of condition 2 and if the logFC of condition 2 is greater than the upper bound of condition 2 (i.e. gene expression output of condition 1 < gene expression output of condition 2)
  assign -1 to element of outcome vector
3. Else (i.e. gene expression output of condition 1 > gene expression output of condition 2)
  assign 1 to the element of outcome vector

```r
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
```

VERY VERY WEAK POINT
Now that we have our pairwise comparison outcome vectors, we create a matrix where the first 6 columns represent the outcome vectors and columns 7:10 represent the fold changes for each condition

```r
fold_changes <- cbind(top_tables[[2]]$logFC, top_tables[[3]]$logFC,
                        top_tables[[1]]$logFC, top_tables[[4]]$logFC)

  limma_outcome_vectors <- cbind(validateOutcomes(limma_outcome_vectors), fold_changes)
```


#### Interaction Mode Classification & Summary
WEAK POINT HACKY CODE
We can then classify each gene into an interaction class and mode. Finally, we create an a summary table with all genes that were classified into a non-null outcome vector.

```r
limma_modes <- apply(limma_outcome_vectors, 1, classifyByThOutcomeVector, "limma") %>%
    unlist

  interaction_top_tables2 <- exciseTopTables2(top_tables, v, limma_modes)
  log_fold_changes2 <- getLogFC(interaction_top_tables2)

  # summary table ---------------------------------------------------------------
  limma_interaction_tbl <- interactionGeneTable(limma_modes, log_fold_changes2,
                                                limma_pw_s)
```


#### Example Workflow: Interaction Mode Classification with Limma

```r
limmaModeClassification <- function(v, metadata) {

  # make design matrix, fit linear model and empirical bayes --------------------
  if (DRUG_COMBINATION == TRUE) {
    design <- model.matrix(~0 + treatment_name + donor_id, data = metadata)
    colnames(design) <- c(levels(factor(metadata$treatment_name)), levels(factor(donors[-1])))
  }

  if (DRUG_COMBINATION == FALSE) {
    design <- model.matrix(~0 + treatment_name:donor_sex, data = metadata)
    colnames(design) <- sub("treatment_name", "", colnames(design))
    colnames(design) <- sub("donor_sex", "", colnames(design))
  }

  fit <- lmFit(v, design)
  eb <- eBayes(fit)

  # differentially expressed genes ----------------------------------------------
  top_tables <- list()

  for (i in 1:ncol(design)) {
    top_tables[[i]] <- topTable(eb, coef = i, number = nrow(v$E), confint = ci,
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

  print("Classifying modes with confidence intervals generated by limma:topTable(confint = TRUE)")

  # Interaction Mode Classification ---------------------------------------------
  # logFC from limma::topTable (control, signal x, signal y, signal xy)
  fold_changes <- cbind(top_tables[[2]]$logFC, top_tables[[3]]$logFC,
                        top_tables[[1]]$logFC, top_tables[[4]]$logFC)

  limma_outcome_vectors <- cbind(validateOutcomes(limma_outcome_vectors), fold_changes)
  limma_modes <- apply(limma_outcome_vectors, 1, classifyByThOutcomeVector, "limma") %>%
    unlist

  interaction_top_tables2 <- exciseTopTables2(top_tables, v, limma_modes)
  log_fold_changes2 <- getLogFC(interaction_top_tables2)

  # summary table ---------------------------------------------------------------
  limma_interaction_tbl <- interactionGeneTable(limma_modes, log_fold_changes2,
                                                limma_pw_s)
  return(list(interaction.tbl = limma_interaction_tbl,
              modes = limma_modes,
              limma_pw_s))
}
```

### Example Workflow: Classification with limma and edgeR

#### Interaction Mode Classification

We are primarily interested in interaction effects between Sitagliptin and Metformin. Our four conditions will therefore be:

1. Vehicle
2. Sitagliptin
3. Metformin
4. Sitagliptin / Metformin

We explicitly state with `DRUG_COMBINATION = TRUE` that condition 4 is a drug combination of signal x and signal y.

First, we will retrieve some libraries and metadata we are interested in. Next, we will code our outcome vectors with confidence intervals or likelihood ratio tests with limma or edgeR respectively. Then we classify interaction classes, modes, and generate a summary table.

`limmaModeClassifier` and `edgeRModeClassifier` output a list with outcome vectors, interaction modes, summary tables and other objects output from limma and edgeR.



```r
library(magrittr)
library(RMySQL)
library(dplyr)
library(stringr)
library(edgeR)
library(sitt)

T0 <- "None"
T1 <- "Sitagliptin"
T2 <- "Metformin"
T12 <- "Sitagliptin / Metformin"

# enumerate theoretical outcome vectors
th_vectors <- outcomeVectorsByMode()

# connect to database
CON <- dbConnect(MySQL(), user = "eld", password = '5thStreet',
                 dbname = 'eld_v4')

# metadata --------------------------------------------------------------------
combo_exists <- sitt::isValidCombination(T1, T2,
                                   DRUG_COMBINATION = TRUE)
metadata <- getInteraction(CON, T1, T2,
                           T12, "EC") %>%
            sitt::filterMetadata(conc = "Mid", sex = "M", cell.type = "EC")
metadata <- getVehicles(CON, "EC", metadata)

counts <- getExpectedCounts(CON, metadata, level = "gene")
genes <- getGeneMetadata(CON, counts, level = "gene")
counts <- asCountsMatrix(counts)

# interaction mode classification with limma -------------------------------
sit_met_limma <- limmaModeClassifier(counts, genes, metadata, DRUG_COMBINATION = TRUE,
                            ci = 0.99)
null_limma <- limmaModeClassifier(counts, genes, metadata, DRUG_COMBINATION = TRUE,
                            ci = 0.99, null = TRUE)

sit_met_edgeR <- edgeRModeClassifier(counts, genes, metadata,
                             DRUG_COMBINATION = TRUE,
                             alpha = 0.01)
null_edgeR <- edgeRModeClassifier(counts, genes, metadata, null = TRUE,
                                  DRUG_COMBINATION = TRUE,
                                  alpha = 0.01)
```



#### Gene Ontology
These objects are actually saved in the /data folder of this package, so we can continue with the analysis by using `base::data()` to load the R objects into RAM.

```r
data(sit_met_edgeR)
```

```
## Warning in data(sit_met_edgeR): data set 'sit_met_edgeR' not found
```

```r
data(null_edgeR)
```

```
## Warning in data(null_edgeR): data set 'null_edgeR' not found
```




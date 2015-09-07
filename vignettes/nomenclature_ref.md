Outcome Vectors, Interaction Modes & Classes Reference
========================================================

```r
library(magrittr)
```

### Enumeration of Outcome Vectors

For 4 gene expression outputs ($e_{0}$, $e_{X}$, $e_{Y}$, $e_{X + Y}$), there are 75 possible outcome vectors if each element of the outcome vector can assume the value of -1, 0, or 1 depending on the outcome of the following 6 pairwise comparisons.

1. $e_{x} vs. e_{0}$
2. $e_{y} vs. e_{0}$
3. $e_{x + y} vs. e_{0}$
4. $e_{y} vs. e_{x}$
5. $e_{x + y} vs. e_{x}$
6. $e_{x + y} vs. e_{y}$


#### No-Interaction Outcome Vectors
$5 / 75$ enumerated vectors are consistent with a lack of interaction. `isNI` returns `TRUE` if an outcome vector matches one of the 5 enumerated. Note that the current implementation of the main classification function will result in all null vectors being classified as "NI", whereas the other 4 are divided into 2 further conjugate pairs (Symmetric Left/Right and Step Up/Down).


```r
ni_list <- list(
sym.left = c(4, 2, 4, 2, "Symmetric Left"),
sym.right = c(2, 4, 2, 4, "Symmetric Right"),
step.up = c(2, 2, 4, 4, "Step Up"),
step.down = c(4, 4, 2, 2, "Step Down"),
null.vector = c(2, 2, 2, 2, "Null Vector")
)
```

![Barplots of theoretically enumerated 'non-interaction' profiles](figure/unnamed-chunk-3-1.png) ![Barplots of theoretically enumerated 'non-interaction' profiles](figure/unnamed-chunk-3-2.png) 

## Interaction Classes

### Interaction Class
Every gene can be classified into an *interaction class*. This is a ternary classification that describes the relationship between the sum of the fractional gene expression outputs of signal x and signal y (condition 2 - condition 1 + condition 3 - condition 1) and the fractional gene expression output of the signal combination (condition 4 - condition 1). More specifically, a negative interaction class exists if the sum of the fractional gene expression outputs of signal x and y are greater than the fractional gene expression of the signal combination and vice versa for a positive interaction class.

An interaction is categorized into a "no-interaction" class when the sum of fractional gene expression outputs is equal to the fractional gene expression output of the signal combination.

Thus, interaction class is present if one of the following inequalities are satisfied:

7. $e_{x+y} - e_{x} - e_{y} + e_{0} > 0$  *Positive*
8. $e_{x+y} - e_{x} - e_{y} + e_{0} < 0$  *Negative*
9. $e_{x+y} - e_{x} - e_{y} + e_{0} = 0$  *No-Interaction*

### Interaction Mode
Each interaction class can be further divided into mathematically and biologically defined interaction modes. There are 5 interaction modes associated with each interaction class and every mode has a conjugate pair in the opposite class. For example, there is a *positive synergy* mode associated with the positive interaction class. The gene expression outputs of a gene exhibiting the positive synergy mode satisfies the following inequalities.

$e_{x + y} >= max(e{0}, e_{x}, e_{y})$ and
$e_{x} - e_{0} > 0$ or $e_{y} - e_{0} > 0$

This interaction mode has a conjugate pair in the negative class, *negative synergy*. The gene expression outputs of a gene exhibiting the negative synergy mode satisfies the following inequalities.

$e_{x + y} <= min(e{0}, e_{x}, e_{y})$ and
$e_{x} - e_{0} < 0$ or $e_{y} - e_{0} < 0$

Note the inequalities that define each mode are the inverse of one another. These two modes, however, are commonly used classifications for drug interactions. Only 2/10 interaction modes were mentioned, suggesting that these interaction classes are heterogenous and allow for multiple interaction profile classifications.

41 outcome vectors are associated with a positive interaction mode and an equal number are associated with a negative interaction mode. There are 75 unique interaction outcome vectors, 5 of which are associated with a lack of interaction.$82 - 70 = 12$, thus 12 outcome vectors are associated with both interaction classes (positive or negative). We see this in [CITE MY SPREADSHEET].

Below, the outcome vectors associated with each interaction class and mode will be hard-coded. The outcome vector - mode associations were constructed from [CITE MY SPREADSHEET] and these results agreed with the full enumeration of interaction profiles in [CITE THEIR BARPLOTS]. Here, we implement a function to take a matrix of outcome vectors associated with an interaction class/mode as input and output the matrix of vectors that are associated with the opposite interaction class/mode.


```r
#' generate outcome vectors for the opposite class/mode
#'
#' @param a_matrix matrix of n rows where each row is a pairwise comparison outcome vector
#' with ncol(a_matrix) == 6
#' @return vecs matrix of n rows where each row maps to the opposite outcome vectors in
#' a_matrix
#'
#'  pwt_n <- enumerateThOutcomeVectors()
#'  oppositeVectors(pwt[10,])
oppositeVectors <- function(a_matrix) {
  # returns opposite sign
  reverseSign <- function(val) {
    if (val == 1) {
      return(-1)
    }
    if (val == -1) {
      return(1)
    }
    else {
      return(0)
    }
  }
  reversed <- apply(a_matrix, 1, function(x) (lapply(x, reverseSign) %>% unlist))
  vecs <- matrix(NaN, nrow = nrow(a_matrix), ncol = 6)
  for (i in 1:ncol(reversed)) {
    vecs[i,] <- reversed[,i]
  }
  return(vecs)
}
```

## Outcome vectors associated with the positive interaction class
There are 41 outcome vectors associated with 5 modes in the positive interaction class: low stabilization, x restores y, y restores x, positive synergy and emergent positive synergy.

### 1. Low Stabilization

$e_{0} > max(e_{x}, e_{y}, e_{x + y})$

The gene expression output of the vehicle condition (condition 1) is greater than the gene expression output of signal x, signal y, and signal x + y (conditions 1, 2, 3)

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
low_stab
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6]
##  [1,]   -1   -1   -1   -1   -1   -1
##  [2,]   -1   -1   -1   -1   -1    0
##  [3,]   -1   -1   -1   -1   -1    1
##  [4,]   -1   -1   -1   -1    0    1
##  [5,]   -1   -1   -1   -1    1    1
##  [6,]   -1   -1   -1    0   -1   -1
##  [7,]   -1   -1   -1    0    0    0
##  [8,]   -1   -1   -1    0    1    1
##  [9,]   -1   -1   -1    1   -1   -1
## [10,]   -1   -1   -1    1    0    1
## [11,]   -1   -1   -1    1    1   -1
## [12,]   -1   -1   -1    1    1    0
## [13,]   -1   -1   -1    1    1    1
```

### 2. X Restores Y

$e_{x} >= e_{x + y}$ and
$e_{x} > e_{y}$ and
$e_{x} >= e_{0}$

The gene expression output of signal x is greater than or equal to the output of the signal combination (condition 4), the output of signal x (condition 2) is greater than the output of signal y (condition 3), and the output of signal x (condition 2) is greater than or equal to the output of the vehicle (condition 1).

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
x_restores_y
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    0   -1   -1   -1   -1    1
## [2,]    0   -1    0   -1    0    1
## [3,]    1   -1   -1   -1   -1    1
## [4,]    1   -1    0   -1   -1    1
## [5,]    1   -1    1   -1   -1    1
## [6,]    1   -1    1   -1    0    1
```

### 3. Y Restores X
$e_{y} >= e_{x + y}$ and
$e_{y} > e_{x}$ and
$e_{y} >= e_{0}$

The gene expression output of signal y is greater than or equal to the gene expression output of the signal combination (condition 4) and the gene expression output of signal y is greater than the output of signal x, and the output of signal y (condition 3) is greater than or equal to the output of the vechicle (condition 1)

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
y_restores_x
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]   -1    0   -1    1    1   -1
## [2,]   -1    0    0    1    1    0
## [3,]   -1    1   -1    1    1   -1
## [4,]   -1    1    0    1    1   -1
## [5,]   -1    1    1    1    1   -1
## [6,]   -1    1    1    1    1    0
```
### 4. Positive Synergy

$e_{x + y} >= max(e{0}, e_{x}, e_{y})$ and
$e_{x} - e_{0} > 0$ or $e_{y} - e_{0} > 0$

The gene expression output of the signal combination (condition 4) is greater than or equal to the maximum output of the vehicle, signal x, and signal y (conditions 1, 2, 3) and the fractional output of signal x (condition 2 - condition 1) is greater than 0 or the fractional output of signal y (condition 3 - condition 1) is greater than 0.


```r
pos_syn
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]   -1    1    1    1    1    1
## [2,]    0    1    1    1    1    1
## [3,]    1   -1    1   -1    1    1
## [4,]    1    0    1   -1    1    1
## [5,]    1    1    1   -1    1    1
## [6,]    1    1    1    0    1    1
## [7,]    1    1    1    1    1    1
```
### 5. Emergent Positive Synergy
$e_{x + y} - e{0} > 0$ and
$sign(0) = +1$ and
$e_{x} - e_{0} <= 0$ & $e_{y} - e_{0} <= 0$

The fractional gene expression output of the signal combination (condition 4 - condition 1) is greater than 0 and the sign of 0 is 1 and the fractional gene expression outputs of both signal x and signal y are less than 0.

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
emer_pos_syn
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6]
##  [1,]   -1   -1    0   -1    1    1
##  [2,]   -1   -1    0    0    1    1
##  [3,]   -1   -1    0    1    1    1
##  [4,]   -1   -1    1   -1    1    1
##  [5,]   -1   -1    1    0    1    1
##  [6,]   -1   -1    1    1    1    1
##  [7,]   -1    0    1    1    1    1
##  [8,]    0   -1    1   -1    1    1
##  [9,]    0    0    1    0    1    1
```

### Positive Outcome Vectors

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
pos_ov
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6]
##  [1,]   -1   -1   -1   -1   -1   -1
##  [2,]   -1   -1   -1   -1   -1    0
##  [3,]   -1   -1   -1   -1   -1    1
##  [4,]   -1   -1   -1   -1    0    1
##  [5,]   -1   -1   -1   -1    1    1
##  [6,]   -1   -1   -1    0   -1   -1
##  [7,]   -1   -1   -1    0    0    0
##  [8,]   -1   -1   -1    0    1    1
##  [9,]   -1   -1   -1    1   -1   -1
## [10,]   -1   -1   -1    1    0    1
## [11,]   -1   -1   -1    1    1   -1
## [12,]   -1   -1   -1    1    1    0
## [13,]   -1   -1   -1    1    1    1
## [14,]    0   -1   -1   -1   -1    1
## [15,]    0   -1    0   -1    0    1
## [16,]    1   -1   -1   -1   -1    1
## [17,]    1   -1    0   -1   -1    1
## [18,]    1   -1    1   -1   -1    1
## [19,]    1   -1    1   -1    0    1
## [20,]   -1    0   -1    1    1   -1
## [21,]   -1    0    0    1    1    0
## [22,]   -1    1   -1    1    1   -1
## [23,]   -1    1    0    1    1   -1
## [24,]   -1    1    1    1    1   -1
## [25,]   -1    1    1    1    1    0
## [26,]   -1    1    1    1    1    1
## [27,]    0    1    1    1    1    1
## [28,]    1   -1    1   -1    1    1
## [29,]    1    0    1   -1    1    1
## [30,]    1    1    1   -1    1    1
## [31,]    1    1    1    0    1    1
## [32,]    1    1    1    1    1    1
## [33,]   -1   -1    0   -1    1    1
## [34,]   -1   -1    0    0    1    1
## [35,]   -1   -1    0    1    1    1
## [36,]   -1   -1    1   -1    1    1
## [37,]   -1   -1    1    0    1    1
## [38,]   -1   -1    1    1    1    1
## [39,]   -1    0    1    1    1    1
## [40,]    0   -1    1   -1    1    1
## [41,]    0    0    1    0    1    1
```

## Outcome vectors associated with the negative interaction class
There are 41 outcome vectors associated with 5 modes in the negative interaction class: high stabilization, x inhibits y, y inhibits x, negative synergy and emergent negative synergy.

### 6. High Stabilization

$e_{0} < max(e_{x}, e_{y}, e_{x + y})$

The gene expression output of the vehicle condition (condition 1) is less than the gene expression output of signal x, signal y, and signal x + y (conditions 1, 2, 3)

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
high_stab
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6]
##  [1,]    1    1    1    1    1    1
##  [2,]    1    1    1    1    1    0
##  [3,]    1    1    1    1    1   -1
##  [4,]    1    1    1    1    0   -1
##  [5,]    1    1    1    1   -1   -1
##  [6,]    1    1    1    0    1    1
##  [7,]    1    1    1    0    0    0
##  [8,]    1    1    1    0   -1   -1
##  [9,]    1    1    1   -1    1    1
## [10,]    1    1    1   -1    0   -1
## [11,]    1    1    1   -1   -1    1
## [12,]    1    1    1   -1   -1    0
## [13,]    1    1    1   -1   -1   -1
```

### 7. X Inhibits Y

$e_{x} <= e_{x + y}$ and
$e_{x} < e_{y}$ and
$e_{x} <= e_{0}$

The gene expression output of signal x is less than or equal to the output of the signal combination (condition 4), the output of signal x (condition 2) is less than the output of signal y (condition 3), and the output of signal x (condition 2) is less than or equal to the output of the vehicle (condition 1).

```r
# X Inhibits Y ----------------------------------------------------------------
x_inhibits_y <- oppositeVectors(x_restores_y)
x_inhibits_y
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    0    1    1    1    1   -1
## [2,]    0    1    0    1    0   -1
## [3,]   -1    1    1    1    1   -1
## [4,]   -1    1    0    1    1   -1
## [5,]   -1    1   -1    1    1   -1
## [6,]   -1    1   -1    1    0   -1
```

```r
x_inhibits_y_s <- apply(x_inhibits_y, 1, stringCoerce)
```

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

### 8. Y Inhibits X
$e_{y} <= e_{x + y}$ and
$e_{y} < e_{x}$ and
$e_{y} <= e_{0}$

The gene expression output of signal y is less than or equal to the gene expression output of the signal combination (condition 4) and the gene expression output of signal y is less than the output of signal x, and the output of signal y (condition 3) is less than or equal to the output of the vechicle (condition 1)

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
y_inhibits_x
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    1    0    1   -1   -1    1
## [2,]    1    0    0   -1   -1    0
## [3,]    1   -1    1   -1   -1    1
## [4,]    1   -1    0   -1   -1    1
## [5,]    1   -1   -1   -1   -1    1
## [6,]    1   -1   -1   -1   -1    0
```

### 9. Negative Synergy

$e_{x + y} <= min(e{0}, e_{x}, e_{y})$ and
$e_{x} - e_{0} < 0$ or $e_{y} - e_{0} < 0$

The gene expression output of the signal combination (condition 4) is less than or equal to the maximum output of the vehicle, signal x, and signal y (conditions 1, 2, 3) and the fractional output of signal x (condition 2 - condition 1) is less than 0 or the fractional output of signal y (condition 3 - condition 1) is less than 0.



```r
neg_syn
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    1   -1   -1   -1   -1   -1
## [2,]    0   -1   -1   -1   -1   -1
## [3,]   -1    1   -1    1   -1   -1
## [4,]   -1    0   -1    1   -1   -1
## [5,]   -1   -1   -1    1   -1   -1
## [6,]   -1   -1   -1    0   -1   -1
## [7,]   -1   -1   -1   -1   -1   -1
```

### 10. Emergent Negative Synergy
$e_{x + y} - e{0} < 0$ and
$sign(0) = -1$ and
$e_{x} - e_{0} >= 0$ & $e_{y} - e_{0} >= 0$

The fractional gene expression output of the signal combination (condition 4 - condition 1) is less than 0 and the sign of 0 is -1 and the fractional gene expression outputs of both signal x and signal y are greater than 0.

```r
emer_neg_syn <- oppositeVectors(emer_pos_syn)
emer_neg_syn
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6]
##  [1,]    1    1    0    1   -1   -1
##  [2,]    1    1    0    0   -1   -1
##  [3,]    1    1    0   -1   -1   -1
##  [4,]    1    1   -1    1   -1   -1
##  [5,]    1    1   -1    0   -1   -1
##  [6,]    1    1   -1   -1   -1   -1
##  [7,]    1    0   -1   -1   -1   -1
##  [8,]    0    1   -1    1   -1   -1
##  [9,]    0    0   -1    0   -1   -1
```

```r
emer_neg_syn_s <- apply(emer_neg_syn, 1, stringCoerce)
```

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

### Negative Outcome Vectors

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
neg_ov
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6]
##  [1,]    1    1    1    1    1    1
##  [2,]    1    1    1    1    1    0
##  [3,]    1    1    1    1    1   -1
##  [4,]    1    1    1    1    0   -1
##  [5,]    1    1    1    1   -1   -1
##  [6,]    1    1    1    0    1    1
##  [7,]    1    1    1    0    0    0
##  [8,]    1    1    1    0   -1   -1
##  [9,]    1    1    1   -1    1    1
## [10,]    1    1    1   -1    0   -1
## [11,]    1    1    1   -1   -1    1
## [12,]    1    1    1   -1   -1    0
## [13,]    1    1    1   -1   -1   -1
## [14,]    0    1    1    1    1   -1
## [15,]    0    1    0    1    0   -1
## [16,]   -1    1    1    1    1   -1
## [17,]   -1    1    0    1    1   -1
## [18,]   -1    1   -1    1    1   -1
## [19,]   -1    1   -1    1    0   -1
## [20,]    1    0    1   -1   -1    1
## [21,]    1    0    0   -1   -1    0
## [22,]    1   -1    1   -1   -1    1
## [23,]    1   -1    0   -1   -1    1
## [24,]    1   -1   -1   -1   -1    1
## [25,]    1   -1   -1   -1   -1    0
## [26,]    1   -1   -1   -1   -1   -1
## [27,]    0   -1   -1   -1   -1   -1
## [28,]   -1    1   -1    1   -1   -1
## [29,]   -1    0   -1    1   -1   -1
## [30,]   -1   -1   -1    1   -1   -1
## [31,]   -1   -1   -1    0   -1   -1
## [32,]   -1   -1   -1   -1   -1   -1
## [33,]    1    1    0    1   -1   -1
## [34,]    1    1    0    0   -1   -1
## [35,]    1    1    0   -1   -1   -1
## [36,]    1    1   -1    1   -1   -1
## [37,]    1    1   -1    0   -1   -1
## [38,]    1    1   -1   -1   -1   -1
## [39,]    1    0   -1   -1   -1   -1
## [40,]    0    1   -1    1   -1   -1
## [41,]    0    0   -1    0   -1   -1
```

### No-Interaction Vectors

```
## Error in match.fun(FUN): object 'stringCoerce' not found
```

```r
ni_ov
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    0    0    0    0    0    0
## [2,]    1    0    1   -1    0    1
## [3,]    0    1    1    1    1    0
## [4,]   -1    0   -1    1    0   -1
## [5,]    0   -1   -1   -1   -1    0
```

### Null Vector

```
## Error in eval(expr, envir, enclos): could not find function "stringCoerce"
```

### No Interaction Vectors

```r
sym_right <- c(1, 0, 1, -1, 0, 1)
sym_left <- c(-1, 0, -1, 1, 0, -1)

# step up and step down
# interaction profiles look like increasing or decreasing step functions
step_up <- c(0, 1, 1, 1, 1, 0)
step_down <- c(0, -1, -1, -1, -1, 0)
```




```r
sessionInfo()
```

```
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_1.0.1 magrittr_1.5  knitr_1.10.5 
## 
## loaded via a namespace (and not attached):
##  [1] MASS_7.3-40      Rcpp_0.11.6      colorspace_1.2-6 digest_0.6.8    
##  [5] evaluate_0.7     formatR_1.2      grid_3.1.1       gtable_0.1.2    
##  [9] munsell_0.4.2    plyr_1.8.2       proto_0.3-10     reshape2_1.4.1  
## [13] scales_0.2.4     stringi_0.4-1    stringr_1.0.0    tools_3.1.1
```


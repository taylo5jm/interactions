Enumeration of Theoretical Outcome Vectors
====================================================


```
## Error in file(filename, "r", encoding = encoding): cannot open the connection
```
## Enumerate Theoretical Outcome Vectors with n conditions and k ranks

Suppose a set of 4 gene expression outputs can be coded as a 6-tuple with repetition, where the elements of the 6-tuple can assume an element in the set $s = (-1, 0, 1)$. Each element in the set of 4 gene expression outputs can assume an element in the set $k = (1, 2, 3, 4)$ where $k$ is the set of possible ranks that a gene expression output can assume.

$4^{4} = 256$ permutations with repitition exist if these assumptions are satisfied. We can then take the ranks of each permutation and average the ties. A gene expression set $g = (1, 1, 1, 1)$ is the same as the set $h = (3, 3, 3, 3)$, when we are considering interaction effects present between the gene expression outputs. By taking the ranks, the number of unique interaction profiles is reduced to 75.  Informally, we can summarize this process as:
$4^{4} = 256$ permutations --> ranks --> remove all duplicated permutations with repetition.

#### Pairwise Comparison Outcome Vectors
Unique interaction profiles can be coded into *pairwise comparison outcome vectors* by comparing the values of the gene expression outputs in 6 pairwise comparisons.

1. $e_{x} vs. e_{0}$
2. $e_{y} vs. e_{0}$
3. $e_{x + y} vs. e_{0}$
4. $e_{y} vs. e_{x}$
5. $e_{x + y} vs. e_{x}$
6. $e_{x + y} vs. e_{y}$

where $e_{0}$, $e_{x}$, $e_{y}$, $e_{x + y}$ are the gene expression outputs for the null, signal $x$, signal $y$, and signal $x + y$. The elements of the 6-tuple can be equal to $-1$, $0$, or $1$ depending on whether the gene expression output of the first condition is greater than the the gene expression output of the second condition.

In this implementation, some functions are defined to take the row ranks of matrix and make the 6 pairwise comparisons for the k = 4 conditions we are interested in.


```r
# rank tuples by row and allow ties
rankTuples <- function(x) {
  ranks <- matrix(0, nrow(x), ncol(x))
  for (i in 1:nrow(x)) {
    ranks[i,] <- rank(x[i,], ties.method = "average")
  }
  return(ranks)
}

# make pairwise comparisons by directly comparing values
pwCompare <- function(v) {
    compare <- function(val1, val2) {
      if (val1 < val2) {
        return(-1)
      }
      if (val1 == val2) {
        return(0)
      }
      if (val1 > val2) {
        return(1)
      }
    }

    pwcomp <- matrix(0, nrow(v), ncol = 6)
    for (i in 1:nrow(v)) {

      e0 <- v[i,1]
      ex <- v[i,2]
      ey <- v[i,3]
      exy <- v[i,4]

      # ex vs. e0
      pwcomp[i,1] <- compare(ex, e0)
      # ey vs. e0
      pwcomp[i,2] <- compare(ey, e0)
      # exy vs. e0
      pwcomp[i,3] <- compare(exy, e0)
      # ey vs. ex
      pwcomp[i,4] <- compare(ey, ex)
      # exy vs. ex
      pwcomp[i,5] <- compare(exy, ex)
      # exy vs. ey
      pwcomp[i,6] <- compare(exy, ey)

    }
    return(pwcomp)
}
```

We then enumerate the theoretical outcome vectors, as described above, when n = 4 conditions are present and can assume k = (1, 2, 3, 4).

```r
#' Enumerate theoretical pairwise comparison outcome vectors given n conditions and k ranks
#'
#' @return matrix where each row represents a pairwise comparison outcome vector
#' @examples
#' enumerateThOutcomeVectors() %>% nrow
enumerateThOutcomeVectors <- function() {
  pwt_n <- gtools::permutations(4, 4, repeats.allowed=TRUE) %>%
    rankTuples() %>% unique() %>%
    pwCompare()
  # pairwise comparison outcome vectors as character vector
  pwt_s <<- apply(pwt_n, 1, stringCoerce)
  return(pwt_n)
}

pwt_n <- enumerateThOutcomeVectors()
```


```r
pheatmap(pwt_n, color = c("#000099", "#FFFFFF", "#990000"), cluster_rows = FALSE, cluster_cols = FALSE)
```

<img src="figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

There should be 75 unique pairwise comparison outcome vectors, when 6 comparisons can take on values in the set $s = (-1, 0, 1)$


```r
dim(pwt_n)
```

```
## [1] 75  6
```

```r
pwt_n
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6]
##  [1,]    0    0    0    0    0    0
##  [2,]    0    0    1    0    1    1
##  [3,]    0    1    0    1    0   -1
##  [4,]    0    1    1    1    1    0
##  [5,]    0    1    1    1    1    1
##  [6,]    0    1    1    1    1   -1
##  [7,]    1    0    0   -1   -1    0
##  [8,]    1    0    1   -1    0    1
##  [9,]    1    0    1   -1    1    1
## [10,]    1    1    0    0   -1   -1
## [11,]    1    1    1    0    0    0
## [12,]    1    1    1    0    1    1
## [13,]    1    1    0    1   -1   -1
## [14,]    1    1    1    1    0   -1
## [15,]    1    1    1    1    1    0
## [16,]    1    1    1    1    1    1
## [17,]    1    1    1    1    1   -1
## [18,]    1    0    1   -1   -1    1
## [19,]    1    1    0   -1   -1   -1
## [20,]    1    1    1   -1   -1    0
## [21,]    1    1    1   -1    0    1
## [22,]    1    1    1   -1    1    1
## [23,]    1    1    1    0   -1   -1
## [24,]    1    1    1    1   -1   -1
## [25,]    1    1    1   -1   -1    1
## [26,]    1    1    1   -1   -1   -1
## [27,]   -1   -1   -1    0    0    0
## [28,]   -1   -1    0    0    1    1
## [29,]   -1   -1    1    0    1    1
## [30,]   -1    0   -1    1    0   -1
## [31,]   -1    0    0    1    1    0
## [32,]   -1    0    1    1    1    1
## [33,]   -1    1   -1    1    0   -1
## [34,]   -1    1    0    1    1   -1
## [35,]   -1    1    1    1    1    0
## [36,]   -1    1    1    1    1    1
## [37,]   -1    1    1    1    1   -1
## [38,]    0   -1   -1   -1   -1    0
## [39,]    0   -1    0   -1    0    1
## [40,]    0   -1    1   -1    1    1
## [41,]    0    0   -1    0   -1   -1
## [42,]    0    1   -1    1   -1   -1
## [43,]    1   -1   -1   -1   -1    0
## [44,]    1   -1    0   -1   -1    1
## [45,]    1   -1    1   -1    0    1
## [46,]    1   -1    1   -1    1    1
## [47,]    1    0   -1   -1   -1   -1
## [48,]    1    1   -1    0   -1   -1
## [49,]    1    1   -1    1   -1   -1
## [50,]    1   -1    1   -1   -1    1
## [51,]    1    1   -1   -1   -1   -1
## [52,]   -1   -1   -1    0    1    1
## [53,]   -1   -1   -1    1    0   -1
## [54,]   -1   -1   -1    1    1    0
## [55,]   -1   -1    0    1    1    1
## [56,]   -1   -1    1    1    1    1
## [57,]   -1    0   -1    1    1   -1
## [58,]   -1    1   -1    1    1   -1
## [59,]   -1   -1   -1   -1   -1    0
## [60,]   -1   -1   -1   -1    0    1
## [61,]   -1   -1    0   -1    1    1
## [62,]   -1   -1    1   -1    1    1
## [63,]   -1   -1   -1    0   -1   -1
## [64,]   -1    0   -1    1   -1   -1
## [65,]   -1    1   -1    1   -1   -1
## [66,]    0   -1   -1   -1   -1    1
## [67,]    0   -1   -1   -1   -1   -1
## [68,]    1   -1   -1   -1   -1    1
## [69,]    1   -1   -1   -1   -1   -1
## [70,]   -1   -1   -1    1    1    1
## [71,]   -1   -1   -1    1    1   -1
## [72,]   -1   -1   -1   -1    1    1
## [73,]   -1   -1   -1    1   -1   -1
## [74,]   -1   -1   -1   -1   -1    1
## [75,]   -1   -1   -1   -1   -1   -1
```

This result is also supported by an implementation of the closed-form recurrence relation solution found in "CITE ME" and implemented by Dr. Stephen Hoang. In this implementation, the recurrence relation solution is used to compute the number of possible vectors when there are n = 4 conditions and k = (1, 2, 3, 4)

```r
# Stephen Hoang
A <- function(n, k) {
  if (n < 1 | k < 1) {
    return(0)
  }
  if (k > n) {
    return(0)
  }
  if (n == 1 & k == 1) {
    return(1)
  }
  res <- k * (A(n - 1, k) + A(n - 1, k - 1))
  return(res)
}

k <- 1:4
lapply(k, FUN = A, n = 4) %>% Reduce(f = "+")
```

```
## [1] 75
```

## Null Distribution

Suppose we are interested in the null distribution of interaction classes and modes. First, we will generate a pseudorandom data set in the form of a 100,000 x 4 matrix, where each row represents a "gene" and each column represents a condition. We will actually generate three of these matrices that follow the Poisson, uniform and normal distribution. Each element in the pseudorandom matrices will be rounded to simulate a discrete distribution of RNA-seq counts.

The rounding also provides an element of statistics in the simulation. Recall there are 3 interaction classes, one of which being the null class.

7. $e_{x+y} - e_{x} - e_{y} + e_{0} > 0$  *Positive*
8. $e_{x+y} - e_{x} - e_{y} + e_{0} < 0$  *Negative*
9. $e_{x+y} - e_{x} - e_{y} + e_{0} = 0$  *No-Interaction*

If we were to continue this simulation without rounding, then the "no-interaction" class would essentially be non-existent, since the probability that any set of gene expression outputs satisfies the "no-interaction" inequality is essentially 0.


```r
library(magrittr)
set.seed(1983)
rlist <- list(
  unif_matrix = runif(100000, min = 2, max = 16) %>% matrix(nrow = 100000, ncol = 4, byrow = TRUE) %>% round,
  norm_matrix = rnorm(100000, 8, 2) %>% matrix(nrow = 100000, ncol = 4, byrow = TRUE) %>% round,
  pois_matrix = rpois(100000, lambda = 8) %>% matrix(nrow = 100000, ncol = 4, byrow = TRUE) %>% round
)
```

We can code the outcome vectors with `pwCompare()`

```r
outcome_vectors <- lapply(rlist, pwCompare)
outcome_vector_tbl <- lapply(outcome_vectors, ovTable, get.class = FALSE)
```



```r
o <- outcome_vector_tbl
sum_table <- data.frame(o[[1]]$Outcome.Vector, o[[1]]$Freq, o[[2]]$Freq, o[[3]]$Freq,
                        o[[1]]$Log.Freq, o[[2]]$Log.Freq, o[[3]]$Log.Freq,
                        o[[1]]$Prop, o[[2]]$Prop, o[[3]]$Prop)
dists <- c("Unif", "Norm", "Pois")
vars <- c("Freq", "Log.Freq", "Prop")
cols <- list()
for (i in 1:length(vars)) {
  cols[[i]] <- lapply(dists, paste, vars[i], sep = ".") %>% unlist
}
colnames(sum_table) <- c("Outcome.Vector", unlist(cols))
sum_table[,5:10] <- signif(sum_table[,5:10], digits = 3)
sum_table
```

```
##       Outcome.Vector Unif.Freq Norm.Freq Pois.Freq Unif.Log.Freq
## 1  -1,-1,-1,-1,-1,-1      2656      1620      2168          7.88
## 2   -1,-1,-1,-1,-1,0       920      1232      1228          6.82
## 3   -1,-1,-1,-1,-1,1      2640      1516      1996          7.88
## 4    -1,-1,-1,-1,0,1       972      1556      1324          6.88
## 5    -1,-1,-1,-1,1,1      2556      1656      2140          7.85
## 6   -1,-1,-1,0,-1,-1       916      1700      1352          6.82
## 7     -1,-1,-1,0,0,0       244       988       528          5.50
## 8     -1,-1,-1,0,1,1       860      1236      1084          6.76
## 9   -1,-1,-1,1,-1,-1      2712      1616      2120          7.91
## 10   -1,-1,-1,1,0,-1       932      1228      1104          6.84
## 11   -1,-1,-1,1,1,-1      2660      1548      2144          7.89
## 12    -1,-1,-1,1,1,0       832      1640      1496          6.72
## 13    -1,-1,-1,1,1,1      2504      1628      2084          7.83
## 14    -1,-1,0,-1,1,1       836      1308       928          6.73
## 15     -1,-1,0,0,1,1       224       788       428          5.41
## 16     -1,-1,0,1,1,1       832      1252       916          6.72
## 17    -1,-1,1,-1,1,1      2664      1644      2212          7.89
## 18     -1,-1,1,0,1,1       920      1184      1108          6.82
## 19     -1,-1,1,1,1,1      2544      1576      2140          7.84
## 20   -1,0,-1,1,-1,-1       840      1232       996          6.73
## 21    -1,0,-1,1,0,-1       232       824       444          5.45
## 22    -1,0,-1,1,1,-1       916      1124       940          6.82
## 23      -1,0,0,1,1,0       204       928       556          5.32
## 24      -1,0,1,1,1,1      1004      1832      1304          6.91
## 25   -1,1,-1,1,-1,-1      2592      1444      2124          7.86
## 26    -1,1,-1,1,0,-1       860      1308      1076          6.76
## 27    -1,1,-1,1,1,-1      2704      1668      2232          7.90
## 28     -1,1,0,1,1,-1      1052      1912      1336          6.96
## 29     -1,1,1,1,1,-1      2844      1504      2208          7.95
## 30      -1,1,1,1,1,0       916      1140       992          6.82
## 31      -1,1,1,1,1,1      2632      1592      2228          7.88
## 32  0,-1,-1,-1,-1,-1       932      1084      1004          6.84
## 33   0,-1,-1,-1,-1,0       216       764       392          5.38
## 34   0,-1,-1,-1,-1,1       916       988       928          6.82
## 35     0,-1,0,-1,0,1       232       908       540          5.45
## 36     0,-1,1,-1,1,1       920      1804      1468          6.82
## 37    0,0,-1,0,-1,-1       216       940       492          5.38
## 38       0,0,0,0,0,0        36       340       168          3.58
## 39       0,0,1,0,1,1       252       884       576          5.53
## 40    0,1,-1,1,-1,-1       896      1640      1372          6.80
## 41      0,1,0,1,0,-1       248      1008       504          5.51
## 42      0,1,1,1,1,-1       932      1144      1076          6.84
## 43       0,1,1,1,1,0       208       832       420          5.34
## 44       0,1,1,1,1,1       832      1256      1116          6.72
## 45  1,-1,-1,-1,-1,-1      2736      1476      2000          7.91
## 46   1,-1,-1,-1,-1,0       904      1284      1148          6.81
## 47   1,-1,-1,-1,-1,1      2684      1612      2128          7.90
## 48    1,-1,0,-1,-1,1       892      1704      1448          6.79
## 49    1,-1,1,-1,-1,1      2676      1612      2256          7.89
## 50     1,-1,1,-1,0,1       940      1248       980          6.85
## 51     1,-1,1,-1,1,1      2592      1628      2180          7.86
## 52   1,0,-1,-1,-1,-1      1096      1788      1384          7.00
## 53     1,0,0,-1,-1,0       196       896       500          5.28
## 54     1,0,1,-1,-1,1       884      1136      1100          6.78
## 55      1,0,1,-1,0,1       224       832       460          5.41
## 56      1,0,1,-1,1,1       812      1236      1136          6.70
## 57   1,1,-1,-1,-1,-1      2608      1404      2176          7.87
## 58    1,1,-1,0,-1,-1       900      1192       988          6.80
## 59    1,1,-1,1,-1,-1      2668      1564      1912          7.89
## 60    1,1,0,-1,-1,-1       884      1224      1016          6.78
## 61     1,1,0,0,-1,-1       212       872       476          5.36
## 62     1,1,0,1,-1,-1       968      1104      1184          6.88
## 63    1,1,1,-1,-1,-1      2632      1644      2124          7.88
## 64     1,1,1,-1,-1,0       920      1812      1540          6.82
## 65     1,1,1,-1,-1,1      2736      1588      2028          7.91
## 66      1,1,1,-1,0,1       888      1304       928          6.79
## 67      1,1,1,-1,1,1      2712      1588      2172          7.91
## 68     1,1,1,0,-1,-1       972      1236       924          6.88
## 69       1,1,1,0,0,0       272       876       468          5.61
## 70       1,1,1,0,1,1       964      1468      1488          6.87
## 71     1,1,1,1,-1,-1      2716      1592      2140          7.91
## 72      1,1,1,1,0,-1       932      1672      1396          6.84
## 73      1,1,1,1,1,-1      2836      1556      2144          7.95
## 74       1,1,1,1,1,0       856      1264       992          6.75
## 75       1,1,1,1,1,1      2632      1572      2192          7.88
##    Norm.Log.Freq Pois.Log.Freq Unif.Prop Norm.Prop Pois.Prop
## 1           7.39          7.68   0.02660   0.01620   0.02170
## 2           7.12          7.11   0.00920   0.01230   0.01230
## 3           7.32          7.60   0.02640   0.01520   0.02000
## 4           7.35          7.19   0.00972   0.01560   0.01320
## 5           7.41          7.67   0.02560   0.01660   0.02140
## 6           7.44          7.21   0.00916   0.01700   0.01350
## 7           6.90          6.27   0.00244   0.00988   0.00528
## 8           7.12          6.99   0.00860   0.01240   0.01080
## 9           7.39          7.66   0.02710   0.01620   0.02120
## 10          7.11          7.01   0.00932   0.01230   0.01100
## 11          7.34          7.67   0.02660   0.01550   0.02140
## 12          7.40          7.31   0.00832   0.01640   0.01500
## 13          7.40          7.64   0.02500   0.01630   0.02080
## 14          7.18          6.83   0.00836   0.01310   0.00928
## 15          6.67          6.06   0.00224   0.00788   0.00428
## 16          7.13          6.82   0.00832   0.01250   0.00916
## 17          7.40          7.70   0.02660   0.01640   0.02210
## 18          7.08          7.01   0.00920   0.01180   0.01110
## 19          7.36          7.67   0.02540   0.01580   0.02140
## 20          7.12          6.90   0.00840   0.01230   0.00996
## 21          6.71          6.10   0.00232   0.00824   0.00444
## 22          7.02          6.85   0.00916   0.01120   0.00940
## 23          6.83          6.32   0.00204   0.00928   0.00556
## 24          7.51          7.17   0.01000   0.01830   0.01300
## 25          7.28          7.66   0.02590   0.01440   0.02120
## 26          7.18          6.98   0.00860   0.01310   0.01080
## 27          7.42          7.71   0.02700   0.01670   0.02230
## 28          7.56          7.20   0.01050   0.01910   0.01340
## 29          7.32          7.70   0.02840   0.01500   0.02210
## 30          7.04          6.90   0.00916   0.01140   0.00992
## 31          7.37          7.71   0.02630   0.01590   0.02230
## 32          6.99          6.91   0.00932   0.01080   0.01000
## 33          6.64          5.97   0.00216   0.00764   0.00392
## 34          6.90          6.83   0.00916   0.00988   0.00928
## 35          6.81          6.29   0.00232   0.00908   0.00540
## 36          7.50          7.29   0.00920   0.01800   0.01470
## 37          6.85          6.20   0.00216   0.00940   0.00492
## 38          5.83          5.12   0.00036   0.00340   0.00168
## 39          6.78          6.36   0.00252   0.00884   0.00576
## 40          7.40          7.22   0.00896   0.01640   0.01370
## 41          6.92          6.22   0.00248   0.01010   0.00504
## 42          7.04          6.98   0.00932   0.01140   0.01080
## 43          6.72          6.04   0.00208   0.00832   0.00420
## 44          7.14          7.02   0.00832   0.01260   0.01120
## 45          7.30          7.60   0.02740   0.01480   0.02000
## 46          7.16          7.05   0.00904   0.01280   0.01150
## 47          7.39          7.66   0.02680   0.01610   0.02130
## 48          7.44          7.28   0.00892   0.01700   0.01450
## 49          7.39          7.72   0.02680   0.01610   0.02260
## 50          7.13          6.89   0.00940   0.01250   0.00980
## 51          7.40          7.69   0.02590   0.01630   0.02180
## 52          7.49          7.23   0.01100   0.01790   0.01380
## 53          6.80          6.21   0.00196   0.00896   0.00500
## 54          7.04          7.00   0.00884   0.01140   0.01100
## 55          6.72          6.13   0.00224   0.00832   0.00460
## 56          7.12          7.04   0.00812   0.01240   0.01140
## 57          7.25          7.69   0.02610   0.01400   0.02180
## 58          7.08          6.90   0.00900   0.01190   0.00988
## 59          7.36          7.56   0.02670   0.01560   0.01910
## 60          7.11          6.92   0.00884   0.01220   0.01020
## 61          6.77          6.17   0.00212   0.00872   0.00476
## 62          7.01          7.08   0.00968   0.01100   0.01180
## 63          7.40          7.66   0.02630   0.01640   0.02120
## 64          7.50          7.34   0.00920   0.01810   0.01540
## 65          7.37          7.61   0.02740   0.01590   0.02030
## 66          7.17          6.83   0.00888   0.01300   0.00928
## 67          7.37          7.68   0.02710   0.01590   0.02170
## 68          7.12          6.83   0.00972   0.01240   0.00924
## 69          6.78          6.15   0.00272   0.00876   0.00468
## 70          7.29          7.31   0.00964   0.01470   0.01490
## 71          7.37          7.67   0.02720   0.01590   0.02140
## 72          7.42          7.24   0.00932   0.01670   0.01400
## 73          7.35          7.67   0.02840   0.01560   0.02140
## 74          7.14          6.90   0.00856   0.01260   0.00992
## 75          7.36          7.69   0.02630   0.01570   0.02190
```


```
## Warning: package 'ggplot2' was built under R version 3.1.3
```

```
## Warning: package 'dplyr' was built under R version 3.1.2
```

```
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
```

<img src="figure/unnamed-chunk-10-1.png" title="Frequency distributions of outcome vectors coded from Uniform, Normal, and Poisson distribution pseudorandom matrices." alt="Frequency distributions of outcome vectors coded from Uniform, Normal, and Poisson distribution pseudorandom matrices." style="display: block; margin: auto;" />

## Null Distribution: Class

We then implement a function to classify an interaction by an interaction class.

```r
classifyByClass <- function(v) {
  class_classifier <- function(v) {
    e0 <- v[1]
    ex <- v[2]
    ey <- v[3]
    exy <- v[4]
    if (exy - ex - ey + e0 > 0) {
      return("Positive")
    }
    if (exy - ex - ey + e0 < 0) {
      return("Negative")
    }
    else {
      return('NI')
    }
  }

  if (is.numeric(v)) {
    class_classifier(v)
  }
  if (is.matrix(v)) {
    apply(v, 1, class_classifier)
  }
}

rlist_classes <- lapply(rlist, classifyByClass)
class_tbl <- lapply(rlist_classes, ovTable, get.class = FALSE)
```

PLOT CLASS DISTRIBUTION HERE

```r
p1 <- qplot(x = 1:nrow(class_tbl[[1]]), y = Freq, data = class_tbl[[1]],
           geom = "bar", stat = "identity", xlab = "Outcome Vector", ylab = "Frequency",
           main = "Frequency Distribution of Interaction Modes") + scale_x_discrete(breaks = 1:10)
p2 <- qplot(x = 1:nrow(class_tbl[[2]]), y = Freq, data = class_tbl[[2]],
            geom = "bar", stat = "identity", xlab = "Outcome Vector", ylab = "Frequency",
            main = "Frequency Distribution of Interaction Modes") +
           scale_x_discrete(breaks = 1:10)
p3 <- qplot(x = 1:nrow(class_tbl[[3]]), y = Freq, data = class_tbl[[3]],
            geom = "bar", stat = "identity", xlab = "Outcome Vector", ylab = "Frequency",
            main = "Frequency Distribution of Interaction Modes") + scale_x_discrete(breaks = 1:10)
multiplot(p1, p2, p3, cols = 1)
```

<img src="figure/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" style="display: block; margin: auto;" />

Finally, we classify each gene into an interaction mode. This function classifies each gene into an interaction class based on the inequalities that define each interaction mode.


```r
library(magrittr)
classifyByMode <- function(v) {

  e0 <- v[1]
  ex <- v[2]
  ey <- v[3]
  exy <- v[4]

  # positive -----------------------------------------------------------
  if (exy - ex - ey + e0 > 0) {
    if (e0 > max(ex, ey, exy)) {
      return("Low.Stab")
    }
    if (ex >= exy & ex > ey & ex >= e0) {
      return("X.Restores.Y")
    }
    if (ey >= exy & ey > ex & ey >= e0) {
      return("Y.Restores.X")
    }
    if (exy >= max(e0, ex, ey) & ex - e0 > 0 | ey - e0 > 0) {
        return("Pos.Syn")
    }
    if (sign(exy - e0) == 1 & sign(e0) == 1 & ex - e0 <= 0 & ey - e0 <= 0) {
      return("Emer.Pos.Syn")
    }

    else {
      return('unsorted')
    }
  }

  # negative ------------------------------------------------------------------
  if (exy - ex - ey + e0 < 0) {

    if (e0 < min(ex, ey, exy)) {
      return("High.Stab")
    }
    if (ex <= exy & ex < ey & ex <= e0) {
      return("X.Inhibits.Y")
    }
    if (ey <= exy & ey < ex & ey <= e0) {
      return("Y.Inhibits.X")
    }
    if (exy < min(e0, ex, ey) & ex - e0 < 0 | ey - e0 < 0) {
      return("Neg.Syn")
    }
    if (sign(exy - e0) == -1 & sign(e0) == 1 & ex - e0 >= 0 & ey - e0 >= 0) {
      return("Emer.Neg.Syn")
    }

    else {
      return("unsorted")
    }
  }

  # no interaction ---------------------------------------------------------------
  if (exy - ex - ey + e0 == 0) {
    return('NI')
  }

  # never executed
  else {
    return('error')
  }
}
```


```r
rclasses <- lapply(rlist, classifyByMode)
```


```r
# PLOT MODES DISTRIBUTION
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
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] dplyr_0.4.1     ggplot2_1.0.1   d3heatmap_0.6.1 pheatmap_1.0.2 
## [5] magrittr_1.5    knitr_1.10.5   
## 
## loaded via a namespace (and not attached):
##  [1] DBI_0.3.1          MASS_7.3-40        RColorBrewer_1.1-2
##  [4] Rcpp_0.11.6        assertthat_0.1     base64enc_0.1-2   
##  [7] colorspace_1.2-6   digest_0.6.8       evaluate_0.7      
## [10] formatR_1.2        gtable_0.1.2       htmltools_0.2.6   
## [13] htmlwidgets_0.4    labeling_0.3       lazyeval_0.1.10   
## [16] munsell_0.4.2      parallel_3.1.1     plyr_1.8.2        
## [19] png_0.1-7          proto_0.3-10       reshape2_1.4.1    
## [22] scales_0.2.4       stringi_0.4-1      stringr_1.0.0     
## [25] tools_3.1.1
```

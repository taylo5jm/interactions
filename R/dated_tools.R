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

pwCompare2 <- function(mean_matrix, sd_matrix) {

  compare <- function(val1, val2, sd) {

    if (val1 < (val2 - (2*sd))) {
      return(-1)
    }

    if (val1 > (val2 + (2*sd))) {
      return(1)
    }

    if (val1 > val2 - (2*sd) & val1 < val2 + (2*sd)) {
      return(0)
    }
  }

  pwcomp <- matrix(-2, nrow(mean_matrix), ncol = 6)

  for (i in 1:nrow(mean_matrix)) {

    e0 <- mean_matrix[i,1]
    ex <- mean_matrix[i,2]
    ey <- mean_matrix[i,3]
    exy <- mean_matrix[i,4]

    pwcomp[i,1] <- compare(ex, e0, sd_matrix[i, 1])

    # ey vs. e0
    pwcomp[i,2] <- compare(ey, e0, sd_matrix[i, 1])

    # exy vs. e0
    pwcomp[i,3] <- compare(exy, e0, sd_matrix[i, 1])

    # ey vs. ex
    pwcomp[i,4] <- compare(ey, ex, sd_matrix[i, 2])

    # exy vs. ex
    pwcomp[i,5] <- compare(exy, ex, sd_matrix[i, 2])

    # exy vs. ey
    pwcomp[i,6] <- compare(exy, ey, sd_matrix[i, 3])
  }

  return(pwcomp)
}

stringCoerce <- function(v) {
  as.character(v) %>%
    paste(collapse = ",", sep = " ") %>%
    return()
}

rankTuples <- function(x) {

  ranks <- matrix(0, nrow(x), ncol(x))

  for (i in 1:nrow(x)) {
    ranks[i,] <- rank(x[i,], ties.method = "average")
  }

  return(ranks)
}

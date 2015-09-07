#' Coerce a numeric vector to a string the wrong way
#'
#' @param v numeric vector
#' @return character vector of length n = 1
#' @export stringCoerce

stringCoerce <- function(v) {
    as.character(v) %>%
      paste(collapse = ",", sep = " ") %>%
      return()
}

#' Compute row ranks and average ties (should be deprecated)
#'
#' @param x matrix of dimensions i x j where the row ranks should be computed
#' @return matrix where each row are the ranks of the original values
#' @export rankTuples
#' @examples
#' library(magrittr)
#' rpois(16, 3) %>% matrix(nrow = 4, ncol = 4) %>% rankTuples

rankTuples <- function(x) {

  ranks <- matrix(0, nrow(x), ncol(x))
  for (i in 1:nrow(x)) {
    ranks[i,] <- rank(x[i,], ties.method = "average")
  }

  return(ranks)
}


#' Determine if an outcome vector is enumerated theoretically
#'
#' @param v outcome vector (numeric) of length n = 6
#' @param th.vectors matrix where each row is an enumerated outcome vector
#' @return Boolean TRUE when outcome vector is not enumerated as part of the 75 theoretical outcome vectors
#' @export isAnomaly
isAnomaly <- function(v, th.vectors) {

  if (any(apply(th.vectors$all_ov, 1, function(x) (all(v == x))) == TRUE)) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}

#' check outcome vector to see if it is consistent with 1/5 no-interaction vectors
#'
#' @param v 6-dimensional outcome vector (numeric vector where length(v) == 6)
#' @return Boolean TRUE if outcome vector matches 1/5 outcome vectors that are consistent
#' with a lack of interaction
#' @export isNi
isNi <- function(v) {

  if (is.numeric(v)) {
    if (all(v == rep(0, 6)) | all(v == c(1, 0, 1, -1, 0, 1)) |
          all(v == c(0, 1, 1, 1, 1, 0)) | all(v == c(-1, 0, -1, 1, 0, -1)) |
          all(v == c(0, -1, -1, -1, -1, 0))) {
      return(TRUE)
    }
    else {
      return(FALSE)
    }
  }
}


#' Check to see if an outcome vector is associated with a particular interaction mode
#'
#' @param ov numeric or character vector of length n = 1 indicating outcome vector
#' @param pwt_n matrix where each row is an outcome vector associated with a particular classification (mode, class)
#' @return Boolean TRUE when outcome vector is associated with a particular mode
#' @export isMode

isMode <- function(ov, pwt_n) {
  # ov is vector coerced by stringCoerce()
  if (is.character(ov)) {
    if (ov %in% pwt_n) {
      return(TRUE)
    }
    else {
      return(FALSE)
    }
  }
  # interaction mode w/ > 1 associated outcome vectors
  if (is.matrix(pwt_n)) {
    apply(pwt_n, 1, function(x) (all(ov == x))) %>% any
  }
  # interaction mode with 1 associated outcome vector
  else {
    return(all(ov == pwt_n))
  }
}

#' Code elements of pairwise outcome vectors with confidence intervals from limma
#'
#' @param condition1 matrix where conditions 1:3 are the log fold changes
#' @param condition2 ""
#' @return numeric vector of length n where each element represents the outcome of the ith
#' pairwise comparison
#' @export codeOVElement

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


classFromOV <- function(ov) {
  if (is.character(ov) | is.factor(ov)) {
    if (ov %in% pos_ov_s & ov %in% neg_ov_s) {
        return("Ambig")
    }

    if (ov %in% pos_ov_s & (!(ov %in% neg_ov_s))) {
        return("Positive")
      }

    if (!(ov %in% pos_ov_s) & ov %in% neg_ov_s) {
      return("Negative")
    }

   if (ov %in% null_ov_s) {
     return("Null")
   }

    if (!(ov %in% pos_ov_s) & !(ov %in% neg_ov_s)) {
      return("Anomaly")
    }

    else {
      return("NI")
    }
  }


  if (is.vector(ov)) {
  # Ambiguous
    if (any(apply(pos_ov, 1, function(x) (all(v == x))) == TRUE)) {
      if(any(apply(neg_ov, 1, function(x) (all(v == x))) == TRUE)) {
        return("Ambig")
      }
    }
    # Postive
    if (any(apply(pos_ov, 1, function(x) (all(v == x))) == TRUE)) {
      if(any(apply(neg_ov, 1, function(x) (all(v == x))) == FALSE)) {
        return("Pos")
      }
    }
    # Negative
    if (any(apply(pos_ov, 1, function(x) (all(v == x))) == FALSE)) {
      if(any(apply(neg_ov, 1, function(x) (all(v == x))) == TRUE)) {
        return("Neg")
      }
    }
    else {
      return("Anomaly")
    }
  }

}

#' Enumerate theoretical pairwise comparison outcome vectors given n conditions and k ranks
#'
#' @return matrix where each row represents a pairwise comparison outcome vector
#' @export enumerateThOutcomeVectors
enumerateThOutcomeVectors <- function() {
  library(magrittr)
  pwt_n <- gtools::permutations(4, 4, repeats.allowed=TRUE) %>%
    rankTuples() %>% unique() %>%
    pwCompare()
  # pairwise comparison outcome vectors as character vector
  pwt_s <<- apply(pwt_n, 1, stringCoerce)
  return(pwt_n)
}

#' Generate a list where each element is a matrix/vector of outcome vectors associated with a
#' known mode/class
#'
#' @return named list where each element contains a vector or matrix with the
#' outcome vector(s) associated with a particular class or mode
#' @export outcomeVectorsByMode
#' @examples
#' library(magrittr)
#' pw <- outcomeVectorsByMode()

outcomeVectorsByMode <- function() {
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
  oppositeVectors <- function(a_matrix) {
    reversed <- apply(a_matrix, 1, function(x) (lapply(x, reverseSign) %>% unlist))
    neg_vecs <- matrix(NaN, nrow = nrow(a_matrix), ncol = 6)
    for (i in 1:ncol(reversed)) {
      neg_vecs[i,] <- reversed[,i]
    }
    return(neg_vecs)
  }

  low_stab <- matrix(c(-1, -1, -1, -1, -1, -1,
                       -1, -1, -1, -1, -1, 0,
                       -1, -1, -1, -1, -1, 1,
                       -1, -1, -1, -1, 0, 1,
                       -1, -1, -1, -1, 1, 1,
                       -1, -1, -1, 0, -1, -1,
                       -1, -1, -1, 0, 0, 0,
                       -1, -1, -1, 0, 1, 1,
                       -1, -1, -1, 1, -1, -1,
                       -1, -1, -1, 1, 0, 1,
                       -1, -1, -1, 1, 1, -1,
                       -1, -1, -1, 1, 1, 0,
                       -1, -1, -1, 1, 1, 1
  ),
  nrow = 13, ncol = 6, byrow = TRUE
  )

  high_stab <- oppositeVectors(low_stab)

  emer_pos_syn <- matrix(c(-1, -1, 0, -1, 1, 1,
                           -1, -1, 0, 0, 1, 1,
                           -1, -1, 0, 1, 1, 1,
                           -1, -1, 1, -1, 1, 1,
                           -1, -1, 1, 0, 1, 1,
                           -1, -1, 1, 1, 1, 1,
                           -1, 0, 1, 1, 1, 1,
                           0, -1, 1, -1, 1, 1,
                           0, 0, 1, 0, 1, 1),
                         nrow = 9, ncol = 6, byrow = TRUE)
  emer_neg_syn <- oppositeVectors(emer_pos_syn)

  y_restores_x <- matrix(c(-1, 0, -1, 1, 1, -1,
                           -1, 0, 0, 1, 1, 0,
                           -1, 1, -1, 1, 1, -1,
                           -1, 1, 0, 1, 1, -1,
                           -1, 1, 1, 1, 1, -1,
                           -1, 1, 1, 1, 1, 0),
                         nrow = 6, ncol = 6, byrow = TRUE)

  y_inhibits_x <- oppositeVectors(y_restores_x)

  x_restores_y <- matrix(c(0, -1, -1, -1, -1, 1,
                           0, -1, 0, -1, 0, 1,
                           1, -1, -1, -1, -1, 1,
                           1, -1, 0, -1, -1, 1,
                           1, -1, 1, -1, -1, 1,
                           1, -1, 1, -1, 0, 1),
                         nrow = 6, ncol = 6, byrow = TRUE)

  x_inhibits_y <- oppositeVectors(x_restores_y)

  pos_syn <- matrix(c(-1, 1, 1, 1, 1, 1,
                      0, 1, 1, 1, 1, 1,
                      1, -1, 1, -1, 1, 1,
                      1, 0, 1, -1, 1, 1,
                      1, 1, 1, -1, 1, 1,
                      1, 1, 1, 0, 1, 1,
                      1, 1, 1, 1, 1, 1
  ), nrow = 7, ncol = 6, byrow = TRUE)

  ni_ov <- matrix(c(rep(0, 6),
                    1, 0, 1, -1, 0, 1,
                    0, 1, 1, 1, 1, 0,
                    -1, 0, -1, 1, 0, -1,
                    0, -1, -1, -1, -1, 0),
                  nrow = 5, ncol = 6, byrow = TRUE)

  sym_right <- c(1, 0, 1, -1, 0, 1)
  sym_left <- c(-1, 0, -1, 1, 0, -1)

  step_up <- c(0, 1, 1, 1, 1, 0)
  step_down <- c(0, -1, -1, -1, -1, 0)
  neg_syn <- oppositeVectors(pos_syn)
  pos_ov <- rbind(low_stab, x_restores_y, y_restores_x, pos_syn, emer_pos_syn)
  neg_ov <- rbind(high_stab, x_inhibits_y, y_inhibits_x, neg_syn, emer_neg_syn)
  all_ov <- rbind(pos_ov, neg_ov)
  ov_list <- list(low_stab = low_stab,
                  x_restores_y = x_restores_y,
                  y_inhibits_x = y_inhibits_x,
                  pos_syn = pos_syn,
                  emer_pos_syn = emer_pos_syn,
                  high_stab = high_stab,
                  x_inhibits_y = x_inhibits_y,
                  y_inhibits_x = y_inhibits_x,
                  neg_syn = neg_syn,
                  emer_neg_syn = emer_neg_syn,
                  pos_ov = pos_ov,
                  neg_ov = neg_ov,
                  all_ov = all_ov,
                  step_up = step_up,
                  step_down = step_down,
                  sym_left = sym_left,
                  sym_right = sym_right)
  return(ov_list)
}




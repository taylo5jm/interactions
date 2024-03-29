% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Interaction.R
\docType{class}
\name{Interaction-class}
\alias{Interaction-class}
\title{Interaction class to hold results from mode classification with the edgeR
limma pipelines}
\description{
Interaction class to hold results from mode classification with the edgeR
limma pipelines
}
\section{Slots}{

 \describe{
   \item{\code{interaction.tbl}:}{Data.frame containing mode classification
   and other data from edgeR or limma}
   \item{\code{modes}:}{Factor vector of length n where each element is an
   interaction mode classification that can mapped to a gene.}
   \item{\code{outcome.vectors}:}{Matrix of dimensions i x 6 where each row
   represents a 6-dimensional outcome vector corresponding to the possible
   comparisons when gene expression outputs from 4 conditions are compared}
   \item{\code{outcome.vectors.s}:}{Character vector of length n where each
   element is an outcome vector as a string}
   \item{\code{voom}:}{EList object}
   \item{\code{d}:}{DGEList object}
   \item{\code{design}:}{list with design matrix and contrasts from
   makeContrasts()}
   }
}


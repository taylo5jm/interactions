% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/edgeR_generics.R
\name{codeEdgeROutcomeVectors}
\alias{codeEdgeROutcomeVectors}
\title{Code pairwise comparison outcome vectors with edgeR likelilood ratio tests p-values
and log fold changes}
\usage{
codeEdgeROutcomeVectors(pairwise, alpha)
}
\arguments{
\item{pairwise}{list of results from glmFitAndLRT()}

\item{alpha}{numeric vector of length n = 1 indicating alpha}
}
\value{
res matrix where res[,1:6] are the elements of the pairwise comparison outcome vectors
and res[,7:12] are the log-fold-change ratios from the 6 contrasts
}
\description{
Code pairwise comparison outcome vectors with edgeR likelilood ratio tests p-values
and log fold changes
}
\examples{
library(edgeR)
library(limma)
}


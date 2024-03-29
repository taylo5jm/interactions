% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/interaction_tools.R
\name{reportnNI}
\alias{reportnNI}
\title{Report number of null interaction vectors coded}
\usage{
reportnNI(modes)
}
\arguments{
\item{character}{vector of length n indicating the interaction mode classifications
of each gene in v$genes}
}
\value{
number of null outcome vectors (0, 0, 0, 0, 0, 0) coded
}
\description{
Report number of null interaction vectors coded
}
\examples{
\dontrun{
reportnNI(limma_modes)
}
}


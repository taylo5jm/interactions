% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/barplots.R
\name{colorBarplot}
\alias{colorBarplot}
\title{Color bars in geom_bar plot given an interaction class and ggplot2}
\usage{
colorBarplot(int_class, gg_plot)
}
\arguments{
\item{int_class}{character vector of length n = 1 indicating interaction class ("Positive", "Negative", "NI",
"Null", "Ambiguous)}

\item{gg_plot}{ggplot2 object}
}
\value{
ggplot2 object with bars colored with respect to interaction class
}
\description{
Color bars in geom_bar plot given an interaction class and ggplot2
}


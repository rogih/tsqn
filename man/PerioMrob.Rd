% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/per_rob_M.r
\name{PerioMrob}
\alias{PerioMrob}
\title{Robust M-periodogram}
\usage{
PerioMrob(series)
}
\arguments{
\item{series}{univariate time series}
}
\value{
a numeric vector containing the robust estimates of the spectral density
}
\description{
This function computes the robust M-periodogram proposed by Reisen et al. (2017).
}
\examples{
PerioMrob(ldeaths)
}
\references{
Reisen, V. A. and Lévy-Leduc, C. and Taqqu, M. (2017) An M-estimator for the long-memory parameter. \emph{To appear in Journal of Statistical Planning and Inference}.

Geweke, J. and Porter-Hudak, S. (1983) The estimation and application of long memory time series models. \emph{Journal of Time Series Analysis}, \bold{4}, 221--238.
}
\author{
Valderio Reisen, Céline Lévy-Leduc and Higor Cotta.
}

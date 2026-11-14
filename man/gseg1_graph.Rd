\name{gseg1_graph}
\alias{gseg1_graph}
\title{Single-Change Scan with a Supplied Similarity Graph}
\description{Applies the single-change graph scan to a supplied edge matrix.}
\usage{
gseg1_graph(n, E, statistics = "m", n0 = 0.05 * n,
  n1 = 0.95 * n, pval.appr = TRUE, skew.corr = TRUE,
  pval.perm = FALSE, B = 100)
}
\arguments{
\item{n}{Number of observations in the sequence.}
\item{E}{Two-column edge matrix for the similarity graph.}
\item{statistics}{Scan statistic or statistics; see \code{\link{gseg1}}.}
\item{n0}{First candidate change-point location.}
\item{n1}{Last candidate change-point location.}
\item{pval.appr}{Whether to calculate the analytic p-value approximation.}
\item{skew.corr}{Whether to include its skewness correction.}
\item{pval.perm}{Whether to calculate a permutation p-value.}
\item{B}{Number of permutations.}
}
\value{A list containing scan statistics and requested p-values.}
\seealso{\code{\link{gseg1}}, \code{\link{gseg_kmst}}}
\examples{
set.seed(1)
x <- matrix(rnorm(80), nrow = 20)
E <- gseg_kmst(x)
ans <- gseg1_graph(nrow(x), E, pval.appr = FALSE)
}

\details{
The requested statistic must have positive permutation variance at each
candidate location. Graphs for which a candidate score is undefined are
rejected rather than silently omitted. For a single candidate location,
set \code{pval.appr = FALSE} and use permutation calibration instead of
the continuous-range analytic scan approximation.
}

\note{
Analytic integrations are checked for convergence and finite probability
outputs in \code{[0,1]}. An unavailable approximation raises an error rather
than returning an invalid probability. Setting \code{pval.appr = FALSE} allows
justified permutation calibration to be used instead. Scan components must be
well-defined on the requested domain; undefined components are not replaced
by zero.
}

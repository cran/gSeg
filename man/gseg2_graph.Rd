\name{gseg2_graph}
\alias{gseg2_graph}
\title{Changed-Interval Scan with a Supplied Similarity Graph}
\description{Applies the changed-interval graph scan to a supplied edge matrix.}
\usage{
gseg2_graph(n, E, statistics = "m", l0 = 0.05 * n,
  l1 = 0.95 * n, pval.appr = TRUE, skew.corr = TRUE,
  pval.perm = FALSE, B = 100)
}
\arguments{
\item{n}{Number of observations in the sequence.}
\item{E}{Two-column edge matrix for the similarity graph.}
\item{statistics}{Scan statistic or statistics; see \code{\link{gseg2}}.}
\item{l0}{Minimum candidate interval length.}
\item{l1}{Maximum candidate interval length.}
\item{pval.appr}{Whether to calculate the analytic p-value approximation.}
\item{skew.corr}{Whether to include its skewness correction.}
\item{pval.perm}{Whether to calculate a permutation p-value.}
\item{B}{Number of permutations.}
}
\value{A list containing scan statistics and requested p-values.}
\seealso{\code{\link{gseg2}}, \code{\link{gseg_kmst}}}
\examples{
set.seed(1)
x <- matrix(rnorm(80), nrow = 20)
E <- gseg_kmst(x)
ans <- gseg2_graph(nrow(x), E, pval.appr = FALSE)
}

\details{
An interval with boundaries \code{(i, j)} contains observations
\code{i + 1, ..., j}. When there are \code{r} tied maximizing intervals,
\code{tauhat} (or \code{tauhat_a}/\code{tauhat_u}) retains the earlier
layout: its first \code{r} entries are the left boundaries and its last
\code{r} entries are the corresponding right boundaries. All returned
intervals satisfy the requested length range after the documented endpoint
adjustments. The analytic approximation requires a nonzero length range;
use permutation calibration for a single permitted interval length.
}

\note{
Analytic integrations are checked for convergence and finite probability
outputs in \code{[0,1]}. An unavailable approximation raises an error rather
than returning an invalid probability. Setting \code{pval.appr = FALSE} allows
justified permutation calibration to be used instead. Scan components must be
well-defined on the requested domain; undefined components are not replaced
by zero.
}

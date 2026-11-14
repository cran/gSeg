\name{gseg2}
\alias{gseg2}
\title{Graph-Based Detection of a Changed Interval}
\description{
Detects an interval whose distribution differs from the remainder of a
sequence. Given observations or dissimilarities, the function constructs a
k-MST and applies the max-type edge-count scan by default. A user-supplied
similarity graph is also supported.
}
\usage{
gseg2(x = NULL, dissimilarity = NULL, statistics = "m",
  l0 = NULL, l1 = NULL, pval.appr = TRUE, skew.corr = TRUE,
  pval.perm = FALSE, B = 100, k = NULL, n = NULL, E = NULL)
}
\arguments{
\item{x}{A numeric matrix with observations in rows. For compatibility, the
first positional argument may instead be the number of observations when the
second positional argument is an edge matrix.}
\item{dissimilarity}{A square symmetric dissimilarity matrix or a \code{dist}
object. For compatibility, the second positional argument may instead be a
two-column edge matrix. Supply exactly one of \code{x} and
\code{dissimilarity} for a data-level call.}
\item{statistics}{The scan statistic or statistics to compute. The default
\code{"m"} reports only the max-type edge-count test (MET). Use \code{"g"}
for GET, \code{"o"} for the original edge-count test, \code{"w"} for the
weighted edge-count test, or \code{"all"} for all four. Long-form names are
also accepted.}
\item{l0}{Minimum candidate interval length. The default is \code{0.05 * N}.}
\item{l1}{Maximum candidate interval length. The default is \code{0.95 * N}.}
\item{pval.appr}{Whether to calculate the analytic p-value approximation.
This approximation requires at least nine observations.}
\item{skew.corr}{Whether to include the skewness correction in the analytic
p-value approximation.}
\item{pval.perm}{Whether to calculate a permutation p-value.}
\item{B}{Number of permutations.}
\item{k}{Number of MST layers for a data-level call. The default is
\code{floor(sqrt(N))}; the largest allowed value is \code{floor(N/2)}.}
\item{n}{Number of observations in a graph-level call. Both named and
positional legacy graph arguments are supported.}
\item{E}{Two-column edge matrix in a graph-level call. Do not supply
\code{k} together with a supplied graph.}
}
\value{
A list with \code{scanZ} and any requested p-values. Each requested component
of \code{scanZ} contains the estimated interval \code{tauhat}, maximum scan
statistic \code{Zmax}, and scan surface. Data-level calls additionally return
the constructed \code{graph} and \code{k}.
}
\seealso{
\code{\link{gseg1}}, \code{\link{gseg2_graph}},
\code{\link{gseg2_discrete}}, \code{\link{gseg_kmst}}
}
\examples{
set.seed(2)
x <- matrix(rnorm(80), nrow = 20)
x[8:13, ] <- x[8:13, ] + 1

# Main interface: build the default k-MST and report MET.
ans <- gseg2(x, pval.appr = FALSE)
names(ans$scanZ)
ans$k

# Earlier graph-level calls remain supported.
ans.graph <- gseg2(nrow(x), ans$graph, pval.appr = FALSE)
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

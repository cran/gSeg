\name{gseg1}
\alias{gseg1}
\title{Graph-Based Detection of a Single Change-Point}
\description{
Detects a single distributional change in a sequence. Given observations or
dissimilarities, the function constructs a k-MST and applies the max-type
edge-count scan by default. A user-supplied similarity graph is also supported.
}
\usage{
gseg1(x = NULL, dissimilarity = NULL, statistics = "m",
  n0 = NULL, n1 = NULL, pval.appr = TRUE, skew.corr = TRUE,
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
for the generalized edge-count test (GET), \code{"o"} for the original
edge-count test, \code{"w"} for the weighted edge-count test, or
\code{"all"} for all four. Long-form names are also accepted.}
\item{n0}{First candidate change-point location. The default is
\code{0.05 * N}.}
\item{n1}{Last candidate change-point location. The default is
\code{0.95 * N}.}
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
of \code{scanZ} contains the estimated change-point \code{tauhat}, maximum
scan statistic \code{Zmax}, and scan path. Data-level calls additionally
return the constructed \code{graph} and \code{k}.
}
\seealso{
\code{\link{gseg2}}, \code{\link{gseg1_graph}},
\code{\link{gseg1_discrete}}, \code{\link{gseg_kmst}}
}
\examples{
set.seed(1)
x <- rbind(matrix(rnorm(40), 10), matrix(rnorm(40, 1), 10))

# Main interface: build the default k-MST and report MET.
ans <- gseg1(x, pval.appr = FALSE)
names(ans$scanZ)
ans$k

# A dissimilarity matrix may be supplied instead.
ans.d <- gseg1(dissimilarity = as.matrix(dist(x)), pval.appr = FALSE)

# Request GET explicitly.
ans.get <- gseg1(x, statistics = "g", pval.appr = FALSE)

# Earlier graph-level calls remain supported.
ans.graph <- gseg1(nrow(x), ans$graph, pval.appr = FALSE)
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

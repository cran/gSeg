\name{gseg2}
\alias{gseg2}
\title{Graph-Based Change-Point Detection for Changed Interval}
\description{This function finds an interval in the sequence where their underlying distribution differs from the rest of the sequence. It provides four graph-based test statistics.
}
\usage{
gseg2(n, E, statistics=c("all","o","w","g","m"), l0=0.05*n, l1=0.95*n, pval.appr=TRUE,
 skew.corr=TRUE, pval.perm=FALSE, B=100)
}
\arguments{
  \item{n}{The number of observations in the sequence.}
  \item{E}{The edge matrix (a "number of edges" by 2 matrix) for the similarity graph.  Each row contains the node indices of an edge.}
   \item{statistics}{The scan statistic to be computed. A character indicating the type of of scan statistic desired. The default is \code{"all"}.

   \code{"all"}: specifies to compute \bold{a}ll of the scan statistics: original, weighted, generalized, and max-type;

   \code{"o", "ori"} or \code{"original"}: specifies the \bold{o}riginal edge-count scan statistic;

   \code{"w"} or \code{"weighted"}: specifies the \bold{w}eighted edge-count scan statistic;

   \code{"g"} or \code{"generalized"}: specifies the \bold{g}eneralized edge-count scan statistic; and

   \code{"m"} or \code{"max"}: specifies the \bold{m}ax-type edge-count scan statistic.

   }
  \item{l0}{The minimum length of the interval to be considered as a changed interval.}
  \item{l1}{The maximum length of the interval to be considered as a changed interval.}
  \item{pval.appr}{If it is TRUE, the function outputs p-value approximation based on asymptotic properties.}
  \item{skew.corr}{This argument is useful only when pval.appr=TRUE.  If skew.corr is TRUE, the p-value approximation would incorporate skewness correction.}
  \item{pval.perm}{If it is TRUE, the function outputs p-value from doing B permutations, where B is another argument that you can specify. Doing permutation could be time consuming, so use this argument with caution as it may take a long time to finish the permutation.}
  \item{B}{This argument is useful only when pval.perm=TRUE.  The default value for B is 100.}
}
\value{
Returns a list \code{scanZ} with \code{tauhat}, \code{Zmax}, and a matrix of the scan statistics for each type of scan statistic specified. See below for more details.
  \item{tauhat}{An estimate of the two ends of the changed interval.}
  \item{Zmax}{The test statistic (maximum of the scan statistics).}
  \item{Z}{A matrix of the original scan statistics (standardized counts) if statistic specified is "all" or "o".}
  \item{Zw}{A matrix of the weighted scan statistics (standardized counts) if statistic specified is "all" or "w".}
  \item{S}{A matrix of the generalized scan statistics (standardized counts) if statistic specified is "all" or "g".}
  \item{M}{A matrix of the max-type scan statistics (standardized counts) if statistic specified is "all" or "m".}
  \item{R}{A matrix of raw counts of the original scan statistic. This output only exists if the statistic specified is "all" or "o".}
  \item{Rw}{A matrix of raw counts of the weighted scan statistic. This output only exists if statistic specified is "all" or "w".}
  \item{pval.appr}{The approximated p-value based on asymptotic theory for each type of statistic specified.}
  \item{pval.perm}{This output exists only when the argument pval.perm is TRUE . It is the permutation p-value from B permutations and appears for each type of statistic specified (same for perm.curve, perm.maxZs, and perm.Z).}
  \item{perm.curve}{A B by 2 matrix with the first column being critical values corresponding to the p-values in the second column.}
  \item{perm.maxZs}{A sorted vector recording the test statistics in the B permutaitons.}
  \item{perm.Z}{A B by n-squared matrix with each row being the vectorized scan statistics from each permutaiton run.}
}
\seealso{
  \code{\link{gSeg}}, \code{\link{gseg1}}, \code{\link{gseg2_discrete}}
}
\examples{
data(Example)
# Five examples, each example is a n-length sequence.
# Ei (i=1,...,5): an edge matrix representing a similarity graph constructed on the
# observations in the ith sequence.  
# Check '?gSeg' to see how the Ei's were constructed.

## E5 is an edge matrix representing a similarity graph.
# It is constructed on a sequence of length n=200 with a change in both mean
# and variance on an interval (tau1 = 155, tau2 = 185).
r5=gseg2(n,E5,statistics="all")

}


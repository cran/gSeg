\name{gseg1}
\alias{gseg1}
\title{Graph-Based Change-Point Detection for Single Change-Point}
\description{This function finds a break point in the sequence where the underlying distribution changes. It provides four graph-based test statistics.}
\usage{
gseg1(n, E, statistics=c("all","o","w","g","m"), n0=0.05*n, n1=0.95*n, pval.appr=TRUE,
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
  \item{n0}{The starting index to be considered as a candidate for the change-point.}
  \item{n1}{The ending index to be considered as a candidate for the change-point.}
  \item{pval.appr}{If it is TRUE, the function outputs p-value approximation based on asymptotic properties.}
  \item{skew.corr}{This argument is useful only when pval.appr=TRUE.  If skew.corr is TRUE, the p-value approximation would incorporate skewness correction.}
  \item{pval.perm}{If it is TRUE, the function outputs p-value from doing B permutations, where B is another argument that you can specify. Doing permutation could be time consuming, so use this argument with caution as it may take a long time to finish the permutation.}
  \item{B}{This argument is useful only when pval.perm=TRUE.  The default value for B is 100.}
}

\value{
  Returns a list \code{scanZ} with \code{tauhat}, \code{Zmax}, and a vector of the scan statistics for each type of scan statistic specified. See below for more details.
  \item{tauhat}{An estimate of the location of the change-point.}
  \item{Zmax}{The test statistic (maximum of the scan statistics).}
  \item{Z}{A vector of the original scan statistics (standardized counts) if statistic specified is "all" or "o".}
  \item{Zw}{A vector of the weighted scan statistics (standardized counts) if statistic specified is "all" or "w".}
  \item{S}{A vector of the generalized scan statistics (standardized counts) if statistic specified is "all" or "g".}
  \item{M}{A vector of the max-type scan statistics (standardized counts) if statistic specified is "all" or "m".}
  \item{R}{A vector of raw counts of the original scan statistic. This output only exists if the statistic specified is "all" or "o".}
  \item{Rw}{A vector of raw counts of the weighted scan statistic. This output only exists if statistic specified is "all" or "w".}
  \item{pval.appr}{The approximated p-value based on asymptotic theory for each type of statistic specified.}
  \item{pval.perm}{This output exists only when the argument pval.perm is TRUE . It is the permutation p-value from B permutations and appears for each type of statistic specified (same for perm.curve, perm.maxZs, and perm.Z).}
  \item{perm.curve}{A B by 2 matrix with the first column being critical values corresponding to the p-values in the second column.}
  \item{perm.maxZs}{A sorted vector recording the test statistics in the B permutaitons.}
  \item{perm.Z}{A B by n matrix with each row being the scan statistics from each permutaiton run.}
}
\seealso{
  \code{\link{gSeg}}, \code{\link{gseg2}}, \code{\link{gseg1_discrete}}
}
\examples{

data(Example)
# Five examples, each example is a n-length sequence.
# Ei (i=1,...,5): an edge matrix representing a similarity graph constructed on the
# observations in the ith sequence.  
# Check '?gSeg' to see how the Ei's were constructed.


## E1 is an edge matrix representing a similarity graph.
# It is constructed on a sequence of length n=200 with a change in mean 
# in the middle of the sequence (tau = 100).

r1 = gseg1(n,E1, statistics="all")  
# output results based on all four statistics
# the scan statistics can be found in r1$scanZ
r1_a = gseg1(n,E1, statistics="w")  
# output results based on the weighted edge-count statistic
r1_b = gseg1(n,E1, statistics=c("w","g"))  
# output results based on the weighted edge-count statistic 
# and generalized edge-count statistic


## E2 is an edge matrix representing a similarity graph.
# It is constructed on a sequence of length n=200 with a change in mean 
# away from the middle of the sequence (tau=45).
r2 = gseg1(n,E2,statistic="all")


## E3 is an edge matrix representing a similarity graph.
# It is constructed on a sequence of length n=200 with a change in both mean 
# and variance away from the middle of the sequence (tau = 145).
r3=gseg1(n,E3,statistic="all")


## E4 is an edge matrix representing a similarity graph.
# It is constructed on a sequence of length n=200 with a change in both mean 
# and variance away from the middle of the sequence (tau = 50).
r4=gseg1(n,E4,statistic="all")

}


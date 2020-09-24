\name{gseg2_discrete}
\alias{gseg2_discrete}
\title{Graph-Based Change-Point Detection for Changed Interval for Data with Repeated Observations}
\description{This function finds an interval in the sequence where their underlying distribution differs from the rest of the sequence when data has repeated observations. It provides four graph-based test statistics.}
\usage{
gseg2_discrete(n, E, id, statistics=c("all","o","w","g","m"), l0=0.05*n, l1=0.95*n, 
   pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100)
}
\arguments{
  \item{n}{The number of observations in the sequence.}
  \item{E}{The edge matrix (a "number of edges" by 2 matrix) for the similarity graph. Each row contains the node indices of an edge.}
  \item{id}{The index of observations (order of observations).}
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
  Returns a list \code{scanZ} with \code{tauhat}, \code{Zmax}, and a vector of the scan statistics for each type of scan statistic specified. See below for more details.
  \item{tauhat_a}{An estimate of the two ends of the changed interval for averaging approach.}
  \item{tauhat_u}{An estimate of the two ends of the changed interval for union approach.}
  \item{Z_a_max}{The test statistic (maximum of the scan statistics) for averaging approach.}
  \item{Z_u_max}{The test statistic (maximum of the scan statistics) for union approach.}
  \item{Zo_a}{A matrix of the original scan statistics (standardized counts) for averaging approach if statistic specified is "all" or "o".}
  \item{Zo_u}{A matrix of the original scan statistics (standardized counts) for union approach if statistic specified is "all" or "o".}
  \item{Zw_a}{A matrix of the weighted scan statistics (standardized counts) for averaging approach if statistic specified is "all" or "w".}
  \item{Zw_u}{A matrix of the weighted scan statistics (standardized counts) for union approach if statistic specified is "all" or "w".}
  \item{S_a}{A matrix of the generalized scan statistics (standardized counts) for averaging approach if statistic specified is "all" or "g".}
  \item{S_u}{A matrix of the generalized scan statistics (standardized counts) for union approach if statistic specified is "all" or "g".}
  \item{M_a}{A matrix of the max-type scan statistics (standardized counts) for averaging approach if statistic specified is "all" or "m".}
  \item{M_u}{A matrix of the max-type scan statistics (standardized counts) for union approach if statistic specified is "all" or "m".}
  \item{Ro_a}{A matrix of raw counts of the original scan statistic for averaging approach. This output only exists if the statistic specified is "all" or "o".}
  \item{Ro_u}{A matrix of raw counts of the original scan statistic for union approach. This output only exists if the statistic specified is "all" or "o".}
  \item{Rw_a}{A matrix of raw counts of the weighted scan statistic for averaging approach. This output only exists if statistic specified is "all" or "w".}
  \item{Rw_a}{A matrix of raw counts of the weighted scan statistic for union approach. This output only exists if statistic specified is "all" or "w".}
  \item{pval.appr}{The approximated p-value based on asymptotic theory for each type of statistic specified.}
  \item{pval.perm}{This output exists only when the argument pval.perm is TRUE . It is the permutation p-value from B permutations and appears for each type of statistic specified (same for perm.curve, perm.maxZs, and perm.Z).}
  \item{perm.curve}{A B by 2 matrix with the first column being critical values corresponding to the p-values in the second column.}
  \item{perm.maxZs}{A sorted vector recording the test statistics in the B permutaitons.}
  \item{perm.Z}{A B by n-squared matrix with each row being the vectorized scan statistics from each permutaiton run.}
}
\seealso{
  \code{\link{gSeg}}, \code{\link{gseg2}}, \code{\link{nnl}}, \code{\link{gseg1_discrete}}
}
\examples{
d = 50
mu = 2
tau = 100
n = 200

set.seed(500)
y1_temp = matrix(rnorm(d*tau),tau)
sam1 = sample(1:tau, replace = TRUE)
y1 = y1_temp[sam1,] 
y2_temp = matrix(rnorm(d*(n-tau),mu/sqrt(d)), n-tau)
sam2 = sample(1:tau, replace = TRUE)
y2 = y2_temp[sam2,] 

y = rbind(y1, y2)

# This data y has repeated observations
y_uni = unique(y)
E = nnl(dist(y_uni), 1)

cha = do.call(paste, as.data.frame(y))    
id = match(cha, unique(cha))

r1 = gseg2_discrete(n, E, id, statistics="all")
}
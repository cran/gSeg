\name{gseg2}
\alias{gseg2}
\title{Graph-Based Change-Point Detection for Changed Interval}
\description{This function finds an interval in the sequence where their underlying distribution differs from the rest of the sequence.
}
\usage{
gseg2(n, E, l0=0.05*n, l1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100)
}
\arguments{
  \item{n}{The number of observations in the sequence.}
  \item{E}{The edge matrix (a "number of edges" by 2 matrix) for the similarity graph.  Each row contains the node indices of an edge.}
  \item{l0}{The minimum length of the interval to be considered as a changed interval.}
  \item{l1}{The maximum length of the interval to be considered as a changed interval.}
  \item{pval.appr}{If it is TRUE, the function outputs p-value approximation based on asymptotic properties.}
  \item{skew.corr}{This argument is useful only when pval.appr=TURE.  If skew.corr is TRUE, the p-value approximation would incorporate skewness correction.}
  \item{pval.perm}{If it is TRUE, the function outputs p-value from doing B permutations, where B is another argument that you can specify.  Use this argument with caution that it may take long time to finish the permutation.}
  \item{B}{This argument is useful only when pval.perm=TRUE.  The default value for B is 100.}
}
\value{
  \item{tauhat}{An estimate of the two ends of the changed interval.}
  \item{Zmax}{The test statistic (maximum of the scan statistics).}
  \item{Z}{A matrix of scan statistics (standardized counts).}
  \item{R}{A matrix of raw counts.}
  \item{pval.appr}{This output exists only when the argument pval.appr is TRUE.  The approximated p-value.}
  \item{pval.perm}{This output exists only when the argument pval.perm is TRUE (same for perm.curve, perm.maxZs, and perm.Z).  p-value calculated from B permutations.}
  \item{perm.curve}{A B by 2 matrix with the first column being critical values corresponding to the p-values in the second column.}
  \item{perm.maxZs}{A sorted vector recording the test statistics in the B permutaitons.}
  \item{perm.Z}{A B by n-squared matrix with each row being the vectorized scan statistics from each permutaiton run.}
}
\seealso{
  \code{\link{gseg1}}
}
\examples{
data(Example) # this example data has: n (the number of observations) and E (an edge matrix)
r2 = gseg2(n,E)
}


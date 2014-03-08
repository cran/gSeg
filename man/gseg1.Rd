\name{gseg1}
\alias{gseg1}
\title{Graph-Based Change-Point Detection for Single Change-Point}
\description{This function finds a break point in the sequence where the underlying distribution changes.
}
\usage{
gseg1(n, E, n0=0.05*n, n1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100)
}
\arguments{
  \item{n}{The number of observations in the sequence.}
  \item{E}{The edge matrix (a "number of edges" by 2 matrix) for the similarity graph.  Each row contains the node indices of an edge.}
  \item{n0}{The starting index to be considered as a candidate for the change-point.}
  \item{n1}{The ending index to be considered as a candidate for the change-point.}
  \item{pval.appr}{If it is TRUE, the function outputs p-value approximation based on asymptotic properties.}
  \item{skew.corr}{This argument is useful only when pval.appr=TURE.  If skew.corr is TRUE, the p-value approximation would incorporate skewness correction.}
  \item{pval.perm}{If it is TRUE, the function outputs p-value from doing B permutations, where B is another argument that you can specify.  Use this argument with caution that it may take long time to finish the permutation.}
  \item{B}{This argument is useful only when pval.perm=TRUE.  The default value for B is 100.}
}
\value{
  \item{tauhat}{An estimate of the location of the change-point.}
  \item{Zmax}{The test statistic (maximum of the scan statistics).}
  \item{Z}{A vector of scan statistics (standardized counts).}
  \item{R}{A vector of raw counts.}
  \item{pval.appr}{This output exists only when the argument pval.appr is TRUE.  The approximated p-value.}
  \item{pval.perm}{This output exists only when the argument pval.perm is TRUE (same for perm.curve, perm.maxZs, and perm.Z).  p-value calculated from B permutations.}
  \item{perm.curve}{A B by 2 matrix with the first column being critical values corresponding to the p-values in the second column.}
  \item{perm.maxZs}{A sorted vector recording the test statistics in the B permutaitons.}
  \item{perm.Z}{A B by n matrix with each row being the scan statistics from each permutaiton run.}
}
\seealso{
  \code{\link{gseg2}}
}
\examples{
data(Example) # this example data has: n (the number of observations) and E (an edge matrix)
r1 = gseg1(n,E)
}


\name{gSeg}
\alias{gseg}
\title{Graph-Based Change-Point Detection}
\description{This package can be used to estimate change-points in a sequence of observations, where the observation can be a vector or a data object, e.g., a network.  A similarity graph is required.  It can be a minimum spanning tree, a minimum distance pairing, a nearest neighbor graph, or a graph based on domain knowledge.  

If you believe the sequence has at most one change point, the function \code{\link{gseg1}} should be used.  If you believe an interval of the sequence has a changed distribution, the function \code{\link{gseg2}} should be used.  If you feel the sequence has multiple change-points, you can use \code{\link{gseg1}} and \code{\link{gseg2}} multiple times. See \code{\link{gseg1}} and \code{\link{gseg2}} for the details of these two function.
}
\author{
Hao Chen and Nancy R. Zhang

Maintainer: Hao Chen (hxchen@ucdavis.edu)
}
\references{
Chen, H. and Zhang, N.R. (2014) Graph-Based Change-Point Detection.

Online access: http://arxiv.org/abs/1209.1625
}

\seealso{
  \code{\link{gseg1}}, \code{\link{gseg2}}
}
\examples{
data(Example) # this example data has: n (the number of observations) and E (an edge matrix)
r1 = gseg1(n,E)
r2 = gseg2(n,E)
}


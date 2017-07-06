\name{gSeg}
\alias{gSeg}
\title{Graph-Based Change-Point Detection}
\description{This package can be used to estimate change-points in a sequence of observations, where the observation can be a vector or a data object, e.g., a network.  A similarity graph is required.  It can be a minimum spanning tree, a minimum distance pairing, a nearest neighbor graph, or a graph based on domain knowledge.

If you believe the sequence has at most one change point, the function \code{\link{gseg1}} should be used.  If you believe an interval of the sequence has a changed distribution, the function \code{\link{gseg2}} should be used.  If you feel the sequence has multiple change-points, you can use \code{\link{gseg1}} and \code{\link{gseg2}} multiple times. See \code{\link{gseg1}} and \code{\link{gseg2}} for the details of these two function.
}
\author{
Hao Chen, Nancy R. Zhang, and Lynna Chu

Maintainer: Hao Chen (hxchen@ucdavis.edu)
}
\references{
	Chen, Hao, and Nancy Zhang. (2015). Graph-based change-point detection. The Annals of Statistics, 43(1), 139-176.
	
	Chu, Lynna, and Hao Chen. (2017). Asymptotic distribution-free change-point detection for modern data. arXiv:1707.00167
}

\seealso{
  \code{\link{gseg1}}, \code{\link{gseg2}}
}
\examples{

data(Example)
# Five examples, each example is a n-length sequence.
# Ei (i=1,...,5): an edge matrix representing a similarity graph constructed on the
# observations in the ith sequence.  
# The following code shows how the Ei's were constructed.

require(ade4) 
# For illustration, we use 'mstree' in this package to construct the similarity graph.
# You can use other ways to construct the graph.

## Sequence 1: change in mean in the middle of the sequence.
d = 50
mu = 2
tau = 100
n = 200
set.seed(500)
y = rbind(matrix(rnorm(d*tau),tau), matrix(rnorm(d*(n-tau),mu/sqrt(d)), n-tau))
y.dist = dist(y)
E1 = mstree(y.dist)
# For illustration, we constructed the minimum spanning tree.
# You can use other ways to construct the graph.

r1 = gseg1(n,E1, statistics="all")  
# output results based on all four statistics
# the scan statistics can be found in r1$scanZ

r1_a = gseg1(n,E1, statistics="w")  
# output results based on the weighted edge-count statistic

r1_b = gseg1(n,E1, statistics=c("w","g"))  
# output results based on the weighted edge-count statistic 
# and generalized edge-count statistic

# image(as.matrix(y.dist))  
# run this if you would like to have some idea on the pairwise distance

## Sequence 2: change in mean away from the middle of the sequence.
d = 50
mu = 2
tau = 45
n = 200
set.seed(500)
y = rbind(matrix(rnorm(d*tau),tau), matrix(rnorm(d*(n-tau),mu/sqrt(d)), n-tau))
y.dist = dist(y)
E2 = mstree(y.dist)
r2 = gseg1(n,E2,statistic="all")
# image(as.matrix(y.dist))  


## Sequence 3: change in both mean and variance away from the middle of the sequence.
d = 50
mu = 2
sigma=0.7
tau = 145
n = 200
set.seed(500)
y = rbind(matrix(rnorm(d*tau),tau), matrix(rnorm(d*(n-tau),mu/sqrt(d),sigma), n-tau))
y.dist = dist(y)
E3 = mstree(y.dist)
r3=gseg1(n,E3,statistic="all")
# image(as.matrix(y.dist)) 


## Sequence 4: change in both mean and variance away from the middle of the sequence.
d = 50
mu = 2
sigma=1.2
tau = 50
n = 200
set.seed(500)
y = rbind(matrix(rnorm(d*tau),tau), matrix(rnorm(d*(n-tau),mu/sqrt(d),sigma), n-tau))
y.dist = dist(y)
E4 = mstree(y.dist)
r4=gseg1(n,E4,statistic="all")
# image(as.matrix(y.dist))  


## Sequence 5: change in both mean and variance happens on an interval.
d = 50
mu = 2
sigma=0.5
tau1 = 155
tau2 = 185
n = 200
set.seed(500)
y1 = matrix(rnorm(d*tau1),tau1)
y2 = matrix(rnorm(d*(tau2-tau1),mu/sqrt(d),sigma), tau2-tau1)
y3 = matrix(rnorm(d*(n-tau2)), n-tau2)
y = rbind(y1, y2, y3)
y.dist = dist(y)
E5 = mstree(y.dist)
r5=gseg2(n,E5,statistics="all")
# image(as.matrix(y.dist))  

}


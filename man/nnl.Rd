\name{nnl}
\alias{nnl}
\title{Construct the Nearest Neighbor Link (NNL)}
\description{This function provides a method to construct the NNL.}
\usage{
nnl(distance, K)
}
\arguments{
  \item{distance}{The distance matrix on the distinct values (a "number of unique observations" by "number of unique observations" matrix).}
  \item{K}{The value of k in "k-MST" or "k-NNL" to construct the similarity graph.}
}
\value{
  \item{E}{The edge matrix representing the similarity graph on the distinct values with the number of edges in the similarity graph being the number of rows and 2 columns. Each row records the subject indices of the two ends of an edge in the similarity graph.}
}
\seealso{
  \code{\link{gSeg}}, \code{\link{gseg1_discrete}}, \code{\link{gseg2_discrete}}
}
\examples{
n = 50
d = 10
dat = matrix(rnorm(d*n),n)
sam = sample(1:n, replace = TRUE)
dat = dat[sam,] 

# This data has repeated observations
dat_uni = unique(dat)
E = nnl(dist(dat_uni), 1)
}

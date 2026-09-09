\name{gseg_kmst}
\alias{gseg_kmst}
\title{Construct a k-MST Similarity Graph}
\usage{
gseg_kmst(x = NULL, dissimilarity = NULL, k = NULL)
}
\arguments{
\item{x}{Optional data matrix with observations in rows.}
\item{dissimilarity}{Optional square dissimilarity matrix or a \code{dist} object. Supply exactly one of \code{x} and \code{dissimilarity}.}
\item{k}{Number of complete MST layers. The default is \code{floor(sqrt(N))}; the largest allowed value is \code{floor(N/2)}.}
}
\value{A two-column edge matrix containing the union of the first \code{k} MST layers.}
\description{Constructs the k-MST used by default in the data-level gSeg interface.}
\examples{
set.seed(1)
x <- matrix(rnorm(60), nrow = 20)
E <- gseg_kmst(x = x)
attr(E, "k")
}

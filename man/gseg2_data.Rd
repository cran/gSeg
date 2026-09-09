\name{gseg2_data}
\alias{gseg2_data}
\title{Changed-Interval Graph Scan from Data or Dissimilarities}
\usage{
gseg2_data(x = NULL, dissimilarity = NULL, k = NULL,
  statistics = "m", ...)
}
\arguments{
\item{x}{Optional data matrix with observations in rows.}
\item{dissimilarity}{Optional square dissimilarity matrix or a \code{dist} object. Supply exactly one of \code{x} and \code{dissimilarity}.}
\item{k}{Number of complete MST layers. The default is \code{floor(sqrt(N))}; the largest allowed value is \code{floor(N/2)}.}
\item{statistics}{Statistics passed to \code{gseg2}. The default reports MET only; use \code{"g"} to request GET or \code{"all"} for all historical statistics.}
\item{...}{Further arguments passed to \code{gseg2}.}
}
\value{The result from \code{gseg2}, with the constructed graph and realized \code{k} added.}
\description{Constructs a k-MST and applies the changed-interval scan.}
\seealso{\code{\link{gseg2}}, \code{\link{gseg1_data}}, \code{\link{gseg_kmst}}}
\examples{
set.seed(1)
x <- matrix(rnorm(60), nrow = 20)
ans <- gseg2_data(x = x, pval.appr = FALSE)
names(ans$scanZ)
}

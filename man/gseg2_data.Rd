\name{gseg2_data}
\alias{gseg2_data}
\title{Data-Level Alias for gseg2}
\description{A compatibility alias for the data-level \code{\link{gseg2}}
interface. New code may call \code{gseg2} directly.}
\usage{
gseg2_data(x = NULL, dissimilarity = NULL, k = NULL,
  statistics = "m", ...)
}
\arguments{
\item{x}{A numeric matrix with observations in rows.}
\item{dissimilarity}{A square symmetric dissimilarity matrix or a \code{dist}
object. Supply exactly one of \code{x} and \code{dissimilarity}.}
\item{k}{Number of MST layers. The default is \code{floor(sqrt(N))}.}
\item{statistics}{Scan statistic or statistics; see \code{\link{gseg2}}.}
\item{...}{Additional arguments passed to \code{gseg2}.}
}
\value{The result from \code{gseg2}.}
\seealso{\code{\link{gseg2}}, \code{\link{gseg_kmst}}}
\examples{
set.seed(1)
x <- matrix(rnorm(60), nrow = 20)
ans <- gseg2_data(x, pval.appr = FALSE)
}

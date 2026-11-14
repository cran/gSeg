\name{gseg1_data}
\alias{gseg1_data}
\title{Data-Level Alias for gseg1}
\description{A compatibility alias for the data-level \code{\link{gseg1}}
interface. New code may call \code{gseg1} directly.}
\usage{
gseg1_data(x = NULL, dissimilarity = NULL, k = NULL,
  statistics = "m", ...)
}
\arguments{
\item{x}{A numeric matrix with observations in rows.}
\item{dissimilarity}{A square symmetric dissimilarity matrix or a \code{dist}
object. Supply exactly one of \code{x} and \code{dissimilarity}.}
\item{k}{Number of MST layers. The default is \code{floor(sqrt(N))}.}
\item{statistics}{Scan statistic or statistics; see \code{\link{gseg1}}.}
\item{...}{Additional arguments passed to \code{gseg1}.}
}
\value{The result from \code{gseg1}.}
\seealso{\code{\link{gseg1}}, \code{\link{gseg_kmst}}}
\examples{
set.seed(1)
x <- matrix(rnorm(60), nrow = 20)
ans <- gseg1_data(x, pval.appr = FALSE)
}

\name{gSeg}
\alias{gSeg}
\title{Graph-Based Change-Point Detection}
\description{
The package detects distributional changes in sequences of multivariate or
object-valued observations. The main functions accept either a numeric data
matrix or a pairwise dissimilarity matrix. They construct a k-MST with
\code{k=floor(sqrt(N))} by default and report the max-type edge-count test
(MET). Users may instead provide their own similarity graph or explicitly
request GET, OET, WET, or all four statistics.

Use \code{\link{gseg1}} when the sequence may contain one change-point and
\code{\link{gseg2}} when it may contain one changed interval. For repeated
observations, use \code{\link{gseg1_discrete}} or
\code{\link{gseg2_discrete}}.
}
\author{
Hao Chen, Nancy R. Zhang, Lynna Chu, and Hoseung Song

Maintainer: Hao Chen (hxchen@ucdavis.edu)
}
\references{
Chen, Hao, and Nancy Zhang. (2015). Graph-based change-point detection.
The Annals of Statistics, 43(1), 139--176.

Chu, Lynna, and Hao Chen. (2019). Asymptotic distribution-free change-point
detection for modern data. The Annals of Statistics, 47(1), 382--414.

Song, Hoseung, and Hao Chen. (2022). Asymptotic distribution-free
change-point detection for data with repeated observations. Biometrika,
109(3), 783--798.
}
\seealso{
\code{\link{gseg1}}, \code{\link{gseg2}},
\code{\link{gseg1_graph}}, \code{\link{gseg2_graph}},
\code{\link{gseg1_discrete}}, \code{\link{gseg2_discrete}},
\code{\link{gseg_kmst}}, \code{\link{nnl}}
}
\examples{
set.seed(500)
N <- 40
d <- 5
x <- rbind(matrix(rnorm(d * 20), 20),
           matrix(rnorm(d * 20, 1 / sqrt(d)), 20))

# The main call accepts the observations directly. It constructs the default
# k-MST and reports MET.
result <- gseg1(x, pval.appr = FALSE)
names(result$scanZ)
result$k

# Pairwise dissimilarities can be supplied when observations are not vectors.
result.d <- gseg1(dissimilarity = as.matrix(dist(x)), pval.appr = FALSE)

# GET is available by explicit request.
result.get <- gseg1(x, statistics = "g", pval.appr = FALSE)

# A user-supplied graph can be analyzed explicitly.
result.graph <- gseg1_graph(N, result$graph, pval.appr = FALSE)
}

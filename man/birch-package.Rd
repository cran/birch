\name{birch-package}
\alias{birch-package}
\docType{package}
\title{
  Working with very large data sets using BIRCH.
}
\description{
  The functions in this package are designed for working with very large
  data sets, but pre-processing them with an algorithm called BIRCH
  (Balanced Iterative Reducing and Clustering using Hierarchies), which
  represents a data set by compact subclusters, each with summary
  statistics (called clustering features). Then, instead of using the
  full data set, these subclusters can be used.

  Commonly, the most advantage with this approach can be obtained with
  combinatorial problems, since (a) the whole data set need not be held
  in memory, and (b) the dimension of the solution space is vastly
  reduced, making it easier to maximise/minimise the objective function.

  A complete explanation of this package is given in \cite{Harrington
  and Salibian-Barrera (2008)}, and discussion of the underlying
  algorithm can be found in \cite{Harrington
  and Salibian-Barrera (2007)}.

  Documentation for developers can be found in the doc subdirectory of
  this package.
}
\details{
\tabular{ll}{
Package: \tab birch\cr
Type: \tab Package\cr
Version: \tab 1.1-2\cr
Date: \tab 2008-02-22\cr
License: \tab GPL\cr
}
The main functions are \code{\link{birch}}, \code{\link{covMcd.birch}},
\code{\link{lts.birch}}, \code{\link{rlga.birch}}, and
\code{\link{kmeans.birch}}. Various generic methods are present,
including those for \code{print}, \code{summary} and \code{plot}.
}
\author{
Justin Harrington \email{harringt@stat.ubc.ca} and Matias Salibian-Barrera \email{matias@stat.ubc.ca}
}
\references{
    Harrington, J and Salibian-Barrera, M (2007) \dQuote{Finding
      Approximate Solutions to Combinatorial Problems with Very Large
      Datasets using BIRCH}, submitted to \emph{Statistical Algorithms
      and Software, 2nd Special Issue} Computational Statistics and Data
    Analysis. A draft can be found at \url{http://www.stat.ubc.ca/~harringt/birch/birch.pdf}.
    
    Harrington, J and Salibian-Barrera, M (2008) \dQuote{birch: Working
    with very large data sets}, submitted to \emph{Journal of
    Statistical Software}. A draft can be found at \url{http://www.stat.ubc.ca/~harringt/birch/birch-jss.pdf}.
}
\keyword{ package }
\seealso{
  \code{\link{birch}}, \code{\link{covMcd.birch}},
  \code{\link{lts.birch}}, \code{\link{rlga.birch}}

}

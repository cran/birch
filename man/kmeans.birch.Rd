\name{kmeans.birch}
\alias{kmeans.birch}
\title{ K-means with BIRCH}
\description{
  Perform the k-means clustering algorithm using a birch object.
}
\usage{
kmeans.birch(birchObject, centers, nstart = 1)
}
\arguments{
  \item{birchObject}{ An object created by the function \code{\link{birch}}.} 
  \item{centers}{ Either the number of clusters or a set of initial (distinct)
    cluster centres.  If a number, a random set of (distinct)
    subclusters in \sQuote{birchObject} is chosen as the initial centres.}
  \item{nstart}{ If 'centers' is a number, how many random sets should be
    chosen?}
}

\details{
  The birch object given by \sQuote{birchObject} is clustered by the k-means
  method, adjusted for dealing with birch objects. The aim is to
  partition the subclusters into k groups such that the sum of squares
  of all the points in each subcluster to the assigned cluster centers
  is minimized.
  
  The result should be approximately similar to that found by performing
  k-means on the original data set. However, this approximation depends
  on the \dQuote{coarseness} of the underlying tree, and the size of the
  combinatorial problem. These aspect as discussed in detail in the
  references.
}
\note{
  In order for this algorithm to produce meaningful results, the number
  of subclusters in the birch object should number in the hundreds, and even better, thousands.
  }

\value{
  Returns a list with components:
  \item{RSS }{The total residual sum-of-squares of the clustering.}
  \item{clust }{A list containing a vector of which subclusters make up the
    clustering (sub) and a vector with the underlying observations that make up the
    clusters (obs)}
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
\author{
  Justin Harrington \email{harringt@stat.ubc.ca} and Matias Salibian-Barrera \email{matias@stat.ubc.ca}
}
\seealso{ \code{\link{birch}}, \code{\link{dist.birch}}}
\examples{
## Load a demo birch Object
data(birchObj)

## Perform k-means, specifying the number of clusters
kOut <- kmeans.birch(birchObj, 2, nstart=10)

## Perform k-means, specifying the initial cluster centers
## See dist.clust for one method of initial cluster centers
kOut <- kmeans.birch(birchObj, matrix(c(0,10), ncol=5, nrow=2))

## To plot using the birch object
plot(birchObj, col=kOut$clust$sub)

## To plot using the underlying data (if available)
\dontrun{plot(x, col=kOut$clust$obs)}
}
\keyword{multivariate}
\keyword{cluster}


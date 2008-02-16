\name{generic methods}
\alias{print.birch}
\alias{summary.birch}
\alias{[.birch}
\alias{dim.birch}
\alias{length.birch}
\alias{cbind.birch}
\alias{plot.birch}
\alias{pairs.birch}
\alias{points.birch}
\alias{lines.birch}
\title{Generic methods for birch objects}
\description{
  A description of the generic methods that are compatible with birch objects.
}

\usage{
\method{print}{birch}(x, \dots)
\method{summary}{birch}(object, \dots)
\method{[}{birch}(x, i, j)
\method{dim}{birch}(x)
\method{length}{birch}(x)
\method{cbind}{birch}(scalar, x)
\method{plot}{birch}(x, centers=FALSE, \dots)
\method{pairs}{birch}(x, centers=FALSE, \dots)
\method{points}{birch}(x, \dots)
\method{lines}{birch}(x, col, \dots)
}

\arguments{
  \item{x, object}{A birch object}
  \item{i, j}{Indices specifying which subclusters and dimensions to extract.}
 \item{scalar}{A scalar to add as a column to the birch object.}
\item{centers}{A Boolean. If TRUE, then just plots the centers of the
subclusters. See the details for more information.}
\item{col}{Color argument for lines.}
\item{\dots}{Pass additional information to the respective methods.}
}

\value{
The \code{dim} method returns a vector of length three
  containing the total number of observations in the tree, the number of
  columns of the data set, and the number of subclusters. The
\code{length} method just returns the number of subclusters. Note that
both of these commands operate on the submitted object, and do not check
if the tree has been updated (e.g. by adding data with the
\code{birch.addToTree} command).

The \code{summary} method returns the mean and covariance of the whole
birch object - equivalently, the whole underlying data set.

The \sQuote{[i,j]} method selects the i-th subcluster and the j-th
column of the birch object, and has functionality similar to that of the
usual indexing of matrices. Similarly, the \code{cbind} method
effectively inserts a column containing a scalar (recycled to the
appropriate length) in front of the data set. 

  The \code{print}, \code{plot}, \code{pairs}, \code{points} and
  \code{lines} methods return  nothing.
}

\details{
 Most of these methods do not require further explanation, as their
 behaviour is fairly standard. The plot and pairs methods, does something
 a little different. 

While they accept all the usual arguments that plot has, instead of
producing plots of the values, which is not possible with a
birch object, instead these methods produce ellipses for each
subclusters, based on its covariance structure, and making the
assumption that the underlying is multivariate normal. The ellipses are
then formed (using the ellipse package) at the 95\% confidence level.

However, if there are a lot of subclusters, then this plot gets
\dQuote{messy}, and so it is possible to plot just the centers of each
subcluster by setting the \code{centers} argument to true. The
\code{points} and \code{lines} methods add the centers and ellipses
respectively to the existing plot.
}
\note{
Care should be taken when using the indexing, as
\preformatted{
birch(z, 1)[,1:2]
}
will not produce the same result as
\preformatted{
birch(z[,1:2], 1)
}
Similarly, there are no guarantees for \code{cbind}.
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
    Statistical
    Software}. A draft can be found at \url{http://www.stat.ubc.ca/~harringt/birch/birch-jss.pdf}. 
}
\seealso{\code{\link{birch}}}
\examples{
## Load demo birch object
data(birchObj)
dim(birchObj)
newObj <- cbind(1, birchObj)
dim(newObj)

dim(birchObj[-1, -2])
}

\keyword{multivariate}
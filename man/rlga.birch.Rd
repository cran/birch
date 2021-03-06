\name{rlga.birch}
\alias{rlga.birch}
\alias{lga.birch}
\title{ Finding the Robust/Non-Robust LGA solution using BIRCH }
\description{
  Performs Linear Grouping Analysis (LGA) and Robust Linear Grouping
  Analysis (RLGA) using a BIRCH object.
}
\usage{
lga.birch(birchObject, k, nsamp=100)
rlga.birch(birchObject, k, alpha=0.5, nsamp=100)
}
\arguments{
  \item{birchObject}{an object created by the function \code{\link{birch}}.}
  \item{k}{the number of clusters}
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the orthogonal residuals are minimized, i.e., alpha*n observations
    are used when calculating orthogonal residuals in
    the hyperplane calculations. Allowed values are between
    0.5 and 1 and the default is 0.5.}
  \item{nsamp}{number of subsets used for initial estimates}
}
\details{
  Robust Linear Grouping (\cite{Garcia-Escudero et al, 2010}) is the
  robust implementation of LGA (\cite{Van Aelst et al, 2006}), and is
  concerned with clustering around hyperplanes. The non-birch versions can
  be found in the package \code{lga}, which is also available on CRAN.
  
  This algorithm is the equivalent design for use with a BIRCH object. For
  further details, please see \cite{Harrington and Salibian-Barrera (2010)}.
}
\value{
  Returns a list containing:
  \item{clust}{A list containing a vector of which subclusters make up the
    clustering (sub) and a vector with the underlying observations that make up the
    clusters (obs). For the robust algorithm, a value of zero indicates it does not belong to the
    best h-subset.}
  \item{ROSS}{the residual sum of squares of orthogonal distances to the
    fitted hyperplanes based on the best data set.}
}
\note{
  In order for this algorithm to produce meaningful results, the number
  of subclusters in the birch object should be in the hundreds, and even better, thousands.
}
\references{ 
  Garcia-Escudero, L.A. and Gordaliza, A. and San Martin, R. and Van
  Aelst, S. and Zamar, R. (2007) \dQuote{Robust Linear Clustering},
  \emph{Unpublished Manuscript}.
  
  Harrington, J. and Salibian-Barrera, M. (2010), \dQuote{Finding
    Approximate Solutions to Combinatorial Problems with Very Large
    Datasets using BIRCH}, \emph{Computational Statistics and Data
    Analysis} 54, 655-667 . 
  
  Harrington, J. and Salibian-Barrera, M. (2008), \dQuote{birch: Working
    with very large data sets}, working paper. 
  
  Van Aelst, S. and Wang, X. and Zamar, R.H. and Zhu, R. (2006)
  \dQuote{Linear grouping using orthogonal regression},
  \emph{Computational Statistics \& Data Analysis} \bold{50}, 1287--1312.
}
\author{
Justin Harrington \email{harringt@stat.ubc.ca} and Matias Salibian-Barrera \email{matias@stat.ubc.ca}
}
\seealso{\code{\link{birch}}, and the non-birch algorithms
  \code{\link[lga:lga]{lga}} and \code{\link[lga:rlga]{rlga}}} 
\examples{
library(MASS) ## for mvrnorm
library(birch)

## Create new data set (that is more applicable to RLGA and LGA
set.seed(1234) 
x <- mvrnorm(1e4, mu=rep(0,2), Sigma=diag(c(0.25,1),2))
x <- rbind(x, mvrnorm(1e4, mu=rep(10,2),
                      Sigma=diag(c(5,0.5),2)))

## Create birch object, and save it
birchObj <- birch(x, 0.5)
length(birchObj)

library(birch)
rlgaOut <- rlga.birch(birchObj, k=2, 0.5)
plot(birchObj, col=rlgaOut$clust$sub+1)

lgaOut <- lga.birch(birchObj, k=2)
plot(birchObj, col=lgaOut$clust$sub)
}
\keyword{robust}
\keyword{multivariate}
\keyword{cluster}

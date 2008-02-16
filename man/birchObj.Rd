\name{birchObj}
\alias{birchObj}
\docType{data}
\title{A birch object for demonstration purposes}

\description{A birch object created from a simulated data for the
  purposes of demonstration}

\usage{data(birchObj)}

\format{
  A birch object with 318 subclusters. The object contains a list of
  \describe{
    \item{\code{sumXi}}{a numeric matrix containing the column sum of
      the observations in each subcluster.}
    \item{\code{sumXisq}}{a numeric array containing the sum-of-squares
      of the observations in each subcluster.}
    \item{\code{N}}{a numeric vector containing the number of
      observations in each subcluster.}
    \item{\code{members}}{a list, each element of which contains a
      vector of the observations numbers present in that subcluster.}
  }
}

\details{
  This birch object was created from a simulated data set of two
  multivariate normals. Its primary use is for demonstrating functions
  in the help files. The sequence of commands for creating the object
  was as follows:
  \preformatted{
    library(MASS) ## for mvrnorm
    library(birch)

    ## Create data set
    set.seed(1234) 
    x <- mvrnorm(1.5e5, mu=rep(0,5), Sigma=diag(1,5))
    x <- rbind(x, mvrnorm(1e5, mu=rep(10,5), Sigma=diag(0.1,5)+0.9))

    ## Create birch object, and save it
    birchObj <- birch(x, 5)
    save(birchObj, file="birchObj.rda")
  }
}

\references{
  Harrington, J and Salibian-Barrera, M (2007) \dQuote{Finding
    Approximate Solutions to Combinatorial Problems with Very Large
    Datasets using BIRCH}, submitted to \emph{Statistical Algorithms
    and Software, 2nd Special Issue} Computational Statistics and Data
  Analysis. \url{http://www.stat.ubc.ca/~harringt/birch.pdf}.
  
  Harrington, J and Salibian-Barrera, M (2008) \dQuote{birch: Working
    with very large data sets}, submitted to \emph{Journal of
    Statistical Software}. \url{http://www.stat.ubc.ca/~harringt/birch-jss.pdf}.
  }
\examples{
data(birchObj)
}
\keyword{datasets}
\keyword{multivariate}
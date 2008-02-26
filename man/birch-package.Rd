\name{birch-package}
\alias{birch-package}
\docType{package}
\title{
  Working with very large data sets using BIRCH.
}
\description{
  
  The functions in this package are designed for working with very large
  data sets by pre-processing the data set with an algorithm called BIRCH
  (Balanced Iterative Reducing and Clustering using Hierarchies), which
  transforms the data set into compact, locally similar subclusters,
  each with summary statistics attached (called clustering
  features). Then, instead of using the full data set, these summary
  statistics can be used.

  This approach is most advangeous in two situations: when the data
  cannot be loaded into memory due to its size; and/or when some form of
  combinatorial optimization is required and the size of the solution
  space makes finding global maximums/minimums difficult.
  
  A complete explanation of this package is given in \cite{Harrington
  and Salibian-Barrera (2008)}, and discussion of the underlying
  algorithms can be found in \cite{Harrington
  and Salibian-Barrera (2007)}.

  Documentation for developers can be found in the doc directory of
  this package. Also, the source code contains dOxygen tags for further
  information.
  
}
\details{
  The main function is
  \tabular{ll}{
    \code{\link{birch}} \tab takes a data set (an R object, text file, etc),
    and creates a birch object\cr
  }
  Various generic methods are present, including
  \tabular{l}{
    \code{print} \cr
    \code{summary}\cr
    \code{plot}\cr
  }
  Finally, some combinatorial-style problems have been
  implemented. These include:
  \tabular{ll}{
  \code{\link{covMcd.birch}} \tab Minimum Covariance Determinant
  (robust estimator for location and dispersion)\cr
  \code{\link{lts.birch}} \tab Least Trimmed Squares (robust regression estimator) \cr
  \code{\link{rlga.birch}} \tab Robust Linear Grouping Analysis (robust
  clustering about hyperplanes)\cr
  \code{\link{kmeans.birch}} \tab k-means\cr
  }
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
  \code{\link{lts.birch}}, \code{\link{rlga.birch}},
  \code{\link{plot.birch}}.
}


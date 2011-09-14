\name{birch}
\alias{birch}
\alias{birch.addToTree}
\alias{birch.getTree}
\alias{birch.killTree}
\title{Create a birch object }
\description{
  This function creates a birch object using the algorithm BIRCH. 
}

\usage{
birch(x, radius, compact=radius, keeptree=FALSE, \dots)
birch.addToTree(x, birchObject, updateDIM = TRUE, ...)
birch.getTree(birchObject)
birch.killTree(birchObject)
}
\arguments{
  \item{x}{a numeric matrix of at least two columns, a file name or a
    connection that is compatible with \code{read.table}.}
  \item{radius}{the closeness criterion}
  \item{compact}{the compactness criterion}
  \item{keeptree}{A Boolean, whether to keep the CF tree in memory.}
  \item{\dots}{Arguments to be passed to \code{read.table} for loading a
    file or connection}
  \item{updateDIM}{Update the dimension of the object? Defaults to TRUE
  (which is desirable!).}
  \item{birchObject}{The output from birch.}
}

\value{
  For \code{birch} and \code{birch.getTree}, a list for each subcluster containing
  \item{N}{the number of observations in the subcluster}
  \item{sumXi}{the linear sum of the observations in the subcluster}
  \item{sumXisq}{the sum of squares of the observations in the subcluster}
  \item{members}{a list, each element containing a observation index
    (membership vector) for the subcluster}
  Both of \code{birch.addToTree} and \code{birch.killTree} return NULL.
}

\details{
  This function creates a CF-Tree not unlike that described in \cite{Zhang
    et al. (1997)}, and used in \cite{Harrington and Salibian-Barrera
    (2010)}. A complete explanation of this package is given in
  \cite{Harrington and Salibian-Barrera (2008)}.

  A full tree structure is used, as is the splitting of nodes
  (as described in the original article). However, the \sQuote{Merging
    Refinement} as described on page 149 is not currently implemented. The
  automatic rebuilding based on page size is not implemented.

  The argument \code{keeptree} allows for the tree to be kept in memory
  after the initial processing has been completed. This allows for
  additional information to be added at a latter stage with
  \code{birch.addToTree}, without needing to process the whole data set
  a second time. However, it should be noted that the CF data (the
  summary statistics of each subcluster, as used by the subsequent
  algorithms) is not returned , and so this command should be followed
  by \code{birch.getTree} command in order to ensure that the correct
  information is used. In other words,
  \preformatted{
    ## Create the tree
    myobject <- birch(x, 0.1, keeptree=TRUE)
    
    ## myobject has no information - let's get some
    myobject <- birch.getTree(myobject)

    ## add some data
    birch.addToTree(y, mybirch)

    ## myobject is now out of date! Until
    mybirch <- birch.getData(mybirch)
  }

  A birch object without summary information is referred to as a
  \dQuote{birch skeleton}. Most algorithms will check for the presence
  of summary data, and request the information if required, but it is
  better for the user to request the data directly.
  
  The birch object produced by the algorithm has a number of compatible
  generic methods, as well as clustering and robustness
  algorithms. These are given in the \sQuote{See Also} section.

  The selection of the parameters for forming the birch object is left
  up to the user, though some guidance is given in \cite{Harrington and
  Salibian-Barrera (2010)}. One consideration is the purpose - for
  example, if simple summary statistics are required (means, covariances
  etc), then a large compactness and radius can be selected, as
  granularity is not required. In the case of clustering and robust
  methods when a refinement step is being done after the birch algorithm
  is completed, then once again the radius and compactness can be
  larger. If there is no refinement step, however, then the selection of
  these parameters is much more important.
  
}
\note{
  
  In this implementation, a limit of matrices with 30 columns has been
  set. It is the intention that this limit will be eventually removed in
  future versions of the package. If the user desires, this limit can be
  increased in the file \sQuote{ll.h}. Other than the matrix column
  limit, there is no practical limit (excluding those imposed by the
  operating system and R) for how many observations can be placed in a
  birch object. In particular, since we retain the observation numbers
  belonging to each subcluster, the maximum length of a vector
  containing integers will determine how large the birch tree can get.

  Some information for developers is provided in the doc directory of
  this package.
}
  
\author{
  Justin Harrington \email{harringt@stat.ubc.ca} and Matias Salibian-Barrera \email{matias@stat.ubc.ca}
}

\references{
  Harrington, J. and Salibian-Barrera, M. (2010), \dQuote{Finding
    Approximate Solutions to Combinatorial Problems with Very Large
    Datasets using BIRCH}, \emph{Computational Statistics and Data
  Analysis} 54, 655-667.
  
  Harrington, J. and Salibian-Barrera, M. (2008), \dQuote{birch: Working
    with very large data sets}, working paper.
  
  Zhang, T. and Ramakrishnan, R. and Livny, M. (1997) \dQuote{BIRCH: A New
    Data Clustering Algorithm and Its Applications}, \emph{Data Mining and
    Knowledge Discovery} \bold{1}, 141--182.
  
}
\seealso{\code{\link{lts.birch}},  \code{\link{rlga.birch}},
  \code{\link{covMcd.birch}}, \code{\link{kmeans.birch}},
  \code{\link{lga.birch}}, and \code{\link{plot.birch}} which also
  contains the other generic methods}
\examples{
## an example
## Create a data set
library(MASS)
set.seed(1234) 
x <- mvrnorm(1e5, mu=rep(0,5), Sigma=diag(1,5))
x <- rbind(x, mvrnorm(1e5, mu=rep(10,5), Sigma=diag(0.1,5)+0.9))

## Create birch object
birchObj <- birch(x, 5, keeptree = TRUE)

## To load directly from a file or connection
\dontrun{birchObj <- birch("myfile.csv", 1, sep=",", header=TRUE)}
\dontrun{birchObj <- birch("http://www.dot.com/myfile.csv", 1, sep=",",
header=TRUE)}

## Leaving a tree in memory
birchObj <- birch(x, 5, keeptree=TRUE)
birch.addToTree(x, birchObj)
birchObj <- birch.getTree(birchObj)
## And don't forget to killTree...
#birch.killTree(birchObj)
}

\keyword{multivariate}

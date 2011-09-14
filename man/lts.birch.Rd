\name{lts.birch}
\alias{lts.birch}
\alias{ltsBirch.refinement}
\title{ Finding the Least Trimmed Squares (LTS) regression estimate using BIRCH}
\description{
  This algorithm searches for the Least Trimmed Squares (LTS) solution using a BIRCH object.
}
\usage{
lts.birch(birchObject, alpha=0.5, intercept = FALSE, nsamp=100)
ltsBirch.refinement(ltsOut, x, y, alpha=0.5, intercept = FALSE)
}
\arguments{
  \item{birchObject}{an object created by the function
    \code{birch}. See details for information on how to specify the
    exploratory and response.} 
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the trimmed residuals are minimized, i.e., alpha*n observations
    are used when computing the trimmed residual sum of squares . Allowed values are between
    0.5 and 1 and the default is 0.5.} 
  \item{intercept}{a Boolean - is there a intercept? }
  \item{nsamp}{number of subsets used for initial estimates}
  \item{ltsOut}{the output from \code{lts.birch}.}
  \item{x, y}{a data set of explanatory and response variables on which
  to perform a set of concentration steps.}
}
\details{
  The algorithm is very similar to the \code{ltsRef} function from the robustbase
  package from \cite{Rousseeuw and Van Driessen (2006)}, except it uses
  a BIRCH object instead. A complete description is given in \cite{Harrington and
  Salibian-Barrera (2010)} and \cite{Harrington and Salibian-Barrera (2008)}

  The algorithm assumes that the last column of the birch object
  contains the response variable, and that all the other columns are
  explanatories. While it is possible to select columns using the usual \preformatted{[,j]},
  it is recommended that the birch object be rebuilt from the underlying
  data set with just the explanatories and response variables selected.
  
  If an intercept is required in the model, either the \code{intercept}
  argument can be set to true, or a column of \sQuote{ones}
  should be column-wise appended to the data (prior to building the
  birch object).

  A summary method is available for the output of this command.
}
\note{
  In order for this algorithm to produce meaningful results, the number
  of subclusters in the birch  object should number in the hundreds, and even better, thousands.
}
\value{
  Returns a list containing:
  \item{best}{A list containing a vector of which subclusters make up the
    clustering (sub) and a vector with the underlying observations that make up the
    clusters (obs).}
  \item{raw.coefficients}{the fitted LTS regression line.}
  \item{Resids}{A list containing the sum of squared residuals for the
    best subset, as well as the sum of squared residuals for the whole
    data set (based on the LTS regression equation).}
}
\references{
  Harrington, J. and Salibian-Barrera, M. (2010), \dQuote{Finding
    Approximate Solutions to Combinatorial Problems with Very Large
    Datasets using BIRCH}, \emph{Computational Statistics and Data
  Analysis} 54, 655-667.
  
  Harrington, J. and Salibian-Barrera, M. (2008), \dQuote{birch: Working
    with very large data sets}, working paper. 
  
  Rousseeuw, P.J. and Van Driessen, K. (2006) \dQuote{Computing LTS
    Regression for Large Data Sets}, \emph{Data Mining and Knowledge
    Discovery} \bold{12}, 29--45.
}
\author{
  Justin Harrington \email{harringt@stat.ubc.ca} and Matias Salibian-Barrera \email{matias@stat.ubc.ca}
}
\seealso{ \code{\link{birch}}, and the original algorithm  \code{\link[robustbase:ltsReg]{ltsReg}}}  
\examples{
data(birchObj)
ltsOut <- lts.birch(birchObj, 0.5)
ltsOut2 <- lts.birch(birchObj, 0.5, intercept=TRUE)
summary(ltsOut2)

## If the original data set was available
\dontrun{refOut <- ltsBirch.refinement(ltsOut2, x, y, 0.5, intercept=TRUE)}
}
\keyword{robust}
\keyword{multivariate}

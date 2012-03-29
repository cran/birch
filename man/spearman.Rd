\name{spearman}
\alias{spearman}
\title{Estimating Spearman's correlation matrix using BIRCH}
\description{
This function builds an estimate of the correlation matrix from a birch 
object or from a dataset. In the last case, it uses birch algorithm to pre-process
the data. 

The Spearman's rho is used as the correlation measure.
}
\usage{
spearman(x, radius, compact = radius, columns = NULL,  estimate , ... )
}
\arguments{
  \item{x}{a birch object, a numeric matrix of at least two columns, a file name or a
    connection that is compatible with \code{read.table}.}
  \item{radius}{the closeness criterion (not used if x is a birch object)}
  \item{compact}{the compactness criterion (not used if x is a birch object)}
  \item{columns}{the choice of variables to use if x is a file name or a connection (default is all variables)}
  \item{estimate}{the choice of corrections to use for the estimation ("Student" and "Woodbury")}
  \item{\dots}{arguments to be passed to \code{read.table} for loading a
    file or connection}
 
}
\details{
  There are two versions of the function. When x is a matrix, a file name or a connection, the 
  first version is used. When x is a birch object, the second version is used.
  
  First version:
  For each pair of variables, one birch object is created using the radius and compactness 
  criteria. The estimations are based on the summary statistics resulting of the usage 
  of birch algorithm that regroups data in compact clusters. 
  The advantage of this approach is that it does not require the full data set to be held 
  in memory. Further details can be found in
  \cite{Charest (2012)}.
  
  The result of the algorithm is a correlation matrix with estimated correlations for each 
  pair of variables.
  
  Second version:
  For each pair of variables, the correlation is estimated based on the summary statistics 
  corresponding to them. 
  
  Estimates:
  \describe{
    \item{"Student", "student", '1'}{Correction using the product-moment correlation with mid-ranks (Kendall (1945), Student (1921))}
    \item{"Woodbury", "woodbury", '2'}{Basic formula without correction (Woodbury, 1945)}
  }	   
  

}
\value{
  Returns a list with components:
    \item{nClus}{matrix showing the number of clusters used for the estimation of each coefficient}
    \item{estimate[,,i]}{one correlation matrix estimate for each selected correction}
}

  
\references{
  Charest (2012)
  
 
}

\author{
Lysiane Charest \email{lysiane.charest@hec.ca} 
}
\seealso{\code{\link{birch}}, \code{\link{kendall}}}

\examples{
## Create a data set
library(MASS)
set.seed(1234) 
x <- mvrnorm(1e5, mu=rep(0,5), Sigma=diag(1,5))
x <- rbind(x, mvrnorm(1e5, mu=rep(10,5), Sigma=diag(0.1,5)+0.9))

##Create csv for the example
write.csv(x,"csvExample.csv", row.names=FALSE)

## Create a birch object
birchObj <- birch(x, radius = 5, keeptree = TRUE)
birchObj <- birch.getTree(birchObj)

##With x: birch object. 
##In this example, altough the birch object was created with V1-V5, 
##we only want the estimates for V2, V4 and V5
fromBirchObj <- spearman(birchObj[,c(2,4,5)], estimate = c("student", "woodbury"))


##With x: dataset. 
##In this example, we requiere estimates for every variables
fromDataset <- spearman(x, radius = 5, estimate = c("student", "woodbury"))

##With x: file name or connection
##In this example, we requiere estimates for variables V2, V4 and V5 only
fromFile <- spearman("csvExample.csv", 1, sep=",", header=TRUE, radius = 5, columns = c(2,4,5), 
estimate = c("student", "woodbury"))
\dontrun{fromConnection <- spearman("http://www.dot.com/myfile.csv", 1, sep=",",
header=TRUE, radius = 5, columns = c(2,4,5), estimate = c("student", "woodbury"))}

}



\keyword{nonparametric}
\keyword{robust}

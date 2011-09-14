
 ## Imprimer info objet birch
print.birch <- function(x, ...){
  cat("Birch Tree\n")
  cat("Built with radius =", attr(x, "radius"), "and compact =",  attr(x, "compact"), "\n")
  cat("Number of leaves:", length(x), "\n")
  cat("Number of underlying observations:", dim(x)[1], "\n")
  if (is.null(x$sumXi)) cat("NB: x is a birch 'skeleton'. No data currently loaded into the object.\n See ?birch for more details.\n")
}

summary.birch <- function(object, ...){
  if (is.null(object$sumXi))
    object <- birch.getTree(object)

  N <- dim(object)[1]
  means <- colSums(object$sumXi)/N
  names(means) <- attr(object, "xcolnames")
  covs <- (rowSums(object$sumXisq, dim=2) - N*means %o% means)/(N-1)
  dimnames(covs) <- list(attr(object, "xcolnames"), attr(object, "xcolnames"))
  outputs <- list(means=means, covs=covs)
  class(outputs) <- "summary.birch"
  return(outputs)
}

print.summary.birch <- function(x, ...){
  cat("Summary data for underlying data\n")
  cat("Column Means:\n")
  print(x$means, ...)
  cat("\nCovariance Matrix:\n")
  print(x$covs, ...)
}


#"<-.birch" <- function() {
#    return(value[])
#}

"[.birch" <- function(x, i= seq(attr(x, "xdim")[3]), j= seq(attr(x, "xdim")[2]) ){
    #print("In crochet")
    iFull= seq(attr(x, "xdim")[3])
    jFull= seq(attr(x, "xdim")[2])
    i= iFull[i]
    j=jFull[j]

  if (is.null(x$sumXi))
    y <- birch.getTree(x)
  else
    y <- x

    if (is.null(attr(x,"internal"))) #No tree to copy. We only update info.
    {
        y$sumXi <- y$sumXi[i, , drop=FALSE]
        y$sumXisq <- y$sumXisq[ , , i, drop=FALSE]
        y$N <- y$N[i]
        attr(y, "xdim")[c(1,3)] <- c(sum(y$N), length(y$N))

        y$sumXi <- y$sumXi[, j, drop=FALSE]
        y$sumXisq <- y$sumXisq[j, j, , drop=FALSE]
        attr(y, "xcolnames") <- attr(x, "xcolnames")[j]
        attr(y, "xdim")[2] <- ncol(y$sumXi)
        attr(y, "radius") <- attr(y, "compact") <- NA
        #print("No tree behind the birch object. Be careful.")
        return(y)
    }
    else
    {
        lea <-  as.vector(i)
        var <- as.vector(j)
        oldD <-attr(x, "xdim")[2]

        attr(y, "internal") <- .Call("LL_assign", attr(x,"internal"),  var, lea, oldD)

        if (!missing(j)) {
           attr(y, "xcolnames") <- attr(x, "xcolnames")[j]
           attr(y, "radius") <- attr(y, "compact") <- NA
        }
        return(birch.getTree(y))
    }
}




dim.birch <- function(x)
  return(attr(x, "xdim"))

length.birch <- function(x)
  return(attr(x, "xdim")[3])

 ## Ajout colonne en avant des autres
cbind.birch <- function(scalar, x){
  if (is.null(x$sumXi))
    y <- birch.getTree(x)
  else
    y <- x
  d <- dim(x)[2] + 1
  #print(d)

  y$sumXi <- cbind(scalar*x$N, x$sumXi)
  y$sumXisq <- array(0, dim=c(d, d, length(y)))
  y$sumXisq[1, 2:d,] <- y$sumXisq[2:d, 1, ] <- t(scalar*x$sumXi)
  y$sumXisq[1,1,] <- scalar^2*x$N
  y$sumXisq[2:d, 2:d, ] <- x$sumXisq
   ## New nb of columns
  attr(y, "xdim")[2] <- d
  attr(y, "radius") <- attr(y, "compact") <- NA
   ## Add name of new column
  if (!is.null(attr(y, "xcolnames")))
    attr(y, "xcolnames") <- c("intercept", attr(y, "xcolnames"))
  attr(y,"internal") <- NULL
  return(y)
}

plot.birch <- function(x, centers=FALSE, ...){
  if (is.null(x$sumXi))
    x <- birch.getTree(x)

  if (dim(x)[2] > 2)
    x <- x[,1:2]
  if (centers)
    plot.birch.centers(x, ...)
  else
    plot.birch.ellipse(x, ...)
}

plot.birch.ellipse <- function(x, col, ...){
  d <- dim(x)[2]
  ellipsedata <- list()
  pointdata <- NULL
  comp <- NULL
  counter <- 1
  ells <- which(x$N > d)
  pointdata <- x$sumXi[x$N < d+1, ]/ x$N[x$N < d+1]
  Nleafs <- length(x$N)
  radius <- attr(x, "radius")
  dotrim <- ifelse(is.na(radius), FALSE, TRUE)

  if (length(ells) > 0)
    for (i in ells){
      edata <- list(means=x$sumXi[i,]/x$N[i],
                    covs=(x$sumXisq[,,i] - x$sumXi[i,] %o% (x$sumXi[i,]/x$N[i]))/(x$N[i]-1))
      ## Confirm that diagonal elements are positive
      diag(edata$covs) <- abs(diag(edata$covs))
      eldata <- ellipse(edata$covs, npoints=50)
      if (dotrim){
        dists <- rep(NA, 50)
        for (i in 1:50)
          dists[i] <- crossprod(eldata[i,])
        eldata[dists > radius , ] <-  sqrt(radius/rowSums(eldata[dists > radius,]^2)) * eldata[dists > radius, ]
      }

      ellipsedata[[counter]] <- sweep(eldata, 2, edata$means[1:2], "+")

      comp <- rbind(comp, ellipsedata[[counter]])
      counter <- counter+1
    }
  comp <- rbind(comp, pointdata)
  dimnames(comp) <- list(NULL, attr(x, "xcolnames"))

  if (missing(col)){
    pcols <- 1
    ecols <- rep(1, length(ells))
  } else if (length(col) > 1){
    pcols <- col[x$N < d+1]
    ecols <- col[x$N > d]
  } else {
    pcols <- col
    ecols <- rep(col, length(ells))
  }

  plot(comp, type="n", ...)

  if (length(ellipsedata) > 0)
    for (i in 1:length(ellipsedata))
      lines(ellipsedata[[i]], col=ecols[i], ...)
  points(pointdata, col=pcols, ...)
}

plot.birch.centers <- function(x, ...) {
  y <- x$sumXi/x$N
  dimnames(y) <- list(NULL, attr(x, "xcolnames"))
  plot(y, ...)
}


pairs.birch <- function(x, centers=FALSE, ...){
  d <- dim(x)[2]
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(mfcol=c(d,d), oma=c(2,2,0,0), mar=c(0,0,0,0)+0.1)
  for (i in 1:d)
    for (j in 1:d){
      if (i == j){
        plot(0, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="", axes=FALSE, type="n")
        text(0,0,ifelse(is.null(attr(x,"colnames")), paste("Var",i), attr(x,"colnames")[i]), fon=2)
        box()
      }
      else if (i > j)
        plot(0, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="", axes=FALSE, type="n")
      else {
        plot(x[,c(i, j)], centers=centers, xlab="", ylab="", axes=FALSE, ...)
        box()
        if (i == 1) ## y-axis
          axis(2)
        if (j == d) ## x-axis
          axis(1)
      }
    }
}

points.birch <- function(x, ...){
  if (is.null(x$sumXi))
    x <- birch.getTree(x)

  y <- x$sumXi/x$N
  plot.xy(xy.coords(y), type="p", ...)
}

lines.birch <- function(x, col=1, ...){
  if (is.null(x$sumXi))
    x <- birch.getTree(x)

  d <- dim(x)[2]
  ellipsedata <- list()
  counter <- 1
  ells <- which(x$N > d)
  radius <- attr(x, "radius")
  dotrim <- ifelse(is.na(radius), FALSE, TRUE)

  for (i in ells){
    edata <- list(means=x$sumXi[i,]/x$N[i],
                  covs=(x$sumXisq[,,i] - x$sumXi[i,] %o% (x$sumXi[i,]/x$N[i]))/(x$N[i]-1))
    ## Confirm that diagonal elements are positive
    diag(edata$covs) <- abs(diag(edata$covs))
    eldata <- ellipse(edata$covs, npoints=50)
    if (dotrim){
      dists <- rep(NA, 50)
      for (i in 1:50)
        dists[i] <- crossprod(eldata[i,])
      eldata[dists > radius , ] <-  sqrt(radius/rowSums(eldata[dists > radius,]^2)) * eldata[dists > radius, ]
    }

    ellipsedata[[counter]] <- sweep(eldata, 2, edata$means[1:2], "+")
    counter <- counter+1
  }

  ecols <- rep(0, length(ells))
  ecols[] <- col ## Will recycle if necessary

  for (i in 1:length(ellipsedata))
      plot.xy(xy.coords(ellipsedata[[i]]), type="l", col=ecols[i], ...)
}

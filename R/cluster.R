dist.birch <- function(birchObject, ...){
  if (is.null(birchObject$sumXi))
    birchObject <- birch.getTree(birchObject)
  return(dist(birchObject$sumXi/birchObject$N, ...))
}

kmeans.birch <- function(birchObject, centers, nstart=1){
  if (class(birchObject) != "birch")
    stop("Not a birch object.")
 if (is.null(birchObject$sumXi))
    birchObject <- birch.getTree(birchObject)
 
  ## Either centers is a (k*d) matrix (in which case it stipulates the starting centers)
  ## or an integer, giving the number of starts
  k <- ifelse(is.matrix(centers), nrow(centers), as.integer(abs(centers)))
  d <- dim(birchObject)[2]
  N <- dim(birchObject)[1]
  Nleafs <- length(birchObject)

  ## Validate arguments
  if (is.matrix(centers)){
    if(nstart > 1)
      warning("Center matrix provided, and nstart > 1. If providing a starting seed, you should only have one start.")
    if (ncol(centers) != d)
      stop("Number of columns in centers does not match number of columns in birch object.")
  }
  if (nstart*k > sum(birchObject$N > 20))
    warning("nstart is too large - not enough subclusters. nstart has been reduced.")
  if ( k < 2 || (!is.matrix(centers) & k - centers != 0))
    stop("Number of centers not an integer > 1.")
  
  ## create data - the trace of each sum-of-squares
  birchObject$trsumXisq <- colSums(apply(birchObject$sumXisq, 3, diag))
  
  if (!is.null(dim(centers))){ ## Centers have been provided
    results <- kmeansBirch.iterate(birchObject, centers, Nleafs, niter=20)
    resids <- results$RSS
    clust <- results$clust
  }
  else {
    index <- which(birchObject$N > 20)
    if (nstart*k >length(index))
      nstart <- floor(length(index)/k)
    startCenters <- matrix(sample(index, nstart*k, replace=F), ncol=k)
    results <- list()
    for (i in 1:nstart)
      results[[i]] <- kmeansBirch.iterate(birchObject, startCenters[i,], Nleafs, niter=20)
    resids <- sapply(results, function(xx) xx$RSS)
    clust <- results[[which.min(resids)]]$clust
  }
  outputs <- list()
  outputs$RSS <- resids[which.min(resids)]
  outputs$clust <- list(sub=clust, obs=rep(NA, N))
  for (i in unique(clust))
    outputs$clust$obs[getMembers(birchObject, which(clust == i))] <- i
  return(outputs)
}


kmeansBirch.iterate <- function(fulltree, centers, Nleafs, niter){
  if (!is.matrix(centers)) ## Have been given index
    centers <- fulltree$sumXi[centers,] /fulltree$N[centers]
  
  dists <- kmeansBirch.calcDist(fulltree, centers, Nleafs)
  for (i in 1:niter){
    membership <- max.col(-dists)
    centers <- kmeansBirch.calcCenters(fulltree, membership)
    dists <- kmeansBirch.calcDist(fulltree, centers, Nleafs)
  }
  return(list(clust=membership, RSS=sum(dists[cbind(1:Nleafs, membership)] )))
}

kmeansBirch.calcDist <- function(x, centers, Nleafs){
  k <- nrow(centers)
  dists <- matrix(x$trsumXisq, ncol=k, nrow=Nleafs)
  for (i in 1:k)
    dists[,i] <- dists[,i] - 2*crossprod(t(x$sumXi), centers[i,]) + x$N *drop(crossprod(centers[i,]))
  return(dists)
}
  

kmeansBirch.calcCenters <- function(x, members){
  k <- length(unique(members))
  centers <- matrix(0, nrow=k, ncol=ncol(x$sumXi))
  for (i in 1:k)
    centers[i,] <- colSums(x$sumXi[members == i,])/sum(x$N[members==i])
  return(centers)
}

`rlga.birch` <- function(birchObject, k, alpha=0.5, nsamp=100){
  if (class(birchObject) != "birch")
    stop("Not a birch object.")
  if (is.null(birchObject$sumXi))
    birchObject <- birch.getTree(birchObject)
  
  d <- dim(birchObject)[2]
  N <- dim(birchObject)[1]
  Nleafs <- length(birchObject)

  birchObject$xtxlong <- birchObject$sumXisq
  attributes(birchObject$xtxlong)$dim <- c(d, d * Nleafs)
  birchObject$xtxlong <- t(birchObject$xtxlong) ## Dimension = (n*d) * d
  
  h <- alpha * sum(birchObject$N) ## h is the number required in the subset
  nList <-  covMcdBirch.startlist(birchObject, Nleafs, d, nsamp*k)

  ## Validate arguments
  if (Nleafs < k*10)
    stop("Not enough subclusters. Make the radius/compact smaller, then try again.")
  if (length(nList) < k)
    stop("Not enough subclusters of full rank for even one start.")
  if (length(nList) < k*nsamp)
    warning("Not enough subclusters of full rank. nsamp has been reduced")
  
  nstarts <- min(nsamp, floor(length(nList)/k))
  starts <- matrix(sample(nList, nstarts*k, replace=FALSE), ncol=k)
  totResids <- rep(NA, nstarts)
  CstepOut <- list()
  
  cat("Progress: ")
  progressChars <- c("\174", "\057", "\055", "\134")
  for (i in 1:nstarts){
    cat(progressChars[(i-1) %% 4 + 1])
    Betas <- rlgaBirch.getBetas(birchObject, split(starts[i,], factor(1:k)))
    CstepOut[[i]] <- rlgaBirch.Cstep(birchObject, params=Betas, h, counter=0)
    totResids[i] <- CstepOut[[i]]$ROSS
    cat("\b")
  }
  best <- which.min(totResids)
  bestData <- CstepOut[[ best ]]
  cat(" done\n")

  outputs <- list()
  outputs$ROSS <- totResids[best]
  outputs$clust <- list()
  outputs$clust$obs <- rep(0, N)
  outputs$clust$sub <- rep(0, Nleafs)
  for (i in 1:k){
    outputs$clust$sub[ bestData$H[[i]] ] <- i
    outputs$clust$obs[ getMembers(birchObject, bestData$H[[i]]) ] <- i
  }
  return(outputs)  
}


`rlgaBirch.Cstep` <-   function(fulltree, params, h, counter, Hold=NULL){
  distances <- rlgaBirch.calcAvgResids(fulltree, params)
  H <- rlgaBirch.getH(fulltree$N, h, distances)
  if (counter > 20 || any(sapply(H, function(xx, yy) sum(yy$N[xx]), yy=fulltree) < 10))
    return(list(params=params,H=H, ROSS=rlgaBirch.getROSS(distances, H, fulltree)))
  else {
    params <- rlgaBirch.getBetas(fulltree, H)
    return(rlgaBirch.Cstep(fulltree, params, h, counter+1, H))
  }
}

`rlgaBirch.calcAvgResids` <- function(fulltree, params){
  return(rlgaBirch.calcResids(fulltree, params)/fulltree$N)
}

`rlgaBirch.calcResids` <- function(fulltree, params){
  d <- ncol(params); k <- nrow(params); n <- length(fulltree$N)
  distances <- matrix(NA, nrow=n, ncol=k)
  for (i in 1:k){
    distances[,i] <- matrix(fulltree$xtxlong %*% params[i, -d], ncol=length(params[i,-1]), byrow=TRUE) %*% params[i, -d] -
      2 * fulltree$sumXi  %*% params[i,-d] * params[i,d] + fulltree$N*crossprod(params[i,d])
  }
  return(distances)
}
`rlgaBirch.getBetas` <- function(fulltree, whichones){
  outputs <- matrix(0, nrow=length(whichones), ncol=(ncol(fulltree$sumXi)+1))
  for (i in 1:length(whichones))
    outputs[i,] <- rlgaBirch.orthreg(getSZbar(fulltree, whichones[[i]]))
  return(outputs)
}

`rlgaBirch.getH` <- function(N, h, distances){
  n <- nrow(distances)
  whichmin <- cbind(1:n, max.col(-distances), N)
  whichmin2 <- whichmin[order(distances[ whichmin[,1:2] ]), , drop=FALSE]
  whichmin3 <- whichmin2[cumsum(whichmin2[,3]) < h, , drop=FALSE]
  whichones <- split(whichmin3[,1], as.factor(whichmin3[,2]))
  return(whichones)
}

`rlgaBirch.getROSS` <- function(distances, best, fulltree){
  totresids <- 0
  for (i in 1:length(best)){
    totresids <- totresids + sum(distances[best[[i]],i]* fulltree$N[ best[[i]] ])
  }
  return(totresids)
}

`rlgaBirch.orthreg` <- function(x){
  es <- svd(x$Sz)$v[,ncol(x$Sz)]
  return(c(es, es %*% x$zbar))
}

`lga.birch` <- function(birchObject, k, nsamp=100)
  return(rlga.birch(birchObject=birchObject, k=k, alpha=1, nsamp=nsamp))

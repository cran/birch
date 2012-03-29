`covMcd.birch` <- function(birchObject, alpha=0.5, nsamp=100){
  if (class(birchObject) != "birch")
    stop("Not a birch object.")
  if (is.null(birchObject$sumXi))
    birchObject <- birch.getTree(birchObject)

  d <- dim(birchObject)[2]
  N <- dim(birchObject)[1]
  Nleafs <- length(birchObject)

  h <- alpha * N  ## h is the number required in the subset
  whichones <- covMcdBirch.startlist(birchObject, Nleafs, d, nsamp)

  ## Validate arguments
  if (length(whichones) < nsamp)
    warning("Not enough subclusters of full rank. nsamp has been reduced.")
  if (alpha < 0.5 || alpha > 1)
    stop("alpha must be in (0.5, 1].")

  detS <- rep(NA, length(whichones))
  CstepOut <- list()

  cat("Progress: ")
  progressChars <- c("\174", "\057", "\055", "\134")
  for (i in seq(along=whichones)){
    cat(progressChars[(i-1) %% 4 + 1])
    P <- getSZbar(birchObject, whichones[i])
    CstepOut[[i]] <- covMcdBirch.Cstep(birchObject, P, h, 0)
    detS[i] <- det(CstepOut[[i]]$params$Sz, log=FALSE)
    cat("\b")
  }
  best <- which.min(detS)
  bestData <- CstepOut[[ best ]]
  cat(" done\n")
  mysubset <- getMembers(birchObject, bestData$H)
  outputs <- list(zbar=bestData$params$zbar, Sz=bestData$params$Sz, Det=detS[best],
              best=list(obs=mysubset, sub=bestData$H))
  class(outputs) <- "covMcd.birch"
  return(outputs)
}

summary.covMcd.birch <- function(object, ...){
  cat("CovMcd using birch\n")
  cat("MCD = ", object$Det)
  cat("\nEstimate of center\n")
  print(object$zbar, ...)
  cat("Estimate of dispersion\n")
  print(object$Sz, ...)
}

covMcdBirch.startlist <- function(fulltree, Nleafs, d, nsamp){
  targetsamp <- min(Nleafs, nsamp)
  startlist <- NULL
  stepthrough <- sample(which(fulltree$N > d + 1), replace=FALSE)
  for (i in stepthrough){
    if (qr(getSZbar(fulltree, i)$Sz)$rank == d)
      startlist <- c(startlist, i)
    if (length(startlist) == targetsamp)
      break
  }
  return(startlist)
}


`covMcdBirch.calcDist` <-
function(fulltree, params){
  return(mahalanobis(fulltree$sumXi/fulltree$N, params$zbar, params$Sz))
}

`covMcdBirch.Cstep` <-
  function(fulltree, params, h, counter){
    counter <- 1
    H <- 1
    Hold <- 0
    while (counter <= 20 && length(intersect(H, Hold)) != max(length(H),length(Hold))){
      D0 <- covMcdBirch.calcDist(fulltree, params)
      Hold <- H
      H <- covMcdBirch.getH(fulltree$N, h, order(D0))
      params <- getSZbar(fulltree, H)
      counter <- counter+1
    }
    if (counter > 20)
      cat("x")
    return(list(params=params,H=H))
  }


`covMcdBirch.getH` <-
  function(N, h,Dorder){
    return(Dorder[ 1:min( which(cumsum(N[Dorder]) >= h)) ])
  }


covMcdBirch.refinement <- function(covOut, x, alpha=0.5){
  ## Validate arguments
  if (!is.matrix(x))
    stop("x must be a matrix")
  if (alpha < 0.5 || alpha > 1)
    stop("alpha must be between 0.5 and 1.")
  if (length(covOut$zbar) != ncol(x))
    stop("Dimensions of covOut and x do not match.")

  best <- covOut$best$obs
  counter <- 0
  oldbest <- NULL
  while (counter <=20 & !identical(oldbest,best)){
    dist <- mahalanobis(x, colMeans(x[best,]), cov(x[best,]))
    oldbest <- best
    best <- which(dist < quantile(dist, alpha))
    counter <- counter+1
  }
  print(counter)
  outputs <- list(Det=det(cov(x[best,])), best=best, zbar=colMeans(x[best,]), Sz=cov(x[best,]))
  class(outputs) <- "covMcd.birch"
  invisible(outputs)
}

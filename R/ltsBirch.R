`lts.birch` <-  function(birchObject, alpha=0.5, intercept=FALSE, nsamp=100){
  if (class(birchObject) != "birch")
    stop("Not a birch object.")
  if (is.null(birchObject$sumXi))
    birchObject <- birch.getTree(birchObject)


  if (intercept) birchObject <- cbind.birch(1, birchObject)
  d <- dim(birchObject)[2]
  N <- dim(birchObject)[1]
  Nleafs <- length(birchObject)

  ## Validate arguments
  if (alpha < 0.5 || alpha > 1)
    stop("alpha must be in (0.5, 1].")

  ## create/amend additional data
  ## Much care is taken if x is a vector, as solve kicks up a fuss
  birchObject$xtx <- birchObject$sumXisq[1:(d-1), 1:(d-1), , drop=FALSE] ## Dimension = (d-1) * (d-1) * n
  birchObject$yty <- birchObject$sumXisq[d, d, ] ## Dimension = n * 1
  birchObject$ytx <- birchObject$sumXisq[1:(d-1), d, ]  ## Dimension = (d-1) * n
  if (is.null(dim(birchObject$ytx)))
    attributes(birchObject$ytx)$dim <- c(1, Nleafs)

  birchObject$xtxlong <- birchObject$xtx
  attributes(birchObject$xtxlong)$dim <- c(d-1, (d-1) * Nleafs)
  birchObject$xtxlong <- t(birchObject$xtxlong) ## Dimension = (n*d) * d
  ## Let's tidy some memory
  birchObject$sumXisq <- birchObject$sumXi <- NULL
  gc()

  h <- alpha * N ## h is the number required in the subset
  whichones <- ltsBirch.startlist(birchObject, Nleafs, d, nsamp)
  if (length(whichones) == 0)
    stop("No sums-of-squares of full rank.Try setting a larger radius/compact")
  if (length(whichones) < nsamp)
    warning("Not enough subclusters of full rank (xtx). nsamp has been reduced.")

  CstepOut <- list()
  totResids <- rep(NA, length(whichones))
  cat("Progress: ")
  progressChars <- c("\174", "\057", "\055", "\134")
  for (i in seq(along=whichones)){
    cat(progressChars[(i-1) %% 4 + 1])
    Betas <- ltsBirch.getBetas(birchObject, whichones[i])
    CstepOut[[i]] <- ltsBirch.Cstep(birchObject, Betas, h, 0)
    totResids[i] <- sum(ltsBirch.calcResids(birchObject, CstepOut[[i]]$params)[ CstepOut[[i]]$H ])
    cat("\b")
  }
  best <- which.min(totResids)
  bestData <- CstepOut[[ best ]]
  cat(" done\n")
  mysubset <- getMembers(birchObject, bestData$H)
  outputs <- list(best=list(leafs=bestData$H, obs=mysubset), raw.coefficients=as.vector(bestData$params), Resids=list(best=totResids[ best ],
                                                                                             whole=sum(ltsBirch.calcResids(birchObject, bestData$params))))
  class(outputs) <- "lts.birch"
  return(outputs)
}

summary.lts.birch <- function(object,...){
  cat("LTS using birch\n")
  cat("LTS = ", object$Resids$best)
  cat("\nCoefficients\n")
  print(object$raw.coefficients, ...)
}


ltsBirch.startlist <- function(fulltree, Nleafs, d, nsamp){
  if (Nleafs < nsamp)
    startlist <- which(apply(fulltree$xtx,3, function(xx) qr(xx)$rank) == (d-1))
  else {
    startlist <- NULL
    stepthrough <- sample(which(fulltree$N > 2), replace=FALSE)
    for (i in stepthrough){
      if (qr(fulltree$xtx[,,i])$rank ==(d-1))
        startlist <- c(startlist, i)
      if (length(startlist) == nsamp)
        break
    }
  }
  return(startlist)
}

`ltsBirch.calcAvgResids` <- function(fulltree, params){
  return(ltsBirch.calcResids(fulltree, params)/fulltree$N)
}

`ltsBirch.calcResids` <- function(fulltree, params){
  return(drop(fulltree$yty - 2* t(params) %*% fulltree$ytx +
              t(params) %*% matrix(fulltree$xtxlong %*% params, nrow=length(params))))
}

`ltsBirch.Cstep` <- function(fulltree, params, h, counter, Hold=NULL){
  D0 <- ltsBirch.calcAvgResids(fulltree, params)
  best <- ltsBirch.getH(fulltree$N, h, order(D0))
  params <- ltsBirch.getBetas(fulltree, best)
  if (counter > 20 || length(intersect(best, Hold)) == max(length(best),length(Hold)))
    return(list(params=params,H=best))
  else
    return(ltsBirch.Cstep(fulltree, params, h, counter+1, best))
}

`ltsBirch.getBetas` <- function(fulltree, whichones){
  if (length(whichones) == 1)
    return(solve(as.matrix(fulltree$xtx[,,whichones])) %*% fulltree$ytx[, whichones])
    else
      return(solve(rowSums(fulltree$xtx[,,whichones, drop=FALSE], dims=2) ) %*% rowSums(fulltree$ytx[, whichones, drop=FALSE]))
}

`ltsBirch.getH` <- function(N, h, Dorder){
  return(Dorder[ 1:min( which(cumsum(N[Dorder]) >= h)) ])
}

ltsBirch.refinement <- function(ltsOut, x, y, alpha=0.5, intercept=FALSE){
  if (is.null(dim(x)))
    x <- as.matrix(x)
  if (intercept)
    x <- cbind(1,x)

  ## Validate arguments
  if (alpha < 0.5 || alpha > 1)
    stop("alpha must be between 0.5 and 1.")
  if (length(ltsOut$raw.coefficients) != ncol(x))
    stop("Dimensions of ltsOut and x do not match.")

  betas <- ltsOut$raw.coefficients
  counter <- 0
  oldbeta <- NULL
  while (counter <= 20 & !identical(oldbeta, betas)){
    res <- (y - x %*% betas)^2
    whichones <- (res < quantile(res, alpha))
    oldbeta <- betas
    betas <- lm.fit(x[whichones, , drop=FALSE], y[whichones])$coeff
    counter <- counter + 1
  }
  print(counter)
  outputs <- list(Resids=list(best=sum(res[whichones])), best=which(whichones), raw.coefficients=betas)
  class(outputs) <- "lts.birch"
  invisible(outputs)
}

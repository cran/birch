birch <- function (x, radius, compact=radius, keeptree=FALSE, ...) {
  ## Argument validation
  if (!is.character(x) & !inherits(x, "connection") & !is.matrix(x))
    stop("Must specify either a matrix, file name or connection.")
  if (radius <= 0 || compact <= 0)
    stop("Compact and/or radius must be greater than zero.")
  if (is.matrix(x))
    if(ncol(x) == 1)
      stop("Must be a matrix of at least two columns.")
  
  birchObject <- list()
  attributes(birchObject) <- list(class="birch", radius=radius, compact=compact)
  if (is.matrix(x)){
    ## Getting data from an R matrix
    if (ncol(x) > 30)
      stop("Matrix can only have up to 30 columns")
    attr(birchObject, "internal") <- .Call("LL_main",  matrix(as.double(x), ncol = ncol(x)),
                                           as.double(radius), as.double(compact), as.integer(1))
    attr(birchObject, "xcolnames") <- colnames(x)
    attr(birchObject, "xdim") <- c(dim(x), .Call("LL_getdim", attr(birchObject, "internal"))[1])
  }
  else {
    ## Getting data from a file/connection
    tmp <- birch.file(x, birchObject, keeptree, ...)
    attr(birchObject, "internal") <- tmp[[1]]
    attr(birchObject, "xcolnames") <- tmp[[2]] ## Gets the column names
    xdim <- .Call("LL_getdim", attr(birchObject, "internal"))
    attr(birchObject, "xdim") <- c(xdim[2], tmp[[3]], xdim[1])
  }
  
  if (!keeptree){
    ## get back the data
    birchObject <- birch.getTree(birchObject)

    ## Kill the tree
    birch.killTree(birchObject)
    attr(birchObject,"internal") <- NULL
  }

  return(birchObject)
}

birch.addToTree <- function(x, birchObject, updateDIM=TRUE, ...){
  ## Argument validation
  if (is.null(attr(birchObject,"internal")))
    stop("No pointer in object. Was tree made with keeptree set to TRUE?")
  if (class(birchObject) != "birch")
    stop("Not a birch object")
  if (!is.character(x) & !inherits(x, "connection") & !is.matrix(x))
    stop("Must specify either a matrix, file name or connection.")
  ## send the data to the tree
  if (is.matrix(x)){
      if (ncol(x) != dim(birchObject)[2])
        stop("Columns in data set do not match those in tree.")
      .Call("LL_adddata", attr(birchObject,"internal"), matrix(as.double(x), ncol = ncol(x)))
  }
  else {
    tmp <- birch.file(x, birchObject, keeptree=TRUE, ...)
  }
  ## In case the output is assigned
  outputs <- birchObject
  outputs$sumXi <- outputs$sumXisq <- outputs$N <- NULL
  if (updateDIM){
    xdim <- .Call("LL_getdim", attr(birchObject, "internal"))
    attr(outputs,"xdim") <- c(xdim[2], dim(birchObject)[2], xdim[1])
  }
  else
    attr(outputs,"xdim") <- c(NA, dim(birchObject)[2], NA)
  return(outputs)
}

birch.getTree <- function(birchObject){
  ## Argument validation
  if (is.null(attr(birchObject,"internal")))
    stop("No pointer in object. Was tree made with keeptree set to TRUE?")
  if (class(birchObject) != "birch")
    stop("Not a birch object")

  ## get the data back from the tree
  birchObject[1:4] <- .Call("LL_getdata", attr(birchObject, "internal"))
  birchObject[[2]] <- t(birchObject[[2]])
  attr(birchObject, "names") <- c("N","sumXi","sumXisq", "members")
  attr(birchObject, "xdim") <- c(sum(birchObject$N), ncol(birchObject$sumXi), length(birchObject$N))
  return(birchObject)
}

birch.killTree <- function(birchObject){
  ## Argument validation
  if (is.null(attr(birchObject,"internal")))
    stop("No pointer in object. Was tree made with keeptree set to TRUE?")
  if (class(birchObject) != "birch")
    stop("Not a birch object")

  ## Kill tree
  .Call("LL_killtree", attr(birchObject,"internal"))
}

birch.file <- function(file, birchObj, keeptree, header, ...){
  ## Get the attributes
  radius <- attr(birchObj, "radius")
  compact <- attr(birchObj, "compact")
  
  ## If the tree exists, then add to it
  addToTree <- !is.null(attr(birchObj, "internal"))

  ## Open file connection
  if (is.character(file)){
    file <- file(file)
    on.exit(close(file))
  }
  if (!inherits(file, "connection"))
    stop("'file' must be a character string or a connection")

  
  cat("Progress: ")
  progressChars <- c("\174", "\057", "\055", "\134")
  counter <- 1
  open(file, "r")

  ## Load first 10 lines
  first <- as.matrix(read.table(file, nrows=10,
                                header=ifelse(missing(header), FALSE, header), ...))
  if (ncol(first) > 30)
    stop("Matrix can only have up to 30 columns")

  if (addToTree){
    pointer <- birchObj
    data <- as.matrix(matrix(as.double(first), ncol = ncol(first)))
    birch.addToTree(data, pointer, updateDIM=FALSE)
  }
  else {
    pointer <- NULL
    attributes(pointer) <- list(class="birch",
                                internal=.Call("LL_main",  matrix(as.double(first), ncol = ncol(first)),
                                  as.double(radius), as.double(compact), as.integer(1)),
                                xdim=c(NA, ncol(first), NA))
  }
  
  ## Then the remainder
  BLOCKSIZE <- floor(1e5/ncol(first)) ## Reads this many rows in at one go
  incomplete <- isDataLeft(file)
  while (incomplete){
    cat(progressChars[(counter-1) %% 4 + 1])
    data <- as.matrix(read.table(file, nrows=BLOCKSIZE, header=FALSE, ...))
    birch.addToTree(data, pointer)
    incomplete <- isDataLeft(file)
    cat("\b"); counter <- counter+1
  }
  cat("done\n")

  return(list(attr(pointer, "internal"), colnames(first), ncol(first)))
}

`getMembers` <- function(fulltree, whichones) {
  return(unlist(fulltree$members[whichones]))
}

`getSZbar` <-
  function(fulltree, whichones){
    n <- sum(fulltree$N[whichones])
    ## Calc new mean
    zbar <- colSums(fulltree$sumXi[whichones,, drop=FALSE])/n
    
    ## New Sd
    Sz <- 1/(n - 1) * (rowSums(fulltree$sumXisq[,, whichones, drop=FALSE], dims=2) - n * zbar %*% t(zbar))
    ## Calc new std dev
    return(list(zbar=zbar, Sz=Sz))
  }

isDataLeft <- function(con){
  oneline <-  readLines(con, n=1)
  if (length(oneline)==0)
    incomplete <- FALSE
  else {
    incomplete <- TRUE
    pushBack(oneline, con)
  }
  return(incomplete)
}

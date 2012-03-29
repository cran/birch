spearman <- function (x, radius, compact=radius, columns = NULL, estimate, ... ) {
  ##Make sure that estimates exist
  possibleEstimate = list('1','2','4', "Student", "student", "Woodbury", "woodbury")
  if (any(!(estimate %in% possibleEstimate)))
    stop ( paste( "\n Invalid choice(s) of estimate(s). \"", as.character(estimate[!(estimate %in% possibleEstimate)]),
                  "\" is not a valid selection. See ?spearman for a list of valid estimates.", sep = ""))
if (class(x) != "birch"){
    ## Argument validation
  if (!is.character(x) & !inherits(x, "connection") & !is.matrix(x))
    stop("Must specify either a matrix, file name or connection.")
  if (radius <= 0 || compact <= 0)
    stop("Compact and/or radius must be greater than zero.")
  if (is.matrix(x))
    if(ncol(x) == 1)
      stop("Must be a matrix of at least two columns.")


  nbEstime <- length(estimate)
  if (is.matrix(x)) {
    #If columns not NULL
    if (!(is.null(columns)))
      stop("Option 'columns' can only be used with file or connection. Use [] instead.")

    #Initialise outputs
    obs = as.double(nrow(x))
    resultMat <- list(Estimate = NULL, nClus = NULL)
    resultMat$estimate <- array(diag(ncol(x)), dim=c(ncol(x),ncol(x),nbEstime),
                                               dimnames = list(colnames(x),colnames(x),paste("Estimate ",estimate,sep = "") ))
    resultMat$nClus <- array(diag(ncol(x)), dim=c(ncol(x),ncol(x)), dimnames = list(colnames(x),colnames(x)) )

    #One birch object for each pair of variables
    for (k in 1:((ncol(x))-1)){
      for (j in (k+1):ncol(x)){
        birch_kj <- birch(x[,c(k,j)],radius,compact,keeptree = TRUE)
        birch_kj <-birch.getTree(birch_kj)
        y <- array(c((birch_kj$sumXi[] / birch_kj$N[]), birch_kj$N),dim = c(dim(birch_kj)[3],  3))


        resultMat$nClus[k,j] <- resultMat$nClus[j,k] <- dim(birch_kj)[3]
        birch.killTree(birch_kj)
        gc()

        #Calculer nb. observations dans les rangs inférieurs pour chaque x
        y <- calculrisi(y)

        for ( l in 1:nbEstime){
            if (1 ==  estimate[l])
              resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime1(y,obs)
            else if ((2 ==  estimate[l]) || ("Student" ==  estimate[l]) || ("student" ==  estimate[l]))
              resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime2(y,obs)
            else if ((4 ==  estimate[l]) || ("Woodbury" ==  estimate[l]) || ("woodbury" ==  estimate[l]))
              resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime4(y,obs)
        }
      }
    }
  }
  else { #x is a file or connection
    ## Need to obtain number of columns and names
    ## Open file connection
    if (is.character(x)){
      file <- file(x)
      on.exit(close(file))
    }
    if (!inherits(file, "connection"))
      stop("'file' must be a character string or a connection")

    #cat("Progress: ")
    progressChars <- c("\174", "\057", "\055", "\134")
    open(file, "r")

    ## Load first 2 lines
    first <- as.matrix(read.table(file, nrows=2, ...))
#    first <- as.matrix(read.table(file, nrows=2,
#                                  header=ifelse(missing(header), FALSE, header), ...))
    ##Retreive columns
    columnsNames <- colnames(first)
    if (!(is.null(columns))){
      columns <- seq(ncol(first))[columns]
      columnsNames <- columnsNames[columns]
    }
    ncolumns <- length(columns)

    ##Initialise outputs
    resultMat <- list(Estimate = NULL, nClus = NULL)
    resultMat$estimate <- array(diag(ncolumns), dim=c(ncolumns,ncolumns,nbEstime),
                                                dimnames = list(columnsNames,columnsNames,paste("Estimate ",estimate,sep = "") ))
    resultMat$nClus <- array(diag(ncolumns), dim=c(ncolumns,ncolumns), dimnames = list(columnsNames,columnsNames) )

    #One birch object for each pair of variables
    for (k in 1:(ncolumns-1)){
      for (j in (k+1):ncolumns){
        if (!(is.null(columns))) birch_kj <- birch(x, columns = columns[c(k,j)],radius,compact,keeptree = TRUE,...)
        else birch_kj <- birch(x, columns = c(k,j),radius,compact,keeptree = TRUE,...)
        birch_kj <-birch.getTree(birch_kj)
        y <- array(c((birch_kj$sumXi[] / birch_kj$N[]), birch_kj$N),dim = c(dim(birch_kj)[3],  3))

        resultMat$nClus[k,j] <- resultMat$nClus[j,k] <- dim(birch_kj)[3]
        obs = as.double(dim(birch_kj)[1])

        birch.killTree(birch_kj)
        gc()

        #Calculer nb. observations dans les rangs inférieurs pour chaque x
        y <- calculrisi(y)

        for ( l in 1:nbEstime){
          if (1 ==  estimate[l])
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime1(y,obs)
          else if ((2 ==  estimate[l]) || ("Student" ==  estimate[l]) || ("student" ==  estimate[l]))
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime2(y,obs)
          else if ((4 ==  estimate[l]) || ("Woodbury" ==  estimate[l]) || ("woodbury" ==  estimate[l]))
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime4(y,obs)
        }
      }
    }
  }

} #else birch object
else{
  #If columns not NULL
  if (!(is.null(columns)))
    stop("Option 'columns' can only be used with file or connection. Use [] instead.")
  resultMat <- birch.spearman(x, estimate,...)
}
return(resultMat)
}


##Calculate the correlation matrix from a birchobj (Spearman's rho)
birch.spearman <- function (birchObj, estimate, ... ) {
  nbEstime = length(estimate)

  #Retreive columns
  columnsNames <- attr(birchObj, "xcolnames")
  ncolumns <- dim(birchObj)[2]
  print("ncolumns")
  print(ncolumns)

  #Initialise outputs
  obs <- as.double(dim(birchObj)[1])

  resultMat <- list(estimate = NULL, nClus = NULL)
  resultMat$estimate <- array(diag(ncolumns), dim=c(ncolumns,ncolumns,nbEstime),
                                               dimnames = list(columnsNames,columnsNames,paste("Estimate ",estimate,sep = "") ))
  resultMat$nClus <- dim(birchObj)[3]

  #One estimate for each pair of variables
  y <- array(c((birchObj$sumXi[] / birchObj$N[]), birchObj$N),dim = c(dim(birchObj)[3],  (ncolumns +1)))
  y <- birch.calculri(y,1, (ncolumns +1))
  for (k in 2:ncolumns){
    ##Calculate ranks
    y <- birch.calculri(y,k,(ncolumns+1))
    ##Calculate estimates
    for ( l in 1:nbEstime){
      if (1 ==  estimate[l])
        resultMat$estimate[k,1,l] <- resultMat$estimate[1,k,l] <-spearman.estime1(y[,c(1,k,(ncolumns+1))],obs)
      else if ((2 ==  estimate[l]) || ("Student" ==  estimate[l]) || ("student" ==  estimate[l]))
        resultMat$estimate[k,1,l] <- resultMat$estimate[1,k,l] <-spearman.estime2(y[,c(1,k,(ncolumns+1))],obs)
      else if ((4 ==  estimate[l]) || ("Woodbury" ==  estimate[l]) || ("woodbury" ==  estimate[l]))
        resultMat$estimate[k,1,l] <- resultMat$estimate[1,k,l] <-spearman.estime4(y[,c(1,k,(ncolumns+1))],obs)
    }
  }
  for (j in 2:(ncolumns-1)){
    for (k in (j+1):ncolumns){
      ##Ranks already exists
      ##Calculate estimates
      for ( l in 1:nbEstime){
        if (1 ==  estimate[l])
          resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime1(y[,c(j,k,(ncolumns+1))],obs)
        else if ((2 ==  estimate[l]) || ("Student" ==  estimate[l]) || ("student" ==  estimate[l]))
          resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime2(y[,c(j,k,(ncolumns+1))],obs)
        else if ((4 ==  estimate[l]) || ("Woodbury" ==  estimate[l]) || ("woodbury" ==  estimate[l]))
          resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-spearman.estime4(y[,c(j,k,(ncolumns+1))],obs)
      }
    }
  }
  return(resultMat)
}




calculrisi <-function(y)
{
    #Calculer les "rangs" (ri) pour x1
    y=y[sort.list(y[,1], method = "quick", na.last=NA), ]
    y[,1] = .Call("LL_boucleSpearman", matrix(as.double(y[,3]),ncol=1))

    #Calculer les "rangs" (si) pour x2
    y=y[sort.list(y[,2], method = "quick", na.last=NA), ]
    y[,2]= .Call("LL_boucleSpearman", matrix(as.double(y[,3]),ncol=1))
    invisible(y)
}

birch.calculri <-function(y,k,max)
{
    #Calculer les "rangs" (ri) pour tous xi
    y=y[sort.list(y[,k], method = "quick", na.last=NA), ]
    y[,k] = .Call("LL_boucleSpearman", matrix(as.double(y[,max]),ncol=1))
    invisible(y)
}

spearman.estime1 <-function(y,obs)
{
    print(dim(y))
    #Option 1 : Toutes les paires
    rho1 = 1-(6*(sum(( y[,3] * (y[,1]-y[,2])^2)))/(obs * (obs ^2 -1)))
    invisible(rho1)
}

#Student
spearman.estime2 <-function(y,obs)
{
    #Option 2 : Juste les paires inter-clusters
    rho2 = 1-(6*(sum(( y[,3] * (y[,1]-y[,2])^2)))/( (obs*(obs ^2 -1))-sum(y[,3]*( y[,3]^2-1)) ))
    invisible(rho2)
}


spearman.estime4 <-function(y,obs)
{
    #Option 4 :
    sumNi = sum(y[,3]*( y[,3]^2-1))
    sumN = (obs * (obs ^2 -1))
    rho4 = 1 - (((6*(sum(( y[,3] * (y[,1]-y[,2])^2)))) + sumNi) / sumN)
    invisible(rho4)
}



kendall <- function (x, radius, compact=radius, columns = NULL, estimate,  ...) {
 ##Make sure that estimates exist
  possibleEstimate = list('1','2','3', "Student", "student", "Woodbury", "woodbury", "Stuart", "stuart")
  if (any(!(estimate %in% possibleEstimate)))
    stop ( paste( "\n Invalid choice(s) of estimate(s). \"", as.character(estimate[!(estimate %in% possibleEstimate)]),
                  "\" is not a valid selection. See ?kendall for a list of valid estimates.", sep = ""))

 if (class(x) != "birch"){

  ## Argument validation
  if (!is.character(x) & !inherits(x, "connection") & !is.matrix(x))
    stop("Must specify either a matrix, file name or connection.")
  if (radius <= 0 || compact <= 0)
    stop("Compact and/or radius must be greater than zero.")
  if (is.matrix(x))
    if(ncol(x) == 1)
      stop("Must be a matrix of at least two columns.")

  nbEstime = length(estimate)

  if (is.matrix(x)) {
    #If columns not NULL
    if (!(is.null(columns)))
      stop("Option 'columns' can only be used with file or connection. Use [] instead.")

    #Initialise outputs
    obs = as.double(nrow(x))
    resultMat <- list(estimate = NULL, nClus = NULL)
    resultMat$estimate <- array(diag(ncol(x)), dim=c(ncol(x),ncol(x),nbEstime),
                                               dimnames = list(colnames(x),colnames(x),paste("Estimate ",estimate,sep = "") ))
    resultMat$nClus <- array(diag(ncol(x)), dim=c(ncol(x),ncol(x)), dimnames = list(colnames(x),colnames(x)) )


    #One birch object for each pair of variables
    for (k in 1:(ncol(x)-1)){
      for (j in (k+1):ncol(x)){
        birch_kj <- birch(x[,c(k,j)],radius,compact,keeptree = TRUE)
        birch_kj <-birch.getTree(birch_kj)
        y <- array(c((birch_kj$sumXi[] / birch_kj$N[]), birch_kj$N),dim = c(dim(birch_kj)[3],  3))

        resultMat$nClus[k,j] <- resultMat$nClus[j,k] <- dim(birch_kj)[3]
        birch.killTree(birch_kj)
        gc()

        y<-y[sort.list(y[,1], method = "quick", na.last=NA), ]

        for ( l in 1:nbEstime){
          if ((1 ==  estimate[l]) || ("Student" ==  estimate[l]) || ("student" ==  estimate[l]))
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime1(y[,-1])
          else if ((2 ==  estimate[l]) || ("Woodbury" ==  estimate[l]) || ("woodbury" ==  estimate[l]))
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime2(y[,-1],obs)
          else if ((3 ==  estimate[l]) || ("Stuart" ==  estimate[l]) || ("stuart" ==  estimate[l]))
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime3(y[,-1],obs)
        }
      }
    }
  }
  else { #x is a file or connection
    ## Need to obtain number of columns and names
    ## Open file connection
    if (is.character(x)){
      file <- file(x)
      on.exit(close(file))
    }
    if (!inherits(file, "connection"))
      stop("'file' must be a character string or a connection")

    #cat("Progress: ")
    progressChars <- c("\174", "\057", "\055", "\134")
    open(file, "r")

    ## Load first 2 lines
    first <- as.matrix(read.table(file, nrows=2, ...))
#    first <- as.matrix(read.table(file, nrows=2,
#                                  header=ifelse(missing(header), FALSE, header), ...))

    ##Retreive columns
    columnsNames <- colnames(first)
    if (!(is.null(columns))){
      columns <- seq(ncol(first))[columns]
      columnsNames <- columnsNames[columns]
    }
    ncolumns <- length(columns)

    resultMat <- list(estimate = NULL, nClus = NULL)
    resultMat$estimate <- array(diag(ncolumns), dim=c(ncolumns,ncolumns,nbEstime),
                                                dimnames = list(columnsNames,columnsNames,paste("Estimate ",estimate,sep = "") ))
    resultMat$nClus <- array(diag(ncolumns), dim=c(ncolumns,ncolumns), dimnames = list(columnsNames,columnsNames) )

    ##Create first birch object to obtain obs
    birch_kj <- birch(x, columns = c(1,2),radius,compact,keeptree = TRUE,...)
    birch_kj <-birch.getTree(birch_kj)
    obs = as.double(dim(birch_kj)[1])

    #One birch object for each pair of variables
    for (k in 1:(ncolumns-1)){
      for (j in (k+1):ncolumns){
        if (!(is.null(columns))) birch_kj <- birch(x, columns = columns[c(k,j)],radius,compact,keeptree = TRUE,...)
        else birch_kj <- birch(x, columns = c(k,j),radius,compact,keeptree = TRUE,...)
        birch_kj <-birch.getTree(birch_kj)
        y <- array(c((birch_kj$sumXi[] / birch_kj$N[]), birch_kj$N),dim = c(dim(birch_kj)[3],  3))

        resultMat$nClus[k,j] <- resultMat$nClus[j,k] <- dim(birch_kj)[3]
        birch.killTree(birch_kj)
        gc()
        y<-y[sort.list(y[,1], method = "quick", na.last=NA), ]

        for ( l in 1:nbEstime){
          if ((1 ==  estimate[l]) || ("Student" ==  estimate[l]) || ("student" ==  estimate[l]))
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime1(y[,-1])
          else if ((2 ==  estimate[l]) || ("Woodbury" ==  estimate[l]) || ("woodbury" ==  estimate[l]))
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime2(y[,-1],obs)
         else if ((3 ==  estimate[l]) || ("Stuart" ==  estimate[l]) || ("stuart" ==  estimate[l]))
            resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime3(y[,-1],obs)
        }
      }
    }
  }
}  #else birch object
else{
  #If columns not NULL
  if (!(is.null(columns)))
    stop("Option 'columns' can only be used with file or connection. Use [] instead.")
  resultMat <- birch.kendall(x, estimate,...)
}
return(resultMat)
}



##Calculate the correlation matrix from a birchobj (Kendall's tau)
birch.kendall <- function (birchObj, estimate, ... ) {
  nbEstime = length(estimate)

  #Retreive columns
  columnsNames <- attr(birchObj, "xcolnames")
  ncolumns <- dim(birchObj)[2]

  #Initialise outputs
  obs <- as.double(dim(birchObj)[1])

  resultMat <- list(estimate = NULL, nClus = NULL)
  resultMat$estimate <- array(diag(ncolumns), dim=c(ncolumns,ncolumns,nbEstime),
                                              dimnames = list(columnsNames,columnsNames,paste("Estimate ",estimate,sep = "") ))
  resultMat$nClus <- dim(birchObj)[3]

  #One estimate for each pair of variables
  y <- array(c((birchObj$sumXi[] / birchObj$N[]), birchObj$N),dim = c(dim(birchObj)[3],  (ncolumns +1)))

  for (k in 1:(ncolumns-1)){
    y<-y[sort.list(y[,k], method = "quick", na.last=NA), ]
    for (j in (k+1):ncolumns){
      for ( l in 1:nbEstime){
        if ((1 ==  estimate[l]) || ("Student" ==  estimate[l]) || ("student" ==  estimate[l]))
          resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime1(y[,c(j,(ncolumns+1))])
        else if ((2 ==  estimate[l]) || ("Woodbury" ==  estimate[l]) || ("woodbury" ==  estimate[l]))
          resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime2(y[,c(j,(ncolumns+1))],obs)
        else if ((3 ==  estimate[l]) || ("Stuart" ==  estimate[l]) || ("stuart" ==  estimate[l]))
          resultMat$estimate[k,j,l] <- resultMat$estimate[j,k,l] <-kendall.estime3(y[,c(j,(ncolumns+1))],obs)
      }
    }
  }
return(resultMat)
}


#Student
kendall.estime1 <-function(y)
{
    resultCall = .Call("LL_boucleKendall",matrix(as.double(y[,1]),ncol=1),matrix(as.double(y[,2]),ncol=1))
    estime1 = (( 2* resultCall[1]) / (resultCall[2])) - 1
    invisible(estime1)
}

#Woodbury
kendall.estime2 <-function(y,obs)
{
    resultCall = .Call("LL_boucleKendall",matrix(as.double(y[,1]),ncol=1),matrix(as.double(y[,2]),ncol=1))
    estime2 = 2* (( 2* resultCall[1]) - resultCall[2]) / (obs*( obs-1))
    invisible(estime2)
}


#Stuart
kendall.estime3 <-function(y,obs)
{
    resultCall <- .Call("LL_boucleKendall",matrix(as.double(y[,1]),ncol=1),matrix(as.double(y[,2]),ncol=1))
    est4 <- (2*dim(y)[1]*( 2*resultCall[1] - resultCall[2]))/((obs*obs)*(dim(y)[1]-1))
    invisible(est4)
}

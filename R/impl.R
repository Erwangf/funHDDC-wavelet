#' funHDDCwavelet : time series clustering using wavelet and gaussian mixture models
#'
#' This package provide an implementation of the funHDDCwavelet algorithm.
#' This algorithm allow the clustering of time series by representing them in a
#' parcimonious and multi-resolution way.
#'
#' @docType package
#' @name funHDDCwavelet
#'
#' @import wavethresh
#' @import stats
#' @import graphics
NULL


weightedCovariance <- function(D,ratio,mu){
  if (is.null(mu)) {
    mu = colMeans(ratio * D)
  }

  ones = rep(1, nrow(D))
  tmp = sqrt(ratio) * (D - ones %*% t(mu))
#  result <- crossprod(tmp) / sum(ratio)
  return(crossprod(tmp) / sum(ratio))
}

.onLoad = function(libname, pkgname){
  packageStartupMessage(
  cat("   __             _   _ ____  ____   ____                         _      _   \n  / _|_   _ _ __ | | | |  _ \\|  _ \\ / ___|_      ____ ___   _____| | ___| |_ \n | |_| | | | '_ \\| |_| | | | | | | | |   \\ \\ /\\ / / _` \\ \\ / / _ \\ |/ _ \\ __|\n |  _| |_| | | | |  _  | |_| | |_| | |___ \\ V  V / (_| |\\ V /  __/ |  __/ |_ \n |_|  \\__,_|_| |_|_| |_|____/|____/ \\____| \\_/\\_/ \\__,_| \\_/ \\___|_|\\___|\\__|\n
                                                                             \n")
  )
}

getBeginEnd <- function(l){
  # level 1 = c0 => 1:1
  if (l == 1) {
    begin = 1
    end = 1
  } else {
    # level 2 = d0 => 2:2
    # level 3 = d1 => 3:4 etc...
    begin = 2^(l - 2) + 1
    end = 2^(l - 1)
  }
  return(list(begin = begin, end = end))
}

decomposeByBlock <- function(Q){
  Qd = list()
  Ad = list()
  size = ncol(Q)
  lvls = 1:(log2(size + 1) + 1)


  for (l in lvls) {
    bounds = getBeginEnd(l)
    begin = bounds$begin
    end = bounds$end

    levelData = Q[begin:end,begin:end]
    eigenData = eigen(levelData)

    Qd[[l]] <- eigenData$vectors
    Ad[[l]] <- eigenData$values
  }
  return(list(Qd = as.matrix(Qd), Ad = Ad))
}

splitCovarByBlock <- function(W){
  Wd <- list()
  size = ncol(W)
  lvls = 1:(log2(size + 1) + 1)


  for (l in lvls) {
    bounds = getBeginEnd(l)
    begin = bounds$begin
    end = bounds$end

    Wd[[l]] <- W[begin:end,begin:end]
  }
  return(Wd)
}

Pr = function(x,A){
  return(t(x) %*% A %*% x)
}

toWavelet <- function(x,wavelet.family,wavelet.filter.number,max.level=round(log2(length(x)))){
  size = length(x)
  if (round(log2(size)) != log2(size)) {
    size = 2^round(log2(size))
    x = spline(x,n = size)$y
  }

  xw = wd(x,family = wavelet.family,filter.number = wavelet.filter.number)
  return(c(xw$C[1],xw$D[1:(2^max.level - 1)]))
}

kCost <- function(x,mu,a,b,d,Pi,QtQ,A){
  if (sum(a <= 0) > 0) {
    # print("A<0")
    # a[a<=0]=0.0000001
  }

  P = length(x)

  Px <- QtQ %*% (x - mu ) + mu

  K = Pr(mu - Px,A) +  sum((x - Px)^2)/b + sum(log(a[1:d])) + (P - d) * log(b) - 2 * log(Pi)

  return(K)
}

kCostWithoutB <- function(x,mu,a,Pi,QtQ,A){
  if (sum(a <= 0) > 0) {
    print("A<0")
  }
#  P = length(x)
  Px <- QtQ %*% (x - mu ) + mu
  K = Pr(mu - Px,A) +  sum(log(a))  - 2 * log(Pi)
  return(as.numeric(K))
}

computeW <- function(X,tm,K,mu,Nis){
  W <- list()
#  N <- nrow(X)
  for (g in 1:K) {
    W[[g]] <- weightedCovariance(X,tm[,g],mu[g,])
  }
  return(W)
}

computeLogLikelihood <- function(Q_all,A_all,Pi,n,W_full,D){

  ll <- 0 # logLikelihood
  k = length(A_all) # number of cluster
  maxLevel = ncol(D)
  for (i in 1:k) { # for each cluster :
    ll_i <- 0 # current cluster likelihood
    Wi <- splitCovarByBlock(W_full[[i]])
    Qi <- Q_all[[i]] # Qi is a list containing for each level a matrix (eigenvectors)
    Ai <- A_all[[i]] # Ai is a list containing for each level a vector of eigenvalues

    Pi_i <- Pi[i] # cluster repartition

    # for each level in the cluster :
    for (lvl in 1:maxLevel) {
 #     sumA <- 0
      Q <- Qi[[lvl]] # Matrix (eigenvectors)
      W <- Wi[[lvl]] # Matrix (empirical covariance)
      a <- Ai[[lvl]] # Vector (eigenvalues)


      # for the first d_{i,lvl} columns of the level
      for (l in 1:D[i,lvl]) {
        ll_i <- ll_i + log(a[l]) + (t(Q[,l]) %*% W %*% Q[,l]) / a[l]
      }
      # for the other columns of the level
      if (D[i,lvl] + 1 <= ncol(Q)) {
        for (l in (D[i,lvl] + 1):ncol(Q)) {
          ll_i <- ll_i + log(mean(a)) + (t(Q[,l]) %*% W %*% Q[,l]) / mean(a)
        }
      }

    }
    ll <- ll + n * Pi_i * ll_i
  }
  return(-ll / 2)
}

getIntrinsecDim <- function(lambda,n,test="scree", value = 0.1){
  lambda <- lambda[lambda >= 0]
  if (length(lambda) == 1) {
    return(1)

  }  else if (test == "mean") {
    meanLambda = mean(lambda)
    return(length(lambda[lambda > meanLambda]))

  } else if (test == "screeBIS") {
    # Scree Test de Cattel
    scree_mat <- matrix(0, ncol = length(lambda), nrow = length(lambda))
    scree_mat[,1] <- lambda
    for (i in 2:length(lambda)) {
      for (j in 2:i) {
        scree_mat[i,j] <- scree_mat[i - 1, j - 1] - scree_mat[i,j - 1]
        if (scree_mat[i,j] < 0) {
          return(i - 1)
        }
      }
    }
    return(length(lambda))

  } else if (test == "scree") {
    lambda <-  (lambda - min(lambda) )/(max(lambda) - min(lambda))
    dimMax = length(lambda)
    i <- 1
    found = FALSE

    while (i < dimMax && found == FALSE) {

      if ((lambda[i] - lambda[i + 1]) <= value) {
        found <- TRUE
      }
      else {
        i <- i + 1
      }
    }

    return(i)
  }

  else if (test == "kss") {

    lambda <- lambda / mean(lambda)

    thr <- 1 + 2 * sqrt((length(lambda) - 1) / (n - 1))
    return(length(lambda[lambda > thr]))


  } else {
    stop(paste("Error",test,"is not a valid test."))
  }
}

getBalancedCluster <- function(clusters,K,minPerClass){
  cl <- list()
  nbPerClass <- c()
  for (g in 1:K) {
    cl[[g]] <- which(clusters == g)
    nbPerClass <- c(nbPerClass,length(cl[[g]]))
  }
  targetClasses <- which(nbPerClass < minPerClass)

  for (tc in targetClasses) {

    while (nbPerClass[[tc]] < minPerClass) {
      # find classes with nbPerClass[[tc]] > minPerClass
      cls <- which(nbPerClass > minPerClass)

      # select one class randomly
      giver <- cls[sample(1:length(cls),1)]

      # add one individual from giver class to target class
      nbToMove <- cl[[giver]][sample(1:length(cl[[giver]]),1)]
      cl[[giver]] <- setdiff(cl[[giver]],nbToMove)
      cl[[tc]] <- c(cl[[tc]],nbToMove)
      nbPerClass[[giver]] <- nbPerClass[[giver]] - 1
      nbPerClass[[tc]] <- nbPerClass[[tc]] + 1
    }
  }

  # convert list[[cluster]]{lineNumber} to c(cluster)
  result <- c()
  for (i in 1:length(clusters)) {
    cl_val <- 0
    j = 0
    while (cl_val == 0 && j < K) {
      j <- j + 1
      if (i %in% cl[[j]]) cl_val <- j
    }
    result <- c(result,cl_val)
  }
  return(result)
}


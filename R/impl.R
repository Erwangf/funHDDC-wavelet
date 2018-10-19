

#' funHDDCwavelet : time series clustering using wavelet and gaussian mixture models
#'
#' This package provide an implementation of the funHDDCwavelet algorithm. This algorithm allow the clustering of time series
#' by representing them in a parcimonious and multi-resolution way.
#'
#' @docType package
#' @name funHDDCwavelet
#'
#' @import wavethresh
#' @import stats
#' @import graphics
NULL




weightedCovariance <- function(D,ratio,mu){
  if(is.null(mu)){
    mu = colMeans(ratio * D)
  }

  ones = rep(1, nrow(D))
  tmp = sqrt(ratio) * (D - ones %*% t(mu))
  result <- crossprod(tmp) / sum(ratio)
  return(crossprod(tmp) / sum(ratio))
}

.onLoad = function(libname, pkgname){
  cat("   __             _   _ ____  ____   ____                         _      _   \n  / _|_   _ _ __ | | | |  _ \\|  _ \\ / ___|_      ____ ___   _____| | ___| |_ \n | |_| | | | '_ \\| |_| | | | | | | | |   \\ \\ /\\ / / _` \\ \\ / / _ \\ |/ _ \\ __|\n |  _| |_| | | | |  _  | |_| | |_| | |___ \\ V  V / (_| |\\ V /  __/ |  __/ |_ \n |_|  \\__,_|_| |_|_| |_|____/|____/ \\____| \\_/\\_/ \\__,_| \\_/ \\___|_|\\___|\\__|\n
                                                                             \n")
}


getBeginEnd <- function(l){
  # level 1 = c0 => 1:1
  if(l == 1){
    begin = 1
    end = 1
  } else {
    # level 2 = d0 => 2:2
    # level 3 = d1 => 3:4 etc...
    begin = 2^(l - 2) +1
    end = 2^(l - 1)
  }
  return(list(begin=begin,end=end))
}

decomposeByBlock <- function(Q){
  Qd = list()
  Ad = list()
  size = ncol(Q)
  lvls = 1:(log2(size+1)+1)


  for(l in lvls){
    bounds = getBeginEnd(l)
    begin = bounds$begin
    end = bounds$end

    levelData = Q[begin:end,begin:end]
    eigenData = eigen(levelData)

    Qd[[l]] <- eigenData$vectors
    Ad[[l]] <- eigenData$values
  }
  return(list(Qd=as.matrix(Qd),Ad=Ad))
}

splitCovarByBlock <- function(W){
  Wd <- list()
  size = ncol(W)
  lvls = 1:(log2(size+1)+1)


  for(l in lvls){
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

toWavelet <- function(x,wavelet.family,wavelet.filter.number){
  size = length(x)
  if(round(log2(size)) != log2(size)){
    size = 2^round(log2(size))
    x = spline(x,n=size)$y
  }

  xw = wd(x,family = wavelet.family,filter.number = wavelet.filter.number)
  return(c(xw$C[1],xw$D))
}

kCost <- function(x,mu,a,b,d,Pi,QtQ,A){
  P = length(x)

  Px <- QtQ %*% ( x - mu ) + mu

  K = Pr(mu-Px,A) +  sum((x-Px)^2)/b + sum(log(a[1:d])) + (P - d) * log(b) - 2 * log(Pi)

  return(K)
}

kCostWithoutB <- function(x,mu,a,Pi,QtQ,A){
  if(sum(a<=0)>0){
    print("HELLO")
  }
  P = length(x)
  Px <- QtQ %*% ( x - mu ) + mu
  K = Pr(mu-Px,A) +  sum(log(a))  - 2 * log(Pi)
  return(as.numeric(K))
}

computeW <- function(X,tm,K,mu,Nis){
  W <- list()
  N <- nrow(X)
  for(g in 1:K){
    W[[g]] <- weightedCovariance(X,tm[,g],mu[g,])
  }
  return(W)
}

computeLogLikelihood <- function(Q_all,A_all,Pi,n,W_full,D){

  ll <- 0 # logLikelihood
  k = length(A_all) # number of cluster
  maxLevel = ncol(D)
  for(i in 1:k){ # for each cluster :
    ll_i <- 0 # current cluster likelihood
    Wi <- splitCovarByBlock(W_full[[i]])
    Qi <- Q_all[[i]] # Qi is a list containing for each level a matrix (eigenvectors)
    Ai <- A_all[[i]] # Ai is a list containing for each level a vector of eigenvalues

    Pi_i <- Pi[i] # cluster repartition

    # for each level in the cluster :
    for(lvl in 1:maxLevel){
      sumA <- 0
      Q <- Qi[[lvl]] # Matrix (eigenvectors)
      W <- Wi[[lvl]] # Matrix (empirical covariance)
      a <- Ai[[lvl]] # Vector (eigenvalues)


      # for the first d_{i,lvl} columns of the level
      for(l in 1:D[i,lvl]){
        ll_i <- ll_i + log(a[l]) +(t(Q[,l]) %*% W %*% Q[,l]) / a[l]
      }
      # for the other columns of the level
      if(D[i,lvl]+1<=ncol(Q)){
        for(l in (D[i,lvl]+1):ncol(Q)){
          ll_i <- ll_i + log(mean(a)) +(t(Q[,l]) %*% W %*% Q[,l]) / mean(a)
        }
      }

    }
    ll <- ll + n * Pi_i * ll_i
  }
  return(- ll / 2)
}

getIntrinsecDim <- function(lambda,n,test="scree"){
  lambda <- lambda[lambda>=0]
  if(length(lambda)==1){
    return(1)

  } else if(substr(test, nchar(test), nchar(test))=="%"){
    # format : 80%, 78.2%, 85.52555%
    maxPercentage= as.numeric(substr(test, 1, nchar(test)-1)) / 100
    total = 0
    d=0
    while(total<maxPercentage){
      d <- d + 1
      total <- total + lambda[d]/sum(lambda)
    }
    return(d)

  } else if(test=="mean"){
    meanLambda = mean(lambda)
    return(length(lambda[lambda>meanLambda]))

  } else if(test=="scree"){
    # Scree Test de Cattel
    scree_mat <- matrix(0,ncol=length(lambda),nrow=length(lambda))
    scree_mat[,1] <- lambda
    for(i in 2:length(lambda)){
      for(j in 2:i){
        scree_mat[i,j] <- scree_mat[i-1,j-1] - scree_mat[i,j-1]
        if(scree_mat[i,j]<0){
          return(i-1)
        }
      }
    }
    return(length(lambda))

  } else if(test=="kss"){

    lambda <- lambda / mean(lambda)

    thr <- 1 + 2 * sqrt((length(lambda)-1)/(n-1))
    return(length(lambda[lambda>thr]))


  } else {
    stop(paste("Error",test,"is not a valid test."))
  }
}

getBalancedCluster <- function(clusters,K,minPerClass){
  cl <- list()
  nbPerClass <- c()
  for(g in 1:K){
    cl[[g]] <- which(clusters==g)
    nbPerClass <- c(nbPerClass,length(cl[[g]]))
  }
  targetClasses <- which(nbPerClass<minPerClass)

  for(tc in targetClasses){

    while(nbPerClass[[tc]]<minPerClass){
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
  for(i in 1:length(clusters)){
    cl_val <- 0
    j = 0
    while(cl_val == 0 && j<K){
      j <- j + 1
      if(i %in% cl[[j]]) cl_val <- j
    }
    result <- c(result,cl_val)
  }
  return(result)
}

#' Adjusted Rand Index
#'
#' Compute the Adjusted Rand Index for two set of classes
#'
#' @param x a vector of classes
#' @param y a vector of classes
#'
#' @examples
#' A = c(1,1,1,2,2,2)
#' B = c(1,2,1,2,2,1)
#' adjustedRandIndex(A,B)
#' @export adjustedRandIndex
adjustedRandIndex <- function(x,y){

    x <- as.vector(x)
    y <- as.vector(y)
    # error control : length of vectors x & y
    if (length(x) != length(y)) {
      stop("Error : arguments must be vectors of the same length")
    }

    tab <- table(x, y)
    if (all(dim(tab) == c(1, 1))){
      return(1)
    }

    a <- sum(choose(tab, 2))
    b <- sum(choose(rowSums(tab), 2)) - a
    c <- sum(choose(colSums(tab), 2)) - a
    d <- choose(sum(tab), 2) - a - b - c
    ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
    return(ARI)

}


#' Generate a time series dataset
#'
#' The dataset composed by 300 curves of 256 time points, from 3 different groups.
#'
#' @return \describe{
#' \item{X}{a matrix of 300 time series (300 rows and 256 columns)}
#' \item{col}{a vector, describing the class of each row (i.e. curve)}
#' }
#'
#' @export generateDataset
generateDataset = function(){
  cluster1.aMean = 0.9
  cluster1.aSd = 0.035
  cluster1.bMean = 0.4
  cluster1.bSd = 0.035

  cluster2.aMean = 0.7
  cluster2.aSd = 0.035
  cluster2.bMean = 0.3
  cluster2.bSd = 0.035

  cluster3.aMean = 0.5
  cluster3.aSd = 0.035
  cluster3.bMean = 0.5
  cluster3.bSd = 0.035
  As = c(rnorm(100,mean = cluster1.aMean, sd=cluster1.aSd ),
         rnorm(100,mean = cluster2.aMean, sd=cluster2.aSd ),
         rnorm(100,mean = cluster3.aMean, sd=cluster3.aSd )
  )
  Bs = c(rnorm(100,mean = cluster1.bMean, sd=cluster1.bSd ),
         rnorm(100,mean = cluster2.bMean, sd=cluster2.bSd ),
         rnorm(100,mean = cluster3.bMean, sd=cluster3.bSd )
  )

  generateSignalHard <- function(noise = 0.1, length=256, a = 1, b=1){
    result = noise * rnorm(length)
    x = 1:length / 10
    result <- result + sin(0.4 * a * x)  + cos(0.2 * a * x )^2 + sin(0.7 * a * x) - cos(0.5 *a * a * x) * sin(0.8 * a * x ) * cos (0.2 * a * a * a * x)
    result <- result + sin(2 * b * x)  + cos(0.8 * b * x )^3 + sin(2 * b * x) - cos(0.8 * b * x) * sin(0.8 * b * x )^2 * cos (4 * b * b * x)^6
    return(result)
  }
  plot(As,Bs,col=c(rep(1,100),rep(2,100),rep(3,100)))

  X = matrix(0,ncol=256,nrow=300)
  for(i in 1:300){
    X[i,] = generateSignalHard(a=As[i],b=Bs[i])
  }
  return(list(X=X,col=c(rep(1,100),rep(2,100),rep(3,100))))
}


#' Plot a dataset of curves in a user friendly way
#'
#' Plot the curves of the dataset \emph{X}, according to the curve's class, provided by the \emph{col} parameter
#' With the \emph{ratio} parameter, you can reduce the number of plotted curves (ex : ratio = 0.1 plot only 10% of the curves).
#'
#' @param X a matrix of time series (one time serie per line of the matrix)
#' @param col the class vector of the data (must be numbers)
#' @param ratio the ratio of curves plotted (1 : all curves plotted, 0.25 : 25\% of the curves plotted)
#' @examples
#' \dontrun{
#' dataset = generateDataset()
#' plot.curve.dataset(X=dataset$X,col=dataset$col,ratio=0.1)
#' }
#'
#' @export plot.curve.dataset
plot.curve.dataset <- function(X,col=rep(1,nrow(X)),ratio=1){
  plot(0,xlim=c(0,ncol(X)),ylim=c(min(X),max(X)))
  indices = sample(1:nrow(X),nrow(X)*ratio)
  for(i in indices){
    lines(X[i,],col=col[i])
  }
}

#' Plot the mean curves of the dataset
#'
#' For each group of curves, plot the mean of the curve, and 10 curves from the plot.
#'
#' @param X a matrix of time series (one time serie per line of the matrix)
#' @param col the class vector of the data (must be numbers)
#'
#' @export plot.mean.curve.dataset
plot.mean.curve.dataset <- function(X,col=rep(1,nrow(X))){
  classes = as.numeric(names(table(clusters)))
  ylim = c(min(X),max(X))
  for(g in classes){
    dataToPlot = X[which(clusters==g),]
    plot(colMeans(dataToPlot),col=g,type="l",lwd=3,ylim=ylim)
    title(paste("Class ",g))
    for(i in sample(1:nrow(dataToPlot), 10,replace=TRUE)){
      lines(dataToPlot[i,])
    }
  }
}


#' Perform the funHDDCwavelet algorithm on a matrix of time series
#'
#' @param X A matrix of time series (each line correspond to a time serie, each column to an time point)
#' @param K The number of classes to find
#' @param minD The minimum dimension of the subspace where a level can be represented (default : 1)
#' @param maxD The maximum dimension of the subspace where a level can be represented (default : 1)
#' @param max.level The maximum wavelet transform level used in the algorithm (default : all levels, depending of the time serie length)
#' @param dimTest The intrinsec dimension estimation method. Available : "scree", "kss", "mean", "XX\%" (ex "85\%") (default: "scree")
#' @param wavelet.family The wavelet family used for the discrete wavelet transform (see wavethresh package). Default : "DaubExPhase"
#' @param wavelet.filter.number The filter number used in for the discrete wavelet transform (see wavethresh package). Default : 4
#' @param init Type of initialization method. Available : "kmeans", "random" (default : "kmeans")
#' @param minPerClass The minimal size of the initial classes
#' @param threshold The threshold used to determine that the log-likelihood converged. Default : 0.01
#' @param minIter The minimal number of iterations of the EM algorithm
#' @param maxIter The maximal number of iterations of the EM algorithm
#' @param poolSize The number of run of the algorithm (the algorithm is executed poolSize times, and the best model is selected via BIC) Default : 5
#' @param verbose if TRUE, prompt some information of the algorithm status over time
#' @param viz if TRUE, plot the first 2 axes of the data after principal component analysis, with current cluster repartitions and colors
#'
#' @examples
#' \dontrun{
#' dataset = generateDataset()
#' X = dataset$X
#' real_cluster = dataset$col
#'
#' # Haar wavelet
#' result = funHDDCwavelet(X,K=3,minD=1,maxD=1,wavelet.family="DaubExPhase",wavelet.filter.number=1,viz=TRUE,minIter=10)
#' clusters = apply(result$tm,1,which.max)
#' adjustedRandIndex(clusters,real_cluster)
#' plot.curve.dataset(X,col=clusters,ratio=0.1)
#' }
#'
#' @export funHDDCwavelet
funHDDCwavelet = function(X, K, minD=1, maxD=1, max.level = round(log2(ncol(X)))+1,dimTest="scree",
                                   wavelet.family="DaubExPhase", wavelet.filter.number=4, init="kmeans",
                                   minPerClass=10,threshold=0.01,minIter=5,maxIter=10, poolSize=5,
                                   verbose=FALSE, viz=F){

  modelCounter = 0
  bestModel = list(BIC=NULL)
  for(tryNumber in 1:poolSize){
    for(K_current in as.vector(K)){
      for(maxD_current in as.vector(maxD)){
        result = internal.funHDDCwavelet(X,K=K_current,minD = minD, maxD = maxD_current, dimTest=dimTest,
                                                wavelet.family = wavelet.family,wavelet.filter.number=wavelet.filter.number,
                                                init=init,minPerClass = minPerClass,threshold=threshold,minIter=minIter,
                                                maxIter=maxIter, verbose=verbose,viz=viz)
        result$K <- K_current
        result$maxD <- maxD_current
        modelCounter <- modelCounter + 1
        if(is.null(bestModel$BIC) || result$BIC < bestModel$BIC){
          bestModel <- result
        }
      }
    }
  }

  cat("----------------------\n")
  cat("Best Model : \n")
  cat(paste("K =",bestModel$K,"\n"))
  cat(paste("maxD =",bestModel$maxD,"\n"))
  cat(paste("BIC =",bestModel$BIC,"\n"))
  cat(paste("Models tested :",modelCounter,"\n"))
  cat("----------------------\n")

  return(bestModel)

}



internal.funHDDCwavelet = function(X, K, minD=1, maxD=1, max.level = round(log2(ncol(X)))+1,dimTest="scree",
                          wavelet.family="DaubExPhase", wavelet.filter.number=4, init="kmeans",
                          minPerClass=10,threshold=0.01,minIter=5,maxIter=10,
                          verbose=FALSE, viz=FALSE
                 ){

  X = as.matrix(X)

  # log likelihood
  ll <- NULL

  sigma = list()





  P = ncol(X)
  N = nrow(X)
  lvls = 1:(min(max.level,log2(P)) + 1)


  tm = matrix(0,ncol=K,nrow=N)

  sigma$Ad = list()
  sigma$Qd = list()

  sigma$b = c()
  sigma$mu = matrix(0,ncol=P,nrow=K)
  sigma$Pi = c()

  #intrinsec dimension initialisation
  sigma$D = matrix(0,ncol=length(lvls),nrow=K)
  for(l in lvls){
    tmp <- getBeginEnd(l)
    maxSize = tmp$end-tmp$begin + 1
    sigma$D[,l] <-  rep(min(maxD,maxSize),K)
  }

  # clusters initialisation
  if(init=="kmeans"){
    initialClusters = kmeans(X,K)$cluster
  }
  else{
    initialClusters = sample(1:K,nrow(X),replace = TRUE)
  }

  # check if some classes are near empty
  initialClusters <- getBalancedCluster(initialClusters,K,minPerClass)

  # parameters initialisation
  for(j in 1:K){
    group <- matrix(X[which(initialClusters==j),],ncol=ncol(X))

    sigma$mu[j,] = colMeans(group)
    sigma$Pi[j] = nrow(group) / N
    V = var(group)

    dec = decomposeByBlock(V)

    sigma$Qd[[j]] = dec$Qd
    sigma$Ad[[j]] = dec$Ad
  }
  # print(sigma$Pi)

  for(iter in 1:maxIter){
    # E step
    Kc = matrix(0,ncol=K,nrow=N)
    for(j in 1:K){

      Qd = sigma$Qd[[j]]
      Ad = sigma$Ad[[j]]
      D = sigma$D[j,] # vector, dimensions for each level (for the cluster j)
      mu = sigma$mu[j,]
      Pi = sigma$Pi[j]

      for(l in lvls){
        d <- D[l]
        bounds = getBeginEnd(l)
        begin = bounds$begin
        end = bounds$end
        Ql = Qd[[l]]
        maxCol = ncol(Ql)
        lambda = Ad[[l]]
        a = Ad[[l]]
        b = 0

        reduction = d<maxCol
        if(reduction){
          Ql[,(d+1) : maxCol] <- 0
          b = mean(lambda[(d+1) : maxCol])
          if(b == 0){

            b = 0.000001
          }
          if(b<0){
            b=abs(b)
          }

          lambda[(d+1) : maxCol] <- b
        }

        Xl = as.matrix(X[,begin:end])
        mul = mu[begin:end]

        Al = tcrossprod(Ql %*% (diag(x=1/lambda,ncol=length(lambda),nrow=length(lambda))), t(Ql))

        QtQ = tcrossprod(Ql)

        for(i in 1:N){
          if(reduction){
            Kc[i,j] = Kc[i,j] + kCost(Xl[i,],mul,a,b,d,Pi,QtQ,Al)
          } else {
            Kc[i,j] = Kc[i,j] + kCostWithoutB(Xl[i,],mul,a,Pi,QtQ,Al)
          }
        }
      }
    }
    # updating tm with cost result
    for(i in 1:N){
      for(j in 1:K){
        tmp = 0

        for(l in 1:K){
          tmp = tmp + exp((Kc[i,j]-Kc[i,l])/2)
        }

        tm[i,j] <- 1/tmp
      }
    }

    # M Step

    # Pi
    sigma$Pi = colMeans(tm)

    # mu
    Nis = colSums(tm)
    sigma$mu = matrix(0,ncol=P,nrow=K)
    for(g in 1:K){
      tmp = 0
      for(i in 1:N){
        tmp = tmp + tm[i,g] * X[i,]
      }
      sigma$mu[g,] = tmp/Nis[g]
    }

    # W
    W <-  computeW(X,tm,K,sigma$mu,Nis)

    for(g in 1:K){
      res = decomposeByBlock(W[[g]])

      sigma$Ad[[g]] = res$Ad
      sigma$Qd[[g]] = res$Qd
    }

    #d
    sigma$D = sigma$D
    for(g in 1:K){
      for(l in lvls){
        a = sigma$Ad[[g]][[l]]
        tmp = getBeginEnd(l)
        levelSize = tmp$end - tmp$begin + 1
        currentClusters = apply(tm, 1, which.max)
        nbInCluster <- sum(currentClusters==g)
        intrinsecDim <- getIntrinsecDim(lambda = a, n=nbInCluster, test = dimTest)

        # condition 1 : intrinsecDim <= levelSize
        if(intrinsecDim>levelSize){
          intrinsecDim <- levelSize
        }

        # condition 2 : minD <= intrinsecDim (si minD<levelSize)
        if(minD < levelSize && intrinsecDim < minD){
          intrinsecDim <- minD
        }

        # condition 3 : intrinsecDim <= maxD
        if(intrinsecDim>maxD){
          intrinsecDim <- maxD
        }

        sigma$D[g,l] <- intrinsecDim
      }
    }

    oldLL <- ll

    # log likelihood
    ll <- computeLogLikelihood(Q_all = sigma$Qd,A_all = sigma$Ad,
                               Pi = sigma$Pi,n = nrow(X),W_full = W,D = sigma$D
    )

    nbParam = 0

    nbParam <-nbParam + K * ncol(X) + K - 1 # number of parameters for mean and pi estimation

    for(l in 1:ncol(sigma$D)){ # for each level
      nbQ <- 0 # number of parameters for Qi matrix estimation
      tmp <- getBeginEnd(l)
      p <- tmp$end - tmp$begin + 1
      for(g in 1:K){ # for each cluster
        nbQ <- nbQ + sigma$D[g,l] * (p - (sigma$D[g,l]+1)/2)
      }
      nbQ <- nbQ / K
      d_barre <- mean(sigma$D[,l])
      nbParam <- nbParam + K * (nbQ + 2 + d_barre)
    }

    bic <- -2 * ll + nbParam * log(nrow(X))
    if(is.nan(bic)){
      print("prout")
    }

    if(verbose){
      print("Pi :")
      print(sigma$Pi)
      print("D :")
      print(sigma$D)
      cat(paste("BIC =",bic,"\n"))
    }


    if(viz){
      par(mfrow=c(2,1))
      clusters = apply(tm,1,which.max)
      plot(prcomp(X)$x,col=clusters)
      title(paste("Iteration :",iter))
      barplot(sigma$Pi,col=1:K)
      title("cluster repartition")
    }
    if(!is.null(oldLL)){

      diffLL <- ll - oldLL
      if(diffLL<0 && abs(diffLL)<threshold && iter>=minIter){

        sigma$tm = tm
        sigma$BIC = bic
        sigma$nbIter = iter
        return(sigma)
      }
    }

  }
  par(mfrow=c(1,1))
  sigma$tm = tm
  sigma$BIC = bic
  sigma$nbIter = iter
  return(sigma)
}

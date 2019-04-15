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
funHDDCwavelet = function(X, K, minD=1, maxD=1, max.level = round(log2(ncol(X))),dimTest="scree",
                          wavelet.family="DaubExPhase", wavelet.filter.number=4, init="kmeans",
                          minPerClass=10,threshold=0.01,minIter=10,maxIter=50, poolSize=5,
                          verbose=FALSE, viz=F){

  modelCounter = 0
  bestModel = list(BIC=NULL)
  for(tryNumber in 1:poolSize){
    for(K_current in as.vector(K)){
      for(maxD_current in as.vector(maxD)){
        for(wavelet.filter.number_current in as.vector(wavelet.filter.number)){


          result = internal.funHDDCwavelet(X,K=K_current,minD = minD, maxD = maxD_current, dimTest=dimTest,max.level = max.level,
                                           wavelet.family = wavelet.family,wavelet.filter.number=wavelet.filter.number_current,
                                           init=init,minPerClass = minPerClass,threshold=threshold,minIter=minIter,
                                           maxIter=maxIter, verbose=verbose,viz=viz)
          result$K <- K_current
          result$maxD <- maxD_current
          result$wavelet.filter.number <- wavelet.filter.number_current
          modelCounter <- modelCounter + 1
          print(result$BIC)
          if(is.null(bestModel$BIC)){
            bestModel <- result
          } else if(result$BIC < bestModel$BIC){
            bestModel <- result
          }
        }
      }
    }
  }

  cat("----------------------\n")
  cat("Best Model : \n")
  cat(paste("K =",bestModel$K,"\n"))
  cat(paste("maxD =",bestModel$maxD,"\n"))
  cat(paste("BIC =",bestModel$BIC,"\n"))
  cat(paste("Filter =",bestModel$wavelet.filter.number,"\n"))
  cat(paste("Models tested :",modelCounter,"\n"))
  cat("----------------------\n")

  return(bestModel)

}



internal.funHDDCwavelet = function(X, K, minD=1, maxD=1, max.level = round(log2(ncol(X))),dimTest="scree",
                                   wavelet.family="DaubExPhase", wavelet.filter.number=4, init="kmeans",
                                   minPerClass=10,threshold=0.01,minIter=5,maxIter=10,
                                   verbose=FALSE, viz=FALSE
){
#  require(wavethresh)

  X = as.matrix(X)

  # clusters initialisation
  if(init=="kmeans"){
    initialClusters = kmeans(X,K)$cluster
  }
  else{
    initialClusters = sample(1:K,nrow(X),replace = TRUE)
  }


  # log likelihood
  ll <- NULL

  sigma = list()

  realMaxLevel = min(max.level,round(log2(ncol(X))))

  # transform : raw ==> wavelets
  X = t(apply(X,1,function(line){
    return(toWavelet(line,wavelet.family,wavelet.filter.number,max.level=realMaxLevel))
  }))




  P = ncol(X)
  N = nrow(X)
  print(P)
  lvls = 1:(realMaxLevel + 1)



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
          if(b <= 0){

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
    ll <- computeLogLikelihood(Q_all = sigma$Qd,
                               A_all = sigma$Ad,
                               Pi = sigma$Pi,
                               n = nrow(X),
                               W_full = W,
                               D = sigma$D
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
      print(paste("Step : ",iter))
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
      par(mfrow=c(1,1))
    }
    if(!is.null(oldLL)){

      diffLL <- ll - oldLL
      if((diffLL<0 || abs(diffLL)<threshold )&& iter>=minIter){

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

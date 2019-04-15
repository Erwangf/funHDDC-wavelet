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

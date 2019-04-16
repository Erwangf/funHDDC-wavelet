#' Plot the mean curves of the dataset
#'
#' For each group of curves, plot the mean of the curve, and 10 curves from the plot.
#'
#' @param X a matrix of time series (one time serie per line of the matrix)
#' @param col the class vector of the data (must be numbers)
#'
#' @export plot.mean.curve.dataset
plot.mean.curve.dataset <- function(X,col=rep(1,nrow(X))){
  classes = as.numeric(names(table(col)))
  ylim = c(min(X),max(X))
  for(g in classes){
    dataToPlot = X[which(col==g),]
    plot(colMeans(dataToPlot),col=g,type="l",lwd=3,ylim=ylim)
    title(paste("Class ",g))
    for(i in sample(1:nrow(dataToPlot), 10,replace=TRUE)){
      lines(dataToPlot[i,])
    }
  }
}

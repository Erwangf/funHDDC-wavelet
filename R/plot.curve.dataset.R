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

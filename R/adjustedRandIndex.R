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
  if (all(dim(tab) == c(1, 1))) {
    return(1)
  }

  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

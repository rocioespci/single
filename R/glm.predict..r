#' Computes prior probability of mutations
#'
#' This is an auxiliary function in single package. It evaluates the sigmoidal function given by the parameters slope and intercept on x.
#' @param x Numeric. Values to evaluate.
#' @param slope Slope of the sigmoidal function to evaluate.
#' @param intercept Intercept of the sigmoidal function to evaluate.
#' @return Numeric.
#' @export glm.predict.
#' @examples
#' x = c(-10:10)
#' y = glm.predict.(x,1,2)
#' plot(x,y)
glm.predict. <- function(x,slope,intercept){
      y <- 1/(1+exp( -( slope*x + intercept ) ))
      return(y)
}

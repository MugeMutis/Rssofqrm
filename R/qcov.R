#' @importFrom quantreg rq
qcov <- function(y, x, tau)
{
  n <- dim(x)[1]
  p <- dim(x)[2]

  cf <- numeric()
  for(i in 1:p){
    mod.i <- rq(y~x[,i], tau, method = "fn")
    cf[i] <- mod.i$coefficients[2]
  }

  output <- diag(cov(x)) * cf

  return(output)
}

#' @importFrom quantreg rq
kim_fun <- function(y, x, w, tau, weight){

  nvar <- dim(x)[2]
  wy <- w%*%y
  wy <- wy *weight
  y <- y *weight
  wx <- w%*%x
  w2x <- w %*% wx

  fit1 <- rq(wy~x+wx+w2x,tau=tau)
  wyhat <- fitted(fit1)
  wyhat <- wyhat *weight

  fit2 <- rq(y~x+wyhat, tau=tau)
  coefs <- fit2$coefficients
  b0 <- coefs[1]
  rho <- coefs[length(coefs)]
  b <- coefs[c(-1, -(nvar+2))]

  return(list(b0 = b0, rho = rho, b = b))
}

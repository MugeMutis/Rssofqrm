#' @importFrom quantreg rq
chernozhukov_fun <- function(y, x, w, tau, weight){

  nvar <- dim(x)[2]
  rho_candidates <- seq(-0.95, 0.95, 0.001)
  nrho <- length(rho_candidates)

  wy <- w%*% y
  wy <- wy * weight
  y <- y * weight
  wx <- w%*%x
  w2x <- w %*% wx

  fit.init <- rq(wy~x+wx+w2x, tau = 0.5)
  wyhat <- fitted(fit.init)
  wyhat <- wyhat * weight


  rho_hat <- numeric()
  for(i in 1:nrho) {
    newy <- y - rho_candidates[i]*wy
    fit <- rq(newy~x+wyhat, tau=tau)
    rho_hat[i] = fit$coef[length(fit$coef)]
  }

  rho_min <- rho_candidates[which(abs(rho_hat)==min(abs(rho_hat)))]

  newy2 <- y - rho_min*wy
  fit.final <- rq(newy2~x, tau=tau)
  coefs <- fit.final$coef

  b0 <- coefs[1]
  rho <- rho_min
  b <- coefs[-1]

  return(list(b0 = b0, rho = rho, b = b))
}

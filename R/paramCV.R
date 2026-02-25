#' @import cvTools
paramCV <- function(y, x, w, tau, gp=NULL, nbasis=NULL, model = c("KM","CH")){

  w <- norm_wei(w)
  model <- match.arg(model)
  n <- dim(x)[1]
  p <- dim(x)[2]

  if(is.null(gp))
    gp <- seq(0, 1, length.out = p)
  if(is.null(nbasis))
    nbasis <- min(10, p/4)

  h = c(1,2,3,4,5,6)

  CVmat <- as.matrix(expand.grid(h))
  CVmat <- cbind(CVmat, NA)
  colnames(CVmat) <- c("h", "BIC")

  for(i in 1:dim(CVmat)[1]){
    h <- CVmat[,1][i]

    fpqr_per <- getRFPQR(y = y, x = x, tau = tau, h = h, nbasis = nbasis, gp = gp,
                         probp1 = 0.95, hampelp2 = 0.975, hampelp3 = 0.999,
                         maxit = 1000, conv = 0.01)
    fsco <- fpqr_per$T
    V <- fpqr_per$V
    details <- fpqr_per
    weights_y <- fpqr_per$weights_y

    if(model == "KM")
      qmodel <- kim_fun(y, fsco, w, tau, weights_y)
    if(model == "CH")
      qmodel <- chernozhukov_fun(y, fsco, w, tau, weights_y)


    b0 <- qmodel$b0
    rho <- qmodel$rho
    b <- qmodel$b

    fits <- solve(diag(n) - rho*w) %*% as.matrix(rep(1,n)) * b0 +
      solve(diag(n) - rho*w) %*% fsco %*% b

    CVmat[,2][i] <- BIC_fun(y, fits, tau, h)
  }

  optima <- CVmat[which.min(CVmat[,2]),][1]
  optim_h <- optima[1]

  return(list(h = optim_h, nbasis=nbasis))
}

#' @importFrom quantreg rq
pqr <- function(y, x, tau, h)
{

  x0 <- x
  y0 <- y

  n <- dim(x)[1]
  p <- dim(x)[2]

  x <- scale(x, scale = F)
  y <- as.matrix(scale(y, scale = F))

  xw <- matrix(, p, h)
  yw <- numeric()
  xs <- matrix(, n, h)
  ys <- matrix(, n, h)
  xl <- matrix(, p, h)
  yl <- numeric()

  for(i in 1:h)
  {
    qcv <- qcov(y, x, tau)
    qsvd <- svd(qcv)

    xw. <- qsvd$u[,1]
    yw. <- qsvd$v[1,]

    buff <- svd_flip(xw., yw.)
    xw. <- buff$x
    yw. <- buff$y

    xs. <- x %*% xw.
    ys. <- y %*% yw.

    xl. <- c(t(x) %*% xs.) / c(t(xs.) %*% xs.)
    x <- x - xs. %*% xl.
    yl. <- c(t(xs.) %*% y) / c(t(xs.) %*% xs.)
    y <- y - xs. %*% yl.

    xw[,i] <- xw.
    yw[i] <- yw.
    xs[,i] <- xs.
    ys[,i] <- ys.
    xl[,i] <- xl.
    yl[i] <- yl.
  }

  xr <- xw %*% ginv(t(xl) %*% xw)
  yr <- yw %*% ginv(t(yl) %*% yw)

  fin.model <- rq(y0~xs, tau)
  pqr.cf <- as.matrix(fin.model$coefficients)
  fin.cf <- xr %*% as.matrix(fin.model$coefficients[-1])

  fits <- cbind(1, xs) %*% pqr.cf
  resds <- y0 - fits

  return(list(x = x0, y = y0, T = xs, R = xr, P = xl, W = xw, d.coef = fin.cf,
              pqr.coef = pqr.cf, fitted.values = fits, residuals = resds))
}

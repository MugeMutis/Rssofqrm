#' Robust Functional Quantile Regression Model with Spatial Autocorrelation
#'
#' Fits a spatial scalar-on-function quantile regression model using either the Kim and Muller (2004)
#' or Chernozhukov and Hansen (2006) method.
#'
#' @param y Numeric vector of length \code{n}. Response variable.
#' @param x Numeric matrix of size \code{n x p}. Functional predictor observed on a common grid.
#' @param w Numeric matrix of size \code{n x n}. Spatial weight matrix.
#' @param tau Numeric. Quantile level to estimate (between 0 and 1).
#' @param tau_comp Numeric. Quantile level to estimate for extracted components (between 0 and 1).
#' @param gp Numeric vector of length \code{p}. Grid of argument values for the functional predictor. If \code{NULL}, defaults to a regular grid between 0 and 1.
#' @param nbasis Integer. Number of basis functions to use in RFPQR. If \code{NULL}, defaults to \code{min(10, p / 4)}.
#' @param h Integer. Number of extracted components. If \code{NULL}, it is optimally determined when CV is True.
#' @param model Character string. Estimation method to use. One of \code{"KM"} (Kim et al.) or \code{"CH"} (Chernozhukov et al.).
#' @param CV Character string. Determine the number of components to extract based on the BIC. One of \code{"TRUE"} or \code{"FALSE"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{b}{Estimated coefficient vector in RFPQR space.}
#'   \item{b0}{Estimated intercept.}
#'   \item{bhat}{Estimated coefficient function in the original function space.}
#'   \item{rho}{Estimated spatial autoregressive parameter.}
#'   \item{opt.h}{Number of extracted optimal components.}
#'   \item{opt.nbasis}{Number of basis functions.}
#'   \item{fitted.values}{Fitted values from the model.}
#'   \item{residuals}{Residuals from the model.}
#'   \item{tau}{Quantile level used.}
#'   \item{details}{Output from \code{getRFPQR}, including RFPQR scores and basis functions.}
#' }
#'
#' @importFrom stats as.formula
#' @examples
#' \dontrun{
#' sim_data <- data_generation(n=250, j=101, rho=0.5, mean.e=1, sig.e=1, out.p = 0)
#' y <- sim_data$y
#' x <- sim_data$x
#' w <- sim_data$w
#' fit_kim <- rssofqr(y=y, x=x, w=w, tau=tau, model = "KM", CV="TRUE")
#' fit_ch <- rssofqr(y=y, x=x, w=w, tau=tau, model = "CH", CV="TRUE")
#' }
#'
#' @export


rssofqr <- function(y, x, w, tau, tau_comp=NULL, gp=NULL, nbasis=NULL, h=NULL, model = c("KM","CH"),CV=c("TRUE","FALSE")){

  w <- norm_wei(w)
  model <- match.arg(model)
  n <- dim(x)[1]
  p <- dim(x)[2]


  if(is.null(gp))
    gp <- seq(0, 1, length.out = p)
  if(is.null(tau_comp))
    tau_comp <- 0.5
  if(is.null(nbasis))
    nbasis <- min(10, p/4)
  if(is.null(h))
    h = 3

  if(CV == TRUE){
    optimod <- paramCV(y, x, w, tau, gp=NULL, nbasis=NULL, model)
    h <- optimod$h
    nbasis <- optimod$nbasis
  }

  fpqr_per <- getRFPQR(y = y, x = x, tau = tau_comp, h = h, nbasis = nbasis, gp = gp,
                      probp1 = 0.95, hampelp2 = 0.975, hampelp3 = 0.999,
                      maxit = 1000, conv = 0.01)

  Amat <- fpqr_per$model.details$Amat
  R <- fpqr_per$R
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

  bhat <- details$model.details$evalbase %*%
      ((t(solve(details$model.details$sinp_mat)) %*% V) %*% b)


  Ainv = solve(diag(n) - rho*w)
  fits <- Ainv %*% ((as.matrix(rep(1,n)) * b0) + (Amat %*% (R %*% as.matrix(b))))
  resids <- y - fits


  return(list(b = b, b0 = b0,
              bhat = bhat,
              rho = rho,
              opt.h=h,
              opt.nbasis=nbasis,
              fitted.values = fits,
              residuals = resids,
              tau = tau,
              details = details))
}

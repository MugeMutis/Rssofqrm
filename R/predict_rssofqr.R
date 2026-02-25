#' Predict Method for Robust Functional Quantile Regression with Spatial Effects
#'
#' Makes predictions for new observations using a fitted robust spatial functional quantile regression model.
#'
#' @param object A fitted model object returned by \code{\link{rssofqr}}.
#' @param xnew Numeric matrix of size \code{n x p}. New functional predictor data (on the same grid used in training).
#' @param wnew Numeric matrix of size \code{n x n}. Spatial weight matrix for the new observations.
#'
#' @return A numeric vector of predicted values of length \code{n}.
#' @examples
#' \dontrun{
#' sim_data <- data_generation(n=250, j=101, rho=0.5, mean.e=1, sig.e=1, out.p = 0)
#' y <- sim_data$y
#' x <- sim_data$x
#' w <- sim_data$w
#' fit_kim <- rssofqr(y=y, x=x, w=w, tau=tau, model = "KM", CV="TRUE")
#' fit_ch <- rssofqr(y=y, x=x, w=w, tau=tau, model = "CH", CV="TRUE")
#' sim_test <- data_generation(n=250, j=101, rho=0.5, mean.e=1, sig.e=1, out.p = 0)
#' y_test <- sim_test$y
#' x_test <- sim_test$x
#' w_test <- sim_test$w
#' predict_kim <- predict_rssofqr(object = fit_kim, xnew = x_test, wnew = w_test)
#' predict_ch <- predict_rssofqr(object = fit_ch, xnew = x_test, wnew = w_test)
#' }
#'
#' @export

predict_rssofqr <- function(object, xnew, wnew)
{
  n <- dim(xnew)[1]
  nbasis <- object$details$model.details$nbasis
  gp <- object$details$model.details$gp
  rho <- object$rho
  b0 <- object$b0
  b <- object$b
  R <- object$details$R

  BS.sol.test <- getAmat(data = xnew, nbf = nbasis, gp = gp)
  Amat.test <- BS.sol.test$Amat
  Amat.test <- scale_fun(data = Amat.test)

  Ainv = solve(diag(n) - rho*wnew)
  predicted.values <- Ainv %*% ((as.matrix(rep(1,n)) * b0) + (Amat.test %*% (R %*% as.matrix(b))))

  return(predicted.values)
}

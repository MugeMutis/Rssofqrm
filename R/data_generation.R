#' Data Generating Process for Spatial Scalar-on-Function Model
#'
#' Generates synthetic data from a spatial scalar-on-function regression model.
#'
#' @param n Integer. Sample size.
#' @param j Integer. Number of grid points.
#' @param rho Numeric. Spatial autoregressive parameter.
#' @param mean.e Numeric. Mean of error terms.
#' @param sig.e Numeric. Standard deviation of error Terms.
#' @param out.p Numeric. Proportion of outliers to introduce (default is 0).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{y}{Response vector of length \code{n}.}
#'   \item{x}{Matrix of functional predictor values.}
#'   \item{w}{Spatial weight matrix of size \code{n x n}.}
#'   \item{tcoefs}{True coefficient function used to generate the data.}
#'   \item{out.index}{Indices of outlier observations, if any.}
#' }
#'
#' @import fda.usc
#' @importFrom stats rnorm
#' @importFrom MASS ginv
#' @examples
#' \dontrun{
#' sim_data <- data_generation(n=250, j=101, rho=0.5, mean.e=1, sig.e=1, out.p = 0)
#' y <- sim_data$y
#' x <- sim_data$x
#' w <- sim_data$w
#' tcoef <- sim_data$tcoefs
#' out.index <- sim_data$out.index
#' }
#'
#' @export
data_generation = function(n, j, rho, mean.e, sig.e, out.p = 0){

  s = seq(0, 1, length.out = j)

  # Define cosine basis
  X.phi1 <- function(nphi, gpx) {
    phi <- matrix(0, nphi, length(gpx))
    for (j in 1:nphi) {
      phi[j, ] <- (j^-1.5) * sqrt(2) * cos(j * pi * gpx)
    }
    return(phi)
  }

  # Define sine basis
  X.phi2 <- function(nphi, gpx) {
    phi <- matrix(0, nphi, length(gpx))
    for (j in 1:nphi) {
      phi[j, ] <- (j^-1.5) * sqrt(2) * sin(j * pi * gpx)
    }
    return(phi)
  }

  # Generate one functional covariate realization
  rX.s <- function(nphi, gpx) {
    xsi <- rnorm(2 * nphi)
    phi_all <- rbind(X.phi1(nphi, gpx), X.phi2(nphi, gpx))
    X <- xsi * phi_all
    Xs <- colSums(X)
    return(Xs)
  }

  # Generate spatial weights matrix
  generate_W <- function(n) {
    wei <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          wei[i, j] <- 1 / (1 + abs(i - j))
        }
      }
    }
    W <- matrix(0, n, n)
    for (i in 1:n) {
      W[i, ] <- wei[i, ] / sum(wei[i, ])
    }
    return(W)
  }

  # Generate functional predictor X (n x length(gpx))
  fX <- t(replicate(n, rX.s(10, s)))

  vBeta = sin(2*pi * s)

  fX = fdata(fX, argvals = s)
  vBeta = fdata(vBeta, argvals = s)

  err = rnorm(n, mean=0, sd=sig.e)

  argx = inprod.fdata(fX, vBeta)

  # Compute spatial weight matrix
  W <- generate_W(n)


  fYe = ginv(diag(n) - rho*W) %*% argx + ginv(diag(n) - rho*W) %*% err
  out.index <- NULL

  if(out.p > 0){

    nout <- round(n * out.p)
    out.index <- sample(1:n, nout)
    err.out <- rnorm(n, mean=mean.e, sd=sig.e)

    # Define cosine basis
    X.phi1.out <- function(nphi, gpx) {
      phi <- matrix(0, nphi, length(gpx))
      for (j in 1:nphi) {
        phi[j, ] <- (j^-.5) * 2 * cos(j * pi * gpx)
      }
      return(phi)
    }

    # Define sine basis
    X.phi2.out <- function(nphi, gpx) {
      phi <- matrix(0, nphi, length(gpx))
      for (j in 1:nphi) {
        phi[j, ] <- (j^-.5) * 2 * sin(j * pi * gpx)
      }
      return(phi)
    }

    # Generate one functional covariate realization
    rX.s.out <- function(nphi, gpx) {
      xsi <- rnorm(2 * nphi)
      phi_all <- rbind(X.phi1.out(nphi, gpx), X.phi2.out(nphi, gpx))
      X <- xsi * phi_all
      Xs <- colSums(X)
      return(Xs)
    }

    fX.out <- t(replicate(n, rX.s.out(10, s)))
    vBeta.out = sin(2*pi * s)

    fX.out = fdata(fX.out, argvals = s)
    vBeta.out = fdata(vBeta.out, argvals = s)

    argx.out = inprod.fdata(fX.out, vBeta.out)

    fYe.out = ginv(diag(n) - rho*W) %*% argx.out + ginv(diag(n) - rho*W) %*% err.out

    fX$data[out.index,] <- fX.out$data[out.index,]
    fYe[out.index,] <- fYe.out[out.index,]
  }

  return(list("y" = fYe, "x" = fX$data, w = W, tcoefs = vBeta$data, out.index = out.index))

}

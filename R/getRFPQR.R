getRFPQR <- function(y, x, tau, h=NULL, nbasis, gp, 
                 probp1 = 0.95, hampelp2 = 0.975, hampelp3 = 0.999,
                 maxit = 1000, conv = 0.01)
{

  x0 <- x
  y0 <- y
  
  BS.sol <- getAmat(data = x, nbf = nbasis, gp = gp)
  x <- BS.sol$Amat
  n <- dim(x)[2]
  
  yx.s <- scale_fun(data = cbind(y, x))
  xs <- yx.s[,-1]
  ys <- as.matrix(yx.s[,1])
  
  
    wx <- sqrt(apply(yx.s[,-1]^2, 1, sum))
    wx <- wx / median(wx)
    wy <- abs(yx.s[,1])
    zerows <- vector(length=0)
    
    if(length(wy) / 2 > sum(wy == 0)){
      wy <- wy / (1.4826 * median(wy))
    }else{
      wy <- wy / (1.4826 * median(wy[wy != 0]))
    }
    
    probct <- qnorm(probp1)
    hampelb <- qnorm(hampelp2)
    hampelr <- qnorm(hampelp3)
    wx[which(wx <= probct)] <- 1
    wx[which(wx > probct & wx <= hampelb)] <- probct / abs(wx[which(wx > probct & wx <= hampelb)])
    wx[which(wx > hampelb & wx <= hampelr)] <- probct * (hampelr - abs(wx[which(wx > hampelb & wx <= hampelr)]))/
      (hampelr - hampelb) * 1 / abs(wx[which(wx > hampelb & wx <= hampelr)])
    wx[which(wx > hampelr)] <- 0
    wy[which(wy <= probct)] <- 1
    wy[which(wy > probct & wy <= hampelb)] <- probct/abs(wy[which(wy > probct & wy <= hampelb)])
    wy[which(wy > hampelb & wy <= hampelr)] <- probct * (hampelr - abs(wy[which(wy > hampelb & wy <= hampelr)]))/
      (hampelr - hampelb) * 1 / abs(wy[which(wy > hampelb & wy <= hampelr)])
    wy[which(wy > hampelr)] <- 0
    
    w <- wx * wy
    
    if(any(w < 1e-6)){
      w0 <- which(w < 1e-6)
      w <- replace(w, list=w0, values = 1e-6)
      we <- w
    } else {
      wxe <- wx
      wye <- wy
      we <- w
    }
    
    yx.w <- as.data.frame(yx.s) * sqrt(we)
    loops <- 1
    rold <- 10^-5
    difference <- 1
    
    while((difference > conv) && loops < maxit)
    {
      m.pqr <- pqr(y = as.matrix(yx.w[,1]), x = as.matrix(yx.w[,-1]),
                   tau = tau, h = h)
      yp <- m.pqr$fitted.values
      r <- yx.s[,1] - yp
      b <- m.pqr$d.coef
      Tpls <- as.data.frame(m.pqr$T) / sqrt(we)
      
      if (length(r) / 2 > sum(r == 0)){
        r <- abs(r) / (1.4826 * median(abs(r)))
      }else{
        r <- abs(r) / (1.4826 * median(abs(r[r != 0])))
      }
      
      dt <- scale_fun(data = Tpls)
      wtn <- sqrt(apply(dt^2, 1, sum))
      wtn <- wtn/median(wtn)
      
      probct <- qnorm(probp1)
      hampelb <- qnorm(hampelp2)
      hampelr <- qnorm(hampelp3)
      wye <- r
      wye[which(r <= probct)] <- 1
      wye[which(r > probct & r <= hampelb)] <- probct / abs(r[which(r > probct & r <= hampelb)])
      wye[which(r > hampelb & r <= hampelr)] <- probct * (hampelr-abs(r[which(r > hampelb & r <= hampelr)]))/
        (hampelr - hampelb) * 1 / abs(r[which(r > hampelb & r <= hampelr)])
      wye[which(r > hampelr)] <- 0
      wye <- as.numeric(wye)
      
      probct <- qchisq(probp1, h)
      hampelb <- qchisq(hampelp2, h)
      hampelr <- qchisq(hampelp3, h)
      wte <- wtn
      wte[which(wtn <= probct)] <- 1
      wte[which(wtn > probct & wtn <= hampelb)] <- probct / abs(wtn[which(wtn > probct & wtn <= hampelb)])
      wte[which(wtn > hampelb & wtn <= hampelr)] <- probct * (hampelr-abs(wtn[which(wtn > hampelb & wtn <= hampelr)])) /
        (hampelr - hampelb) * 1 / abs(wtn[which(wtn > hampelb & wtn <= hampelr)])
      wte[which(wtn > hampelr)] <- 0
      
      difference <- abs(sum(b^2) - rold)/rold
      rold <- sum(b^2)
      we <- wye * wte
      
      if(any(we < 1e-6)){
        w0 <- which(we < 1e-6)
        we <- replace(we, list = w0, values = 1e-6)
        zerows <- unique(c(zerows, as.numeric(names(w0))))
      }
      
      if(length(zerows) >= (n/2)){
        break
      }
      
      yx.w <- as.data.frame(yx.s) * sqrt(we)
      
      loops <- loops + 1
    }
    
    if (difference > conv){
      warning(paste("Convergence was not achieved. The scaled difference between norms of the coefficient vectors is ",
                    round(difference, digits=4)))
    }
    
    w <- we
    w[zerows] <- 0
    wt <- wte
    wt[zerows] <- 0
    wy <- wye
    wy[zerows] <- 0
  
  P <- m.pqr$P
  R <- m.pqr$R
  W <- m.pqr$W
  T <- m.pqr$T
  
 
  model.details <- list()
  model.details$Amat <- BS.sol$Amat
  model.details$bs_basis <- BS.sol$bs_basis
  model.details$inp_mat <- BS.sol$inp_mat
  model.details$sinp_mat <- BS.sol$sinp_mat
  model.details$evalbase <- BS.sol$evalbase
  model.details$nbasis <- nbasis
  model.details$gp <- gp

  
  
  return(list(V=W, T=T, R=R, P=P, model.details = model.details, weights=w, weights_y=wy, weights_x=wt))
  
}
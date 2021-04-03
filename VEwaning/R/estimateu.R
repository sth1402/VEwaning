#' @include gFunctions.R
estimateu <- function(data, 
                      times, 
                      preu, 
                      lag, 
                      gFunc,
                      theta, 
                      v, 
                      dObj) {

  n <- nrow(x = data)
  m <- length(x = theta)
  nTimes <- length(x = times)
 
  ##  Calculate dtilNu(t) and dtilYu(t) at U >=TP

  # {t_m - R_i - lag}
  tRl <- -outer(X = data[,dObj$R] + lag, Y = times, FUN = "-")
  gFuncR <- gFunction(gFunc = gFunc, u = tRl, theta = theta, knts = v)
  expgR <- exp(x = gFuncR$gu)
  
  # {t_m - E_i = lag}
  tEl <- -outer(X = data[,dObj$E] + lag, Y = times, FUN = "-")
  gFuncE <- gFunction(gFunc = gFunc, u = tEl, theta = theta, knts = v)
  expgE <- exp(x = gFuncE$gu)

  Yut <- preu$Y*{preu$wt0 * expgR + preu$wt1 * expgE}

  dNut <- preu$dNtTilde

  ##   Get Zu, ZuminusZbar, and unblinded part of estimating equation

  estu <- rep(x = 0.0, times = m)
#  ZmZbar <- array(data = 0.0, dim = c(n, nTimes, m))
  influmat <- matrix(data = 0.0, nrow = n, ncol = m)
  ZmZbar <- list()
  Yu.sum <- colSums(x = Yut)

  # {nt}
  dLambdaHat = colSums(x = dNut) / Yu.sum

  # {n x nt}
  YdLambdaHat = sweep(x = Yut, 
                      MARGIN = 2L, 
                      STATS = dLambdaHat, 
                      FUN = "*")

  # {n x nt} dN - dLambda Y
  dNmYdLambdaHat <- dNut - YdLambdaHat

  for (d in 2L:m) {
    Z <- data[,dObj$A] * gFuncE$gutheta[[ d ]] + 
         {1L - data[,dObj$A]} * gFuncR$gutheta[[ d ]]
    Zbar <- colSums(x = Yut*Z)/Yu.sum
    ZmZbar[[ d ]] <- t(x = t(x = Z) - Zbar)

    estu[d] <- sum(ZmZbar[[ d ]] * dNut)

    influmat[,d] <- rowSums(x = ZmZbar[[ d ]] * dNmYdLambdaHat)
  }

  meatu <- crossprod(x = influmat)

  dNut <- t(x = dNut)

  gradu <- matrix(data = 0.0, nrow = m, ncol = m)
  for (d1 in 2L:m) {
    for (d2 in d1:m) {
      res <- colSums(x = ZmZbar[[ d1 ]]*ZmZbar[[ d2 ]]*Yut) / Yu.sum

      gradu[d1,d2] <- sum(dNut*res)
      gradu[d2,d1] <- gradu[d1,d2]
    }
  }      

  return( list("estu" = estu, "gradu" = gradu, "meatu" = meatu))
}

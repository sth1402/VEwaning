#' @include gFunctions.R
estimateb <- function(data, 
                      times, 
                      preb, 
                      lag, 
                      theta,
                      gFunc, 
                      v, 
                      dObj) {

  n <- nrow(x = data)
  m <- length(x = theta)
  nTimes <- length(x = times)
    
  #  exponential term for tilYb(t)
  # (t_m - E_i - lag)
  tEl <- -outer(X = data[,dObj$E] + lag, Y = times, FUN = "-") 
  gFuncE <- gFunction(gFunc = gFunc, u = tEl, theta = theta, knts = v)
  expgE <- exp(x = theta[1L] + gFuncE$gu)

  Ybt <- preb$Y * {preb$wt0 + preb$wt1*expgE}

  # {n x nt}
  dNbt <- preb$dNtTilde   

  #   Get Zb, ZbminusZbar, and blinded part of estimating equation

  #{nt}
  Yb.sum <- colSums(x = Ybt)


  # {nt}
  dLambdaHat = colSums(x = dNbt) / Yb.sum

  # {n x nt}
  YdLambdaHat = sweep(x = Ybt, 
                      MARGIN = 2L, 
                      STATS = dLambdaHat, 
                      FUN = "*")

  # {n x nt} dN - dLambda Y
  dNmYdLambdaHat <- dNbt - YdLambdaHat

  estb <- rep(x = 0.0, times = m)
  influmat <- matrix(data = 0.0, nrow = n, ncol = m)
  ZbmZbar <- list()

  for (d in 1L:m) {
    # {n x nt}
    Z <- gFuncE$gutheta[[ d ]]
    Z[data[,dObj$A] == 0L,] <- 0.0

    # {nt}
    Zbar <- colSums(x = Ybt*Z) / Yb.sum

    # {n x nt}
    ZbmZbar[[ d ]] <- t(x = t(x = Z) - Zbar)

    estb[d] <- sum( ZbmZbar[[ d ]]*dNbt )

    influmat[,d] <- rowSums(x = ZbmZbar[[ d ]] * {dNmYdLambdaHat})
  }


  meatb <- crossprod(x = influmat)
  dNbt <- t(x = dNbt)

  gradb <- matrix(data = 0.0, nrow = m, ncol = m)
  for (d1 in 1L:m) {
    for (d2 in d1:m) {
      res <- colSums(x = ZbmZbar[[ d1 ]]*ZbmZbar[[ d2 ]]*Ybt) / Yb.sum
      gradb[d1,d2] <- sum(dNbt*res)
      gradb[d2,d1] <- gradb[d1,d2]
    }
  }

  return( list("estb" = estb, "gradb" = gradb, "meatb" = meatb) )
}
